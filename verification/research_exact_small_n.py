#!/usr/bin/env python3
"""Certified exhaustive computation of F(n) for small orders.

The search uses two exact invariances.  Switching by a diagonal sign matrix
first makes every edge incident with vertex 0 positive.  Permuting the other
n-1 vertices then shows that it is enough to inspect one representative of
each unlabeled graph on n-1 vertices, where graph edges encode the remaining
negative coefficients.  Brendan McKay's ``geng`` supplies exactly those
representatives.

No floating-point arithmetic is used to determine F(n).  For a negative-edge
graph G, a spin set S not containing vertex 0 has energy

    binom(n, 2) - 2|E(G)| - 2|S|(n-|S|) + 4|delta_G(S)|.

The script checks this formula against direct quadratic-form evaluation,
checks the graph6 decoder against NetworkX when available, and independently
enumerates all 2^binom(n,2) signings through n=6 without switching or graph
isomorphism reduction.

Expected exact values:

    n       2  3  4  5  6  7  8   9  10
    F(n)    1  3  4  4  5  9  10  12  13

The normal default stops at n=9.  Passing ``--max-n 10`` checks 274,668
unlabeled graphs on nine vertices and takes about one minute on an Apple M4.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import random
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


# Number of unlabeled simple graphs on k vertices, OEIS A000088.  These are
# completeness assertions for the external generator, not inputs to the
# optimization.
UNLABELED_GRAPH_COUNTS = {
    1: 1,
    2: 2,
    3: 4,
    4: 11,
    5: 34,
    6: 156,
    7: 1044,
    8: 12346,
    9: 274668,
    10: 12005168,
}

EXPECTED_F = {
    2: 1,
    3: 3,
    4: 4,
    5: 4,
    6: 5,
    7: 9,
    8: 10,
    9: 12,
    10: 13,
}

EXPECTED_OPTIMAL_SWITCHING_CLASS_COUNTS = {7: 6, 8: 2, 9: 15, 10: 2}

# These hashes commit to the complete newline-delimited graph6 stream emitted
# by nauty 2.9.3.  They are printed for auditability; --strict-stream-digests
# makes them mandatory because canonical labels could in principle differ
# across nauty versions without changing the exhaustive result.
EXPECTED_GENG_STREAM_SHA256 = {
    2: "ecf5de1a2ecc66a1876a832804c64f6b5125784e94c82285d9720621c613ab46",
    3: "b7cd2a004ade86133158ffa94292f1d79a1fa154874706bf33b9e841cd3fa4cb",
    4: "aefbaa12a956ed1f415fa897c455185134275a89a57ce1ef7d38f771c0d9129e",
    5: "b809c405cd3b8fb3cc836a9e7471658a8cf02cbc258029bddbecc9f2caf2ea32",
    6: "db56d888b80afedd817dc7ee048204866e01a1044bb6090156a0285ef4ffa67e",
    7: "5ee3aab11e44b22b9f71922c7e51b5fe050a0de8ee525e08754233464c1a4172",
    8: "3c9d236c155206869769a6fdec67a879f31a8e1b5b772d765d384449fa6b5882",
    9: "6b740e1c1ec4f6c7d5539e2e236da0f1ad6aa3120d534590b0ea1f09ddc0b345",
    10: "ce9c5d4d27c8e55de5f0c6348ec781a650382e16bdff26b6c3418fa00a9cfcf9",
}


@dataclass(frozen=True)
class SearchResult:
    n: int
    value: int
    graph6: str
    maximizing_masks: tuple[int, ...]
    representative_energy_spectrum: tuple[int, ...]
    optimal_representatives: int
    optimal_representatives_sha256: str
    optimal_graph6: tuple[str, ...]
    graphs_checked: int
    generator_sha256: str


def locate_geng(explicit: str | None) -> str:
    """Locate nauty's geng without changing PATH."""
    candidates = [
        explicit,
        shutil.which("geng"),
        "/opt/homebrew/opt/nauty/bin/geng",
        "/usr/local/bin/geng",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return candidate
    raise FileNotFoundError(
        "nauty geng was not found; install nauty or pass --geng /absolute/path"
    )


def graph6_adjacency(line: bytes, expected_order: int) -> tuple[int, ...]:
    """Decode a small graph6 record (orders at most 62) into bit adjacencies."""
    record = line.strip()
    if not record or record.startswith(b">>"):
        raise ValueError(f"unexpected graph6 record: {record!r}")
    order = record[0] - 63
    if order != expected_order or not 0 <= order <= 62:
        raise ValueError((order, expected_order, record))

    bits: list[int] = []
    for byte in record[1:]:
        value = byte - 63
        if not 0 <= value < 64:
            raise ValueError(record)
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))

    needed = order * (order - 1) // 2
    if len(bits) < needed or any(bits[needed:]):
        raise ValueError(("noncanonical graph6 padding", record))

    adjacency = [0] * order
    index = 0
    for column in range(1, order):
        for row in range(column):
            if bits[index]:
                adjacency[row] |= 1 << column
                adjacency[column] |= 1 << row
            index += 1
    return tuple(adjacency)


def edge_count(adjacency: tuple[int, ...]) -> int:
    return sum(mask.bit_count() for mask in adjacency) // 2


def cut_size(adjacency: tuple[int, ...], subset: int) -> int:
    return sum(
        (adjacency[v] & ~subset).bit_count()
        for v in range(len(adjacency))
        if subset >> v & 1
    )


def energy_from_cut(adjacency: tuple[int, ...], subset: int) -> int:
    """Energy when root spin is +1 and precisely ``subset`` has spin -1."""
    residual_order = len(adjacency)
    n = residual_order + 1
    total_edges = n * (n - 1) // 2
    subset_size = subset.bit_count()
    return (
        total_edges
        - 2 * edge_count(adjacency)
        - 2 * subset_size * (n - subset_size)
        + 4 * cut_size(adjacency, subset)
    )


def direct_energy(adjacency: tuple[int, ...], subset: int) -> int:
    """Independent O(n^2) evaluation of the signed quadratic form."""
    spins = [1] + [-1 if subset >> v & 1 else 1 for v in range(len(adjacency))]
    value = sum(spins[0] * spins[v + 1] for v in range(len(adjacency)))
    for row in range(len(adjacency)):
        for column in range(row + 1, len(adjacency)):
            coefficient = -1 if adjacency[row] >> column & 1 else 1
            value += coefficient * spins[row + 1] * spins[column + 1]
    return value


def maximum_absolute_energy(
    adjacency: tuple[int, ...], stop_at: int | None = None
) -> tuple[int, tuple[int, ...]]:
    """Return M(A), optionally stopping once the given lower bound is reached."""
    maximum = -1
    masks: list[int] = []
    for subset in range(1 << len(adjacency)):
        value = abs(energy_from_cut(adjacency, subset))
        if value > maximum:
            maximum = value
            masks = [subset]
        elif value == maximum:
            masks.append(subset)
        if stop_at is not None and maximum >= stop_at:
            return maximum, tuple(masks)
    return maximum, tuple(masks)


def geng_records(geng: str, order: int) -> Iterator[bytes]:
    process = subprocess.Popen(
        [geng, "-q", str(order)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout is not None
    try:
        for line in process.stdout:
            if line.strip():
                yield line.strip()
    finally:
        stderr = (
            process.stderr.read().decode("utf-8", errors="replace")
            if process.stderr
            else ""
        )
        return_code = process.wait()
        if return_code:
            raise RuntimeError(f"geng exited {return_code}: {stderr}")


def exact_search(n: int, geng: str, networkx_crosscheck: bool) -> SearchResult:
    residual_order = n - 1
    expected_count = UNLABELED_GRAPH_COUNTS[residual_order]
    digest = hashlib.sha256()
    count = 0
    best = n * (n - 1) // 2 + 1
    best_record = ""
    best_masks: tuple[int, ...] = ()
    best_spectrum: tuple[int, ...] = ()
    optimal_count = 0
    optimal_digest = hashlib.sha256()
    optimal_records: list[str] = []

    nx = None
    if networkx_crosscheck:
        try:
            import networkx as nx_module

            nx = nx_module
        except ImportError as error:
            raise RuntimeError("--networkx-crosscheck requires networkx") from error

    for raw_record in geng_records(geng, residual_order):
        digest.update(raw_record + b"\n")
        adjacency = graph6_adjacency(raw_record, residual_order)
        count += 1

        # A deterministic sparse sample catches a second graph6-decoder error
        # without making NetworkX part of the proof search.
        if nx is not None and (count <= 3 or count % 997 == 0):
            graph = nx.from_graph6_bytes(raw_record)
            nx_edges = {
                (min(row, column), max(row, column)) for row, column in graph.edges()
            }
            decoded_edges = {
                (row, column)
                for row in range(residual_order)
                for column in range(row + 1, residual_order)
                if adjacency[row] >> column & 1
            }
            if nx_edges != decoded_edges:
                raise AssertionError((raw_record, nx_edges, decoded_edges))

        # Energies all have the parity of binom(n,2), so a value larger than
        # ``best`` is at least ``best + 2``.  Stopping at best+1 rejects such a
        # graph early while still fully evaluating and counting every graph
        # that actually attains the optimum.
        maximum, masks = maximum_absolute_energy(adjacency, stop_at=best + 1)
        if maximum < best:
            # Recompute without early termination and save a full exact
            # spectrum; this is the upper-bound certificate.
            maximum, masks = maximum_absolute_energy(adjacency)
            spectrum = tuple(
                sorted(
                    energy_from_cut(adjacency, subset)
                    for subset in range(1 << residual_order)
                )
            )
            if maximum != max(abs(value) for value in spectrum):
                raise AssertionError((n, raw_record, maximum))
            best = maximum
            best_record = raw_record.decode("ascii")
            best_masks = masks
            best_spectrum = spectrum
            optimal_count = 1
            optimal_digest = hashlib.sha256(raw_record + b"\n")
            optimal_records = [raw_record.decode("ascii")]
        elif maximum == best:
            optimal_count += 1
            optimal_digest.update(raw_record + b"\n")
            optimal_records.append(raw_record.decode("ascii"))

    if count != expected_count:
        raise AssertionError(
            ("incomplete geng stream", residual_order, count, expected_count)
        )
    if not best_record:
        raise AssertionError(("no signing was evaluated", n))

    adjacency = graph6_adjacency(best_record.encode("ascii"), residual_order)
    for subset in range(1 << residual_order):
        formula_value = energy_from_cut(adjacency, subset)
        direct_value = direct_energy(adjacency, subset)
        if formula_value != direct_value:
            raise AssertionError((n, best_record, subset, formula_value, direct_value))

    return SearchResult(
        n=n,
        value=best,
        graph6=best_record,
        maximizing_masks=best_masks,
        representative_energy_spectrum=best_spectrum,
        optimal_representatives=optimal_count,
        optimal_representatives_sha256=optimal_digest.hexdigest(),
        optimal_graph6=tuple(optimal_records),
        graphs_checked=count,
        generator_sha256=digest.hexdigest(),
    )


def adjacency_from_labeled_edge_mask(order: int, mask: int) -> tuple[int, ...]:
    adjacency = [0] * order
    bit = 0
    for row in range(order):
        for column in range(row + 1, order):
            if mask >> bit & 1:
                adjacency[row] |= 1 << column
                adjacency[column] |= 1 << row
            bit += 1
    return tuple(adjacency)


def unreduced_direct_search(n: int) -> int:
    """Independent search over every signing, using only the definition."""
    edges = tuple(itertools.combinations(range(n), 2))
    best = len(edges) + 1
    for coefficient_mask in range(1 << len(edges)):
        maximum = -1
        for spin_mask in range(1 << (n - 1)):
            spins = (1,) + tuple(
                -1 if spin_mask >> (vertex - 1) & 1 else 1 for vertex in range(1, n)
            )
            value = 0
            for edge_index, (row, column) in enumerate(edges):
                coefficient = -1 if coefficient_mask >> edge_index & 1 else 1
                value += coefficient * spins[row] * spins[column]
            maximum = max(maximum, abs(value))
            if maximum >= best:
                break
        best = min(best, maximum)
    return best


def signed_matrix(adjacency: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    n = len(adjacency) + 1
    matrix = [[0] * n for _ in range(n)]
    for vertex in range(1, n):
        matrix[0][vertex] = matrix[vertex][0] = 1
    for row in range(n - 1):
        for column in range(row + 1, n - 1):
            value = -1 if adjacency[row] >> column & 1 else 1
            matrix[row + 1][column + 1] = value
            matrix[column + 1][row + 1] = value
    return tuple(tuple(row) for row in matrix)


def matrix_product(
    left: tuple[tuple[int, ...], ...], right: tuple[tuple[int, ...], ...]
) -> tuple[tuple[int, ...], ...]:
    order = len(left)
    return tuple(
        tuple(
            sum(left[row][index] * right[index][column] for index in range(order))
            for column in range(order)
        )
        for row in range(order)
    )


def characteristic_polynomial_coefficients(
    matrix: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    """Exact Faddeev-LeVerrier coefficients of det(lambda I - A)."""
    order = len(matrix)
    identity = tuple(
        tuple(1 if row == column else 0 for column in range(order))
        for row in range(order)
    )
    auxiliary = identity
    coefficients = [1]
    for degree in range(1, order + 1):
        product = matrix_product(matrix, auxiliary)
        negative_trace = -sum(product[index][index] for index in range(order))
        quotient, remainder = divmod(negative_trace, degree)
        if remainder:
            raise AssertionError(("nonintegral characteristic coefficient", degree))
        coefficients.append(quotient)
        auxiliary = tuple(
            tuple(
                product[row][column] + (quotient if row == column else 0)
                for column in range(order)
            )
            for row in range(order)
        )
    return tuple(coefficients)


def optimizer_invariants(result: SearchResult) -> dict[str, object]:
    adjacency = graph6_adjacency(result.graph6.encode("ascii"), result.n - 1)
    matrix = signed_matrix(adjacency)
    square = matrix_product(matrix, matrix)
    conference = all(
        square[row][column] == (result.n - 1 if row == column else 0)
        for row in range(result.n)
        for column in range(result.n)
    )
    coefficients = characteristic_polynomial_coefficients(matrix)
    histogram: dict[str, int] = {}
    for value in result.representative_energy_spectrum:
        key = str(value)
        histogram[key] = histogram.get(key, 0) + 1
    return {
        "optimizer_graph6": result.graph6,
        "negative_residual_edges": edge_count(adjacency),
        "negative_residual_degree_sequence": sorted(
            mask.bit_count() for mask in adjacency
        ),
        "signed_row_sums": sorted(sum(row) for row in matrix),
        "characteristic_polynomial_coefficients": list(coefficients),
        "signed_determinant": (-1) ** result.n * coefficients[-1],
        "conference_identity": conference,
        "energy_histogram": histogram,
    }


def optimizer_catalogue_invariants(result: SearchResult) -> dict[str, object]:
    characteristic_polynomials: set[tuple[int, ...]] = set()
    determinants: set[int] = set()
    negative_edge_counts: set[int] = set()
    conference_records: list[str] = []
    for record in result.optimal_graph6:
        adjacency = graph6_adjacency(record.encode("ascii"), result.n - 1)
        matrix = signed_matrix(adjacency)
        square = matrix_product(matrix, matrix)
        coefficients = characteristic_polynomial_coefficients(matrix)
        characteristic_polynomials.add(coefficients)
        determinants.add((-1) ** result.n * coefficients[-1])
        negative_edge_counts.add(edge_count(adjacency))
        if all(
            square[row][column] == (result.n - 1 if row == column else 0)
            for row in range(result.n)
            for column in range(result.n)
        ):
            conference_records.append(record)
    return {
        "conference_representatives": conference_records,
        "distinct_characteristic_polynomial_coefficients": [
            list(coefficients) for coefficients in sorted(characteristic_polynomials)
        ],
        "signed_determinants": sorted(determinants),
        "negative_residual_edge_counts": sorted(negative_edge_counts),
    }


def normalized_graph_at_root(
    matrix: tuple[tuple[int, ...], ...], root: int
) -> tuple[int, ...]:
    remaining = [vertex for vertex in range(len(matrix)) if vertex != root]
    adjacency = [0] * len(remaining)
    for row_index, row in enumerate(remaining):
        for column_index in range(row_index + 1, len(remaining)):
            column = remaining[column_index]
            normalized_coefficient = (
                matrix[root][row] * matrix[row][column] * matrix[root][column]
            )
            if normalized_coefficient == -1:
                adjacency[row_index] |= 1 << column_index
                adjacency[column_index] |= 1 << row_index
    return tuple(adjacency)


def optimizer_switching_classes(result: SearchResult) -> list[list[str]]:
    """Classify optimal root-normalized graphs up to switching and permutation."""
    try:
        import networkx as nx
    except ImportError as error:
        raise RuntimeError("--classify-switching-optima requires networkx") from error

    families: list[tuple[str, list[object]]] = []
    for record in result.optimal_graph6:
        adjacency = graph6_adjacency(record.encode("ascii"), result.n - 1)
        matrix = signed_matrix(adjacency)
        graphs: list[object] = []
        for root in range(result.n):
            rooted = normalized_graph_at_root(matrix, root)
            graph = nx.Graph()
            graph.add_nodes_from(range(result.n - 1))
            graph.add_edges_from(
                (row, column)
                for row in range(result.n - 1)
                for column in range(row + 1, result.n - 1)
                if rooted[row] >> column & 1
            )
            graphs.append(graph)
        families.append((record, graphs))

    classes: list[list[tuple[str, list[object]]]] = []
    for record, family in families:
        for equivalence_class in classes:
            representative_family = equivalence_class[0][1]
            equivalent = any(
                sorted(dict(left.degree()).values())
                == sorted(dict(right.degree()).values())
                and nx.is_isomorphic(left, right)
                for left in family
                for right in representative_family
            )
            if equivalent:
                equivalence_class.append((record, family))
                break
        else:
            classes.append([(record, family)])

    output = [
        [record for record, _ in equivalence_class] for equivalence_class in classes
    ]
    if sorted(record for group in output for record in group) != sorted(
        result.optimal_graph6
    ):
        raise AssertionError(
            ("switching classification lost a representative", result.n)
        )
    expected = EXPECTED_OPTIMAL_SWITCHING_CLASS_COUNTS.get(result.n)
    if expected is not None and len(output) != expected:
        raise AssertionError(("switching class count", result.n, len(output), expected))
    return output


def verify_gf9_paley_conference() -> dict[str, object]:
    """Construct the order-10 Paley matrix over F_3[t]/(t^2+1)."""
    field = tuple(itertools.product(range(3), repeat=2))

    def subtract(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        return ((left[0] - right[0]) % 3, (left[1] - right[1]) % 3)

    def multiply(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        # t^2 = -1 in F_3[t]/(t^2+1).
        return (
            (left[0] * right[0] - left[1] * right[1]) % 3,
            (left[0] * right[1] + left[1] * right[0]) % 3,
        )

    nonzero_squares = {multiply(value, value) for value in field if value != (0, 0)}

    def quadratic_character(value: tuple[int, int]) -> int:
        if value == (0, 0):
            return 0
        return 1 if value in nonzero_squares else -1

    order = 10
    matrix = [[0] * order for _ in range(order)]
    for vertex in range(1, order):
        matrix[0][vertex] = matrix[vertex][0] = 1
    for row, left in enumerate(field, start=1):
        for column, right in enumerate(field, start=1):
            if row != column:
                matrix[row][column] = quadratic_character(subtract(left, right))
    immutable = tuple(tuple(row) for row in matrix)
    if immutable != tuple(
        tuple(matrix[column][row] for column in range(order)) for row in range(order)
    ):
        raise AssertionError("the GF(9) Paley matrix is not symmetric")
    square = matrix_product(immutable, immutable)
    expected_square = tuple(
        tuple(9 if row == column else 0 for column in range(order))
        for row in range(order)
    )
    if square != expected_square:
        raise AssertionError("the GF(9) Paley conference identity failed")

    histogram: dict[int, int] = {}
    maximum = 0
    for spin_mask in range(1 << (order - 1)):
        spins = (1,) + tuple(
            -1 if spin_mask >> (vertex - 1) & 1 else 1 for vertex in range(1, order)
        )
        energy = sum(
            immutable[row][column] * spins[row] * spins[column]
            for row in range(order)
            for column in range(row + 1, order)
        )
        histogram[energy] = histogram.get(energy, 0) + 1
        maximum = max(maximum, abs(energy))
    if maximum != 15 or sum(histogram.values()) != 1 << (order - 1):
        raise AssertionError(("GF(9) Paley Boolean maximum", maximum, histogram))
    return {
        "order": order,
        "maximum_absolute_energy": maximum,
        "energy_histogram": {str(key): histogram[key] for key in sorted(histogram)},
        "conference_identity": "PASSED",
    }


def verify_switching_invariance(seed: int = 413935) -> None:
    generator = random.Random(seed)
    for n in range(3, 9):
        residual_order = n - 1
        edge_total = residual_order * (residual_order - 1) // 2
        for _ in range(25):
            adjacency = adjacency_from_labeled_edge_mask(
                residual_order, generator.randrange(1 << edge_total)
            )
            original = max(
                abs(direct_energy(adjacency, subset))
                for subset in range(1 << residual_order)
            )

            # Relabel all n vertices and switch by arbitrary vertex signs.
            coefficients = [[0] * n for _ in range(n)]
            for vertex in range(1, n):
                coefficients[0][vertex] = coefficients[vertex][0] = 1
            for row in range(residual_order):
                for column in range(row + 1, residual_order):
                    value = -1 if adjacency[row] >> column & 1 else 1
                    coefficients[row + 1][column + 1] = value
                    coefficients[column + 1][row + 1] = value
            permutation = list(range(n))
            generator.shuffle(permutation)
            switches = [generator.choice((-1, 1)) for _ in range(n)]
            transformed = [
                [
                    switches[row]
                    * switches[column]
                    * coefficients[permutation[row]][permutation[column]]
                    for column in range(n)
                ]
                for row in range(n)
            ]
            transformed_maximum = 0
            for subset in range(1 << (n - 1)):
                spins = [1] + [-1 if subset >> v & 1 else 1 for v in range(n - 1)]
                value = sum(
                    transformed[row][column] * spins[row] * spins[column]
                    for row in range(n)
                    for column in range(row + 1, n)
                )
                transformed_maximum = max(transformed_maximum, abs(value))
            if transformed_maximum != original:
                raise AssertionError((n, original, transformed_maximum))


def verify_corruption_controls(results: Iterable[SearchResult]) -> None:
    checked = list(results)
    if not checked:
        raise AssertionError("no results supplied to corruption controls")

    # Removing absolute values is a tempting but invalid simplification.  The
    # frustrated triangle (one residual negative edge) has energies
    # (-3, 1, 1, 1): its one-sided and absolute maxima differ.
    frustrated_triangle = (0b10, 0b01)
    energies = [energy_from_cut(frustrated_triangle, subset) for subset in range(4)]
    if max(energies) == max(abs(value) for value in energies):
        raise AssertionError("dropping absolute values was not detected")

    # Flipping a stored spectrum entry beyond the certified maximum must be
    # rejected by the same certificate predicate.
    result = checked[-1]
    corrupted = result.representative_energy_spectrum[:-1] + (result.value + 2,)
    if max(abs(value) for value in corrupted) == result.value:
        raise AssertionError("corrupted spectrum was not detected")


def result_as_json(
    result: SearchResult,
    emit_all_optima: bool,
    classify_switching_optima: bool,
) -> dict[str, object]:
    output: dict[str, object] = {
        "n": result.n,
        "F": result.value,
        "normalized": result.value / result.n**1.5,
        "optimizer_graph6": result.graph6,
        "maximizing_masks": list(result.maximizing_masks),
        "optimizer_energy_spectrum_sha256": hashlib.sha256(
            json.dumps(
                result.representative_energy_spectrum, separators=(",", ":")
            ).encode()
        ).hexdigest(),
        "optimal_unlabeled_rooted_representatives": result.optimal_representatives,
        "optimal_representatives_sha256": result.optimal_representatives_sha256,
        "graphs_checked": result.graphs_checked,
        "geng_stream_sha256": result.generator_sha256,
    }
    if result.n >= 7:
        output["optimizer_invariants"] = optimizer_invariants(result)
        output["optimizer_catalogue_invariants"] = optimizer_catalogue_invariants(
            result
        )
    if emit_all_optima:
        output["all_optimal_graph6"] = list(result.optimal_graph6)
    if classify_switching_optima:
        output["optimal_switching_classes"] = optimizer_switching_classes(result)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=2)
    parser.add_argument("--max-n", type=int, default=9)
    parser.add_argument("--geng")
    parser.add_argument("--networkx-crosscheck", action="store_true")
    parser.add_argument("--emit-all-optima", action="store_true")
    parser.add_argument("--classify-switching-optima", action="store_true")
    parser.add_argument("--strict-stream-digests", action="store_true")
    parser.add_argument(
        "--labeled-crosscheck-max-n",
        type=int,
        default=6,
        help="independently enumerate every unreduced signing through this n",
    )
    arguments = parser.parse_args()
    if not 2 <= arguments.min_n <= arguments.max_n <= 11:
        raise SystemExit("require 2 <= min-n <= max-n <= 11")

    geng = locate_geng(arguments.geng)
    verify_switching_invariance()
    results: list[SearchResult] = []
    for n in range(arguments.min_n, arguments.max_n + 1):
        result = exact_search(n, geng, arguments.networkx_crosscheck)
        if n <= arguments.labeled_crosscheck_max_n:
            independent = unreduced_direct_search(n)
            if independent != result.value:
                raise AssertionError(
                    ("unreduced crosscheck", n, independent, result.value)
                )
        expected = EXPECTED_F.get(n)
        if expected is not None and result.value != expected:
            raise AssertionError(("expected F(n)", n, result.value, expected))
        if (
            arguments.strict_stream_digests
            and result.generator_sha256 != EXPECTED_GENG_STREAM_SHA256[n]
        ):
            raise AssertionError(
                (
                    "unexpected geng stream",
                    n,
                    result.generator_sha256,
                    EXPECTED_GENG_STREAM_SHA256[n],
                )
            )
        results.append(result)
        print(
            json.dumps(
                result_as_json(
                    result,
                    arguments.emit_all_optima,
                    arguments.classify_switching_optima,
                ),
                sort_keys=True,
            ),
            flush=True,
        )

    verify_corruption_controls(results)
    paley_certificate = verify_gf9_paley_conference()
    print(
        json.dumps(
            {
                "corruption_controls": "PASSED",
                "deterministic_seed": 413935,
                "exact_integer_arithmetic": True,
                "expected_F_values": "PASSED",
                "independent_unreduced_search_through_n": min(
                    arguments.max_n, arguments.labeled_crosscheck_max_n
                ),
                "paley_gf9_conference": paley_certificate,
                "switching_invariance": "PASSED",
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Complete order-9 weighted Bellman and two-point geometry experiment.

Nauty ``geng`` supplies one representative of each unlabeled graph on eight
vertices.  Root switching turns these 12,346 records into a complete set of
root-normalized order-9 signings.  All energy, distance, covering-radius, and
extension calculations below use exact integer NumPy arrays.

The completeness boundary is explicit: the script trusts nauty's generator,
while asserting the known record count and the repository's committed stream
digest.  NetworkX independently decodes a deterministic sample and performs
the switching-permutation classification of the optimum records.
"""

from __future__ import annotations

import hashlib
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import research_exact_small_n as small  # noqa: E402

EXPECTED_PAIR_COUNTS = {
    (12, 0): 20,
    (12, 1): 35,
    (14, 0): 452,
    (14, 1): 612,
    (16, 0): 3049,
    (16, 1): 36,
    (18, 0): 3676,
    (18, 1): 2,
    (20, 0): 2456,
    (22, 0): 1204,
    (24, 0): 496,
    (26, 0): 190,
    (28, 0): 74,
    (30, 0): 28,
    (32, 0): 10,
    (34, 0): 4,
    (36, 0): 2,
}
EXPECTED_COLLISION = {
    "GHOgmo": (4, 4, 15),
    "Gxd?Dc": (3, 3, 17),
}
EXPECTED_COLLISION_HISTOGRAM = ((2, 124), (6, 85), (10, 37), (14, 10))
EXPECTED_COLLISION_PAIR_DISTANCE = (10, 16, 14, 20, 40)


def projective_spins(order: int) -> np.ndarray:
    masks = np.arange(1 << (order - 1), dtype=np.uint16)
    output = np.ones((len(masks), order), dtype=np.int16)
    for index in range(order - 1):
        output[:, index + 1] = 1 - 2 * ((masks >> index) & 1)
    return output


def coefficients(adjacency: tuple[int, ...], order: int) -> np.ndarray:
    values: list[int] = []
    for row in range(order):
        for column in range(row + 1, order):
            if row == 0:
                values.append(1)
            else:
                values.append(-1 if adjacency[row - 1] >> (column - 1) & 1 else 1)
    return np.asarray(values, dtype=np.int16)


def sparse_key(counts: np.ndarray) -> tuple[tuple[int, int], ...]:
    nonzero = np.flatnonzero(counts)
    return tuple((int(index), int(counts[index])) for index in nonzero)


def group_add(
    groups: dict[object, tuple[set[int], int]], key: object, deficit: int
) -> None:
    if key not in groups:
        groups[key] = ({deficit}, 1)
        return
    deficits, count = groups[key]
    deficits.add(deficit)
    groups[key] = (deficits, count + 1)


def mixed_summary(
    groups: dict[object, tuple[set[int], int]],
) -> tuple[int, int, list[int]]:
    sizes = sorted(
        (count for deficits, count in groups.values() if len(deficits) > 1),
        reverse=True,
    )
    return len(sizes), sum(sizes), sizes


def record_geometry(
    adjacency: tuple[int, ...],
    pair_products: np.ndarray,
    distance: np.ndarray,
    dot_absolute: np.ndarray,
) -> dict[str, object]:
    coefficient_vector = coefficients(adjacency, 9)
    energies = np.abs(pair_products @ coefficient_vector).astype(np.int16)
    maximum = int(energies.max())
    if np.any((maximum - energies) % 2):
        raise AssertionError("quadratic energy deficits are not even")
    weights = ((maximum - energies) // 2).astype(np.int16)
    extremizers = np.flatnonzero(energies == maximum)

    rho_ext = int(distance[:, extremizers].min(axis=1).max())
    rho_weighted = int((distance + weights[None, :]).min(axis=1).max())
    extension_from_radius = maximum + 9 - 2 * rho_weighted
    extension_direct = int((energies[None, :] + dot_absolute).max(axis=1).min())
    if extension_direct != extension_from_radius:
        raise AssertionError(
            ("weighted Bellman identity", extension_direct, extension_from_radius)
        )

    histogram = sparse_key(np.bincount(energies, minlength=37))
    maximizer_pair_distance = tuple(
        int(value)
        for value in np.bincount(
            distance[np.ix_(extremizers, extremizers)].ravel(), minlength=5
        )
    )
    energy_indices = (energies // 2).astype(np.int16)
    colored_codes = (
        ((energy_indices[:, None] * 19 + energy_indices[None, :]) * 5 + distance)
        .astype(np.int32)
        .ravel()
    )
    colored_two_point = sparse_key(np.bincount(colored_codes, minlength=19 * 19 * 5))
    return {
        "maximum": maximum,
        "deficit": 4 - rho_weighted,
        "rho_ext": rho_ext,
        "rho_weighted": rho_weighted,
        "extension": extension_direct,
        "histogram": histogram,
        "maximizer_pair_distance": maximizer_pair_distance,
        "colored_two_point": colored_two_point,
    }


def networkx_decode_check(raw_record: bytes, adjacency: tuple[int, ...]) -> None:
    try:
        import networkx as nx
    except ImportError as error:
        raise RuntimeError("this verifier requires NetworkX") from error
    graph = nx.from_graph6_bytes(raw_record)
    decoded = {
        (row, column)
        for row in range(8)
        for column in range(row + 1, 8)
        if adjacency[row] >> column & 1
    }
    actual = {(min(row, column), max(row, column)) for row, column in graph.edges()}
    if actual != decoded:
        raise AssertionError(("independent graph6 decoder", raw_record))


def main() -> None:
    geng = small.locate_geng(None)
    states = projective_spins(9)
    edges = tuple((row, column) for row in range(9) for column in range(row + 1, 9))
    pair_products = np.column_stack(
        [states[:, row] * states[:, column] for row, column in edges]
    ).astype(np.int16)
    dot_absolute = np.abs(states @ states.T).astype(np.int16)
    distance = ((9 - dot_absolute) // 2).astype(np.int16)
    if not np.array_equal(dot_absolute, 9 - 2 * distance):
        raise AssertionError("projective distance normalization failed")

    stream_digest = hashlib.sha256()
    pair_counts: Counter[tuple[int, int]] = Counter()
    histogram_groups: dict[object, tuple[set[int], int]] = {}
    maximizer_groups: dict[object, tuple[set[int], int]] = {}
    colored_groups: dict[object, tuple[set[int], int]] = {}
    optimum_records: list[str] = []
    record_deficits: dict[str, int] = {}
    count = 0

    for raw_record in small.geng_records(geng, 8):
        count += 1
        stream_digest.update(raw_record + b"\n")
        adjacency = small.graph6_adjacency(raw_record, 8)
        if count <= 3 or count % 997 == 0:
            networkx_decode_check(raw_record, adjacency)
        geometry = record_geometry(adjacency, pair_products, distance, dot_absolute)
        maximum = int(geometry["maximum"])
        deficit = int(geometry["deficit"])
        pair_counts[(maximum, deficit)] += 1
        histogram = geometry["histogram"]
        maximizer_key = (histogram, geometry["maximizer_pair_distance"])
        group_add(histogram_groups, histogram, deficit)
        group_add(maximizer_groups, maximizer_key, deficit)
        group_add(colored_groups, geometry["colored_two_point"], deficit)
        if maximum == 12:
            record = raw_record.decode("ascii")
            optimum_records.append(record)
            record_deficits[record] = deficit

    if count != small.UNLABELED_GRAPH_COUNTS[8]:
        raise AssertionError(("incomplete order-9 catalogue", count))
    expected_digest = small.EXPECTED_GENG_STREAM_SHA256[9]
    if stream_digest.hexdigest() != expected_digest:
        raise AssertionError(
            ("order-9 stream digest", stream_digest.hexdigest(), expected_digest)
        )
    if dict(pair_counts) != EXPECTED_PAIR_COUNTS:
        raise AssertionError(("full order-9 Bellman pairs", dict(pair_counts)))

    pareto = {
        pair
        for pair in pair_counts
        if not any(
            other != pair and other[0] <= pair[0] and other[1] <= pair[1]
            for other in pair_counts
        )
    }
    if pareto != {(12, 0)}:
        raise AssertionError(("order-9 Pareto frontier", pareto))

    dummy_result = small.SearchResult(
        n=9,
        value=12,
        graph6=optimum_records[0],
        maximizing_masks=(),
        representative_energy_spectrum=(),
        optimal_representatives=len(optimum_records),
        optimal_representatives_sha256="",
        optimal_graph6=tuple(optimum_records),
        graphs_checked=count,
        generator_sha256=stream_digest.hexdigest(),
    )
    switching_classes = small.optimizer_switching_classes(dummy_result)
    class_deficits: Counter[int] = Counter()
    for equivalence_class in switching_classes:
        deficits = {record_deficits[record] for record in equivalence_class}
        if len(deficits) != 1:
            raise AssertionError(("weighted deficit not switching-invariant", deficits))
        class_deficits[next(iter(deficits))] += 1
    if dict(class_deficits) != {0: 4, 1: 11}:
        raise AssertionError(("order-9 class deficits", dict(class_deficits)))

    histogram_summary = mixed_summary(histogram_groups)
    maximizer_summary = mixed_summary(maximizer_groups)
    colored_summary = mixed_summary(colored_groups)
    if histogram_summary != (10, 874, [176, 154, 144, 108, 84, 80, 58, 24, 24, 22]):
        raise AssertionError(("histogram collision summary", histogram_summary))
    if maximizer_summary != (3, 112, [50, 36, 26]):
        raise AssertionError(("maximizer geometry summary", maximizer_summary))
    if colored_summary != (0, 0, []):
        raise AssertionError(("colored two-point invariant", colored_summary))

    collision_data: dict[str, dict[str, object]] = {}
    for record, expected_profile in EXPECTED_COLLISION.items():
        raw_record = record.encode("ascii")
        adjacency = small.graph6_adjacency(raw_record, 8)
        networkx_decode_check(raw_record, adjacency)
        geometry = record_geometry(adjacency, pair_products, distance, dot_absolute)
        profile = (
            int(geometry["rho_ext"]),
            int(geometry["rho_weighted"]),
            int(geometry["extension"]),
        )
        if profile != expected_profile:
            raise AssertionError(("collision profile", record, profile))
        if geometry["histogram"] != EXPECTED_COLLISION_HISTOGRAM:
            raise AssertionError(("collision histogram", record, geometry["histogram"]))
        if geometry["maximizer_pair_distance"] != EXPECTED_COLLISION_PAIR_DISTANCE:
            raise AssertionError(("collision pair distance", record))
        collision_data[record] = geometry

    left, right = (collision_data[record] for record in EXPECTED_COLLISION)
    if left["colored_two_point"] == right["colored_two_point"]:
        raise AssertionError("colored two-point corruption control went undetected")

    print(f"order_9_root_records={count}")
    print("order_9_pareto={(12,0)}")
    print("order_9_M12_root_pairs={(12,0):20,(12,1):35}")
    print("order_9_M12_switching_classes={delta0:4,delta1:11}")
    print("histogram_mixed_groups=10,records=874")
    print("histogram_plus_maximizer_distance_mixed_groups=3,records=112")
    print("energy_colored_two_point_mixed_groups=0")
    print("collision=GHOgmo:(4,4,15),Gxd?Dc:(3,3,17)")
    print("collision_histogram={2:124,6:85,10:37,14:10}")
    print("collision_maximizer_pair_distance=(10,16,14,20,40)")
    print(f"geng_stream_sha256={stream_digest.hexdigest()}")
    print("deterministic_seed=413935")
    print("order_9_weighted_geometry=PASSED")


if __name__ == "__main__":
    main()

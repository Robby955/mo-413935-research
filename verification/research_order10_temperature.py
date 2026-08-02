#!/usr/bin/env python3
"""Exhaustive order-10 finite-temperature phase verification.

Root switching identifies order-10 signings with simple graphs on the nine
non-root vertices.  Nauty ``geng`` supplies one representative of every
unlabeled graph of order nine.  The verifier evaluates the complete catalogue
and proves, in exact integer polynomial arithmetic, that three absolute-energy
histograms give the lower envelope of the augmented partition function.

The only floating-point operations are the displayed numerical phase
thresholds.  Catalogue completeness trusts nauty's generator, fail-closed
against both the known graph count and the committed newline-delimited stream
digest.  A deterministic sample of records and every phase representative are
independently recomputed with the cut-energy formula.
"""

from __future__ import annotations

import hashlib
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import research_exact_small_n as small  # noqa: E402

Polynomial = tuple[int, ...]
Histogram = tuple[int, ...]

ORDER = 10
PROJECTIVE_STATES = 1 << (ORDER - 1)
ENERGY_LEVELS = tuple(range(1, ORDER * (ORDER - 1) // 2 + 1, 2))

EXPECTED_HISTOGRAM_COUNT = 6012
EXPECTED_PHASE_RECORDS = (
    ("HEhbtjK",),
    ("HCRbdO{", "HCrb`qi", "HEhuTxm", "HEhutxm"),
    ("HCRbczQ", "HCpdehU", "HCrbdxz", "HCZbeyz"),
)
EXPECTED_PHASE_SPARSE_HISTOGRAMS = (
    {1: 180, 7: 180, 9: 140, 15: 12},
    {1: 108, 3: 88, 5: 96, 7: 84, 9: 60, 11: 48, 13: 24, 15: 4},
    {3: 200, 5: 192, 11: 80, 13: 40},
)
EXPECTED_PHASE_POLYNOMIALS = (
    (512, 2816, 4512, 3680, 2120, 792, 156, 12),
    (512, 2816, 5280, 4704, 2232, 576, 76, 4),
    (512, 2816, 5792, 5600, 2520, 520, 40),
)
EXPECTED_PHASE_MAXIMA = (15, 15, 13)


def trim(polynomial: Iterable[int]) -> Polynomial:
    coefficients = list(polynomial)
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(coefficients or [0])


def polynomial_add(left: Polynomial, right: Polynomial) -> Polynomial:
    output = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        output[index] += value
    for index, value in enumerate(right):
        output[index] += value
    return trim(output)


def polynomial_subtract(left: Polynomial, right: Polynomial) -> Polynomial:
    return polynomial_add(left, tuple(-value for value in right))


def polynomial_multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    output = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            output[left_index + right_index] += left_value * right_value
    return trim(output)


def polynomial_scale(scale: int, polynomial: Polynomial) -> Polynomial:
    return trim(scale * value for value in polynomial)


def odd_cosh_polynomials() -> tuple[Polynomial, ...]:
    """Return R_j with cosh((2j+1)t) = cosh(t) R_j(4 sinh(t)^2)."""
    output: list[Polynomial] = [(1,), (1, 1)]
    for index in range(1, len(ENERGY_LEVELS) - 1):
        output.append(
            polynomial_subtract(
                polynomial_multiply((2, 1), output[index]), output[index - 1]
            )
        )
    return tuple(output)


ODD_COSH_POLYNOMIALS = odd_cosh_polynomials()


def histogram_from_sparse(sparse: dict[int, int]) -> Histogram:
    if any(energy not in ENERGY_LEVELS for energy in sparse):
        raise AssertionError(("invalid odd energy", sparse))
    return tuple(sparse.get(energy, 0) for energy in ENERGY_LEVELS)


def histogram_polynomial(histogram: Histogram) -> Polynomial:
    output: Polynomial = (0,)
    for count, basis in zip(histogram, ODD_COSH_POLYNOMIALS, strict=True):
        if count:
            output = polynomial_add(output, polynomial_scale(count, basis))
    return output


def coefficientwise_dominates(polynomial: Polynomial, candidate: Polynomial) -> bool:
    length = max(len(polynomial), len(candidate))
    return all(
        (polynomial[index] if index < len(polynomial) else 0)
        >= (candidate[index] if index < len(candidate) else 0)
        for index in range(length)
    )


def projective_pair_products() -> np.ndarray:
    masks = np.arange(PROJECTIVE_STATES, dtype=np.uint16)
    spins = np.ones((PROJECTIVE_STATES, ORDER), dtype=np.int16)
    for index in range(ORDER - 1):
        spins[:, index + 1] = 1 - 2 * ((masks >> index) & 1)
    edges = tuple(
        (row, column) for row in range(ORDER) for column in range(row + 1, ORDER)
    )
    return np.column_stack(
        [spins[:, row] * spins[:, column] for row, column in edges]
    ).astype(np.int16)


def coefficient_vector(adjacency: tuple[int, ...]) -> list[int]:
    output: list[int] = []
    for row in range(ORDER):
        for column in range(row + 1, ORDER):
            output.append(
                1 if row == 0 or not (adjacency[row - 1] >> (column - 1) & 1) else -1
            )
    return output


def exact_histogram(absolute_energies: np.ndarray) -> Histogram:
    if absolute_energies.shape != (PROJECTIVE_STATES,):
        raise AssertionError(("energy vector shape", absolute_energies.shape))
    counts = np.bincount(absolute_energies, minlength=ENERGY_LEVELS[-1] + 1)
    if int(counts.sum()) != PROJECTIVE_STATES:
        raise AssertionError("histogram lost projective states")
    if np.any(counts[0::2]):
        raise AssertionError("order-10 energy had even parity")
    return tuple(int(counts[energy]) for energy in ENERGY_LEVELS)


def independent_energy_check(
    adjacency: tuple[int, ...], signed_energies: np.ndarray
) -> None:
    expected = np.asarray(
        [
            small.energy_from_cut(adjacency, subset)
            for subset in range(PROJECTIVE_STATES)
        ],
        dtype=np.int16,
    )
    if not np.array_equal(signed_energies, expected):
        mismatch = int(np.flatnonzero(signed_energies != expected)[0])
        raise AssertionError(
            (
                "independent cut-energy recomputation",
                mismatch,
                int(signed_energies[mismatch]),
                int(expected[mismatch]),
            )
        )


def verify_conference_record(adjacency: tuple[int, ...]) -> None:
    coefficients = coefficient_vector(adjacency)
    matrix = np.zeros((ORDER, ORDER), dtype=np.int16)
    edge_index = 0
    for row in range(ORDER):
        for column in range(row + 1, ORDER):
            matrix[row, column] = matrix[column, row] = coefficients[edge_index]
            edge_index += 1
    gram = matrix @ matrix
    if not np.array_equal(gram, (ORDER - 1) * np.eye(ORDER, dtype=np.int16)):
        raise AssertionError("high-temperature phase record is not conference")


def evaluate_batch(
    batch: list[tuple[int, bytes, tuple[int, ...]]],
    pair_products: np.ndarray,
    histograms: set[Histogram],
    phase_records: dict[Histogram, list[str]],
    phase_histograms: set[Histogram],
    independently_checked: set[str],
) -> None:
    coefficient_matrix = np.asarray(
        [coefficient_vector(adjacency) for _, _, adjacency in batch], dtype=np.int16
    )
    signed_energy_matrix = coefficient_matrix @ pair_products.T
    if signed_energy_matrix.dtype != np.int16:
        raise AssertionError(
            ("unexpected exact matrix dtype", signed_energy_matrix.dtype)
        )

    for (record_number, raw_record, adjacency), signed_energies in zip(
        batch, signed_energy_matrix, strict=True
    ):
        absolute_energies = np.abs(signed_energies)
        histogram = exact_histogram(absolute_energies)
        histograms.add(histogram)
        record = raw_record.decode("ascii")
        if histogram in phase_histograms:
            phase_records[histogram].append(record)

        # Deterministic sparse checks catch coefficient ordering, graph6
        # decoding, and matrix-orientation errors independently of NumPy.
        if (
            record_number <= 3
            or record_number % 9973 == 0
            or any(record in records for records in EXPECTED_PHASE_RECORDS)
        ):
            independent_energy_check(adjacency, signed_energies)
            independently_checked.add(record)
        if record == EXPECTED_PHASE_RECORDS[0][0]:
            verify_conference_record(adjacency)


def exact_catalogue() -> tuple[set[Histogram], dict[Histogram, list[str]], str, int]:
    geng = small.locate_geng(None)
    pair_products = projective_pair_products()
    phase_histograms = tuple(
        histogram_from_sparse(sparse) for sparse in EXPECTED_PHASE_SPARSE_HISTOGRAMS
    )
    phase_histogram_set = set(phase_histograms)
    histogram_set: set[Histogram] = set()
    phase_records: dict[Histogram, list[str]] = defaultdict(list)
    independently_checked: set[str] = set()
    stream_digest = hashlib.sha256()
    batch: list[tuple[int, bytes, tuple[int, ...]]] = []
    count = 0

    for raw_record in small.geng_records(geng, ORDER - 1):
        stream_digest.update(raw_record + b"\n")
        adjacency = small.graph6_adjacency(raw_record, ORDER - 1)
        count += 1
        batch.append((count, raw_record, adjacency))
        if len(batch) == 4096:
            evaluate_batch(
                batch,
                pair_products,
                histogram_set,
                phase_records,
                phase_histogram_set,
                independently_checked,
            )
            batch.clear()
    if batch:
        evaluate_batch(
            batch,
            pair_products,
            histogram_set,
            phase_records,
            phase_histogram_set,
            independently_checked,
        )

    expected_count = small.UNLABELED_GRAPH_COUNTS[ORDER - 1]
    if count != expected_count:
        raise AssertionError(("incomplete rooted order-10 catalogue", count))
    digest = stream_digest.hexdigest()
    expected_digest = small.EXPECTED_GENG_STREAM_SHA256[ORDER]
    if digest != expected_digest:
        raise AssertionError(("order-10 geng stream digest", digest, expected_digest))
    if len(histogram_set) != EXPECTED_HISTOGRAM_COUNT:
        raise AssertionError(("absolute-energy histogram count", len(histogram_set)))
    expected_independent = {
        record for records in EXPECTED_PHASE_RECORDS for record in records
    }
    if not expected_independent.issubset(independently_checked):
        raise AssertionError(
            ("phase representative missed independent check", independently_checked)
        )
    return histogram_set, phase_records, digest, count


def prove_phase_envelope(
    histograms: set[Histogram], phase_records: dict[Histogram, list[str]]
) -> tuple[Polynomial, Polynomial, Polynomial]:
    candidate_histograms = tuple(
        histogram_from_sparse(sparse) for sparse in EXPECTED_PHASE_SPARSE_HISTOGRAMS
    )
    candidate_polynomials = tuple(
        histogram_polynomial(histogram) for histogram in candidate_histograms
    )
    if candidate_polynomials != EXPECTED_PHASE_POLYNOMIALS:
        raise AssertionError(("candidate polynomial", candidate_polynomials))
    candidate_maxima = tuple(
        max(energy for energy, count in zip(ENERGY_LEVELS, histogram) if count)
        for histogram in candidate_histograms
    )
    if candidate_maxima != EXPECTED_PHASE_MAXIMA:
        raise AssertionError(("candidate maxima", candidate_maxima))

    for histogram, expected_records in zip(
        candidate_histograms, EXPECTED_PHASE_RECORDS, strict=True
    ):
        records = tuple(phase_records[histogram])
        if records != expected_records:
            raise AssertionError(("candidate record list", records, expected_records))

    for histogram in histograms:
        polynomial = histogram_polynomial(histogram)
        if not any(
            coefficientwise_dominates(polynomial, candidate)
            for candidate in candidate_polynomials
        ):
            raise AssertionError(
                ("undominated order-10 histogram", histogram, polynomial)
            )

    first, second, third = candidate_polynomials
    first_factor = polynomial_scale(
        8,
        polynomial_multiply(
            (0, 0, 1),
            polynomial_multiply(
                (-2, 1),
                polynomial_multiply(
                    (1, 1),
                    polynomial_multiply((3, 1), polynomial_multiply((4, 1), (4, 1))),
                ),
            ),
        ),
    )
    second_factor = polynomial_scale(
        4,
        polynomial_multiply(
            (0, 0, 1),
            polynomial_multiply(polynomial_multiply((4, 1), (4, 1)), (-8, -10, 1, 1)),
        ),
    )
    if polynomial_subtract(first, second) != first_factor:
        raise AssertionError("P0-P1 factorization")
    if polynomial_subtract(second, third) != second_factor:
        raise AssertionError("P1-P2 factorization")

    # The cubic z^3+z^2-10z-8 has exactly one positive root by Descartes'
    # rule, and the exact signs at 2 and 4 locate it.  Together with the two
    # factorizations this proves the three stated open phase intervals.
    def cubic(z: int) -> int:
        return z**3 + z**2 - 10 * z - 8

    if cubic(2) >= 0 or cubic(4) <= 0:
        raise AssertionError("positive cubic root bracket")
    coefficient_signs = (-1, -1, 1, 1)
    sign_changes = sum(
        left != right for left, right in zip(coefficient_signs, coefficient_signs[1:])
    )
    if sign_changes != 1:
        raise AssertionError("Descartes sign-change certificate")
    return first, second, third


def positive_cubic_root() -> float:
    lower, upper = 2.0, 4.0
    for _ in range(100):
        middle = (lower + upper) / 2.0
        value = middle**3 + middle**2 - 10.0 * middle - 8.0
        if value < 0.0:
            lower = middle
        else:
            upper = middle
    root = (lower + upper) / 2.0
    if not 3.0838723594 < root < 3.0838723595:
        raise AssertionError(("cubic threshold", root))
    return root


def corruption_controls(candidates: tuple[Polynomial, ...]) -> None:
    corrupted_histogram = dict(EXPECTED_PHASE_SPARSE_HISTOGRAMS[0])
    corrupted_histogram[1] -= 1
    if (
        histogram_polynomial(histogram_from_sparse(corrupted_histogram))
        == candidates[0]
    ):
        raise AssertionError("histogram corruption was not detected")

    wrong_factor = polynomial_scale(
        8,
        polynomial_multiply(
            (0, 0, 1),
            polynomial_multiply(
                (-1, 1),
                polynomial_multiply(
                    (1, 1),
                    polynomial_multiply((3, 1), polynomial_multiply((4, 1), (4, 1))),
                ),
            ),
        ),
    )
    if polynomial_subtract(candidates[0], candidates[1]) == wrong_factor:
        raise AssertionError("phase-factor corruption was not detected")

    impossible_candidates = tuple(
        polynomial_add(candidate, (1,)) for candidate in candidates
    )
    if any(
        coefficientwise_dominates(candidates[0], candidate)
        for candidate in impossible_candidates
    ):
        raise AssertionError("coefficientwise-direction corruption was not detected")


def main() -> None:
    histograms, phase_records, digest, count = exact_catalogue()
    candidates = prove_phase_envelope(histograms, phase_records)
    corruption_controls(candidates)

    zeta = positive_cubic_root()
    first_threshold = math.asinh(1.0 / math.sqrt(2.0))
    second_threshold = math.asinh(math.sqrt(zeta) / 2.0)
    first_extensive = math.sqrt(ORDER) * first_threshold
    second_extensive = math.sqrt(ORDER) * second_threshold

    print(f"order_10_root_records={count}")
    print(f"order_10_absolute_energy_histograms={len(histograms)}")
    print(
        "phase_records="
        + ",".join(str(len(records)) for records in EXPECTED_PHASE_RECORDS)
    )
    print("phase_maxima=" + ",".join(map(str, EXPECTED_PHASE_MAXIMA)))
    print("phase_transition_z=2,root(z^3+z^2-10z-8)")
    print(f"positive_cubic_root={zeta:.12f}")
    print(f"temperature_thresholds={first_threshold:.12f},{second_threshold:.12f}")
    print(f"extensive_beta_thresholds={first_extensive:.12f},{second_extensive:.12f}")
    print(f"geng_stream_sha256={digest}")
    print("independent_recomputation=PASSED")
    print("corruption_controls=histogram,coefficientwise_direction,phase_factorization")
    print("order_10_temperature_phases=PASSED")


if __name__ == "__main__":
    main()

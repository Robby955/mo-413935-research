#!/usr/bin/env python3
"""Exact finite calibration of the microcanonical composition theorem.

The script evaluates every balanced split of banked optimal signings at
orders 12, 13, and 14.  For each split it computes the three complete local
deficit histograms, their product-convolution order statistic Lambda, and the
true relative-gauge gain.  The latter is exact because the identity gauge is
globally optimal and every relative gauge is another signing of the same
order.

All pass/fail decisions use integer arithmetic.  The displayed
near-subadditive targets use floating point only as descriptive output.
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter
from dataclasses import dataclass
from itertools import combinations

import check_conference_examples as conference
import research_order13_certify as order13


Matrix = list[list[int]]
Histogram = Counter[int]
ProfileRecord = tuple[int, int, int, int, int, int, int, int]


@dataclass(frozen=True)
class ExpectedCase:
    label: str
    order: int
    left_order: int
    optimum: int
    left_optimum: int
    right_optimum: int
    partitions: int
    profile_types: int
    digest: str
    lambda_range: tuple[int, int]
    scalar_bound_range: tuple[int, int]
    geometry_gap_range: tuple[int, int]
    best_bound_count: int
    child_optimal_count: int
    child_optimal_best_bound: int
    target_lattice_energy: int
    target_count_range: tuple[int, int]


EXPECTED = {
    order13.SURVIVOR_RECORDS[0]: ExpectedCase(
        order13.SURVIVOR_RECORDS[0],
        12,
        6,
        18,
        5,
        5,
        462,
        12,
        "2c29479f48ff9644bfb860640ba35450af837502371b17f59153653b708fbc01",
        (0, 10),
        (24, 26),
        (6, 8),
        7,
        7,
        24,
        16,
        (213312, 269122),
    ),
    order13.SURVIVOR_RECORDS[1]: ExpectedCase(
        order13.SURVIVOR_RECORDS[1],
        12,
        6,
        18,
        5,
        5,
        462,
        12,
        "5ff5291913cd0cf4db0e3cb2648934e99607fc59f184ff734abb9a765c31d444",
        (0, 10),
        (24, 26),
        (6, 8),
        7,
        7,
        24,
        16,
        (213312, 269122),
    ),
    "C14-minus-infinity": ExpectedCase(
        "C14-minus-infinity",
        13,
        6,
        20,
        5,
        9,
        1716,
        9,
        "b9ce312b08962293639003979cac33a0c6a852dd40209a00bd23653eb437d8aa",
        (2, 10),
        (28, 30),
        (8, 10),
        364,
        52,
        28,
        20,
        (436864, 495632),
    ),
    "C14": ExpectedCase(
        "C14",
        14,
        7,
        21,
        9,
        9,
        1716,
        3,
        "e96d8070d53fcdd2a47ce95eddb865d666c74dd02557df53ae8ec993d3fbd4b7",
        (6, 10),
        (31, 33),
        (10, 12),
        364,
        624,
        31,
        27,
        (231581, 305465),
    ),
}

EXPECTED_GAIN_PAIRS = {
    order13.SURVIVOR_RECORDS[0]: (
        (0, 6, 6),
        (2, 8, 1),
        (2, 10, 70),
        (4, 12, 135),
        (8, 16, 120),
        (10, 18, 130),
    ),
    order13.SURVIVOR_RECORDS[1]: (
        (0, 6, 6),
        (2, 8, 1),
        (2, 10, 70),
        (4, 12, 135),
        (8, 16, 120),
        (10, 18, 130),
    ),
    "C14-minus-infinity": (
        (2, 10, 52),
        (6, 14, 156),
        (6, 16, 338),
        (8, 16, 156),
        (8, 18, 546),
        (10, 20, 468),
    ),
    "C14": (
        (6, 16, 364),
        (6, 18, 260),
        (10, 22, 1092),
    ),
}


# (independent ceiling, profile gain, number of product triples above the
# exact lattice threshold, number of balanced splits).  The threshold is the
# smallest energy of the correct parity strictly above
# (F(left)^(2/3) + F(right)^(2/3))^(3/2).
EXPECTED_TARGET_COUNTS = {
    order13.SURVIVOR_RECORDS[0]: (
        (24, 0, 253376, 6),
        (26, 2, 213312, 1),
        (28, 2, 254464, 10),
        (28, 2, 257792, 30),
        (28, 2, 261184, 30),
        (30, 4, 262400, 60),
        (30, 4, 265536, 15),
        (30, 4, 266752, 60),
        (34, 8, 261208, 60),
        (34, 8, 264982, 60),
        (36, 10, 261538, 60),
        (36, 10, 265752, 30),
        (36, 10, 269122, 40),
    ),
    order13.SURVIVOR_RECORDS[1]: (
        (24, 0, 253376, 6),
        (26, 2, 213312, 1),
        (28, 2, 254464, 10),
        (28, 2, 257792, 30),
        (28, 2, 261184, 30),
        (30, 4, 262400, 60),
        (30, 4, 265536, 15),
        (30, 4, 266752, 60),
        (34, 8, 261208, 60),
        (34, 8, 264982, 60),
        (36, 10, 261538, 60),
        (36, 10, 265752, 30),
        (36, 10, 269122, 40),
    ),
    "C14-minus-infinity": (
        (30, 2, 436864, 52),
        (34, 6, 454144, 156),
        (36, 6, 460370, 26),
        (36, 6, 465164, 156),
        (36, 6, 476592, 156),
        (36, 8, 453052, 156),
        (38, 8, 474548, 156),
        (38, 8, 484960, 78),
        (38, 8, 485036, 156),
        (38, 8, 495632, 156),
        (40, 10, 483552, 78),
        (40, 10, 484084, 78),
        (40, 10, 485426, 156),
        (40, 10, 494272, 156),
    ),
    "C14": (
        (37, 6, 231581, 364),
        (39, 6, 231807, 182),
        (39, 6, 303122, 78),
        (43, 10, 304908, 546),
        (43, 10, 305465, 546),
    ),
}


SPINS = {order: order13.projective_spins(order) for order in range(1, 15)}


def principal(matrix: Matrix, vertices: tuple[int, ...]) -> Matrix:
    return [[matrix[row][column] for column in vertices] for row in vertices]


def graph_profile(matrix: Matrix) -> tuple[int, Histogram]:
    energy_counts = Counter(
        order13.energy(matrix, spin) for spin in SPINS[len(matrix)]
    )
    maximum = max(abs(value) for value in energy_counts)
    deficits: Histogram = Counter()
    for value, count in energy_counts.items():
        deficits[maximum - value] += count
        deficits[maximum + value] += count
    assert sum(deficits.values()) == 2 ** len(matrix)
    return maximum, deficits


def rectangular_profile(
    matrix: Matrix,
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> tuple[int, Histogram]:
    absolute_energies: Histogram = Counter()
    for right_spin in SPINS[len(right)]:
        fields = [
            sum(
                matrix[row][column] * right_spin[index]
                for index, column in enumerate(right)
            )
            for row in left
        ]
        for left_spin in SPINS[len(left)]:
            value = abs(
                sum(sign * field for sign, field in zip(left_spin, fields))
            )
            absolute_energies[value] += 1
    maximum = max(absolute_energies)
    deficits = Counter(
        {maximum - value: count for value, count in absolute_energies.items()}
    )
    assert sum(deficits.values()) == 2 ** (len(left) + len(right) - 2)
    return maximum, deficits


def convolve(left: Histogram, right: Histogram) -> Histogram:
    result: Histogram = Counter()
    for left_value, left_count in left.items():
        for right_value, right_count in right.items():
            result[left_value + right_value] += left_count * right_count
    return result


def order_statistic(histogram: Histogram, rank: int) -> int:
    cumulative = 0
    for value in sorted(histogram):
        cumulative += histogram[value]
        if cumulative >= rank:
            return value
    raise AssertionError(("rank exceeds histogram", rank, cumulative))


def profile_record(
    matrix: Matrix,
    left: tuple[int, ...],
    right: tuple[int, ...],
    optimum: int,
    target_lattice_energy: int,
    graph_cache: dict[tuple[int, ...], tuple[int, Histogram]],
) -> tuple[ProfileRecord, int]:
    for vertices in (left, right):
        if vertices not in graph_cache:
            graph_cache[vertices] = graph_profile(principal(matrix, vertices))
    left_maximum, left_profile = graph_cache[left]
    right_maximum, right_profile = graph_cache[right]
    cross_maximum, cross_profile = rectangular_profile(matrix, left, right)
    independent_ceiling = left_maximum + right_maximum + cross_maximum
    product_profile = convolve(convolve(left_profile, right_profile), cross_profile)
    order = len(matrix)
    assert sum(product_profile.values()) == 2 ** (2 * order - 2)
    profile_gain = order_statistic(product_profile, 2 ** (order - 1))
    scalar_bound = independent_ceiling - profile_gain
    true_gain = independent_ceiling - optimum
    target_count = sum(
        count
        for deficit, count in product_profile.items()
        if independent_ceiling - deficit >= target_lattice_energy
    )
    return (
        (
        left_maximum,
        right_maximum,
        cross_maximum,
        independent_ceiling,
        profile_gain,
        scalar_bound,
        true_gain,
        true_gain - profile_gain,
        ),
        target_count,
    )


def target_is_strictly_between(
    left_optimum: int,
    right_optimum: int,
    lower: int,
    upper: int,
) -> bool:
    """Compare (a^(2/3)+b^(2/3))^(3/2) with integers exactly.

    If y is the square of the target, then
    (y-a^2-b^2)^3 = 27a^2b^2 y.  The corresponding cubic is strictly
    increasing for y>a^2+b^2+3ab, which contains all comparison points used
    here.  No floating-point arithmetic enters this pass/fail decision.
    """

    a = left_optimum
    b = right_optimum
    baseline = a * a + b * b

    def comparison_polynomial(z: int) -> int:
        y = z * z
        assert y > baseline + 3 * a * b
        return (y - baseline) ** 3 - 27 * a * a * b * b * y

    return comparison_polynomial(lower) < 0 < comparison_polynomial(upper)


def scan_case(matrix: Matrix, expected: ExpectedCase) -> Counter[ProfileRecord]:
    full_maximum, _ = graph_profile(matrix)
    assert len(matrix) == expected.order
    assert full_maximum == expected.optimum

    vertices = tuple(range(expected.order))
    graph_cache: dict[tuple[int, ...], tuple[int, Histogram]] = {}
    records: Counter[ProfileRecord] = Counter()
    target_counts: Counter[tuple[int, int, int]] = Counter()
    digest = hashlib.sha256()
    for left in combinations(vertices, expected.left_order):
        # Balanced complements describe the same unordered split.  Requiring
        # vertex zero selects exactly one of the two orientations.
        if 2 * expected.left_order == expected.order and 0 not in left:
            continue
        right = tuple(vertex for vertex in vertices if vertex not in left)
        record, target_count = profile_record(
            matrix,
            left,
            right,
            expected.optimum,
            expected.target_lattice_energy,
            graph_cache,
        )
        records[record] += 1
        target_counts[(record[3], record[4], target_count)] += 1
        digest.update(
            (
                ",".join(map(str, left))
                + "|"
                + ",".join(map(str, record))
                + "\n"
            ).encode("ascii")
        )

    assert sum(records.values()) == expected.partitions
    assert len(records) == expected.profile_types
    assert digest.hexdigest() == expected.digest

    lambdas = [record[4] for record in records]
    scalar_bounds = [record[5] for record in records]
    geometry_gaps = [record[7] for record in records]
    assert (min(lambdas), max(lambdas)) == expected.lambda_range
    assert (min(scalar_bounds), max(scalar_bounds)) == expected.scalar_bound_range
    assert (min(geometry_gaps), max(geometry_gaps)) == expected.geometry_gap_range
    assert min(geometry_gaps) > 0

    gain_pairs = Counter()
    for record, count in records.items():
        gain_pairs[(record[4], record[6])] += count
    expected_gain_pairs = Counter(
        {
            (profile_gain, true_gain): count
            for profile_gain, true_gain, count in EXPECTED_GAIN_PAIRS[expected.label]
        }
    )
    assert gain_pairs == expected_gain_pairs

    expected_target_counts = Counter(
        {
            (ceiling, profile_gain, target_count): count
            for ceiling, profile_gain, target_count, count in (
                EXPECTED_TARGET_COUNTS[expected.label]
            )
        }
    )
    assert target_counts == expected_target_counts
    observed_target_counts = [key[2] for key in target_counts]
    assert (
        min(observed_target_counts),
        max(observed_target_counts),
    ) == expected.target_count_range
    fiber_count = 2 ** (expected.order - 1)
    assert min(observed_target_counts) >= fiber_count

    best_bound = min(scalar_bounds)
    best_bound_count = sum(
        count for record, count in records.items() if record[5] == best_bound
    )
    assert best_bound_count == expected.best_bound_count

    child_optimal_count = sum(
        count
        for record, count in records.items()
        if record[0] == expected.left_optimum
        and record[1] == expected.right_optimum
    )
    child_optimal_best_bound = min(
        record[5]
        for record in records
        if record[0] == expected.left_optimum
        and record[1] == expected.right_optimum
    )
    assert child_optimal_count == expected.child_optimal_count
    assert child_optimal_best_bound == expected.child_optimal_best_bound

    assert target_is_strictly_between(
        expected.left_optimum,
        expected.right_optimum,
        expected.target_lattice_energy - 2,
        expected.target_lattice_energy,
    )

    near_subadditive_target = (
        expected.left_optimum ** (2.0 / 3.0)
        + expected.right_optimum ** (2.0 / 3.0)
    ) ** 1.5
    print(
        f"case={expected.label} order={expected.order} "
        f"partitions={expected.partitions} profile_types={expected.profile_types} "
        f"lambda={expected.lambda_range[0]}..{expected.lambda_range[1]} "
        f"scalar_bound={expected.scalar_bound_range[0]}..{expected.scalar_bound_range[1]} "
        f"geometry_gap={expected.geometry_gap_range[0]}..{expected.geometry_gap_range[1]} "
        f"child_optimal_splits={child_optimal_count} "
        f"best_bound={best_bound} target={near_subadditive_target:.12f} "
        f"excess={best_bound - near_subadditive_target:.12f} "
        f"target_lattice_energy={expected.target_lattice_energy} "
        f"target_triples={expected.target_count_range[0]}..{expected.target_count_range[1]} "
        f"fiber_count={fiber_count} "
        f"sha256={expected.digest}"
    )
    print(
        "gain_pairs="
        + ",".join(
            f"{profile_gain}/{true_gain}:{count}"
            for (profile_gain, true_gain), count in sorted(gain_pairs.items())
        )
    )
    return records


def verify_corruption_controls(matrix: Matrix) -> None:
    # Omitting the balanced orientation guard doubles the split count.
    all_balanced = math.comb(14, 7)
    guarded_balanced = sum(1 for subset in combinations(range(14), 7) if 0 in subset)
    assert all_balanced == 2 * guarded_balanced == 3432

    left = tuple(range(7))
    right = tuple(range(7, 14))
    correct_maximum, correct_profile = rectangular_profile(matrix, left, right)
    signed_values = []
    for right_spin in SPINS[7]:
        fields = [
            sum(
                matrix[row][column] * right_spin[index]
                for index, column in enumerate(right)
            )
            for row in left
        ]
        signed_values.extend(
            sum(sign * field for sign, field in zip(left_spin, fields))
            for left_spin in SPINS[7]
        )
    corrupted_maximum = max(signed_values)
    corrupted_profile = Counter(
        corrupted_maximum - value for value in signed_values
    )
    assert (correct_maximum, correct_profile) != (
        corrupted_maximum,
        corrupted_profile,
    )
    print("corruption_controls=balanced_double_count,cross_absolute_value")


def main() -> None:
    survivor_records = []
    for record in order13.SURVIVOR_RECORDS:
        survivor_records.append(scan_case(order13.graph6_signing(record), EXPECTED[record]))
    # The two rooted survivor records represent the same unrooted optimum;
    # their aggregate profile distributions must agree even though their
    # subset-labelled stream digests differ.
    assert survivor_records[0] == survivor_records[1]

    c14 = conference.paley_conference(13)
    conference.check_conference_identity(c14)
    c13 = conference.principal_submatrix(c14, 0)
    scan_case(c13, EXPECTED["C14-minus-infinity"])
    scan_case(c14, EXPECTED["C14"])
    verify_corruption_controls(c14)
    print("relative_profile_calibration=PASSED")


if __name__ == "__main__":
    main()

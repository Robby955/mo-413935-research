#!/usr/bin/env python3
"""Exact checks for the framed common-Hadamard order-16 lift.

This verifier separates two statements:

* the displayed order-16 sign matrix has Boolean quadratic maximum 30, so
  F(16) <= 30;
* among lifts with the pinned oriented Hadamard cross frame and arbitrary
  order-four signings on the diagonal, the minimum maximum is exactly 30.

The restricted lower bound is a six-state certificate, not an exhaustive
search over the 64^4 choices of diagonal blocks.
"""

from __future__ import annotations

import hashlib
import itertools
from collections import Counter


if not __debug__:
    raise RuntimeError("assertions must be enabled")


Matrix = tuple[tuple[int, ...], ...]
Spins = tuple[int, ...]

VERTEX_PAIRS = tuple(itertools.combinations(range(4), 2))
CLOUD_PAIRS = VERTEX_PAIRS

ORIENTED_H: Matrix = (
    (1, 1, 1, 1),
    (1, 1, -1, -1),
    (1, -1, 1, -1),
    (-1, 1, 1, -1),
)
BASE_SIGNS = (1, 1, 1, -1, 1, -1)

# Edge order: 01, 02, 03, 12, 13, 23.
P_EDGES = (1, 1, 1, -1, 1, -1)
R_EDGES = (-1, -1, 1, 1, 1, -1)
WITNESS_INTERNAL_EDGES = (P_EDGES, R_EDGES, P_EDGES, R_EDGES)

# Each entry is a (+28 state, -28 state).  The other three clouds have the
# same projective spins.  On cloud 1, the three pairs give the directed cycle
# 0 -> 3 -> 2 -> 0 in PROJECTIVE_REPS.
OBSTRUCTION_PAIRS: tuple[tuple[tuple[Spins, ...], tuple[Spins, ...]], ...] = (
    (
        (
            (1, 1, 1, 1),
            (1, 1, 1, 1),
            (-1, -1, -1, -1),
            (1, 1, 1, -1),
        ),
        (
            (1, 1, 1, 1),
            (-1, -1, 1, 1),
            (-1, -1, -1, -1),
            (-1, -1, -1, 1),
        ),
    ),
    (
        (
            (1, 1, -1, -1),
            (1, 1, -1, -1),
            (-1, -1, 1, 1),
            (1, 1, -1, 1),
        ),
        (
            (1, 1, -1, -1),
            (-1, -1, 1, -1),
            (-1, -1, 1, 1),
            (-1, -1, 1, -1),
        ),
    ),
    (
        (
            (1, 1, -1, -1),
            (1, 1, -1, 1),
            (-1, -1, 1, 1),
            (1, 1, -1, 1),
        ),
        (
            (1, 1, -1, -1),
            (-1, -1, -1, -1),
            (-1, -1, 1, 1),
            (-1, -1, 1, -1),
        ),
    ),
)

EXPECTED_FULL_HISTOGRAM = (
    (-30, 38),
    (-28, 88),
    (-26, 186),
    (-24, 456),
    (-22, 670),
    (-20, 1056),
    (-18, 1634),
    (-16, 1904),
    (-14, 2518),
    (-12, 2960),
    (-10, 3050),
    (-8, 3640),
    (-6, 3774),
    (-4, 4088),
    (-2, 4514),
    (0, 4384),
    (2, 4514),
    (4, 4088),
    (6, 3774),
    (8, 3640),
    (10, 3050),
    (12, 2960),
    (14, 2518),
    (16, 1904),
    (18, 1634),
    (20, 1056),
    (22, 670),
    (24, 456),
    (26, 186),
    (28, 88),
    (30, 38),
)
EXPECTED_HISTOGRAM_SHA256 = (
    "875e4f931630501ec7730abad70df7ab029602f7a48b4b10d56df9de1e319388"
)
EXPECTED_MATRIX_SHA256 = (
    "352392c57458568ddbf2920d4cb487f67d21fb491ab8c75c6862e9c7fc6a9181"
)

PROJECTIVE_REPS = tuple(
    (1,) + tail for tail in itertools.product((1, -1), repeat=3)
)


def edge_matrix(edges: tuple[int, ...]) -> Matrix:
    matrix = [[0] * 4 for _ in range(4)]
    for value, (left, right) in zip(edges, VERTEX_PAIRS):
        matrix[left][right] = value
        matrix[right][left] = value
    return tuple(tuple(row) for row in matrix)


def transpose(matrix: Matrix) -> Matrix:
    return tuple(
        tuple(matrix[row][column] for row in range(len(matrix)))
        for column in range(len(matrix[0]))
    )


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum(
                left[row][inner] * right[inner][column]
                for inner in range(len(right))
            )
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def is_hadamard(matrix: Matrix) -> bool:
    product = matmul(matrix, transpose(matrix))
    return product == tuple(
        tuple(4 if row == column else 0 for column in range(4))
        for row in range(4)
    )


def projective_spins(order: int):
    for mask in range(1 << (order - 1)):
        yield (1,) + tuple(
            -1 if (mask >> (coordinate - 1)) & 1 else 1
            for coordinate in range(1, order)
        )


def projective_key(spins: Spins) -> Spins:
    sign = spins[0]
    return tuple(sign * value for value in spins)


def quadratic_energy(matrix: Matrix, spins: Spins) -> int:
    return sum(
        matrix[left][right] * spins[left] * spins[right]
        for left in range(len(matrix))
        for right in range(left + 1, len(matrix))
    )


def internal_energy(edges: tuple[int, ...], spins: Spins) -> int:
    return sum(
        value * spins[left] * spins[right]
        for value, (left, right) in zip(edges, VERTEX_PAIRS)
    )


def cross_energy(clouds: tuple[Spins, ...]) -> int:
    return sum(
        BASE_SIGNS[pair_index]
        * sum(
            clouds[left_cloud][row]
            * ORIENTED_H[row][column]
            * clouds[right_cloud][column]
            for row in range(4)
            for column in range(4)
        )
        for pair_index, (left_cloud, right_cloud) in enumerate(CLOUD_PAIRS)
    )


def block_energy(clouds: tuple[Spins, ...]) -> int:
    return cross_energy(clouds) + sum(
        internal_energy(edges, clouds[cloud])
        for cloud, edges in enumerate(WITNESS_INTERNAL_EDGES)
    )


def assemble_matrix(internal_edges: tuple[tuple[int, ...], ...]) -> Matrix:
    matrix = [[0] * 16 for _ in range(16)]
    for cloud, edges in enumerate(internal_edges):
        block = edge_matrix(edges)
        for row in range(4):
            for column in range(4):
                matrix[4 * cloud + row][4 * cloud + column] = block[row][column]
    for pair_index, (left_cloud, right_cloud) in enumerate(CLOUD_PAIRS):
        sign = BASE_SIGNS[pair_index]
        for row in range(4):
            for column in range(4):
                value = sign * ORIENTED_H[row][column]
                matrix[4 * left_cloud + row][4 * right_cloud + column] = value
                matrix[4 * right_cloud + column][4 * left_cloud + row] = value
    return tuple(tuple(row) for row in matrix)


def verify_structure(matrix: Matrix) -> None:
    if not is_hadamard(ORIENTED_H):
        raise AssertionError("the pinned cross frame is not Hadamard")
    if ORIENTED_H == transpose(ORIENTED_H):
        raise AssertionError("the pinned cross frame unexpectedly became symmetric")
    if any(len(row) != 16 for row in matrix):
        raise AssertionError("assembled matrix is not 16 by 16")
    for row in range(16):
        for column in range(16):
            if matrix[row][column] != matrix[column][row]:
                raise AssertionError("assembled matrix is not symmetric")
            expected_values = (0,) if row == column else (-1, 1)
            if matrix[row][column] not in expected_values:
                raise AssertionError("assembled matrix is not an admissible signing")


def verify_order_four_blocks() -> None:
    optimum = min(
        max(
            abs(internal_energy(edges, spins))
            for spins in projective_spins(4)
        )
        for edges in itertools.product((1, -1), repeat=6)
    )
    if optimum != 4:
        raise AssertionError(("recomputed F(4)", optimum))
    for label, edges in (("P", P_EDGES), ("R", R_EDGES)):
        maximum = max(
            abs(internal_energy(edges, spins))
            for spins in projective_spins(4)
        )
        if maximum != optimum:
            raise AssertionError((label, maximum, optimum))

    # R_ij = s_i s_j P_{perm(i),perm(j)}.
    permutation = (1, 0, 3, 2)
    switching = (1, -1, -1, -1)
    p_matrix = edge_matrix(P_EDGES)
    r_matrix = edge_matrix(R_EDGES)
    transformed = tuple(
        tuple(
            switching[row]
            * switching[column]
            * p_matrix[permutation[row]][permutation[column]]
            for column in range(4)
        )
        for row in range(4)
    )
    if transformed != r_matrix:
        raise AssertionError("the P-to-R switching-permutation certificate failed")


def verify_six_state_obstruction() -> None:
    directed_edges: list[tuple[int, int]] = []
    for positive, negative in OBSTRUCTION_PAIRS:
        if cross_energy(positive) != 28 or cross_energy(negative) != -28:
            raise AssertionError("an obstruction state has the wrong cross energy")
        positive_keys = tuple(projective_key(cloud) for cloud in positive)
        negative_keys = tuple(projective_key(cloud) for cloud in negative)
        if any(
            positive_keys[cloud] != negative_keys[cloud]
            for cloud in (0, 2, 3)
        ):
            raise AssertionError("the paired states do not cancel other blocks")
        directed_edges.append(
            (
                PROJECTIVE_REPS.index(positive_keys[1]),
                PROJECTIVE_REPS.index(negative_keys[1]),
            )
        )
    if tuple(directed_edges) != ((0, 3), (3, 2), (2, 0)):
        raise AssertionError(("obstruction cycle", directed_edges))

    # If a full lift had maximum at most 28, each pair would imply
    # q_D(source) <= q_D(target).  The directed cycle would make the three
    # response values equal.  But q_D(rep_0)-q_D(rep_2) is twice a sum of
    # three signs, hence cannot be zero.
    for edges in itertools.product((1, -1), repeat=6):
        q0 = internal_energy(edges, PROJECTIVE_REPS[0])
        q2 = internal_energy(edges, PROJECTIVE_REPS[2])
        q3 = internal_energy(edges, PROJECTIVE_REPS[3])
        odd_sum_formula = 2 * (edges[1] + edges[3] + edges[5])
        if q0 - q2 != odd_sum_formula:
            raise AssertionError("the odd-sum response identity failed")
        if odd_sum_formula == 0:
            raise AssertionError("a sum of three signs vanished")
        if q0 == q3 == q2:
            raise AssertionError("the impossible response equality occurred")


def common_internal_profile() -> Counter[int]:
    states: list[tuple[tuple[int, ...], int]] = []
    for spins in projective_spins(16):
        clouds = tuple(spins[4 * cloud : 4 * cloud + 4] for cloud in range(4))
        response_indices = tuple(
            PROJECTIVE_REPS.index(projective_key(cloud)) for cloud in clouds
        )
        states.append((response_indices, cross_energy(clouds)))

    profile: Counter[int] = Counter()
    for edges in itertools.product((1, -1), repeat=6):
        responses = tuple(
            internal_energy(edges, representative)
            for representative in PROJECTIVE_REPS
        )
        maximum = max(
            abs(cross + sum(responses[index] for index in indices))
            for indices, cross in states
        )
        profile[maximum] += 1
    return profile


def verify_witness(matrix: Matrix) -> tuple[str, int]:
    matrix_digest = hashlib.sha256(repr(matrix).encode()).hexdigest()
    if matrix_digest != EXPECTED_MATRIX_SHA256:
        raise AssertionError(("matrix digest", matrix_digest))
    projective_histogram: Counter[int] = Counter()
    cross_maximum = 0
    for spins in projective_spins(16):
        clouds = tuple(spins[4 * cloud : 4 * cloud + 4] for cloud in range(4))
        direct = quadratic_energy(matrix, spins)
        blocked = block_energy(clouds)
        if direct != blocked:
            raise AssertionError("direct and blockwise energy computations disagree")
        projective_histogram[direct] += 1
        cross_maximum = max(cross_maximum, abs(cross_energy(clouds)))

    full_histogram = tuple(
        (energy, 2 * count) for energy, count in sorted(projective_histogram.items())
    )
    digest = hashlib.sha256(repr(list(full_histogram)).encode()).hexdigest()
    if full_histogram != EXPECTED_FULL_HISTOGRAM:
        raise AssertionError("the complete energy histogram changed")
    if digest != EXPECTED_HISTOGRAM_SHA256:
        raise AssertionError(("energy histogram digest", digest))
    maximum = max(abs(energy) for energy in projective_histogram)
    projective_maximizers = sum(
        count
        for energy, count in projective_histogram.items()
        if abs(energy) == maximum
    )
    if maximum != 30 or projective_maximizers != 38 or cross_maximum != 28:
        raise AssertionError(
            ("witness maxima", maximum, projective_maximizers, cross_maximum)
        )
    return digest, projective_maximizers


def verify_corruption_controls(matrix: Matrix) -> None:
    corrupted_h = [list(row) for row in ORIENTED_H]
    corrupted_h[0][0] *= -1
    if is_hadamard(tuple(tuple(row) for row in corrupted_h)):
        raise AssertionError("Hadamard corruption was not detected")

    corrupted_matrix = [list(row) for row in matrix]
    corrupted_matrix[0][4] *= -1
    if all(
        corrupted_matrix[row][column] == corrupted_matrix[column][row]
        for row in range(16)
        for column in range(16)
    ):
        raise AssertionError("transpose corruption was not detected")

    corrupted_pair = [
        [list(cloud) for cloud in state]
        for state in OBSTRUCTION_PAIRS[0]
    ]
    corrupted_pair[0][0][1] *= -1
    positive = tuple(tuple(cloud) for cloud in corrupted_pair[0])
    negative = tuple(tuple(cloud) for cloud in corrupted_pair[1])
    if cross_energy(positive) == 28 and all(
        projective_key(positive[cloud]) == projective_key(negative[cloud])
        for cloud in (0, 2, 3)
    ):
        raise AssertionError("six-state corruption was not detected")


def main() -> None:
    matrix = assemble_matrix(WITNESS_INTERNAL_EDGES)
    verify_structure(matrix)
    verify_order_four_blocks()
    verify_six_state_obstruction()
    histogram_digest, projective_maximizers = verify_witness(matrix)

    uniform_profile = common_internal_profile()
    expected_uniform_profile = Counter({42: 32, 40: 18, 46: 8, 38: 6})
    if uniform_profile != expected_uniform_profile:
        raise AssertionError(("common-internal profile", uniform_profile))
    verify_corruption_controls(matrix)

    print("order=16")
    print("cross_frame=common_oriented_H4")
    print("base_signs=+,+,+,-,+,-")
    print("internal_pattern=P,R,P,R")
    print("internal_maxima=4,4,4,4")
    print("P_R_switching_permutation_equivalent=TRUE")
    print("projective_spins_checked=32768")
    print("cross_only_maximum=28")
    print("lift_maximum=30")
    print(f"projective_maximizers={projective_maximizers}")
    print("full_maximizers=76")
    print(f"matrix_sha256={EXPECTED_MATRIX_SHA256}")
    print(f"energy_histogram_sha256={histogram_digest}")
    print("six_state_obstruction_cycle=0->3->2->0")
    print("fixed_H_arbitrary_internal_minimum=30")
    print("common_literal_internal_minimum=38")
    print("common_literal_internal_profile=M38:6,M40:18,M42:32,M46:8")
    print("corruption_controls=hadamard_entry,transpose,six_state_pair")
    print("framed_hadamard_lift_30_verification=PASSED")


if __name__ == "__main__":
    main()

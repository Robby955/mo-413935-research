#!/usr/bin/env python3
"""Exact verifier for a non-common-H, quotient-free order-16 cloud signing."""

from __future__ import annotations

import hashlib
import itertools
from functools import cache
from collections import Counter


if not __debug__:
    raise RuntimeError("assertions must be enabled")


Matrix = tuple[tuple[int, ...], ...]

HADAMARD_BLOCKS: tuple[Matrix, ...] = (
    ((-1, -1, 1, -1), (1, -1, 1, 1), (-1, -1, -1, 1), (-1, 1, 1, 1)),
    ((-1, 1, 1, 1), (-1, 1, -1, -1), (1, 1, 1, -1), (-1, -1, 1, -1)),
    ((1, 1, -1, 1), (-1, 1, 1, 1), (-1, -1, -1, 1), (-1, 1, -1, -1)),
    ((1, -1, -1, -1), (1, -1, 1, 1), (1, 1, -1, 1), (-1, -1, -1, 1)),
    ((1, 1, -1, 1), (-1, 1, 1, 1), (-1, 1, -1, -1), (-1, -1, -1, 1)),
    ((-1, -1, -1, 1), (-1, -1, 1, -1), (-1, 1, 1, 1), (1, -1, 1, 1)),
)

INTERNAL_BLOCKS: tuple[Matrix, ...] = (
    ((0, -1, 1, 1), (-1, 0, -1, -1), (1, -1, 0, 1), (1, -1, 1, 0)),
    ((0, -1, 1, -1), (-1, 0, 1, -1), (1, 1, 0, 1), (-1, -1, 1, 0)),
    ((0, 1, 1, -1), (1, 0, -1, 1), (1, -1, 0, -1), (-1, 1, -1, 0)),
    ((0, -1, -1, 1), (-1, 0, 1, 1), (-1, 1, 0, -1), (1, 1, -1, 0)),
)

CLOUD_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
HADAMARD_SEED: Matrix = (
    (1, 1, 1, 1),
    (1, -1, 1, -1),
    (1, 1, -1, -1),
    (1, -1, -1, 1),
)

CANONICAL_CLOUD_ORDER = (0, 1, 3, 2)
CANONICAL_COMMON_H: Matrix = (
    (1, 1, 1, -1),
    (-1, 1, -1, -1),
    (1, 1, -1, 1),
    (1, -1, -1, -1),
)
CANONICAL_GAUGES: dict[int, Matrix] = {
    0: ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)),
    1: ((-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, 0, -1), (0, 0, -1, 0)),
    3: ((0, 0, 1, 0), (0, 0, 0, -1), (-1, 0, 0, 0), (0, 1, 0, 0)),
    2: ((0, 0, 0, 1), (0, 1, 0, 0), (1, 0, 0, 0), (0, 0, 1, 0)),
}
CANONICAL_BASE_SIGNS = (1, 1, 1, -1, 1, -1)


def transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[row][column] for row in range(len(matrix))) for column in range(len(matrix[0])))


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum(left[row][inner] * right[inner][column] for inner in range(len(right)))
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def matrix_equal_up_to_sign(left: Matrix, right: Matrix) -> bool:
    return left == right or left == tuple(tuple(-value for value in row) for row in right)


def is_hadamard(matrix: Matrix) -> bool:
    product = matmul(matrix, transpose(matrix))
    return product == tuple(
        tuple(4 if row == column else 0 for column in range(4))
        for row in range(4)
    )


def is_signed_permutation(matrix: Matrix) -> bool:
    return all(
        sum(abs(value) for value in row) == 1 for row in matrix
    ) and all(
        sum(abs(matrix[row][column]) for row in range(4)) == 1
        for column in range(4)
    )


def scale_exact(matrix: Matrix, divisor: int) -> Matrix | None:
    if any(value % divisor for row in matrix for value in row):
        return None
    return tuple(tuple(value // divisor for value in row) for row in matrix)


@cache
def signed_hadamards() -> tuple[Matrix, ...]:
    matrices: set[Matrix] = set()
    for row_permutation in itertools.permutations(range(4)):
        for column_permutation in itertools.permutations(range(4)):
            permuted = tuple(
                tuple(
                    HADAMARD_SEED[row_permutation[row]][column_permutation[column]]
                    for column in range(4)
                )
                for row in range(4)
            )
            for row_signs in itertools.product((-1, 1), repeat=4):
                for column_signs in itertools.product((-1, 1), repeat=4):
                    matrices.add(
                        tuple(
                            tuple(
                                row_signs[row]
                                * permuted[row][column]
                                * column_signs[column]
                                for column in range(4)
                            )
                            for row in range(4)
                        )
                    )
    result = tuple(sorted(matrices))
    if len(result) != 768 or not all(is_hadamard(matrix) for matrix in result):
        raise AssertionError("order-four Hadamard catalogue is corrupted")
    return result


def assemble_matrix() -> Matrix:
    matrix = [[0 for _ in range(16)] for _ in range(16)]
    for cloud, internal in enumerate(INTERNAL_BLOCKS):
        for row in range(4):
            for column in range(4):
                matrix[4 * cloud + row][4 * cloud + column] = internal[row][column]
    for block, (left_cloud, right_cloud) in zip(HADAMARD_BLOCKS, CLOUD_PAIRS):
        for row in range(4):
            for column in range(4):
                matrix[4 * left_cloud + row][4 * right_cloud + column] = block[row][column]
                matrix[4 * right_cloud + column][4 * left_cloud + row] = block[row][column]
    return tuple(tuple(row) for row in matrix)


def quadratic_energy(matrix: Matrix, spins: tuple[int, ...]) -> int:
    return sum(
        matrix[row][column] * spins[row] * spins[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def block_energy(spins: tuple[int, ...], include_internal: bool = True) -> int:
    clouds = tuple(spins[4 * cloud : 4 * cloud + 4] for cloud in range(4))
    energy = 0
    if include_internal:
        for cloud, block in enumerate(INTERNAL_BLOCKS):
            energy += quadratic_energy(block, clouds[cloud])
    for block, (left_cloud, right_cloud) in zip(HADAMARD_BLOCKS, CLOUD_PAIRS):
        energy += sum(
            clouds[left_cloud][row]
            * block[row][column]
            * clouds[right_cloud][column]
            for row in range(4)
            for column in range(4)
        )
    return energy


def cross_edge_energies(spins: tuple[int, ...]) -> tuple[int, ...]:
    clouds = tuple(spins[4 * cloud : 4 * cloud + 4] for cloud in range(4))
    return tuple(
        sum(
            clouds[left_cloud][row]
            * block[row][column]
            * clouds[right_cloud][column]
            for row in range(4)
            for column in range(4)
        )
        for block, (left_cloud, right_cloud) in zip(HADAMARD_BLOCKS, CLOUD_PAIRS)
    )


def projective_spins(order: int):
    for mask in range(1 << (order - 1)):
        yield (1,) + tuple(
            -1 if (mask >> (coordinate - 1)) & 1 else 1
            for coordinate in range(1, order)
        )


def recompute_f4() -> int:
    edges = tuple(itertools.combinations(range(4), 2))
    optimum = 100
    for signing in itertools.product((-1, 1), repeat=len(edges)):
        maximum = 0
        for spins in projective_spins(4):
            energy = sum(
                coefficient * spins[left] * spins[right]
                for coefficient, (left, right) in zip(signing, edges)
            )
            maximum = max(maximum, abs(energy))
        optimum = min(optimum, maximum)
    return optimum


def fixed_template_profile(
    state_terms: tuple[tuple[int, tuple[int, ...]], ...]
) -> Counter[tuple[int, int]]:
    """Test the non-adaptive operation that only flips whole cross blocks."""

    profile: Counter[tuple[int, int]] = Counter()
    for coefficients in itertools.product((-1, 1), repeat=6):
        base_maximum = max(
            abs(
                sum(
                    coefficient * spins[left] * spins[right]
                    for coefficient, (left, right) in zip(
                        coefficients, CLOUD_PAIRS
                    )
                )
            )
            for spins in projective_spins(4)
        )
        lift_maximum = max(
            abs(
                internal
                + sum(
                    coefficient * term
                    for coefficient, term in zip(coefficients, terms)
                )
            )
            for internal, terms in state_terms
        )
        profile[(base_maximum, lift_maximum)] += 1
    return profile


def directed_block(
    by_pair: dict[tuple[int, int], Matrix], left: int, right: int
) -> Matrix:
    if left < right:
        return by_pair[(left, right)]
    return transpose(by_pair[(right, left)])


def common_h_solution_count_for_order(
    blocks: tuple[Matrix, ...],
    cloud_order: tuple[int, ...],
    *,
    require_symmetric_h: bool = False,
    require_common_internal: bool = False,
) -> int:
    """Count common-H gauges for one ordering of the four clouds.

    Reordering matters when H is nonsymmetric: every inverted edge uses H^T.
    Independent base-edge signs are allowed.
    """

    by_pair = dict(zip(CLOUD_PAIRS, blocks))
    identity = tuple(
        tuple(1 if row == column else 0 for column in range(4))
        for row in range(4)
    )
    solutions = 0
    for common in signed_hadamards():
        if require_symmetric_h and common != transpose(common):
            continue
        gauges: dict[int, Matrix] = {cloud_order[0]: identity}
        valid = True
        for cloud in cloud_order[1:]:
            # Star-edge signs may be made positive by replacing a whole cloud
            # gauge R_cloud with -R_cloud.
            star = directed_block(by_pair, cloud_order[0], cloud)
            numerator = matmul(transpose(star), common)
            gauge = scale_exact(numerator, 4)
            if gauge is None or not is_signed_permutation(gauge):
                valid = False
                break
            gauges[cloud] = gauge
        if not valid:
            continue
        for left_position in range(1, 4):
            for right_position in range(left_position + 1, 4):
                left = cloud_order[left_position]
                right = cloud_order[right_position]
                transformed = matmul(
                    matmul(
                        transpose(gauges[left]),
                        directed_block(by_pair, left, right),
                    ),
                    gauges[right],
                )
                if not matrix_equal_up_to_sign(transformed, common):
                    valid = False
                    break
            if not valid:
                break
        if valid and require_common_internal:
            common_internal = INTERNAL_BLOCKS[cloud_order[0]]
            for cloud in cloud_order[1:]:
                transformed_internal = matmul(
                    matmul(transpose(gauges[cloud]), INTERNAL_BLOCKS[cloud]),
                    gauges[cloud],
                )
                if not matrix_equal_up_to_sign(
                    transformed_internal, common_internal
                ):
                    valid = False
                    break
        if valid:
            solutions += 1
    return solutions


def common_h_solution_profile(
    blocks: tuple[Matrix, ...],
    *,
    require_symmetric_h: bool = False,
    require_common_internal: bool = False,
) -> dict[tuple[int, ...], int]:
    return {
        cloud_order: common_h_solution_count_for_order(
            blocks,
            cloud_order,
            require_symmetric_h=require_symmetric_h,
            require_common_internal=require_common_internal,
        )
        for cloud_order in itertools.permutations(range(4))
    }


def base_signing_maximum(coefficients: tuple[int, ...]) -> int:
    return max(
        abs(
            sum(
                coefficient * spins[left] * spins[right]
                for coefficient, (left, right) in zip(
                    coefficients, CLOUD_PAIRS
                )
            )
        )
        for spins in projective_spins(4)
    )


def verify_canonical_common_h_representation() -> None:
    by_pair = dict(zip(CLOUD_PAIRS, HADAMARD_BLOCKS))
    observed_signs: list[int] = []
    for left_position in range(4):
        for right_position in range(left_position + 1, 4):
            left = CANONICAL_CLOUD_ORDER[left_position]
            right = CANONICAL_CLOUD_ORDER[right_position]
            transformed = matmul(
                matmul(
                    transpose(CANONICAL_GAUGES[left]),
                    directed_block(by_pair, left, right),
                ),
                CANONICAL_GAUGES[right],
            )
            if transformed == CANONICAL_COMMON_H:
                observed_signs.append(1)
            elif transformed == tuple(
                tuple(-value for value in row) for row in CANONICAL_COMMON_H
            ):
                observed_signs.append(-1)
            else:
                raise AssertionError("canonical common-H representation failed")
    if tuple(observed_signs) != CANONICAL_BASE_SIGNS:
        raise AssertionError(("canonical base signs", tuple(observed_signs)))
    if not is_hadamard(CANONICAL_COMMON_H):
        raise AssertionError("canonical common H is not Hadamard")
    if CANONICAL_COMMON_H == transpose(CANONICAL_COMMON_H):
        raise AssertionError("canonical common H unexpectedly became symmetric")
    if base_signing_maximum(CANONICAL_BASE_SIGNS) != 4:
        raise AssertionError("canonical base signing is not F(4)-optimal")


def verify_structure(matrix: Matrix) -> None:
    if any(len(row) != 16 for row in matrix):
        raise AssertionError("assembled matrix is not square")
    for row in range(16):
        if matrix[row][row] != 0:
            raise AssertionError("assembled matrix has nonzero diagonal")
        for column in range(16):
            if matrix[row][column] != matrix[column][row]:
                raise AssertionError("assembled matrix is not symmetric")
            if row != column and matrix[row][column] not in (-1, 1):
                raise AssertionError("assembled matrix has a nonsign entry")
    if not all(is_hadamard(block) for block in HADAMARD_BLOCKS):
        raise AssertionError("a cross block is not Hadamard")
    if any(sum(map(sum, block)) != 0 for block in HADAMARD_BLOCKS):
        raise AssertionError("a cross block has nonzero quotient sum")
    if any(
        sum(block[row][column] for row in range(4) for column in range(row + 1, 4))
        != 0
        for block in INTERNAL_BLOCKS
    ):
        raise AssertionError("an internal block is not fixed-half")
    if sum(matrix[row][column] for row in range(16) for column in range(row + 1, 16)) != 0:
        raise AssertionError("the order-16 signing is not fixed-half")


def verify_corruption_controls(matrix: Matrix) -> None:
    corrupted_block = [list(row) for row in HADAMARD_BLOCKS[0]]
    corrupted_block[0][0] *= -1
    if is_hadamard(tuple(tuple(row) for row in corrupted_block)):
        raise AssertionError("Hadamard-entry corruption was not detected")

    corrupted_matrix = [list(row) for row in matrix]
    corrupted_matrix[4][0] *= -1
    if all(
        corrupted_matrix[row][column] == corrupted_matrix[column][row]
        for row in range(16)
        for column in range(16)
    ):
        raise AssertionError("lower-block transpose corruption was not detected")

    synthetic = tuple(
        tuple(
            tuple(
                (1 if (left + right) % 2 == 0 else -1)
                * HADAMARD_SEED[row][column]
                for column in range(4)
            )
            for row in range(4)
        )
        for left, right in CLOUD_PAIRS
    )
    if not any(common_h_solution_profile(synthetic).values()):
        raise AssertionError("common-H detector rejected a signed common system")


def main() -> None:
    matrix = assemble_matrix()
    verify_structure(matrix)

    histogram: Counter[int] = Counter()
    cross_maximum = 0
    independent_checks = 0
    state_terms: list[tuple[int, tuple[int, ...]]] = []
    for spins in projective_spins(16):
        direct = quadratic_energy(matrix, spins)
        blocked = block_energy(spins)
        if direct != blocked:
            raise AssertionError("direct and block energy computations disagree")
        histogram[direct] += 1
        terms = cross_edge_energies(spins)
        cross_energy = sum(terms)
        cross_maximum = max(cross_maximum, abs(cross_energy))
        state_terms.append((blocked - cross_energy, terms))
        independent_checks += 1

    maximum = max(abs(energy) for energy in histogram)
    maximizers = sum(count for energy, count in histogram.items() if abs(energy) == maximum)
    histogram_digest = hashlib.sha256(repr(sorted(histogram.items())).encode()).hexdigest()
    expected_digest = "b4950bbfbac3a6223c6260572ce698fd54168777da04adfa317c95456baa38ca"
    if maximum != 32 or maximizers != 14 or cross_maximum != 28:
        raise AssertionError(("unexpected lift maxima", maximum, maximizers, cross_maximum))
    if histogram_digest != expected_digest:
        raise AssertionError(("energy histogram digest", histogram_digest))

    f4 = recompute_f4()
    if f4 != 4 or maximum != 8 * f4:
        raise AssertionError(("lossless-scale benchmark", f4, maximum))
    internal_maxima = tuple(
        max(
            abs(quadratic_energy(INTERNAL_BLOCKS[cloud], spins))
            for spins in projective_spins(4)
        )
        for cloud in CANONICAL_CLOUD_ORDER
    )
    if internal_maxima != (6, 6, 4, 4):
        raise AssertionError(("canonical internal maxima", internal_maxima))
    common_profile = common_h_solution_profile(HADAMARD_BLOCKS)
    nonzero_common_profile = {
        order: count for order, count in common_profile.items() if count
    }
    expected_common_profile = {
        (0, 1, 3, 2): 2,
        (0, 2, 3, 1): 2,
        (1, 0, 2, 3): 2,
        (1, 3, 2, 0): 2,
        (2, 0, 1, 3): 2,
        (2, 3, 1, 0): 2,
        (3, 1, 0, 2): 2,
        (3, 2, 0, 1): 2,
    }
    if nonzero_common_profile != expected_common_profile:
        raise AssertionError(("common-H cloud-order profile", nonzero_common_profile))
    common_h_common_d_profile = common_h_solution_profile(
        HADAMARD_BLOCKS,
        require_common_internal=True,
    )
    if any(common_h_common_d_profile.values()):
        raise AssertionError(
            (
                "a common-H gauge also made the internal blocks common up to sign",
                common_h_common_d_profile,
            )
        )
    verify_canonical_common_h_representation()
    template_profile = fixed_template_profile(tuple(state_terms))
    expected_profile = Counter({(4, 40): 40, (4, 38): 8, (6, 40): 8, (6, 32): 8})
    if template_profile != expected_profile:
        raise AssertionError(("fixed-template profile", template_profile))
    verify_corruption_controls(matrix)

    print("clouds=4x4")
    print("quotient_free_hadamard_blocks=6")
    print("fixed_half_internal_blocks=4")
    print(f"projective_spins_checked={independent_checks}")
    print("base_F4=4")
    print("cross_only_maximum=28")
    print("lift_maximum=32")
    print("projective_maximizers=14")
    print(f"energy_histogram_sha256={histogram_digest}")
    print("order4_hadamards_checked=768")
    print(
        "common_H_nonzero_cloud_orders="
        + ";".join(
            "".join(map(str, order)) + f":{count}"
            for order, count in nonzero_common_profile.items()
        )
    )
    print("common_H_total_gauge_solutions=16")
    print("canonical_common_H_cloud_order=0132")
    print("canonical_common_H_is_symmetric=FALSE")
    print("canonical_base_signs=+,+,+,-,+,-")
    print("canonical_base_maximum=4")
    print("canonical_internal_maxima=6,6,4,4")
    print("common_H_common_D_up_to_sign_gauge_solutions=0")
    print(
        "stored_orientation_fixed_template_profile="
        "M4_to_38:8,M4_to_40:40,M6_to_32:8,M6_to_40:8"
    )
    print("corruption_controls=hadamard_entry,lower_block_transpose,synthetic_common_H")
    print("cloud_dependent_hadamard_lift_verification=PASSED")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact small-case checks for relative-switching profile composition.

Every signing in block sizes 2+3 is exhausted.  The script checks the
balanced product-to-gauge map, the exact fiber maximum identity, the
microcanonical order-statistic bound, and its exponential-profile corollary.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product


Spin = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


def projective_spins(order: int) -> tuple[Spin, ...]:
    return tuple((1, *tail) for tail in product((-1, 1), repeat=order - 1))


def symmetric_signings(order: int) -> tuple[Matrix, ...]:
    edges = [(left, right) for left in range(order) for right in range(left + 1, order)]
    result = []
    for signs in product((-1, 1), repeat=len(edges)):
        matrix = [[0] * order for _ in range(order)]
        for (left, right), sign in zip(edges, signs):
            matrix[left][right] = sign
            matrix[right][left] = sign
        result.append(tuple(tuple(row) for row in matrix))
    return tuple(result)


def symmetric_matrix(order: int, signs: tuple[int, ...]) -> Matrix:
    edges = [(left, right) for left in range(order) for right in range(left + 1, order)]
    assert len(edges) == len(signs)
    matrix = [[0] * order for _ in range(order)]
    for (left, right), sign in zip(edges, signs):
        matrix[left][right] = sign
        matrix[right][left] = sign
    return tuple(tuple(row) for row in matrix)


def rectangular_signings(rows: int, columns: int) -> tuple[Matrix, ...]:
    result = []
    for signs in product((-1, 1), repeat=rows * columns):
        matrix = tuple(
            tuple(signs[row * columns + column] for column in range(columns))
            for row in range(rows)
        )
        result.append(matrix)
    return tuple(result)


def quadratic_energy(matrix: Matrix, spin: Spin) -> int:
    return sum(
        matrix[left][right] * spin[left] * spin[right]
        for left in range(len(spin))
        for right in range(left + 1, len(spin))
    )


def rectangular_energy(matrix: Matrix, left: Spin, right: Spin) -> int:
    return sum(
        matrix[row][column] * left[row] * right[column]
        for row in range(len(left))
        for column in range(len(right))
    )


def augmented_states(matrix: Matrix) -> tuple[tuple[int, Spin, int], ...]:
    return tuple(
        (orientation, spin, orientation * quadratic_energy(matrix, spin))
        for orientation in (-1, 1)
        for spin in projective_spins(len(matrix))
    )


def cross_states(matrix: Matrix) -> tuple[tuple[Spin, Spin, int], ...]:
    return tuple(
        (left, right, abs(rectangular_energy(matrix, left, right)))
        for left in projective_spins(len(matrix))
        for right in projective_spins(len(matrix[0]))
    )


def coordinate_product(left: Spin, right: Spin) -> Spin:
    result = tuple(a * b for a, b in zip(left, right))
    assert result[0] == 1
    return result


def aligned_full_matrix(
    left_matrix: Matrix,
    right_matrix: Matrix,
    cross_matrix: Matrix,
    left_gauge: Spin,
    right_gauge: Spin,
    relative_sign: int,
) -> Matrix:
    left_order = len(left_matrix)
    right_order = len(right_matrix)
    order = left_order + right_order
    result = [[0] * order for _ in range(order)]
    for row in range(left_order):
        for column in range(left_order):
            result[row][column] = (
                left_gauge[row]
                * left_matrix[row][column]
                * left_gauge[column]
            )
    for row in range(right_order):
        for column in range(right_order):
            result[left_order + row][left_order + column] = (
                relative_sign
                * right_gauge[row]
                * right_matrix[row][column]
                * right_gauge[column]
            )
    for row in range(left_order):
        for column in range(right_order):
            value = cross_matrix[row][column]
            result[row][left_order + column] = value
            result[left_order + column][row] = value
    return tuple(tuple(row) for row in result)


def maximum_energy(matrix: Matrix) -> int:
    return max(
        abs(quadratic_energy(matrix, spin))
        for spin in projective_spins(len(matrix))
    )


def powers_of_two_profile(deficits: list[int]) -> Fraction:
    return sum((Fraction(1, 2**deficit) for deficit in deficits), Fraction())


def check_case(
    left_matrix: Matrix,
    right_matrix: Matrix,
    cross_matrix: Matrix,
    run_profile_checks: bool,
) -> tuple[int, int, int, bool, int, int]:
    left_states = augmented_states(left_matrix)
    right_states = augmented_states(right_matrix)
    rectangle_states = cross_states(cross_matrix)

    left_maximum = max(state[2] for state in left_states)
    right_maximum = max(state[2] for state in right_states)
    rectangle_maximum = max(state[2] for state in rectangle_states)
    independent_ceiling = left_maximum + right_maximum + rectangle_maximum

    fibers: dict[tuple[Spin, Spin, int], list[int]] = {}
    all_deficits = []
    for left_orientation, left_spin, left_energy in left_states:
        for right_orientation, right_spin, right_energy in right_states:
            relative_sign = left_orientation * right_orientation
            for rectangle_left, rectangle_right, rectangle_energy in rectangle_states:
                gauge = (
                    coordinate_product(left_spin, rectangle_left),
                    coordinate_product(right_spin, rectangle_right),
                    relative_sign,
                )
                total_energy = left_energy + right_energy + rectangle_energy
                fibers.setdefault(gauge, []).append(total_energy)
                all_deficits.append(independent_ceiling - total_energy)

    total_order = len(left_matrix) + len(right_matrix)
    gauge_count = 2 ** (total_order - 1)
    fiber_size = 2 ** (total_order - 1)
    assert len(fibers) == gauge_count
    assert all(len(values) == fiber_size for values in fibers.values())

    fiber_identity_checks = 0
    fiber_deficits = []
    for (left_gauge, right_gauge, relative_sign), energies in fibers.items():
        full_matrix = aligned_full_matrix(
            left_matrix,
            right_matrix,
            cross_matrix,
            left_gauge,
            right_gauge,
            relative_sign,
        )
        full_maximum = maximum_energy(full_matrix)
        assert max(energies) == full_maximum
        fiber_deficits.append(independent_ceiling - full_maximum)
        fiber_identity_checks += 1

    all_deficits.sort()
    order_statistic = all_deficits[gauge_count - 1]
    best_fiber_deficit = max(fiber_deficits)
    assert best_fiber_deficit >= order_statistic
    assert (
        min(max(energies) for energies in fibers.values())
        <= independent_ceiling - order_statistic
    )

    profile_checks = 0
    if run_profile_checks:
        left_deficits = [left_maximum - state[2] for state in left_states]
        right_deficits = [right_maximum - state[2] for state in right_states]
        rectangle_deficits = [
            rectangle_maximum - state[2] for state in rectangle_states
        ]
        profile_product = (
            powers_of_two_profile(left_deficits)
            * powers_of_two_profile(right_deficits)
            * powers_of_two_profile(rectangle_deficits)
        )
        for threshold in range(independent_ceiling + 2):
            count = sum(deficit < threshold for deficit in all_deficits)
            assert Fraction(count) <= (2**threshold) * profile_product
            profile_checks += 1

    return (
        fiber_identity_checks,
        len(all_deficits),
        profile_checks,
        best_fiber_deficit > order_statistic,
        order_statistic,
        best_fiber_deficit,
    )


def verify_exhaustive_two_plus_three() -> None:
    cases = 0
    fiber_checks = 0
    triple_checks = 0
    profile_checks = 0
    strict_geometry_gains = 0
    for left_matrix in symmetric_signings(2):
        for right_matrix in symmetric_signings(3):
            for cross_index, cross_matrix in enumerate(rectangular_signings(2, 3)):
                checked = check_case(
                    left_matrix,
                    right_matrix,
                    cross_matrix,
                    run_profile_checks=(cross_index % 17 == 0),
                )
                (
                    case_fibers,
                    case_triples,
                    case_profiles,
                    strict_gain,
                    _order_statistic,
                    _best_fiber_deficit,
                ) = checked
                cases += 1
                fiber_checks += case_fibers
                triple_checks += case_triples
                profile_checks += case_profiles
                strict_geometry_gains += int(strict_gain)

    assert cases == 1024
    assert fiber_checks == 16384
    assert triple_checks == 262144
    assert strict_geometry_gains > 0
    print(f"block_cases_checked={cases}")
    print(f"balanced_fibers_checked={fiber_checks}")
    print(f"product_triples_checked={triple_checks}")
    print(f"exponential_profile_checks={profile_checks}")
    print(f"strict_geometry_gains={strict_geometry_gains}")


def verify_two_plus_four_collision() -> None:
    left_matrix = symmetric_matrix(2, (1,))
    right_matrix = symmetric_matrix(4, (-1, -1, -1, -1, -1, 1))
    cross_high = (
        (-1, -1, -1, -1),
        (-1, -1, 1, 1),
    )
    cross_low = (
        (-1, -1, -1, -1),
        (-1, 1, -1, 1),
    )
    high = check_case(left_matrix, right_matrix, cross_high, False)
    low = check_case(left_matrix, right_matrix, cross_low, False)
    assert high[-2:] == (0, 4)
    assert low[-2:] == (0, 2)
    print("scalar_collision=lambda:0,0 true_gains:4,2")


def verify_corruption_control() -> None:
    left_matrix = symmetric_signings(2)[0]
    right_matrix = symmetric_signings(3)[0]
    cross_matrix = rectangular_signings(2, 3)[0]
    left_states = augmented_states(left_matrix)
    right_states = augmented_states(right_matrix)
    rectangle_states = cross_states(cross_matrix)

    corrupted: dict[tuple[Spin, Spin, int], list[int]] = {}
    for left_orientation, left_spin, left_energy in left_states:
        for _right_orientation, right_spin, right_energy in right_states:
            for rectangle_left, rectangle_right, rectangle_energy in rectangle_states:
                gauge = (
                    coordinate_product(left_spin, rectangle_left),
                    coordinate_product(right_spin, rectangle_right),
                    left_orientation,  # Wrong: the right orientation is omitted.
                )
                corrupted.setdefault(gauge, []).append(
                    left_energy + right_energy + rectangle_energy
                )

    mismatch_detected = False
    for (left_gauge, right_gauge, relative_sign), energies in corrupted.items():
        full_matrix = aligned_full_matrix(
            left_matrix,
            right_matrix,
            cross_matrix,
            left_gauge,
            right_gauge,
            relative_sign,
        )
        if max(energies) != maximum_energy(full_matrix):
            mismatch_detected = True
            break
    assert mismatch_detected
    print("corruption_control=relative_orientation_omission_detected")


def main() -> None:
    verify_exhaustive_two_plus_three()
    verify_two_plus_four_collision()
    verify_corruption_control()
    print("relative_profile_composition_verification=PASSED")


if __name__ == "__main__":
    main()

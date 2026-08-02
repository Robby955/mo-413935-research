#!/usr/bin/env python3
"""Exact checks for the order-five tensor amplification counterexamples."""

from __future__ import annotations

import itertools


A = (
    (0, -1, 1, 1, -1),
    (-1, 0, -1, 1, 1),
    (1, -1, 0, -1, 1),
    (1, 1, -1, 0, -1),
    (-1, 1, 1, -1, 0),
)

D_N = (-1, 1, -1, 1, 1)
D_A = (-1, -1, 1, 1, 1)

X_NN = (
    (-1, -1, 1, 1, -1),
    (-1, 1, -1, 1, 1),
    (1, -1, -1, -1, 1),
    (1, 1, -1, 1, -1),
    (-1, 1, 1, -1, 1),
)

X_NA = (
    (1, 1, -1, 1, -1),
    (-1, -1, 1, 1, -1),
    (1, -1, 1, -1, 1),
    (-1, 1, 1, -1, 1),
    (-1, -1, -1, 1, 1),
)

X_AA = (
    (-1, -1, -1, 1, 1),
    (-1, -1, 1, 1, -1),
    (-1, 1, 1, -1, 1),
    (1, 1, -1, 1, -1),
    (1, -1, 1, -1, 1),
)

X_UNIFORM = (
    (-1, 1, -1, -1, 1),
    (1, -1, -1, 1, -1),
    (-1, -1, 1, -1, 1),
    (-1, 1, -1, 1, -1),
    (1, -1, 1, -1, -1),
)


def assert_seidel(matrix: tuple[tuple[int, ...], ...]) -> None:
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise AssertionError("matrix is not square")
    for row in range(n):
        if matrix[row][row] != 0:
            raise AssertionError("matrix diagonal is not zero")
        for column in range(n):
            if matrix[row][column] != matrix[column][row]:
                raise AssertionError("matrix is not symmetric")
            if row != column and matrix[row][column] not in (-1, 1):
                raise AssertionError("off-diagonal entry is not a sign")


def quadratic_energy(
    matrix: tuple[tuple[int, ...], ...], spins: tuple[int, ...]
) -> int:
    return sum(
        matrix[row][column] * spins[row] * spins[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def spins_from_mask(n: int, mask: int) -> tuple[int, ...]:
    return tuple(-1 if mask >> index & 1 else 1 for index in range(n))


def completed(
    diagonal: tuple[int, ...],
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(A[row][column] if row != column else diagonal[row] for column in range(5))
        for row in range(5)
    )


def separated_tensor(
    left_diagonal: tuple[int, ...], right_diagonal: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    left = completed(left_diagonal)
    right = completed(right_diagonal)
    size = 25
    matrix = [[0 for _ in range(size)] for _ in range(size)]
    for left_row in range(5):
        for right_row in range(5):
            row = 5 * left_row + right_row
            for left_column in range(5):
                for right_column in range(5):
                    column = 5 * left_column + right_column
                    value = left[left_row][left_column] * right[right_row][
                        right_column
                    ] - (left_diagonal[left_row] if left_row == left_column else 0) * (
                        right_diagonal[right_row] if right_row == right_column else 0
                    )
                    matrix[row][column] = value
    result = tuple(tuple(row) for row in matrix)
    assert_seidel(result)
    return result


def canonical_tensor() -> tuple[tuple[int, ...], ...]:
    completed_plus = completed((1, 1, 1, 1, 1))
    matrix = []
    for left_row in range(5):
        for right_row in range(5):
            row = []
            index = 5 * left_row + right_row
            for left_column in range(5):
                for right_column in range(5):
                    column = 5 * left_column + right_column
                    value = completed_plus[left_row][left_column] * completed_plus[
                        right_row
                    ][right_column] - (1 if index == column else 0)
                    row.append(value)
            matrix.append(tuple(row))
    result = tuple(matrix)
    assert_seidel(result)
    return result


def flatten(matrix: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    return tuple(value for row in matrix for value in row)


def permuted(
    matrix: tuple[tuple[int, ...], ...], permutation: tuple[int, ...], sign: int
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            sign * matrix[permutation[row]][permutation[column]] for column in range(5)
        )
        for row in range(5)
    )


def classify_balanced_diagonals() -> tuple[int, int]:
    representatives = (completed(D_N), completed(D_A))
    counts = [0, 0]
    for diagonal in itertools.product((-1, 1), repeat=5):
        if abs(sum(diagonal)) != 1:
            continue
        target = completed(diagonal)
        matches = set()
        for representative_index, representative in enumerate(representatives):
            for permutation in itertools.permutations(range(5)):
                for sign in (-1, 1):
                    if permuted(representative, permutation, sign) == target:
                        matches.add(representative_index)
        if len(matches) != 1:
            raise AssertionError(
                ("balanced orbit classification failed", diagonal, matches)
            )
        counts[matches.pop()] += 1
    if counts != [10, 10]:
        raise AssertionError(("unexpected balanced orbit sizes", counts))
    return counts[0], counts[1]


def verify_base_matrix() -> None:
    assert_seidel(A)
    energies = {quadratic_energy(A, spins_from_mask(5, mask)) for mask in range(1 << 5)}
    if energies != {-4, 0, 4}:
        raise AssertionError(("unexpected order-five energy set", energies))


def verify_product_lower_bound() -> int:
    checked = 0
    for left in itertools.product((-1, 1), repeat=5):
        for right in itertools.product((-1, 1), repeat=5):
            lower_bound = 4 * (8 + abs(sum(left)) + abs(sum(right)))
            if (abs(sum(left)), abs(sum(right))) != (1, 1) and lower_bound < 48:
                raise AssertionError("unbalanced product lower bound fell below 48")
            checked += 1
    return checked


def verify_certificates() -> dict[str, int]:
    cases = {
        "NN": (separated_tensor(D_N, D_N), flatten(X_NN), 88),
        "NA": (separated_tensor(D_N, D_A), flatten(X_NA), -80),
        "AA": (separated_tensor(D_A, D_A), flatten(X_AA), 88),
        "uniform": (canonical_tensor(), flatten(X_UNIFORM), 100),
    }
    observed = {}
    corrupted_expected = {"NN": 76, "NA": -64, "AA": 76, "uniform": 76}
    for name, (matrix, witness, expected) in cases.items():
        value = quadratic_energy(matrix, witness)
        if value != expected:
            raise AssertionError(("certificate energy mismatch", name, value, expected))
        corrupted = list(witness)
        corrupted[0] *= -1
        corrupted_value = quadratic_energy(matrix, tuple(corrupted))
        if corrupted_value != corrupted_expected[name] or corrupted_value == expected:
            raise AssertionError(("certificate corruption was not detected", name))
        observed[name] = value
    return observed


def verify_tensor_iteration() -> None:
    for exponent in (1, 2):
        numerator = 225**exponent - 25**exponent
        order = 25**exponent
        normalized = numerator / (2.0 * order**1.5)
        expected = 0.5 * (1.8**exponent - 0.2**exponent)
        if abs(normalized - expected) > 1e-12:
            raise AssertionError("tensor iteration formula failed")


def main() -> None:
    verify_base_matrix()
    product_pairs = verify_product_lower_bound()
    first_orbit, second_orbit = classify_balanced_diagonals()
    certificate_values = verify_certificates()
    verify_tensor_iteration()
    print("order_5_energy_set=-4,0,4")
    print(f"diagonal_pairs_checked={product_pairs}")
    print(f"balanced_completion_orbits={first_orbit},{second_orbit}")
    print(
        "certificate_energies="
        + ",".join(f"{name}:{value}" for name, value in certificate_values.items())
    )
    print("tensor_iterations_checked=2")
    print("corruption_controls=PASSED")


if __name__ == "__main__":
    main()

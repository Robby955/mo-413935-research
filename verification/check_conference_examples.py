#!/usr/bin/env python3
"""Exact checks for the small Paley conference examples used in STATUS.md."""

from __future__ import annotations

import math


def legendre_symbol(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def paley_conference(prime: int) -> list[list[int]]:
    if prime % 4 != 1:
        raise ValueError("prime must be 1 modulo 4")
    order = prime + 1
    matrix = [[0] * order for _ in range(order)]
    for column in range(1, order):
        matrix[0][column] = matrix[column][0] = 1
    for row in range(1, order):
        for column in range(row + 1, order):
            value = legendre_symbol((row - 1) - (column - 1), prime)
            matrix[row][column] = matrix[column][row] = value
    return matrix


def energy(matrix: list[list[int]], signs: list[int]) -> int:
    return sum(
        matrix[row][column] * signs[row] * signs[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def maximum_absolute_energy(matrix: list[list[int]]) -> int:
    order = len(matrix)
    best = 0
    for bits in range(1 << (order - 1)):
        signs = [1] + [
            1 if (bits >> (index - 1)) & 1 else -1
            for index in range(1, order)
        ]
        best = max(best, abs(energy(matrix, signs)))
    return best


def check_conference_identity(matrix: list[list[int]]) -> None:
    order = len(matrix)
    for row in range(order):
        for column in range(order):
            product = sum(
                matrix[row][index] * matrix[index][column]
                for index in range(order)
            )
            expected = order - 1 if row == column else 0
            if product != expected:
                raise AssertionError((row, column, product, expected))


def main() -> None:
    expected_maxima = {5: 5, 13: 21, 17: 33}
    for prime, expected in expected_maxima.items():
        matrix = paley_conference(prime)
        check_conference_identity(matrix)
        maximum = maximum_absolute_energy(matrix)
        if maximum != expected:
            raise AssertionError((prime, maximum, expected))
        order = prime + 1
        ratio = maximum / order**1.5
        ceiling = math.sqrt((order - 1) / order) / 2
        print(
            f"order={order} maximum={maximum} "
            f"ratio={ratio:.9f} spectral_ceiling={ceiling:.9f}"
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact finite checks for the square-order Paley Boolean eigenvector.

For an odd prime ``p``, this script constructs GF(p^2), the symmetric Paley
conference matrix of order p^2+1, and the sign eigenvector obtained by making
(p+1)/2 additive GF(p)-cosets positive.  It checks

    C^2 = p^2 I,       C x = p x,
    Q_C(x) = p(p^2+1)/2.

The proof in the note works for every odd prime power m by replacing GF(p)
with GF(m).  This script deliberately checks only prime m so that it remains
self-contained and does not depend on a finite-field package.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class QuadraticField:
    prime: int
    nonsquare: int

    def add(self, left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        return (
            (left[0] + right[0]) % self.prime,
            (left[1] + right[1]) % self.prime,
        )

    def negate(self, value: tuple[int, int]) -> tuple[int, int]:
        return (-value[0] % self.prime, -value[1] % self.prime)

    def multiply(
        self, left: tuple[int, int], right: tuple[int, int]
    ) -> tuple[int, int]:
        return (
            (left[0] * right[0] + self.nonsquare * left[1] * right[1])
            % self.prime,
            (left[0] * right[1] + left[1] * right[0]) % self.prime,
        )

    def power(self, value: tuple[int, int], exponent: int) -> tuple[int, int]:
        result = (1, 0)
        base = value
        while exponent:
            if exponent & 1:
                result = self.multiply(result, base)
            base = self.multiply(base, base)
            exponent >>= 1
        return result

    def quadratic_character(self, value: tuple[int, int]) -> int:
        if value == (0, 0):
            return 0
        result = self.power(value, (self.prime**2 - 1) // 2)
        if result == (1, 0):
            return 1
        if result == (self.prime - 1, 0):
            return -1
        raise AssertionError(("Euler criterion outside {+-1}", value, result))


def legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def least_nonsquare(prime: int) -> int:
    for value in range(2, prime):
        if legendre(value, prime) == -1:
            return value
    raise AssertionError(("no nonsquare", prime))


def paley_matrix(prime: int) -> tuple[list[list[int]], tuple[tuple[int, int], ...]]:
    field = QuadraticField(prime, least_nonsquare(prime))
    elements = tuple((a, b) for b in range(prime) for a in range(prime))
    order = prime**2 + 1
    matrix = [[0] * order for _ in range(order)]
    for index in range(1, order):
        matrix[0][index] = matrix[index][0] = 1
    for row, left in enumerate(elements, start=1):
        for column in range(row + 1, order):
            right = elements[column - 1]
            difference = field.add(left, field.negate(right))
            character = field.quadratic_character(difference)
            if character not in (-1, 1):
                raise AssertionError(("off-diagonal zero character", row, column))
            matrix[row][column] = matrix[column][row] = character
    return matrix, elements


def subfield_coset_vector(
    prime: int, elements: tuple[tuple[int, int], ...]
) -> list[int]:
    positive_cosets = (prime + 1) // 2
    quotient_signs = [1 if b < positive_cosets else -1 for b in range(prime)]
    if sum(quotient_signs) != 1:
        raise AssertionError("unbalanced quotient signs")
    return [1] + [quotient_signs[b] for _, b in elements]


def matrix_vector(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [
        sum(coefficient * value for coefficient, value in zip(row, vector, strict=True))
        for row in matrix
    ]


def check_conference_identity(matrix: list[list[int]], prime: int) -> None:
    order = len(matrix)
    for row in range(order):
        for column in range(order):
            product = sum(
                matrix[row][index] * matrix[index][column] for index in range(order)
            )
            expected = prime**2 if row == column else 0
            if product != expected:
                raise AssertionError(("conference identity", row, column, product))


def energy(matrix: list[list[int]], vector: list[int]) -> int:
    return sum(
        matrix[row][column] * vector[row] * vector[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def verify_prime(prime: int) -> tuple[int, int]:
    matrix, elements = paley_matrix(prime)
    vector = subfield_coset_vector(prime, elements)
    check_conference_identity(matrix, prime)
    image = matrix_vector(matrix, vector)
    expected_image = [prime * value for value in vector]
    if image != expected_image:
        raise AssertionError(("Boolean eigenvector", prime))
    witnessed_energy = energy(matrix, vector)
    spectral_ceiling = prime * (prime**2 + 1) // 2
    if witnessed_energy != spectral_ceiling:
        raise AssertionError(
            ("spectral ceiling witness", prime, witnessed_energy, spectral_ceiling)
        )

    corrupted = list(vector)
    corrupted[1] *= -1
    if matrix_vector(matrix, corrupted) == [prime * value for value in corrupted]:
        raise AssertionError("coset-constant corruption was not detected")
    return len(matrix), witnessed_energy


def main() -> None:
    records = []
    for prime in (3, 5, 7):
        order, witnessed_energy = verify_prime(prime)
        records.append(f"m={prime}:N={order}:M={witnessed_energy}")
    print("paley_subfield_witnesses=" + ",".join(records))
    print("corruption_controls=coset_constancy")
    print("paley_subfield_verification=PASSED")


if __name__ == "__main__":
    main()

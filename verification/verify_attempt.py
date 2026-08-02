#!/usr/bin/env python3
"""Independent finite checks for the MO 413935 partial-attempt note."""

from __future__ import annotations

import math
import random


def sign_vectors(length: int):
    for bits in range(1 << length):
        yield tuple(1 if (bits >> index) & 1 else -1 for index in range(length))


def edge_characters(signs: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        signs[row] * signs[column]
        for row in range(len(signs))
        for column in range(row + 1, len(signs))
    )


def hamming(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(a != b for a, b in zip(left, right, strict=True))


def verify_covering_radius_identity(max_order: int = 5) -> None:
    corruption_detected = False
    for order in range(2, max_order + 1):
        edge_count = order * (order - 1) // 2
        cut_words = {edge_characters(signs) for signs in sign_vectors(order)}
        augmented = cut_words | {tuple(-entry for entry in word) for word in cut_words}

        direct = math.inf
        covering_radius = 0
        unaugmented_radius = 0
        for coefficients in sign_vectors(edge_count):
            direct_norm = max(
                abs(sum(a * c for a, c in zip(coefficients, word, strict=True)))
                for word in cut_words
            )
            direct = min(direct, direct_norm)
            covering_radius = max(
                covering_radius,
                min(hamming(coefficients, word) for word in augmented),
            )
            unaugmented_radius = max(
                unaugmented_radius,
                min(hamming(coefficients, word) for word in cut_words),
            )

        reconstructed = edge_count - 2 * covering_radius
        if reconstructed != direct:
            raise AssertionError((order, direct, reconstructed))
        if edge_count - 2 * unaugmented_radius != direct:
            corruption_detected = True

    if not corruption_detected:
        raise AssertionError("dropping the negative codewords was not detected")


def matrix_square(matrix: list[list[int]]) -> list[list[int]]:
    order = len(matrix)
    return [
        [
            sum(matrix[row][index] * matrix[index][column] for index in range(order))
            for column in range(order)
        ]
        for row in range(order)
    ]


def verify_gaussian_correlation_inequality() -> None:
    generator = random.Random(413935)
    checked_edges = 0
    for order in range(3, 13):
        for _ in range(30):
            matrix = [[0] * order for _ in range(order)]
            for row in range(order):
                for column in range(row + 1, order):
                    value = generator.choice((-1, 1))
                    matrix[row][column] = matrix[column][row] = value

            square = matrix_square(matrix)
            parameter = math.sqrt(order / (order - 1))
            variance = 1 + parameter * parameter * (order - 1) / order
            v = 2 * parameter / (math.sqrt(order) * variance)

            for row in range(order):
                for column in range(row + 1, order):
                    coefficient = matrix[row][column]
                    u = (
                        coefficient
                        * parameter
                        * parameter
                        * square[row][column]
                        / (order * variance)
                    )
                    left = u - v
                    right = u + v
                    tolerance = 1e-12
                    if left < -1 - tolerance or right > 1 + tolerance:
                        raise AssertionError((order, row, column, left, right))
                    left = max(-1.0, left)
                    right = min(1.0, right)
                    difference = math.asin(right) - math.asin(left)
                    if difference + tolerance < 2 * v:
                        raise AssertionError((order, row, column, difference, 2 * v))
                    checked_edges += 1

    small_v = 1e-4
    if math.asin(small_v) - math.asin(-small_v) >= 2.01 * small_v:
        raise AssertionError("strengthened derivative corruption was not detected")
    print(f"gaussian_edges_checked={checked_edges}")


def verify_optimized_bound() -> None:
    for order in range(2, 10_000):
        edge_count = order * (order - 1) / 2
        parameter = math.sqrt(order / (order - 1))
        variance = 1 + parameter * parameter * (order - 1) / order
        derived = 4 * parameter * edge_count / (math.pi * math.sqrt(order) * variance)
        target = order * math.sqrt(order - 1) / math.pi
        if not math.isclose(derived, target, rel_tol=1e-13, abs_tol=1e-13):
            raise AssertionError((order, derived, target))


def main() -> None:
    verify_covering_radius_identity()
    verify_gaussian_correlation_inequality()
    verify_optimized_bound()
    print("covering_radius_identity=VERIFIED_orders_2_through_5")
    print("optimized_lower_bound=VERIFIED_orders_2_through_9999")
    print("corruption_controls=PASSED")


if __name__ == "__main__":
    main()

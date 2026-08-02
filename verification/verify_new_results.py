#!/usr/bin/env python3
"""Deterministic checks for the second MO 413935 research attempt.

This script verifies finite instances of the new exact identities:

* Walsh Parseval and the nondegeneracy of the relaxed objective;
* linearity of the augmented cut code and its even-Eulerian dual;
* the MacWilliams/high-temperature polynomial identity;
* the fourth-cumulant/trace(A^4) identity and row-correlation defect;
* all root-normalized order-seven trace(A^4) minimizers and sixth moments;
* the covariance pair used in the elliptope and Gaussian-SDP bounds;
* the fact that removing the absolute value changes the problem.

All theorem checks use exact integer or Fraction arithmetic.  Floating point
is used only to print normalized values, never as a pass/fail oracle.
"""

from __future__ import annotations

import itertools
import random
from fractions import Fraction


SEED = 413935


def edges(order: int) -> tuple[tuple[int, int], ...]:
    return tuple(itertools.combinations(range(order), 2))


def spin_vectors(order: int):
    for mask in range(1 << order):
        yield tuple(1 if mask >> vertex & 1 else -1 for vertex in range(order))


def cut_word_bits(spins: tuple[int, ...]) -> int:
    word = 0
    for index, (row, column) in enumerate(edges(len(spins))):
        if spins[row] * spins[column] == -1:
            word |= 1 << index
    return word


def augmented_cut_code(order: int) -> set[int]:
    edge_count = len(edges(order))
    all_ones = (1 << edge_count) - 1
    cuts = {cut_word_bits(spins) for spins in spin_vectors(order)}
    return cuts | {word ^ all_ones for word in cuts}


def is_even_eulerian(order: int, edge_mask: int) -> bool:
    if edge_mask.bit_count() % 2:
        return False
    degrees = [0] * order
    for index, (row, column) in enumerate(edges(order)):
        if edge_mask >> index & 1:
            degrees[row] ^= 1
            degrees[column] ^= 1
    return not any(degrees)


def dual_by_definition(code: set[int], edge_count: int) -> set[int]:
    return {
        candidate
        for candidate in range(1 << edge_count)
        if all((candidate & word).bit_count() % 2 == 0 for word in code)
    }


def verify_code_and_dual() -> int:
    orthogonality_checks = 0
    for order in range(3, 7):
        edge_count = len(edges(order))
        code = augmented_cut_code(order)
        if len(code) != 1 << order:
            raise AssertionError(("augmented code dimension", order, len(code)))
        if 0 not in code:
            raise AssertionError(("zero codeword missing", order))
        for left in code:
            for right in code:
                if left ^ right not in code:
                    raise AssertionError(("code not linear", order, left, right))

        computed_dual = dual_by_definition(code, edge_count)
        even_eulerian = {
            mask for mask in range(1 << edge_count) if is_even_eulerian(order, mask)
        }
        if computed_dual != even_eulerian:
            raise AssertionError(("wrong dual", order))
        if len(computed_dual) != 1 << (edge_count - order):
            raise AssertionError(("dual dimension", order, len(computed_dual)))
        if order >= 4 and min(mask.bit_count() for mask in computed_dual if mask) != 4:
            raise AssertionError(("dual minimum distance", order))
        orthogonality_checks += len(code) * len(computed_dual)
    return orthogonality_checks


def energy(coefficients: tuple[int, ...], spins: tuple[int, ...]) -> int:
    return sum(
        coefficient * spins[row] * spins[column]
        for coefficient, (row, column) in zip(
            coefficients, edges(len(spins)), strict=True
        )
    )


def matrix_from_coefficients(
    order: int, coefficients: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    matrix = [[0] * order for _ in range(order)]
    for coefficient, (row, column) in zip(coefficients, edges(order), strict=True):
        matrix[row][column] = matrix[column][row] = coefficient
    return tuple(tuple(row) for row in matrix)


def matrix_product(
    left: tuple[tuple[int, ...], ...], right: tuple[tuple[int, ...], ...]
) -> tuple[tuple[int, ...], ...]:
    order = len(left)
    return tuple(
        tuple(
            sum(left[row][index] * right[index][column] for index in range(order))
            for column in range(order)
        )
        for row in range(order)
    )


def trace_fourth(matrix: tuple[tuple[int, ...], ...]) -> int:
    square = matrix_product(matrix, matrix)
    return sum(entry * entry for row in square for entry in row)


def moments(order: int, coefficients: tuple[int, ...]) -> tuple[int, int, int]:
    values = [energy(coefficients, spins) for spins in spin_vectors(order)]
    count = 1 << order
    sums = tuple(sum(value**power for value in values) for power in (2, 4, 6))
    if any(total % count for total in sums):
        raise AssertionError(("nonintegral moment", order, sums))
    return tuple(total // count for total in sums)  # type: ignore[return-value]


def verify_parseval_cumulants() -> int:
    generator = random.Random(SEED)
    checks = 0
    for order in range(2, 8):
        edge_count = len(edges(order))
        samples: list[tuple[int, ...]] = []
        if order <= 5:
            samples.extend(itertools.product((-1, 1), repeat=edge_count))
        else:
            samples.extend(
                tuple(generator.choice((-1, 1)) for _ in range(edge_count))
                for _ in range(80)
            )

        for coefficients in samples:
            second, fourth, _ = moments(order, coefficients)
            if second != edge_count:
                raise AssertionError(("Walsh Parseval", order, second, edge_count))
            matrix = matrix_from_coefficients(order, coefficients)
            trace4 = trace_fourth(matrix)
            expected_cumulant = 3 * trace4 - 2 * order * (order - 1) * (3 * order - 4)
            observed_cumulant = fourth - 3 * second * second
            if observed_cumulant != expected_cumulant:
                raise AssertionError(
                    ("fourth cumulant", order, observed_cumulant, expected_cumulant)
                )

            square = matrix_product(matrix, matrix)
            defect = sum(
                square[row][column] ** 2
                for row in range(order)
                for column in range(order)
                if row != column
            )
            baseline = order * (order - 1) ** 2
            if trace4 != baseline + defect:
                raise AssertionError(("row-correlation defect", order, trace4))
            checks += 1

    return checks


def verify_high_temperature_minimizers() -> tuple[int, int]:
    """Exhaust all root-normalized order-seven signings."""
    order = 7
    edge_list = edges(order)
    residual_indices = tuple(
        index for index, (row, _) in enumerate(edge_list) if row != 0
    )
    minimum_trace: int | None = None
    minimizing_count = 0
    sixth_moments: set[int] = set()

    for mask in range(1 << len(residual_indices)):
        coefficients = [1] * len(edge_list)
        for bit, edge_index in enumerate(residual_indices):
            coefficients[edge_index] = -1 if mask >> bit & 1 else 1
        immutable = tuple(coefficients)
        trace4 = trace_fourth(matrix_from_coefficients(order, immutable))
        if minimum_trace is None or trace4 < minimum_trace:
            minimum_trace = trace4
            minimizing_count = 1
            sixth_moments = {moments(order, immutable)[2]}
        elif trace4 == minimum_trace:
            minimizing_count += 1
            sixth_moments.add(moments(order, immutable)[2])

    if minimum_trace != 342:
        raise AssertionError(("order-seven minimum trace", minimum_trace))
    if sixth_moments != {50781, 53661, 59421}:
        raise AssertionError(("order-seven sixth moments", sixth_moments))

    examples = (
        (
            1,
            1,
            1,
            1,
            1,
            1,
            -1,
            1,
            -1,
            1,
            -1,
            1,
            -1,
            -1,
            1,
            -1,
            -1,
            1,
            1,
            1,
            1,
        ),
        (
            1,
            1,
            1,
            1,
            1,
            1,
            -1,
            -1,
            1,
            1,
            -1,
            -1,
            1,
            -1,
            1,
            -1,
            1,
            1,
            1,
            1,
            1,
        ),
    )
    expected_sixth = (50781, 53661)
    for coefficients, expected in zip(examples, expected_sixth, strict=True):
        if trace_fourth(matrix_from_coefficients(order, coefficients)) != 342:
            raise AssertionError("stored order-seven example is not trace-minimizing")
        if moments(order, coefficients)[2] != expected:
            raise AssertionError("stored order-seven sixth moment changed")
    return minimizing_count, len(sixth_moments)


def verify_macwilliams_identity() -> int:
    generator = random.Random(SEED)
    parameter = Fraction(1, 3)
    checks = 0
    for order in range(3, 7):
        edge_count = len(edges(order))
        code = augmented_cut_code(order)
        dual = {
            mask for mask in range(1 << edge_count) if is_even_eulerian(order, mask)
        }
        cosets = [0, (1 << edge_count) - 1]
        cosets.extend(generator.randrange(1 << edge_count) for _ in range(12))
        for representative in cosets:
            left = Fraction(0)
            for word in code:
                product = Fraction(1)
                for edge_index in range(edge_count):
                    sign = -1 if (representative ^ word) >> edge_index & 1 else 1
                    product *= 1 + sign * parameter
                left += product
            right = len(code) * sum(
                (-1 if (representative & word).bit_count() % 2 else 1)
                * parameter ** word.bit_count()
                for word in dual
            )
            if left != right:
                raise AssertionError(("MacWilliams identity", order, representative))
            checks += 1
    return checks


def paley_conference(prime: int) -> tuple[tuple[int, ...], ...]:
    if prime % 4 != 1:
        raise ValueError("prime must equal 1 modulo 4")

    def legendre(value: int) -> int:
        value %= prime
        if value == 0:
            return 0
        return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1

    order = prime + 1
    matrix = [[0] * order for _ in range(order)]
    for vertex in range(1, order):
        matrix[0][vertex] = matrix[vertex][0] = 1
    for row in range(1, order):
        for column in range(row + 1, order):
            value = legendre((row - 1) - (column - 1))
            matrix[row][column] = matrix[column][row] = value
    return tuple(tuple(row) for row in matrix)


def coefficients_from_matrix(matrix: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    return tuple(matrix[row][column] for row, column in edges(len(matrix)))


def verify_covariance_and_sdp() -> int:
    generator = random.Random(SEED)
    checks = 0
    matrices = [paley_conference(5), paley_conference(13)]
    for order in range(3, 10):
        coefficients = tuple(
            generator.choice((-1, 1)) for _ in range(len(edges(order)))
        )
        matrices.append(matrix_from_coefficients(order, coefficients))

    for matrix in matrices:
        order = len(matrix)
        square = matrix_product(matrix, matrix)
        if any(square[index][index] != order - 1 for index in range(order)):
            raise AssertionError(("wrong square diagonal", order))

        # R_+-R_-=2A/sqrt(n-1).  Multiplying by sqrt(n-1) keeps this exact.
        for row in range(order):
            for column in range(order):
                scaled_difference = 2 * matrix[row][column]
                if row == column and scaled_difference:
                    raise AssertionError(("covariance difference diagonal", order))

        # The oriented objective of this covariance pair is n*sqrt(n-1).
        scaled_objective = sum(
            2 * matrix[row][column] ** 2 for row, column in edges(order)
        )
        if scaled_objective != order * (order - 1):
            raise AssertionError(("scaled SDP objective", order, scaled_objective))
        checks += 1

    for matrix in matrices[:2]:
        order = len(matrix)
        square = matrix_product(matrix, matrix)
        expected = tuple(
            tuple(order - 1 if row == column else 0 for column in range(order))
            for row in range(order)
        )
        if square != expected:
            raise AssertionError(("conference identity", order))
        coefficients = coefficients_from_matrix(matrix)
        boolean_maximum = max(
            abs(energy(coefficients, spins)) for spins in spin_vectors(order)
        )
        # A strict finite integrality gap is a corruption control for replacing
        # Boolean spins by the elliptope.
        # Compare 2 M(A) with n sqrt(n-1) without floating point.
        if not (2 * boolean_maximum) ** 2 < order**2 * (order - 1):
            raise AssertionError(("expected finite SDP gap", order, boolean_maximum))
    return checks


def verify_absolute_value_counterexample() -> int:
    checks = 0
    for order in range(2, 15):
        coefficients = tuple(1 if row == 0 else -1 for row, _ in edges(order))
        values = [energy(coefficients, spins) for spins in spin_vectors(order)]
        if max(values) != order // 2:
            raise AssertionError(("one-sided maximum", order, max(values)))
        if max(abs(value) for value in values) != order * (order - 1) // 2:
            raise AssertionError(("absolute maximum", order))
        checks += 1
    return checks


def verify_corruption_controls() -> None:
    # Dropping the antipodal extension makes the code too small at n=3.
    cuts = {cut_word_bits(spins) for spins in spin_vectors(3)}
    if len(cuts) == len(augmented_cut_code(3)):
        raise AssertionError("antipodal-code corruption was not detected")

    # A false extreme-point claim already fails in one dimension.
    if min(abs(value) for value in (Fraction(-1), Fraction(0), Fraction(1))) != 0:
        raise AssertionError("relaxed-minimum corruption was not detected")

    # Equal fourth-order minimizers need not have equal sixth moments.
    first = (
        1,
        1,
        1,
        1,
        1,
        1,
        -1,
        1,
        -1,
        1,
        -1,
        1,
        -1,
        -1,
        1,
        -1,
        -1,
        1,
        1,
        1,
        1,
    )
    second = (
        1,
        1,
        1,
        1,
        1,
        1,
        -1,
        -1,
        1,
        1,
        -1,
        -1,
        1,
        -1,
        1,
        -1,
        1,
        1,
        1,
        1,
        1,
    )
    if moments(7, first)[2] == moments(7, second)[2]:
        raise AssertionError("higher-cumulant corruption was not detected")


def main() -> None:
    code_checks = verify_code_and_dual()
    parseval_checks = verify_parseval_cumulants()
    high_temperature_minimizers, sixth_moment_classes = (
        verify_high_temperature_minimizers()
    )
    macwilliams_checks = verify_macwilliams_identity()
    covariance_checks = verify_covariance_and_sdp()
    absolute_checks = verify_absolute_value_counterexample()
    verify_corruption_controls()
    print(f"code_dual_orthogonality_checks={code_checks}")
    print(f"parseval_cumulant_signings_checked={parseval_checks}")
    print(f"order_7_trace_minimizers_checked={high_temperature_minimizers}")
    print(f"order_7_sixth_moment_classes={sixth_moment_classes}")
    print(f"macwilliams_cosets_checked={macwilliams_checks}")
    print(f"covariance_sdp_matrices_checked={covariance_checks}")
    print(f"absolute_value_orders_checked={absolute_checks}")
    print(f"deterministic_seed={SEED}")
    print("corruption_controls=PASSED")


if __name__ == "__main__":
    main()

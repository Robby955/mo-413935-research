#!/usr/bin/env python3
"""Exact degree-19 moment certificate for the C14 labeled-shell vacancy.

This checker reconstructs the balanced 7+7 occupancy directly in gauge
space.  It does not use the Walsh inversion in
``verify_labeled_shell_parseval.py``.  The local graph and rectangular states
are bucketed by deficit, and every subthreshold product triple is sent to its
gauge by XOR.

The vacancy certificate uses the integer polynomial

    P(x) = (1-x) product_a ((x-a)(x-a-1))/(a(a+1)).

For integer x >= 1 every adjacent-root factor is nonnegative, so P(x) <= 0,
while P(0)=1.  Therefore E P(b)>0 forces an empty fiber.  All arithmetic in
the certificate is integer or Fraction arithmetic.
"""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict
from fractions import Fraction

if not __debug__:
    raise RuntimeError("verification requires Python assertions")


Matrix = tuple[tuple[int, ...], ...]
Spin = tuple[int, ...]

GROUP_SIZE = 1 << 13
TARGET_ENERGY = 27
ADJACENT_ROOTS = (9, 17, 26, 36, 46, 56, 67, 76, 86)
EXPECTED_CERTIFICATE = Fraction(
    1_707_454_816_960_049_615,
    99_244_391_564_512_637_853_696,
)
EXPECTED_LOCALIZING_NUMERATOR = 584_163_517_696_745_929_254_421_003_286_532
EXPECTED_CHEBYSHEV_ORDER = 14


def projective_spins(order: int) -> tuple[Spin, ...]:
    return tuple(
        (1,) + tail for tail in itertools.product((-1, 1), repeat=order - 1)
    )


def spin_mask(spin: Spin) -> int:
    return sum(1 << index for index, value in enumerate(spin[1:]) if value < 0)


def legendre_symbol(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def paley_conference(prime: int) -> Matrix:
    if prime % 4 != 1:
        raise ValueError("prime must be 1 modulo 4")
    order = prime + 1
    matrix = [[0] * order for _ in range(order)]
    for column in range(1, order):
        matrix[0][column] = matrix[column][0] = 1
    for row in range(1, order):
        for column in range(row + 1, order):
            sign = legendre_symbol(row - column, prime)
            matrix[row][column] = matrix[column][row] = sign
    return tuple(tuple(row) for row in matrix)


def assert_conference(matrix: Matrix) -> None:
    order = len(matrix)
    for row in range(order):
        for column in range(order):
            value = sum(matrix[row][index] * matrix[index][column] for index in range(order))
            expected = order - 1 if row == column else 0
            if value != expected:
                raise AssertionError(("conference identity", row, column, value))


def principal(matrix: Matrix, vertices: tuple[int, ...]) -> Matrix:
    return tuple(
        tuple(matrix[row][column] for column in vertices) for row in vertices
    )


def graph_energy(matrix: Matrix, spin: Spin) -> int:
    return sum(
        matrix[row][column] * spin[row] * spin[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def rectangular_energy(matrix: Matrix, left: Spin, right: Spin) -> int:
    return sum(
        matrix[row][column] * left[row] * right[column]
        for row in range(len(left))
        for column in range(len(right))
    )


def local_state_buckets() -> tuple[
    dict[int, list[int]],
    dict[int, list[int]],
    dict[int, list[int]],
    tuple[int, int, int],
]:
    conference = paley_conference(13)
    assert_conference(conference)
    left_vertices = tuple(range(7))
    right_vertices = tuple(range(7, 14))
    left = principal(conference, left_vertices)
    right = principal(conference, right_vertices)
    cross = tuple(
        tuple(conference[row][column] for column in right_vertices)
        for row in left_vertices
    )
    spins = projective_spins(7)

    left_maximum = max(abs(graph_energy(left, spin)) for spin in spins)
    right_maximum = max(abs(graph_energy(right, spin)) for spin in spins)
    cross_maximum = max(
        abs(rectangular_energy(cross, x, y)) for x in spins for y in spins
    )
    maxima = (left_maximum, right_maximum, cross_maximum)
    if maxima != (11, 11, 21):
        raise AssertionError(("local maxima", maxima))

    left_states: defaultdict[int, list[int]] = defaultdict(list)
    right_states: defaultdict[int, list[int]] = defaultdict(list)
    cross_states: defaultdict[int, list[int]] = defaultdict(list)

    for spin in spins:
        value = graph_energy(left, spin)
        for orientation in (-1, 1):
            deficit = left_maximum - orientation * value
            gauge = spin_mask(spin) | (int(orientation < 0) << 12)
            left_states[deficit].append(gauge)

    for spin in spins:
        value = graph_energy(right, spin)
        for orientation in (-1, 1):
            deficit = right_maximum - orientation * value
            gauge = (spin_mask(spin) << 6) | (int(orientation < 0) << 12)
            right_states[deficit].append(gauge)

    for x in spins:
        for y in spins:
            deficit = cross_maximum - abs(rectangular_energy(cross, x, y))
            gauge = spin_mask(x) | (spin_mask(y) << 6)
            cross_states[deficit].append(gauge)

    if sum(map(len, left_states.values())) != 128:
        raise AssertionError("left state count")
    if sum(map(len, right_states.values())) != 128:
        raise AssertionError("right state count")
    if sum(map(len, cross_states.values())) != 4096:
        raise AssertionError("cross state count")

    return dict(left_states), dict(right_states), dict(cross_states), maxima


def direct_occupancy(cutoff: int) -> tuple[list[int], int]:
    left_states, right_states, cross_states, _ = local_state_buckets()
    occupancy = [0] * GROUP_SIZE
    total = 0
    for left_deficit, left_gauges in left_states.items():
        for right_deficit, right_gauges in right_states.items():
            for cross_deficit, cross_gauges in cross_states.items():
                if left_deficit + right_deficit + cross_deficit > cutoff:
                    continue
                for left_gauge in left_gauges:
                    for right_gauge in right_gauges:
                        partial_gauge = left_gauge ^ right_gauge
                        for cross_gauge in cross_gauges:
                            occupancy[partial_gauge ^ cross_gauge] += 1
                            total += 1
    return occupancy, total


def multiply_polynomials(
    left: list[Fraction], right: list[Fraction]
) -> list[Fraction]:
    output = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            output[left_degree + right_degree] += left_value * right_value
    return output


def certificate_coefficients() -> list[Fraction]:
    coefficients = [Fraction(1), Fraction(-1)]
    for root in ADJACENT_ROOTS:
        denominator = root * (root + 1)
        factor = [
            Fraction(1),
            Fraction(-(2 * root + 1), denominator),
            Fraction(1, denominator),
        ]
        coefficients = multiply_polynomials(coefficients, factor)
    if len(coefficients) != 20:
        raise AssertionError(("certificate degree", len(coefficients) - 1))
    return coefficients


def evaluate_polynomial(coefficients: list[Fraction], value: int) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def localizing_polynomial(value: int) -> int:
    result = 1
    for root in ADJACENT_ROOTS:
        result *= 2 * value - (2 * root + 1)
    return result


def chebyshev_value(order: int, argument: Fraction) -> Fraction:
    if order == 0:
        return Fraction(1)
    previous = Fraction(1)
    current = argument
    for _ in range(1, order):
        previous, current = current, 2 * argument * current - previous
    return current


def verify() -> None:
    maxima = (11, 11, 21)
    cutoff = sum(maxima) - TARGET_ENERGY
    if cutoff != 16:
        raise AssertionError(("target cutoff", cutoff))
    occupancy, total = direct_occupancy(cutoff)
    histogram = Counter(occupancy)
    positive_values = [value for value in occupancy if value]
    if total != 304_908 or sum(occupancy) != total:
        raise AssertionError(("target triple count", total, sum(occupancy)))
    if histogram[0] != 1 or min(positive_values) != 6 or max(occupancy) != 87:
        raise AssertionError(
            ("occupancy extrema", histogram[0], min(positive_values), max(occupancy))
        )

    coefficients = certificate_coefficients()
    if evaluate_polynomial(coefficients, 0) != 1:
        raise AssertionError("P(0) != 1")
    if any(
        evaluate_polynomial(coefficients, value) > 0
        for value in range(1, max(occupancy) + 1)
    ):
        raise AssertionError("P is positive at a positive occupancy")

    direct_expectation = sum(
        evaluate_polynomial(coefficients, value) for value in occupancy
    ) / GROUP_SIZE
    moments = [
        Fraction(sum(value**degree for value in occupancy), GROUP_SIZE)
        for degree in range(len(coefficients))
    ]
    moment_expectation = sum(
        coefficient * moment
        for coefficient, moment in zip(coefficients, moments, strict=True)
    )
    if direct_expectation != moment_expectation:
        raise AssertionError(("moment reconstruction", direct_expectation, moment_expectation))
    if direct_expectation != EXPECTED_CERTIFICATE or direct_expectation <= 0:
        raise AssertionError(("certificate value", direct_expectation))

    certified_empty_fibers = (
        (GROUP_SIZE * direct_expectation).numerator
        + (GROUP_SIZE * direct_expectation).denominator
        - 1
    ) // (GROUP_SIZE * direct_expectation).denominator
    if certified_empty_fibers != 1:
        raise AssertionError(("certified empty fibers", certified_empty_fibers))

    localizing_numerator = sum(
        (1 - value) * localizing_polynomial(value) ** 2
        for value in occupancy
    )
    if localizing_numerator != EXPECTED_LOCALIZING_NUMERATOR:
        raise AssertionError(("localizing numerator", localizing_numerator))
    if localizing_numerator <= 0:
        raise AssertionError("localizing matrix certificate has wrong sign")

    minimum_positive = min(positive_values)
    maximum = max(occupancy)
    chebyshev_argument = Fraction(
        maximum + minimum_positive, maximum - minimum_positive
    )
    chebyshev_threshold = GROUP_SIZE * (maximum - 1)
    generic_order = next(
        order
        for order in range(1, 100)
        if chebyshev_value(order, chebyshev_argument) ** 2
        > chebyshev_threshold
    )
    if generic_order != EXPECTED_CHEBYSHEV_ORDER:
        raise AssertionError(("Chebyshev order", generic_order))
    if (
        chebyshev_value(generic_order - 1, chebyshev_argument) ** 2
        > chebyshev_threshold
    ):
        raise AssertionError("Chebyshev minimality")

    # Corruption control 1: replacing 1-x by 1+x destroys the required sign.
    leading_factors = [Fraction(1), Fraction(1)]
    remaining_factors = [Fraction(1)]
    for root in ADJACENT_ROOTS:
        denominator = root * (root + 1)
        remaining_factors = multiply_polynomials(
            remaining_factors,
            [
                Fraction(1),
                Fraction(-(2 * root + 1), denominator),
                Fraction(1, denominator),
            ],
        )
    wrong_leading = multiply_polynomials(leading_factors, remaining_factors)
    if not any(
        evaluate_polynomial(wrong_leading, value) > 0
        for value in range(1, max(occupancy) + 1)
    ):
        raise AssertionError("leading-sign corruption was not detected")

    # Corruption control 2: replacing adjacent roots 9,10 by 9,11 makes the
    # factor negative at x=10 and invalidates the integer sign certificate.
    corrupted = [Fraction(1), Fraction(-1)]
    corrupted = multiply_polynomials(
        corrupted,
        [Fraction(1), Fraction(-20, 99), Fraction(1, 99)],
    )
    for root in ADJACENT_ROOTS[1:]:
        denominator = root * (root + 1)
        corrupted = multiply_polynomials(
            corrupted,
            [
                Fraction(1),
                Fraction(-(2 * root + 1), denominator),
                Fraction(1, denominator),
            ],
        )
    if evaluate_polynomial(corrupted, 10) <= 0:
        raise AssertionError("nonadjacent-factor corruption was not detected")

    wrong_cutoff_total = sum(direct_occupancy(cutoff - 2)[0])
    if wrong_cutoff_total == total:
        raise AssertionError("target-cutoff corruption was not detected")

    roots_text = ",".join(map(str, ADJACENT_ROOTS))
    print(
        "c14_direct="
        f"maxima:{maxima[0]},{maxima[1]},{maxima[2]} "
        f"B:{total} K:{GROUP_SIZE} "
        f"occupancy:0x{histogram[0]},min_positive:{min(positive_values)},max:{max(occupancy)}"
    )
    print(
        f"degree19=adjacent_roots:{roots_text} "
        f"expectation:{direct_expectation} certified_empty_fibers:{certified_empty_fibers}"
    )
    print(
        "localizing_order:9 "
        f"positive_one_minus_b_numerator:{localizing_numerator}"
    )
    print(
        f"chebyshev_generic_order:{generic_order} "
        f"range:{minimum_positive}..{maximum}"
    )
    print("corruption_controls=leading_sign,nonconsecutive_factor,target_cutoff")
    print("arithmetic=integer,fraction")
    print("labeled_shell_moment_certificate=PASSED")


if __name__ == "__main__":
    verify()

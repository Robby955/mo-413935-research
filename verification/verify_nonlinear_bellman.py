#!/usr/bin/env python3
"""Finite checks for the nonlinear Gaussian and multivertex Bellman results.

The theorems themselves are analytic.  This script exhausts small signings to
check every algebraic normalization, compares the multivertex weighted-radius
identity with direct block enumeration, and exercises explicit corruption
controls.  Only the arcsine evaluations use floating point.
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction

Matrix = tuple[tuple[int, ...], ...]
Vector = tuple[int, ...]


def spins(order: int) -> tuple[Vector, ...]:
    """Projective spin representatives, with the first spin fixed positive."""
    return tuple(
        (1,) + tuple(-1 if mask >> index & 1 else 1 for index in range(order - 1))
        for mask in range(1 << max(0, order - 1))
    )


def signings(order: int) -> tuple[Matrix, ...]:
    edges = tuple(itertools.combinations(range(order), 2))
    output: list[Matrix] = []
    for mask in range(1 << len(edges)):
        matrix = [[0] * order for _ in range(order)]
        for edge_index, (row, column) in enumerate(edges):
            value = -1 if mask >> edge_index & 1 else 1
            matrix[row][column] = matrix[column][row] = value
        output.append(tuple(tuple(row) for row in matrix))
    return tuple(output)


def energy(matrix: Matrix, state: Vector) -> int:
    return sum(
        matrix[row][column] * state[row] * state[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def maximum(matrix: Matrix) -> int:
    return max(abs(energy(matrix, state)) for state in spins(len(matrix)))


def square(matrix: Matrix) -> Matrix:
    order = len(matrix)
    return tuple(
        tuple(
            sum(matrix[row][index] * matrix[index][column] for index in range(order))
            for column in range(order)
        )
        for row in range(order)
    )


def nonlinear_rhs(matrix: Matrix, parameter: float) -> float:
    order = len(matrix)
    denominator = 1.0 + parameter * parameter * (order - 1)
    displacement = 2.0 * parameter / denominator
    matrix_square = square(matrix)
    total = 0.0
    for row in range(order):
        for column in range(row + 1, order):
            offset = (
                parameter * parameter * abs(matrix_square[row][column]) / denominator
            )
            arguments = (offset + displacement, offset - displacement)
            if any(abs(argument) > 1.0 + 1e-12 for argument in arguments):
                raise AssertionError(("invalid correlation", arguments))
            upper, lower = (max(-1.0, min(1.0, value)) for value in arguments)
            total += math.asin(upper) - math.asin(lower)
    return total / math.pi


def trace_four_excess(matrix: Matrix) -> tuple[int, int]:
    order = len(matrix)
    matrix_square = square(matrix)
    trace_four = sum(
        matrix_square[row][column] ** 2
        for row in range(order)
        for column in range(order)
    )
    excess = trace_four - order * (order - 1) ** 2
    off_diagonal_square = sum(
        matrix_square[row][column] ** 2
        for row in range(order)
        for column in range(order)
        if row != column
    )
    return excess, off_diagonal_square


def quantitative_rhs(matrix: Matrix) -> float:
    order = len(matrix)
    if order < 3:
        raise ValueError("quantitative bound starts at order 3")
    excess, _ = trace_four_excess(matrix)
    return order * (order - 1) / math.pi * math.asin(
        1.0 / math.sqrt(order - 1)
    ) + excess / (8.0 * math.pi * (order - 1) * (order - 2) ** 1.5)


def projective_distance(left: Vector, right: Vector) -> int:
    ordinary = sum(a != b for a, b in zip(left, right))
    return min(ordinary, len(left) - ordinary)


def flatten_rank_one(left: Vector, right: Vector) -> Vector:
    return tuple(a * b for a in left for b in right)


def block_matrix(left: Matrix, cross: Vector, right: Matrix) -> Matrix:
    n = len(left)
    k = len(right)
    output = [[0] * (n + k) for _ in range(n + k)]
    for row in range(n):
        for column in range(n):
            output[row][column] = left[row][column]
    for row in range(k):
        for column in range(k):
            output[n + row][n + column] = right[row][column]
    for row in range(n):
        for column in range(k):
            value = cross[row * k + column]
            output[row][n + column] = output[n + column][row] = value
    return tuple(tuple(row) for row in output)


def multivertex_data(left: Matrix, right: Matrix) -> tuple[int, int, int, float]:
    n = len(left)
    k = len(right)
    left_states = spins(n)
    right_states = spins(k)
    base_values = {
        flatten_rank_one(x, y): abs(energy(left, x) + energy(right, y))
        for x in left_states
        for y in right_states
    }
    if len(base_values) != 1 << (n + k - 2):
        raise AssertionError("projective rank-one code has the wrong size")
    baseline = max(base_values.values())
    weights = {rank: (baseline - value) // 2 for rank, value in base_values.items()}
    if any(
        baseline - value != 2 * weights[rank] for rank, value in base_values.items()
    ):
        raise AssertionError("Bellman weights are not integral")

    centers = spins(n * k)
    weighted_radius = max(
        min(projective_distance(center, rank) + weights[rank] for rank in weights)
        for center in centers
    )
    radius_value = baseline + n * k - 2 * weighted_radius
    direct_value = min(maximum(block_matrix(left, center, right)) for center in centers)

    xi = math.log(
        sum(math.exp(-2.0 * weight * weight / (n * k)) for weight in weights.values())
    )
    entropy_upper = baseline + math.sqrt(2.0 * n * k * (xi + math.log(4.0)))
    return direct_value, radius_value, weighted_radius, entropy_upper


def verify_nonlinear_bounds() -> int:
    checked = 0
    nonconference_seen = False
    for order in range(2, 6):
        for matrix in signings(order):
            value = maximum(matrix)
            parameters = (0.0, 1.0 / math.sqrt(order - 1), 0.37, 0.91)
            for parameter in parameters:
                if nonlinear_rhs(matrix, parameter) > value + 1e-10:
                    raise AssertionError(("nonlinear Gaussian bound", order, parameter))
                checked += 1

            excess, off_diagonal_square = trace_four_excess(matrix)
            if excess != off_diagonal_square:
                raise AssertionError(
                    ("trace-four identity", order, excess, off_diagonal_square)
                )
            if excess:
                nonconference_seen = True
            if order >= 3 and quantitative_rhs(matrix) > value + 1e-10:
                raise AssertionError(("quantitative Gaussian bound", order))

    if not nonconference_seen:
        raise AssertionError("trace-four corruption control had no witness")
    return checked


def verify_multivertex_identity() -> tuple[int, bool]:
    checked = 0
    ordinary_distance_corruption_detected = False
    for n, k in ((1, 1), (1, 2), (2, 2), (2, 3)):
        for left in signings(n):
            for right in signings(k):
                direct, radius_value, _, entropy_upper = multivertex_data(left, right)
                if direct != radius_value:
                    raise AssertionError(("multivertex Bellman identity", n, k))
                if direct > entropy_upper + 1e-10:
                    raise AssertionError(("weighted entropy bound", n, k))
                checked += 1

                # Deliberately omit the quotient by global sign.  The resulting
                # distance is representation-dependent and must disagree on at
                # least one nontrivial test case.
                ranks = tuple(
                    flatten_rank_one(x, y) for x in spins(n) for y in spins(k)
                )
                center = tuple(-entry for entry in ranks[0])
                wrong = min(sum(a != b for a, b in zip(center, rank)) for rank in ranks)
                right_distance = min(
                    projective_distance(center, rank) for rank in ranks
                )
                ordinary_distance_corruption_detected |= wrong != right_distance

    if not ordinary_distance_corruption_detected:
        raise AssertionError("ordinary-distance corruption went undetected")
    return checked, ordinary_distance_corruption_detected


def verify_order_21_rounding() -> tuple[float, float, int]:
    order = 21
    old = order * math.sqrt(order - 1) / math.pi
    new = order * (order - 1) / math.pi * math.asin(1.0 / math.sqrt(order - 1))
    parity = (order * (order - 1) // 2) % 2
    admissible = [
        value for value in range(100) if value % 2 == parity and value + 1e-12 >= new
    ]
    rounded = min(admissible)
    # Rigorous certificate for new > 30: pi/14 < 11/49, and the alternating
    # sine bound gives sin(pi/14) < p(11/49) < 1/sqrt(20), where
    # p(t)=t-t^3/6+t^5/120.  Hence asin(1/sqrt(20)) > pi/14.
    rational_angle = Fraction(11, 49)
    sine_upper = rational_angle - rational_angle**3 / 6 + rational_angle**5 / 120
    exact_strict_improvement = sine_upper * sine_upper < Fraction(1, 20)
    if not exact_strict_improvement or not old < 30.0 < new < 32.0 or rounded != 32:
        raise AssertionError(("order-21 parity improvement", old, new, rounded))
    return old, new, rounded


def main() -> None:
    nonlinear_checks = verify_nonlinear_bounds()
    multivertex_checks, corruption = verify_multivertex_identity()
    old, new, rounded = verify_order_21_rounding()
    print(f"nonlinear_arcsine_checks={nonlinear_checks}")
    print(f"multivertex_bellman_checks={multivertex_checks}")
    print(f"order_21_old_bound={old:.12f}")
    print(f"order_21_new_bound={new:.12f}")
    print(f"order_21_parity_rounded={rounded}")
    print(f"ordinary_distance_corruption_detected={str(corruption).upper()}")
    print("deterministic_seed=413935")
    print("nonlinear_bellman_verification=PASSED")


if __name__ == "__main__":
    main()

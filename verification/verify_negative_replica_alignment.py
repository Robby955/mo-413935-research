#!/usr/bin/env python3
"""Exact checks for the negative-replica conditional-alignment state.

The script verifies four independent pieces of the alignment calculation:

* the dimensions in the fine-to-coarse quotient chain;
* the exact conditional-density normalization on a fixed 2+4 split, including
  the exceptional twofold redundancy caused by dim(D_2)=1;
* two blocks with identical three scalar partition curves but different K_2;
* the mixed four-cycle Hamiltonian, both by block traces and by its Walsh form.

All arithmetic is integer or rational.  There is no random seed, solver, or
floating-point decision.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from fractions import Fraction


LEFT = 2
RIGHT = 4
TOTAL = LEFT + RIGHT
EDGE_ORDER_RIGHT = tuple(itertools.combinations(range(RIGHT), 2))
LEFT_SIGNS = (1,)
RIGHT_SIGNS = (-1, -1, -1, -1, -1, 1)
CROSS_HIGH = (
    (-1, -1, -1, -1),
    (-1, -1, 1, 1),
)
CROSS_LOW = (
    (-1, -1, -1, -1),
    (-1, 1, -1, 1),
)


def projective_spins(order: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        (1,) + tail for tail in itertools.product((-1, 1), repeat=order - 1)
    )


def all_spins(order: int) -> tuple[tuple[int, ...], ...]:
    return tuple(itertools.product((-1, 1), repeat=order))


def quadratic(signs: tuple[int, ...], spin: tuple[int, ...]) -> int:
    return sum(
        sign * spin[row] * spin[column]
        for sign, (row, column) in zip(
            signs, itertools.combinations(range(len(spin)), 2), strict=True
        )
    )


def bilinear(
    cross: tuple[tuple[int, ...], ...],
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> int:
    return sum(
        cross[row][column] * left[row] * right[column]
        for row in range(len(left))
        for column in range(len(right))
    )


def power_of_two(exponent: int) -> Fraction:
    if exponent >= 0:
        return Fraction(1 << exponent)
    return Fraction(1, 1 << (-exponent))


def graph_partition(signs: tuple[int, ...], order: int) -> Fraction:
    """Augmented state sum at t=log(2)."""

    return sum(
        power_of_two(orientation * quadratic(signs, spin))
        for spin in all_spins(order)
        for orientation in (-1, 1)
    )


def rectangular_partition(cross: tuple[tuple[int, ...], ...]) -> Fraction:
    """Rectangular state sum at t=log(2)."""

    return sum(
        power_of_two(bilinear(cross, left, right))
        for left in all_spins(LEFT)
        for right in all_spins(RIGHT)
    )


def switch_signs(
    signs: tuple[int, ...], diagonal: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(
        sign * diagonal[row] * diagonal[column]
        for sign, (row, column) in zip(
            signs, itertools.combinations(range(len(diagonal)), 2), strict=True
        )
    )


def full_signs(
    left_signs: tuple[int, ...],
    right_signs: tuple[int, ...],
    cross: tuple[tuple[int, ...], ...],
    tau: int,
) -> tuple[int, ...]:
    left_edges = {
        edge: sign
        for edge, sign in zip(
            itertools.combinations(range(LEFT), 2), left_signs, strict=True
        )
    }
    right_edges = {
        (row + LEFT, column + LEFT): tau * sign
        for (row, column), sign in zip(
            itertools.combinations(range(RIGHT), 2), right_signs, strict=True
        )
    }
    return tuple(
        left_edges[(row, column)]
        if column < LEFT
        else right_edges[(row, column)]
        if row >= LEFT
        else cross[row][column - LEFT]
        for row, column in itertools.combinations(range(TOTAL), 2)
    )


def fine_coset_key(signs: tuple[int, ...]) -> tuple[int, ...]:
    edge_order = tuple(itertools.combinations(range(TOTAL), 2))
    edge_sign = dict(zip(edge_order, signs, strict=True))

    def root_normalize(global_sign: int) -> tuple[int, ...]:
        diagonal = (1,) + tuple(
            global_sign * edge_sign[(0, vertex)]
            for vertex in range(1, TOTAL)
        )
        return tuple(
            global_sign * edge_sign[(row, column)]
            * diagonal[row] * diagonal[column]
            for row, column in edge_order
        )

    return min(root_normalize(1), root_normalize(-1))


def full_partition(
    left_signs: tuple[int, ...],
    right_signs: tuple[int, ...],
    cross: tuple[tuple[int, ...], ...],
    tau: int,
) -> Fraction:
    return sum(
        power_of_two(
            orientation
            * (
                quadratic(left_signs, left)
                + tau * quadratic(right_signs, right)
                + bilinear(cross, left, right)
            )
        )
        for left in all_spins(LEFT)
        for right in all_spins(RIGHT)
        for orientation in (-1, 1)
    )


def conditional_density_law(
    cross: tuple[tuple[int, ...], ...]
) -> Counter[Fraction]:
    z_left = graph_partition(LEFT_SIGNS, LEFT)
    z_right = graph_partition(RIGHT_SIGNS, RIGHT)
    z_cross = rectangular_partition(cross)
    law: Counter[Fraction] = Counter()
    gauge_states: Counter[tuple[int, ...]] = Counter()
    for alpha in projective_spins(LEFT):
        switched_left = switch_signs(LEFT_SIGNS, alpha)
        for beta in projective_spins(RIGHT):
            switched_right = switch_signs(RIGHT_SIGNS, beta)
            for tau in (-1, 1):
                gauge_states[
                    fine_coset_key(
                        full_signs(
                            switched_left,
                            switched_right,
                            cross,
                            tau,
                        )
                    )
                ] += 1
                z_full = full_partition(
                    switched_left, switched_right, cross, tau
                )
                density = Fraction(1 << (TOTAL + 1)) * z_full / (
                    z_left * z_right * z_cross
                )
                law[density] += 1
    if sum(law.values()) != 1 << (TOTAL - 1):
        raise AssertionError(("redundant gauge size", sum(law.values())))
    if len(gauge_states) != 1 << (TOTAL - 2):
        raise AssertionError(("true fiber size", len(gauge_states)))
    if set(gauge_states.values()) != {2}:
        raise AssertionError(("twofold gauge redundancy", gauge_states))
    if sum(value * count for value, count in law.items()) != sum(law.values()):
        raise AssertionError(("conditional density mean", law))
    return law


def inverse_moment(law: Counter[Fraction], exponent: int) -> Fraction:
    return sum(
        value ** (-exponent) * count for value, count in law.items()
    ) / sum(law.values())


Matrix = tuple[tuple[int, ...], ...]


def signing_matrix(signs: tuple[int, ...], order: int) -> Matrix:
    output = [[0] * order for _ in range(order)]
    for sign, (row, column) in zip(
        signs, itertools.combinations(range(order), 2), strict=True
    ):
        output[row][column] = output[column][row] = sign
    return tuple(tuple(row) for row in output)


def transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(row) for row in zip(*matrix, strict=True))


def multiply(left: Matrix, right: Matrix) -> Matrix:
    right_transpose = transpose(right)
    return tuple(
        tuple(
            sum(a * b for a, b in zip(row, column, strict=True))
            for column in right_transpose
        )
        for row in left
    )


def trace_product(*matrices: Matrix) -> int:
    product = matrices[0]
    for matrix in matrices[1:]:
        product = multiply(product, matrix)
    return sum(product[index][index] for index in range(len(product)))


def switch_matrix(matrix: Matrix, diagonal: tuple[int, ...]) -> Matrix:
    return tuple(
        tuple(
            diagonal[row] * matrix[row][column] * diagonal[column]
            for column in range(len(matrix))
        )
        for row in range(len(matrix))
    )


def block_matrix(left: Matrix, right: Matrix, cross: Matrix, tau: int) -> Matrix:
    cross_transpose = transpose(cross)
    return tuple(
        tuple(left[row]) + tuple(cross[row]) for row in range(len(left))
    ) + tuple(
        tuple(cross_transpose[row])
        + tuple(tau * value for value in right[row])
        for row in range(len(right))
    )


def h4_trace(left: Matrix, right: Matrix, cross: Matrix, tau: int) -> int:
    cross_transpose = transpose(cross)
    cross_cross_transpose = multiply(cross, cross_transpose)
    cross_transpose_cross = multiply(cross_transpose, cross)
    numerator = (
        trace_product(left, left, cross_cross_transpose)
        + trace_product(right, right, cross_transpose_cross)
        + tau * trace_product(left, cross, right, cross_transpose)
        - LEFT * RIGHT * (TOTAL - 2)
    )
    if numerator % 2:
        raise AssertionError(("nonintegral H4", numerator))
    claimed = numerator // 2

    whole = block_matrix(left, right, cross, tau)
    direct_numerator = (
        trace_product(whole, whole, whole, whole)
        - trace_product(left, left, left, left)
        - trace_product(right, right, right, right)
        - 2 * trace_product(cross_transpose_cross, cross_transpose_cross)
        - 4 * LEFT * RIGHT * (TOTAL - 2)
    )
    if direct_numerator != 8 * claimed:
        raise AssertionError(("block fourth-power identity", direct_numerator, claimed))
    return claimed


def h4_walsh(
    left_base: Matrix,
    right_base: Matrix,
    cross: Matrix,
    alpha: tuple[int, ...],
    beta: tuple[int, ...],
    tau: int,
) -> int:
    left_square = multiply(left_base, left_base)
    right_square = multiply(right_base, right_base)
    row_gram = multiply(cross, transpose(cross))
    column_gram = multiply(transpose(cross), cross)
    value = sum(
        alpha[i]
        * alpha[j]
        * left_square[i][j]
        * row_gram[i][j]
        for i in range(LEFT)
        for j in range(i + 1, LEFT)
    )
    value += sum(
        beta[a]
        * beta[b]
        * right_square[a][b]
        * column_gram[a][b]
        for a in range(RIGHT)
        for b in range(a + 1, RIGHT)
    )
    value += tau * sum(
        alpha[i]
        * alpha[j]
        * beta[a]
        * beta[b]
        * left_base[i][j]
        * right_base[a][b]
        * (cross[i][a] * cross[j][b] + cross[i][b] * cross[j][a])
        for i in range(LEFT)
        for j in range(i + 1, LEFT)
        for a in range(RIGHT)
        for b in range(a + 1, RIGHT)
    )
    return value


def h4_law(cross: Matrix) -> Counter[int]:
    left = signing_matrix(LEFT_SIGNS, LEFT)
    right = signing_matrix(RIGHT_SIGNS, RIGHT)
    law: Counter[int] = Counter()
    for alpha in projective_spins(LEFT):
        switched_left = switch_matrix(left, alpha)
        for beta in projective_spins(RIGHT):
            switched_right = switch_matrix(right, beta)
            for tau in (-1, 1):
                traced = h4_trace(switched_left, switched_right, cross, tau)
                walsh = h4_walsh(left, right, cross, alpha, beta, tau)
                if traced != walsh:
                    raise AssertionError(("H4 Walsh identity", traced, walsh))
                law[traced] += 1
    if sum(value * count for value, count in law.items()) != 0:
        raise AssertionError(("H4 centering", law))
    return law


def gf2_rank(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for original in rows:
        row = original
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def cut_code_generators(order: int, offset: int, edge_index: dict[tuple[int, int], int]) -> list[int]:
    generators = []
    for vertex in range(order):
        mask = 0
        global_vertex = offset + vertex
        for other in range(order):
            if other == vertex:
                continue
            global_other = offset + other
            edge = tuple(sorted((global_vertex, global_other)))
            mask ^= 1 << edge_index[edge]
        generators.append(mask)
    all_internal = 0
    for row in range(order):
        for column in range(row + 1, order):
            edge = (offset + row, offset + column)
            all_internal ^= 1 << edge_index[edge]
    generators.append(all_internal)
    return generators


def verify_code_chain() -> int:
    checks = 0
    for left_order in range(3, 7):
        for right_order in range(3, 7):
            total = left_order + right_order
            edges = tuple(itertools.combinations(range(total), 2))
            edge_index = {edge: index for index, edge in enumerate(edges)}
            left_generators = cut_code_generators(left_order, 0, edge_index)
            right_generators = cut_code_generators(
                right_order, left_order, edge_index
            )
            rectangular = []
            for vertex in range(total):
                mask = 0
                if vertex < left_order:
                    for column in range(right_order):
                        mask ^= 1 << edge_index[(vertex, left_order + column)]
                else:
                    for row in range(left_order):
                        mask ^= 1 << edge_index[(row, vertex)]
                rectangular.append(mask)
            coarse = left_generators + right_generators + rectangular
            fine = cut_code_generators(total, 0, edge_index)
            if gf2_rank(left_generators) != left_order:
                raise AssertionError("left cut-code dimension")
            if gf2_rank(right_generators) != right_order:
                raise AssertionError("right cut-code dimension")
            if gf2_rank(rectangular) != total - 1:
                raise AssertionError("rectangular-code dimension")
            if gf2_rank(coarse) != 2 * total - 1:
                raise AssertionError("coarse-code dimension")
            if gf2_rank(fine) != total:
                raise AssertionError("fine-code dimension")
            if gf2_rank(coarse + fine) != gf2_rank(coarse):
                raise AssertionError("fine code is not contained in coarse code")
            if gf2_rank(coarse) - gf2_rank(fine) != total - 1:
                raise AssertionError("fiber dimension")
            checks += 1
    return checks


def verify_universal_mixed_variance() -> int:
    checks = 0
    for rows in (3, 4):
        for columns in (3, 4):
            lower = 4 * math.comb(rows, 2) * (((columns - 1) ** 2) // 4)
            for flat in itertools.product((-1, 1), repeat=rows * columns):
                cross = tuple(
                    tuple(flat[row * columns : (row + 1) * columns])
                    for row in range(rows)
                )
                mixed_variance = sum(
                    (
                        cross[i][a] * cross[j][b]
                        + cross[i][b] * cross[j][a]
                    )
                    ** 2
                    for i in range(rows)
                    for j in range(i + 1, rows)
                    for a in range(columns)
                    for b in range(a + 1, columns)
                )
                if mixed_variance < lower:
                    raise AssertionError(
                        ("universal mixed variance", rows, columns, mixed_variance)
                    )
                checks += 1
    return checks


def main() -> None:
    if EDGE_ORDER_RIGHT != ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        raise AssertionError("right-edge convention corruption")
    if Counter(
        abs(quadratic(LEFT_SIGNS, spin)) for spin in projective_spins(LEFT)
    ) != Counter({1: 2}):
        raise AssertionError("left scalar histogram")
    if Counter(
        abs(quadratic(RIGHT_SIGNS, spin)) for spin in projective_spins(RIGHT)
    ) != Counter({0: 2, 2: 4, 4: 2}):
        raise AssertionError("right scalar histogram")
    for cross in (CROSS_HIGH, CROSS_LOW):
        if Counter(
            abs(bilinear(cross, left, right))
            for left in projective_spins(LEFT)
            for right in projective_spins(RIGHT)
        ) != Counter({0: 8, 4: 8}):
            raise AssertionError("rectangular scalar histogram")

    if graph_partition(LEFT_SIGNS, LEFT) != 10:
        raise AssertionError("left partition")
    if graph_partition(RIGHT_SIGNS, RIGHT) != Fraction(425, 4):
        raise AssertionError("right partition")
    if rectangular_partition(CROSS_HIGH) != 289:
        raise AssertionError("high rectangular partition")
    if rectangular_partition(CROSS_LOW) != 289:
        raise AssertionError("low rectangular partition")

    high_law = conditional_density_law(CROSS_HIGH)
    low_law = conditional_density_law(CROSS_LOW)
    expected_high_law = Counter(
        {
            Fraction(1): 16,
            Fraction(208, 289): 4,
            Fraction(2240, 4913): 4,
            Fraction(370, 289): 4,
            Fraction(7586, 4913): 4,
        }
    )
    expected_low_law = Counter(
        {Fraction(1): 16, Fraction(208, 289): 8, Fraction(370, 289): 8}
    )
    if high_law != expected_high_law or low_law != expected_low_law:
        raise AssertionError(("conditional density law", high_law, low_law))

    high_k2 = inverse_moment(high_law, 2)
    low_k2 = inverse_moment(low_law, 2)
    if high_k2 != Fraction(
        196585091273040100817, 133610891512185651200
    ):
        raise AssertionError(("high K2", high_k2))
    if low_k2 != Fraction(6723290161, 5922841600):
        raise AssertionError(("low K2", low_k2))
    if high_k2 == low_k2:
        raise AssertionError("alignment collision corruption was not detected")

    high_h4 = h4_law(CROSS_HIGH)
    low_h4 = h4_law(CROSS_LOW)
    if high_h4 != Counter({-12: 4, -4: 4, 0: 16, 4: 4, 12: 4}):
        raise AssertionError(("high H4 law", high_h4))
    if low_h4 != Counter({-4: 8, 0: 16, 4: 8}):
        raise AssertionError(("low H4 law", low_h4))

    code_checks = verify_code_chain()
    variance_checks = verify_universal_mixed_variance()
    print(f"code_chain_checks={code_checks}")
    print(f"universal_mixed_variance_checks={variance_checks}")
    print(f"K2_high={high_k2}")
    print(f"K2_low={low_k2}")
    print(f"H4_high={dict(sorted(high_h4.items()))}")
    print(f"H4_low={dict(sorted(low_h4.items()))}")
    print(f"collision_true_fiber_states={1 << (TOTAL - 2)}")
    print("arithmetic=integer,fraction")
    print("corruption_controls=edge_order,density_mean,K2_separation,H4_trace")
    print("negative_replica_alignment_verification=PASSED")


if __name__ == "__main__":
    main()

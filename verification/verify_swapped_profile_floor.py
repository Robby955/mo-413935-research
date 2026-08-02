#!/usr/bin/env python3
"""Exact checks for the swapped-profile scalar floor.

The swapped profile uses projective absolute graph states and signed full
spin pairs for the rectangular block.  For a balanced ``r+r`` split, the
proof checked here has two parts.

1.  If ``mu_r = E|eps_1 + ... + eps_r|``, every ``r x r`` sign matrix has
    at least ``2**r`` signed full-spin pairs with rectangular energy at
    least ``r*mu_r - 3*r``.
2.  Consequently, the ``2**(2*r-1)``-st largest swapped product energy is
    at least ``max(M(A), M(B)) + r*mu_r - 3*r``.

All arithmetic used in the asserted inequalities is integer or Fraction
arithmetic.  The exhaustive rectangular check covers every sign matrix
through order four.  The raw product-floor check covers every
switching-normalized balanced block triple through order three.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction
from itertools import product
from math import comb

if not __debug__:
    raise RuntimeError("verification requires Python assertions")


Spin = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


EXPECTED_RECTANGULAR_MATRICES = 66_066
EXPECTED_RECTANGULAR_PAIRS = 16_810_248
EXPECTED_RAW_CASES = 66
EXPECTED_RAW_PRODUCT_STATES = 65_664
EXPECTED_STREAM_SHA256 = (
    "9552e95262c5cbc683e64dfea427da2952cc3b686442af3146385f67c7567b32"
)


def full_spins(order: int) -> tuple[Spin, ...]:
    return tuple(product((-1, 1), repeat=order))


def projective_spins(order: int) -> tuple[Spin, ...]:
    return tuple((1, *tail) for tail in product((-1, 1), repeat=order - 1))


def sign_matrices(order: int) -> tuple[Matrix, ...]:
    return tuple(
        tuple(
            tuple(signs[row * order + column] for column in range(order))
            for row in range(order)
        )
        for signs in product((-1, 1), repeat=order * order)
    )


def switching_normalized_signings(order: int) -> tuple[Matrix, ...]:
    residual_edges = tuple(
        (left, right)
        for left in range(1, order)
        for right in range(left + 1, order)
    )
    result = []
    for signs in product((-1, 1), repeat=len(residual_edges)):
        matrix = [[0] * order for _ in range(order)]
        for vertex in range(1, order):
            matrix[0][vertex] = matrix[vertex][0] = 1
        for (left, right), sign in zip(residual_edges, signs):
            matrix[left][right] = matrix[right][left] = sign
        result.append(tuple(tuple(row) for row in matrix))
    return tuple(result)


def switching_normalized_rectangles(order: int) -> tuple[Matrix, ...]:
    residual_entries = tuple(
        (row, column)
        for row in range(1, order)
        for column in range(1, order)
    )
    result = []
    for signs in product((-1, 1), repeat=len(residual_entries)):
        matrix = [[1] * order for _ in range(order)]
        for (row, column), sign in zip(residual_entries, signs):
            matrix[row][column] = sign
        result.append(tuple(tuple(row) for row in matrix))
    return tuple(result)


def rademacher_absolute_mean(order: int) -> Fraction:
    return Fraction(
        sum(comb(order, minus) * abs(order - 2 * minus) for minus in range(order + 1)),
        1 << order,
    )


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


def row_fields(matrix: Matrix, right: Spin) -> tuple[int, ...]:
    return tuple(
        sum(entry * sign for entry, sign in zip(row, right)) for row in matrix
    )


def maximizing_left_spin(fields: tuple[int, ...]) -> Spin:
    return tuple(1 if field >= 0 else -1 for field in fields)


def flip(spin: Spin, coordinate: int) -> Spin:
    return tuple(-entry if index == coordinate else entry for index, entry in enumerate(spin))


def cross_multiplicity_witnesses(matrix: Matrix) -> tuple[set[tuple[Spin, Spin]], Fraction, int]:
    order = len(matrix)
    mu = rademacher_absolute_mean(order)
    good_cutoff = order * mu - order
    energy_cutoff = order * mu - 3 * order
    witnesses: set[tuple[Spin, Spin]] = set()
    good_count = 0
    absolute_sum = 0
    for right in full_spins(order):
        fields = row_fields(matrix, right)
        row_maximum = sum(abs(field) for field in fields)
        absolute_sum += row_maximum
        if row_maximum < good_cutoff:
            continue
        good_count += 1
        maximizing = maximizing_left_spin(fields)
        candidates = (maximizing,) + tuple(
            flip(maximizing, coordinate) for coordinate in range(order)
        )
        assert len(set(candidates)) == order + 1
        for left in candidates:
            assert rectangular_energy(matrix, left, right) >= energy_cutoff
            witnesses.add((left, right))

    assert Fraction(absolute_sum, 1 << order) == order * mu
    assert (order + 1) * good_count >= 1 << order
    assert len(witnesses) == (order + 1) * good_count
    assert len(witnesses) >= 1 << order
    return witnesses, energy_cutoff, good_count


def verify_rectangular_floor() -> tuple[int, int, str]:
    matrix_count = 0
    pair_count = 0
    digest = hashlib.sha256()
    omitted_flip_failure_detected = False
    for order in range(1, 5):
        rank = 1 << order
        for matrix_index, matrix in enumerate(sign_matrices(order)):
            witnesses, cutoff, good_count = cross_multiplicity_witnesses(matrix)
            energies = sorted(
                (
                    rectangular_energy(matrix, left, right)
                    for left in full_spins(order)
                    for right in full_spins(order)
                ),
                reverse=True,
            )
            assert Fraction(energies[rank - 1]) >= cutoff
            omitted_flip_failure_detected |= good_count < rank
            matrix_count += 1
            pair_count += len(energies)
            digest.update(
                (
                    f"cross,{order},{matrix_index},{cutoff},"
                    f"{good_count},{len(witnesses)},{energies[rank - 1]}\n"
                ).encode("ascii")
            )

    assert omitted_flip_failure_detected
    assert matrix_count == EXPECTED_RECTANGULAR_MATRICES
    assert pair_count == EXPECTED_RECTANGULAR_PAIRS
    return matrix_count, pair_count, digest.hexdigest()


def raw_swapped_energy_floor(
    graph_left: Matrix,
    graph_right: Matrix,
    rectangle: Matrix,
) -> tuple[Fraction, int, int]:
    order = len(graph_left)
    assert len(graph_right) == len(rectangle) == order
    graph_left_energies = tuple(
        abs(quadratic_energy(graph_left, spin)) for spin in projective_spins(order)
    )
    graph_right_energies = tuple(
        abs(quadratic_energy(graph_right, spin)) for spin in projective_spins(order)
    )
    cross_energies = tuple(
        rectangular_energy(rectangle, left, right)
        for left in full_spins(order)
        for right in full_spins(order)
    )
    product_energies = sorted(
        (
            left_energy + right_energy + cross_energy
            for left_energy in graph_left_energies
            for right_energy in graph_right_energies
            for cross_energy in cross_energies
        ),
        reverse=True,
    )
    rank = 1 << (2 * order - 1)
    assert len(product_energies) == rank * rank
    rank_energy = product_energies[rank - 1]
    floor = (
        max(max(graph_left_energies), max(graph_right_energies))
        + order * rademacher_absolute_mean(order)
        - 3 * order
    )
    assert Fraction(rank_energy) >= floor
    return floor, rank_energy, len(product_energies)


def verify_raw_product_floor() -> tuple[int, int, str]:
    case_count = 0
    state_count = 0
    digest = hashlib.sha256()
    incompatible_double_maximum_rejected = False
    for order in (2, 3):
        graphs = switching_normalized_signings(order)
        rectangles = switching_normalized_rectangles(order)
        for left_index, graph_left in enumerate(graphs):
            for right_index, graph_right in enumerate(graphs):
                for rectangle_index, rectangle in enumerate(rectangles):
                    floor, rank_energy, product_states = raw_swapped_energy_floor(
                        graph_left, graph_right, rectangle
                    )
                    left_maximum = max(
                        abs(quadratic_energy(graph_left, spin))
                        for spin in projective_spins(order)
                    )
                    right_maximum = max(
                        abs(quadratic_energy(graph_right, spin))
                        for spin in projective_spins(order)
                    )
                    corrupted_floor = (
                        left_maximum
                        + right_maximum
                        + order * rademacher_absolute_mean(order)
                    )
                    incompatible_double_maximum_rejected |= (
                        Fraction(rank_energy) < corrupted_floor
                    )
                    case_count += 1
                    state_count += product_states
                    digest.update(
                        (
                            f"raw,{order},{left_index},{right_index},"
                            f"{rectangle_index},{floor},{rank_energy},{product_states}\n"
                        ).encode("ascii")
                    )

    assert incompatible_double_maximum_rejected
    assert case_count == EXPECTED_RAW_CASES
    assert state_count == EXPECTED_RAW_PRODUCT_STATES
    return case_count, state_count, digest.hexdigest()


def main() -> None:
    rectangular_matrices, rectangular_pairs, cross_digest = verify_rectangular_floor()
    raw_cases, raw_states, raw_digest = verify_raw_product_floor()
    stream_digest = hashlib.sha256(f"{cross_digest}\n{raw_digest}\n".encode("ascii")).hexdigest()
    assert stream_digest == EXPECTED_STREAM_SHA256
    print(f"rectangular_matrices_checked={rectangular_matrices}")
    print(f"rectangular_pairs_checked={rectangular_pairs}")
    print(f"raw_balanced_cases_checked={raw_cases}")
    print(f"raw_product_states_checked={raw_states}")
    print(f"stream_sha256={stream_digest}")
    print("corruption_controls=radius_one_multiplicity,incompatible_double_maximum")
    print("swapped_profile_floor=PASSED")


if __name__ == "__main__":
    main()

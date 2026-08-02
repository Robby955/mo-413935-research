#!/usr/bin/env python3
r"""Exact finite checks for the cavity and hereditary inequalities.

For a signing ``A`` and a vertex partition ``S \sqcup T``, write

    Q_A(x, y) = q(x) + h(x) . y + r(y).

The script exhausts every signing, every vertex partition, and every fixed
spin vector on ``S`` through order five.  It checks the exact Walsh
second-moment decomposition

    E_y Q_A(x, y)^2 = q(x)^2 + ||h(x)||_2^2 + binom(|T|, 2)

and the exact pairing identity obtained from ``y`` and ``-y``.  It also
checks the parity and integer-lattice consequences, independently recomputes
the small values of ``F``, records the order-four failure of exact
subadditivity for ``H(n) = F(n)^(2/3)``, and exercises deliberate corruption
controls.

All pass/fail decisions use integer or Fraction arithmetic.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Iterator
from fractions import Fraction


MAX_ORDER = 5
EXPECTED_F = {1: 0, 2: 1, 3: 3, 4: 4, 5: 4}


def edges(order: int) -> tuple[tuple[int, int], ...]:
    return tuple(itertools.combinations(range(order), 2))


def sign_vectors(length: int) -> Iterator[tuple[int, ...]]:
    return itertools.product((-1, 1), repeat=length)


def energy(
    edge_list: tuple[tuple[int, int], ...],
    coefficients: tuple[int, ...],
    spins: tuple[int, ...],
) -> int:
    return sum(
        coefficient * spins[row] * spins[column]
        for coefficient, (row, column) in zip(
            coefficients, edge_list, strict=True
        )
    )


def spins_on_partition(
    order: int,
    subset: tuple[int, ...],
    complement: tuple[int, ...],
    subset_spins: tuple[int, ...],
    complement_spins: tuple[int, ...],
) -> tuple[int, ...]:
    spins = [0] * order
    for vertex, value in zip(subset, subset_spins, strict=True):
        spins[vertex] = value
    for vertex, value in zip(complement, complement_spins, strict=True):
        spins[vertex] = value
    if any(value == 0 for value in spins):
        raise AssertionError("partition did not assign every spin")
    return tuple(spins)


def block_data(
    edge_list: tuple[tuple[int, int], ...],
    coefficients: tuple[int, ...],
    subset: tuple[int, ...],
    complement: tuple[int, ...],
    subset_spins: tuple[int, ...],
) -> tuple[int, tuple[int, ...], tuple[tuple[int, int, int], ...]]:
    """Return q, the cross field h, and indexed internal T edges."""
    subset_positions = {vertex: index for index, vertex in enumerate(subset)}
    complement_positions = {
        vertex: index for index, vertex in enumerate(complement)
    }
    q = 0
    field = [0] * len(complement)
    internal_complement: list[tuple[int, int, int]] = []

    for coefficient, (row, column) in zip(
        coefficients, edge_list, strict=True
    ):
        row_in_subset = row in subset_positions
        column_in_subset = column in subset_positions
        if row_in_subset and column_in_subset:
            q += (
                coefficient
                * subset_spins[subset_positions[row]]
                * subset_spins[subset_positions[column]]
            )
        elif row_in_subset:
            field[complement_positions[column]] += (
                coefficient * subset_spins[subset_positions[row]]
            )
        elif column_in_subset:
            field[complement_positions[row]] += (
                coefficient * subset_spins[subset_positions[column]]
            )
        else:
            internal_complement.append(
                (
                    complement_positions[row],
                    complement_positions[column],
                    coefficient,
                )
            )

    return q, tuple(field), tuple(internal_complement)


def complement_energy(
    internal_edges: tuple[tuple[int, int, int], ...],
    spins: tuple[int, ...],
) -> int:
    return sum(
        coefficient * spins[row] * spins[column]
        for row, column, coefficient in internal_edges
    )


def mu_numerator(length: int) -> int:
    """Return 2^length times E|epsilon_1 + ... + epsilon_length|."""
    return sum(abs(sum(spins)) for spins in sign_vectors(length))


def mu(length: int) -> Fraction:
    return Fraction(mu_numerator(length), 1 << length)


def verify_mu_formula() -> int:
    checks = 0
    for length in range(MAX_ORDER + 1):
        observed = mu(length)
        if length == 0:
            expected = Fraction(0)
        else:
            expected = Fraction(
                length
                * math.comb(length - 1, (length - 1) // 2),
                1 << (length - 1),
            )
        if observed != expected:
            raise AssertionError(("mu formula", length, observed, expected))
        checks += 1
    return checks


def verify_cavity_inequalities() -> tuple[dict[int, int], dict[str, int]]:
    minimum_maximum: dict[int, int] = {}
    counts = {
        "signings": 0,
        "partitions": 0,
        "fixed_x": 0,
        "paired_y": 0,
    }

    for order in range(1, MAX_ORDER + 1):
        edge_list = edges(order)
        order_minimum: int | None = None
        all_full_spins = tuple(sign_vectors(order))

        for coefficients in sign_vectors(len(edge_list)):
            counts["signings"] += 1
            maximum = max(
                abs(energy(edge_list, coefficients, spins))
                for spins in all_full_spins
            )
            if order_minimum is None or maximum < order_minimum:
                order_minimum = maximum

            for subset_mask in range(1 << order):
                subset = tuple(
                    vertex
                    for vertex in range(order)
                    if subset_mask >> vertex & 1
                )
                complement = tuple(
                    vertex
                    for vertex in range(order)
                    if not (subset_mask >> vertex & 1)
                )
                subset_size = len(subset)
                complement_size = len(complement)
                complement_spins = tuple(sign_vectors(complement_size))
                counts["partitions"] += 1

                fixed_x_records: list[tuple[int, int]] = []
                for subset_spins in sign_vectors(subset_size):
                    counts["fixed_x"] += 1
                    q, field, internal_edges = block_data(
                        edge_list,
                        coefficients,
                        subset,
                        complement,
                        subset_spins,
                    )
                    field_square = sum(value * value for value in field)
                    direct_second_sum = 0
                    base_absolute_sum = 0
                    field_absolute_sum = 0

                    for y in complement_spins:
                        counts["paired_y"] += 1
                        r = complement_energy(internal_edges, y)
                        linear = sum(
                            coefficient * value
                            for coefficient, value in zip(field, y, strict=True)
                        )
                        full_spins = spins_on_partition(
                            order, subset, complement, subset_spins, y
                        )
                        negative_y_spins = spins_on_partition(
                            order,
                            subset,
                            complement,
                            subset_spins,
                            tuple(-value for value in y),
                        )
                        plus_energy = energy(
                            edge_list, coefficients, full_spins
                        )
                        minus_energy = energy(
                            edge_list, coefficients, negative_y_spins
                        )
                        if plus_energy != q + r + linear:
                            raise AssertionError(("block expansion", order))
                        if minus_energy != q + r - linear:
                            raise AssertionError(("global T flip", order))

                        paired_maximum = max(
                            abs(plus_energy), abs(minus_energy)
                        )
                        paired_sum = abs(q + r) + abs(linear)
                        if paired_maximum != paired_sum:
                            raise AssertionError(("pairing identity", order))
                        if paired_maximum > maximum:
                            raise AssertionError(("maximum domination", order))

                        direct_second_sum += plus_energy * plus_energy
                        base_absolute_sum += abs(q + r)
                        field_absolute_sum += abs(linear)

                    state_count = 1 << complement_size
                    expected_second_sum = state_count * (
                        q * q
                        + field_square
                        + math.comb(complement_size, 2)
                    )
                    if direct_second_sum != expected_second_sum:
                        raise AssertionError(
                            (
                                "Walsh second moment",
                                order,
                                subset_mask,
                                direct_second_sum,
                                expected_second_sum,
                            )
                        )
                    if expected_second_sum > state_count * maximum * maximum:
                        raise AssertionError(("exact square inequality", order))

                    expected_parity = subset_size % 2
                    if any(abs(value) % 2 != expected_parity for value in field):
                        raise AssertionError(("cross-field parity", order))
                    parity_floor = complement_size * expected_parity
                    if field_square < parity_floor:
                        raise AssertionError(("cross-field square floor", order))

                    if base_absolute_sum < state_count * abs(q):
                        raise AssertionError(("Jensen step", order))
                    if (
                        state_count * maximum
                        < state_count * abs(q) + field_absolute_sum
                    ):
                        raise AssertionError(("averaged pairing inequality", order))

                    # The p=1 sharp Khintchine bound, squared to stay exact.
                    if (
                        2 * field_absolute_sum * field_absolute_sum
                        < field_square * state_count * state_count
                    ):
                        raise AssertionError(("Khintchine consequence", order))

                    # Odd |S| forces every field coordinate to be a nonzero
                    # odd integer. Coordinatewise contraction to magnitude one
                    # gives the exact integer-lattice lower bound mu_k.
                    if expected_parity and field_absolute_sum < mu_numerator(
                        complement_size
                    ):
                        raise AssertionError(("mu consequence", order))

                    fixed_x_records.append((abs(q), field_square))

                submaximum = max(value for value, _ in fixed_x_records)
                for absolute_q, field_square in fixed_x_records:
                    if absolute_q != submaximum:
                        continue
                    exact_lower_square = (
                        submaximum * submaximum
                        + field_square
                        + math.comb(complement_size, 2)
                    )
                    parity_lower_square = (
                        submaximum * submaximum
                        + math.comb(complement_size, 2)
                        + complement_size * (subset_size % 2)
                    )
                    if maximum * maximum < exact_lower_square:
                        raise AssertionError(("hereditary field bound", order))
                    if maximum * maximum < parity_lower_square:
                        raise AssertionError(("hereditary parity bound", order))

        if order_minimum is None:
            raise AssertionError(("no signings", order))
        minimum_maximum[order] = order_minimum

    if minimum_maximum != EXPECTED_F:
        raise AssertionError(
            ("independent small F values", minimum_maximum, EXPECTED_F)
        )
    return minimum_maximum, counts


def verify_f_consequences(values: dict[int, int]) -> int:
    checks = 0
    for total in range(1, MAX_ORDER + 1):
        edge_parity = math.comb(total, 2) % 2
        if values[total] % 2 != edge_parity:
            raise AssertionError(("F parity", total, values[total]))

        for subset_size in range(1, total + 1):
            complement_size = total - subset_size
            square_rhs = (
                values[subset_size] ** 2
                + math.comb(complement_size, 2)
                + complement_size * (subset_size % 2)
            )
            if values[total] ** 2 < square_rhs:
                raise AssertionError(
                    ("F hereditary square consequence", subset_size, total)
                )

            increment_numerator = (
                mu_numerator(complement_size) if subset_size % 2 else 0
            )
            if (
                (values[total] - values[subset_size])
                * (1 << complement_size)
                < increment_numerator
            ):
                raise AssertionError(
                    ("F mu consequence", subset_size, total)
                )

            if subset_size % 2:
                difference = values[total] - values[subset_size]
                if 2 * difference * difference < complement_size:
                    raise AssertionError(
                        ("F Khintchine consequence", subset_size, total)
                    )
            checks += 1
    return checks


def verify_h_subadditivity_failure(values: dict[int, int]) -> tuple[int, int, int]:
    """Give an exact counterexample to H(n+k) <= H(n) + H(k)."""
    left_order = 2
    right_order = 2
    total = left_order + right_order
    if values[left_order] != values[right_order]:
        raise AssertionError("the exact comparison expects equal child values")

    # Here H(2) + H(2) = 2 because F(2)=1. Cubing both positive
    # sides reduces H(4) > 2 to the integer comparison F(4)^2 > 8.
    child = values[left_order]
    if child != 1 or values[total] ** 2 <= 8 * child * child:
        raise AssertionError(("missing H subadditivity counterexample", values))
    return left_order, right_order, total


def expect_corruption_detected(
    label: str, corrupted_check: Callable[[], None]
) -> str:
    try:
        corrupted_check()
    except AssertionError:
        return label
    raise AssertionError(("corruption was not detected", label))


def verify_corruption_controls(values: dict[int, int]) -> tuple[str, ...]:
    def doubled_internal_walsh_term() -> None:
        complement_size = 2
        correct = math.comb(complement_size, 2)
        corrupted = complement_size * (complement_size - 1)
        if correct != corrupted:
            raise AssertionError("double-counted symmetric internal edges")

    def quadratic_term_treated_as_odd() -> None:
        edge_list = edges(4)
        coefficients = (1,) * len(edge_list)
        subset = (0, 1)
        complement = (2, 3)
        subset_spins = (1, 1)
        y = (1, 1)
        q, field, internal_edges = block_data(
            edge_list, coefficients, subset, complement, subset_spins
        )
        r = complement_energy(internal_edges, y)
        linear = sum(a * b for a, b in zip(field, y, strict=True))
        actual = energy(
            edge_list,
            coefficients,
            spins_on_partition(
                4, subset, complement, subset_spins, tuple(-v for v in y)
            ),
        )
        corrupted = q - r - linear
        if actual != corrupted:
            raise AssertionError("quadratic T term was incorrectly negated")

    def strengthened_mu_constant() -> None:
        length = 3
        all_one_field_sum = sum(
            abs(sum(y)) for y in sign_vectors(length)
        )
        corrupted_lower_sum = mu_numerator(length) + 1
        if all_one_field_sum < corrupted_lower_sum:
            raise AssertionError("strictly strengthened mu bound")

    def asserted_exact_h_subadditivity() -> None:
        child = values[2]
        if values[4] ** 2 > 8 * child * child:
            raise AssertionError("H exact subadditivity fails at 2+2")

    controls = (
        ("walsh_edge_double_count", doubled_internal_walsh_term),
        ("quadratic_flip_parity", quadratic_term_treated_as_odd),
        ("mu_strengthening", strengthened_mu_constant),
        ("h_exact_subadditivity", asserted_exact_h_subadditivity),
    )
    return tuple(
        expect_corruption_detected(label, check) for label, check in controls
    )


def format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    mu_checks = verify_mu_formula()
    values, counts = verify_cavity_inequalities()
    consequence_checks = verify_f_consequences(values)
    left, right, total = verify_h_subadditivity_failure(values)
    corruption_controls = verify_corruption_controls(values)

    f_receipt = ",".join(f"F({order})={values[order]}" for order in values)
    mu_receipt = ",".join(
        format_fraction(mu(length)) for length in range(MAX_ORDER + 1)
    )
    print(f"orders_exhausted=1..{MAX_ORDER}")
    print(f"signings_checked={counts['signings']}")
    print(f"partitions_checked={counts['partitions']}")
    print(f"fixed_x_checks={counts['fixed_x']}")
    print(f"paired_y_checks={counts['paired_y']}")
    print(f"mu_formula_checks={mu_checks}")
    print(f"f_consequence_checks={consequence_checks}")
    print(f"small_F_values={f_receipt}")
    print(f"mu_values_k_0_to_{MAX_ORDER}={mu_receipt}")
    print(f"h_exact_subadditivity_counterexample={left}+{right}->{total}")
    print(f"corruption_controls={','.join(corruption_controls)}")
    print("cavity_hereditary_verification=PASSED")


if __name__ == "__main__":
    main()

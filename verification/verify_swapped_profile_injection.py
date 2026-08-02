#!/usr/bin/env python3
"""Exact checks for the swapped-profile injection theorem.

The swapped profile uses projective absolute graph energies but signed full
spin pairs for the rectangular block.  It is not the augmented/projective
profile used by the exact max-plus theorem.  Nevertheless, a balanced map
and a one-sided gauge injection give a valid (generally weaker) composition
bound.  This script checks that construction exhaustively on every
switching-normalized signing at block sizes 2+2, 2+3, 2+4, and 3+3.

Only Python integer arithmetic is used.  No solver or random search enters.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from itertools import product

if not __debug__:
    raise RuntimeError("verification requires Python assertions")


Spin = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]
Gauge = tuple[Spin, Spin, int]
State = tuple[Spin, Spin, Spin, Spin]


SPLITS = ((2, 2), (2, 3), (2, 4), (3, 3))

# normalized cases, represented labelled cases, gauges, states, strict fiber
# dominations, and the complete (Lambda_swap, true gain) distribution.
EXPECTED_SPLITS = {
    (2, 2): (2, 64, 16, 128, 4, ((0, 0, 1), (0, 2, 1))),
    (2, 3): (8, 1_024, 128, 2_048, 64, ((2, 4, 8),)),
    (2, 4): (
        64,
        32_768,
        2_048,
        65_536,
        688,
        ((0, 2, 12), (0, 4, 12), (2, 4, 24), (4, 4, 14), (4, 6, 2)),
    ),
    (3, 3): (
        64,
        32_768,
        2_048,
        65_536,
        628,
        ((2, 6, 24), (4, 6, 40)),
    ),
}
EXPECTED_STREAM_SHA256 = (
    "638daefed306506cac5f7a724a64b4717d902601a80b2a5b7e73a9da768fbec9"
)


def projectivize(spin: Spin) -> Spin:
    return tuple(spin[0] * entry for entry in spin)


def projective_spins(order: int) -> tuple[Spin, ...]:
    return tuple((1, *tail) for tail in product((-1, 1), repeat=order - 1))


def full_spins(order: int) -> tuple[Spin, ...]:
    return tuple(product((-1, 1), repeat=order))


def switching_normalized_signings(order: int) -> tuple[Matrix, ...]:
    """One representative after making every edge from vertex zero positive."""

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


def switching_normalized_rectangles(rows: int, columns: int) -> tuple[Matrix, ...]:
    """One representative after making the first row and column positive."""

    residual_entries = tuple(
        (row, column)
        for row in range(1, rows)
        for column in range(1, columns)
    )
    result = []
    for signs in product((-1, 1), repeat=len(residual_entries)):
        matrix = [[1] * columns for _ in range(rows)]
        for (row, column), sign in zip(residual_entries, signs):
            matrix[row][column] = sign
        result.append(tuple(tuple(row) for row in matrix))
    return tuple(result)


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


def maximum_energy(matrix: Matrix) -> int:
    return max(
        abs(quadratic_energy(matrix, spin))
        for spin in projective_spins(len(matrix))
    )


def maximum_energy_witness(matrix: Matrix) -> tuple[int, Spin, int]:
    best = -1
    best_spin: Spin | None = None
    best_orientation = 1
    for spin in projective_spins(len(matrix)):
        value = quadratic_energy(matrix, spin)
        if abs(value) > best:
            best = abs(value)
            best_spin = spin
            best_orientation = 1 if value >= 0 else -1
    assert best_spin is not None
    return best, best_spin, best_orientation


def coordinate_class(left: Spin, right: Spin) -> Spin:
    assert len(left) == len(right)
    return projectivize(tuple(a * b for a, b in zip(left, right)))


def swapped_gauge(state: State) -> Gauge:
    graph_left, graph_right, cross_left, cross_right = state
    return (
        coordinate_class(graph_left, cross_left),
        coordinate_class(graph_right, cross_right),
        cross_left[0],
    )


def aligned_full_matrix(
    graph_left: Matrix,
    graph_right: Matrix,
    rectangle: Matrix,
    gauge: Gauge,
) -> Matrix:
    left_gauge, right_gauge, relative_sign = gauge
    left_order = len(graph_left)
    right_order = len(graph_right)
    order = left_order + right_order
    result = [[0] * order for _ in range(order)]
    for row in range(left_order):
        for column in range(left_order):
            result[row][column] = (
                left_gauge[row]
                * graph_left[row][column]
                * left_gauge[column]
            )
    for row in range(right_order):
        for column in range(right_order):
            result[left_order + row][left_order + column] = (
                relative_sign
                * right_gauge[row]
                * graph_right[row][column]
                * right_gauge[column]
            )
    for row in range(left_order):
        for column in range(right_order):
            value = rectangle[row][column]
            result[row][left_order + column] = value
            result[left_order + column][row] = value
    return tuple(tuple(row) for row in result)


def swapped_energy(
    graph_left: Matrix,
    graph_right: Matrix,
    rectangle: Matrix,
    state: State,
) -> int:
    left_spin, right_spin, cross_left, cross_right = state
    return (
        abs(quadratic_energy(graph_left, left_spin))
        + abs(quadratic_energy(graph_right, right_spin))
        + rectangular_energy(rectangle, cross_left, cross_right)
    )


def injected_state(
    gauge: Gauge,
    maximizing_spin: Spin,
    orientation: int,
    left_order: int,
) -> State:
    """Construct the swapped state that dominates one gauge maximum."""

    left_gauge, right_gauge, relative_sign = gauge
    left_spin = maximizing_spin[:left_order]
    right_spin = maximizing_spin[left_order:]
    graph_left = coordinate_class(left_gauge, left_spin)
    graph_right = coordinate_class(right_gauge, right_spin)
    common_sign = relative_sign * orientation * left_spin[0]
    cross_left = tuple(
        common_sign * orientation * entry for entry in left_spin
    )
    cross_right = tuple(common_sign * entry for entry in right_spin)
    assert cross_left[0] == relative_sign
    return graph_left, graph_right, cross_left, cross_right


def order_statistic(values: list[int], rank: int) -> int:
    assert 1 <= rank <= len(values)
    return sorted(values)[rank - 1]


def graph_maximum(matrix: Matrix) -> int:
    return max(
        abs(quadratic_energy(matrix, spin))
        for spin in projective_spins(len(matrix))
    )


def rectangle_maximum(matrix: Matrix) -> int:
    return max(
        abs(rectangular_energy(matrix, left, right))
        for left in projective_spins(len(matrix))
        for right in projective_spins(len(matrix[0]))
    )


def check_case(
    graph_left: Matrix,
    graph_right: Matrix,
    rectangle: Matrix,
) -> tuple[int, int, int, int, int, int]:
    left_order = len(graph_left)
    right_order = len(graph_right)
    total_order = left_order + right_order
    graph_left_states = projective_spins(left_order)
    graph_right_states = projective_spins(right_order)
    cross_left_states = full_spins(left_order)
    cross_right_states = full_spins(right_order)

    independent_ceiling = (
        graph_maximum(graph_left)
        + graph_maximum(graph_right)
        + rectangle_maximum(rectangle)
    )
    fibers: dict[Gauge, list[tuple[State, int]]] = {}
    deficits = []
    for left_spin in graph_left_states:
        for right_spin in graph_right_states:
            for cross_left in cross_left_states:
                for cross_right in cross_right_states:
                    state = (left_spin, right_spin, cross_left, cross_right)
                    energy = swapped_energy(
                        graph_left, graph_right, rectangle, state
                    )
                    assert energy <= independent_ceiling
                    fibers.setdefault(swapped_gauge(state), []).append((state, energy))
                    deficits.append(independent_ceiling - energy)

    gauge_count = 1 << (total_order - 1)
    fiber_size = 1 << (total_order - 1)
    assert len(fibers) == gauge_count
    assert all(len(entries) == fiber_size for entries in fibers.values())
    assert len(deficits) == 1 << (2 * total_order - 2)

    gauge_maxima: dict[Gauge, int] = {}
    injected_states: set[State] = set()
    strict_fiber_dominations = 0
    for gauge, entries in fibers.items():
        aligned = aligned_full_matrix(
            graph_left, graph_right, rectangle, gauge
        )
        gauge_maximum, spin, orientation = maximum_energy_witness(aligned)
        gauge_maxima[gauge] = gauge_maximum
        state = injected_state(gauge, spin, orientation, left_order)
        assert swapped_gauge(state) == gauge
        injected_energy = swapped_energy(
            graph_left, graph_right, rectangle, state
        )
        assert injected_energy >= gauge_maximum
        assert state in {candidate for candidate, _energy in entries}
        assert state not in injected_states
        injected_states.add(state)
        fiber_maximum = max(energy for _state, energy in entries)
        assert fiber_maximum >= injected_energy >= gauge_maximum
        strict_fiber_dominations += int(fiber_maximum > gauge_maximum)

    assert len(injected_states) == gauge_count
    true_gain = independent_ceiling - min(gauge_maxima.values())
    swapped_lambda = order_statistic(deficits, gauge_count)
    assert swapped_lambda <= true_gain

    # Corruption control: deleting the cross-left sign from the gauge loses
    # the relative-sign coordinate, halves the number of fibers, and doubles
    # their size.  It cannot support the rank used above.
    corrupted = Counter(
        (gauge[0], gauge[1])
        for gauge, entries in fibers.items()
        for _entry in entries
    )
    assert len(corrupted) == gauge_count // 2
    assert set(corrupted.values()) == {2 * fiber_size}

    return (
        gauge_count,
        len(deficits),
        swapped_lambda,
        true_gain,
        strict_fiber_dominations,
        independent_ceiling,
    )


def represented_labelled_cases(left_order: int, right_order: int) -> int:
    edge_count = (
        left_order * (left_order - 1) // 2
        + right_order * (right_order - 1) // 2
        + left_order * right_order
    )
    return 1 << edge_count


def main() -> None:
    total_normalized_cases = 0
    total_labelled_cases = 0
    total_gauges = 0
    total_states = 0
    total_strict_fiber_dominations = 0
    digest = hashlib.sha256()

    for left_order, right_order in SPLITS:
        left_signings = switching_normalized_signings(left_order)
        right_signings = switching_normalized_signings(right_order)
        rectangles = switching_normalized_rectangles(left_order, right_order)
        case_count = 0
        split_gauges = 0
        split_states = 0
        split_strict = 0
        split_pairs: Counter[tuple[int, int]] = Counter()
        for left_index, graph_left in enumerate(left_signings):
            for right_index, graph_right in enumerate(right_signings):
                for rectangle_index, rectangle in enumerate(rectangles):
                    checked = check_case(graph_left, graph_right, rectangle)
                    (
                        gauge_count,
                        state_count,
                        swapped_lambda,
                        true_gain,
                        strict_dominations,
                        independent_ceiling,
                    ) = checked
                    case_count += 1
                    split_gauges += gauge_count
                    split_states += state_count
                    split_strict += strict_dominations
                    split_pairs[(swapped_lambda, true_gain)] += 1
                    digest.update(
                        (
                            f"{left_order},{right_order},{left_index},"
                            f"{right_index},{rectangle_index},"
                            f"{swapped_lambda},{true_gain},"
                            f"{strict_dominations},{independent_ceiling}\n"
                        ).encode("ascii")
                    )

        expected_normalized = 1 << (
            (left_order - 1) * (left_order - 2) // 2
            + (right_order - 1) * (right_order - 2) // 2
            + (left_order - 1) * (right_order - 1)
        )
        labelled_cases = represented_labelled_cases(left_order, right_order)
        assert case_count == expected_normalized
        assert labelled_cases == case_count * (1 << (2 * (left_order + right_order) - 3))
        expected = EXPECTED_SPLITS[(left_order, right_order)]
        assert (
            case_count,
            labelled_cases,
            split_gauges,
            split_states,
            split_strict,
        ) == expected[:5]
        assert split_pairs == Counter(
            {(profile_gain, true_gain): count for profile_gain, true_gain, count in expected[5]}
        )
        total_normalized_cases += case_count
        total_labelled_cases += labelled_cases
        total_gauges += split_gauges
        total_states += split_states
        total_strict_fiber_dominations += split_strict
        print(
            f"split={left_order}+{right_order} "
            f"normalized_cases={case_count} represented_labelled={labelled_cases} "
            f"gauges={split_gauges} states={split_states} "
            f"strict_fiber_dominations={split_strict} "
            "lambda_gain_pairs="
            + ",".join(
                f"{profile_gain}/{true_gain}:{count}"
                for (profile_gain, true_gain), count in sorted(split_pairs.items())
            )
        )

    # These assertions turn truncation, convention changes, and enumeration
    # corruption into hard failures.  The digest is filled after an
    # independent first run of this exact deterministic stream.
    assert total_normalized_cases == 138
    assert total_labelled_cases == 66_624
    assert total_gauges == 4_240
    assert total_states == 133_248
    assert total_strict_fiber_dominations == 1_384
    assert digest.hexdigest() == EXPECTED_STREAM_SHA256
    print(f"normalized_cases_checked={total_normalized_cases}")
    print(f"represented_labelled_cases={total_labelled_cases}")
    print(f"balanced_gauges_checked={total_gauges}")
    print(f"swapped_states_checked={total_states}")
    print(f"strict_fiber_dominations={total_strict_fiber_dominations}")
    print(f"stream_sha256={digest.hexdigest()}")
    print("corruption_controls=relative_sign_omission,max_plus_strictness")
    print("swapped_profile_injection=PASSED")


if __name__ == "__main__":
    main()

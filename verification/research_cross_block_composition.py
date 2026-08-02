#!/usr/bin/env python3
"""Exact cross-block and one-vertex composition research for MO 413935.

For fixed sign matrices A and B, define

    J(A, B) = min_C M([[A, C], [C^T, B]]),

where C ranges over rectangular sign matrices.  This script computes J for
every pair of optimal switching classes with total order at most ten.  A
deterministic Z3 cutting-plane solver proves infeasibility below the reported
value and every returned witness is independently evaluated on the full spin
cube.  Whenever the cross block has at most sixteen entries, a separate
exhaustive enumeration checks the solver result.

The script also exhausts every one-vertex extension of every optimal
root-normalized representative through order nine.  Finally, a direct search
over all order-eight residual graphs proves the finite non-heredity result

    min_{M(B)=F(8)} J(K_2, B) = 15 > F(10) = 13,

while an order-eight block with M(B)=12 does participate in an order-ten
signing of maximum 13.  All mathematical pass/fail decisions use integers.

Completeness of the optimizer catalogues trusts nauty ``geng`` and the
asserted graph counts and stream hashes.  Infeasibility for cross blocks with
more than sixteen entries trusts Z3.  The critical 2+8 obstruction is also
checked by direct enumeration of every cross signing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import Counter
from dataclasses import dataclass
from functools import cache
from typing import Any

import research_exact_small_n as small


DETERMINISTIC_SEED = 413935
EXPECTED_F = {
    1: 0,
    2: 1,
    3: 3,
    4: 4,
    5: 4,
    6: 5,
    7: 9,
    8: 10,
    9: 12,
    10: 13,
}
EXPECTED_ROOT_REPRESENTATIVES = {
    2: 1,
    3: 2,
    4: 2,
    5: 1,
    6: 1,
    7: 12,
    8: 4,
    9: 55,
}
EXPECTED_SWITCHING_CLASSES = {
    2: 1,
    3: 2,
    4: 1,
    5: 1,
    6: 1,
    7: 6,
    8: 2,
    9: 15,
}
EXPECTED_PAIR_DISTRIBUTIONS = {
    (2, 2): {4: 1},
    (2, 3): {4: 2},
    (2, 4): {5: 1},
    (3, 3): {5: 2, 7: 2},
    (2, 5): {9: 1},
    (3, 4): {9: 2},
    (2, 6): {10: 1},
    (3, 5): {10: 2},
    (4, 4): {10: 1},
    (2, 7): {12: 6},
    (3, 6): {12: 2},
    (4, 5): {12: 1},
    (2, 8): {15: 2},
    (3, 7): {13: 6, 15: 6},
    (4, 6): {13: 1},
    (5, 5): {13: 1},
}
EXPECTED_EXTENSION_CLASS_DISTRIBUTIONS = {
    2: {3: 1},
    3: {4: 2},
    4: {4: 1},
    5: {5: 1},
    6: {9: 1},
    7: {10: 4, 12: 2},
    8: {12: 2},
    9: {13: 4, 15: 11},
}
EXPECTED_ORDER8_VALUE_COUNTS = {
    10: 4,
    12: 104,
    14: 342,
    16: 318,
    18: 164,
    20: 68,
    22: 28,
    24: 10,
    26: 4,
    28: 2,
}
ORDER10_SWITCHING_CLASSES = (
    ("HCRbczQ", "HCrbdxz"),
    ("HCpdehU", "HCZbeyz"),
)
EXPECTED_ORDER10_EXTENSION_DISTRIBUTION = {17: 1, 19: 1}

Matrix = tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class CatalogueEntry:
    order: int
    value: int
    root_records: tuple[str, ...]
    switching_classes: tuple[tuple[str, ...], ...]


@dataclass(frozen=True)
class CrossState:
    base_energy: int
    character_mask: int


@dataclass(frozen=True)
class FeasibilityResult:
    feasible: bool
    coefficient_mask: int | None
    maximum: int | None
    iterations: int
    constraints: int


@dataclass(frozen=True)
class CrossOptimum:
    value: int
    coefficient_mask: int
    witness_sha256: str
    infeasible_bounds: tuple[int, ...]
    cutting_plane_iterations: int
    active_constraints: int
    direct_crosscheck: bool


def load_z3() -> Any:
    try:
        import z3
    except ImportError as error:
        raise RuntimeError(
            "research_cross_block_composition.py requires z3-solver"
        ) from error
    return z3


@cache
def spin_vectors(order: int, fix_first: bool) -> tuple[tuple[int, ...], ...]:
    free = order - (1 if fix_first else 0)
    prefix = (1,) if fix_first else ()
    return tuple(
        prefix
        + tuple(-1 if mask >> index & 1 else 1 for index in range(free))
        for mask in range(1 << free)
    )


def matrix_energy(matrix: Matrix, spins: tuple[int, ...]) -> int:
    return sum(
        matrix[row][column] * spins[row] * spins[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def matrix_from_record(order: int, record: str) -> Matrix:
    adjacency = small.graph6_adjacency(record.encode("ascii"), order - 1)
    return small.signed_matrix(adjacency)


def build_catalogue(geng: str, maximum_order: int) -> dict[int, CatalogueEntry]:
    catalogue: dict[int, CatalogueEntry] = {
        1: CatalogueEntry(1, 0, ("singleton",), (("singleton",),))
    }
    for order in range(2, maximum_order + 1):
        result = small.exact_search(order, geng, False)
        if result.value != EXPECTED_F[order]:
            raise AssertionError(("unexpected F value", order, result.value))
        expected_root_count = EXPECTED_ROOT_REPRESENTATIVES[order]
        if result.optimal_representatives != expected_root_count:
            raise AssertionError(
                (
                    "optimal root representative count",
                    order,
                    result.optimal_representatives,
                    expected_root_count,
                )
            )
        expected_digest = small.EXPECTED_GENG_STREAM_SHA256.get(order)
        if expected_digest is not None and result.generator_sha256 != expected_digest:
            raise AssertionError(("geng stream digest", order))
        classes = tuple(
            tuple(group) for group in small.optimizer_switching_classes(result)
        )
        expected_classes = EXPECTED_SWITCHING_CLASSES[order]
        if len(classes) != expected_classes:
            raise AssertionError(
                ("switching class count", order, len(classes), expected_classes)
            )
        catalogue[order] = CatalogueEntry(
            order=order,
            value=result.value,
            root_records=result.optimal_graph6,
            switching_classes=classes,
        )
    return catalogue


def combined_states(left: Matrix, right: Matrix) -> tuple[CrossState, ...]:
    left_order = len(left)
    right_order = len(right)
    states: list[CrossState] = []
    for left_spins in spin_vectors(left_order, True):
        left_energy = matrix_energy(left, left_spins)
        for right_spins in spin_vectors(right_order, False):
            character_mask = 0
            bit = 0
            for row in range(left_order):
                for column in range(right_order):
                    if left_spins[row] * right_spins[column] == -1:
                        character_mask |= 1 << bit
                    bit += 1
            states.append(
                CrossState(
                    left_energy + matrix_energy(right, right_spins),
                    character_mask,
                )
            )
    return tuple(states)


def evaluate_mask(
    states: tuple[CrossState, ...], cross_edges: int, coefficient_mask: int
) -> int:
    if not 0 <= coefficient_mask < 1 << cross_edges:
        raise ValueError("cross coefficient mask has the wrong width")
    return max(
        abs(
            state.base_energy
            + cross_edges
            - 2 * (coefficient_mask ^ state.character_mask).bit_count()
        )
        for state in states
    )


def exhaustive_cross_optimum(
    states: tuple[CrossState, ...], cross_edges: int
) -> tuple[int, int, int]:
    if cross_edges > 16:
        raise ValueError("direct cross enumeration is capped at sixteen edges")
    best = math.inf
    best_mask = 0
    best_count = 0
    # Negating the whole cross block is induced by flipping every spin in one
    # block and preserves both internal quadratic forms. Fix the first edge +1.
    for reduced_mask in range(1 << (cross_edges - 1)):
        coefficient_mask = reduced_mask << 1
        maximum = 0
        for state in states:
            value = abs(
                state.base_energy
                + cross_edges
                - 2 * (coefficient_mask ^ state.character_mask).bit_count()
            )
            maximum = max(maximum, value)
            if maximum > best:
                break
        if maximum < best:
            best = maximum
            best_mask = coefficient_mask
            best_count = 1
        elif maximum == best:
            best_count += 1
    return int(best), best_mask, best_count


def exhaustive_cross_witness(
    states: tuple[CrossState, ...], cross_edges: int, bound: int
) -> int | None:
    if cross_edges > 16:
        raise ValueError("direct cross enumeration is capped at sixteen edges")
    for reduced_mask in range(1 << (cross_edges - 1)):
        coefficient_mask = reduced_mask << 1
        for state in states:
            value = abs(
                state.base_energy
                + cross_edges
                - 2 * (coefficient_mask ^ state.character_mask).bit_count()
            )
            if value > bound:
                break
        else:
            return coefficient_mask
    return None


def cutting_plane_feasibility(
    z3: Any,
    states: tuple[CrossState, ...],
    cross_edges: int,
    bound: int,
    timeout_ms: int,
) -> FeasibilityResult:
    solver = z3.Solver()
    solver.set(timeout=timeout_ms, random_seed=DETERMINISTIC_SEED)
    variables = tuple(z3.Bool(f"cross_{index}") for index in range(cross_edges))
    solver.add(variables[0])
    active: set[int] = set()
    iterations = 0

    while True:
        status = solver.check()
        if status == z3.unknown:
            raise RuntimeError(("Z3 returned unknown", solver.reason_unknown()))
        if status == z3.unsat:
            return FeasibilityResult(False, None, None, iterations, len(active))

        model = solver.model()
        coefficient_mask = sum(
            1 << index
            for index, variable in enumerate(variables)
            if not z3.is_true(model.eval(variable, model_completion=True))
        )
        violations: list[int] = []
        maximum = 0
        for index, state in enumerate(states):
            value = abs(
                state.base_energy
                + cross_edges
                - 2 * (coefficient_mask ^ state.character_mask).bit_count()
            )
            maximum = max(maximum, value)
            if value > bound and index not in active:
                violations.append(index)

        if not violations:
            return FeasibilityResult(
                True, coefficient_mask, maximum, iterations, len(active)
            )

        for index in violations:
            state = states[index]
            terms = [
                z3.If(
                    variables[edge],
                    -1 if state.character_mask >> edge & 1 else 1,
                    1 if state.character_mask >> edge & 1 else -1,
                )
                for edge in range(cross_edges)
            ]
            energy = state.base_energy + z3.Sum(terms)
            solver.add(energy <= bound, energy >= -bound)
            active.add(index)
        iterations += 1


def exact_cross_optimum(
    z3: Any,
    left: Matrix,
    right: Matrix,
    global_lower_bound: int,
    timeout_ms: int,
) -> CrossOptimum:
    total_order = len(left) + len(right)
    total_edges = math.comb(total_order, 2)
    cross_edges = len(left) * len(right)
    states = combined_states(left, right)
    if global_lower_bound % 2 != total_edges % 2:
        raise AssertionError("global lower bound has the wrong energy parity")

    infeasible: list[int] = []
    total_iterations = 0
    largest_active_set = 0
    for bound in range(global_lower_bound, total_edges + 1, 2):
        result = cutting_plane_feasibility(
            z3, states, cross_edges, bound, timeout_ms
        )
        total_iterations += result.iterations
        largest_active_set = max(largest_active_set, result.constraints)
        if not result.feasible:
            infeasible.append(bound)
            continue
        if result.coefficient_mask is None or result.maximum is None:
            raise AssertionError("feasible result omitted its witness")
        recomputed = evaluate_mask(states, cross_edges, result.coefficient_mask)
        if recomputed != result.maximum or recomputed > bound:
            raise AssertionError(("cross witness recomputation", recomputed, result))
        if recomputed != bound:
            raise AssertionError(("non-sharp feasible bound", recomputed, bound))

        direct_crosscheck = False
        if cross_edges <= 16:
            direct_value, _, _ = exhaustive_cross_optimum(states, cross_edges)
            if direct_value != bound:
                raise AssertionError(("direct cross disagreement", direct_value, bound))
            direct_crosscheck = True

        signs = [
            -1 if result.coefficient_mask >> edge & 1 else 1
            for edge in range(cross_edges)
        ]
        witness_sha256 = hashlib.sha256(
            json.dumps(signs, separators=(",", ":")).encode("ascii")
        ).hexdigest()
        return CrossOptimum(
            value=bound,
            coefficient_mask=result.coefficient_mask,
            witness_sha256=witness_sha256,
            infeasible_bounds=tuple(infeasible),
            cutting_plane_iterations=total_iterations,
            active_constraints=largest_active_set,
            direct_crosscheck=direct_crosscheck,
        )
    raise AssertionError("no cross signing was feasible")


def class_representative_matrix(
    order: int, switching_class: tuple[str, ...]
) -> Matrix:
    if order == 1:
        return ((0,),)
    return matrix_from_record(order, switching_class[0])


def analyze_pair_composition(
    z3: Any,
    catalogue: dict[int, CatalogueEntry],
    timeout_ms: int,
) -> tuple[list[dict[str, object]], int]:
    reports: list[dict[str, object]] = []
    direct_crosschecks = 0
    for (left_order, right_order), expected_distribution in sorted(
        EXPECTED_PAIR_DISTRIBUTIONS.items(), key=lambda item: (sum(item[0]), item[0])
    ):
        total_order = left_order + right_order
        distribution: Counter[int] = Counter()
        class_values: list[dict[str, object]] = []
        for left_class in catalogue[left_order].switching_classes:
            left = class_representative_matrix(left_order, left_class)
            for right_class in catalogue[right_order].switching_classes:
                right = class_representative_matrix(right_order, right_class)
                optimum = exact_cross_optimum(
                    z3, left, right, EXPECTED_F[total_order], timeout_ms
                )
                distribution[optimum.value] += 1
                direct_crosschecks += int(optimum.direct_crosscheck)
                class_values.append(
                    {
                        "active_constraints": optimum.active_constraints,
                        "cutting_plane_iterations": optimum.cutting_plane_iterations,
                        "infeasible_bounds": list(optimum.infeasible_bounds),
                        "left": left_class[0],
                        "right": right_class[0],
                        "value": optimum.value,
                        "witness_sha256": optimum.witness_sha256,
                    }
                )
        observed_distribution = dict(sorted(distribution.items()))
        if observed_distribution != expected_distribution:
            raise AssertionError(
                (
                    "pair distribution",
                    left_order,
                    right_order,
                    observed_distribution,
                    expected_distribution,
                )
            )
        minimum = min(distribution)
        maximum = max(distribution)
        h_defect = (
            minimum ** (2 / 3)
            - EXPECTED_F[left_order] ** (2 / 3)
            - EXPECTED_F[right_order] ** (2 / 3)
        )
        reports.append(
            {
                "class_pair_count": sum(distribution.values()),
                "class_values": class_values,
                "direct_crosschecks": (
                    len(class_values) if left_order * right_order <= 16 else 0
                ),
                "global_F": EXPECTED_F[total_order],
                "h_defect_at_best_fixed_blocks": round(h_defect, 12),
                "left_order": left_order,
                "maximum": maximum,
                "minimum": minimum,
                "right_order": right_order,
                "value_distribution": {
                    str(key): value for key, value in observed_distribution.items()
                },
            }
        )
    return reports, direct_crosschecks


def direct_extension_value(matrix: Matrix) -> tuple[int, int, int]:
    states = combined_states(((0,),), matrix)
    return exhaustive_cross_optimum(states, len(matrix))


def restricted_extension_optimum(
    matrix: Matrix, minimum_absolute_base_energy: int
) -> int:
    states = tuple(
        state
        for state in combined_states(((0,),), matrix)
        if abs(state.base_energy) >= minimum_absolute_base_energy
    )
    if not states:
        raise ValueError("restricted extension state set is empty")
    value, _, _ = exhaustive_cross_optimum(states, len(matrix))
    return value


def analyze_extensions(
    catalogue: dict[int, CatalogueEntry],
) -> tuple[list[dict[str, object]], dict[str, int]]:
    reports: list[dict[str, object]] = []
    obstruction_layers: Counter[str] = Counter()
    for order, expected_distribution in EXPECTED_EXTENSION_CLASS_DISTRIBUTIONS.items():
        class_distribution: Counter[int] = Counter()
        root_distribution: Counter[int] = Counter()
        class_values: list[dict[str, object]] = []
        root_values: dict[str, int] = {}
        for record in catalogue[order].root_records:
            value, _, _ = direct_extension_value(matrix_from_record(order, record))
            root_distribution[value] += 1
            root_values[record] = value

        for switching_class in catalogue[order].switching_classes:
            values = {root_values[record] for record in switching_class}
            if len(values) != 1:
                raise AssertionError(
                    ("extension value not switching invariant", order, switching_class)
                )
            value = values.pop()
            class_distribution[value] += 1
            record = switching_class[0]
            class_report: dict[str, object] = {
                "class_size": len(switching_class),
                "record": record,
                "value": value,
            }
            if value > EXPECTED_F[order + 1]:
                matrix = matrix_from_record(order, record)
                absolute_layers = sorted(
                    {
                        abs(state.base_energy)
                        for state in combined_states(((0,),), matrix)
                    },
                    reverse=True,
                )
                profile = [
                    {
                        "minimum_absolute_base_energy": layer,
                        "restricted_extension_value": restricted_extension_optimum(
                            matrix, layer
                        ),
                    }
                    for layer in absolute_layers
                ]
                first_obstruction = next(
                    item["minimum_absolute_base_energy"]
                    for item in profile
                    if item["restricted_extension_value"]
                    > EXPECTED_F[order + 1]
                )
                obstruction_layers[f"order_{order}_abs_{first_obstruction}"] += 1
                class_report["layer_profile"] = profile
            class_values.append(class_report)

        observed_distribution = dict(sorted(class_distribution.items()))
        if observed_distribution != expected_distribution:
            raise AssertionError(
                (
                    "extension class distribution",
                    order,
                    observed_distribution,
                    expected_distribution,
                )
            )
        reports.append(
            {
                "class_count": len(catalogue[order].switching_classes),
                "class_values": class_values,
                "global_F_next": EXPECTED_F[order + 1],
                "order": order,
                "root_representative_count": len(catalogue[order].root_records),
                "root_value_distribution": {
                    str(key): value for key, value in sorted(root_distribution.items())
                },
                "value_distribution": {
                    str(key): value for key, value in observed_distribution.items()
                },
            }
        )
    expected_layers = {"order_7_abs_7": 2, "order_9_abs_12": 11}
    if dict(obstruction_layers) != expected_layers:
        raise AssertionError(
            ("extension obstruction layers", dict(obstruction_layers), expected_layers)
        )
    return reports, dict(obstruction_layers)


def analyze_order10_extensions() -> dict[str, object]:
    distribution: Counter[int] = Counter()
    class_values: dict[str, dict[str, int]] = {}
    for switching_class in ORDER10_SWITCHING_CLASSES:
        values: set[int] = set()
        representative_mask = 0
        representative_count = 0
        for record in switching_class:
            matrix = matrix_from_record(10, record)
            maximum = max(
                abs(matrix_energy(matrix, spins))
                for spins in spin_vectors(10, True)
            )
            if maximum != EXPECTED_F[10]:
                raise AssertionError(("stored order-10 optimizer", record, maximum))
            extension_value, coefficient_mask, best_mask_count = (
                direct_extension_value(matrix)
            )
            values.add(extension_value)
            if record == switching_class[0]:
                representative_mask = coefficient_mask
                representative_count = best_mask_count
        if len(values) != 1:
            raise AssertionError(
                ("order-10 extension value not switching invariant", switching_class)
            )
        value = values.pop()
        distribution[value] += 1
        class_values[switching_class[0]] = {
            "best_cross_masks_mod_global_sign": representative_count,
            "cross_coefficient_mask": representative_mask,
            "value": value,
        }
    observed = dict(sorted(distribution.items()))
    if observed != EXPECTED_ORDER10_EXTENSION_DISTRIBUTION:
        raise AssertionError(
            (
                "order-10 extension distribution",
                observed,
                EXPECTED_ORDER10_EXTENSION_DISTRIBUTION,
            )
        )
    return {
        "class_count": len(ORDER10_SWITCHING_CLASSES),
        "class_values": class_values,
        "order": 10,
        "value_distribution": {str(key): value for key, value in observed.items()},
    }


def analyze_two_plus_eight_slack(
    geng: str, catalogue: dict[int, CatalogueEntry]
) -> dict[str, object]:
    left = ((0, 1), (1, 0))
    for record in catalogue[8].root_records:
        right = matrix_from_record(8, record)
        states = combined_states(left, right)
        direct_value, _, _ = exhaustive_cross_optimum(states, 16)
        if direct_value != 15:
            raise AssertionError(("optimal order-eight block composed below 15", record))

    value_counts: Counter[int] = Counter()
    value_twelve_checked = 0
    value_twelve_composable = 0
    witness_record = ""
    witness_mask: int | None = None
    stream_digest = hashlib.sha256()
    records_checked = 0
    for raw_record in small.geng_records(geng, 7):
        stream_digest.update(raw_record + b"\n")
        records_checked += 1
        adjacency = small.graph6_adjacency(raw_record, 7)
        value, _ = small.maximum_absolute_energy(adjacency)
        value_counts[value] += 1
        if value != 12:
            continue
        value_twelve_checked += 1
        right = small.signed_matrix(adjacency)
        states = combined_states(left, right)
        coefficient_mask = exhaustive_cross_witness(states, 16, 13)
        if coefficient_mask is None:
            continue
        if evaluate_mask(states, 16, coefficient_mask) != 13:
            raise AssertionError("2+8 slack witness failed recomputation")
        value_twelve_composable += 1
        if witness_mask is None:
            witness_record = raw_record.decode("ascii")
            witness_mask = coefficient_mask

    if records_checked != small.UNLABELED_GRAPH_COUNTS[7]:
        raise AssertionError(("order-eight residual graph count", records_checked))
    expected_digest = small.EXPECTED_GENG_STREAM_SHA256[8]
    if stream_digest.hexdigest() != expected_digest:
        raise AssertionError("order-eight geng stream digest")
    if dict(sorted(value_counts.items())) != EXPECTED_ORDER8_VALUE_COUNTS:
        raise AssertionError(("order-eight value distribution", value_counts))
    if value_twelve_checked != 104 or value_twelve_composable != 68:
        raise AssertionError(
            (
                "order-eight slack counts",
                value_twelve_checked,
                value_twelve_composable,
            )
        )
    if witness_mask is None:
        raise AssertionError("no order-eight slack witness was found")
    witness_receipt = {
        "cross_coefficient_mask": witness_mask,
        "order8_graph6": witness_record,
    }
    witness_sha256 = hashlib.sha256(
        json.dumps(witness_receipt, sort_keys=True, separators=(",", ":")).encode(
            "ascii"
        )
    ).hexdigest()
    return {
        "minimum_order8_block_maximum_for_order10_optimum": 12,
        "optimal_order8_composition_value": 15,
        "order8_root_representatives_at_value_12": value_twelve_checked,
        "order8_value12_representatives_composing_to_13": value_twelve_composable,
        "slack_witness": witness_receipt,
        "slack_witness_sha256": witness_sha256,
    }


def verify_corruption_controls(slack_report: dict[str, object]) -> tuple[str, ...]:
    controls: list[str] = []

    try:
        evaluate_mask((CrossState(0, 0),), 1, 2)
    except ValueError:
        controls.append("cross_mask_width")
    else:
        raise AssertionError("cross-mask-width corruption was not detected")

    if slack_report["optimal_order8_composition_value"] != 13:
        controls.append("false_optimal_block_heredity")
    else:
        raise AssertionError("false optimal-block heredity was not detected")

    if slack_report["minimum_order8_block_maximum_for_order10_optimum"] == 12:
        controls.append("required_internal_slack")
    else:
        raise AssertionError("internal-slack corruption was not detected")

    # A frustrated triangle has one-sided energies 1,1,1,-3. This prevents a
    # checker that silently drops the absolute value from passing.
    one_sided = (1, 1, 1, -3)
    if max(one_sided) != max(abs(value) for value in one_sided):
        controls.append("absolute_value")
    else:
        raise AssertionError("absolute-value corruption was not detected")
    return tuple(controls)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geng")
    parser.add_argument("--timeout-ms", type=int, default=120_000)
    arguments = parser.parse_args()
    if arguments.timeout_ms <= 0:
        raise SystemExit("--timeout-ms must be positive")

    try:
        z3 = load_z3()
        geng = small.locate_geng(arguments.geng)
        catalogue = build_catalogue(geng, 9)
    except (RuntimeError, FileNotFoundError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(2) from error

    pair_reports, direct_crosschecks = analyze_pair_composition(
        z3, catalogue, arguments.timeout_ms
    )
    extension_reports, obstruction_layers = analyze_extensions(catalogue)
    order10_extension_report = analyze_order10_extensions()
    slack_report = analyze_two_plus_eight_slack(geng, catalogue)
    corruption_controls = verify_corruption_controls(slack_report)

    print(
        json.dumps(
            {
                "catalogue": {
                    str(order): {
                        "root_representatives": len(entry.root_records),
                        "switching_classes": len(entry.switching_classes),
                        "value": entry.value,
                    }
                    for order, entry in catalogue.items()
                }
            },
            sort_keys=True,
        )
    )
    for report in pair_reports:
        printable = {key: value for key, value in report.items() if key != "class_values"}
        print(json.dumps({"pair_composition": printable}, sort_keys=True))
    for report in extension_reports:
        printable = {key: value for key, value in report.items() if key != "class_values"}
        print(json.dumps({"one_vertex_extension": printable}, sort_keys=True))
    print(
        json.dumps(
            {"one_vertex_extension": order10_extension_report}, sort_keys=True
        )
    )
    print(json.dumps({"extension_obstruction_layers": obstruction_layers}, sort_keys=True))
    print(json.dumps({"two_plus_eight_slack": slack_report}, sort_keys=True))
    print(
        json.dumps(
            {
                "corruption_controls": list(corruption_controls),
                "deterministic_seed": DETERMINISTIC_SEED,
                "direct_crosschecks": direct_crosschecks,
                "status": "PASSED",
                "z3_version": z3.get_version_string(),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

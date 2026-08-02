#!/usr/bin/env python3
"""Finite checks for RESEARCH_CONTINUATION.md.

The analytic entropy step uses the data-processing inequality and is not
proved by this script.  These checks independently enumerate the cavity and
block partition identities, and test every algebraic part of the scaled
Gaussian covariance and determinant bound.
"""

from __future__ import annotations

import itertools
import math
import random
from collections.abc import Iterable


DETERMINISTIC_SEED = 413935
TOLERANCE = 2e-10


def edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(itertools.combinations(range(n), 2))


def signing_from_mask(n: int, mask: int) -> dict[tuple[int, int], int]:
    return {
        edge: (-1 if mask >> index & 1 else 1) for index, edge in enumerate(edges(n))
    }


def signings(n: int) -> Iterable[dict[tuple[int, int], int]]:
    for mask in range(1 << len(edges(n))):
        yield signing_from_mask(n, mask)


def spin(mask: int, vertex: int) -> int:
    return -1 if mask >> vertex & 1 else 1


def energy(n: int, coefficients: dict[tuple[int, int], int], spin_mask: int) -> int:
    return sum(
        value * spin(spin_mask, row) * spin(spin_mask, column)
        for (row, column), value in coefficients.items()
    )


def partition(
    n: int, coefficients: dict[tuple[int, int], int], temperature: float
) -> float:
    return sum(
        2.0 * math.cosh(temperature * energy(n, coefficients, spin_mask))
        for spin_mask in range(1 << n)
    )


def minimum_log_partition(n: int, temperature: float) -> float:
    return min(
        math.log(partition(n, coefficients, temperature))
        for coefficients in signings(n)
    )


def extend_signing(
    n: int,
    coefficients: dict[tuple[int, int], int],
    incident_signs: tuple[int, ...],
) -> dict[tuple[int, int], int]:
    if len(incident_signs) != n:
        raise ValueError("wrong number of incident signs")
    extended = dict(coefficients)
    for vertex, value in enumerate(incident_signs):
        extended[(vertex, n)] = value
    return extended


def check_cavity_identities() -> int:
    checked = 0
    for n in range(1, 5):
        for coefficients in signings(n):
            for temperature in (0.0, 0.17, 0.61):
                base_partition = partition(n, coefficients, temperature)
                extension_sum = 0.0
                for incident_signs in itertools.product((-1, 1), repeat=n):
                    extended = extend_signing(n, coefficients, incident_signs)
                    extension_sum += partition(n + 1, extended, temperature)
                    checked += 1
                extension_average = extension_sum / (1 << n)
                expected = 2.0 * base_partition * math.cosh(temperature) ** n
                if not math.isclose(
                    extension_average,
                    expected,
                    rel_tol=TOLERANCE,
                    abs_tol=TOLERANCE,
                ):
                    raise AssertionError(("cavity average mismatch", n, temperature))
    return checked


def check_minimum_partition_inequalities() -> int:
    values: dict[tuple[int, float], float] = {}
    temperatures = (0.0, 0.19, 0.53)
    for n in range(1, 6):
        for temperature in temperatures:
            values[(n, temperature)] = minimum_log_partition(n, temperature)

    checked = 0
    for temperature in temperatures:
        for n, k in ((1, 1), (1, 2), (2, 2), (2, 3)):
            combined = values[(n + k, temperature)]
            block_upper = (
                values[(n, temperature)]
                + values[(k, temperature)]
                - math.log(2.0)
                + n * k * math.log(math.cosh(temperature))
            )
            if combined > block_upper + TOLERANCE:
                raise AssertionError(("block inequality failed", n, k, temperature))
            checked += 1

        for n in range(1, 5):
            increment = values[(n + 1, temperature)] - values[(n, temperature)]
            if increment < math.log(2.0) - TOLERANCE:
                raise AssertionError(("cavity lower bound failed", n, temperature))
            if increment > (
                math.log(2.0) + n * math.log(math.cosh(temperature)) + TOLERANCE
            ):
                raise AssertionError(("cavity upper bound failed", n, temperature))
            checked += 1
    return checked


def matrix_from_signing(
    n: int, coefficients: dict[tuple[int, int], int]
) -> list[list[float]]:
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    for (row, column), value in coefficients.items():
        matrix[row][column] = float(value)
        matrix[column][row] = float(value)
    return matrix


def matrix_square(matrix: list[list[float]]) -> list[list[float]]:
    n = len(matrix)
    return [
        [
            sum(matrix[row][middle] * matrix[middle][column] for middle in range(n))
            for column in range(n)
        ]
        for row in range(n)
    ]


def determinant(matrix: list[list[float]]) -> float:
    work = [row[:] for row in matrix]
    result = 1.0
    n = len(work)
    for column in range(n):
        pivot_row = max(range(column, n), key=lambda row: abs(work[row][column]))
        pivot = work[pivot_row][column]
        if abs(pivot) < 1e-14:
            return 0.0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            result = -result
        pivot = work[column][column]
        result *= pivot
        for row in range(column + 1, n):
            multiplier = work[row][column] / pivot
            for target_column in range(column + 1, n):
                work[row][target_column] -= multiplier * work[column][target_column]
    return result


def determinant_floor(damping: float) -> float:
    if not 0.0 <= damping < 1.0:
        raise ValueError("damping must be in [0,1)")
    if damping <= 2.0 / 3.0:
        return (1.0 - damping / 2.0) ** 2
    return 2.0 * damping * (1.0 - damping)


def covariance_pair(
    matrix: list[list[float]], damping: float
) -> tuple[list[list[float]], list[list[float]]]:
    n = len(matrix)
    q = n - 1
    square = matrix_square(matrix)
    covariances = []
    for orientation in (1.0, -1.0):
        covariance = []
        for row in range(n):
            covariance_row = []
            for column in range(n):
                identity = 1.0 if row == column else 0.0
                base = 0.5 * (
                    identity
                    + square[row][column] / q
                    + orientation * 2.0 * matrix[row][column] / math.sqrt(q)
                )
                covariance_row.append((1.0 - damping) * identity + damping * base)
            covariance.append(covariance_row)
        covariances.append(covariance)
    return covariances[0], covariances[1]


def covariance_samples() -> Iterable[dict[tuple[int, int], int]]:
    for n in range(2, 6):
        yield from signings(n)
    generator = random.Random(DETERMINISTIC_SEED)
    n = 6
    edge_count = len(edges(n))
    for _ in range(128):
        yield signing_from_mask(n, generator.randrange(1 << edge_count))


def check_covariance_algebra() -> tuple[int, int, bool]:
    edge_checks = 0
    determinant_checks = 0
    false_strengthening_detected = False
    for coefficients in covariance_samples():
        n = max(max(edge) for edge in coefficients) + 1
        matrix = matrix_from_signing(n, coefficients)
        for damping in (0.0, 0.2, 0.67, 0.9):
            plus, minus = covariance_pair(matrix, damping)
            for index in range(n):
                if not math.isclose(plus[index][index], 1.0, abs_tol=TOLERANCE):
                    raise AssertionError("plus covariance diagonal corruption")
                if not math.isclose(minus[index][index], 1.0, abs_tol=TOLERANCE):
                    raise AssertionError("minus covariance diagonal corruption")

            q = n - 1
            for (row, column), value in coefficients.items():
                oriented_difference = value * (
                    math.asin(plus[row][column]) - math.asin(minus[row][column])
                )
                expected_floor = 2.0 * damping / math.sqrt(q)
                if oriented_difference < expected_floor - TOLERANCE:
                    raise AssertionError(
                        ("arcsine orientation failed", n, damping, row, column)
                    )
                edge_checks += 1

            determinant_product = determinant(plus) * determinant(minus)
            certified_floor = determinant_floor(damping) ** n
            if determinant_product < certified_floor - 2e-8:
                raise AssertionError(
                    (
                        "determinant floor failed",
                        n,
                        damping,
                        determinant_product,
                        certified_floor,
                    )
                )
            if determinant_product < (1.0 - damping**2) ** n - 1e-8:
                false_strengthening_detected = True
            determinant_checks += 1

    if not false_strengthening_detected:
        raise AssertionError("false conference-style determinant bound survived")
    return edge_checks, determinant_checks, false_strengthening_detected


def log_mean_negative_power(values: list[float], power: float) -> float:
    exponents = [-power * math.log(value) for value in values]
    maximum = max(exponents)
    return maximum + math.log(
        sum(math.exp(exponent - maximum) for exponent in exponents) / len(values)
    )


def check_negative_replica_sandwich() -> int:
    checks = 0
    for n in range(2, 6):
        coefficients_list = list(signings(n))
        signing_count = len(coefficients_list)
        expected_count = 1 << len(edges(n))
        if signing_count != expected_count:
            raise AssertionError("signing count corruption")
        for beta in (0.7, 2.0):
            temperature = beta / math.sqrt(n)
            partitions = [
                partition(n, coefficients, temperature)
                for coefficients in coefficients_list
            ]
            minimax = math.log(min(partitions)) / (beta * n)
            for power in (0.5, float(n)):
                soft_minimum = -log_mean_negative_power(partitions, power) / (
                    power * beta * n
                )
                error = len(edges(n)) * math.log(2.0) / (power * beta * n)
                if soft_minimum < minimax - TOLERANCE:
                    raise AssertionError("negative-replica lower sandwich failed")
                if soft_minimum > minimax + error + TOLERANCE:
                    raise AssertionError("negative-replica upper sandwich failed")
                checks += 1
    return checks


def sylvester(order: int) -> tuple[tuple[int, ...], ...]:
    if order == 1:
        return ((1,),)
    if order % 2:
        raise ValueError("Sylvester order must be a power of two")
    smaller = sylvester(order // 2)
    top = tuple(tuple(row) + tuple(row) for row in smaller)
    bottom = tuple(tuple(row) + tuple(-value for value in row) for row in smaller)
    return top + bottom


def check_hadamard_cavity_example() -> int:
    order = 16
    hadamard = sylvester(order)
    for row in range(order):
        for column in range(order):
            inner_product = sum(
                hadamard[index][row] * hadamard[index][column] for index in range(order)
            )
            expected = order if row == column else 0
            if inner_product != expected:
                raise AssertionError("Hadamard covariance corruption")

    bent = []
    for point in range(order):
        bits = tuple(point >> index & 1 for index in range(4))
        phase = (bits[0] & bits[1]) ^ (bits[2] & bits[3])
        bent.append(-1 if phase else 1)
    transform = [
        sum(row[column] * bent[column] for column in range(order)) for row in hadamard
    ]
    if {abs(value) for value in transform} != {4}:
        raise AssertionError("bent Walsh spectrum is not flat")

    checks = 0
    for temperature in (0.11, 0.37):
        hadamard_cost = (
            sum(math.cosh(temperature * value) for value in transform) / order
        )
        expected_hadamard = math.cosh(temperature * math.sqrt(order))
        uniform_cost = math.cosh(temperature) ** order
        if not math.isclose(
            hadamard_cost,
            expected_hadamard,
            rel_tol=TOLERANCE,
            abs_tol=TOLERANCE,
        ):
            raise AssertionError("Hadamard cavity equality failed")
        if uniform_cost <= 0.0:
            raise AssertionError("uniform cavity cost corruption")
        checks += 1
    return checks


def abstract_b(index: int) -> float:
    if index == 1:
        return 0.0
    return index**2 * (1.0 + 0.1 * math.sin(math.log(index)))


def check_abstract_countermodel() -> int:
    checks = 0
    for n in range(1, 301):
        increment = 0.1 * (abstract_b(n + 1) - abstract_b(n))
        if not 0.0 <= increment <= n + TOLERANCE:
            raise AssertionError(("abstract cavity coefficient failed", n, increment))
        for k in (1, 2, 7, 29):
            if abstract_b(n + k) + TOLERANCE < abstract_b(n) + abstract_b(k):
                raise AssertionError(("abstract superadditivity failed", n, k))
            checks += 1

    for damping_index in range(1, 1000):
        damping = damping_index / 1000.0
        if determinant_floor(damping) > math.exp(-damping) + TOLERANCE:
            raise AssertionError("abstract entropy-barrier comparison failed")
        checks += 1
    return checks


def verify_corruption_controls() -> None:
    zero_temperature_increment = minimum_log_partition(2, 0.0) - minimum_log_partition(
        1, 0.0
    )
    if zero_temperature_increment >= 2.0 * math.log(2.0) - TOLERANCE:
        raise AssertionError("incorrect cavity factor was not detected")

    for damping in (0.2, 0.8):
        minimizing_u = 0.0 if damping <= 2.0 / 3.0 else 3.0 - 2.0 / damping
        value = (
            1.0 - damping / 2.0 + damping * minimizing_u / 2.0
        ) ** 2 - damping**2 * minimizing_u
        if not math.isclose(
            value,
            determinant_floor(damping),
            rel_tol=TOLERANCE,
            abs_tol=TOLERANCE,
        ):
            raise AssertionError("determinant minimizer corruption was not detected")
        if value >= determinant_floor(damping) * (1.0 + 1e-6):
            raise AssertionError("inflated determinant floor was not detected")


def main() -> None:
    cavity_checks = check_cavity_identities()
    partition_checks = check_minimum_partition_inequalities()
    edge_checks, determinant_checks, false_bound = check_covariance_algebra()
    replica_checks = check_negative_replica_sandwich()
    hadamard_checks = check_hadamard_cavity_example()
    countermodel_checks = check_abstract_countermodel()
    verify_corruption_controls()
    print(f"cavity_extensions_checked={cavity_checks}")
    print(f"partition_inequalities_checked={partition_checks}")
    print(f"covariance_edges_checked={edge_checks}")
    print(f"determinant_products_checked={determinant_checks}")
    print(f"negative_replica_checks={replica_checks}")
    print(f"hadamard_cavity_checks={hadamard_checks}")
    print(f"abstract_countermodel_checks={countermodel_checks}")
    print(f"false_conference_determinant_bound_detected={str(false_bound).upper()}")
    print(f"deterministic_seed={DETERMINISTIC_SEED}")
    print("corruption_controls=PASSED")


if __name__ == "__main__":
    main()

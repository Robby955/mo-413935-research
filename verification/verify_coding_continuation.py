#!/usr/bin/env python3
"""Exact checks for the continued coding and multiaffine results."""

from __future__ import annotations

import itertools
from collections import defaultdict
from fractions import Fraction


def complete_edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(itertools.combinations(range(n), 2))


def toggle_boundary(boundary: frozenset[int], edge: tuple[int, int]) -> frozenset[int]:
    result = set(boundary)
    row, column = edge
    if row != column:
        for vertex in (row, column):
            if vertex in result:
                result.remove(vertex)
            else:
                result.add(vertex)
    return frozenset(result)


def transfer_sectors(
    graph_edges: tuple[tuple[int, int], ...], weights: tuple[Fraction, ...]
) -> dict[tuple[frozenset[int], int], Fraction]:
    if len(graph_edges) != len(weights):
        raise ValueError("edge and weight counts differ")
    sectors: dict[tuple[frozenset[int], int], Fraction] = defaultdict(Fraction)
    for mask in range(1 << len(graph_edges)):
        boundary: frozenset[int] = frozenset()
        parity = 0
        term = Fraction(1)
        for index, edge in enumerate(graph_edges):
            if mask >> index & 1:
                boundary = toggle_boundary(boundary, edge)
                parity ^= 1
                term *= weights[index]
        sectors[(boundary, parity)] += term
    return dict(sectors)


def p_and_r(
    graph_edges: tuple[tuple[int, int], ...], weights: tuple[Fraction, ...]
) -> tuple[Fraction, Fraction]:
    sectors = transfer_sectors(graph_edges, weights)
    return (
        sectors.get((frozenset(), 0), Fraction(0)),
        sectors.get((frozenset(), 1), Fraction(0)),
    )


def contract_graph(
    vertex_count: int,
    graph_edges: tuple[tuple[int, int], ...],
    weights: tuple[Fraction, ...],
    edge_index: int,
) -> tuple[int, tuple[tuple[int, int], ...], tuple[Fraction, ...]]:
    row, column = graph_edges[edge_index]
    if row == column:
        raise ValueError("only nonloop contraction is supported")
    representative = {
        vertex: (row if vertex == column else vertex) for vertex in range(vertex_count)
    }
    surviving_vertices = sorted(set(representative.values()))
    compressed = {vertex: index for index, vertex in enumerate(surviving_vertices)}
    contracted_edges = []
    contracted_weights = []
    for index, edge in enumerate(graph_edges):
        if index == edge_index:
            continue
        contracted_edges.append(
            (
                compressed[representative[edge[0]]],
                compressed[representative[edge[1]]],
            )
        )
        contracted_weights.append(weights[index])
    return (
        len(surviving_vertices),
        tuple(contracted_edges),
        tuple(contracted_weights),
    )


def check_deletion_contraction() -> tuple[int, bool]:
    checks = 0
    omitted_term_detected = False
    for vertex_count in (3, 4):
        graph_edges = complete_edges(vertex_count)
        weights = tuple(
            Fraction(((-1) ** index) * (index % 3 + 1), index % 4 + 5)
            for index in range(len(graph_edges))
        )
        p_graph, r_graph = p_and_r(graph_edges, weights)
        for edge_index, edge_weight in enumerate(weights):
            deleted_edges = graph_edges[:edge_index] + graph_edges[edge_index + 1 :]
            deleted_weights = weights[:edge_index] + weights[edge_index + 1 :]
            p_deleted, r_deleted = p_and_r(deleted_edges, deleted_weights)
            _, contracted_edges, contracted_weights = contract_graph(
                vertex_count, graph_edges, weights, edge_index
            )
            p_contracted, r_contracted = p_and_r(contracted_edges, contracted_weights)
            expected_p = p_deleted + edge_weight * (r_contracted - r_deleted)
            expected_r = r_deleted + edge_weight * (p_contracted - p_deleted)
            if p_graph != expected_p or r_graph != expected_r:
                raise AssertionError(
                    ("deletion-contraction failed", vertex_count, edge_index)
                )
            if p_graph != p_deleted + edge_weight * r_contracted:
                omitted_term_detected = True
            checks += 1
    if not omitted_term_detected:
        raise AssertionError("omitted deletion term was not detected")
    return checks, omitted_term_detected


def check_transfer_updates() -> int:
    graph_edges = complete_edges(4)
    weights = tuple(Fraction(index + 1, index + 7) for index in range(len(graph_edges)))
    recursive: dict[tuple[frozenset[int], int], Fraction] = {
        (frozenset(), 0): Fraction(1)
    }
    checks = 0
    for edge, weight in zip(graph_edges, weights, strict=True):
        keys = {
            (boundary, parity)
            for boundary_size in range(5)
            for boundary in map(
                frozenset, itertools.combinations(range(4), boundary_size)
            )
            for parity in (0, 1)
        }
        updated = {}
        for boundary, parity in keys:
            updated[(boundary, parity)] = recursive.get(
                (boundary, parity), Fraction(0)
            ) + weight * recursive.get(
                (toggle_boundary(boundary, edge), parity ^ 1), Fraction(0)
            )
        recursive = updated
        direct = transfer_sectors(
            graph_edges[: checks // 32 + 1], weights[: checks // 32 + 1]
        )
        for key in keys:
            if recursive.get(key, Fraction(0)) != direct.get(key, Fraction(0)):
                raise AssertionError(("transfer update failed", edge, key))
            checks += 1
    return checks


def augmented_cut_code(n: int) -> set[int]:
    graph_edges = complete_edges(n)
    words = set()
    for vertex_mask in range(1 << n):
        for complement in (0, 1):
            word = 0
            for index, (row, column) in enumerate(graph_edges):
                bit = (
                    (vertex_mask >> row & 1) ^ (vertex_mask >> column & 1) ^ complement
                )
                word |= bit << index
            words.add(word)
    return words


def coset_probabilities(n: int, bias: Fraction) -> tuple[dict[int, Fraction], set[int]]:
    graph_edges = complete_edges(n)
    code = augmented_cut_code(n)
    probabilities: dict[int, Fraction] = defaultdict(Fraction)
    for word in range(1 << len(graph_edges)):
        ones = word.bit_count()
        probability = ((1 - bias) / 2) ** ones * ((1 + bias) / 2) ** (
            len(graph_edges) - ones
        )
        representative = min(word ^ codeword for codeword in code)
        probabilities[representative] += probability
    return dict(probabilities), code


def signed_polynomial(n: int, word: int, bias: Fraction) -> Fraction:
    graph_edges = complete_edges(n)
    weights = tuple(
        (-bias if word >> index & 1 else bias) for index in range(len(graph_edges))
    )
    return p_and_r(graph_edges, weights)[0]


def check_noise_cosets() -> tuple[int, bool]:
    checks = 0
    strict_monotonicity_detected = False
    for n in (3, 4):
        minima = []
        for bias in (Fraction(1, 4), Fraction(1, 2)):
            probabilities, code = coset_probabilities(n, bias)
            graph_edges = complete_edges(n)
            scale = Fraction(2**n, 2 ** len(graph_edges))
            polynomial_values = []
            for word in range(1 << len(graph_edges)):
                representative = min(word ^ codeword for codeword in code)
                polynomial = signed_polynomial(n, word, bias)
                if probabilities[representative] != scale * polynomial:
                    raise AssertionError(("coset-noise identity failed", n, bias, word))
                polynomial_values.append(polynomial)
                checks += 1
            minima.append(min(polynomial_values))
        if minima[0] < minima[1]:
            raise AssertionError(("noise monotonicity failed", n, minima))
        if minima[0] > minima[1]:
            strict_monotonicity_detected = True
    if not strict_monotonicity_detected:
        raise AssertionError("strict noise monotonicity example was not found")
    return checks, strict_monotonicity_detected


def check_multiaffine_box() -> tuple[int, Fraction, Fraction]:
    graph_edges = complete_edges(4)
    radius = Fraction(1, 3)
    endpoint_values = []
    for signs in itertools.product((-1, 1), repeat=len(graph_edges)):
        weights = tuple(radius * sign for sign in signs)
        endpoint_values.append(p_and_r(graph_edges, weights)[0])
    endpoint_minimum = min(endpoint_values)

    grid_minimum = None
    grid_points = 0
    for values in itertools.product(
        (-radius, Fraction(0), radius), repeat=len(graph_edges)
    ):
        polynomial = p_and_r(graph_edges, values)[0]
        grid_minimum = (
            polynomial if grid_minimum is None else min(grid_minimum, polynomial)
        )
        grid_points += 1
    if grid_minimum != endpoint_minimum:
        raise AssertionError(
            ("box endpoint property failed", grid_minimum, endpoint_minimum)
        )

    triangle_edges = complete_edges(3)
    triangle_weights = (radius, radius, radius)
    if p_and_r(triangle_edges, triangle_weights) != (Fraction(1), radius**3):
        raise AssertionError("triangle P/R corruption")

    plus_star = triangle_edges + ((0, 3), (1, 3), (2, 3))
    plus_weights = triangle_weights + (radius, radius, radius)
    mixed_weights = triangle_weights + (-radius, radius, radius)
    plus_value = p_and_r(plus_star, plus_weights)[0]
    mixed_value = p_and_r(plus_star, mixed_weights)[0]
    if plus_value != 1 + 3 * radius**4 or mixed_value != 1 - radius**4:
        raise AssertionError("K3 vertex-extension counterexample failed")
    return grid_points, plus_value, mixed_value


def main() -> None:
    deletion_checks, deletion_corruption = check_deletion_contraction()
    transfer_checks = check_transfer_updates()
    noise_checks, strict_noise = check_noise_cosets()
    grid_points, plus_value, mixed_value = check_multiaffine_box()
    print(f"deletion_contraction_checks={deletion_checks}")
    print(f"transfer_sector_checks={transfer_checks}")
    print(f"noise_cosets_checked={noise_checks}")
    print(f"multiaffine_grid_points={grid_points}")
    print(f"k3_extension_values={plus_value},{mixed_value}")
    print(f"strict_noise_monotonicity_detected={str(strict_noise).upper()}")
    print(f"omitted_deletion_term_detected={str(deletion_corruption).upper()}")
    print("corruption_controls=PASSED")


if __name__ == "__main__":
    main()

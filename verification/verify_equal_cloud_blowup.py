#!/usr/bin/env python3
"""Exact tiny checks for fixed-half equal-cloud blow-up obstructions."""

from __future__ import annotations

import hashlib
import itertools
from collections.abc import Iterable

if not __debug__:
    raise RuntimeError("verification requires Python assertions")


Edge = tuple[int, int]


def edges_of_order(order: int) -> tuple[Edge, ...]:
    return tuple(itertools.combinations(range(order), 2))


def doubled_cut_discrepancy(order: int, graph: frozenset[Edge]) -> int:
    """Return max_S |2 e(S,S^c) - |S||S^c||."""

    best = 0
    for mask in range(1 << order):
        size = mask.bit_count()
        crossing = sum(((mask >> i) ^ (mask >> j)) & 1 for i, j in graph)
        best = max(best, abs(2 * crossing - size * (order - size)))
    return best


def cut_deviation(
    order: int, graph: frozenset[Edge], mask: int
) -> int:
    """Return twice the signed density-1/2 deviation of one cut."""

    size = mask.bit_count()
    crossing = sum(((mask >> i) ^ (mask >> j)) & 1 for i, j in graph)
    return 2 * crossing - size * (order - size)


def fixed_half_graphs(order: int) -> Iterable[frozenset[Edge]]:
    edges = edges_of_order(order)
    for selected in itertools.combinations(edges, len(edges) // 2):
        yield frozenset(selected)


def equal_cloud_blowup(
    order: int, cloud_size: int, graph: frozenset[Edge]
) -> frozenset[Edge]:
    """Construct a fixed-half complete/empty equal-cloud blow-up."""

    total_order = order * cloud_size
    target = len(edges_of_order(total_order)) // 2
    result: set[Edge] = set()
    for left, right in graph:
        for a in range(cloud_size):
            for b in range(cloud_size):
                result.add((left * cloud_size + a, right * cloud_size + b))

    internal_candidates = tuple(
        (cloud * cloud_size + a, cloud * cloud_size + b)
        for cloud in range(order)
        for a, b in itertools.combinations(range(cloud_size), 2)
    )
    required_internal = target - len(result)
    if not 0 <= required_internal <= len(internal_candidates):
        raise ValueError("fixed-half internal completion is infeasible")
    result.update(internal_candidates[:required_internal])
    if len(result) != target:
        raise AssertionError("fixed-half completion count mismatch")
    return frozenset(result)


def seidel_switch(
    order: int, graph: frozenset[Edge], switching_mask: int
) -> frozenset[Edge]:
    result = set(graph)
    for edge in edges_of_order(order):
        left, right = edge
        if ((switching_mask >> left) ^ (switching_mask >> right)) & 1:
            if edge in result:
                result.remove(edge)
            else:
                result.add(edge)
    return frozenset(result)


def rebalance_to_fixed_half(
    order: int, graph: frozenset[Edge]
) -> tuple[frozenset[Edge], int]:
    all_edges = edges_of_order(order)
    target = len(all_edges) // 2
    result = set(graph)
    edits = abs(len(result) - target)
    if len(result) > target:
        for edge in sorted(result)[:edits]:
            result.remove(edge)
    else:
        missing = (edge for edge in all_edges if edge not in result)
        result.update(itertools.islice(missing, edits))
    if len(result) != target:
        raise AssertionError("rebalancing failed")
    return frozenset(result), edits


def lifted_mask(base_mask: int, order: int, cloud_size: int) -> int:
    result = 0
    for cloud in range(order):
        if (base_mask >> cloud) & 1:
            for coordinate in range(cloud_size):
                result |= 1 << (cloud * cloud_size + coordinate)
    return result


def verify_feasibility_exception() -> None:
    # For n=2, a fixed-half base has either zero or one complete cross block.
    # The internal capacity 2*binom(k,2) leaves the global half-edge target in
    # the gap for every k>1.
    for cloud_size in range(2, 7):
        base = next(iter(fixed_half_graphs(2)))
        try:
            equal_cloud_blowup(2, cloud_size, base)
        except ValueError:
            continue
        raise AssertionError(("n=2 exception was not detected", cloud_size))


def main() -> None:
    verify_feasibility_exception()
    cases = ((3, 2), (3, 3), (4, 2))
    blowups_checked = 0
    cloud_cuts_checked = 0
    switches_checked = 0
    transcript: list[str] = []
    corruption_detected = False

    for order, cloud_size in cases:
        for graph in fixed_half_graphs(order):
            base_discrepancy = doubled_cut_discrepancy(order, graph)
            blowup = equal_cloud_blowup(order, cloud_size, graph)
            total_order = order * cloud_size
            total_edges = len(edges_of_order(total_order))
            if len(blowup) != total_edges // 2:
                raise AssertionError("blow-up is not fixed-half")

            cloud_maximum = 0
            for mask in range(1 << order):
                base_deviation = cut_deviation(order, graph, mask)
                blowup_deviation = cut_deviation(
                    total_order,
                    blowup,
                    lifted_mask(mask, order, cloud_size),
                )
                if blowup_deviation != cloud_size**2 * base_deviation:
                    raise AssertionError("cloud-union scaling identity failed")
                cloud_maximum = max(cloud_maximum, abs(blowup_deviation))
                cloud_cuts_checked += 1
            if cloud_maximum != cloud_size**2 * base_discrepancy:
                raise AssertionError("cloud-union maximum scaling failed")

            if not corruption_detected and base_discrepancy > 0:
                corrupted = set(blowup)
                first_cross = (0, cloud_size)
                if first_cross in corrupted:
                    corrupted.remove(first_cross)
                else:
                    corrupted.add(first_cross)
                corrupted_deviation = cut_deviation(
                    total_order,
                    frozenset(corrupted),
                    lifted_mask(1, order, cloud_size),
                )
                correct_deviation = cut_deviation(
                    total_order, blowup, lifted_mask(1, order, cloud_size)
                )
                corruption_detected = corrupted_deviation != correct_deviation

            target_parity = total_edges & 1
            for switching_mask in range(1 << total_order):
                switched = seidel_switch(total_order, blowup, switching_mask)
                rebalanced, edits = rebalance_to_fixed_half(total_order, switched)
                final_discrepancy = doubled_cut_discrepancy(
                    total_order, rebalanced
                )
                # This is twice
                # d(L') >= k^2 d(G) - (r+t_N)/2.
                if final_discrepancy < (
                    cloud_size**2 * base_discrepancy - edits - target_parity
                ):
                    raise AssertionError("switch-edit lower bound failed")
                transcript.append(
                    f"{order},{cloud_size},{len(graph)},{switching_mask},"
                    f"{edits},{final_discrepancy}"
                )
                switches_checked += 1
            blowups_checked += 1

    if not corruption_detected:
        raise AssertionError("cross-edge corruption was not detected")
    digest = hashlib.sha256("\n".join(transcript).encode()).hexdigest()
    expected_digest = (
        "3e0d47a7b780432841530c2231881d0bddc7169c0abc9f5238bab1e9347ed535"
    )
    if digest != expected_digest:
        raise AssertionError(("switch-profile digest mismatch", digest))
    print("orders_and_clouds=3x2,3x3,4x2")
    print(f"fixed_half_blowups_checked={blowups_checked}")
    print(f"cloud_union_cuts_checked={cloud_cuts_checked}")
    print(f"switch_rebalance_checks={switches_checked}")
    print("n2_infeasible_cloud_sizes=2..6")
    print(f"switch_profile_sha256={digest}")
    print("corruption_controls=cross_edge,n2_feasibility")
    print("equal_cloud_blowup_verification=PASSED")


if __name__ == "__main__":
    main()

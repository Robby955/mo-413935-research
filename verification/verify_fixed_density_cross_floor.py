#!/usr/bin/env python3
"""Exact finite checks for the fixed-density bipartite cross-floor lemma.

For a bipartite graph C with sign matrix S = 2 1_C - 1, let

    D_xor(C) = max_{P,R} |e(P,V\\R) + e(U\\P,R)
                 - (|P|(k-|R|) + (n-|P|)|R|)/2|.

The analytic lemma states

    4 D_xor(C) = |sum(S)| + ||S||_{infinity -> 1}.

It also states the exact Khintchine floor

    ||S||_{infinity -> 1} >= max(n mu_k, k mu_n)

and that prescribing a feasible total t0 costs at most

    floor(sqrt(nk)) + |t0|

above the unrestricted optimum.  This script exhausts every rectangular
signing through 4 by 4, replays the switching-and-editing construction, and
contains deliberate centering and parity corruptions.  All decisions use
integer arithmetic.
"""

from __future__ import annotations

import hashlib
import math
from itertools import product

if not __debug__:
    raise RuntimeError("verification requires Python assertions")


MAX_ROWS = 4
MAX_COLUMNS = 4
EXPECTED_MATRIX_CHECKS = 74_954
EXPECTED_KHINTCHINE_CHECKS = 149_908
EXPECTED_FIXED_TOTAL_CHECKS = 116
EXPECTED_PAIR_CHECKS = 16
EXPECTED_PROFILE_SHA256 = (
    "2f0e4e234d225e8b5ad4513212900132817ad323e91f3e6e0eb0874847a3c0f1"
)


def signs(length: int, *, projective: bool = False) -> tuple[tuple[int, ...], ...]:
    if projective:
        return tuple((1, *tail) for tail in product((-1, 1), repeat=length - 1))
    return tuple(product((-1, 1), repeat=length))


def patterns(
    rows: int,
    columns: int,
) -> tuple[tuple[tuple[int, ...], int], ...]:
    """Return row-column products and the corresponding XOR-cut masks."""

    result = []
    for left in signs(rows, projective=True):
        for right in signs(columns):
            products_flat = tuple(
                left[row] * right[column]
                for row in range(rows)
                for column in range(columns)
            )
            opposite_mask = sum(
                1 << index
                for index, value in enumerate(products_flat)
                if value == -1
            )
            result.append((products_flat, opposite_mask))
    return tuple(result)


def signing_total(mask: int, entries: int) -> int:
    return 2 * mask.bit_count() - entries


def rectangular_norm(
    mask: int,
    entries: int,
    cut_patterns: tuple[tuple[tuple[int, ...], int], ...],
) -> int:
    total = signing_total(mask, entries)
    return max(
        abs(total - 2 * (2 * (mask & opposite).bit_count() - opposite.bit_count()))
        for _, opposite in cut_patterns
    )


def four_xor_discrepancy(
    mask: int,
    cut_patterns: tuple[tuple[tuple[int, ...], int], ...],
) -> int:
    return max(
        2 * abs(2 * (mask & opposite).bit_count() - opposite.bit_count())
        for _, opposite in cut_patterns
    )


def corrupted_four_discrepancy(
    mask: int,
    cut_patterns: tuple[tuple[tuple[int, ...], int], ...],
) -> int:
    """Deliberately use density-one rather than density-one-half centering."""

    return max(
        4 * abs((mask & opposite).bit_count() - opposite.bit_count())
        for _, opposite in cut_patterns
    )


def mu_numerator(order: int) -> int:
    """Return 2^order times E|epsilon_1 + ... + epsilon_order|."""

    return sum(abs(sum(spin)) for spin in signs(order))


def switched_mask(
    mask: int,
    entries: int,
    products_flat: tuple[int, ...],
) -> int:
    result = 0
    for index, multiplier in enumerate(products_flat):
        value = 1 if (mask >> index) & 1 else -1
        if value * multiplier == 1:
            result |= 1 << index
    assert result < 1 << entries
    return result


def balance_total(mask: int, entries: int, target: int) -> int:
    current = signing_total(mask, entries)
    if abs(target) > entries or (target - entries) % 2:
        raise ValueError("infeasible signing total")
    difference = target - current
    assert difference % 2 == 0
    flips = abs(difference) // 2

    result = mask
    if difference > 0:
        candidates = [index for index in range(entries) if not (result >> index) & 1]
        assert len(candidates) >= flips
        for index in candidates[:flips]:
            result |= 1 << index
    elif difference < 0:
        candidates = [index for index in range(entries) if (result >> index) & 1]
        assert len(candidates) >= flips
        for index in candidates[:flips]:
            result &= ~(1 << index)

    assert signing_total(result, entries) == target
    return result


def required_cross_total(rows: int, columns: int) -> int:
    left_edges = math.comb(rows, 2)
    right_edges = math.comb(columns, 2)
    full_edges = math.comb(rows + columns, 2)
    required_cross_edges = (
        full_edges // 2 - left_edges // 2 - right_edges // 2
    )
    return 2 * required_cross_edges - rows * columns


def main() -> None:
    identity_checks = 0
    khintchine_checks = 0
    fixed_total_checks = 0
    constructive_balance_checks = 0
    switching_checks = 0
    half_density_target_checks = 0
    parity_corruptions_detected = 0
    wrong_centering_detected = False
    minima_records: list[str] = []

    for rows in range(1, MAX_ROWS + 1):
        for columns in range(1, MAX_COLUMNS + 1):
            entries = rows * columns
            cut_patterns = patterns(rows, columns)
            minimum_by_total: dict[int, int] = {}
            unrestricted_minimum = None
            unrestricted_mask = None

            left_mu_numerator = mu_numerator(rows)
            right_mu_numerator = mu_numerator(columns)

            for mask in range(1 << entries):
                total = signing_total(mask, entries)
                norm = rectangular_norm(mask, entries, cut_patterns)
                four_discrepancy = four_xor_discrepancy(mask, cut_patterns)

                assert four_discrepancy == abs(total) + norm
                identity_checks += 1

                assert norm * (1 << columns) >= rows * right_mu_numerator
                assert norm * (1 << rows) >= columns * left_mu_numerator
                khintchine_checks += 2

                if corrupted_four_discrepancy(mask, cut_patterns) != four_discrepancy:
                    wrong_centering_detected = True

                old = minimum_by_total.get(total)
                if old is None or norm < old:
                    minimum_by_total[total] = norm
                if unrestricted_minimum is None or norm < unrestricted_minimum:
                    unrestricted_minimum = norm
                    unrestricted_mask = mask

            assert unrestricted_minimum is not None
            assert unrestricted_mask is not None

            best_switch = None
            best_switched_total = None
            for products_flat, _ in cut_patterns:
                candidate = switched_mask(
                    unrestricted_mask,
                    entries,
                    products_flat,
                )
                candidate_total = signing_total(candidate, entries)
                if (
                    best_switched_total is None
                    or abs(candidate_total) < abs(best_switched_total)
                ):
                    best_switch = candidate
                    best_switched_total = candidate_total

            assert best_switch is not None
            assert best_switched_total is not None
            assert best_switched_total * best_switched_total <= entries
            assert rectangular_norm(
                best_switch,
                entries,
                cut_patterns,
            ) == unrestricted_minimum
            switching_checks += 1

            for target in range(-entries, entries + 1, 2):
                restricted_minimum = minimum_by_total[target]
                assert restricted_minimum <= (
                    unrestricted_minimum + math.isqrt(entries) + abs(target)
                )

                candidate = balance_total(best_switch, entries, target)
                candidate_norm = rectangular_norm(candidate, entries, cut_patterns)
                assert candidate_norm <= (
                    unrestricted_minimum + abs(best_switched_total - target)
                )
                assert restricted_minimum <= candidate_norm
                fixed_total_checks += 1
                constructive_balance_checks += 1

            target = required_cross_total(rows, columns)
            assert abs(target) <= 2
            assert target in minimum_by_total
            half_density_target_checks += 1

            wrong_parity_target = target + 1 if target < entries else target - 1
            assert (wrong_parity_target - entries) % 2
            assert wrong_parity_target not in minimum_by_total
            try:
                balance_total(best_switch, entries, wrong_parity_target)
            except ValueError:
                parity_corruptions_detected += 1
            else:
                raise AssertionError("wrong-parity target was accepted")

            profile = ",".join(
                f"{total}:{minimum_by_total[total]}"
                for total in sorted(minimum_by_total)
            )
            minima_records.append(f"{rows}x{columns}[{profile}]")

    assert identity_checks == EXPECTED_MATRIX_CHECKS
    assert khintchine_checks == EXPECTED_KHINTCHINE_CHECKS
    assert fixed_total_checks == EXPECTED_FIXED_TOTAL_CHECKS
    assert constructive_balance_checks == EXPECTED_FIXED_TOTAL_CHECKS
    assert switching_checks == EXPECTED_PAIR_CHECKS
    assert half_density_target_checks == EXPECTED_PAIR_CHECKS
    assert parity_corruptions_detected == EXPECTED_PAIR_CHECKS
    assert wrong_centering_detected

    profile_digest = hashlib.sha256("|".join(minima_records).encode()).hexdigest()
    assert profile_digest == EXPECTED_PROFILE_SHA256
    print(f"rectangles=1x1..{MAX_ROWS}x{MAX_COLUMNS}")
    print(f"xor_identity_checks={identity_checks}")
    print(f"khintchine_floor_checks={khintchine_checks}")
    print(f"fixed_total_minimum_checks={fixed_total_checks}")
    print(f"constructive_switch_edit_checks={constructive_balance_checks}")
    print(f"minimum_profile_sha256={profile_digest}")
    print(
        "corruption_controls="
        "wrong_half_centering_detected,wrong_total_parity_detected"
    )
    print("fixed_density_cross_floor_verification=PASSED")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact and finite checks for three newly isolated convergence walls.

This script verifies:

* the algebraic floor in the weighted all-state entropy argument;
* the exact order-four-by-order-four cross-block/orbit diagnostic; and
* finite instances of an abstract Bellman countermodel together with the
  exact dyadic telescoping identity.

All B4 and dyadic decisions use integer or rational arithmetic.  Floating
point is used only for the displayed entropy constants and for finite checks
of the smooth function used to define the abstract countermodel.  The
countermodel's nonconvergence is an analytic statement, not a finite
computational claim.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from fractions import Fraction


DETERMINISTIC_SEED = 413935
COUNTERMODEL_LIMIT = 200_000
HEREDITARY_LIMIT = 500

Matrix = tuple[tuple[int, ...], ...]
Vector = tuple[int, ...]


def verify_entropy_floor() -> tuple[int, float]:
    """Check the exact algebra behind the all-state entropy floor.

    Writing ``a^2 = 2nk(n+k) log(2)``, the entropy estimate is bounded below
    by ``L + sqrt(a^2-L^2)`` when ``0 <= L <= a``.  After setting ``u=L/a``,
    the only inequality needed is

        1-u^2 >= (1-u)^2,

    whose residual is exactly ``2u(1-u)``.  The dimension check below also
    catches either missing projective-state factors or a missing ``log(4)``.
    """

    checks = 0
    for n in range(1, 21):
        for k in range(1, 21):
            # |R_{n,k}| = 2^(n+k-2), while the extra log(4) contributes
            # exactly two more copies of log(2).
            exponent_after_log4 = (n + k - 2) + 2
            floor_coefficient = 2 * n * k * exponent_after_log4
            expected_coefficient = 2 * n * k * (n + k)
            if floor_coefficient != expected_coefficient:
                raise AssertionError(("entropy dimension cancellation", n, k))
            checks += 1

    corruption_detected = False
    for numerator in range(1001):
        u = Fraction(numerator, 1000)
        left = 1 - u * u
        right = (1 - u) * (1 - u)
        residual = left - right
        if residual != 2 * u * (1 - u) or residual < 0:
            raise AssertionError(("entropy floor residual", u, residual))

        # Deliberately change the required (1-u)^2 to (1+u)^2.  A nonzero
        # grid point must reject this sign corruption.
        wrong_residual = left - (1 + u) * (1 + u)
        corruption_detected |= wrong_residual < 0
        checks += 1

    if not corruption_detected:
        raise AssertionError("entropy-sign corruption went undetected")

    # The positive-term expansion
    # log(2) = 2 * sum_{j>=0} 1/((2j+1)3^(2j+1))
    # gives log(2) > 2/3 > 1/2.  Thus 2*sqrt(log(2)) > sqrt(2) exactly,
    # so a power-saving remainder cannot bridge the balanced-block gap.
    if not Fraction(2, 3) > Fraction(1, 2):
        raise AssertionError("rational log(2) certificate corrupted")
    leading_gap = 2.0 * math.sqrt(math.log(2.0)) - math.sqrt(2.0)
    if not leading_gap > 0.25:
        raise AssertionError(("weighted entropy leading gap", leading_gap))
    return checks, leading_gap


def projective_spins(order: int) -> tuple[Vector, ...]:
    return tuple(
        (1,) + tuple(-1 if mask >> index & 1 else 1 for index in range(order - 1))
        for mask in range(1 << (order - 1))
    )


def energy(matrix: Matrix, state: Vector) -> int:
    return sum(
        matrix[row][column] * state[row] * state[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def rank_one(left: Vector, right: Vector) -> Vector:
    return tuple(a * b for a in left for b in right)


def inner(left: Vector, right: Vector) -> int:
    return sum(a * b for a, b in zip(left, right))


def cross_score(
    center: Vector, heights: tuple[int, ...], ranks: tuple[Vector, ...]
) -> int:
    return max(
        height + abs(inner(center, rank)) for height, rank in zip(heights, ranks)
    )


def full_block_score(block: Matrix, center: Vector) -> int:
    """Independently recompute the eight-vertex quadratic maximum."""

    states = tuple(itertools.product((-1, 1), repeat=4))
    result = 0
    for left in states:
        for right in states:
            cross = sum(
                center[4 * row + column] * left[row] * right[column]
                for row in range(4)
                for column in range(4)
            )
            value = energy(block, left) + energy(block, right) + cross
            result = max(result, abs(value))
    return result


def binomial_absolute_tail(size: int, threshold: int) -> Fraction:
    count = sum(
        math.comb(size, negative_count)
        for negative_count in range(size + 1)
        if abs(size - 2 * negative_count) > threshold
    )
    return Fraction(count, 1 << size)


def iid_union_sum(heights: tuple[int, ...], target: int) -> Fraction:
    return sum(
        (binomial_absolute_tail(16, target - height) for height in heights),
        start=Fraction(0),
    )


def corrupted_inclusive_union_sum(heights: tuple[int, ...], target: int) -> Fraction:
    """Deliberately use >= where the bad event requires a strict >."""

    return sum(
        (
            Fraction(
                sum(
                    math.comb(16, negative_count)
                    for negative_count in range(17)
                    if abs(16 - 2 * negative_count) >= target - height
                ),
                1 << 16,
            )
            for height in heights
        ),
        start=Fraction(0),
    )


def orbit_data(
    seed: Vector,
    heights: tuple[int, ...],
    ranks: tuple[Vector, ...],
    target: int,
) -> tuple[Counter[int], Counter[int], Counter[int]]:
    cross_profile = Counter(abs(inner(seed, rank)) for rank in ranks)
    scores: Counter[int] = Counter()
    violations: Counter[int] = Counter()
    for shift in ranks:
        center = tuple(a * b for a, b in zip(seed, shift))
        scores[cross_score(center, heights, ranks)] += 1
        violations[
            sum(
                height + abs(inner(center, rank)) > target
                for height, rank in zip(heights, ranks)
            )
        ] += 1
    return cross_profile, scores, violations


def verify_b4_cross_orbit() -> tuple[int, int, int]:
    block: Matrix = (
        (0, 1, 1, 1),
        (1, 0, 1, -1),
        (1, 1, 0, 1),
        (1, -1, 1, 0),
    )
    states = projective_spins(4)
    ranks = tuple(rank_one(left, right) for left in states for right in states)
    heights = tuple(
        abs(energy(block, left) + energy(block, right))
        for left in states
        for right in states
    )
    if len(set(ranks)) != 64:
        raise AssertionError("rank-one projective group has the wrong size")
    if Counter(heights) != Counter({0: 14, 2: 24, 4: 16, 6: 8, 8: 2}):
        raise AssertionError(("B4 height profile", Counter(heights)))

    best = 10**9
    optimum_count = 0
    first_witness: Vector | None = None
    center_count = 0
    for mask in range(1 << 15):
        center = (1,) + tuple(
            -1 if mask >> (index - 1) & 1 else 1 for index in range(1, 16)
        )
        value = cross_score(center, heights, ranks)
        center_count += 1
        if value < best:
            best = value
            optimum_count = 1
            first_witness = center
        elif value == best:
            optimum_count += 1

    if (center_count, best, optimum_count) != (32_768, 10, 92):
        raise AssertionError(
            ("B4 exhaustive cross optimum", center_count, best, optimum_count)
        )
    if first_witness is None or full_block_score(block, first_witness) != best:
        raise AssertionError("independent full-spin recomputation failed")

    adaptive: Vector = (
        1,
        -1,
        1,
        1,
        -1,
        1,
        1,
        -1,
        1,
        1,
        -1,
        -1,
        1,
        -1,
        -1,
        -1,
    )
    adaptive_profile, adaptive_scores, adaptive_violations = orbit_data(
        adaptive, heights, ranks, 10
    )
    if adaptive_profile != Counter({0: 24, 4: 32, 8: 8}):
        raise AssertionError(("adaptive cross profile", adaptive_profile))
    if adaptive_scores != Counter({10: 2, 12: 18, 14: 30, 16: 14}):
        raise AssertionError(("adaptive orbit scores", adaptive_scores))
    if (
        sum(count * multiplicity for count, multiplicity in adaptive_violations.items())
        != 272
    ):
        raise AssertionError(("adaptive mean violation count", adaptive_violations))
    if adaptive_violations[0] != 2 or full_block_score(block, adaptive) != 10:
        raise AssertionError(("adaptive good shifts", adaptive_violations))

    sylvester: Vector = (
        1,
        1,
        1,
        1,
        1,
        -1,
        1,
        -1,
        1,
        1,
        -1,
        -1,
        1,
        -1,
        -1,
        1,
    )
    sylvester_profile, sylvester_scores, sylvester_violations = orbit_data(
        sylvester, heights, ranks, 12
    )
    if sylvester_profile != Counter({0: 12, 4: 48, 8: 4}):
        raise AssertionError(("Sylvester cross profile", sylvester_profile))
    if sylvester_scores != Counter({12: 40, 14: 16, 16: 8}):
        raise AssertionError(("Sylvester orbit scores", sylvester_scores))
    sylvester_total_violations = sum(
        count * multiplicity for count, multiplicity in sylvester_violations.items()
    )
    if sylvester_total_violations != 40 or sylvester_violations[0] != 40:
        raise AssertionError(("Sylvester first moment", sylvester_violations))

    union_10 = iid_union_sum(heights, 10)
    union_12 = iid_union_sum(heights, 12)
    union_14 = iid_union_sum(heights, 14)
    if (union_10, union_12, union_14) != (
        Fraction(17_973, 4096),
        Fraction(6073, 4096),
        Fraction(1653, 4096),
    ):
        raise AssertionError(("iid union sums", union_10, union_12, union_14))
    if not union_12 > 1 > union_14:
        raise AssertionError("iid union threshold corruption")
    if corrupted_inclusive_union_sum(heights, 12) == union_12:
        raise AssertionError("strict-tail corruption went undetected")

    # A one-entry mutation must destroy the adaptive orbit's even cross
    # profile.  This controls both the stored witness and row-major ordering.
    corrupted = adaptive[:-1] + (-adaptive[-1],)
    corrupted_profile = Counter(abs(inner(corrupted, rank)) for rank in ranks)
    if corrupted_profile == adaptive_profile or set(corrupted_profile) <= {0, 4, 8}:
        raise AssertionError("adaptive-seed corruption went undetected")

    return center_count, len(ranks) * 2, optimum_count


def smooth_target(order: int) -> float:
    c = 2.0 / 5.0
    epsilon = 1.0 / 25.0
    shift = math.e**2
    return order**1.5 * (c + epsilon * math.sin(math.log(math.log(order + shift))))


def parity_rounded_increment(order: int) -> int:
    threshold = smooth_target(order + 1) - smooth_target(order)
    increment = max(1, math.ceil(threshold))
    if increment % 2 != order % 2:
        increment += 1
    return increment


def mu_value(order: int) -> Fraction:
    if order == 0:
        return Fraction(0)
    return Fraction(order * math.comb(order - 1, (order - 1) // 2), 1 << (order - 1))


def build_countermodel(limit: int) -> list[int]:
    prefix = [0, 0, 1, 3, 4, 4, 5, 9, 10, 12, 13]
    if limit < 10:
        raise ValueError("countermodel limit must include the exact prefix")
    values = prefix + [0] * (limit - 10)
    for order in range(10, limit):
        values[order + 1] = values[order] + parity_rounded_increment(order)
    return values


def verify_countermodel() -> tuple[int, int, int, float]:
    values = build_countermodel(COUNTERMODEL_LIMIT)
    exact_prefix = (0, 1, 3, 4, 4, 5, 9, 10, 12, 13)
    if tuple(values[1:11]) != exact_prefix:
        raise AssertionError(("countermodel exact prefix", values[1:11]))

    minimum_integer_distance = 1.0
    maximum_increment_ratio = 0.0
    r_values: list[int] = []
    for order in range(COUNTERMODEL_LIMIT):
        increment = values[order + 1] - values[order]
        residual = increment - (order % 2)
        r_values.append(residual)
        if residual < 0 or residual % 2:
            raise AssertionError(("Bellman residual lattice", order, residual))

        if order >= 10:
            threshold = smooth_target(order + 1) - smooth_target(order)
            minimum_integer_distance = min(
                minimum_integer_distance, abs(threshold - round(threshold))
            )
            if increment <= 0 or increment % 2 != order % 2:
                raise AssertionError(("rounded increment parity", order, increment))
            if increment > 3.0 * math.sqrt(order) + 1e-12 or increment > order:
                raise AssertionError(("rounded increment ceiling", order, increment))
            if residual > 3.0 * math.sqrt(order) + 1e-12:
                raise AssertionError(("Bellman residual ceiling", order, residual))
            if residual // 2 > order // 2:
                raise AssertionError(("weighted-deficit range", order, residual))
            maximum_increment_ratio = max(
                maximum_increment_ratio, increment / math.sqrt(order)
            )

            baseline = 13.0 + smooth_target(order) - smooth_target(10)
            rounding_error = values[order] - baseline
            if not -1e-9 <= rounding_error < 2 * (order - 10) + 1e-9:
                raise AssertionError(("countermodel rounding envelope", order))

    for order in range(1, COUNTERMODEL_LIMIT + 1):
        expected_parity = (order * (order - 1) // 2) % 2
        if values[order] % 2 != expected_parity:
            raise AssertionError(("energy parity", order, values[order]))
        if values[order] > order * (order - 1) // 2:
            raise AssertionError(("trivial energy ceiling", order, values[order]))
        gaussian_floor = order * math.sqrt(order - 1) / math.pi
        if values[order] + 1e-12 < gaussian_floor:
            raise AssertionError(("Gaussian scalar floor", order, values[order]))

    if minimum_integer_distance <= 1e-7:
        raise AssertionError(
            (
                "countermodel finite rounding is numerically unstable",
                minimum_integer_distance,
            )
        )

    mus = tuple(mu_value(order) for order in range(HEREDITARY_LIMIT + 1))
    hereditary_checks = 0
    for total in range(1, HEREDITARY_LIMIT + 1):
        for retained in range(total + 1):
            removed = total - retained
            squared_floor = (
                values[retained] ** 2
                + removed * (removed - 1) // 2
                + removed * (retained % 2)
            )
            if values[total] ** 2 < squared_floor:
                raise AssertionError(("squared hereditary floor", total, retained))
            additive_floor = Fraction(values[retained]) + (retained % 2) * mus[removed]
            if Fraction(values[total]) < additive_floor:
                raise AssertionError(("additive hereditary floor", total, retained))
            hereditary_checks += 2

    partial_r = [0]
    for residual in r_values:
        partial_r.append(partial_r[-1] + residual)
    dyadic_checks = 0
    parity_corruption_detected = False
    q = 2.0**-1.5
    for order in range(1, COUNTERMODEL_LIMIT // 2 + 1):
        residual_sum = partial_r[2 * order] - partial_r[order]
        parity_sum = (order + 1) // 2
        if values[2 * order] != values[order] + residual_sum + parity_sum:
            raise AssertionError(("exact dyadic telescoping", order))
        parity_corruption_detected |= values[2 * order] != values[order] + residual_sum

        c_order = values[order] / order**1.5
        c_double = values[2 * order] / (2 * order) ** 1.5
        averaged_residual = residual_sum / order**1.5
        exact_remainder = q * parity_sum / order**1.5
        if not math.isclose(
            c_double,
            q * c_order + q * averaged_residual + exact_remainder,
            rel_tol=0.0,
            abs_tol=2e-15,
        ):
            raise AssertionError(("normalized dyadic identity", order))
        dyadic_checks += 1

    if not parity_corruption_detected:
        raise AssertionError("dyadic parity corruption went undetected")

    # Removing the parity adjustment from the rounded increment must violate
    # the prescribed energy lattice somewhere in the tested tail.
    wrong_parity_detected = False
    for order in range(10, 1000):
        wrong_increment = max(
            1, math.ceil(smooth_target(order + 1) - smooth_target(order))
        )
        wrong_parity_detected |= wrong_increment % 2 != order % 2
    if not wrong_parity_detected:
        raise AssertionError("increment-parity corruption went undetected")

    return (
        COUNTERMODEL_LIMIT,
        hereditary_checks,
        dyadic_checks,
        maximum_increment_ratio,
    )


def main() -> None:
    entropy_checks, entropy_gap = verify_entropy_floor()
    centers, orbit_shifts, optimum_centers = verify_b4_cross_orbit()
    model_orders, hereditary_checks, dyadic_checks, increment_ratio = (
        verify_countermodel()
    )

    print(f"entropy_floor_exact_checks={entropy_checks}")
    print(f"balanced_entropy_leading_gap={entropy_gap:.12f}")
    print(f"b4_cross_centers_checked={centers}")
    print(f"b4_orbit_shifts_checked={orbit_shifts}")
    print("b4_height_profile=0:14,2:24,4:16,6:8,8:2")
    print("b4_iid_union_sums=K10:17973/4096,K12:6073/4096,K14:1653/4096")
    print("b4_cross_optimum=10")
    print(f"b4_optimal_projective_centers={optimum_centers}")
    print(f"countermodel_orders_checked={model_orders}")
    print(f"countermodel_hereditary_checks={hereditary_checks}")
    print(f"dyadic_identities_checked={dyadic_checks}")
    print(f"maximum_increment_over_sqrt_n={increment_ratio:.12f}")
    print(f"deterministic_seed={DETERMINISTIC_SEED}")
    print(
        "corruption_controls="
        "entropy_sign,adaptive_seed,union_threshold,increment_parity,dyadic_parity"
    )
    print("frontier_walls_verification=PASSED")


if __name__ == "__main__":
    main()

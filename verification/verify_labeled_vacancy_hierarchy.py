#!/usr/bin/env python3
"""Exact checks for Fourier and collision vacancy criteria.

The input is the labeled balanced 7+7 shell occupancy for the standard
Paley conference matrix of order 14.  This script verifies three facts.

* The negative occupancy moment at inverse temperature ``log(K)`` detects
  the unique vacant gauge exactly.
* Odd factorial-moment truncations are valid lower bounds for that moment;
  in this example degree 87 is the first truncation that certifies vacancy.
* The Fourier--Bochner principal-minor hierarchy cannot give a small witness:
  the observed one-zero, positive-gap-six spectrum rules out every witness
  on at most 6826 characters.

All decisions use integer or Fraction arithmetic.  The shell construction
and Walsh transform are imported from the independent Parseval verifier.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import comb, prod

if not __debug__:
    raise RuntimeError("verification requires Python assertions")

import verify_labeled_shell_parseval as parseval


EXPECTED_OCCUPANCY = (8192, 304_908, 1, 6, 87)
EXPECTED_FIRST_COLLISION_DEGREE = 87
EXPECTED_PSD_SAFE_DIMENSION = 6826
LOCALIZING_ROOTS = (9, 17, 26, 36, 46, 56, 67, 76, 86)
EXPECTED_NORMALIZED_MARGIN = Fraction(
    -36_795_384_082_687_448_302_747_606_657,
    4_816_759_830_492_505_652_837_357_886_720_000,
)


def c14_occupancy() -> tuple[list[int], list[int]]:
    left, right, cross, maxima = parseval.shell_arrays_c14()
    if maxima != (11, 11, 21):
        raise AssertionError(("local maxima", maxima))
    transformed = parseval.labeled_target_transform(left, right, cross, 16)
    group_size = len(transformed)
    inverse = parseval.fwht(transformed)
    if any(value % group_size for value in inverse):
        raise AssertionError("nonintegral Walsh inverse")
    occupancy = [value // group_size for value in inverse]
    if parseval.fwht(occupancy) != transformed:
        raise AssertionError("Walsh round trip")
    return occupancy, transformed


def verify_soft_vacancy(occupancy: list[int]) -> Fraction:
    group_size = len(occupancy)
    law = Counter(occupancy)
    soft_moment = sum(
        count * Fraction(1, group_size) ** value
        for value, count in law.items()
    ) / group_size
    threshold = Fraction(1, group_size)
    if not soft_moment > threshold:
        raise AssertionError("log(K) negative moment missed vacancy")

    # If all occupancies were at least one, every summand K^{-b(g)} would
    # be at most K^{-1}.  Replacing the single zero by one confirms the
    # converse side of the exact threshold on this same finite state.
    filled = occupancy.copy()
    filled[filled.index(0)] = 1
    filled_moment = sum(Fraction(1, group_size) ** value for value in filled)
    filled_moment /= group_size
    if filled_moment > threshold:
        raise AssertionError("filled occupancy passed vacancy threshold")
    return soft_moment


def verify_collision_truncations(
    occupancy: list[int], soft_moment: Fraction
) -> int:
    group_size = len(occupancy)
    law = Counter(occupancy)
    q = Fraction(group_size - 1, group_size)
    partial = Fraction(1)
    first_certifying_degree: int | None = None
    odd_lower_bounds = 0
    even_upper_bounds = 0
    maximum = max(occupancy)
    for degree in range(1, maximum + 1):
        factorial_moment = Fraction(
            sum(
                count * comb(value, degree)
                for value, count in law.items()
                if value >= degree
            ),
            group_size,
        )
        partial += (-q) ** degree * factorial_moment
        if degree % 2:
            if partial > soft_moment:
                raise AssertionError(("odd truncation is not a lower bound", degree))
            odd_lower_bounds += 1
            if partial > Fraction(1, group_size) and first_certifying_degree is None:
                first_certifying_degree = degree
        else:
            if partial < soft_moment:
                raise AssertionError(("even truncation is not an upper bound", degree))
            even_upper_bounds += 1

    if partial != soft_moment:
        raise AssertionError("complete factorial expansion mismatch")
    if first_certifying_degree != EXPECTED_FIRST_COLLISION_DEGREE:
        raise AssertionError(("first collision degree", first_certifying_degree))
    if (odd_lower_bounds, even_upper_bounds) != (44, 43):
        raise AssertionError(("truncation counts", odd_lower_bounds, even_upper_bounds))
    return first_certifying_degree


def verify_fourier_bochner_floor(
    occupancy: list[int], transformed: list[int]
) -> tuple[int, Fraction]:
    group_size = len(occupancy)
    law = Counter(occupancy)
    zero_count = law[0]
    minimum_positive = min(value for value in occupancy if value > 0)
    safe_dimension = (
        (minimum_positive - 1) * group_size
        // (minimum_positive * zero_count)
    )
    if safe_dimension != EXPECTED_PSD_SAFE_DIMENSION:
        raise AssertionError(("PSD-safe dimension", safe_dimension))

    # For c=b-1, the full normalized Fourier kernel has eigenvalues c(g).
    # Its trace coefficient and total Fourier sum are checked without a
    # floating eigensolver.  The unique empty gauge supplies eigenvalue -1.
    constant_coefficient = Fraction(transformed[0], group_size) - 1
    nontrivial_maximum = max(
        Fraction(abs(value), group_size) for value in transformed[1:]
    )
    if nontrivial_maximum > constant_coefficient:
        raise AssertionError("one-character condition unexpectedly certifies")
    kernel_transform = transformed.copy()
    kernel_transform[0] -= group_size
    if sum(kernel_transform) != -group_size:
        raise AssertionError("empty-gauge Fourier eigenvalue corruption")

    # Corruption control: omitting the subtraction of one changes the
    # negative eigenvalue from -1 to zero and cannot certify vacancy.
    if sum(transformed) != 0:
        raise AssertionError("unshifted-kernel corruption not detected")
    return safe_dimension, nontrivial_maximum


def verify_localizing_margin(occupancy: list[int]) -> Fraction:
    def witness(value: int) -> int:
        return prod(2 * value - (2 * root + 1) for root in LOCALIZING_ROOTS)

    group_size = len(occupancy)
    value_at_zero = witness(0)
    signed_sum = sum((value - 1) * witness(value) ** 2 for value in occupancy)
    margin = Fraction(signed_sum, group_size * value_at_zero**2)
    if margin != EXPECTED_NORMALIZED_MARGIN:
        raise AssertionError(("normalized localizing margin", margin))
    if not Fraction(-1, group_size) <= margin < 0:
        raise AssertionError("unique-hole margin bound failed")

    # Filling the unique zero with one changes every normalized localizing
    # form by exactly q(0)^2/K.  For the normalized witness this is 1/K.
    filled = occupancy.copy()
    filled[filled.index(0)] = 1
    filled_sum = sum((value - 1) * witness(value) ** 2 for value in filled)
    filled_margin = Fraction(filled_sum, group_size * value_at_zero**2)
    if filled_margin - margin != Fraction(1, group_size):
        raise AssertionError("fill-one rank-one perturbation failed")
    if filled_margin < 0:
        raise AssertionError("filled occupancy retained a negative witness")
    return margin


def main() -> None:
    occupancy, transformed = c14_occupancy()
    law = Counter(occupancy)
    observed = (
        len(occupancy),
        sum(occupancy),
        law[0],
        min(value for value in occupancy if value > 0),
        max(occupancy),
    )
    if observed != EXPECTED_OCCUPANCY:
        raise AssertionError(("C14 occupancy", observed))
    soft_moment = verify_soft_vacancy(occupancy)
    first_degree = verify_collision_truncations(occupancy, soft_moment)
    safe_dimension, maximum_coefficient = verify_fourier_bochner_floor(
        occupancy, transformed
    )
    normalized_margin = verify_localizing_margin(occupancy)
    print(
        "c14_vacancy="
        f"K:{observed[0]} B:{observed[1]} zeros:{observed[2]} "
        f"positive_range:{observed[3]}..{observed[4]}"
    )
    print(f"soft_inverse_temperature=log(K) first_collision_degree={first_degree}")
    print(
        f"fourier_psd_safe_dimension={safe_dimension} "
        f"maximum_nontrivial_coefficient={maximum_coefficient}"
    )
    print(f"normalized_localizing_margin={normalized_margin}")
    print(
        "corruption_controls="
        "filled_vacancy,unshifted_fourier_kernel,fill_one_localizer"
    )
    print("arithmetic=integer,fraction")
    print("labeled_vacancy_hierarchy=PASSED")


if __name__ == "__main__":
    main()

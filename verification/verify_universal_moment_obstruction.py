#!/usr/bin/env python3
"""Exact finite checks for the universal occupancy-moment obstruction.

The script verifies the explicit six-point vacant/nonvacant pair matching
moments through degree two.  It also checks, over a deterministic finite
range, the exact derivative-vector product used by the pigeonhole proof and
the simpler sufficient upper bound.  It does not claim to enumerate the
asymptotic polynomial collision.
"""

from __future__ import annotations

from math import comb, factorial, prod


if not __debug__:
    raise RuntimeError("verification requires Python assertions")


VACANT = (0, 4, 5, 6, 6, 6)
NONVACANT = (1, 2, 6, 6, 6, 6)
MAXIMUM_K = 128
MAXIMUM_K_ORDER = 32
EXPECTED_COUNT_SUMMARY = (3631, 449, 469, 20)


def raw_moment(values: tuple[int, ...], degree: int) -> int:
    return sum(value**degree for value in values)


def falling(value: int, degree: int) -> int:
    return prod(range(value - degree + 1, value + 1)) if degree else 1


def derivative_range_maximum(order: int, degree: int) -> int:
    return factorial(degree) * comb(order + 1, degree + 1)


def derivative_vector_upper_bound(order: int, multiplicity: int) -> int:
    return prod(
        derivative_range_maximum(order, degree) + 1
        for degree in range(multiplicity)
    )


def simple_vector_upper_bound(order: int, multiplicity: int) -> int:
    exponent = multiplicity * (multiplicity + 1) // 2
    return (order + 1) ** exponent


def verify_explicit_pair() -> None:
    order = 6
    assert len(VACANT) == len(NONVACANT) == order
    assert VACANT.count(0) == 1
    assert all(0 <= value <= order for value in VACANT)
    assert all(1 <= value <= order for value in NONVACANT)

    observed_moments = []
    for degree in range(4):
        left = raw_moment(VACANT, degree)
        right = raw_moment(NONVACANT, degree)
        observed_moments.append((left, right))
    assert observed_moments[:3] == [(6, 6), (27, 27), (149, 149)]
    assert observed_moments[3] == (837, 873)

    # Before the common padding at value six, the signed support polynomial
    # is 1-x-x^2+x^4+x^5-x^6.  Its first three derivatives vanish at one,
    # while the third derivative does not.
    coefficients = (1, -1, -1, 0, 1, 1, -1)
    derivatives = tuple(
        sum(
            coefficient * falling(index, degree)
            for index, coefficient in enumerate(coefficients)
        )
        for degree in range(4)
    )
    assert derivatives[:3] == (0, 0, 0)
    assert derivatives[3] != 0

    # Corruption controls: the zero is essential, and unequal padding loses
    # both cardinality and the first-moment identity.
    zero_corruption = (1, *VACANT[1:])
    assert raw_moment(zero_corruption, 1) != raw_moment(NONVACANT, 1)
    padding_corruption = NONVACANT[:-1]
    assert len(padding_corruption) != len(VACANT)
    assert raw_moment(padding_corruption, 1) != raw_moment(VACANT, 1)


def verify_derivative_counts() -> tuple[int, int, int, int]:
    cases = 0
    simple_collision_cases = 0
    exact_collision_cases = 0
    exact_only_cases = 0
    zeroth_range_corruption_detected = False
    for order in range(1, MAXIMUM_K + 1):
        maximum_multiplicity = min(order + 1, MAXIMUM_K_ORDER)
        for multiplicity in range(1, maximum_multiplicity + 1):
            exact = derivative_vector_upper_bound(order, multiplicity)
            simple = simple_vector_upper_bound(order, multiplicity)
            polynomial_count = 1 << (order + 1)

            # Directly recompute the hockey-stick derivative maximum.
            for degree in range(multiplicity):
                direct = sum(falling(index, degree) for index in range(order + 1))
                assert direct == derivative_range_maximum(order, degree)

            if multiplicity == 1:
                assert exact == order + 2
                assert exact > simple
                assert exact < polynomial_count
                zeroth_range_corruption_detected = True
            else:
                # The apparent K+2 versus K+1 loss at derivative zero is
                # paid for by the sharper first-derivative range.
                paired = (order + 2) * (order * (order + 1) // 2 + 1)
                assert paired <= (order + 1) ** 3
                assert exact <= simple

            simple_collision = simple < polynomial_count
            exact_collision = exact < polynomial_count
            assert not simple_collision or exact_collision
            simple_collision_cases += int(simple_collision)
            exact_collision_cases += int(exact_collision)
            exact_only_cases += int(exact_collision and not simple_collision)
            cases += 1

    assert zeroth_range_corruption_detected
    assert exact_only_cases > 0
    return cases, simple_collision_cases, exact_collision_cases, exact_only_cases


def main() -> None:
    verify_explicit_pair()
    cases, simple_cases, exact_cases, exact_only_cases = verify_derivative_counts()
    assert (
        cases,
        simple_cases,
        exact_cases,
        exact_only_cases,
    ) == EXPECTED_COUNT_SUMMARY
    print("explicit_pair=K:6 matched_moments:0..2 cubic:837/873")
    print(
        f"derivative_count_cases={cases} "
        f"simple_collision_cases={simple_cases} "
        f"exact_collision_cases={exact_cases} "
        f"exact_only_cases={exact_only_cases}"
    )
    print("corruption_controls=zero_removal,unequal_padding,zeroth_range")
    print("arithmetic=integer")
    print("universal_moment_obstruction=PASSED")


if __name__ == "__main__":
    main()

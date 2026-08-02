#!/usr/bin/env python3
"""Deterministic checks for the Paley interval-leakage construction.

The asymptotic theorem is analytic.  This script independently checks its
finite arithmetic inputs, the interval energy formula, representative exact
witnesses, and the Fourier-half estimate in one small case.
"""

from __future__ import annotations

import cmath
import math


CASES = (
    # p, progression level, least nonresidue, exact full witness energy
    (73, 3, 5, 265),
    (241, 5, 7, 1705),
    (2521, 7, 11, 60625),
    (9241, 11, 13, 430729),
)


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    residue = pow(value, (prime - 1) // 2, prime)
    assert residue in (1, prime - 1)
    return 1 if residue == 1 else -1


def primes_through(limit: int) -> list[int]:
    return [value for value in range(3, limit + 1, 2) if is_prime(value)]


def progression_modulus(level: int) -> int:
    result = 8
    for prime in primes_through(level):
        result *= prime
    return result


def least_nonresidue(prime: int) -> int:
    for value in range(1, prime):
        if legendre(value, prime) == -1:
            return value
    raise AssertionError("a prime field must have a quadratic nonresidue")


def interval_signs(prime: int) -> list[int]:
    half = (prime - 1) // 2
    return [-1 if value < half else 1 for value in range(prime)]


def interval_core_energy_fast(prime: int) -> int:
    """Return Q_K(f) using the exact interval autocorrelation formula."""
    half = (prime - 1) // 2
    return 4 * sum(
        (half - difference) * legendre(difference, prime)
        for difference in range(1, half)
    )


def core_energy_direct(signs: list[int], prime: int) -> int:
    result = 0
    for left in range(prime):
        for right in range(left + 1, prime):
            result += (
                legendre(left - right, prime) * signs[left] * signs[right]
            )
    return result


def fourier_half_energies(signs: list[int], prime: int) -> tuple[float, float]:
    energies = [0.0, 0.0]
    scale = math.sqrt(prime)
    for frequency in range(1, prime):
        transform = sum(
            signs[value]
            * cmath.exp(-2j * math.pi * frequency * value / prime)
            for value in range(prime)
        ) / scale
        index = 0 if legendre(frequency, prime) == 1 else 1
        energies[index] += abs(transform) ** 2
    return energies[0], energies[1]


def verify_cases() -> None:
    for prime, level, expected_nonresidue, expected_energy in CASES:
        assert is_prime(prime)
        assert prime % progression_modulus(level) == 1
        assert prime % 4 == 1
        assert all(legendre(value, prime) == 1 for value in range(1, level + 1))

        nonresidue = least_nonresidue(prime)
        assert nonresidue == expected_nonresidue
        signs = interval_signs(prime)
        assert sum(signs) == 1

        core_energy = interval_core_energy_fast(prime)
        full_energy = abs(core_energy) + 1
        assert full_energy == expected_energy

        leakage_bound = 2.0 * prime / (nonresidue - 1)
        ratio_lower_bound = (
            1.0
            - 4.0 * prime / ((prime + 1) * (nonresidue - 1))
            - (math.sqrt(prime) - 1.0) ** 2 / (prime * (prime + 1))
        )
        witness_ratio = 2.0 * full_energy / ((prime + 1) * math.sqrt(prime))
        assert witness_ratio + 1e-12 >= ratio_lower_bound

        print(
            f"p={prime} level={level} least_nonresidue={nonresidue} "
            f"S=1 Q_core={core_energy} Q_full={full_energy} "
            f"leakage_bound={leakage_bound:.12f} "
            f"witness_ratio={witness_ratio:.12f}"
        )


def verify_small_fourier_case() -> None:
    prime = 73
    signs = interval_signs(prime)
    fast_energy = interval_core_energy_fast(prime)
    assert core_energy_direct(signs, prime) == fast_energy

    residue_energy, nonresidue_energy = fourier_half_energies(signs, prime)
    parseval_target = prime - sum(signs) ** 2 / prime
    assert abs(residue_energy + nonresidue_energy - parseval_target) < 1e-9
    spectral_difference = math.sqrt(prime) * (
        residue_energy - nonresidue_energy
    )
    assert abs(spectral_difference - 2 * fast_energy) < 1e-8
    bound = 2.0 * prime / (least_nonresidue(prime) - 1)
    assert min(residue_energy, nonresidue_energy) <= bound + 1e-10
    print(
        "fourier_case=p73 "
        f"E_residue={residue_energy:.12f} "
        f"E_nonresidue={nonresidue_energy:.12f} parseval=PASS"
    )


def verify_corruption_controls() -> None:
    assert not is_prime(72)
    signs = interval_signs(73)
    signs[0] *= -1
    assert sum(signs) != 1
    assert core_energy_direct(signs, 73) != interval_core_energy_fast(73)
    print("corruption_controls=composite_modulus,interval_endpoint")


def main() -> None:
    verify_cases()
    verify_small_fourier_case()
    verify_corruption_controls()
    print("paley_least_nonresidue_verification=PASSED")


if __name__ == "__main__":
    main()

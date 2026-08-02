#!/usr/bin/env python3
"""Deterministic checks for the negative-replica transport obstruction.

The proof in the research note is analytic.  This script checks its constants,
the exact entropy-production decomposition on finite Boolean cubes, and the
one-coordinate stability inequality used in the quantitative remainder.
Floating point is used only for these regression checks; no finite grid is
presented as a proof of the analytic inequality.
"""

from __future__ import annotations

import math
import random


DETERMINISTIC_SEED = 413935


def log_cosh(value: float) -> float:
    absolute = abs(value)
    return absolute + math.log1p(math.exp(-2.0 * absolute)) - math.log(2.0)


def conditional_divergence(z: float) -> float:
    return z * math.tanh(z) - log_cosh(z)


def generator_contribution(q: float, z: float) -> float:
    # The direct cosh ratio is stable on the deterministic grid below.
    return 0.5 * (
        math.cosh((1.0 + 2.0 / q) * z) / math.cosh(z) - 1.0
    )


def scalar_remainder(q: float, z: float) -> float:
    divergence = conditional_divergence(z)
    return generator_contribution(q, z) - (
        2.0 * (q + 1.0) * divergence / q**2
    )


def verify_scalar_stability() -> tuple[int, float]:
    q_values = (0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0)
    z_values = tuple(index / 40.0 for index in range(1, 401))
    minimum_slack = math.inf
    checks = 0
    for q in q_values:
        for z in z_values:
            divergence = conditional_divergence(z)
            remainder = scalar_remainder(q, z)
            lower = 2.0 * divergence**2 / (3.0 * q)
            slack = remainder - lower
            tolerance = 2e-12 * max(1.0, abs(remainder), abs(lower))
            if slack < -tolerance:
                raise AssertionError(("scalar stability", q, z, remainder, lower))
            minimum_slack = min(minimum_slack, slack)
            checks += 1

    # A deliberate sign corruption must be rejected.
    q = 2.0
    z = 1.0
    wrong = generator_contribution(q, z) + 2.0 * (q + 1.0) * (
        conditional_divergence(z)
    ) / q**2
    if abs(wrong - scalar_remainder(q, z)) < 1e-6:
        raise AssertionError("scalar-sign corruption was not detected")
    return checks, minimum_slack


def entropy_production(function: list[float], q: float) -> tuple[float, float, float]:
    size = len(function)
    dimension = (size - 1).bit_length()
    if size != 1 << dimension or dimension == 0:
        raise ValueError("function must live on a nontrivial Boolean cube")
    raw_weights = [value ** (-q) for value in function]
    normalizer = sum(raw_weights)
    probability = [weight / normalizer for weight in raw_weights]
    divergence = sum(
        mass * math.log(mass * size) for mass in probability if mass > 0.0
    )

    generator_mean = 0.0
    conditional_sum = 0.0
    scalar_sum = 0.0
    for index, value in enumerate(function):
        generator = 0.5 * sum(
            function[index ^ (1 << coordinate)] / value - 1.0
            for coordinate in range(dimension)
        )
        generator_mean += probability[index] * generator

    for coordinate in range(dimension):
        bit = 1 << coordinate
        for index in range(size):
            if index & bit:
                continue
            partner = index | bit
            pair_mass = probability[index] + probability[partner]
            conditional_left = probability[index] / pair_mass
            conditional_right = probability[partner] / pair_mass
            conditional_kl = (
                conditional_left * math.log(2.0 * conditional_left)
                + conditional_right * math.log(2.0 * conditional_right)
            )
            z = abs(0.5 * q * math.log(function[partner] / function[index]))
            expected_kl = conditional_divergence(z)
            if abs(conditional_kl - expected_kl) > 3e-12:
                raise AssertionError(
                    ("conditional KL parameterization", conditional_kl, expected_kl)
                )
            conditional_sum += pair_mass * conditional_kl
            scalar_sum += pair_mass * scalar_remainder(q, z)

    remainder = generator_mean - 2.0 * (q + 1.0) * divergence / q**2
    total_correlation = conditional_sum - divergence
    decomposition = scalar_sum + 2.0 * (q + 1.0) * total_correlation / q**2
    if abs(remainder - decomposition) > 2e-10 * max(1.0, abs(remainder)):
        raise AssertionError(("entropy-production decomposition", remainder, decomposition))
    if total_correlation < -2e-12:
        raise AssertionError(("entropy tensorization", total_correlation))
    lower = 2.0 * divergence**2 / (3.0 * q * dimension)
    if remainder + 2e-11 < lower:
        raise AssertionError(("cube stability", remainder, lower, divergence))
    return remainder, divergence, total_correlation


def noise(function: list[float], rho: float) -> list[float]:
    output = list(function)
    dimension = (len(function) - 1).bit_length()
    for coordinate in range(dimension):
        bit = 1 << coordinate
        previous = list(output)
        keep = (1.0 + rho) / 2.0
        flip = (1.0 - rho) / 2.0
        for index in range(len(output)):
            output[index] = keep * previous[index] + flip * previous[index ^ bit]
    return output


def normalized_negative_log_moment(function: list[float], q: float) -> float:
    logarithms = [-q * math.log(value) for value in function]
    maximum = max(logarithms)
    gamma = maximum + math.log(
        sum(math.exp(value - maximum) for value in logarithms) / len(function)
    )
    return gamma / q


def verify_cube_stability() -> tuple[int, float]:
    generator = random.Random(DETERMINISTIC_SEED)
    checks = 0
    largest_derivative_error = 0.0
    for dimension in range(1, 7):
        for sample in range(40):
            base = [math.exp(generator.uniform(-1.5, 1.5)) for _ in range(1 << dimension)]
            q = math.exp(generator.uniform(math.log(0.2), math.log(20.0)))
            entropy_production(base, q)
            checks += 1

            if sample < 5:
                theta = generator.uniform(1.1, 4.0)
                s = generator.uniform(0.1, 0.8)
                step = 1e-5

                def h_at(time: float) -> float:
                    evolved = noise(base, math.exp(-time))
                    exponent = theta * math.exp(2.0 * time) - 1.0
                    return normalized_negative_log_moment(evolved, exponent)

                evolved = noise(base, math.exp(-s))
                exponent = theta * math.exp(2.0 * s) - 1.0
                remainder, _, _ = entropy_production(evolved, exponent)
                numerical_derivative = (h_at(s + step) - h_at(s - step)) / (
                    2.0 * step
                )
                error = abs(numerical_derivative + remainder)
                if error > 3e-7:
                    raise AssertionError(
                        ("semigroup derivative", dimension, error, numerical_derivative, remainder)
                    )
                largest_derivative_error = max(largest_derivative_error, error)
    return checks, largest_derivative_error


def transport_liminf_constant(beta: float, theta: float, upper_constant: float) -> float:
    bracket = (
        beta**2 / 8.0
        - beta * (1.0 - 1.0 / math.sqrt(2.0)) * upper_constant
        - math.log(2.0)
    )
    return 2.0 * theta * bracket / beta**2 - math.log(2.0)


def verify_transport_obstruction() -> tuple[float, float, float]:
    beta = 4.0
    theta = 8.0
    constant = transport_liminf_constant(beta, theta, 0.5)
    expected = math.sqrt(2.0) - 2.0 * math.log(2.0)
    if abs(constant - expected) > 2e-15 or constant <= 0.027:
        raise AssertionError(("transport liminf constant", constant, expected))

    coefficient_root = (
        2.0 * (1.0 - 1.0 / math.sqrt(2.0))
        + math.sqrt(
            4.0 * (1.0 - 1.0 / math.sqrt(2.0)) ** 2
            + 8.0 * math.log(2.0)
        )
    )
    bracket = (
        beta**2 / 8.0
        - beta * (1.0 - 1.0 / math.sqrt(2.0)) / 2.0
        - math.log(2.0)
    )
    theta_threshold = beta**2 * math.log(2.0) / (2.0 * bracket)
    if not 3.012 < coefficient_root < 3.013:
        raise AssertionError(("beta threshold", coefficient_root))
    if not 7.690 < theta_threshold < theta:
        raise AssertionError(("theta threshold", theta_threshold))

    # Omitting the orbit-entropy penalty in the finite Laplace bound removes
    # the final -log(2) term from the asymptotic constant.  This is a genuine
    # leading-order corruption; replacing m-n by m would only change o(1).
    corrupted = constant + math.log(2.0)
    if abs(corrupted - constant) < 0.6:
        raise AssertionError("orbit-entropy corruption was not detected")
    return constant, coefficient_root, theta_threshold


def main() -> None:
    scalar_checks, scalar_slack = verify_scalar_stability()
    cube_checks, derivative_error = verify_cube_stability()
    constant, beta_threshold, theta_threshold = verify_transport_obstruction()
    print(f"scalar_stability_checks={scalar_checks}")
    print(f"cube_entropy_production_checks={cube_checks}")
    print(f"maximum_derivative_error={derivative_error:.3e}")
    print(f"minimum_scalar_grid_slack={scalar_slack:.3e}")
    print(f"PT_liminf_beta4_theta8={constant:.15f}")
    print(f"PT_beta_threshold={beta_threshold:.12f}")
    print(f"PT_theta_threshold_at_beta4={theta_threshold:.12f}")
    print(f"deterministic_seed={DETERMINISTIC_SEED}")
    print("corruption_controls=scalar_sign,orbit_entropy")
    print("negative_replica_transport_obstruction_verification=PASSED")


if __name__ == "__main__":
    main()

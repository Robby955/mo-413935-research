#!/usr/bin/env python3
"""Exact verification of the labeled-shell Parseval composition bound.

The checker has two independent finite targets.

* At split 2+4, two existing alignment examples have scalar order statistic
  Lambda=0.  Their labeled zero-deficit occupancy has nonzero Fourier mass,
  which certifies gain 2.  Direct gauge enumeration gives true gains 4 and 2.
* For one fixed 7+7 split of the Paley conference matrix C14, an exact
  shellwise Walsh convolution reconstructs all 8192 fiber occupancies at the
  near-subadditive target.  It verifies Parseval exactly and demonstrates
  that the generic L2 sufficient bound does not detect the known empty fiber.

All theorem decisions use Python integers or Fraction.  Floating point is
used only for concise descriptive output.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter, defaultdict
from fractions import Fraction

if not __debug__:
    raise RuntimeError("verification requires Python assertions")


Matrix = tuple[tuple[int, ...], ...]
Spin = tuple[int, ...]

LEFT_SIGNS = (1,)
RIGHT_SIGNS = (-1, -1, -1, -1, -1, 1)
CROSS_HIGH = (
    (-1, -1, -1, -1),
    (-1, -1, 1, 1),
)
CROSS_LOW = (
    (-1, -1, -1, -1),
    (-1, 1, -1, 1),
)


def projective_spins(order: int) -> tuple[Spin, ...]:
    return tuple(
        (1,) + tail for tail in itertools.product((-1, 1), repeat=order - 1)
    )


def quadratic(signs: tuple[int, ...], spin: Spin) -> int:
    return sum(
        sign * spin[row] * spin[column]
        for sign, (row, column) in zip(
            signs, itertools.combinations(range(len(spin)), 2), strict=True
        )
    )


def bilinear(cross: Matrix, left: Spin, right: Spin) -> int:
    return sum(
        cross[row][column] * left[row] * right[column]
        for row in range(len(left))
        for column in range(len(right))
    )


def multiply_spins(left: Spin, right: Spin) -> Spin:
    return tuple(a * b for a, b in zip(left, right, strict=True))


def spin_mask(spin: Spin) -> int:
    """Encode the projective coordinates after the fixed first coordinate."""

    return sum(1 << index for index, value in enumerate(spin[1:]) if value < 0)


def fwht(values: list[int]) -> list[int]:
    """Unnormalized Walsh-Hadamard transform."""

    output = values.copy()
    width = 1
    while width < len(output):
        for start in range(0, len(output), 2 * width):
            for offset in range(width):
                low = start + offset
                high = low + width
                a = output[low]
                b = output[high]
                output[low] = a + b
                output[high] = a - b
        width *= 2
    return output


def graph_matrix(signs: tuple[int, ...], order: int) -> Matrix:
    matrix = [[0] * order for _ in range(order)]
    for sign, (row, column) in zip(
        signs, itertools.combinations(range(order), 2), strict=True
    ):
        matrix[row][column] = matrix[column][row] = sign
    return tuple(tuple(row) for row in matrix)


def matrix_energy(matrix: Matrix, spin: Spin) -> int:
    return sum(
        matrix[row][column] * spin[row] * spin[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def switch_signs(
    signs: tuple[int, ...], diagonal: Spin
) -> tuple[int, ...]:
    return tuple(
        sign * diagonal[row] * diagonal[column]
        for sign, (row, column) in zip(
            signs,
            itertools.combinations(range(len(diagonal)), 2),
            strict=True,
        )
    )


def full_2_plus_4_signs(
    left_signs: tuple[int, ...],
    right_signs: tuple[int, ...],
    cross: Matrix,
    tau: int,
) -> tuple[int, ...]:
    edge_order = tuple(itertools.combinations(range(6), 2))
    left_edges = {
        edge: sign
        for edge, sign in zip(
            itertools.combinations(range(2), 2), left_signs, strict=True
        )
    }
    right_edges = {
        (row + 2, column + 2): tau * sign
        for (row, column), sign in zip(
            itertools.combinations(range(4), 2), right_signs, strict=True
        )
    }
    return tuple(
        left_edges[(row, column)]
        if column < 2
        else right_edges[(row, column)]
        if row >= 2
        else cross[row][column - 2]
        for row, column in edge_order
    )


def small_fiber_occupancy(cross: Matrix) -> tuple[list[int], list[int], int]:
    left_spins = projective_spins(2)
    right_spins = projective_spins(4)
    left_maximum = max(abs(quadratic(LEFT_SIGNS, spin)) for spin in left_spins)
    right_maximum = max(
        abs(quadratic(RIGHT_SIGNS, spin)) for spin in right_spins
    )
    cross_maximum = max(
        abs(bilinear(cross, left, right))
        for left in left_spins
        for right in right_spins
    )

    deficits: list[int] = []
    occupancy = [0] * 32
    for sigma in (-1, 1):
        for z in left_spins:
            left_deficit = left_maximum - sigma * quadratic(LEFT_SIGNS, z)
            for eta in (-1, 1):
                for w in right_spins:
                    right_deficit = right_maximum - eta * quadratic(
                        RIGHT_SIGNS, w
                    )
                    for x in left_spins:
                        for y in right_spins:
                            cross_deficit = cross_maximum - abs(
                                bilinear(cross, x, y)
                            )
                            total_deficit = (
                                left_deficit + right_deficit + cross_deficit
                            )
                            deficits.append(total_deficit)
                            if total_deficit != 0:
                                continue
                            alpha = spin_mask(multiply_spins(z, x))
                            beta = spin_mask(multiply_spins(w, y))
                            tau = int(sigma * eta < 0)
                            group_index = alpha | (beta << 1) | (tau << 4)
                            occupancy[group_index] += 1

    return occupancy, sorted(deficits), left_maximum + right_maximum + cross_maximum


def small_true_gain(cross: Matrix, independent_ceiling: int) -> int:
    left_spins = projective_spins(2)
    right_spins = projective_spins(4)
    full_spins = projective_spins(6)
    best = math.inf
    for alpha in left_spins:
        switched_left = switch_signs(LEFT_SIGNS, alpha)
        for beta in right_spins:
            switched_right = switch_signs(RIGHT_SIGNS, beta)
            for tau in (-1, 1):
                signing = full_2_plus_4_signs(
                    switched_left, switched_right, cross, tau
                )
                maximum = max(
                    abs(quadratic(signing, spin)) for spin in full_spins
                )
                best = min(best, maximum)
    return independent_ceiling - int(best)


def verify_small_example(name: str, cross: Matrix) -> None:
    occupancy, deficits, independent_ceiling = small_fiber_occupancy(cross)
    group_size = 32
    if len(deficits) != group_size * group_size:
        raise AssertionError((name, "product size", len(deficits)))
    if any(deficit % 2 for deficit in deficits):
        raise AssertionError((name, "deficit parity"))
    if deficits[group_size - 1] != 0:
        raise AssertionError((name, "scalar Lambda", deficits[group_size - 1]))
    expected_occupancy = Counter({0: 8, 1: 16, 2: 8})
    if Counter(occupancy) != expected_occupancy:
        raise AssertionError((name, "occupancy law", Counter(occupancy)))

    transformed = fwht(occupancy)
    if transformed[0] != group_size:
        raise AssertionError((name, "Fourier mean", transformed[0]))
    expected_coefficients = (
        {2: 16, 12: 16}
        if name == "high"
        else {21: -16, 27: 16}
    )
    observed_coefficients = {
        index: value
        for index, value in enumerate(transformed)
        if index and value
    }
    if observed_coefficients != expected_coefficients:
        raise AssertionError(
            (name, "nontrivial Fourier coefficients", observed_coefficients)
        )

    fourier_variance = Fraction(
        sum(value * value for value in transformed[1:]), group_size * group_size
    )
    direct_variance = sum(
        (Fraction(value) - 1) ** 2 for value in occupancy
    ) / group_size
    if (
        fourier_variance != Fraction(1, 2)
        or direct_variance != Fraction(1, 2)
    ):
        raise AssertionError(
            (name, "Parseval variance", fourier_variance, direct_variance)
        )
    l2_certifies = Fraction(0) < fourier_variance / (group_size - 1)
    if not l2_certifies:
        raise AssertionError((name, "labeled bonus not detected"))

    # Here B_s=K, so mu=1.  Any positive variance makes
    # mu-sqrt(V/(K-1))<1, forcing an empty integer-valued fiber.
    certified_gain = 2 if l2_certifies else 0
    if certified_gain != 2:
        raise AssertionError((name, "certified gain", certified_gain))
    if min(occupancy) != 0:
        raise AssertionError((name, "Parseval prediction has no empty fiber"))
    true_gain = small_true_gain(cross, independent_ceiling)
    expected_true_gain = 4 if name == "high" else 2
    if true_gain != expected_true_gain:
        raise AssertionError((name, "true gain", true_gain))

    # Corruption control: dropping the relative orientation bit merges
    # distinct fibers and destroys the certified occupancy law.
    corrupted = [0] * group_size
    for index, value in enumerate(occupancy):
        corrupted[index & ~(1 << 4)] += value
    if Counter(corrupted) == expected_occupancy:
        raise AssertionError((name, "relative-orientation corruption"))

    coefficient_text = ",".join(
        f"{index}:{value}" for index, value in expected_coefficients.items()
    )
    print(
        f"small_{name}=lambda:0 occupancy:0x8,1x16,2x8 "
        f"fourier:{coefficient_text} V:1/2 "
        f"certified_gain:{certified_gain} true_gain:{true_gain}"
    )


def legendre_symbol(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def paley_conference(prime: int) -> Matrix:
    if prime % 4 != 1:
        raise ValueError("prime must be 1 modulo 4")
    order = prime + 1
    matrix = [[0] * order for _ in range(order)]
    for column in range(1, order):
        matrix[0][column] = matrix[column][0] = 1
    for row in range(1, order):
        for column in range(row + 1, order):
            sign = legendre_symbol((row - 1) - (column - 1), prime)
            matrix[row][column] = matrix[column][row] = sign
    return tuple(tuple(row) for row in matrix)


def transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(row) for row in zip(*matrix, strict=True))


def multiply(left: Matrix, right: Matrix) -> Matrix:
    columns = transpose(right)
    return tuple(
        tuple(
            sum(a * b for a, b in zip(row, column, strict=True))
            for column in columns
        )
        for row in left
    )


def assert_conference(matrix: Matrix) -> None:
    product = multiply(matrix, matrix)
    order = len(matrix)
    expected = tuple(
        tuple(order - 1 if row == column else 0 for column in range(order))
        for row in range(order)
    )
    if product != expected:
        raise AssertionError("C14 conference identity")


def principal(matrix: Matrix, vertices: tuple[int, ...]) -> Matrix:
    return tuple(
        tuple(matrix[row][column] for column in vertices) for row in vertices
    )


def rectangular_energy(cross: Matrix, left: Spin, right: Spin) -> int:
    return sum(
        cross[row][column] * left[row] * right[column]
        for row in range(len(left))
        for column in range(len(right))
    )


def shell_arrays_c14() -> tuple[
    dict[int, list[int]],
    dict[int, list[int]],
    dict[int, list[int]],
    tuple[int, int, int],
]:
    conference = paley_conference(13)
    assert_conference(conference)
    left_vertices = tuple(range(7))
    right_vertices = tuple(range(7, 14))
    left = principal(conference, left_vertices)
    right = principal(conference, right_vertices)
    cross = tuple(
        tuple(conference[row][column] for column in right_vertices)
        for row in left_vertices
    )
    spins = projective_spins(7)
    left_maximum = max(abs(matrix_energy(left, spin)) for spin in spins)
    right_maximum = max(abs(matrix_energy(right, spin)) for spin in spins)
    cross_maximum = max(
        abs(rectangular_energy(cross, x, y)) for x in spins for y in spins
    )
    maxima = (left_maximum, right_maximum, cross_maximum)
    if maxima != (11, 11, 21):
        raise AssertionError(("C14 local maxima", maxima))

    group_size = 1 << 13
    left_shells: defaultdict[int, list[int]] = defaultdict(
        lambda: [0] * group_size
    )
    right_shells: defaultdict[int, list[int]] = defaultdict(
        lambda: [0] * group_size
    )
    cross_shells: defaultdict[int, list[int]] = defaultdict(
        lambda: [0] * group_size
    )

    for spin in spins:
        value = matrix_energy(left, spin)
        for orientation in (-1, 1):
            deficit = left_maximum - orientation * value
            index = spin_mask(spin) | (int(orientation < 0) << 12)
            left_shells[deficit][index] += 1
    for spin in spins:
        value = matrix_energy(right, spin)
        for orientation in (-1, 1):
            deficit = right_maximum - orientation * value
            index = (spin_mask(spin) << 6) | (int(orientation < 0) << 12)
            right_shells[deficit][index] += 1
    for x in spins:
        for y in spins:
            deficit = cross_maximum - abs(rectangular_energy(cross, x, y))
            index = spin_mask(x) | (spin_mask(y) << 6)
            cross_shells[deficit][index] += 1

    return dict(left_shells), dict(right_shells), dict(cross_shells), maxima


def scalar_target_count(
    left_shells: dict[int, list[int]],
    right_shells: dict[int, list[int]],
    cross_shells: dict[int, list[int]],
    cutoff: int,
) -> int:
    total = 0
    for left_deficit, left_values in left_shells.items():
        left_count = sum(left_values)
        for right_deficit, right_values in right_shells.items():
            right_count = sum(right_values)
            for cross_deficit, cross_values in cross_shells.items():
                if left_deficit + right_deficit + cross_deficit <= cutoff:
                    total += left_count * right_count * sum(cross_values)
    return total


def transformed_shells(shells: dict[int, list[int]]) -> dict[int, list[int]]:
    return {deficit: fwht(values) for deficit, values in shells.items()}


def labeled_target_transform(
    left_shells: dict[int, list[int]],
    right_shells: dict[int, list[int]],
    cross_shells: dict[int, list[int]],
    cutoff: int,
) -> list[int]:
    left_fourier = transformed_shells(left_shells)
    right_fourier = transformed_shells(right_shells)
    cross_fourier = transformed_shells(cross_shells)
    group_size = len(next(iter(left_shells.values())))
    output = [0] * group_size
    for left_deficit, left_values in left_fourier.items():
        for right_deficit, right_values in right_fourier.items():
            for cross_deficit, cross_values in cross_fourier.items():
                if left_deficit + right_deficit + cross_deficit > cutoff:
                    continue
                for index in range(group_size):
                    output[index] += (
                        left_values[index]
                        * right_values[index]
                        * cross_values[index]
                    )
    return output


def verify_c14_split() -> None:
    left_shells, right_shells, cross_shells, maxima = shell_arrays_c14()
    group_size = 8192
    # The balanced child target is 18*sqrt(2).  Its square is 648, so exact
    # integer comparison places it strictly between the odd energies 25 and
    # 27.  Since L=43, energy >=27 is exactly deficit <=16.
    target_square = 648
    if not 25 * 25 < target_square < 27 * 27:
        raise AssertionError("C14 target lattice interval")
    if math.comb(14, 2) % 2 != 1 or sum(maxima) != 43:
        raise AssertionError("C14 target parity")
    cutoff = sum(maxima) - 27
    if cutoff != 16:
        raise AssertionError(("C14 target cutoff", cutoff))
    scalar_count = scalar_target_count(
        left_shells, right_shells, cross_shells, cutoff
    )
    if scalar_count != 304_908:
        raise AssertionError(("C14 scalar target count", scalar_count))

    transformed = labeled_target_transform(
        left_shells, right_shells, cross_shells, cutoff
    )
    if transformed[0] != scalar_count:
        raise AssertionError(("C14 Fourier mean", transformed[0], scalar_count))
    nontrivial_count = sum(value != 0 for value in transformed[1:])
    if nontrivial_count != 8159:
        raise AssertionError(("C14 nontrivial coefficient count", nontrivial_count))

    inverse = fwht(transformed)
    if any(value % group_size for value in inverse):
        raise AssertionError("C14 nonintegral inverse transform")
    occupancy = [value // group_size for value in inverse]
    if sum(occupancy) != scalar_count:
        raise AssertionError(("C14 occupancy sum", sum(occupancy)))
    if (min(occupancy), max(occupancy), occupancy.count(0)) != (0, 87, 1):
        raise AssertionError(
            (
                "C14 occupancy range/zeros",
                min(occupancy),
                max(occupancy),
                occupancy.count(0),
            )
        )

    nontrivial_square_sum = sum(value * value for value in transformed[1:])
    fourier_variance = Fraction(
        nontrivial_square_sum, group_size * group_size
    )
    mean = Fraction(scalar_count, group_size)
    direct_variance = sum((Fraction(value) - mean) ** 2 for value in occupancy)
    direct_variance /= group_size
    if fourier_variance != direct_variance:
        raise AssertionError(
            ("C14 Parseval", fourier_variance, direct_variance)
        )

    # The L2 sufficient condition is
    # mu - sqrt(V/(K-1)) < 1.  Its negation can be checked exactly by
    # squaring because mu>1 and V>=0.
    l2_certifies = (mean - 1) ** 2 < fourier_variance / (group_size - 1)
    if l2_certifies:
        raise AssertionError("C14 generic L2 bound unexpectedly certifies")
    if min(occupancy) != 0:
        raise AssertionError("C14 labeled target has no empty fiber")

    # Corruption controls for the strict target convention and Walsh inverse.
    wrong_cutoff_count = scalar_target_count(
        left_shells, right_shells, cross_shells, cutoff - 2
    )
    if wrong_cutoff_count == scalar_count:
        raise AssertionError("C14 target-cutoff corruption")
    test_vector = [index % 7 - 3 for index in range(32)]
    twice_transformed = fwht(fwht(test_vector))
    if twice_transformed != [32 * value for value in test_vector]:
        raise AssertionError("Walsh inverse corruption")

    print(
        "c14_split="
        f"maxima:{maxima[0]},{maxima[1]},{maxima[2]} "
        f"B:{scalar_count} K:{group_size} "
        f"nontrivial_fourier:{nontrivial_count} "
        f"V:{fourier_variance} occupancy:{min(occupancy)}..{max(occupancy)} "
        f"zero_fibers:{occupancy.count(0)} "
        "l2_certifies:false"
    )


def main() -> None:
    verify_small_example("high", CROSS_HIGH)
    verify_small_example("low", CROSS_LOW)
    verify_c14_split()
    print("corruption_controls=relative_orientation,target_cutoff,walsh_inverse")
    print("arithmetic=integer,fraction")
    print("labeled_shell_parseval_verification=PASSED")


if __name__ == "__main__":
    main()

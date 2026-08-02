#!/usr/bin/env python3
"""Finite exact checks for the fixed-half cut-discrepancy equivalence.

The analytic theorem is proved in the research note.  This script exhausts
all signings through order six, independently computes both minimizations,
checks the pointwise fixed-layer identity, and includes a centering corruption
control.  Every pass/fail decision uses integer arithmetic.
"""

from __future__ import annotations

import math
from itertools import product


EXPECTED_F = {2: 1, 3: 3, 4: 4, 5: 4, 6: 5}


def projective_spins(order: int) -> list[tuple[int, ...]]:
    return [(1, *tail) for tail in product((-1, 1), repeat=order - 1)]


def main() -> None:
    rows = []
    corrupted_centering_detected = False
    pointwise_checks = 0
    switching_checks = 0

    for order in range(2, 7):
        edges = [(i, j) for i in range(order) for j in range(i + 1, order)]
        edge_count = len(edges)
        fixed_negative_count = edge_count // 2
        spins = projective_spins(order)

        minimum_maximum = None
        minimum_four_discrepancy = None

        for mask in range(1 << edge_count):
            negative_count = mask.bit_count()
            energies = []
            for spin in spins:
                energy = 0
                for edge_index, (left, right) in enumerate(edges):
                    sign = -1 if (mask >> edge_index) & 1 else 1
                    energy += sign * spin[left] * spin[right]
                energies.append(energy)

            maximum = max(abs(energy) for energy in energies)
            if minimum_maximum is None or maximum < minimum_maximum:
                minimum_maximum = maximum

            # Switching by a spin changes the total signing sum to Q_A(spin).
            assert min(abs(energy) for energy in energies) ** 2 <= edge_count
            switching_checks += 1

            if negative_count != fixed_negative_count:
                continue

            four_discrepancy = 0
            corrupted_four_discrepancy = 0
            for spin in spins:
                side = [index for index, value in enumerate(spin) if value == 1]
                cut_size = len(side) * (order - len(side))
                negative_cut = 0
                for edge_index, (left, right) in enumerate(edges):
                    if spin[left] != spin[right] and (mask >> edge_index) & 1:
                        negative_cut += 1

                # Four times |e_G(S,S^c) - |S||S^c|/2|.
                four_discrepancy = max(
                    four_discrepancy,
                    2 * abs(2 * negative_cut - cut_size),
                )
                # Deliberately wrong density-one centering.
                corrupted_four_discrepancy = max(
                    corrupted_four_discrepancy,
                    4 * abs(negative_cut - cut_size),
                )

            target_total = edge_count % 2
            assert abs(four_discrepancy - maximum) <= target_total
            pointwise_checks += 1
            if abs(corrupted_four_discrepancy - maximum) > target_total:
                corrupted_centering_detected = True

            if (
                minimum_four_discrepancy is None
                or four_discrepancy < minimum_four_discrepancy
            ):
                minimum_four_discrepancy = four_discrepancy

        assert minimum_maximum == EXPECTED_F[order]
        assert minimum_four_discrepancy is not None
        assert minimum_maximum - 1 <= minimum_four_discrepancy
        excess = minimum_four_discrepancy - minimum_maximum - 2
        assert excess <= 0 or excess * excess <= edge_count
        rows.append((order, minimum_maximum, minimum_four_discrepancy))

    assert corrupted_centering_detected
    print(
        "orders="
        + ",".join(
            f"{order}:F={minimum}:4H={four_h}"
            for order, minimum, four_h in rows
        )
    )
    print(f"fixed_layer_pointwise_checks={pointwise_checks}")
    print(f"switching_mean_square_checks={switching_checks}")
    print("corruption_control=wrong_cut_centering_detected")
    print("cut_discrepancy_equivalence_verification=PASSED")


if __name__ == "__main__":
    main()

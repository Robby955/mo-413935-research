#!/usr/bin/env python3
"""Finite checks for the negative-replica composition program.

For a signing ``A`` of order ``n`` put

    f_{n,t}(A) = Z_A(t) / (2^(n+1) cosh(t)^binom(n,2)),
    G_n(q,t)   = log E_A f_{n,t}(A)^(-q).

The accompanying note proves analytically that

    G_{n+k}(q,t) >= G_n(q,t) + G_k(q,t).

This script first checks that inequality by direct enumeration in small
orders.  It then evaluates the reverse-hypercontractive transport defect

    Delta_n(theta,beta)
      = (q_{2n}/q_n) G_n(q_n,beta/sqrt(n))
        - G_n(q_{2n},beta/sqrt(2n)),

where q_r = theta/tanh(beta/sqrt(r))^2 - 1.  The second calculation sums
over every labeled switching-normalized disorder through order nine.  Nauty
``geng`` supplies one representative of each unlabeled residual graph and
``countg`` supplies its automorphism-group order, so the exact labeled
multiplicity is (n-1)!/|Aut(G)|.

All graph, energy, and multiplicity calculations are exact.  Exponentials,
logarithms, and the displayed transport defects use double precision and are
finite evidence only.  In particular, a positive finite defect is not a
disproof of asymptotic transport saturation.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import math
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass

import research_exact_small_n as small


DETERMINISTIC_SEED = 413935
PARAMETERS = ((1.0, 1.0), (2.0, 4.0))
EXPECTED_DEFECT_RATIOS = {
    (1.0, 1.0): {
        4: 0.0049185306,
        5: 0.0087837540,
        6: 0.0116569524,
        7: 0.0138659997,
        8: 0.0156035833,
        9: 0.0170007335,
    },
    (2.0, 4.0): {
        4: 0.1484980782,
        5: 0.2771315131,
        6: 0.3275853637,
        7: 0.2860893478,
        8: 0.3076235038,
        9: 0.3159700451,
    },
}


@dataclass(frozen=True)
class DisorderClass:
    graph6: bytes
    labeled_multiplicity: int
    absolute_energy_histogram: tuple[tuple[int, int], ...]


def logsumexp(values: list[float]) -> float:
    if not values:
        raise ValueError("logsumexp requires a nonempty input")
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def log_cosh(value: float) -> float:
    absolute = abs(value)
    return absolute + math.log1p(math.exp(-2.0 * absolute)) - math.log(2.0)


def root_disorders(order: int):
    """Yield every labeled switching-normalized residual adjacency."""
    residual_order = order - 1
    residual_edges = tuple(itertools.combinations(range(residual_order), 2))
    for edge_mask in range(1 << len(residual_edges)):
        adjacency = [0] * residual_order
        for edge_index, (row, column) in enumerate(residual_edges):
            if edge_mask >> edge_index & 1:
                adjacency[row] |= 1 << column
                adjacency[column] |= 1 << row
        yield tuple(adjacency)


def absolute_energy_histogram(
    adjacency: tuple[int, ...],
) -> tuple[tuple[int, int], ...]:
    histogram = Counter(
        abs(small.energy_from_cut(adjacency, subset))
        for subset in range(1 << len(adjacency))
    )
    return tuple(sorted(histogram.items()))


def log_normalized_partition_from_histogram(
    order: int, histogram: tuple[tuple[int, int], ...], temperature: float
) -> float:
    terms = [
        math.log(count) + log_cosh(temperature * energy)
        for energy, count in histogram
    ]
    return (
        logsumexp(terms)
        - (order - 1) * math.log(2.0)
        - math.comb(order, 2) * log_cosh(temperature)
    )


def direct_gamma(order: int, exponent: float, temperature: float) -> float:
    terms = []
    for adjacency in root_disorders(order):
        histogram = absolute_energy_histogram(adjacency)
        log_partition = log_normalized_partition_from_histogram(
            order, histogram, temperature
        )
        terms.append(-exponent * log_partition)
    disorder_bits = math.comb(order - 1, 2)
    return logsumexp(terms) - disorder_bits * math.log(2.0)


def verify_supermultiplicativity() -> int:
    checks = 0
    for left, right in ((1, 1), (1, 2), (2, 2), (2, 3)):
        total = left + right
        for exponent in (0.5, 2.0):
            for temperature in (0.0, 0.3, 1.1):
                parent = direct_gamma(total, exponent, temperature)
                children = direct_gamma(left, exponent, temperature) + direct_gamma(
                    right, exponent, temperature
                )
                if parent + 2e-12 < children:
                    raise AssertionError(
                        (
                            "negative-replica supermultiplicativity",
                            left,
                            right,
                            exponent,
                            temperature,
                            parent,
                            children,
                        )
                    )
                if temperature == 0.0 and (
                    abs(parent) > 2e-12 or abs(children) > 2e-12
                ):
                    raise AssertionError("t=0 normalization corruption")
                checks += 1
    return checks


def locate_countg() -> str:
    candidates = (
        shutil.which("countg"),
        "/opt/homebrew/opt/nauty/bin/countg",
        "/usr/local/bin/countg",
    )
    for candidate in candidates:
        if candidate:
            return candidate
    raise FileNotFoundError("nauty countg was not found")


def automorphism_orders(countg: str, records: tuple[bytes, ...]) -> tuple[int, ...]:
    process = subprocess.run(
        [countg, "-q", "-V", "--a", "-1"],
        input=b"\n".join(records) + b"\n",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    orders: list[int] = []
    expected_index = 1
    for line in process.stdout.decode("utf-8", errors="strict").splitlines():
        if not line.startswith("Graph "):
            continue
        label, value = line.split(":", 1)
        graph_index = int(label.split()[1])
        if graph_index != expected_index:
            raise AssertionError(("countg graph order", graph_index, expected_index))
        orders.append(int(value.strip()))
        expected_index += 1
    if len(orders) != len(records):
        raise AssertionError(("countg record count", len(orders), len(records)))
    return tuple(orders)


def class_catalogue(order: int, countg: str) -> tuple[DisorderClass, ...]:
    residual_order = order - 1
    records = tuple(small.geng_records(small.locate_geng(None), residual_order))
    expected_count = small.UNLABELED_GRAPH_COUNTS[residual_order]
    if len(records) != expected_count:
        raise AssertionError(("incomplete geng stream", order, len(records)))

    digest = hashlib.sha256()
    for record in records:
        digest.update(record + b"\n")
    expected_digest = small.EXPECTED_GENG_STREAM_SHA256[order]
    if digest.hexdigest() != expected_digest:
        raise AssertionError(("geng stream digest", order, digest.hexdigest()))

    group_orders = automorphism_orders(countg, records)
    factorial = math.factorial(residual_order)
    classes: list[DisorderClass] = []
    for record, group_order in zip(records, group_orders, strict=True):
        if factorial % group_order:
            raise AssertionError(("invalid automorphism order", order, group_order))
        adjacency = small.graph6_adjacency(record, residual_order)
        classes.append(
            DisorderClass(
                graph6=record,
                labeled_multiplicity=factorial // group_order,
                absolute_energy_histogram=absolute_energy_histogram(adjacency),
            )
        )

    multiplicity_total = sum(entry.labeled_multiplicity for entry in classes)
    expected_total = 1 << math.comb(residual_order, 2)
    if multiplicity_total != expected_total:
        raise AssertionError(
            ("labeled multiplicity total", order, multiplicity_total, expected_total)
        )
    return tuple(classes)


def class_gamma(
    order: int,
    classes: tuple[DisorderClass, ...],
    exponent: float,
    temperature: float,
    *,
    ignore_multiplicities: bool = False,
) -> float:
    terms = []
    for entry in classes:
        log_partition = log_normalized_partition_from_histogram(
            order, entry.absolute_energy_histogram, temperature
        )
        log_weight = 0.0 if ignore_multiplicities else math.log(entry.labeled_multiplicity)
        terms.append(log_weight - exponent * log_partition)
    denominator = (
        math.log(len(classes))
        if ignore_multiplicities
        else math.comb(order - 1, 2) * math.log(2.0)
    )
    return logsumexp(terms) - denominator


def transport_defect(
    order: int,
    classes: tuple[DisorderClass, ...],
    beta: float,
    theta: float,
) -> float:
    natural_temperature = beta / math.sqrt(order)
    contracted_temperature = beta / math.sqrt(2 * order)
    natural_u = math.tanh(natural_temperature)
    contracted_u = math.tanh(contracted_temperature)
    natural_q = theta / natural_u**2 - 1.0
    contracted_q = theta / contracted_u**2 - 1.0
    if natural_q <= 0.0:
        raise AssertionError(("nonpositive natural exponent", order, beta, theta))
    invariant_error = abs(
        (1.0 + natural_q) * natural_u**2
        - (1.0 + contracted_q) * contracted_u**2
    )
    if invariant_error > 2e-13:
        raise AssertionError(("invariant parameter curve", invariant_error))

    natural = class_gamma(order, classes, natural_q, natural_temperature)
    contracted = class_gamma(order, classes, contracted_q, contracted_temperature)
    defect = contracted_q / natural_q * natural - contracted
    if defect < -2e-9:
        raise AssertionError(
            ("reverse-hypercontractive defect sign", order, beta, theta, defect)
        )
    return max(0.0, defect)


def verify_transport(max_order: int) -> dict[tuple[float, float], dict[int, float]]:
    countg = locate_countg()
    results = {parameters: {} for parameters in PARAMETERS}
    unweighted_corruption_detected = False
    for order in range(4, max_order + 1):
        classes = class_catalogue(order, countg)
        for beta, theta in PARAMETERS:
            ratio = transport_defect(order, classes, beta, theta) / order**2
            expected = EXPECTED_DEFECT_RATIOS[(beta, theta)][order]
            # This is a floating-point regression check, not an exact
            # certificate.  The stored values are rounded to ten decimals.
            if abs(ratio - expected) > 2e-7:
                raise AssertionError(
                    ("transport defect regression", order, beta, theta, ratio, expected)
                )
            results[(beta, theta)][order] = ratio

        if order == min(8, max_order):
            beta, theta = PARAMETERS[0]
            temperature = beta / math.sqrt(order)
            exponent = theta / math.tanh(temperature) ** 2 - 1.0
            weighted = class_gamma(order, classes, exponent, temperature)
            unweighted = class_gamma(
                order,
                classes,
                exponent,
                temperature,
                ignore_multiplicities=True,
            )
            unweighted_corruption_detected = abs(weighted - unweighted) > 1e-6

    if not unweighted_corruption_detected:
        raise AssertionError("unlabeled-class multiplicity corruption was not detected")
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=9)
    arguments = parser.parse_args()
    if not 8 <= arguments.max_order <= 9:
        raise SystemExit("require 8 <= max-order <= 9")

    supermultiplicativity_checks = verify_supermultiplicativity()
    results = verify_transport(arguments.max_order)
    print(f"negative_replica_supermultiplicativity_checks={supermultiplicativity_checks}")
    for beta, theta in PARAMETERS:
        encoded = ",".join(
            f"{order}:{ratio:.10f}"
            for order, ratio in sorted(results[(beta, theta)].items())
        )
        print(f"transport_beta_{beta:g}_theta_{theta:g}_delta_over_n2={encoded}")
    print(f"transport_orders=4..{arguments.max_order}")
    print(f"deterministic_seed={DETERMINISTIC_SEED}")
    print("corruption_controls=unlabeled_multiplicity,t0_normalization")
    print("negative_replica_transport_verification=PASSED")


if __name__ == "__main__":
    main()

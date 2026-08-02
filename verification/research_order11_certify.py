#!/usr/bin/env python3
"""Exact nauty certificate for F(11)=17 and witness check for F(12)=18.

After switching all edges at a root positive, an order-11 signing is encoded
by a graph G on the ten remaining vertices.  If e=|E(G)| and S is the set of
negative residual spins, then

    Q(S) = 55 - 2e - 2|S|(11-|S|) + 4|delta_G(S)|.

If M<=15, the all-positive state and global-negation/complement symmetry let
us assume 20<=e<=22.  Singleton states then force maximum degree at most six.
The default nauty command therefore generates exactly the 2,153,606
unlabeled graphs in this reduced range.  Every cut energy is evaluated with
integer bit operations and no graph survives.

Passing ``--full-stream`` instead reads all 12,005,168 unlabeled graphs on ten
vertices and applies the reductions internally.  Both modes commit to their
complete newline-delimited graph6 streams with SHA-256 hashes.

Completeness trusts nauty ``geng`` to emit one representative of every
unlabeled graph.  No MILP, SAT, floating point, or timeout status is used.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import math
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass

import research_exact_small_n as small


DETERMINISTIC_SEED = 413935
RESIDUAL_ORDER = 10
EDGE_PAIRS = tuple(
    (row, column)
    for column in range(1, RESIDUAL_ORDER)
    for row in range(column)
)
EDGE_COUNT = len(EDGE_PAIRS)
FULL_RECORDS = 12_005_168
FILTERED_RECORDS = 2_153_606
FULL_SHA256 = "5650c7c979fdffd8c0f99a2f2ee8775938ec2a3dd69aa65be1207936824fc5b3"
FILTERED_SHA256 = "b62da4d7ebfaab4ccd801fd509f2fc85f6c2b815c8c1d2e969de7aa6a82c322d"
ORDER11_WITNESS = b"ICRbczQMo"
ORDER12_WITNESS = b"JWUuDOR\\K{?"


def graph6_edge_mask(record: bytes, order: int) -> int:
    record = record.strip()
    if not record or record[0] - 63 != order or order > 62:
        raise ValueError(("graph6 order", record, order))
    needed = math.comb(order, 2)
    mask = 0
    bit_index = 0
    for byte in record[1:]:
        value = byte - 63
        if not 0 <= value < 64:
            raise ValueError(("graph6 byte", byte))
        for shift in range(5, -1, -1):
            bit = value >> shift & 1
            if bit_index < needed:
                if bit:
                    mask |= 1 << bit_index
            elif bit:
                raise ValueError(("nonzero graph6 padding", record))
            bit_index += 1
    if bit_index < needed:
        raise ValueError(("short graph6 record", record))
    return mask


def incidence_masks() -> tuple[int, ...]:
    return tuple(
        sum(1 << edge_index for edge_index, edge in enumerate(EDGE_PAIRS) if v in edge)
        for v in range(RESIDUAL_ORDER)
    )


def cut_masks() -> tuple[int, ...]:
    return tuple(
        sum(
            1 << edge_index
            for edge_index, (row, column) in enumerate(EDGE_PAIRS)
            if ((subset >> row) & 1) ^ ((subset >> column) & 1)
        )
        for subset in range(1 << RESIDUAL_ORDER)
    )


INCIDENCE_MASKS = incidence_masks()
CUT_MASKS = cut_masks()
SIZE_PRIORITY = (1, 9, 5, 4, 6, 3, 7, 2, 8, 0, 10)
SUBSET_ORDER = tuple(
    sorted(
        range(1 << RESIDUAL_ORDER),
        key=lambda subset: (SIZE_PRIORITY.index(subset.bit_count()), subset),
    )
)


def maximum_degree(edge_mask: int) -> int:
    return max((edge_mask & incident).bit_count() for incident in INCIDENCE_MASKS)


def violates_bound(edge_mask: int, bound: int = 15) -> tuple[bool, int]:
    edge_count = edge_mask.bit_count()
    checks = 0
    for subset in SUBSET_ORDER:
        subset_size = subset.bit_count()
        energy = (
            55
            - 2 * edge_count
            - 2 * subset_size * (11 - subset_size)
            + 4 * (edge_mask & CUT_MASKS[subset]).bit_count()
        )
        checks += 1
        if abs(energy) > bound:
            return True, checks
    return False, checks


@dataclass(frozen=True)
class StreamResult:
    records: int
    eligible_records: int
    subset_checks: int
    survivors: int
    independently_crosschecked_edge_counts: tuple[int, ...]
    sha256: str
    seconds: float


def certify_stream(geng: str, full_stream: bool) -> StreamResult:
    command = [geng, "-q", "10"]
    if not full_stream:
        command = [geng, "-q", "-D6", "10", "20:22"]
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout is not None
    digest = hashlib.sha256()
    records = 0
    eligible = 0
    subset_checks = 0
    survivors = 0
    independently_crosschecked_edge_counts: set[int] = set()
    started = time.monotonic()
    for raw_record in process.stdout:
        record = raw_record.strip()
        if not record:
            continue
        digest.update(record + b"\n")
        records += 1
        edge_mask = graph6_edge_mask(record, RESIDUAL_ORDER)
        edge_count = edge_mask.bit_count()
        if full_stream and not (
            20 <= edge_count <= 22 and maximum_degree(edge_mask) <= 6
        ):
            continue
        eligible += 1
        if edge_count not in independently_crosschecked_edge_counts:
            adjacency = small.graph6_adjacency(record, RESIDUAL_ORDER)
            for subset in range(1 << RESIDUAL_ORDER):
                size = subset.bit_count()
                cut = (edge_mask & CUT_MASKS[subset]).bit_count()
                mask_energy = (
                    55
                    - 2 * edge_count
                    - 2 * size * (11 - size)
                    + 4 * cut
                )
                if mask_energy != direct_energy(adjacency, subset):
                    raise AssertionError(
                        ("eligible sample energy disagreement", edge_count, subset)
                    )
            independently_crosschecked_edge_counts.add(edge_count)
        violated, checks = violates_bound(edge_mask)
        subset_checks += checks
        if not violated:
            survivors += 1
            raise AssertionError(("order-11 M<=15 survivor", record.decode("ascii")))

    stderr = process.stderr.read().decode("utf-8", errors="replace")
    return_code = process.wait()
    seconds = time.monotonic() - started
    if return_code:
        raise RuntimeError(("geng failure", return_code, stderr))

    expected_records = FULL_RECORDS if full_stream else FILTERED_RECORDS
    expected_digest = FULL_SHA256 if full_stream else FILTERED_SHA256
    if records != expected_records:
        raise AssertionError(("geng record count", records, expected_records))
    if eligible != FILTERED_RECORDS:
        raise AssertionError(("eligible record count", eligible, FILTERED_RECORDS))
    if digest.hexdigest() != expected_digest:
        raise AssertionError(("geng stream digest", digest.hexdigest(), expected_digest))
    if survivors:
        raise AssertionError(("unexpected survivors", survivors))
    if independently_crosschecked_edge_counts != {20, 21, 22}:
        raise AssertionError(
            (
                "eligible edge-count crosschecks",
                independently_crosschecked_edge_counts,
                {20, 21, 22},
            )
        )
    return StreamResult(
        records=records,
        eligible_records=eligible,
        subset_checks=subset_checks,
        survivors=survivors,
        independently_crosschecked_edge_counts=tuple(
            sorted(independently_crosschecked_edge_counts)
        ),
        sha256=digest.hexdigest(),
        seconds=seconds,
    )


def direct_energy(adjacency: tuple[int, ...], subset: int) -> int:
    order = len(adjacency) + 1
    spins = [1] + [
        -1 if subset >> vertex & 1 else 1 for vertex in range(order - 1)
    ]
    value = sum(spins[0] * spins[vertex] for vertex in range(1, order))
    for row in range(order - 1):
        for column in range(row + 1, order - 1):
            coefficient = -1 if adjacency[row] >> column & 1 else 1
            value += coefficient * spins[row + 1] * spins[column + 1]
    return value


def verify_witness(record: bytes, order: int, expected: int) -> tuple[int, int]:
    adjacency = small.graph6_adjacency(record, order - 1)
    histogram: Counter[int] = Counter()
    maximum = 0
    for subset in range(1 << (order - 1)):
        cut_energy = small.energy_from_cut(adjacency, subset)
        independent_energy = direct_energy(adjacency, subset)
        if cut_energy != independent_energy:
            raise AssertionError(
                ("witness energy disagreement", order, subset, cut_energy, independent_energy)
            )
        histogram[cut_energy] += 1
        maximum = max(maximum, abs(cut_energy))
    if maximum != expected:
        raise AssertionError(("witness maximum", order, maximum, expected))
    return maximum, sum(count for energy, count in histogram.items() if abs(energy) == maximum)


def verify_corruption_controls() -> None:
    # First check the local graph6 bit convention against the separately
    # implemented adjacency decoder and direct quadratic-form evaluation.
    adjacency = small.graph6_adjacency(ORDER11_WITNESS, 10)
    edge_mask = graph6_edge_mask(ORDER11_WITNESS, 10)
    for subset in range(1 << 10):
        size = subset.bit_count()
        cut = (edge_mask & CUT_MASKS[subset]).bit_count()
        energy = 55 - 2 * edge_mask.bit_count() - 2 * size * (11 - size) + 4 * cut
        if energy != direct_energy(adjacency, subset):
            raise AssertionError("graph6 edge ordering disagrees with direct evaluation")

    # Reversing the graph6 payload convention is an intentional corruption and
    # must be detected by the independent adjacency calculation.
    reversed_mask = sum(
        ((edge_mask >> index) & 1) << (EDGE_COUNT - 1 - index)
        for index in range(EDGE_COUNT)
    )
    reversed_mismatch = False
    for subset in range(1 << 10):
        size = subset.bit_count()
        cut = (reversed_mask & CUT_MASKS[subset]).bit_count()
        energy = 55 - 2 * reversed_mask.bit_count() - 2 * size * (11 - size) + 4 * cut
        if energy != direct_energy(adjacency, subset):
            reversed_mismatch = True
            break
    if not reversed_mismatch:
        raise AssertionError("reversed graph6 ordering corruption was not detected")

    false_bound_violated, _ = violates_bound(edge_mask, bound=15)
    if not false_bound_violated:
        raise AssertionError("false order-11 bound corruption was not detected")

    triangle_adjacency = (0b10, 0b01)
    values = [direct_energy(triangle_adjacency, subset) for subset in range(4)]
    if max(values) == max(abs(value) for value in values):
        raise AssertionError("absolute-value corruption was not detected")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--full-stream",
        action="store_true",
        help="read all 12,005,168 graphs instead of the exact reduced stream",
    )
    parser.add_argument("--geng")
    arguments = parser.parse_args()
    geng = small.locate_geng(arguments.geng)

    result = certify_stream(geng, arguments.full_stream)
    order11_maximum, order11_maximizers = verify_witness(ORDER11_WITNESS, 11, 17)
    order12_maximum, order12_maximizers = verify_witness(ORDER12_WITNESS, 12, 18)
    verify_corruption_controls()

    # F(11)>15 and odd energy parity give F(11)>=17; the first witness gives
    # the reverse inequality.  Monotonicity plus even order-12 energy parity
    # then gives F(12)>=18; the second witness gives equality.
    if math.comb(11, 2) % 2 != 1 or math.comb(12, 2) % 2 != 0:
        raise AssertionError("energy parity corruption")
    print(f"geng_mode={'full' if arguments.full_stream else 'filtered'}")
    print(f"geng_records={result.records}")
    print(f"eligible_order11_records={result.eligible_records}")
    print(f"order11_subset_checks={result.subset_checks}")
    print(f"order11_M_le_15_survivors={result.survivors}")
    print(
        "independently_crosschecked_edge_counts="
        + ",".join(map(str, result.independently_crosschecked_edge_counts))
    )
    print(f"geng_stream_sha256={result.sha256}")
    print(f"order11_witness={ORDER11_WITNESS.decode('ascii')} maximum={order11_maximum} maximizers={order11_maximizers}")
    print(f"order12_witness={ORDER12_WITNESS.decode('ascii')} maximum={order12_maximum} maximizers={order12_maximizers}")
    print("certified_values=F(11)=17,F(12)=18")
    print(f"deterministic_seed={DETERMINISTIC_SEED}")
    print("corruption_controls=graph6_ordering,absolute_value,false_bound,energy_parity")
    print("order11_order12_verification=PASSED")


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)

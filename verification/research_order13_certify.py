#!/usr/bin/env python3
"""Computer-assisted certificate for F(13)=20 and F(14)=21.

Every order-12 signing is switching-equivalent to a root-normalized signing,
whose residual negative edges form a graph on 11 vertices.  In full mode this
script runs nauty ``geng`` in eight disjoint ``res/mod`` shards, hashes and
counts the exact byte stream sent to a separately compiled C threshold
scanner, and checks both subprocess exit statuses.  The committed counts and
SHA-256 digests make truncation, overlap changes, and generator drift fail
closed.

The threshold scanner stops only after witnessing ``|Q| >= 19``.  Order-12
energies are even, so only signings with exact maximum at most 18 survive.
There are two rooted residual survivors.  A separate pure-Python evaluator
checks their complete spin cubes, all projective one-vertex extensions, and
all vertex deletions.  Both survivors have M=18 and extension minimum 24.
All other order-12 signings have M>=20 and therefore extension minimum at
least 20.  The Bellman identity gives F(13)>=20; a principal submatrix of the
order-14 Paley conference matrix gives equality.  Monotonicity, parity, and
the full conference witness then give F(14)=21.

Completeness trusts nauty to emit one representative of every unlabeled graph
in the documented shards.  No SAT/MILP status, floating point, cached scan
log, or supplied executable is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import check_conference_examples as conference


DETERMINISTIC_SEED = 413935
RESIDUAL_ORDER = 11
SIGNING_ORDER = 12
TOTAL_RECORDS = 1_018_997_864
SURVIVOR_RECORDS = ("JCpVdXyxpz?", "JCpdUg{[dM?")


@dataclass(frozen=True)
class ShardExpectation:
    records: int
    bytes: int
    sha256: str
    survivors: tuple[str, ...] = ()


SHARDS = (
    ShardExpectation(
        119_431_209,
        1_433_174_508,
        "0ae80a506fc9ef5aca212fb41eea05a453dd8326dc078fb1e2228b0c96cc2d5c",
    ),
    ShardExpectation(
        128_496_882,
        1_541_962_584,
        "54cd49e34942d1905213b802bff7adc565bd2c3251ac2268414ea179857559c5",
    ),
    ShardExpectation(
        128_472_053,
        1_541_664_636,
        "92ad1076355bd829ed0596dc70b9203dcffd88b5a4bc549de84aa37e6ceea318",
        (SURVIVOR_RECORDS[0],),
    ),
    ShardExpectation(
        121_592_284,
        1_459_107_408,
        "946fdda70e727892256a9e10c6e8d9823f2c1366f867499baa9a451b5d46e2e8",
    ),
    ShardExpectation(
        119_556_409,
        1_434_676_908,
        "b3dbbe8d186b4be67dbb119d41f34cb1637444bfa0207710549b9d04557a2798",
        (SURVIVOR_RECORDS[1],),
    ),
    ShardExpectation(
        134_239_743,
        1_610_876_916,
        "6acefb47e065175ccb85bb157d4cc50d6b1c8519d782146f2c449bc2582c2520",
    ),
    ShardExpectation(
        143_004_566,
        1_716_054_792,
        "6e2d97e700e524f244c14ddd1f0a832404b1661366bbebf4d5324640efc3d275",
    ),
    ShardExpectation(
        124_204_718,
        1_490_456_616,
        "73551d55bb414836abc29cf2700aee4552c9d476210cfba156434a9530dcc522",
    ),
)


@dataclass(frozen=True)
class ShardResult:
    shard: int
    records: int
    bytes: int
    sha256: str
    survivors: tuple[str, ...]


def graph6_signing(record: str) -> list[list[int]]:
    """Decode a short graph6 record and append an all-positive root."""
    if not record:
        raise ValueError("empty graph6 record")
    residual_order = ord(record[0]) - 63
    if residual_order != RESIDUAL_ORDER:
        raise ValueError(("graph6 order", residual_order, RESIDUAL_ORDER))
    edge_count = residual_order * (residual_order - 1) // 2
    data_characters = (edge_count + 5) // 6
    if len(record) != 1 + data_characters:
        raise ValueError(("graph6 length", len(record), 1 + data_characters))

    bits: list[int] = []
    for character in record[1:]:
        value = ord(character) - 63
        if not 0 <= value < 64:
            raise ValueError(("graph6 byte", character))
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))
    if any(bits[edge_count:]):
        raise ValueError("nonzero graph6 padding")

    matrix = [[0] * (residual_order + 1) for _ in range(residual_order + 1)]
    bit_index = 0
    for column in range(1, residual_order):
        for row in range(column):
            coefficient = -1 if bits[bit_index] else 1
            matrix[row][column] = matrix[column][row] = coefficient
            bit_index += 1
    for vertex in range(residual_order):
        matrix[vertex][residual_order] = matrix[residual_order][vertex] = 1
    return matrix


def projective_spins(order: int) -> list[tuple[int, ...]]:
    return [
        (1,)
        + tuple(
            -1 if mask >> (vertex - 1) & 1 else 1
            for vertex in range(1, order)
        )
        for mask in range(1 << (order - 1))
    ]


def energy(matrix: list[list[int]], spins: tuple[int, ...]) -> int:
    return sum(
        matrix[row][column] * spins[row] * spins[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def exact_energies(matrix: list[list[int]]) -> tuple[list[int], list[tuple[int, ...]]]:
    spins = projective_spins(len(matrix))
    return [energy(matrix, state) for state in spins], spins


def maximum_absolute_energy(matrix: list[list[int]]) -> int:
    values, _ = exact_energies(matrix)
    return max(abs(value) for value in values)


def extension_minimum(matrix: list[list[int]]) -> tuple[int, int]:
    """Directly enumerate every projective incident column."""
    values, spins = exact_energies(matrix)
    best = 10**9
    optimal_centers = 0
    for column in projective_spins(len(matrix)):
        candidate = max(
            abs(value)
            + abs(sum(coefficient * spin for coefficient, spin in zip(column, state)))
            for value, state in zip(values, spins)
        )
        if candidate < best:
            best = candidate
            optimal_centers = 1
        elif candidate == best:
            optimal_centers += 1
    return best, optimal_centers


def delete_vertex(matrix: list[list[int]], deleted: int) -> list[list[int]]:
    retained = [index for index in range(len(matrix)) if index != deleted]
    return [[matrix[row][column] for column in retained] for row in retained]


def compile_scanner(compiler: str, output: Path) -> None:
    source = Path(__file__).with_name("order12_threshold_scan.c")
    command = [
        compiler,
        "-std=c11",
        "-O3",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        "-Wconversion",
        "-Wshadow",
        "-Werror",
        str(source),
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True)


def parse_scanner_output(output: str) -> tuple[int, tuple[str, ...]]:
    survivors: list[str] = []
    summary: tuple[int, int] | None = None
    for line in output.splitlines():
        survivor_match = re.fullmatch(r"SURVIVOR (\S+) M=(\d+)", line)
        if survivor_match:
            record, maximum = survivor_match.groups()
            if int(maximum) > 18:
                raise AssertionError(("invalid survivor maximum", line))
            survivors.append(record)
            continue
        summary_match = re.fullmatch(r"SUMMARY records=(\d+) survivors=(\d+)", line)
        if summary_match:
            if summary is not None:
                raise AssertionError(("duplicate summary", output))
            summary = tuple(map(int, summary_match.groups()))
            continue
        raise AssertionError(("unparsed scanner output", line))
    if summary is None or summary[1] != len(survivors):
        raise AssertionError(("scanner summary", summary, survivors))
    return summary[0], tuple(survivors)


def verify_quick(scanner: Path) -> None:
    empty_graph = "J" + "?" * 10
    fixture = (empty_graph,) + SURVIVOR_RECORDS
    process = subprocess.run(
        [str(scanner)],
        input=("\n".join(fixture) + "\n").encode("ascii"),
        capture_output=True,
        check=True,
    )
    records, c_survivors = parse_scanner_output(process.stdout.decode("ascii"))
    if process.stderr or records != len(fixture) or c_survivors != SURVIVOR_RECORDS:
        raise AssertionError(("C scanner fixture", process.stderr, records, c_survivors))
    if maximum_absolute_energy(graph6_signing(empty_graph)) != 66:
        raise AssertionError("empty residual corruption control")

    malformed = SURVIVOR_RECORDS[0][:-1] + "@"
    malformed_run = subprocess.run(
        [str(scanner)], input=(malformed + "\n").encode("ascii"), capture_output=True
    )
    if malformed_run.returncode == 0:
        raise AssertionError("nonzero graph6 padding was accepted")

    for record in SURVIVOR_RECORDS:
        matrix = graph6_signing(record)
        values, _ = exact_energies(matrix)
        maximum = max(abs(value) for value in values)
        maximizers = sum(abs(value) == maximum for value in values)
        extension, optimal_centers = extension_minimum(matrix)
        deletion_maxima = tuple(
            maximum_absolute_energy(delete_vertex(matrix, vertex))
            for vertex in range(SIGNING_ORDER)
        )
        if (maximum, maximizers, extension, optimal_centers) != (18, 20, 24, 772):
            raise AssertionError(
                (record, maximum, maximizers, extension, optimal_centers)
            )
        if deletion_maxima != (17,) * SIGNING_ORDER:
            raise AssertionError((record, "deletion maxima", deletion_maxima))
        print(
            f"order12_survivor={record} M={maximum} maximizers={maximizers} "
            f"extension_minimum={extension} optimal_centers={optimal_centers} "
            "deletion_maxima=17x12"
        )

    corrupted = graph6_signing(SURVIVOR_RECORDS[0])
    corrupted[0][1] *= -1
    corrupted[1][0] *= -1
    if maximum_absolute_energy(corrupted) != 20:
        raise AssertionError("survivor edge-flip corruption was not detected")

    conference_matrix = conference.paley_conference(13)
    conference.check_conference_identity(conference_matrix)
    conference_maximum = conference.maximum_absolute_energy(conference_matrix)
    principal_maxima = tuple(
        conference.maximum_absolute_energy(
            conference.principal_submatrix(conference_matrix, deleted)
        )
        for deleted in range(14)
    )
    if conference_maximum != 21 or principal_maxima != (20,) * 14:
        raise AssertionError((conference_maximum, principal_maxima))
    print("paley_C14_M=21 principal_order13_maxima=20x14")
    print("corruption_controls=graph6_padding,empty_residual,edge_flip")
    print("order13_order14_quick_verification=PASSED")


def run_shard(geng: str, scanner: Path, shard: int) -> ShardResult:
    expectation = SHARDS[shard]
    with tempfile.TemporaryFile() as generator_stderr, tempfile.TemporaryFile() as scanner_stdout, tempfile.TemporaryFile() as scanner_stderr:
        generator = subprocess.Popen(
            [geng, "-q", str(RESIDUAL_ORDER), f"{shard}/8"],
            stdout=subprocess.PIPE,
            stderr=generator_stderr,
        )
        scan = subprocess.Popen(
            [str(scanner)],
            stdin=subprocess.PIPE,
            stdout=scanner_stdout,
            stderr=scanner_stderr,
        )
        assert generator.stdout is not None and scan.stdin is not None
        digest = hashlib.sha256()
        records = 0
        byte_count = 0
        try:
            while True:
                chunk = generator.stdout.read(1 << 20)
                if not chunk:
                    break
                digest.update(chunk)
                records += chunk.count(b"\n")
                byte_count += len(chunk)
                scan.stdin.write(chunk)
            generator.stdout.close()
            scan.stdin.close()
            generator_returncode = generator.wait()
            scanner_returncode = scan.wait()
        except BaseException:
            generator.kill()
            scan.kill()
            generator.wait()
            scan.wait()
            raise

        generator_stderr.seek(0)
        scanner_stdout.seek(0)
        scanner_stderr.seek(0)
        generator_error = generator_stderr.read().decode("utf-8", errors="replace")
        scanner_output = scanner_stdout.read().decode("ascii")
        scanner_error = scanner_stderr.read().decode("utf-8", errors="replace")

    if generator_returncode or scanner_returncode or generator_error or scanner_error:
        raise RuntimeError(
            {
                "shard": shard,
                "generator_returncode": generator_returncode,
                "scanner_returncode": scanner_returncode,
                "generator_stderr": generator_error,
                "scanner_stderr": scanner_error,
            }
        )
    scanner_records, survivors = parse_scanner_output(scanner_output)
    result = ShardResult(shard, records, byte_count, digest.hexdigest(), survivors)
    if scanner_records != records:
        raise AssertionError(("producer/consumer record mismatch", result, scanner_records))
    if (
        result.records != expectation.records
        or result.bytes != expectation.bytes
        or result.sha256 != expectation.sha256
        or result.survivors != expectation.survivors
    ):
        raise AssertionError(("shard receipt mismatch", result, expectation))
    return result


def verify_full(geng: str, scanner: Path, jobs: int) -> None:
    results: list[ShardResult] = []
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = [
            executor.submit(run_shard, geng, scanner, shard) for shard in range(8)
        ]
        for future in as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda result: result.shard)

    for result in results:
        print(
            json.dumps(
                {
                    "bytes": result.bytes,
                    "records": result.records,
                    "sha256": result.sha256,
                    "shard": result.shard,
                    "survivors": list(result.survivors),
                },
                sort_keys=True,
            )
        )
    total = sum(result.records for result in results)
    survivors = tuple(
        record for result in results for record in result.survivors
    )
    if total != TOTAL_RECORDS or survivors != SURVIVOR_RECORDS:
        raise AssertionError(("full scan aggregate", total, survivors))
    print(f"order12_residual_records={total}")
    print("order12_threshold_survivors=" + ",".join(survivors))
    print("certified_values=F(13)=20,F(14)=21")
    print("certificate_boundary=nauty_unlabeled_graph_completeness")
    print("order13_order14_full_verification=PASSED")


def resolve_executable(requested: str | None, alternatives: tuple[str, ...]) -> str:
    if requested is not None:
        path = shutil.which(requested) if os.sep not in requested else requested
        if path is None or not os.access(path, os.X_OK):
            raise FileNotFoundError(requested)
        return str(Path(path).resolve())
    for alternative in alternatives:
        path = shutil.which(alternative)
        if path is not None:
            return str(Path(path).resolve())
    raise FileNotFoundError(" or ".join(alternatives))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-stream", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--geng")
    parser.add_argument("--cc")
    arguments = parser.parse_args()
    if arguments.jobs < 1 or arguments.jobs > 8:
        parser.error("--jobs must be between 1 and 8")

    compiler = resolve_executable(arguments.cc, ("cc", "gcc", "clang"))
    with tempfile.TemporaryDirectory(prefix="mo413935-order13-") as directory:
        scanner = Path(directory) / "order12_threshold_scan"
        compile_scanner(compiler, scanner)
        verify_quick(scanner)
        if arguments.full_stream:
            geng = resolve_executable(arguments.geng, ("geng", "nauty-geng"))
            verify_full(geng, scanner, arguments.jobs)
    print(f"deterministic_seed={DETERMINISTIC_SEED}")


if __name__ == "__main__":
    main()

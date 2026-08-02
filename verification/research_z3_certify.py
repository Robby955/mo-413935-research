#!/usr/bin/env python3
"""Independent Z3 certificates for the exact values F(7), F(8), and F(9).

After switching, every root edge is +1.  The remaining coefficients are Z3
integer variables constrained to {-1,+1}.  For every one of the 2^(n-1)
spin vectors with root spin +1, the script asserts the two exact linear
inequalities ``-T <= Q_A(x) <= T``.

For each n it proves UNSAT at F(n)-2 and SAT at F(n).  Quadratic energies
have the parity of binom(n,2), so these two checks certify the exact optimum.
The SAT model is then evaluated independently with plain Python integer
arithmetic over every spin vector.

Expected status output (timings and witness hashes can vary with Z3 version):

    {"n": 7, "sat_bound": 9,  "unsat_bound": 7,  "status": "PASSED"}
    {"n": 8, "sat_bound": 10, "unsat_bound": 8,  "status": "PASSED"}
    {"n": 9, "sat_bound": 12, "unsat_bound": 10, "status": "PASSED"}
    {"corruption_controls": "PASSED", "deterministic_seed": 413935}

The default n=7..9 run took about 90 seconds with Z3 4.16.0 on an Apple M4.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import time
from typing import Any


EXPECTED_F = {7: 9, 8: 10, 9: 12}
DETERMINISTIC_SEED = 413935


def load_z3() -> Any:
    try:
        import z3
    except ImportError as error:
        raise RuntimeError(
            "z3-solver is required for this optional independent certificate; "
            "install z3-solver and rerun"
        ) from error
    return z3


def residual_edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(itertools.combinations(range(1, n), 2))


def spins_from_mask(n: int, mask: int) -> tuple[int, ...]:
    return (1,) + tuple(-1 if mask >> (vertex - 1) & 1 else 1 for vertex in range(1, n))


def direct_energy(
    n: int, coefficients: dict[tuple[int, int], int], spin_mask: int
) -> int:
    spins = spins_from_mask(n, spin_mask)
    root_edges = sum(spins[vertex] for vertex in range(1, n))
    residual = sum(
        coefficients[(row, column)] * spins[row] * spins[column]
        for row, column in residual_edges(n)
    )
    return root_edges + residual


def validate_witness(
    n: int, coefficients: dict[tuple[int, int], int], claimed_bound: int
) -> int:
    expected_edges = set(residual_edges(n))
    if set(coefficients) != expected_edges:
        raise ValueError("witness has the wrong coefficient keys")
    if any(value not in (-1, 1) for value in coefficients.values()):
        raise ValueError("witness contains a coefficient outside {-1,+1}")
    maximum = max(
        abs(direct_energy(n, coefficients, spin_mask))
        for spin_mask in range(1 << (n - 1))
    )
    if maximum > claimed_bound:
        raise AssertionError(
            ("witness violates claimed bound", n, maximum, claimed_bound)
        )
    return maximum


def build_solver(z3: Any, n: int, bound: int, timeout_ms: int):
    solver = z3.Solver()
    solver.set(timeout=timeout_ms, random_seed=DETERMINISTIC_SEED)
    edges = residual_edges(n)
    variables = {edge: z3.Int(f"a_{n}_{edge[0]}_{edge[1]}") for edge in edges}
    for variable in variables.values():
        solver.add(z3.Or(variable == 1, variable == -1))

    # Negating the full signing preserves max |Q|.  After renormalizing its
    # root edges to +1, every residual coefficient is negated.  We may
    # therefore choose the member of this pair with at most half of the
    # residual coefficients negative.  This is exact symmetry breaking.
    solver.add(
        z3.Sum([z3.If(variable == -1, 1, 0) for variable in variables.values()])
        <= len(edges) // 2
    )

    for spin_mask in range(1 << (n - 1)):
        spins = spins_from_mask(n, spin_mask)
        energy = sum(spins[vertex] for vertex in range(1, n)) + sum(
            variables[(row, column)] * spins[row] * spins[column]
            for row, column in edges
        )
        solver.add(energy <= bound, energy >= -bound)
    return solver, variables


def certify_order(z3: Any, n: int, target: int, timeout_ms: int) -> dict[str, object]:
    lower_bound = target - 2
    started = time.monotonic()
    lower_solver, _ = build_solver(z3, n, lower_bound, timeout_ms)
    lower_status = lower_solver.check()
    lower_seconds = time.monotonic() - started
    if lower_status != z3.unsat:
        reason = lower_solver.reason_unknown() if lower_status == z3.unknown else ""
        raise AssertionError(
            ("lower certificate was not UNSAT", n, lower_bound, lower_status, reason)
        )

    started = time.monotonic()
    upper_solver, variables = build_solver(z3, n, target, timeout_ms)
    upper_status = upper_solver.check()
    upper_seconds = time.monotonic() - started
    if upper_status != z3.sat:
        reason = upper_solver.reason_unknown() if upper_status == z3.unknown else ""
        raise AssertionError(
            ("upper certificate was not SAT", n, target, upper_status, reason)
        )

    model = upper_solver.model()
    coefficients = {
        edge: model.eval(variable).as_long() for edge, variable in variables.items()
    }
    witnessed_maximum = validate_witness(n, coefficients, target)
    if witnessed_maximum != target:
        raise AssertionError(
            ("unexpectedly non-sharp SAT witness", n, witnessed_maximum)
        )

    ordered_coefficients = [coefficients[edge] for edge in residual_edges(n)]
    witness_sha256 = hashlib.sha256(
        json.dumps(ordered_coefficients, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    return {
        "n": n,
        "sat_bound": target,
        "sat_seconds": round(upper_seconds, 3),
        "status": "PASSED",
        "unsat_bound": lower_bound,
        "unsat_seconds": round(lower_seconds, 3),
        "witness_maximum": witnessed_maximum,
        "witness_sha256": witness_sha256,
        "_witness": coefficients,
    }


def verify_corruption_controls(certificates: list[dict[str, object]]) -> None:
    if not certificates:
        raise AssertionError("no certificates supplied to corruption controls")
    certificate = certificates[0]
    n = int(certificate["n"])
    target = int(certificate["sat_bound"])
    coefficients = dict(certificate["_witness"])

    first_edge = residual_edges(n)[0]
    malformed = dict(coefficients)
    malformed[first_edge] = 0
    try:
        validate_witness(n, malformed, target)
    except ValueError:
        pass
    else:
        raise AssertionError("zero-coefficient corruption was not detected")

    try:
        validate_witness(n, coefficients, target - 2)
    except AssertionError:
        pass
    else:
        raise AssertionError("false-bound corruption was not detected")

    # A frustrated triangle has one-sided energies (1,1,1,-3).  This catches
    # the tempting corruption that removes the absolute value.
    triangle = {(1, 2): -1}
    triangle_energies = [direct_energy(3, triangle, mask) for mask in range(4)]
    if max(triangle_energies) == max(abs(value) for value in triangle_energies):
        raise AssertionError("absolute-value corruption was not detected")


def printable_certificate(certificate: dict[str, object]) -> dict[str, object]:
    return {key: value for key, value in certificate.items() if key != "_witness"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=7)
    parser.add_argument("--max-n", type=int, default=9)
    parser.add_argument("--timeout-ms", type=int, default=180_000)
    arguments = parser.parse_args()
    if not 7 <= arguments.min_n <= arguments.max_n <= 9:
        raise SystemExit("require 7 <= min-n <= max-n <= 9")
    if arguments.timeout_ms <= 0:
        raise SystemExit("timeout-ms must be positive")

    try:
        z3 = load_z3()
    except RuntimeError as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(2) from error

    certificates: list[dict[str, object]] = []
    for n in range(arguments.min_n, arguments.max_n + 1):
        certificate = certify_order(z3, n, EXPECTED_F[n], arguments.timeout_ms)
        certificates.append(certificate)
        print(
            json.dumps(printable_certificate(certificate), sort_keys=True), flush=True
        )

    verify_corruption_controls(certificates)
    print(
        json.dumps(
            {
                "corruption_controls": "PASSED",
                "deterministic_seed": DETERMINISTIC_SEED,
                "z3_version": z3.get_version_string(),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

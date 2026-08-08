"""Rigorous upper bounds on F(n) = min_A max_x |sum_{i<j} a_ij x_i x_j|.

Inner maximum is EXACT (full enumeration of 2^(n-1) sign vectors, x_n fixed to +1
since Q(-x) = Q(x)).  Outer minimisation over A is simulated annealing, so every
number produced is a valid upper bound F(n) <= M(A) for the witness A printed.

Incremental trick: flipping a_ij changes every energy by -2 a_ij x_i x_j, so one
annealing step costs O(2^(n-1)) instead of O(2^(n-1) n^2).
"""

import argparse
import json
import math
import sys
import time

import numpy as np


def columns(n):
    """col[i] = +-1 pattern of coordinate i over all 2^(n-1) sign vectors."""
    rows = 1 << (n - 1)
    idx = np.arange(rows, dtype=np.int64)
    cols = np.empty((n, rows), dtype=np.int8)
    for i in range(n - 1):
        cols[i] = np.where((idx >> i) & 1, -1, 1).astype(np.int8)
    cols[n - 1] = 1
    return cols


def energies(a, cols, n):
    """E[x] = sum_{i<j} a_ij x_i x_j, exactly, as int32."""
    rows = cols.shape[1]
    e = np.zeros(rows, dtype=np.int32)
    for i in range(n):
        for j in range(i + 1, n):
            if a[i, j] > 0:
                e += cols[i] * cols[j]
            else:
                e -= cols[i] * cols[j]
    return e


def random_matrix(n, rng):
    a = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(i + 1, n):
            v = 1 if rng.random() < 0.5 else -1
            a[i, j] = v
            a[j, i] = v
    return a


def paley_conference(n):
    """Symmetric conference matrix of order n = q+1, q prime, q = 1 mod 4."""
    q = n - 1
    if q < 5 or q % 4 != 1:
        return None
    for d in range(2, int(q ** 0.5) + 1):
        if q % d == 0:
            return None
    chi = np.zeros(q, dtype=np.int8)
    res = {(k * k) % q for k in range(1, q)}
    for k in range(1, q):
        chi[k] = 1 if k in res else -1
    c = np.zeros((n, n), dtype=np.int8)
    for i in range(1, n):
        c[0, i] = 1
        c[i, 0] = 1
    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                c[i, j] = chi[(i - j) % q]
    for i in range(n):
        c[i, i] = 0
    # off-diagonal zeros cannot occur for i != j here; force +-1
    c[c == 0] = 1
    np.fill_diagonal(c, 0)
    c = np.triu(c) + np.triu(c, 1).T
    return c.astype(np.int8)


def anneal(n, steps, seed, t0, t1, init="random", verbose=False):
    rng = np.random.default_rng(seed)
    cols = columns(n)
    a = None
    if init == "conference":
        a = paley_conference(n)
    if a is None:
        a = random_matrix(n, rng)
    e = energies(a, cols, n)
    cur = int(np.abs(e).max())
    best = cur
    best_a = a.copy()
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    npairs = len(pairs)
    t_start = time.time()
    for step in range(steps):
        temp = t0 * (t1 / t0) ** (step / max(1, steps - 1))
        pi, pj = pairs[rng.integers(npairs)]
        delta = np.int32(-2 * a[pi, pj])
        prod = (cols[pi] * cols[pj]).astype(np.int32)
        e += delta * prod
        new = int(np.abs(e).max())
        if new <= cur or rng.random() < math.exp((cur - new) / temp):
            cur = new
            a[pi, pj] = -a[pi, pj]
            a[pj, pi] = a[pi, pj]
            if new < best:
                best = new
                best_a = a.copy()
                if verbose:
                    print(f"  n={n} step={step} best={best} "
                          f"ratio={best / n ** 1.5:.4f}", flush=True)
        else:
            e -= delta * prod
    return best, best_a, time.time() - t_start


def verify(a, n):
    """Recompute M(A) from scratch, independently of the annealing state."""
    cols = columns(n)
    e = energies(a, cols, n)
    return int(np.abs(e).max())


def basin_hop(n, rounds, steps, seed, kick, t0, t1, cols, a0=None, verbose=False):
    """Iterated local search: kick the incumbent, re-anneal, keep improvements."""
    rng = np.random.default_rng(seed)
    a = random_matrix(n, rng) if a0 is None else a0.copy()
    e = energies(a, cols, n)
    cur = int(np.abs(e).max())
    best, best_a = cur, a.copy()
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    npairs = len(pairs)
    for rd in range(rounds):
        a = best_a.copy()
        for _ in range(kick):
            pi, pj = pairs[rng.integers(npairs)]
            a[pi, pj] = -a[pi, pj]
            a[pj, pi] = a[pi, pj]
        e = energies(a, cols, n)
        cur = int(np.abs(e).max())
        for step in range(steps):
            temp = t0 * (t1 / t0) ** (step / max(1, steps - 1))
            pi, pj = pairs[rng.integers(npairs)]
            delta = np.int32(-2 * a[pi, pj])
            prod = (cols[pi] * cols[pj]).astype(np.int32)
            e += delta * prod
            new = int(np.abs(e).max())
            if new <= cur or rng.random() < math.exp((cur - new) / temp):
                cur = new
                a[pi, pj] = -a[pi, pj]
                a[pj, pi] = a[pi, pj]
                if new < best:
                    best, best_a = new, a.copy()
                    if verbose:
                        print(f"  n={n} round={rd} best={best} "
                              f"ratio={best / n ** 1.5:.4f}", flush=True)
            else:
                e -= delta * prod
    return best, best_a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--rounds", type=int, default=40)
    ap.add_argument("--kick", type=int, default=6)
    ap.add_argument("--restarts", type=int, default=4)
    ap.add_argument("--t0", type=float, default=3.0)
    ap.add_argument("--t1", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=413935)
    args = ap.parse_args()

    n = args.n
    overall = None
    overall_a = None
    cols = columns(n)
    conf = paley_conference(n)
    if conf is not None:
        m_conf = verify(conf, n)
        print(f"n={n} paley_conference M={m_conf} "
              f"ceiling={0.5 * n * math.sqrt(n - 1):.3f}", flush=True)
    for r in range(args.restarts):
        t_start = time.time()
        seed_a = conf if (r == 0 and conf is not None) else None
        best, best_a = basin_hop(
            n, args.rounds, args.steps, args.seed + 1000 * r, args.kick,
            args.t0, args.t1, cols, a0=seed_a, verbose=True)
        checked = verify(best_a, n)
        assert checked == best, f"verification mismatch {checked} != {best}"
        print(f"n={n} restart={r} M={best} ratio={best / n ** 1.5:.4f} "
              f"secs={time.time() - t_start:.1f}", flush=True)
        if overall is None or best < overall:
            overall = best
            overall_a = best_a
    ceiling = 0.5 * n * math.sqrt(n - 1)
    out = {
        "n": n,
        "F_upper_bound": overall,
        "ratio_n32": overall / n ** 1.5,
        "spectral_ceiling": ceiling,
        "fraction_of_ceiling": overall / ceiling,
        "steps": args.steps,
        "restarts": args.restarts,
        "witness_upper_triangle": [
            int(overall_a[i, j]) for i in range(n) for j in range(i + 1, n)
        ],
        "witness_M_reverified": verify(overall_a, n),
        "spectrum_absmax": float(
            np.abs(np.linalg.eigvalsh(overall_a.astype(float))).max()),
        "sqrt_n_minus_1": math.sqrt(n - 1),
    }
    print("RESULT " + json.dumps(out), flush=True)


if __name__ == "__main__":
    sys.exit(main())

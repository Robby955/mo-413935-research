"""Maximum of a Paley conference matrix against its spectral ceiling.

For a symmetric conference matrix C of order n, C^2 = (n-1)I, so every
eigenvalue is +-sqrt(n-1) and Cauchy-Schwarz gives

    M(C) <= n*sqrt(n-1)/2,

the ceiling that produces limsup F(n)/n^(3/2) <= 1/2.  The ceiling is attained
only if some sign vector aligns fully with the top eigenspace.  This script
measures how close the Paley family comes.

Exact by full enumeration for n <= 26.  Beyond that the value reported is a
multi-start local-search maximum, hence a LOWER bound on M(C).  That is the
useful direction here: it can only understate how close Paley gets to its
ceiling, so an observed ratio near 1 is evidence that cannot be an artefact of
an incomplete search.

Usage:
    python3 verification/paley_scan.py [max_q]
"""

import math
import sys

import numpy as np


def is_prime(q):
    if q < 2:
        return False
    for d in range(2, int(q ** 0.5) + 1):
        if q % d == 0:
            return False
    return True


def paley_conference(q):
    """Symmetric conference matrix of order q+1, for prime q = 1 mod 4."""
    residues = {(k * k) % q for k in range(1, q)}
    chi = np.zeros(q, dtype=np.int64)
    for k in range(1, q):
        chi[k] = 1 if k in residues else -1
    n = q + 1
    C = np.ones((n, n), dtype=np.int64)
    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                C[i, j] = chi[(i - j) % q]
    np.fill_diagonal(C, 0)
    return C


def exact_max(C):
    """Exact M(C) by enumerating all 2^(n-1) sign vectors, x_n fixed to +1."""
    n = C.shape[0]
    rows = 1 << (n - 1)
    idx = np.arange(rows, dtype=np.int64)
    best = 0
    chunk = 1 << 22
    Cf = C.astype(np.float32)
    for s in range(0, rows, chunk):
        ii = idx[s:s + chunk]
        X = np.empty((len(ii), n), dtype=np.float32)
        for i in range(n - 1):
            X[:, i] = np.where((ii >> i) & 1, -1.0, 1.0)
        X[:, n - 1] = 1.0
        E = np.einsum("ij,ij->i", X @ Cf, X) / 2.0
        best = max(best, int(np.abs(E).max()))
    return best


def search_max(C, restarts=4000, seed=0, iters=300):
    """Multi-start steepest-ascent maximum; a lower bound on M(C)."""
    n = C.shape[0]
    rng = np.random.default_rng(seed)
    Cf = C.astype(np.float32)
    best = 0.0
    for sign in (1.0, -1.0):
        X = rng.choice(np.array([-1, 1], dtype=np.float32), size=(restarts, n))
        for _ in range(iters):
            gain = -2.0 * sign * X * (X @ Cf)
            k = np.argmax(gain, axis=1)
            g = gain[np.arange(restarts), k]
            act = g > 0
            if not act.any():
                break
            X[np.arange(restarts)[act], k[act]] *= -1
        v = sign * np.einsum("ij,ij->i", X @ Cf, X) / 2.0
        best = max(best, float(v.max()))
    return int(best)


def main(max_q):
    print(f"{'q':>5} {'n':>5} {'M':>9}  {'kind':<16} {'ceiling':>10} "
          f"{'M/ceiling':>10} {'M/n^1.5':>9}")
    for q in range(5, max_q + 1, 4):
        if not is_prime(q):
            continue
        C = paley_conference(q)
        n = q + 1
        assert np.array_equal(C @ C, (n - 1) * np.eye(n, dtype=np.int64)), \
            f"not a conference matrix at q={q}"
        ceiling = 0.5 * n * math.sqrt(n - 1)
        if n <= 26:
            m, kind = exact_max(C), "exact"
        else:
            m, kind = search_max(C, seed=q), "search lower bd"
        print(f"{q:>5} {n:>5} {m:>9}  {kind:<16} {ceiling:>10.2f} "
              f"{m / ceiling:>10.4f} {m / n ** 1.5:>9.4f}", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 113)

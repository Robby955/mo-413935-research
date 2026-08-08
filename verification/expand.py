"""One level of the tower: order n -> order n+1, keeping every extension with M <= T.

Input  : graph6 file of graphs on m = n-1 vertices (order-n matrices, last row +1)
Output : graph6 file of graphs on m+1 = n vertices (order-(n+1) matrices, last row +1)

Extension rule.  Fix x_{n+1} = +1 (legitimate, E is quadratic).  With v in {-1,1}^n
the sign row of the new vertex,

    E(x) = E_B(x) + <v, x>,   x in {-1,1}^n,

so one (2^n x 2^n) block gives M for every v at once.

Normalising the extension.  Let w_i = [v_i = -1] xor [v_n = -1] for i = 1..n, so
w_n = 0.  Switching by v puts the new vertex's row at all +1, and the resulting
graph on n vertices is: G Seidel-switched on W = {i : w_i = 1}, plus a new vertex
joined to exactly W.  v and -v give the same graph, so 2^n rows give 2^(n-1) graphs.

Usage:  python3 expand.py in.g6 m thresh out.g6 [--res R --mod M]
"""

import argparse
import subprocess
import sys
import time

import numpy as np


def pair_order(m):
    return [(i, j) for j in range(1, m) for i in range(j)]


def g6_decode(line, m):
    npairs = m * (m - 1) // 2
    nb = (npairs + 5) // 6
    b = line.strip().encode()
    assert b[0] == m + 63, (b[0], m)
    body = np.frombuffer(b[1:1 + nb], dtype=np.uint8) - 63
    bits = np.unpackbits(body.astype(np.uint8)).reshape(nb, 8)[:, 2:8]
    return bits.reshape(-1)[:npairs].astype(np.uint8)


def g6_encode_many(E, m):
    """E: (K, npairs) uint8 bit array -> list of graph6 strings on m vertices."""
    npairs = m * (m - 1) // 2
    nb = (npairs + 5) // 6
    pad = nb * 6 - npairs
    if pad:
        E = np.concatenate([E, np.zeros((E.shape[0], pad), dtype=np.uint8)], axis=1)
    E = E.reshape(E.shape[0], nb, 6)
    two = np.zeros((E.shape[0], nb, 2), dtype=np.uint8)
    packed = np.packbits(np.concatenate([two, E], axis=2), axis=2)[:, :, 0]
    packed = packed + 63
    head = np.full((E.shape[0], 1), m + 63, dtype=np.uint8)
    out = np.concatenate([head, packed], axis=1)
    return [bytes(r).decode() for r in out]


def build_shared(n):
    """n = order of the SOURCE matrix.  Sign vectors over its n vertices."""
    rows = 1 << n
    idx = np.arange(rows, dtype=np.int64)
    X = np.empty((rows, n), dtype=np.float32)
    for i in range(n):
        X[:, i] = np.where((idx >> i) & 1, -1.0, 1.0)
    m = n - 1
    pairs = pair_order(m)
    Q = np.empty((rows, len(pairs)), dtype=np.float32)
    for k, (i, j) in enumerate(pairs):
        Q[:, k] = X[:, i] * X[:, j]
    base = Q.sum(axis=1) + (X[:, :m] * X[:, [m]]).sum(axis=1)
    return X, Q, base.astype(np.float32)


def new_edges(e, wbits, m):
    """e: (npairs,) uint8 for G on m vertices;  wbits: (K, n) uint8 with n=m+1.
    Returns (K, npairs_new) uint8 for the graphs on n = m+1 vertices."""
    n = m + 1
    pairs_old = pair_order(m)
    pairs_new = pair_order(n)
    K = wbits.shape[0]
    out = np.empty((K, len(pairs_new)), dtype=np.uint8)
    old_index = {p: k for k, p in enumerate(pairs_old)}
    for k, (i, j) in enumerate(pairs_new):
        if j < m:
            out[:, k] = e[old_index[(i, j)]] ^ wbits[:, i] ^ wbits[:, j]
        else:                       # j == m, the freshly added vertex
            out[:, k] = wbits[:, i]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("m", type=int, help="vertices of the input graphs")
    ap.add_argument("thresh", type=float)
    ap.add_argument("outfile")
    ap.add_argument("--res", type=int, default=0)
    ap.add_argument("--mod", type=int, default=1)
    ap.add_argument("--count-only", action="store_true")
    a = ap.parse_args()

    m = a.m
    n = m + 1                                   # order of the source matrix
    X, Q, base = build_shared(n)
    D = X @ X.T                                 # (2^n, 2^n) inner products <v,x>
    lines = [l.strip() for l in open(a.infile) if l.strip()]
    lines = lines[a.res::a.mod]
    half = 1 << (n - 1)                         # v and -v are equivalent

    fh = open(a.outfile, "w")
    kept = 0
    t0 = time.time()
    for gi, line in enumerate(lines):
        e = g6_decode(line, m)
        EB = base - 2.0 * (Q @ e.astype(np.float32))          # (2^n,)
        Mv = np.abs(EB[:, None] + D).max(axis=0)              # (2^n,) per v
        good = np.flatnonzero(Mv[:half] <= a.thresh)          # v with v_n = +1
        if len(good):
            vb = np.empty((len(good), n), dtype=np.uint8)
            for i in range(n):
                vb[:, i] = ((good >> i) & 1).astype(np.uint8)
            wb = vb ^ vb[:, [n - 1]]
            E = new_edges(e, wb, m)
            if not a.count_only:
                fh.write("\n".join(g6_encode_many(E, n)) + "\n")
            kept += len(good)
        if gi and gi % 200 == 0:
            print(f"  {gi}/{len(lines)} raw_kept={kept} "
                  f"{(time.time()-t0)/gi*1000:.1f} ms/graph", flush=True)
    fh.close()
    dt = time.time() - t0
    print(f"source m={m} (order {n}) graphs={len(lines)} thresh={int(a.thresh)} "
          f"raw_extensions_kept={kept} {dt:.1f}s "
          f"({dt/max(len(lines),1)*1000:.1f} ms/graph)", flush=True)


if __name__ == "__main__":
    main()

"""expand.py for the upper levels: same mathematics, blocked over v to keep the
peak memory near one column block instead of a full (2^n x 2^n) temporary.

Usage: python3 expand_big.py in.g6 m thresh out.g6 [--res R --mod M] [--block B]
"""

import argparse
import time

import numpy as np

from expand import g6_decode, g6_encode_many, new_edges, pair_order


def build_shared(n):
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
    base = (Q.sum(axis=1) + (X[:, :m] * X[:, [m]]).sum(axis=1)).astype(np.float32)
    return X, Q, base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("m", type=int)
    ap.add_argument("thresh", type=float)
    ap.add_argument("outfile")
    ap.add_argument("--res", type=int, default=0)
    ap.add_argument("--mod", type=int, default=1)
    ap.add_argument("--block", type=int, default=1024)
    ap.add_argument("--count-only", action="store_true")
    a = ap.parse_args()

    m = a.m
    n = m + 1
    X, Q, base = build_shared(n)
    half = 1 << (n - 1)
    XT = np.ascontiguousarray(X.T)                     # (n, 2^n)
    lines = [l.strip() for l in open(a.infile) if l.strip()][a.res::a.mod]

    fh = open(a.outfile, "w")
    kept = 0
    t0 = time.time()
    for gi, line in enumerate(lines):
        e = g6_decode(line, m)
        EB = base - 2.0 * (Q @ e.astype(np.float32))   # (2^n,)
        good = []
        for s in range(0, half, a.block):
            t = min(s + a.block, half)
            Dblk = X[s:t] @ XT                          # (B, 2^n)
            Mv = np.abs(Dblk + EB[None, :]).max(axis=1)
            good.append(np.flatnonzero(Mv <= a.thresh) + s)
        good = np.concatenate(good) if good else np.array([], dtype=np.int64)
        if len(good):
            vb = np.empty((len(good), n), dtype=np.uint8)
            for i in range(n):
                vb[:, i] = ((good >> i) & 1).astype(np.uint8)
            wb = vb ^ vb[:, [n - 1]]
            E = new_edges(e, wb, m)
            if not a.count_only:
                fh.write("\n".join(g6_encode_many(E, n)) + "\n")
            kept += len(good)
        if gi and gi % 100 == 0:
            print(f"  {gi}/{len(lines)} raw_kept={kept} "
                  f"{(time.time()-t0)/gi*1000:.1f} ms/graph", flush=True)
    fh.close()
    dt = time.time() - t0
    print(f"source m={m} (order {n}) graphs={len(lines)} thresh={int(a.thresh)} "
          f"raw_extensions_kept={kept} {dt:.1f}s "
          f"({dt/max(len(lines),1)*1000:.1f} ms/graph)", flush=True)


if __name__ == "__main__":
    main()

"""Two-stage expand: reject most sign rows against the K most extreme energies
before touching the full (2^n x 2^n) block.

Stage 1 uses the K sign vectors x with the largest |E_B(x)|.  A row v survives
only if max over those K of |E_B(x) + <v,x>| <= T.  Stage 2 does the exact full
maximum on the survivors only.  Stage 1 is a relaxation of stage 2, so nothing
valid is lost.

Usage: python3 expand_prune.py in.g6 m thresh out.g6 [--res R --mod M] [--K 256]
"""

import argparse
import sys
import time

import numpy as np

from expand import g6_decode, g6_encode_many, new_edges, pair_order
from expand_big import build_shared


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("m", type=int)
    ap.add_argument("thresh", type=float)
    ap.add_argument("outfile")
    ap.add_argument("--res", type=int, default=0)
    ap.add_argument("--mod", type=int, default=1)
    ap.add_argument("--K", type=int, default=256)
    ap.add_argument("--block", type=int, default=2048)
    a = ap.parse_args()

    m = a.m
    n = m + 1
    X, Q, base = build_shared(n)
    half = 1 << (n - 1)
    Xh = np.ascontiguousarray(X[:half])          # rows v with v_n = +1
    XT = np.ascontiguousarray(X.T)
    def source_lines():
        with open(a.infile) as fp:
            for i, l in enumerate(fp):
                if i % a.mod == a.res:
                    l = l.strip()
                    if l:
                        yield l
    lines = source_lines()

    fh = open(a.outfile, "w")
    kept = 0
    stage2 = 0
    t0 = time.time()
    for gi, line in enumerate(lines):
        e = g6_decode(line, m)
        EB = base - 2.0 * (Q @ e.astype(np.float32))
        top = np.argpartition(-np.abs(EB), a.K)[:a.K]
        Dk = Xh @ np.ascontiguousarray(X[top].T)          # (half, K)
        Mk = np.abs(Dk + EB[top][None, :]).max(axis=1)
        cand = np.flatnonzero(Mk <= a.thresh)
        stage2 += len(cand)
        good = []
        if len(cand):
            for s in range(0, len(cand), a.block):
                idx = cand[s:s + a.block]
                Mv = np.abs(Xh[idx] @ XT + EB[None, :]).max(axis=1)
                good.append(idx[Mv <= a.thresh])
            good = np.concatenate(good)
        if len(good):
            vb = np.empty((len(good), n), dtype=np.uint8)
            for i in range(n):
                vb[:, i] = ((good >> i) & 1).astype(np.uint8)
            wb = vb ^ vb[:, [n - 1]]
            E = new_edges(e, wb, m)
            fh.write("\n".join(g6_encode_many(E, n)) + "\n")
            kept += len(good)
        if gi and gi % 500 == 0:
            print(f"  {gi} raw_kept={kept} stage2={stage2} "
                  f"{(time.time()-t0)/gi*1000:.2f} ms/graph",
                  file=sys.stderr, flush=True)
    fh.close()
    dt = time.time() - t0
    print(f"source m={m} (order {n}) graphs={gi+1} thresh={int(a.thresh)} "
          f"raw_extensions_kept={kept} stage2_candidates={stage2} {dt:.1f}s "
          f"({dt/max(gi+1,1)*1000:.2f} ms/graph)",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()

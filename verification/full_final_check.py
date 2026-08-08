"""Independent full check of the last tower level: no order-15 matrix has M <= 25.

Reads every order-14 record in T14_thresh25.g6 and, for each, asks whether ANY
one-vertex extension reaches M <= 25.  Written from scratch; shares no code with
the pipeline that produced the file.

Prune: the constraint |E_A(x) + <v,x>| <= 25 must hold for every x, so imposing
it on the K most extreme |E_A(x)| is a necessary condition.  That can only keep
too many candidates, never too few, so survivors get the full 2^14 check.

v and -v give the same maximum, since E_A(-x) = E_A(x), so only half the sign
rows are enumerated.
"""

import argparse

import numpy as np

N = 14                      # order of the source matrices
M_VERT = N - 1              # graphs are on 13 vertices
NPAIRS = M_VERT * (M_VERT - 1) // 2
NB = (NPAIRS + 5) // 6
THRESH = 25
K = 192


def build():
    rows = 1 << N
    idx = np.arange(rows, dtype=np.int64)
    X = np.empty((rows, N), dtype=np.float32)
    for i in range(N):
        X[:, i] = np.where((idx >> i) & 1, -1.0, 1.0)
    pairs = [(i, j) for j in range(1, M_VERT) for i in range(j)]
    Q = np.empty((rows, NPAIRS), dtype=np.float32)
    for k, (i, j) in enumerate(pairs):
        Q[:, k] = X[:, i] * X[:, j]
    base = Q.sum(axis=1) + (X[:, :M_VERT] * X[:, [M_VERT]]).sum(axis=1)
    return X, Q, base.astype(np.float32)


def decode(line):
    b = line.encode()
    body = np.frombuffer(b[1:1 + NB], dtype=np.uint8) - 63
    bits = np.unpackbits(body.astype(np.uint8)).reshape(NB, 8)[:, 2:8]
    return bits.reshape(-1)[:NPAIRS].astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("--res", type=int, default=0)
    ap.add_argument("--mod", type=int, default=1)
    a = ap.parse_args()

    X, Q, base = build()
    half = 1 << (N - 1)
    Vh = np.ascontiguousarray(X[:half])          # sign rows with v_N = +1
    VhT = np.ascontiguousarray(Vh.T)

    seen = 0
    survivors = 0
    hits = []
    best_overall = None
    with open(a.infile) as fh:
        for i, line in enumerate(fh):
            if i % a.mod != a.res:
                continue
            line = line.strip()
            if not line:
                continue
            seen += 1
            E = base - 2.0 * (Q @ decode(line))
            top = np.argpartition(-np.abs(E), K)[:K]
            Dk = np.ascontiguousarray(X[top]) @ VhT          # (K, half)
            ok = (np.abs(E[top][:, None] + Dk) <= THRESH).all(axis=0)
            cand = np.flatnonzero(ok)
            if cand.size == 0:
                continue
            survivors += int(cand.size)
            D = X @ np.ascontiguousarray(Vh[cand].T)         # full check
            Mv = np.abs(E[:, None] + D).max(axis=0)
            mn = float(Mv.min())
            if best_overall is None or mn < best_overall:
                best_overall = mn
            if mn <= THRESH:
                hits.append((line, int(mn)))
    print(f"res={a.res} seen={seen} prune_survivors={survivors} "
          f"best_full_M={best_overall} hits<=25={len(hits)}", flush=True)
    for h in hits[:5]:
        print(f"   HIT {h[0]} M={h[1]}", flush=True)


if __name__ == "__main__":
    main()

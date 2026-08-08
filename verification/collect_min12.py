"""Exhaustive F(12), collecting every minimiser, for the one-vertex lift to F(13).

Same reduction as verify_F_exhaustive: F(n) = min over non-isomorphic graphs on
n-1 vertices of M(A(G)), A(G) having last row all +1 and a_ij = -1 on edges.

Records the minimum AND every graph attaining M <= THRESH, because any order-13
matrix with M <= 18 restricts to an order-12 matrix with M <= 18, i.e. to an
optimal one.  Extending only those settles F(13).
"""

import subprocess
import sys

import numpy as np

THRESH = 18
CAP = 20_000_000


def pair_order(m):
    return [(i, j) for j in range(1, m) for i in range(j)]


def build_tables(n):
    m = n - 1
    rows = 1 << m
    idx = np.arange(rows, dtype=np.int64)
    X = np.empty((rows, m), dtype=np.float32)
    for i in range(m):
        X[:, i] = np.where((idx >> i) & 1, -1.0, 1.0)
    pairs = pair_order(m)
    Q = np.empty((rows, len(pairs)), dtype=np.float32)
    for k, (i, j) in enumerate(pairs):
        Q[:, k] = X[:, i] * X[:, j]
    return Q, (Q.sum(axis=1) + X.sum(axis=1)).astype(np.float32), len(pairs)


def main(n, out_path):
    m = n - 1
    Q, base, npairs = build_tables(n)
    QT = np.ascontiguousarray(Q.T)
    nbytes = (npairs + 5) // 6
    linelen = 1 + nbytes + 1
    chunk = max(2000, 30_000_000 // (1 << m))

    proc = subprocess.Popen(["geng", "-q", str(m)], stdout=subprocess.PIPE)
    best = None
    kept = []
    total = 0
    buf = b""
    fh = open(out_path, "w")
    while True:
        data = proc.stdout.read(linelen * chunk)
        if not data:
            break
        buf += data
        usable = (len(buf) // linelen) * linelen
        block, buf = buf[:usable], buf[usable:]
        if not block:
            continue
        arr = np.frombuffer(block, dtype=np.uint8).reshape(-1, linelen)
        body = (arr[:, 1:1 + nbytes] - 63).astype(np.uint8)
        bits = np.unpackbits(body, axis=1).reshape(-1, nbytes, 8)[:, :, 2:8]
        e = bits.reshape(len(arr), nbytes * 6)[:, :npairs].astype(np.float32)
        M = np.abs(base[None, :] - 2.0 * (e @ QT)).max(axis=1)
        mn = float(M.min())
        if best is None or mn < best:
            best = mn
        hit = np.flatnonzero(M <= THRESH)
        if len(hit) and len(kept) < CAP:
            for k in hit:
                fh.write(bytes(arr[k, :1 + nbytes]).decode() + "\n")
            kept.extend(hit[:1])
            if len(kept) == 1:
                fh.flush()
        total += len(arr)
        if total % 50_000_000 < chunk:
            print(f"  ...{total} graphs, running min {int(best)}", flush=True)
    proc.stdout.close()
    rc = proc.wait()
    fh.close()
    assert rc == 0, f"geng exited {rc}"
    assert not buf, "trailing bytes"
    n_kept = sum(1 for _ in open(out_path))
    print(f"n={n} graphs_on_{m}_vertices={total} F({n})={int(best)} "
          f"minimisers_with_M<={THRESH}: {n_kept}", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2])

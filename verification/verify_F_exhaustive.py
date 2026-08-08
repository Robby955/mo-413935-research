"""Independent exhaustive verification of F(n).

Every sign matrix is switching-equivalent to one whose last row is all +1, and
M is invariant under switching and under relabelling.  So

    F(n) = min over non-isomorphic graphs G on n-1 vertices of M(A(G)),

where A(G) has a_{in} = +1 and, for i<j<n, a_ij = -1 exactly on the edges of G.
Graphs come from nauty geng, so the enumeration is isomorph-free and complete.

With S = sum_{i<j<n} x_i x_j and T = sum_{j<n} x_j x_n,

    E_G(x) = (S + T) - 2 * sum_{(i,j) in E(G)} x_i x_j,

so a chunk of graphs is one matrix product against the pair-product table.
"""

import subprocess
import sys

import numpy as np


def pair_order(m):
    """graph6 bit order: (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),..."""
    return [(i, j) for j in range(1, m) for i in range(j)]


def build_tables(n):
    m = n - 1                      # vertices carried by the graph
    rows = 1 << m                  # sign vectors with x_n = +1
    idx = np.arange(rows, dtype=np.int64)
    X = np.empty((rows, m), dtype=np.float32)
    for i in range(m):
        X[:, i] = np.where((idx >> i) & 1, -1.0, 1.0)
    pairs = pair_order(m)
    Q = np.empty((rows, len(pairs)), dtype=np.float32)
    for k, (i, j) in enumerate(pairs):
        Q[:, k] = X[:, i] * X[:, j]
    S = Q.sum(axis=1)              # sum over all i<j<n
    T = X.sum(axis=1)              # x_n = +1
    return Q, (S + T).astype(np.float32), len(pairs)


def run(n, chunk=None):
    m = n - 1
    Q, base, npairs = build_tables(n)
    if chunk is None:                      # cap the energy block near 120 MB
        chunk = max(2000, 30_000_000 // (1 << m))
    nbytes = (npairs + 5) // 6
    linelen = 1 + nbytes + 1       # order byte + payload + newline
    QT = np.ascontiguousarray(Q.T)

    proc = subprocess.Popen(["geng", "-q", str(m)], stdout=subprocess.PIPE)
    best = None
    best_g6 = None
    total = 0
    buf = b""
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
        assert arr[0, 0] == m + 63, f"unexpected graph6 order byte {arr[0,0]}"
        body = (arr[:, 1:1 + nbytes] - 63).astype(np.uint8)
        bits = np.unpackbits(body, axis=1).reshape(-1, nbytes, 8)[:, :, 2:8]
        e = bits.reshape(len(arr), nbytes * 6)[:, :npairs].astype(np.float32)
        E = base[None, :] - 2.0 * (e @ QT)
        M = np.abs(E).max(axis=1)
        k = int(np.argmin(M))
        if best is None or M[k] < best:
            best = float(M[k])
            best_g6 = bytes(arr[k, :1 + nbytes]).decode()
        total += len(arr)
    proc.stdout.close()
    rc = proc.wait()
    assert rc == 0, f"geng exited {rc}"
    assert not buf, "trailing bytes in geng stream"
    print(f"n={n}  graphs_on_{m}_vertices={total}  F({n})={int(best)}  "
          f"witness_graph6={best_g6}", flush=True)
    return int(best)


if __name__ == "__main__":
    for n in [int(a) for a in sys.argv[1:]]:
        run(n)

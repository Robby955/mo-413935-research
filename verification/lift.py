"""One-vertex lift: rule out low-M matrices at order n+1 from the minimisers at order n.

If A has order n+1 with M(A) <= T, deleting the last vertex leaves B of order n
with M(B) <= T.  When T equals F(n), that forces B to be a minimiser.  So it is
enough to extend every minimiser in all 2^n ways and check none reaches T.

Energy of the extension: with the new vertex's signs v and x_{n+1} = +1,
    E(x) = E_B(x) + <v, x>,
so a whole minimiser is one matrix product over (2^n sign vectors) x (2^n v).
"""

import sys

import numpy as np


def pair_order(m):
    return [(i, j) for j in range(1, m) for i in range(j)]


def energies_order_n(g6_lines, n):
    """E_B over ALL 2^n sign vectors, for each graph on n-1 vertices."""
    m = n - 1
    pairs = pair_order(m)
    npairs = len(pairs)
    nb = (npairs + 5) // 6
    rows = 1 << n                      # all sign vectors, x_i free for i<=n
    idx = np.arange(rows, dtype=np.int64)
    X = np.empty((rows, n), dtype=np.float32)
    for i in range(n):
        X[:, i] = np.where((idx >> i) & 1, -1.0, 1.0)
    # pair products among the first m coordinates, plus the last vertex's edges
    Q = np.empty((rows, npairs), dtype=np.float32)
    for k, (i, j) in enumerate(pairs):
        Q[:, k] = X[:, i] * X[:, j]
    base = Q.sum(axis=1) + (X[:, :m] * X[:, [m]]).sum(axis=1)
    out = []
    for line in g6_lines:
        g6 = line.strip().encode()
        body = np.frombuffer(g6[1:1 + nb], dtype=np.uint8) - 63
        bits = np.unpackbits(body.astype(np.uint8)).reshape(nb, 8)[:, 2:6 + 2]
        e = bits.reshape(-1)[:npairs].astype(np.float32)
        out.append(base - 2.0 * (Q @ e))
    return np.array(out), X


def lift(g6_lines, n, thresh):
    """Return the minimum M over every one-vertex extension of the given order-n
    matrices, and a witness if any extension reaches <= thresh."""
    E, X = energies_order_n(g6_lines, n)      # (G, 2^n)
    D = X @ X.T                               # (2^n, 2^n): <v, x>, shared by all
    best = None
    witness = None
    for gi in range(E.shape[0]):
        allM = np.abs(E[gi][:, None] + D).max(axis=0)   # per v
        k = int(np.argmin(allM))
        if best is None or allM[k] < best:
            best = float(allM[k])
            witness = (gi, k)
        if best <= thresh:
            break
    return best, witness


if __name__ == "__main__":
    path, n, thresh = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    lines = [l for l in open(path) if l.strip()]
    print(f"minimisers read: {len(lines)}  order n={n}  threshold {thresh}",
          flush=True)
    best, w = lift(lines, n, thresh)
    print(f"min M over all one-vertex extensions = {int(best)}", flush=True)
    print("REACHES threshold" if best <= thresh else
          f"NO extension reaches {int(thresh)}", flush=True)

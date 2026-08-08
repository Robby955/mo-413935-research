"""Utilities: M of a graph6 record, and the graph6 normal form of a sign matrix."""

import subprocess
import sys

import numpy as np


def pair_order(m):
    return [(i, j) for j in range(1, m) for i in range(j)]


def g6_bits(line, m):
    npairs = m * (m - 1) // 2
    nb = (npairs + 5) // 6
    b = line.strip().encode()
    assert b[0] == m + 63
    body = np.frombuffer(b[1:1 + nb], dtype=np.uint8) - 63
    bits = np.unpackbits(body.astype(np.uint8)).reshape(nb, 8)[:, 2:8]
    return bits.reshape(-1)[:npairs].astype(np.int8)


def matrix_of(line, m):
    """Order n = m+1 sign matrix with last row/col all +1."""
    n = m + 1
    A = np.ones((n, n), dtype=np.int8)
    e = g6_bits(line, m)
    for k, (i, j) in enumerate(pair_order(m)):
        A[i, j] = A[j, i] = 1 - 2 * e[k]
    np.fill_diagonal(A, 0)
    return A


def M_of(A):
    n = A.shape[0]
    rows = 1 << (n - 1)
    idx = np.arange(rows, dtype=np.int64)
    X = np.ones((rows, n), dtype=np.int32)
    for i in range(n - 1):
        X[:, i] = np.where((idx >> i) & 1, -1, 1)
    E = np.zeros(rows, dtype=np.int32)
    for i in range(n):
        for j in range(i + 1, n):
            E += int(A[i, j]) * X[:, i] * X[:, j]
    return int(np.abs(E).max())


def g6_of_matrix(A):
    """Switch the last row of A to all +1, return graph6 of the resulting graph."""
    n = A.shape[0]
    d = np.array([A[i, n - 1] for i in range(n - 1)] + [1], dtype=np.int8)
    B = (d[:, None] * A * d[None, :])
    assert all(B[i, n - 1] == 1 for i in range(n - 1))
    m = n - 1
    npairs = m * (m - 1) // 2
    bits = np.zeros(npairs, dtype=np.uint8)
    for k, (i, j) in enumerate(pair_order(m)):
        bits[k] = 1 if B[i, j] == -1 else 0
    nb = (npairs + 5) // 6
    pad = np.zeros(nb * 6 - npairs, dtype=np.uint8)
    bb = np.concatenate([bits, pad]).reshape(nb, 6)
    two = np.zeros((nb, 2), dtype=np.uint8)
    packed = np.packbits(np.concatenate([two, bb], axis=1), axis=1)[:, 0] + 63
    return chr(m + 63) + "".join(chr(c) for c in packed)


def canon(line):
    out = subprocess.run(["labelg", "-q"], input=line + "\n", capture_output=True,
                         text=True)
    return out.stdout.strip()


if __name__ == "__main__":
    path, m = sys.argv[1], int(sys.argv[2])
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        A = matrix_of(line, m)
        print(f"{line}  order={m+1}  M={M_of(A)}  canon={canon(line)}")

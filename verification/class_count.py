"""Count switching classes, not rooted records, in a list of graph6 records.

The enumerations in this directory emit one rooted normal form per graph: a
representative of a switching class whose last row and column are +1.  A single
switching class contributes one such record for each orbit of its automorphism
group on the root, so raw record counts overstate the number of classes.

For an order-n signing A and each vertex k, switching by the k-th row makes
every edge at k positive; deleting k then leaves an ordinary graph on n-1
vertices.  Canonicalising those n graphs with nauty's labelg gives a set that
depends only on the switching class, so two records lie in the same class
exactly when their sets agree.

Usage:
    python3 verification/class_count.py FILE.g6 [FILE.g6 ...]

Requires nauty (labelg) on PATH.
"""

import subprocess
import sys

import numpy as np

from minmax_search import verify

LIMIT = 50          # records per file; canonicalising is O(n) labelg calls each


def decode(line):
    """graph6 record on m vertices -> order-(m+1) signing with last row all +1."""
    b = line.strip().encode()
    m = b[0] - 63
    npairs = m * (m - 1) // 2
    nb = (npairs + 5) // 6
    body = np.frombuffer(b[1:1 + nb], dtype=np.uint8) - 63
    bits = np.unpackbits(body.astype(np.uint8)).reshape(nb, 8)[:, 2:8]
    bits = bits.reshape(-1)[:npairs]
    A = np.ones((m + 1, m + 1), dtype=np.int8)
    np.fill_diagonal(A, 0)
    k = 0
    for j in range(1, m):
        for i in range(j):
            if bits[k]:
                A[i, j] = A[j, i] = -1
            k += 1
    return A, m + 1


def g6_encode(adj, m):
    bits = [int(adj[i, j]) for j in range(1, m) for i in range(j)]
    out = chr(m + 63)
    for c in range(0, len(bits), 6):
        v = 0
        for t in range(6):
            v = (v << 1) | (bits[c + t] if c + t < len(bits) else 0)
        out += chr(v + 63)
    return out


def rooting_set(A, n):
    """Canonical forms of all n rootings; an invariant of the switching class."""
    forms = []
    for k in range(n):
        d = A[k].copy()
        d[k] = 1
        B = A * np.outer(d, d)
        keep = [i for i in range(n) if i != k]
        S = B[np.ix_(keep, keep)]
        forms.append(g6_encode((S == -1).astype(np.int8), n - 1))
    done = subprocess.run(["labelg", "-q"], input="\n".join(forms) + "\n",
                          capture_output=True, text=True)
    if done.returncode != 0:
        sys.exit("labelg failed; is nauty on PATH?")
    return frozenset(l.strip() for l in done.stdout.split() if l.strip())


def main(paths):
    for path in paths:
        lines = [l for l in open(path) if l.strip()]
        if len(lines) > LIMIT:
            print(f"{path}: {len(lines)} records, above the {LIMIT} record "
                  f"limit for class counting")
            continue
        classes = {}
        for line in lines:
            A, n = decode(line)
            classes.setdefault(rooting_set(A, n), []).append(
                (line.strip(), verify(A, n)))
        print(f"{path}: {len(lines)} rooted records, {len(classes)} switching "
              f"classes")
        for i, recs in enumerate(classes.values(), 1):
            codes = ", ".join(r for r, _ in recs)
            print(f"   class {i}: M = {recs[0][1]}, records: {codes}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])

"""Sample geng residue classes and count graphs with M(A(G)) <= THRESH.

Usage:  python3 sample_T.py m thresh mod res1 res2 ...

m       vertices of the graph (matrix order is m+1)
thresh  keep graphs with M <= thresh
mod     geng res/mod split
res*    residues to run

Prints one line per residue: residue, graphs seen, hits, wall seconds.
Then a pooled estimate scaled to the full count for that m, when known.
"""

import os
import subprocess
import sys
import time

import numpy as np

TOTAL = {8: 12346, 9: 274668, 10: 12005168, 11: 1018997864}


def pair_order(m):
    return [(i, j) for j in range(1, m) for i in range(j)]


def build_tables(m):
    rows = 1 << m
    idx = np.arange(rows, dtype=np.int64)
    X = np.empty((rows, m), dtype=np.float32)
    for i in range(m):
        X[:, i] = np.where((idx >> i) & 1, -1.0, 1.0)
    pairs = pair_order(m)
    Q = np.empty((rows, len(pairs)), dtype=np.float32)
    for k, (i, j) in enumerate(pairs):
        Q[:, k] = X[:, i] * X[:, j]
    base = (Q.sum(axis=1) + X.sum(axis=1)).astype(np.float32)
    return np.ascontiguousarray(Q.T), base, len(pairs)


def scan(m, thresh, args, chunk=100_000, keep_path=None):
    QT, base, npairs = build_tables(m)
    nbytes = (npairs + 5) // 6
    linelen = 1 + nbytes + 1
    cmd = ["geng", "-q", str(m)] + args
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    total = 0
    hits = 0
    best = None
    buf = b""
    fh = open(keep_path, "w") if keep_path else None
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
        hit = np.flatnonzero(M <= thresh)
        hits += len(hit)
        if fh is not None:
            for k in hit:
                fh.write(bytes(arr[k, :1 + nbytes]).decode() + "\n")
        total += len(arr)
    proc.stdout.close()
    rc = proc.wait()
    if fh is not None:
        fh.close()
    assert rc == 0, f"geng exited {rc}"
    assert not buf, "trailing bytes"
    return total, hits, best


def main():
    m = int(sys.argv[1])
    thresh = float(sys.argv[2])
    mod = int(sys.argv[3])
    residues = [int(a) for a in sys.argv[4:]]
    seen = 0
    hits = 0
    best = None
    t_all = time.time()
    for r in residues:
        t0 = time.time()
        args = [] if mod == 1 else [f"{r}/{mod}"]
        keep = os.environ.get("KEEP_PREFIX")
        kp = f"{keep}_m{m}_t{int(thresh)}_{r}of{mod}.g6" if keep else None
        tot, h, b = scan(m, thresh, args, keep_path=kp)
        dt = time.time() - t0
        seen += tot
        hits += h
        if best is None or (b is not None and b < best):
            best = b
        print(f"res {r}/{mod}: graphs={tot} hits={h} min={int(b)} "
              f"{dt:.1f}s  ({tot/max(dt,1e-9):,.0f} g/s)", flush=True)
    dt = time.time() - t_all
    print(f"POOLED m={m} thresh={int(thresh)} classes={len(residues)}/{mod} "
          f"graphs_seen={seen} hits={hits} min={int(best)} {dt:.1f}s", flush=True)
    if m in TOTAL:
        frac = seen / TOTAL[m]
        print(f"sampling fraction = {frac:.6f} of {TOTAL[m]}", flush=True)
        if hits:
            print(f"extrapolated |T| = {hits / frac:,.0f}", flush=True)


if __name__ == "__main__":
    main()

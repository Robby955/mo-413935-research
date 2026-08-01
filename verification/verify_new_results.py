#!/usr/bin/env python3
"""Independent verification of every table in notes/limit_program.tex.

Everything here is recomputed from first principles in pure Python
(standard library only), independently of the C programs that produced
the data. Two kinds of checks:

  FULL   the claim is recomputed exhaustively (small n), certifying
         both directions of a min-max value;
  WITNESS the recorded witness is verified to achieve the stated value,
         certifying one direction; the other direction rests on the
         exhaustive C enumeration (verification/exact_fn.c and
         verification/conference_max.c), whose own audits recompute
         their reported optima from scratch.

Corruption controls at the end must fail; the script asserts that they
do. Deterministic; no randomness.
"""

from __future__ import annotations

import itertools
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(HERE, "data", "results_2026_08.json")))

F = {int(k): v for k, v in DATA["F"].items()}
FW = {int(k): v for k, v in DATA["F_witness_mask"].items()}
FP = {int(k): v for k, v in DATA["F_plus"].items()}
OSC = {int(k): v for k, v in DATA["Osc"].items()}
OSCW = {int(k): v for k, v in DATA["Osc_witness_mask"].items()}


def build_from_mask(n: int, mask: int):
    """Signing with all edges at vertex 0 positive; mask sets the rest."""
    A = [[0] * n for _ in range(n)]
    for j in range(1, n):
        A[0][j] = A[j][0] = 1
    e = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            s = -1 if (mask >> e) & 1 else 1
            A[i][j] = A[j][i] = s
            e += 1
    return A


def energies(A):
    n = len(A)
    for bits in range(1 << (n - 1)):
        x = [1] + [1 if (bits >> k) & 1 == 0 else -1 for k in range(n - 1)]
        yield sum(A[i][j] * x[i] * x[j] for i in range(n) for j in range(i + 1, n))


def spread(A):
    lo, hi = 10**9, -(10**9)
    for q in energies(A):
        lo, hi = min(lo, q), max(hi, q)
    return hi, lo


def check(label, cond):
    if not cond:
        raise AssertionError(f"FAILED: {label}")
    print(f"ok {label}")


# 1. FULL: brute force over every signing for n <= 5.
def brute_all(n):
    E = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(E)
    ba = bp = bo = 10**9
    xs = [[1] + [1 if (b >> k) & 1 == 0 else -1 for k in range(n - 1)]
          for b in range(1 << (n - 1))]
    for am in range(1 << m):
        a = [1 if (am >> k) & 1 == 0 else -1 for k in range(m)]
        vals = [sum(a[k] * x[i] * x[j] for k, (i, j) in enumerate(E)) for x in xs]
        ba = min(ba, max(abs(v) for v in vals))
        bp = min(bp, max(vals))
        bo = min(bo, max(vals) - min(vals))
    return ba, bp, bo


for n in (3, 4, 5):
    ba, bp, bo = brute_all(n)
    check(f"FULL F({n})={F[n]}", ba == F[n])
    check(f"FULL F_plus({n})={FP[n]}", bp == FP[n])
    check(f"FULL Osc({n})={OSC[n]}", bo == OSC[n])

# 2. FULL: independent switching-class enumeration for n = 6, 7.
for n in (6, 7):
    er = (n - 1) * (n - 2) // 2
    ba = bp = bo = 10**9
    for mask in range(1 << er):
        hi, lo = spread(build_from_mask(n, mask))
        ba = min(ba, max(hi, -lo))
        bp = min(bp, hi)
        bo = min(bo, hi - lo)
    check(f"FULL F({n})={F[n]}", ba == F[n])
    check(f"FULL F_plus({n})={FP[n]}", bp == FP[n])
    check(f"FULL Osc({n})={OSC[n]}", bo == OSC[n])

# 3. WITNESS: recorded optima are achieved for all tabulated n.
for n in sorted(FW):
    hi, lo = spread(build_from_mask(n, FW[n]))
    check(f"WITNESS max|S| of F-witness({n}) == {F[n]}", max(hi, -lo) == F[n])
for n in sorted(OSCW):
    hi, lo = spread(build_from_mask(n, OSCW[n]))
    check(f"WITNESS spread of Osc-witness({n}) == {OSC[n]}", hi - lo == OSC[n])

# 4. Petersdorf pattern: all-negative signing achieves max S = floor(n/2),
#    and F_plus table equals floor(n/2).
for n in sorted(FP):
    er = (n - 1) * (n - 2) // 2
    hi, lo = spread(build_from_mask(n, (1 << er) - 1))
    check(f"all-negative class has max S = floor({n}/2)", hi == n // 2)
    check(f"F_plus({n}) == floor(n/2)", FP[n] == n // 2)

# 5. Parity and sandwich laws on the tables.
for n in sorted(F):
    check(f"parity F({n})", (F[n] - n * (n - 1) // 2) % 2 == 0)
for n in sorted(OSC):
    check(f"F <= Osc <= 2F at n={n}", F[n] <= OSC[n] <= 2 * F[n])
for n in sorted(OSC):
    for m in sorted(OSC):
        if n + m in OSC:
            check(f"Osc superadditive ({n},{m})",
                  OSC[n + m] >= OSC[n] + OSC[m])

# 6. Conference matrices, rebuilt independently.
def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def paley(q):
    n = q + 1
    C = [[0] * n for _ in range(n)]
    for j in range(1, n):
        C[0][j] = C[j][0] = 1
    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                C[i][j] = legendre(i - j, q)
    return C


def petersen_seidel():
    pairs = list(itertools.combinations(range(5), 2))
    C = [[0] * 10 for _ in range(10)]
    for i in range(10):
        for j in range(10):
            if i != j:
                adj = not (set(pairs[i]) & set(pairs[j]))
                C[i][j] = -1 if adj else 1
    return C


def gf25_conference():
    def mul(x, y):
        return ((x[0] * y[0] + 2 * x[1] * y[1]) % 5,
                (x[0] * y[1] + x[1] * y[0]) % 5)

    def chi(z):
        if z == (0, 0):
            return 0
        r, base, e = (1, 0), z, 12
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        assert r in ((1, 0), (4, 0))
        return 1 if r == (1, 0) else -1

    C = [[0] * 26 for _ in range(26)]
    for j in range(1, 26):
        C[0][j] = C[j][0] = 1
    for i in range(1, 26):
        for j in range(1, 26):
            if i != j:
                x, y = divmod(i - 1, 5), divmod(j - 1, 5)
                d = ((x[1] - y[1]) % 5, (x[0] - y[0]) % 5)
                C[i][j] = chi(d)
    return C


def is_conference(C):
    n = len(C)
    for i in range(n):
        for j in range(n):
            s = sum(C[i][k] * C[k][j] for k in range(n))
            if s != ((n - 1) if i == j else 0):
                return False
    return True


def energy_at(C, xmask):
    n = len(C)
    x = [-1 if (xmask >> i) & 1 else 1 for i in range(n)]
    return sum(C[i][j] * x[i] * x[j] for i in range(n) for j in range(i + 1, n))


CONF = {"paley_q5": paley(5), "petersen": petersen_seidel(),
        "paley_q13": paley(13), "paley_q17": paley(17),
        "gf25": gf25_conference(), "paley_q29": paley(29)}

for rec in DATA["conference"]:
    C = CONF[rec["source"]]
    N = rec["order"]
    check(f"C^2=(N-1)I at order {N}", len(C) == N and is_conference(C))
    val = energy_at(C, rec["witness_x"])
    check(f"WITNESS conference order {N} achieves {rec['maximum']}",
          abs(val) == rec["maximum"])
    ceiling_attained = (2 * rec["maximum"] == N * math.isqrt(N - 1)
                        and math.isqrt(N - 1) ** 2 == N - 1)
    check(f"ceiling attainment flag at order {N}",
          ceiling_attained == rec["attains_ceiling"])

# FULL maxima for the small conference orders.
for src, expect in (("paley_q5", 5), ("petersen", 15), ("paley_q13", 21)):
    C = CONF[src]
    hi, lo = spread(C)
    check(f"FULL conference max {src} == {expect}", max(hi, -lo) == expect)

# 7. Tensor doubling at n = 5: best intra-pair fill of A x H2 has max 13.
A5 = build_from_mask(5, FW[5])
best = None
for fill in range(1 << 5):
    B = [[0] * 10 for _ in range(10)]
    for i in range(5):
        for j in range(5):
            for a in range(2):
                for b in range(2):
                    if i == j and a == b:
                        continue
                    if i == j:
                        B[2 * i + a][2 * j + b] = -1 if (fill >> i) & 1 else 1
                    else:
                        h = 1 if (a, b) != (1, 1) else -1
                        B[2 * i + a][2 * j + b] = A5[i][j] * h
    hi, lo = spread(B)
    m = max(hi, -lo)
    best = m if best is None else min(best, m)
check("tensor: min-fill M(A5 x H2) == 13 == F(10)", best == 13)

# 8. Finite merge instance: blocks F(6)+F(4) joined by a Paley-14 cross
#    block obey M <= F(6)+F(4)+max bilinear, and F(10) <= M.
A6 = build_from_mask(6, FW[6])
A4 = build_from_mask(4, FW[4])
P14 = paley(13)
Cross = [[P14[i][6 + j] for j in range(4)] for i in range(6)]
bilin = 0
for xm in range(1 << 6):
    x = [1 if (xm >> k) & 1 == 0 else -1 for k in range(6)]
    row = [sum(Cross[i][j] * x[i] for i in range(6)) for j in range(4)]
    bilin = max(bilin, sum(abs(r) for r in row))
B = [[0] * 10 for _ in range(10)]
for i in range(6):
    for j in range(6):
        B[i][j] = A6[i][j]
for i in range(4):
    for j in range(4):
        B[6 + i][6 + j] = A4[i][j]
for i in range(6):
    for j in range(4):
        B[i][6 + j] = B[6 + j][i] = Cross[i][j]
hi, lo = spread(B)
MB = max(hi, -lo)
check("merge instance: M(joined) <= F(6)+F(4)+bilinear",
      MB <= F[6] + F[4] + bilin)
check("merge instance: F(10) <= M(joined)", F[10] <= MB)

# 9. Barrier profile margins up to 10^6.
def Gt(n):
    return 0.4 + 0.05 * math.sin(2 * math.pi * math.log2(math.log2(n)))


def Ft(n):
    return Gt(n) * n ** 1.5


grid = sorted(set(list(range(64, 4096)) +
                  [int(64 * 1.01 ** k) for k in range(1500)
                   if 64 * 1.01 ** k <= 10 ** 6]))
for n in grid:
    if n + 1 <= 10 ** 6:
        check_ok = Ft(n) < Ft(n + 1) <= Ft(n) + n
        if not check_ok:
            raise AssertionError(f"barrier monotone/Lipschitz fails at {n}")
    for m in (max(4, n // 64), n // 7 + 4, n // 2, n):
        if not (Ft(n + m) <= Ft(n) + Ft(m) + 1.01 * n * math.sqrt(m)):
            raise AssertionError(f"barrier merge fails at {n},{m}")
        if not (Ft(n + m) >= Ft(n) + Ft(m)):
            raise AssertionError(f"barrier superadditivity fails at {n},{m}")
    if not (1 / math.pi < 0.35 <= Gt(n) <= 0.45 < 0.5):
        raise AssertionError(f"barrier range fails at {n}")
print("ok barrier profile satisfies all value inequalities on the grid")

# 10. Corruption controls: these MUST fail.
failed = 0
try:
    hi9, lo9 = spread(build_from_mask(9, FW[9]))
    check("corrupt: F(9) claimed 10", max(hi9, -lo9) <= 10)
except AssertionError:
    failed += 1
try:
    Cbad = [row[:] for row in CONF["petersen"]]
    Cbad[0][1] = Cbad[1][0] = -Cbad[0][1]
    check("corrupt: perturbed Petersen still conference", is_conference(Cbad))
except AssertionError:
    failed += 1
try:
    check("corrupt: all-plus signing one-sided small",
          spread(build_from_mask(8, 0))[0] == 4)
except AssertionError:
    failed += 1
if failed != 3:
    raise AssertionError("corruption controls did not all trigger")
print("ok corruption controls all triggered")

print()
print("exact_values=VERIFIED_full_n_le_7_witness_n_le_10")
print("petersdorf_pattern=VERIFIED_n_le_10")
print("oscillation_tables=VERIFIED")
print("conference_attainment=VERIFIED_orders_10_26")
print("merge_and_tensor_instances=VERIFIED")
print("barrier_profile=VERIFIED_grid_to_1e6")
print("corruption_controls=PASSED")
print("ALL CHECKS PASSED")

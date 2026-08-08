"""Are the optimal matrices flat-spectrum, or do they trade norm for alignment?

For each n with a known/found optimum A:
  rho   = ||A||_op / sqrt(n-1)        1.0 <=> conference (flattest possible)
  m     = M(A) / (0.5 * n * ||A||)    alignment fraction; 2/pi = 0.6366 is the
                                      Nesterov rounding floor for ANY matrix
  ratio = M(A) / n^{3/2}
"""
import math, numpy as np, itertools
from minmax_search import verify, random_matrix, columns, energies, basin_hop

def stats(A, label):
    n = A.shape[0]
    M = verify(A, n)
    rho = float(np.abs(np.linalg.eigvalsh(A.astype(float))).max())
    print(f"{label:<26} n={n:3d} M={M:5d} ratio={M/n**1.5:.4f} "
          f"||A||/sqrt(n-1)={rho/math.sqrt(n-1):.4f} "
          f"m=M/(n||A||/2)={M/(0.5*n*rho):.4f}", flush=True)

def paley_prime(q):
    res = {(k*k) % q for k in range(1, q)}
    chi = np.zeros(q, dtype=np.int8)
    for k in range(1, q): chi[k] = 1 if k in res else -1
    n = q+1
    C = np.ones((n, n), dtype=np.int8)
    for i in range(1, n):
        for j in range(1, n):
            if i != j: C[i, j] = chi[(i-j) % q]
    np.fill_diagonal(C, 0)
    return C

# exact optimum witnesses found earlier in this session
W = {
10: [1,1,1,-1,1,1,1,-1,1,-1,-1,1,1,1,1,1,1,-1,-1,1,1,1,1,-1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,1,1],
14: [1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1,1,-1,1,1,-1,1,1,-1,-1,-1,-1,1,1,-1,1,-1,1,1,-1,-1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,-1,1,1,-1,1,-1,1,1,1,-1,1,1,-1,1],
}
for n, w in W.items():
    A = np.zeros((n, n), dtype=np.int8); k = 0
    for i in range(n):
        for j in range(i+1, n):
            A[i, j] = A[j, i] = w[k]; k += 1
    stats(A, f"exact optimum n={n}")

# n=13 optimum via vertex deletion from the n=14 optimum (M=20, matches repo)
A14 = np.zeros((14, 14), dtype=np.int8); k = 0
for i in range(14):
    for j in range(i+1, 14):
        A14[i, j] = A14[j, i] = W[14][k]; k += 1
best = None
for d in range(14):
    keep = [i for i in range(14) if i != d]
    B = np.ascontiguousarray(A14[np.ix_(keep, keep)])
    v = verify(B, 13)
    if best is None or v < best[0]: best = (v, B)
stats(best[1], "optimum n=13 (deletion)")

for q in (5, 13, 17):
    stats(paley_prime(q), f"Paley conference q={q}")

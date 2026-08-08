"""max_x x^T A x / n^{3/2} for A = Haar-random with FLAT spectrum +-sqrt(n).

Same spectral profile as a conference matrix (A^2 = n I), but no arithmetic
structure.  Compares against:
  2/pi   = 0.6366  Nesterov rounding floor (any such A)
  sqrt(15)/4 = 0.9682  my first-moment/annealed cap for Haar-flat A
  ~0.94-0.95         what Paley conference matrices actually achieve
"""
import math, numpy as np

def cube_max(A, restarts, rng, iters=600):
    n = A.shape[0]
    Af = A.astype(np.float32)
    best = 0.0
    for sign in (1.0, -1.0):
        X = rng.choice(np.array([-1, 1], dtype=np.float32), size=(restarts, n))
        for _ in range(iters):
            G = X @ Af
            gain = -2.0 * sign * X * G
            k = np.argmax(gain, axis=1)
            g = gain[np.arange(restarts), k]
            act = g > 1e-9
            if not act.any(): break
            X[np.arange(restarts)[act], k[act]] *= -1
        v = sign * np.einsum('ij,ij->i', X @ Af, X)
        best = max(best, float(v.max()))
    return best

print(f"{'n':>5} {'m=max/n^1.5':>12}   (2/pi=0.6366, annealed cap=0.9682)")
for n in [64, 128, 256, 512, 1024]:
    rng = np.random.default_rng(413935 + n)
    vals = []
    reps = 3 if n <= 512 else 2
    for r in range(reps):
        Z = rng.standard_normal((n, n)).astype(np.float64)
        Q, _ = np.linalg.qr(Z)
        d = np.ones(n); d[: n // 2] = -1.0
        rng.shuffle(d)
        A = (Q * d) @ Q.T * math.sqrt(n)   # spectrum exactly +-sqrt(n)
        A = (A + A.T) / 2
        np.fill_diagonal(A, 0.0)           # match the zero-diagonal constraint
        restarts = 3000 if n <= 256 else 800
        m = cube_max(A, restarts, rng) / n ** 1.5
        vals.append(m)
    print(f"{n:>5} {np.mean(vals):>12.4f}   spread={np.std(vals):.4f} reps={reps}", flush=True)

# MathOverflow 413935: proof-search ledger

Status date: 2026-08-01 (second pass; the 2026-07-25 ledger is
preserved unchanged below and remains accurate except where the
August audit section says otherwise)

## August 2026 program: summary

Full statements and proofs: `notes/limit_program.tex`. Claim-by-claim
audit of the July material: `AUDIT.md`. Independent recomputation of
every table: `verification/verify_new_results.py` (ALL CHECKS PASSED).

New proved facts.

1. Exact values (exhaustive over switching classes, audited):
   F(2..10) = 1, 3, 4, 4, 5, 9, 10, 12, 13. Not in OEIS. One-sided
   optimum floor(n/2) for n <= 10 (equivalent to Petersdorf 1966).
   F(n) has the parity of C(n,2).
2. F(10) = 13 beats the order-10 conference class (15). Optimal
   spectra are flat only at n = 5, 6. Flatness and optimality diverge.
3. Conference ceiling attainment: symmetric conference matrices attain
   max = (1/2) N sqrt(N-1) exactly at N = 10 (Petersen) and N = 26
   (GF(25) Paley); impossible when sqrt(N-1) is irrational. The July
   "conference spectral shortcut" wall was an artifact of the orders
   sampled. Values at 6, 14, 18, 30: 5, 21, 33, 75 (ratios rising:
   0.340, 0.401, 0.432, 0.456).
4. Even-odd inequality: M(B) >= |even| + |odd| for every bipartition
   and every (x, y); hence every two-block construction pays
   Theta(n sqrt m) in the cross block, both directions pinned.
5. Merge theorem: F(n+m) <= F(n) + F(m) + sqrt(nm q(n)),
   q(n) = (1+o(1)) n. Corollaries: G is continuous on multiplicative
   windows (modulus C sqrt(eps)); the limit set of G is a closed
   interval in [1/pi, 1/2]; the limit exists iff it exists along
   conference orders; a one-scale-to-all-scales amplification
   inequality would imply convergence.
6. Oscillation functional Osc(n) = min_A (max S - min S):
   superadditive (exactly), F <= Osc <= 2F, and the Gaussian-rounding
   method of the July note, optimized over all correlation-matrix
   pairs, proves exactly F >= Osc/2. Data: Osc = 2F at n = 4, 5, 6,
   8, 9, 10 and 2F - 2 at n = 3, 7. The 1/pi bound is the linear
   slice; no spectrally defined slice can beat it (Cauchy-Schwarz),
   and the adversary of every polynomial slice is flat-spectrum.
7. Fixed-temperature equivalence: |Phi_{n,beta} - G(n)| <=
   log 2/(beta sqrt n). The beta-uniformity worry is empty; the open
   problem is pure n -> infinity interpolation at fixed beta.
8. Relaxations collapse: box couplings give 0; sphere couplings give
   exactly sqrt(C(n,2)) = Theta(n) (single-edge concentration);
   near-optimal signings tend to the zero graphon (cut norm
   <= 4 M = O(n^{3/2})). First-order limit theories are blind here.
9. Consistency barrier: G~(n) = 2/5 + (1/20) sin(2 pi log2 log2 n)
   satisfies every inequality proved about the values F(n) in this
   repository (verified with margins to 10^6) and diverges. No local
   estimate, however sharpened, can decide the question.

New certified upper bounds: F(11) <= 19, F(12) <= 18 (tensor of the
exact n = 6 optimizer), F(13) <= 20, F(14) <= 21, F(15) <= 29,
F(16) <= 32, F(17) <= 32, F(20) <= 44.

New wall data.

- Tensor doubling A x H2 (exact, best intra-pair fill) moves ratios
  toward roughly 0.47-0.50 instead of preserving them:
  0.358 -> 0.411, 0.340 -> 0.433, 0.486 -> 0.477, 0.442 -> 0.500,
  0.411 -> 0.492. At n = 5 the doubled matrix is exactly optimal at
  n = 10, the loss being absorbed by parity.
- Restriction anti-decay: principal submatrices of the order-30 Paley
  matrix have ratios decreasing in size from about 0.7 (k = 8) to
  0.456 (k = 30): restrictions of good large signings are relatively
  poor small signings; the reversed amplification route fails on
  conference matrices.
- Greedy or beam one-vertex extension from the exact F(10) optimizer
  stalls near ratio 0.52 by n = 14, worse than conference
  restrictions.

Live proof target, sharpened. The July targets are subsumed by one
statement (Conjecture, amplification lemma): there is eta(n) -> 0
with F(N) <= (N/n)^{3/2} F(n) + eta(n) N^{3/2} for all N >= n^2.
This implies convergence. Comparisons at any fixed ratio (for example
proving G(2n) - G(n) -> 0) provably do NOT suffice: the barrier
profile satisfies them. Nonconvergence would need a lower-bound
mechanism active only at special scales; no candidate mechanism is
visible, all algebraic families exist at log-dense orders, and the
subsequence corollary smooths those out.

Metadata corrections (2026-08-01): the question was asked by Paata
Ivanisvili on 2022-01-16; earlier files credited it to Konrad
Swanepoel (README) and A. Kontorovich, 2021 (paper bibliography).
Both fixed. MathOverflow still shows zero answers as of 2026-08-01.

---

The original ledger follows.

Status date: 2026-07-25

## Problem

For
\[
M_n=\min_{a_{ij}\in\{-1,1\}}
 \max_{x\in\{-1,1\}^n}
\left|\sum_{1\le i<j\le n}a_{ij}x_ix_j\right|,
\]
decide whether \(M_n/n^{3/2}\) has a limit.

No proof of existence or nonexistence has yet survived audit.

## Exact reformulations

Let \(E=\binom n2\), let \(G\) be the graph of negative coefficients, and
let \(\delta(S)\) be a cut of \(K_n\). Switching by the sign vector whose
negative coordinates form \(S\) replaces \(G\) by
\(G\mathbin{\triangle}\delta(S)\). Hence
\[
M_n=\min_G\max_{S\subseteq[n]}
\left|E-2\left|G\mathbin{\triangle}\delta(S)\right|\right|.
\]
Equivalently, \(M_n\) is the least possible maximum deviation from \(E/2\)
of the weights in a coset of the binary cut space of \(K_n\).

If \(B\) is a switching representative with maximum total signed edge
sum \(T\), and \(C\) is its largest signed cut sum, then the extreme
energies in that switching class are \(T\) and \(T-2C\). Thus the class
contributes
\[
\max\{T,\,2C-T\}.
\]

## Proved elementary facts

1. \(M_n\le M_{n+1}\). Restrict a coefficient array on \(n+1\) vertices
   to the first \(n\). For every \(x\in\{-1,1\}^n\), the restricted
   energy is the average of the two extended energies obtained by setting
   the last sign to \(1\) and \(-1\).

2. \(M_{n+1}\le M_n+n\). Extend an optimal array arbitrarily by one
   vertex and use the triangle inequality.

3. The problem is the minimum \(L^\infty\) norm of a homogeneous
   degree-two Walsh polynomial with all level-two Fourier coefficients
   in \(\{-1,1\}\).

4. For a symmetric conference matrix \(A\) of order \(n\),
   \(A^2=(n-1)I\) gives
   \[
   \max_x\left|\frac12x^\mathsf{T}Ax\right|
   \le \frac12n\sqrt{n-1}.
   \]
   This is only an upper bound. It does not show that the normalized
   optimum tends to \(1/2\).

5. A Gaussian-sign argument gives the stronger exact lower bound
   \[
   M_n\ge \frac{n\sqrt{n-1}}{\pi}.
   \]
   The proof and its covariance calculation are given in
   `paper/main.tex` and checked numerically in
   `verify_attempt.py`.

6. A bipartite exposure gives an independent, weaker lower bound. Fix a
   set \(I\) of \(k\) vertices, choose its signs uniformly, and then
   choose every sign outside \(I\) to make its signed sum into \(I\)
   nonnegative. If
   \(\mu_k=\mathbb E|\varepsilon_1+\cdots+\varepsilon_k|\), the resulting
   cross energy has expectation \((n-k)\mu_k\). Flipping every sign
   outside \(I\) reverses the cross energy and preserves both internal
   energies, so one of the two full energies has absolute value at least
   the cross energy. Therefore
   \[
   M_n\ge \max_{0\le k\le n}(n-k)\mu_k,
   \]
   and taking \(k\sim n/3\) gives
   \[
   \liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
   \ge \frac23\sqrt{\frac{2}{3\pi}}.
   \]

7. The Paley construction and the prime number theorem in the progression
   \(1\bmod 4\) give
   \[
   \limsup_{n\to\infty}\frac{M_n}{n^{3/2}}\le\frac12.
   \]
   Indeed, choose a prime \(q\equiv1\bmod4\) with \(q+1\ge n\) and
   \((q+1)/n\to1\), form the symmetric Paley conference matrix of order
   \(q+1\), and take an \(n\)-vertex principal submatrix. Its operator
   norm is at most \(\sqrt q\), so its quadratic sign norm is at most
   \(n\sqrt q/2\).

8. Let
   \[
   D_n=\{(x_ix_j)_{i<j}:x\in\{-1,1\}^n\}
   \cup
   \{-(x_ix_j)_{i<j}:x\in\{-1,1\}^n\}.
   \]
   If \(\rho(D_n)\) is its covering radius in
   \(\{-1,1\}^{\binom n2}\), then
   \[
   M_n=\binom n2-2\rho(D_n).
   \]
   This follows by converting correlation to Hamming distance and
   observing that the added negatives turn absolute correlation into
   nearest-codeword distance.

## Audited walls

### Ordinary block composition

Joining optimal arrays on disjoint vertex blocks introduces a complete
bipartite sign form. Its unavoidable sign norm is of the same
\(N^{3/2}\) order as the target. It cannot be put into a lower-order
error term.

### Rank-one blow-up

Replacing each coefficient by a rank-one \(k\)-by-\(k\) sign block gives
the exact bound
\[
M_{nk}\le k^2M_n+O(nk^2),
\]
whose \(k^2\) scaling is too large.

### Annealed partition function

For a random cross block and inverse temperature \(t\), averaging the
partition function incurs
\((\cosh t)^{nm}\). To make entropy negligible at the \(N^{3/2}\)
ground-state scale requires \(t\sqrt N\to\infty\), while making this
annealed cross cost negligible requires \(t\sqrt N\to0\). The direct
annealed interpolation therefore cannot prove the limit.

### Conference spectral shortcut

Small Paley conference examples do not attain their spectral ceiling:
the exact normalized maxima at orders \(6,14,18\) are approximately
\(0.3402,0.4009,0.4321\), while their spectral ceilings are
\(0.4564,0.4818,0.4859\). A proof equating the two is false.

## Live proof targets

1. Find a quenched interpolation or composition theorem whose cross-block
   cost is absorbed at the same normalized constant, rather than bounded
   independently.
2. Prove an asymptotically lossless tensor or lift theorem for the
   degree-two Walsh sign norm.
3. Alternatively, find two subsequences with separated normalized limits.

The public PDF records the proved partial results and explicitly states that
the limit question remains unresolved.

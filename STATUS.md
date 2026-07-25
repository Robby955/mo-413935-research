# MathOverflow 413935: proof-search ledger

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

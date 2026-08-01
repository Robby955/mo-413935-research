# Independent audit of every claim in this repository

Audit date: 2026-08-01. Method: every proof re-derived from scratch,
every script re-run, every numeric claim recomputed by independent
code, all metadata checked against the primary source.

Verdict keys: PROVED (re-derived, correct), VERIFIED (numeric claim
reproduced by independent computation), CORRECTED (true statement
replaced by a sharper or repaired one), WRONG (factual error, fixed
on this branch).

## Mathematical claims of paper/main.tex (July note)

1. Augmented cut-code identity, F(n) = C(n,2) - 2 rho(D_n)
   (Proposition 1). PROVED. Re-derived via the correlation-distance
   dictionary. Addition: D_n corresponds over GF(2) to the cut code
   extended by the all-ones word, which is again linear; covering
   radii through n = 10 are 0, 0, 1, 3, 5, 6, 9, 12, 16.
   Exhaustive check for n <= 5 in verify_attempt.py re-run: PASSED.

2. Gaussian-sign lower bound F(n) >= n sqrt(n-1)/pi (Theorem 2).
   PROVED. Covariance algebra checked entry by entry:
   variances 1 + t^2(n-1)/n from (A^2)_ii = n-1, covariances
   +-(2t/sqrt n) a_ij + (t^2/n)(A^2)_ij, sign alignment v > 0 constant,
   arcsine derivative bound valid on the closed interval, optimization
   at t = sqrt(n/(n-1)) exact. Numeric verifier re-run: PASSED.
   Context added in the new note: this is the linear slice of a
   variational program whose full value is the oscillation functional.

3. Conference upper bound limsup <= 1/2 (Theorem 3). PROVED.
   Symmetric Paley conference matrices exist at order q+1 for prime
   q = 1 mod 4; C^2 = (q)I checked programmatically at all orders
   used; operator norm of a principal submatrix bounded by the full
   norm; density of usable orders from the prime number theorem in
   arithmetic progressions.

4. Monotonicity F(n) <= F(n+1) <= F(n) + n (Proposition 4). PROVED.

5. Ledger fact 6 (bipartite exposure, liminf >= (2/3) sqrt(2/(3 pi))).
   PROVED, and strengthened: it is the odd part of the even-odd
   inequality of the new note (Lemma 3.1 there), which also carries an
   unused additive even term. Numerically 0.3071, weaker than 1/pi,
   correctly labeled in the ledger.

6. Wall: ordinary block composition loses at leading order for
   balanced blocks. PROVED and now two-sided: the cross cost of any
   two-block signing is between (sqrt(2/pi) - o(1)) n sqrt(m) and
   (1 + o(1)) n sqrt(m); merges are exactly a Theta(n sqrt m)
   mechanism (new note, Sections 3 and 4).

7. Wall: rank-one blow-up gives M_{nk} <= k^2 M_n + O(nk^2).
   VERIFIED as an upper-bound scaling statement. The new note adds
   exact tensor data: A x H_2 doubling moves normalized ratios toward
   roughly 0.47-0.50 instead of preserving them.

8. Wall: annealed partition function window is empty. VERIFIED
   quantitatively: the sandwich needs beta sqrt(N) -> infinity, the
   annealed cross cost beta sqrt(N) -> 0. Both re-derived.

9. Wall: "conference spectral shortcut is false; orders 6, 14, 18 do
   not attain the ceiling." CORRECTED. The nonattainment at those
   orders is forced by sqrt(q) being irrational (the form is an
   integer). At orders with q a perfect square the ceiling IS
   attained: order 10 (Petersen two-graph, value 15) and order 26
   (Paley over GF(25), value 65), both verified exactly with recorded
   witnesses. There is no uniform conference deficit, and the wall as
   stated was an artifact of the sample.

10. Conference maxima 5, 21, 33 at orders 6, 14, 18 and normalized
    ratios/ceilings. VERIFIED by re-run and by an independent
    implementation; extended to orders 10, 26, 30 (15, 65, 75).

## Verification scripts

- verify_attempt.py: re-run, output matches README exactly, corruption
  controls trigger as designed. PASSED.
- check_conference_examples.py: re-run, matches. PASSED.

## Metadata

- README states the problem is "Konrad Swanepoel's". WRONG. The
  question was asked by Paata Ivanisvili (MathOverflow user, question
  413935, dated 2022-01-16). Fixed on this branch.
- paper/main.tex bibliography credits the question to A. Kontorovich,
  2021. WRONG (same correction; the year is 2022). Fixed on this
  branch. The mathematical body of the July note is unaffected.
- MathOverflow state re-checked 2026-08-01: zero posted answers. The
  asker's best lower-bound constant mentioned in comments is 2^{-5/2}
  (via Defant, Mastylo, Perez), which both repository bounds dominate.

## Not previously recorded, now added

- Exact values F(2..10) = 1, 3, 4, 4, 5, 9, 10, 12, 13 (not in OEIS
  as of 2026-08-01).
- One-sided optimum equals floor(n/2) for n <= 10 with the
  all-negative class as unique witness, equivalent to Petersdorf's
  maximum-frustration theorem.
- Oscillation optima Osc(3..10) = 4, 8, 8, 10, 16, 20, 24, 26.
- F(10) = 13 beats the order-10 conference class (15): spectral-
  ceiling attainment and minimax optimality separate at n = 10.

See notes/limit_program.tex for the new results and proofs, and
verification/verify_new_results.py for the independent checks of every
table above.

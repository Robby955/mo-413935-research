# MathOverflow 413935: Min-max quadratic forms of signs

Research notes, exact computations, and reproducible verification for Paata
Ivanisvili's [MathOverflow question](https://mathoverflow.net/questions/413935/min-max-of-a-quadratic-form-of-plus-minus-ones).

For a symmetric zero-diagonal sign matrix `A`, define

```text
Q_A(x) = sum_{1 <= i < j <= n} a_ij x_i x_j,
M(A)   = max_{x in {+-1}^n} |Q_A(x)|,
F(n)   = min_A M(A).
```

The question is whether `F(n) / n^(3/2)` converges. The limit remains open.

## Focused outputs

The project is now split into two narrower manuscripts:

- [Finite structure in a min-max quadratic sign problem](paper/mo-413935-finite-results.pdf)
  covers the exact values through order 15, computational certificate
  boundaries, optimizer non-heredity, the weighted Bellman identity, and the
  pinned framed-Hadamard theorem proving `F(16) <= 30`.
- [Relative-gauge composition for quadratic sign discrepancy](paper/mo-413935-composition-framework.pdf)
  gives the exact labeled block-composition law, its finite-temperature
  alignment form, the audited scalar obstructions, and the precise sufficient
  cross-order theorem that remains open.

The [complete research note](paper/mo-413935-second-attempt.pdf) and its
[source](paper/second_attempt.tex) remain the broad archive. The tag
`research-frontier-2026-08-05` freezes the exploration frontier before this
publication split. Failed approaches have not been deleted.

The active research scope is in [ACTIVE_RESEARCH.md](ACTIVE_RESEARCH.md). A
self-contained prompt for a public frontier model is in
[COMPOSITION_FRONTIER_PROMPT.md](COMPOSITION_FRONTIER_PROMPT.md). The earlier
[broad prompt](FRONTIER_MODEL_PROMPT.md) is retained for provenance.

## Current rigorous frontier

The audited universal bounds are

```text
1/pi <= liminf F(n)/n^(3/2)
     <= limsup F(n)/n^(3/2) <= 1/2.
```

The exact finite sequence is

```text
n:     2  3  4  5  6  7   8   9  10  11  12  13  14  15
F(n):  1  3  4  4  5  9  10  12  13  17  18  20  21  27
```

Direct enumeration stops at order 12, where all 1,018,997,864 graphs on 11
vertices are scored without any pruning filter. Orders 13, 14 and 15 need no
enumeration beyond that. They follow from `F(12) = 18` by a one-vertex lift
and the parity law `F(n) = choose(n,2) mod 2`, which forces every value to
have the parity of the number of terms in the form. Order 15 iterates the lift
three times, through level sizes

```text
|B_12(24)| = 82,502,142   |B_13(24)| = 282,202,131   |B_14(25)| = 1,313,164
```

and an empty level at order 15, so no order-15 signing has maximum at most 25;
the order-14 conference matrix extends to one of maximum 27.

The lower certificates through order 12 are computer-assisted. Their
completeness trusts nauty's isomorph-free graph generation; all stream counts,
digests, producer exits, surviving records, and explicit witnesses are checked
separately. The three enumeration sizes 274,668, 12,005,168 and 1,018,997,864
equal the published counts of non-isomorphic graphs on 9, 10 and 11 vertices.
The order-16 result is only

```text
F(16) <= 30.
```

Combined with `F(15) = 27` this pins order 16 to two values: monotonicity gives
`F(16) >= 27`, and `choose(16,2) = 120` is even, so

```text
F(16) is 28 or 30.
```

Deciding between them needs the complete set of order-15 minimisers, which is
one threshold above the order-15 tower and well beyond its cost.

Up to switching and permutation the minimiser is unique at orders 12, 13 and
14. The unique minimiser at order 14 is the Paley conference matrix, and at
order 13 it is that matrix with one vertex deleted. Order 10 has two
minimisers. Record counts in the raw enumeration output are larger than these,
because one switching class contributes one rooted record per orbit of its
automorphism group on the root.

Other banked results include:

- `F(n) <= F(n+1) <= F(n)+n` and the parity-refined puncturing bound;
- the exact augmented cut-code identity
  `F(n) = choose(n,2) - 2 rho(D_n)`;
- the Gaussian lower bound `F(n) >= n sqrt(n-1) / pi`;
- an exact fixed-density cut-discrepancy equivalence up to `O(n)`;
- an exact energy-weighted covering-radius formula for one-vertex extension;
- finite optimizer non-heredity, including the complete `2+8` obstruction;
- the exact relative-gauge max-plus convolution for arbitrary two-block
  composition;
- an exact finite-temperature conditional-alignment chain whose
  zero-temperature slope is the optimizer-compatible composition gain;
- the framed order-16 construction and matching lower bound 30 inside its
  pinned oriented-Hadamard family.

The broad archive also records the cube and elliptope relaxations, linear cut
code and signed MacWilliams identity, graphon and spectral losses, negative
replica transport obstruction, scalar microcanonical profile, labeled Fourier
occupancy hierarchy, Paley alignment, Hadamard lifts, and all failed
amplification attempts.

## The active wall

Set

```text
H(n) = F(n)^(2/3).
```

A uniform estimate

```text
H(n+k) <= H(n) + H(k) + O((n+k)^(1-delta))
```

for some `delta > 0` would force convergence. The missing theorem is not a
separate bound on the internal blocks and the rectangular cross block. Such a
bound loses the full leading scale because both contributions can cancel on
the same spin assignment.

The surviving exact state is the labeled relative-switching fiber. For a
deficit threshold `s`, let `b_s(g)` count subthreshold product states in gauge
fiber `g`. Then

```text
b_s(g) = 0
```

is exactly the assertion that gauge `g` achieves the desired composition
gain. The next useful result must prove an empty fiber, preferably an
abundance of good fibers, at the power-saving threshold. The complete Fourier
factorization of `b_s` and the mixed four-cycle alignment Hamiltonian are the
current starting points.

Scalar energy profiles, low moments, the weighted Hamming union bound, scalar
negative-replica transport, and graph-plus-rectangular transport have audited
leading-order obstructions. The Paley route still lacks a minimax rigidity
lower bound; on the known dense aligned sequence that lower bound is
equivalent to proving the full limit is `1/2`.

## Repository map

- [Active research program](ACTIVE_RESEARCH.md)
- [Spectral numerics behind the active route](SPECTRAL_NUMERICS.md)
- [Composition-only frontier prompt](COMPOSITION_FRONTIER_PROMPT.md)
- [Finite-results source](paper/finite_results.tex)
- [Composition-framework source](paper/composition_framework.tex)
- [Broad research note](paper/mo-413935-second-attempt.pdf)
- [Independent audit](AUDIT.md)
- [Literature and concept map](LITERATURE.md)
- [Continued proof search](RESEARCH_CONTINUATION.md)
- [Proof-search ledger](RESEARCH_LEDGER.md)
- [Original ledger](STATUS.md)
- [Verification guide](verification/README.md)

## Reproduce the focused checks

Use Python 3.10 or newer. The following checks use exact integer or rational
arithmetic for their pass/fail decisions:

```bash
python3 verification/verify_attempt.py
python3 verification/research_exact_small_n.py --max-n 10
python3 verification/research_order11_certify.py
python3 verification/research_order13_certify.py
python3 verification/research_cross_block_composition.py
python3 verification/verify_nonlinear_bellman.py
python3 verification/verify_relative_profile_composition.py
python3 verification/verify_labeled_shell_parseval.py
python3 verification/verify_negative_replica_alignment.py
python3 verification/verify_framed_hadamard_lift_30.py
```

Some exhaustive checks require nauty, NetworkX, NumPy, or `z3-solver`; the
order-13 full-stream scan is expensive. The verification guide gives exact
dependencies, expected output, deterministic seeds, stream digests,
corruption controls, and trust boundaries.

The exact values through order 15 have a separate, cheaper chain. It needs
nauty and NumPy only:

```bash
python3 verification/verify_F_exhaustive.py 10 11
python3 verification/lift.py verification/order12_minimisers.g6 12 18
python3 verification/class_count.py verification/order12_minimisers.g6
```

The first reruns the unpruned enumeration at orders 10 and 11, in seconds and
minutes respectively; adding `12` costs about 20 minutes and settles
`F(12) = 18`. The second lifts the order-12 minimisers and returns 24, which
is what gives `F(13) = 20`. The third shows the two order-12 records are a
single switching class.

Order 15 is the expensive one. It reruns with

```bash
PASSES=8 JOBS=16 bash verification/tower25.sh
```

which needs about 15 GB of scratch and 15 to 20 hours on 16 cores, and prints
an empty level at order 15. Its last level can be rechecked on its own, in
about 20 minutes single-core, with `verification/full_final_check.py`.

The framed order-16 construction also has an independent strict-C verifier:

```bash
cc -std=c11 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror \
  verification/verify_framed_hadamard_lift_30.c \
  -o /tmp/verify_framed_hadamard_lift_30
/tmp/verify_framed_hadamard_lift_30
```

Build the two focused manuscripts with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -jobname=mo-413935-finite-results \
  -output-directory=paper paper/finite_results.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -jobname=mo-413935-composition-framework \
  -output-directory=paper paper/composition_framework.tex
```

## Provenance

Rob Sneiderman directs and preserves the project. Proofs, computational
certificates, nonclaims, and failed routes are retained for independent audit.

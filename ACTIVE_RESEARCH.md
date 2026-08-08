# Active research program

This file records the deliberately narrowed research program after the broad
exploration phase. The complete proof-search history remains in
`RESEARCH_LEDGER.md`, `RESEARCH_CONTINUATION.md`, and
`paper/second_attempt.tex`. The tag `research-frontier-2026-08-05` freezes that
broad frontier.

## Objective

The long-term question is whether

```text
F(n) / n^(3/2)
```

converges. That remains open, and this phase does not target it directly. The
reason is recorded in the archive: a divergent profile satisfies every value
inequality proved so far, so no local estimate can decide convergence, and the
standard existence machinery for a limit of this kind averages over disorder,
whereas here the minimum over signings selects the disorder adversarially.

The active target is instead the strict inequality

```text
limsup F(n)/n^(3/2) < 1/2.
```

It admits partial credit: any explicit `delta > 0` would be the first movement
on the interval `[1/pi, 1/2]` since the question was posed. It is also the side
where the obstruction is now precisely located, which the previous target was
not.

## Why this target, and where it is blocked

`SPECTRAL_NUMERICS.md` records the measurements. In outline:

- The bound `1/2` is exactly the spectral ceiling `n*sqrt(n-1)/2`, which is
  also the value of the semidefinite relaxation for a conference matrix. No
  spectral or relaxation argument can go below it.
- A flat-spectrum matrix with no arithmetic structure reaches only about
  `0.93` of that ceiling, and a first-moment count caps it at
  `sqrt(15)/4 = 0.9682`. Realising that behaviour with `+-1` entries would give
  a constant near `0.484`.
- A `+-1` matrix with a flat spectrum is a conference matrix, and the
  constructible conference matrices are arithmetic. Paley reaches at least
  `0.978` of its ceiling at order 102 and tends to `1`, so it **violates** the
  generic cap rather than obeying it.

So the quasirandom heuristic that predicts a sub-`1/2` constant fails on the
only family that can be written down. That is the wall.

## Acceptable next results

A new result belongs in the active program only if it does at least one of the
following:

1. Exhibits an explicit infinite family of signings with
   `M(A) <= (1/2 - delta) n^(3/2)` for a fixed `delta > 0`, with the upper
   bound proved rather than searched. Local search cannot supply this: an
   incomplete maximisation understates `M`, which is the wrong direction.
2. Proves an alignment bound `m(A) <= 1 - delta` for an explicit conference
   family other than Paley, for instance a Mathon family, where the
   least-nonresidue mechanism that drives Paley to its ceiling is absent.
3. Proves that near-flat signings, with spectral norm `(1+eps)sqrt(n)` and
   suppressed alignment, beat the flat ones asymptotically. The order-10
   optimum does exactly this at finite order: norm `1.311*sqrt(n-1)`,
   alignment `0.661` against the Nesterov floor `0.6366`, beating the
   conference matrix 13 to 15.
4. Improves the lower bound above `1/pi`, which means beating Nesterov
   rounding on flat-spectrum instances.

Finite computation remains useful when it tests a proposed family or supplies
a counterexample. Growing the exact-value table by itself is not an active
asymptotic objective.

## Banked, not active

The following are retained as proved results, finite models, or closed routes,
but are not independent research programs in this phase:

- the power-saving near-subadditivity target
  `H(n+m) <= H(n)+H(m)+O((n+m)^(1-delta))` for `H(n) = F(n)^(2/3)`, together
  with the labeled relative-gauge composition law and the empty-fiber
  condition `min_g b_s(g) = 0` that would have delivered it. The law is exact
  and is written up in `paper/composition_framework.tex`; what is banked is the
  attempt to close it, after a long sequence of leading-order obstructions;
- exact values through order 15, and `F(16)` in `{28, 30}`;
- the weighted Bellman identities and finite optimizer non-heredity;
- scalar cube, elliptope, graphon, spectral, cumulant, and moment relaxations;
- scalar negative-replica transport and graph-plus-rectangular transport;
- the weighted Hamming-ball union-bound composition certificate;
- isolated Hadamard lifts without an iterable framed state;
- Paley alignment without a minimax rigidity lower bound.

These results are documented in the archive and in the focused finite-results
manuscript. They should not be reproved unless a materially different mechanism
feeds the active target.

## Publication split

- `paper/finite_results.tex` contains the exact finite results, certificate
  boundaries, Bellman state, non-heredity, and framed order-16 theorem.
- `paper/composition_framework.tex` contains the exact relative-gauge
  composition law, its scalar and labeled refinements, the proved
  obstructions, and the settling lemma.
- `paper/second_attempt.tex` remains the complete research archive. It is not
  being shortened or overwritten.

No manuscript claims that the MathOverflow limit problem is solved.

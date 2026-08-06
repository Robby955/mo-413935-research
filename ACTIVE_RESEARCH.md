# Active research program

This file records the deliberately narrowed research program after the broad
exploration phase. The complete proof-search history remains in
`RESEARCH_LEDGER.md`, `RESEARCH_CONTINUATION.md`, and
`paper/second_attempt.tex`. The tag `research-frontier-2026-08-05` freezes that
broad frontier.

## Objective

Prove that

```text
F(n) / n^(3/2)
```

converges, or isolate a comparably decisive obstruction. The limit is not
currently known to exist.

The convergence target used in this phase is

```text
H(n + m) <= H(n) + H(m) + O((n + m)^(1 - delta)),
H(n) = F(n)^(2/3),
```

for some fixed `delta > 0`. This power-saving near-subadditivity would force
the convergence of `H(n) / n`, and hence of `F(n) / n^(3/2)`.

## Only active asymptotic route

The active route is the exact labeled relative-gauge composition law.

For fixed internal signings `A`, `B`, and rectangular cross signing `C`, the
product of their local state spaces maps onto a relative-switching group with
equal fibers. The maximum of the full signing in each gauge is exactly the
maximum local energy in the corresponding fiber. Equivalently, its gain below
the independent ceiling is the minimum total local deficit in that fiber.

At a proposed target deficit `s`, let `b_s(g)` count subthreshold product
triples in the gauge fiber `g`. The exact task is to prove that a suitable
choice of near-optimal blocks and cross seed has

```text
min_g b_s(g) = 0
```

at the deficit needed for the power-saving inequality above. A stronger and
more stable theorem would prove that exponentially many fibers are empty or
good.

The next acceptable theorem must therefore be one of the following:

1. A direct character-sum or Fourier argument proving an empty relative-gauge
   fiber at the required threshold.
2. A quantitative abundance theorem for good fibers, strong enough to survive
   repeated proportional composition with a summable error.
3. A closed variational limit for the full relative-switching distribution
   that implies the same power-saving composition inequality.

The first mixed four-cycle Hamiltonian and the complete labeled-shell Fourier
factorization are concrete starting states. They are useful only if they lead
to a cross-order estimate.

## Acceptance test for new work

A new result belongs in the active program only if it does at least one of the
following:

- proves or quantitatively advances the displayed near-subadditivity bound;
- proves a nontrivial uniform bound on the labeled fiber minimum or on the
  abundance of good gauges;
- supplies an iterable invariant with a proved composition rule and
  `o(N^(3/2))` loss;
- proves nonconvergence by matching infinite-family upper and lower bounds.

Finite computation remains useful when it distinguishes candidate labeled
states or tests a proposed theorem. Growing the exact-value table by itself is
not an active asymptotic objective.

## Banked, not active

The following are retained as proved results, finite models, or closed routes,
but are not independent research programs in this phase:

- exact values through order 14 and the order-16 upper bound `F(16) <= 30`;
- the weighted Bellman identities and finite optimizer non-heredity;
- scalar cube, elliptope, graphon, spectral, cumulant, and moment relaxations;
- scalar negative-replica transport and graph-plus-rectangular transport;
- the weighted Hamming-ball union-bound composition certificate;
- isolated Hadamard lifts without an iterable framed state;
- Paley alignment without a minimax rigidity lower bound.

These results are documented in the archive and in the focused finite-results
manuscript. They should not be reproved unless a materially different mechanism
feeds the active cross-order target.

## Publication split

- `paper/finite_results.tex` contains the exact finite results, certificate
  boundaries, Bellman state, non-heredity, and framed order-16 theorem.
- `paper/composition_framework.tex` contains the exact relative-gauge
  composition law, its scalar and labeled refinements, the proved
  obstructions, and the settling lemma.
- `paper/second_attempt.tex` remains the complete research archive. It is not
  being shortened or overwritten.

No manuscript claims that the MathOverflow limit problem is solved.

# MathOverflow 413935: Min-max quadratic forms of signs

Research notes, exact computations, and reproducible verification for Paata
Ivanisvili's [MathOverflow question](https://mathoverflow.net/questions/413935/min-max-of-a-quadratic-form-of-plus-minus-ones).

## Current frontier

The question asks whether `F(n) / n^(3/2)` converges. The limit remains open,
with the following audited bounds and identities:

```text
1/π ≤ liminf F(n)/n^(3/2) ≤ limsup F(n)/n^(3/2) ≤ 1/2
F(n) = choose(n, 2) - 2ρ(D_n)
F(n) ≤ F(n + 1) ≤ F(n) + n
```

The repository now goes substantially beyond those baseline results:

- The direct coupling-cube relaxation has the unique optimizer zero and a
  leading-order integrality gap. It cannot be rounded into a solution.
- The adversarial elliptope relaxation has the exact normalized limit `1/2`.
- The augmented cut code is linear. Its dual consists of even-cardinality
  Eulerian subgraphs, producing an exact signed MacWilliams formula.
- Ordinary graphon and empirical spectral limits provably erase the relevant
  `n^(3/2)` information.
- Exhaustive computation gives
  `F(2), ..., F(10) = 1, 3, 4, 4, 5, 9, 10, 12, 13`.
  Exact extension of the two optimal order-10 classes gives the additional
  bound `13 ≤ F(11) ≤ 17`; the exact order-11 value remains unknown.
- Exact hereditary cavity inequalities strengthen the induced-submatrix lower
  bounds through Walsh orthogonality and block pairing. Their universal gains
  are subleading, so they do not settle convergence.
- Retaining the nonlinear term in the Gaussian arcsine argument gives an
  exact trace-four bound. A parity theorem for its Gram defect gives the new
  finite consequences `F(20) ≥ 30` and `F(21) ≥ 32`. Optimizing the full
  one-parameter certificate still has universal asymptotic constant exactly
  `1/π`, so this entire Gaussian mechanism cannot settle convergence by
  itself.
- Exact obstruction results rule out several natural composition mechanisms,
  including edge-separable saturation, near-saturated cross blocks,
  bounded-rank local lifts, and canonical Seidel/Kronecker amplification.
- Complete cross-block optimization exposes finite non-heredity. No optimal
  order-10 signing contains an optimal order-8 principal submatrix: completing
  optimal blocks of sizes 2 and 8 forces maximum 15, although `F(10) = 13`.
  Allowing the order-8 block maximum to rise from 10 to 12 is necessary and
  sufficient to recover 13.
- One-vertex extension has an exact energy-weighted projective
  covering-radius formula and Bellman identity. The ordinary covering radius
  of the exact maximizers explains every order-9 and order-10 extension class,
  while two order-7 classes prove that lower energy layers cannot generally be
  discarded. Complete Bellman Pareto frontiers through residual order 8 show
  exact trades between internal energy and covering deficit.
- Two optimal order-9 classes have the same complete absolute-energy
  histogram, hence identical scalar partition functions at every temperature,
  but extension values 13 and 15. Scalar free energy alone is therefore not a
  closed cavity state.
- The one-vertex Bellman identity extends exactly to arbitrary two-block
  composition as a weighted covering radius of the projective rank-one code.
  The associated weighted union bound is now proved incapable of reaching
  the required composition constant. Retaining the switching-orbit
  dependence gives an exact weighted sumset noncoverage formulation instead.
- A complete pass over all 12,346 root-normalized order-9 signings gives the
  full Bellman frontier `{(12, 0)}`. A stronger collision proves that the
  energy histogram plus the pair-distance law of exact maximizers still does
  not determine extension. The complete energy-coloured two-point law does
  separate all deficits at order 9, a finite result only.
- The complete order-10 catalogue has an exact three-phase finite-temperature
  minimizer: a conference signing at high temperature, an intermediate
  non-ground-state signing, and an `F(10)=13` signing at low temperature.
  Thus the adversarial optimizer genuinely changes with temperature.
- Density-one control of Bellman-optimal predecessors is not the decisive
  scalar wall. An explicit nonconvergent countermodel satisfies all current
  scalar cross-order inequalities even with `O(sqrt(n))` Bellman cost at every
  order. The precise one-vertex target is stabilization of the normalized
  dyadic Bellman cost, or a corresponding Bellman-Cesaro law.

## Main open routes

The finite-temperature route is existence of the extensive minimax
free-energy limit, with the negative-replica formulation providing the
sharpest current settling lemma.

The ground-state route is a power-saving near-subadditivity theorem for
`H(n) = F(n)^(2/3)`. The finite non-heredity result shows that the construction
must use a composable family of near-optimal blocks rather than arbitrary
exact minimizers. A uniform defect of order `O((n + m)^(1 - δ))`, for some
positive `δ`, would force `H(n) / n` and therefore `F(n) / n^(3/2)` to
converge.

The present wall is optimizer-sensitive composition. Internal block energy
and the cross field can cancel for the same spin assignment, while bounding
them separately costs the full leading order.

For one-vertex growth, the exact state is now known: the Pareto profile of the
internal maximum and its energy-weighted covering deficit. The next question
is whether that profile has a uniform asymptotic law over the near-optimal
energy window. Exact ground states alone are provably insufficient.

The exact multivertex state sharpens this route, but its weighted-entropy
union bound has a universal balanced floor
`2 sqrt(log 2) N^(3/2)`, strictly above the required asymptotic ceiling
`sqrt(2) N^(3/2)`. The live replacement is a weighted noncoverage theorem for
the row-and-column switching orbit of a low-norm cross seed. This is a
max-plus convolution problem on the projective rank-one group and retains the
higher-order dependence lost by first moment estimates.

## Research package

- [Frontier-model research prompt](FRONTIER_MODEL_PROMPT.md)
- [Continued proof search](RESEARCH_CONTINUATION.md)
- [Main research note](paper/mo-413935-second-attempt.pdf)
- [Main note source](paper/second_attempt.tex)
- [Independent claim audit](AUDIT.md)
- [Literature and concept map](LITERATURE.md)
- [Proof-search ledger](RESEARCH_LEDGER.md)
- [Original proof-search ledger](STATUS.md)
- [Original research note](paper/mo-413935-ai-attempt.pdf)
- [Verification guide](verification/README.md)

Failed routes are retained with their exact obstruction rather than removed
from the research record.

## Reproduce the checks

Use Python 3.10 or newer. The core checks use only the standard library.

```bash
python3 verification/verify_attempt.py
python3 verification/check_conference_examples.py
python3 verification/verify_new_results.py
python3 verification/verify_continuation.py
python3 verification/verify_coding_continuation.py
python3 verification/verify_amplification_obstructions.py
python3 verification/verify_cavity_hereditary.py
python3 verification/verify_nonlinear_bellman.py
python3 verification/verify_frontier_walls.py
```

The exact cross-block research check requires nauty `geng`, NetworkX, and
`z3-solver`. The order-9 geometry check requires nauty, NetworkX, and NumPy;
the order-10 temperature check requires nauty and NumPy.

```bash
python3 verification/research_cross_block_composition.py
python3 verification/research_order9_weighted_geometry.py
python3 verification/research_order10_temperature.py
```

The exhaustive search through order 10 requires nauty `geng`; NetworkX adds
independent decoding and switching-class checks.

```bash
python3 verification/research_exact_small_n.py \
  --max-n 10 \
  --networkx-crosscheck \
  --classify-switching-optima \
  --strict-stream-digests
```

The independent trusted-solver check for orders 7, 8, and 9 requires
`z3-solver`.

```bash
python3 verification/research_z3_certify.py
```

See the [verification guide](verification/README.md) for expected output,
certificate boundaries, dependencies, seeds, and the order-10 Z3 timeout.

## Provenance

Rob Sneiderman directs and preserves the project. The research notes,
proof-search records, and verification code retain exact nonclaims and are
designed for independent replay.

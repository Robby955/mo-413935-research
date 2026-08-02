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
- Exact finite computation gives
  `F(2), ..., F(14) = 1, 3, 4, 4, 5, 9, 10, 12, 13, 17, 18, 20, 21`.
  The order-11 lower certificate scans all 12,005,168 unlabeled residual
  graphs and exactly filters 2,153,606 eligible classes for cut evaluation.
  Direct generation of that reduced stream reproduces the zero-survivor result
  with the same evaluator.
  This is a computer-assisted theorem whose completeness trusts nauty. The
  order-12 upper bound is an explicit signing; its lower bound follows from
  order 11 by puncturing and parity.
- The order-13 lower certificate hashes and scans all 1,018,997,864 unlabeled
  residual graphs on eleven vertices in eight fail-closed nauty shards. Exactly
  two rooted records survive the order-12 threshold; direct enumeration gives
  maximum 18 and one-vertex extension minimum 24 for both. All other records
  have maximum at least 20. A principal submatrix of the order-14 Paley matrix
  supplies the matching upper bound, and parity then gives `F(14) = 21`.
  Completeness has the same explicit nauty trust boundary as the order-11
  certificate. This also makes the order-14 conference matrix optimal, but
  the certificate alone does not establish uniqueness.
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
- The annealed-normalized negative replica is exactly supermultiplicative at
  fixed exponent and temperature. Boolean reverse hypercontractivity follows
  the parameter curve needed by block composition, but the proposed
  power-saving reverse comparison is false: at
  `(beta, theta) = (4, 8)` its defect has rigorous liminf at least
  `sqrt(2) - 2 log(2) = 0.027919...` after division by `n^2`.
- The transport defect has an exact entropy-production integral and a
  quantitative stability remainder. The full block law instead factors
  through a conditional relative-switching alignment term. Its
  zero-temperature slope is exactly
  `F(n) + F(k) + B_square(n,k) - F(n+k)`, so it is provably a leading-order
  state, not an entropy correction that can be dropped.
- The first alignment Hamiltonian is an explicit mixed four-cycle trace.
  Exact `2 + 4` examples have identical graph and rectangular scalar
  partition curves but different conditional alignment moments. Thus even all
  three local scalar free-energy curves do not close the composition state.
  The exceptional twofold gauge redundancy at block size 2 is explicitly
  accounted for in the verifier.
- The relative-gauge law is an exact convolution of three local Gibbs laws.
  At zero temperature this yields a deterministic microcanonical composition
  theorem: the `2^(n+k-1)`-st smallest sum of the three local energy deficits
  is a guaranteed composition gain. This guarantee is optimal if the additive
  labels are discarded. It replaces optimizer-count heuristics by a precise
  near-maximal-profile large-deviation target and has no known analogue of the
  failed weighted union bound's universal leading floor.
- For every odd prime power `m`, the square-order Paley conference matrix of
  order `m^2+1` has a Boolean eigenvector and attains its spectral ceiling.
  This is the known regular-conference construction, not a new result and not
  a statement about `F(m^2+1)`. Already `F(10)=13` while that Paley matrix has
  maximum 15.
- For every symmetric Paley order `q+1`, an exact Fourier-leakage identity
  reduces Boolean spectral alignment to finding an asymptotically balanced
  sign function whose additive Fourier energy has `o(q)` mass in one
  quadratic-character half. For prime `p`, a half-interval sign function has
  leakage at most `2p/(ell(p)-1)`, where `ell(p)` is the least quadratic
  nonresidue. The prime number theorem in fixed progressions produces a
  multiplicatively dense prime sequence with `ell(p) -> infinity`, proving
  asymptotic Paley spectral alignment on that sequence. Exact zero leakage is
  impossible for nonconstant Boolean functions at prime order. Source-built
  exhaustive scans give `M(C_6), M(C_14), M(C_18), M(C_30) = 5, 21, 33, 75`.
  The Paley alignment theorem is still only about a selected matrix family.
  A separate dense-order minimax rigidity lower bound is required to say
  anything decisive about the limit of `F(n)/n^(3/2)`.
- Density-one control of Bellman-optimal predecessors is not the decisive
  scalar wall. An explicit nonconvergent countermodel satisfies all current
  scalar cross-order inequalities even with `O(sqrt(n))` Bellman cost at every
  order. The precise one-vertex target is stabilization of the normalized
  dyadic Bellman cost, or a corresponding Bellman-Cesaro law.

## Main open routes

The finite-temperature route is existence of the extensive minimax
free-energy limit. Scalar negative-replica parameter transport is now closed:
its defect is leading order in a rigorous parameter range. The surviving
state is the conditional distribution on relative block switchings. A useful
next theorem must control that alignment free energy with a summable defect;
simply writing the required almost-superadditivity in alignment coordinates
does not make it weaker.

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

The new scalar intermediate target is the microcanonical profile theorem:
construct near-optimal graph blocks and a cross seed for which fewer than
`2^(n+k-1)` product triples lie below the deficit needed for power-saving
near-subadditivity. Constantly many exact maximizers do not suffice; the
entire leading-scale near-maximal profile must be thin. If this quantile is
too small, the additive labels in the exact relative-gauge convolution are
the remaining geometric fallback.

A separate value-specific route now has only one mathematical wall. The
least-nonresidue interval theorem supplies an asymptotically balanced Paley
sign function with Fourier leakage `o(p)` on a multiplicatively dense prime
sequence. What remains is minimax rigidity
`F(p+1) >= (1-o(1)) M(C_{p+1})` on that sequence. The finite optimality of
`C_6` and `C_14` does not establish this lower bound, and `F(10)<M(C_10)`
shows that exact conference optimality is false at some orders.

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
python3 verification/verify_paley_subfield.py
python3 verification/verify_negative_replica_transport_obstruction.py
python3 verification/verify_negative_replica_alignment.py
python3 verification/research_order13_certify.py
```

The order-30 Paley alignment scan is a self-contained C check:

```bash
cc -std=c11 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror \
  verification/research_paley_alignment.c -lm \
  -o /tmp/research_paley_alignment
/tmp/research_paley_alignment
```

The exact cross-block research check requires nauty `geng`, NetworkX, and
`z3-solver`. The order-9 geometry check requires nauty, NetworkX, and NumPy;
the order-10 temperature check requires nauty and NumPy.

```bash
python3 verification/research_cross_block_composition.py
python3 verification/research_order9_weighted_geometry.py
python3 verification/research_order10_temperature.py
python3 verification/research_negative_replica_transport.py
python3 verification/research_order11_certify.py
```

The full order-13 lower certificate scans 1,018,997,864 residual graphs and
requires nauty `geng`:

```bash
python3 verification/research_order13_certify.py \
  --full-stream --jobs 8 --geng /absolute/path/to/geng
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

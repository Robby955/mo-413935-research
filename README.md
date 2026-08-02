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
- Exact hereditary cavity inequalities strengthen the induced-submatrix lower
  bounds through Walsh orthogonality and block pairing. Their universal gains
  are subleading, so they do not settle convergence.
- Exact obstruction results rule out several natural composition mechanisms,
  including edge-separable saturation, near-saturated cross blocks,
  bounded-rank local lifts, and canonical Seidel/Kronecker amplification.

## Main open routes

The finite-temperature route is existence of the extensive minimax
free-energy limit, with the negative-replica formulation providing the
sharpest current settling lemma.

The ground-state route is a power-saving near-subadditivity theorem for
`H(n) = F(n)^(2/3)`. A uniform defect of order
`O((n + m)^(1 - δ))`, for some positive `δ`, would force
`H(n) / n` and therefore `F(n) / n^(3/2)` to converge.

The present wall is optimizer-sensitive composition. Internal block energy
and the cross field can cancel for the same spin assignment, while bounding
them separately costs the full leading order.

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

Rob Sneiderman directs the project. OpenAI Codex produced the research notes,
proof searches, and verification code, with the finite checks designed for
independent replay.

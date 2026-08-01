# An AI Attempt at MathOverflow 413935

OpenAI Codex attempted Paata Ivanisvili's
[min-max problem for quadratic forms of signs](https://mathoverflow.net/questions/413935/min-max-of-a-quadratic-form-of-plus-minus-ones)
(MathOverflow 413935, January 2022).
This repository contains the resulting research note, its source,
reproducible checks, and a second-pass research program.

## Status: not solved

The original question asks whether

$$
\lim_{n\to\infty}\frac{F(n)}{n^{3/2}}
$$

exists. This attempt does not prove existence or nonexistence, and must not
be cited or submitted as a solution.

## Partial result

The note derives

$$
\frac1\pi
\le
\liminf_{n\to\infty}\frac{F(n)}{n^{3/2}}
\le
\limsup_{n\to\infty}\frac{F(n)}{n^{3/2}}
\le
\frac12
$$

and the exact augmented cut-code identity

$$
F(n)=\binom n2-2\rho(D_n).
$$

These claims have computational consistency checks but have not received
independent human mathematical review.

## August 2026 second pass

A follow-up program audited every claim above (all mathematical claims
stand; two attribution errors were fixed) and added new results:
exact values F(2..10) = 1, 3, 4, 4, 5, 9, 10, 12, 13; an exactly
superadditive oscillation functional that identifies the reach of the
Gaussian-rounding method; a merge theorem with matching block floor
and its continuity corollaries (the limit exists iff it exists along
conference orders); conference spectral-ceiling attainment at orders
10 and 26; the collapse of box, sphere, and graphon relaxations; and
an explicit divergent profile satisfying every proved inequality,
which shows local estimates cannot decide the question. The problem
remains open.

- [Second research note (PDF)](notes/limit_program.pdf)
- [Second note source](notes/limit_program.tex)
- [Claim-by-claim audit](AUDIT.md)
- [New verification suite](verification/verify_new_results.py)
- [Machine-readable results and witnesses](verification/data/results_2026_08.json)
- Exhaustive enumerators: `verification/exact_fn.c`,
  `verification/conference_max.c`, `verification/ext_chain.c`

## Files

- [Research note (PDF)](paper/mo-413935-ai-attempt.pdf)
- [LaTeX source](paper/main.tex)
- [Proof-search ledger](STATUS.md)
- [Main verification script](verification/verify_attempt.py)
- [Conference-matrix checks](verification/check_conference_examples.py)

## Reproduce the checks

The scripts require Python 3 and use only the standard library.

```bash
python3 verification/verify_attempt.py
python3 verification/check_conference_examples.py
python3 verification/verify_new_results.py
```

Expected main-verifier output:

```text
gaussian_edges_checked=8550
covering_radius_identity=VERIFIED_orders_2_through_5
optimized_lower_bound=VERIFIED_orders_2_through_9999
corruption_controls=PASSED
```

These computational checks audit finite instances and algebra used in the
writeup. They are not a substitute for mathematical review.

## Transcript and attribution

The PDF and mathematical argument were produced by OpenAI Codex. Rob
Sneiderman directed the attempt, preserved the artifacts, and is publishing
them without claiming a solution.

The PDF is a distilled research note, not the complete conversation
transcript. The challenge requests the full AI transcript, so an exported
transcript should be shared separately with any contest submission.

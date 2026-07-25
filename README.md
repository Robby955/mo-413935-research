# An AI Attempt at MathOverflow 413935

OpenAI Codex attempted Konrad Swanepoel's
[min-max problem for quadratic forms of signs](https://mathoverflow.net/questions/413935/min-max-of-a-quadratic-form-of-plus-minus-ones).
This repository contains the resulting research note, its source, and
reproducible checks.

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

# An AI attempt at MathOverflow 413935

This repository records an incomplete attempt by OpenAI Codex at Konrad
Swanepoel's [min-max problem for quadratic forms of
signs](https://mathoverflow.net/questions/413935/min-max-of-a-quadratic-form-of-plus-minus-ones).

## Status: not solved

The attempt does not prove that
\[
\lim_{n\to\infty}\frac{F(n)}{n^{3/2}}
\]
exists, and it does not prove that the limit fails to exist. It must not be
cited or submitted as a solution.

The note proves the partial bounds
\[
\frac1\pi
\le
\liminf_{n\to\infty}\frac{F(n)}{n^{3/2}}
\le
\limsup_{n\to\infty}\frac{F(n)}{n^{3/2}}
\le
\frac12
\]
and gives the exact augmented cut-code identity
\[
F(n)=\binom n2-2\rho(D_n).
\]

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

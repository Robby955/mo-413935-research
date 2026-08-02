# Independent audit of the original attempt

Audit date: 2026-08-01
Audited revision: `b5b966a45a62a053cfe90edf3befce6a14984030`

## Scope and conclusion

Every tracked source file, verification script, and the rendered PDF at the
audited revision was read.  The manifest was checked before any edits, the two
original verification scripts were rerun, and the original TeX source was
rebuilt independently.

The central mathematical claims survive audit:

\[
\frac1\pi\le \liminf_n\frac{F(n)}{n^{3/2}}
\le \limsup_n\frac{F(n)}{n^{3/2}}\le\frac12,
\qquad
F(n)=\binom n2-2\rho(D_n),
\]

and \(F(n)\le F(n+1)\le F(n)+n\).  None of them proves convergence.
The repository contained two bibliographic errors and one factor-of-two prose
error, recorded below.  The descriptions of failed proof mechanisms are valid
only for the displayed mechanisms; they are not impossibility theorems.

The status labels in this ledger mean:

- **proved:** a complete proof was independently checked;
- **flawed:** the statement as written is false or factually wrong;
- **incomplete:** useful reasoning is present, but it does not establish the
  broader statement that its wording might suggest;
- **computational only:** exact finite computation, not an asymptotic proof;
- **redundant:** correct but already follows from another audited claim.

## Claim ledger

| Claim | Status | Proof or check location | Audit finding |
|---|---|---|---|
| Switching gives \(Q_A(x)=E-2|G\mathbin\triangle\delta(S)|\). | proved | `STATUS.md`, “Exact reformulations” | Direct expansion is exact. Switching and vertex relabeling preserve \(M(A)\). |
| “\(M_n\) is the least possible maximum deviation from \(E/2\).” | flawed | `STATUS.md`, original line 27 | The energy is \(2(E/2-w)\), so \(M_n\) is **twice** the least possible maximum weight deviation. The formula immediately before it is correct. |
| A switching class with maximum total signed sum \(T\) and maximum signed cut sum \(C\) contributes \(\max\{T,2C-T\}\). | proved | `STATUS.md`, “Exact reformulations” | After choosing a representative attaining \(T\), switching over a cut subtracts twice its signed cut sum. Taking both extrema gives the formula. |
| \(F(n)\le F(n+1)\). | proved | `paper/main.tex`, Proposition “Elementary regularity” | A restricted energy is the average of its two one-vertex extensions; absolute value and minimization are used in the correct order. |
| \(F(n+1)\le F(n)+n\). | proved | same proposition | For the new linear term \(L\), \(\max_{s=\pm1}|Q+sL|=|Q|+|L|\le F(n)+n\). |
| \(F(n)\) is the minimum \(L^\infty\)-norm of a degree-two Walsh polynomial with all level-two coefficients in \(\{\pm1\}\). | proved | `STATUS.md`, item 3 | This is exactly the definition after identifying \(x_ix_j\) with Walsh characters. |
| Conference spectral bound \(M(C)\le n\sqrt{n-1}/2\). | proved | `paper/main.tex`, Theorem “conference-matrix upper bound” | Cauchy--Schwarz and \(C^2=(n-1)I\) give an upper bound only; no Boolean attainment is assumed. |
| Gaussian-sign lower bound \(F(n)\ge n\sqrt{n-1}/\pi\). | proved | `paper/main.tex`, Theorem “Gaussian-sign lower bound” | The covariance calculation, arcsine-law step, factor of two between two expectations, and optimization in \(t\) were all rederived. Endpoint correlations are handled by continuity. |
| Bipartite exposure bound \(F(n)\ge\max_k(n-k)\mathbb E|\sum_{i=1}^k\varepsilon_i|\). | proved | `STATUS.md`, item 6 | Conditioning on the exposed block and flipping every outside spin gives the claimed absolute-energy lower bound. Its asymptotic constant is weaker than \(1/\pi\). |
| Paley principal submatrices imply \(\limsup F(n)/n^{3/2}\le1/2\). | proved | `paper/main.tex`, Theorem “conference-matrix upper bound” | Symmetric Paley conference matrices at prime orders \(q+1\), compression of operator norm, and primes \(q\equiv1\pmod4\) with next-order ratio tending to one suffice. |
| \(F(n)=\binom n2-2\rho(D_n)\) for the augmented cut code. | proved | `paper/main.tex`, Proposition “Augmented cut-code identity” | Correlation equals \(m-2d_H\); antipodal augmentation turns maximum absolute correlation into nearest-codeword correlation. |
| The interval \(1/\pi\le\liminf\le\limsup\le1/2\). | redundant | combination of the two audited theorems | Correct, but it is only a collation of the lower and upper bounds. |
| Direct random cross-block composition has a leading \(N^{3/2}\) union-bound cost. | proved | `paper/main.tex`, “composition wall” | The displayed subgaussian union bound has size \(\sqrt{2\log2\,nm(n+m)}\), which is leading order for comparable blocks. |
| Ordinary block composition cannot prove convergence. | incomplete | `STATUS.md`, “Ordinary block composition” | The calculation rules out the stated independent-cross-block estimate. It does not rule out correlated, optimized, or quenched composition mechanisms. |
| Rank-one blow-up has \(F(nk)\le k^2F(n)+O(nk^2)\). | proved | `STATUS.md`, “Rank-one blow-up” | More precisely, choosing arbitrary internal clone-block signs gives \(F(nk)\le k^2F(n)+n\binom k2\). This specific lift is not asymptotically lossless. |
| Rank-one blow-ups cannot yield amplification. | incomplete | same section | The displayed naive lift fails; the audit does not exclude a different lift with cancellation or additional structure. |
| Direct annealed interpolation has incompatible temperature requirements. | proved | `STATUS.md`, “Annealed partition function” | For the stated independent cross block, \((\cosh t)^{nm}\) gives a leading cross cost when entropy is suppressed. |
| Annealed methods cannot prove convergence. | incomplete | same section | Only the direct factorized annealed estimate is excluded. A quenched interpolation or a different variational representation remains possible. |
| The Paley conference examples of orders \(6,14,18\) have maxima \(5,21,33\). | computational only | `verification/check_conference_examples.py` | Exhaustive Boolean evaluation reproduces the values and confirms strict gaps below the spectral ceilings. |
| The conference spectral ceiling is attained by Boolean spins. | flawed | explicitly rejected in `STATUS.md` | The exact examples above are counterexamples to this tempting shortcut. The repository correctly rejects it. |
| Covering-radius identity verified through \(n=5\). | computational only | `verification/verify_attempt.py` | The exhaustive loops do check both sides, but only for the printed finite range. |
| Gaussian covariance algebra verified on 8,550 edges and optimized formula through \(n=9999\). | computational only | `verification/verify_attempt.py` | Deterministic seed and corruption controls reproduce the output. Floating-point sampling supports, but is not needed for, the analytic proof. |
| The question was posed by Konrad Swanepoel / A. Kontorovich in 2021. | flawed | original `README.md` and `paper/main.tex` bibliography | The live MathOverflow record identifies Paata Ivanisvili and 16 January 2022. The attribution and year are corrected in the current files. |

## Verification boundaries

The original scripts use deterministic inputs where randomness occurs and
contain useful negative controls.  They do not constitute machine-checked
proofs.  In particular:

- the original covering-radius enumeration stops at \(n=5\);
- the Gaussian verifier samples finite matrices and uses floating point;
- the conference verifier proves only the three finite maxima it enumerates;
- `zip(..., strict=True)` means the documented runtime should be Python 3.10
  or newer, not an unspecified Python 3 release;
- the manifest authenticated the original tracked bytes, not the truth of the
  mathematical claims.

The second-attempt scripts add exact-arithmetic identity checks and exhaustive
small-order searches.  Their own scope and dependencies are documented in
`verification/README.md`.

## Post-audit claim additions

These claims were developed after the fixed-revision audit above. They are
listed separately so the original audit boundary remains explicit.

| Claim | Status | Proof or check location | Exact boundary |
|---|---|---|---|
| The weighted projective radius satisfies \(E(B)=M(B)+n-2\rho_{\mathrm w}(B)\), and minimizing it gives the exact Bellman identity for \(F(n+1)\). | proved | `RESEARCH_CONTINUATION.md`, Theorem 15; `paper/second_attempt.tex` | Elementary pairing and projective Hamming-distance proof; no computation is used. |
| Only configurations with \(\lvert Q_B(x)\rvert\ge M(B)-2\lfloor n/2\rfloor\) can affect the weighted radius. | proved | immediately after Theorem 15 | An exact maximizer always supplies weighted distance at most \(\lfloor n/2\rfloor\); deeper configurations cannot minimize. |
| Exact maximizers alone determine every optimal-class extension at orders 9 and 10, but not at order 7. | computational only | `verification/research_cross_block_composition.py` | Complete optimizer catalogues trust nauty stream completeness; all incident signs are directly enumerated. |
| Two optimal order-9 classes have identical absolute-energy histograms, hence identical scalar partition functions at every temperature, but extension values 13 and 15. | computational only | `RESEARCH_CONTINUATION.md`, Proposition 16; same script | The two records are in distinct catalogue classes; histograms, covering radii, and incident-sign extensions are directly enumerated. |
| The full Pareto frontiers are \(\{(5,2),(7,1),(9,0)\}\) at residual order 6 and \(\{(10,1),(12,0)\}\) at residual order 8. | computational only | same script, `full_extension_pareto` records | Every root-normalized residual graph through order 8 is exhausted; graph counts and stream hashes are asserted, and extension values are independently recomputed. |
| Optimal blocks of orders 2 and 8 force combined maximum 15, while order-8 internal maximum 12 is the least budget permitting combined maximum 13. | computational only | `RESEARCH_CONTINUATION.md`, Theorem 14; same script | The critical cross blocks are directly exhausted; completeness of residual representatives trusts nauty. |
| The nonlinear square-covariance formula (17.1), its trace-four stability refinement, and the necessary condition \(\|A_n^2-(n-1)I\|_F=o(n^2)\) at asymptotic constant \(1/\pi\). | proved | `RESEARCH_CONTINUATION.md`, Theorem 17; `paper/second_attempt.tex` | The two Gaussian correlation matrices, exact arcsine difference, convexity, and trace identity were rederived. Finite normalization checks exhaust signings through order 5. |
| The finite improvement \(F(21)\ge32\). | proved | consequence of Theorem 17; `verification/verify_nonlinear_bellman.py` | A rational sine bound proves strictness above 30; energy parity forces the next admissible value 32. Floating-point output is not used for strictness. |
| The multivertex weighted rank-one covering identity for \(J(B,D)\) and the resulting exact two-block Bellman formula for \(F(n+k)\). | proved | `RESEARCH_CONTINUATION.md`, Theorem 18; `paper/second_attempt.tex` | Elementary pairing and projective Hamming distance; 23 small block pairs are independently enumerated. |
| The weighted entropy upper bound \(J\le L+\sqrt{2nk(\Xi+\log4)}\). | proved | `RESEARCH_CONTINUATION.md`, Proposition 19 | Hoeffding plus a finite union bound. Its generic specialization still has leading-order cost. The constant-matching profile law stated after it remains conjectural. |
| Bellman-optimal predecessors have \(\sigma_n,\delta_n=O(g(N)\sqrt N)\) for all but \(O(N/g(N))\) indices in every dyadic interval. | proved | `RESEARCH_CONTINUATION.md`, Proposition 20 | Exact Bellman increment, telescoping, the random-sign upper bound, and Markov. Sparse exceptional orders remain uncontrolled. |
| The complete order-9 Bellman frontier is \(\mathcal B_9=\{(12,0)\}\), with 20/35 root representatives and 4/11 switching classes by deficit. | computational only | `verification/research_order9_weighted_geometry.py` | All 12,346 root-normalized records are evaluated with exact integer arrays. Completeness trusts nauty count/digest; class grouping trusts NetworkX isomorphism. |
| The records `GHOgmo` and `Gxd?Dc` have identical energy histograms and exact-maximizer pair-distance laws but extension values 15 and 17. | computational only | same order-9 verifier | Both records are directly decoded; energies, radii, pair enumerators, and incident-sign extensions are recomputed exactly. |
| The complete energy-coloured two-point law separates every weighted deficit at order 9. | computational only | same order-9 verifier | Exact finite catalogue statement only; no sufficiency is claimed at other orders or asymptotically. |

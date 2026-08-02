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
| The nonlinear square-covariance formula (17.1), its trace-four refinement, and the necessary condition \(\|A_n^2-(n-1)I\|_F=o(n^2)\) at asymptotic constant \(1/\pi\). | proved | `RESEARCH_CONTINUATION.md`, Theorems 17 and 22; `paper/second_attempt.tex` | The two Gaussian correlation matrices, exact arcsine difference, Jensen compression, and trace identity were rederived. The norm condition is only vanishing normalized Gram defect, not conference rigidity. |
| The finite improvements \(F(20)\ge30\) and \(F(21)\ge32\). | proved | Theorems 17 and 21; `verification/verify_nonlinear_bellman.py` | Exact rational inequalities prove both strict thresholds; the energy lattice supplies the parity rounding. Floating-point output is not used for strictness. |
| Optimizing the complete one-parameter nonlinear square-covariance certificate has minimax asymptotic constant exactly \(1/\pi\). | proved | `RESEARCH_CONTINUATION.md`, Theorem 22 | Paley conference principal blocks of order \(n+o(n)\) give a uniform-in-parameter upper barrier; the universal parameter gives the matching lower bound. |
| The multivertex weighted rank-one covering identity for \(J(B,D)\) and the resulting exact two-block Bellman formula for \(F(n+k)\). | proved | `RESEARCH_CONTINUATION.md`, Theorem 18; `paper/second_attempt.tex` | Elementary pairing and projective Hamming distance; 23 small block pairs are independently enumerated. |
| The weighted entropy upper bound \(J\le L+\sqrt{2nk(\Xi+\log4)}\). | proved | `RESEARCH_CONTINUATION.md`, Proposition 19 | Hoeffding plus a finite union bound. The bound itself is correct. |
| The proposed constant-matching use of that entropy bound, equation (19.2). | flawed | `RESEARCH_CONTINUATION.md`, Proposition 23 | The left side has universal floor \(\sqrt{2nk(n+k)\log2}\); in balanced blocks its constant \(2\sqrt{\log2}\) exceeds the required \(\sqrt2\) by a leading amount. |
| The switching-orbit bad set is exactly \(\mathcal B_K=\bigcup_R R C_{K-h(R)}\). | proved | `RESEARCH_CONTINUATION.md`, Proposition 23; `verification/verify_frontier_walls.py` | Elementary group action. It is an exact dependence-sensitive reformulation, not yet a noncoverage theorem. |
| Bellman-optimal predecessors have \(\sigma_n,\delta_n=O(g(N)\sqrt N)\) for all but \(O(N/g(N))\) indices in every dyadic interval. | proved | `RESEARCH_CONTINUATION.md`, Proposition 20 | Exact Bellman increment, telescoping, the random-sign upper bound, and Markov. Theorem 24 shows that even uniform scalar \(O(\sqrt n)\) control would not force convergence. |
| The normalized sequence converges iff its normalized dyadic Bellman cost converges; a scalar countermodel satisfies all current scalar inequalities with \(O(\sqrt n)\) cost but oscillates. | proved | `RESEARCH_CONTINUATION.md`, Theorem 24; `verification/verify_frontier_walls.py` | The equivalence is an exact contraction argument. The countermodel is an abstract integer sequence, not a family of sign matrices. |
| The complete order-9 Bellman frontier is \(\mathcal B_9=\{(12,0)\}\), with 20/35 root representatives and 4/11 switching classes by deficit. | computational only | `verification/research_order9_weighted_geometry.py` | All 12,346 root-normalized records are evaluated with exact integer arrays. Completeness trusts nauty count/digest; class grouping trusts NetworkX isomorphism. |
| The records `GHOgmo` and `Gxd?Dc` have identical energy histograms and exact-maximizer pair-distance laws but extension values 15 and 17. | computational only | same order-9 verifier | Both records are directly decoded; energies, radii, pair enumerators, and incident-sign extensions are recomputed exactly. |
| The complete energy-coloured two-point law separates every weighted deficit at order 9. | computational only | same order-9 verifier | Exact finite catalogue statement only; no sufficiency is claimed at other orders or asymptotically. |
| The order-10 minimax partition function has exactly three positive-temperature histogram phases with thresholds \(0.658478948\ldots\) and \(0.792460762\ldots\). | computational only | `verification/research_order10_temperature.py`; `paper/second_attempt.tex` | All 274,668 rooted records and 6,012 histograms are exhausted in exact polynomial arithmetic. Completeness trusts the asserted nauty count and stream digest. |
| The order-10 optimum-extension catalogue alone implies \(F(11)\ge15\). | computational only | `RESEARCH_CONTINUATION.md`, exact orders 11 and 12 | Conditional on the complete nauty catalogue, a hypothetical order-11 maximum 13 would make every order-10 principal submatrix optimal, contradicting the extension minima 17 and 19; parity gives 15. |
| \(F(11)=17\) and \(F(12)=18\). | computational only | `verification/research_order11_certify.py`; `paper/second_attempt.tex` | The full nauty stream of 12,005,168 residual classes is hashed and exactly filtered to 2,153,606 analytically eligible classes; exact cut evaluation finds no maximum-15 survivor. Direct generation of the reduced stream reproduces the same result with the same decoder and evaluator, while deterministic samples at each eligible edge count are recomputed by a separate adjacency formula. Completeness trusts nauty. Both upper witnesses are evaluated over their full projective spin cubes by two formulas; order 12 then uses puncturing and parity. No solver status is used. |
| The annealed-normalized negative moment satisfies \(\mathcal G_{n+k}(q,t)\ge\mathcal G_n(q,t)+\mathcal G_k(q,t)\). | proved | `RESEARCH_CONTINUATION.md`, Theorem 25; `paper/second_attempt.tex` | Jensen over the cross disorder and the two global orientations gives the exact factors, which cancel under annealed normalization. The absolute value is retained by the augmented spin. |
| Reverse hypercontractivity transports \(\mathcal G_n/q\) along \((1+q)\tanh^2t=\theta\). | proved | same locations; Borell citation in `LITERATURE.md` | The parameter curve and inequality direction are exact. The direction is opposite to the missing lower comparison. Require positive exponents and sufficiently large orders in the extensive parametrization. |
| The power-saving transport estimate (PT). | flawed | RESEARCH_CONTINUATION.md, Theorem 27; verification/verify_negative_replica_transport_obstruction.py | PT is false. The finite-disorder Laplace bound and the audited limsup constant \(1/2\) give \(\liminf \Delta_n/n^2\ge\sqrt2-2\log2>0\) at \((\beta,\theta)=(4,8)\). The old finite data were only suggestive; the new analytic estimate is the disproof. |
| The reverse-hypercontractive transport defect has an exact entropy-production integral and satisfies \(R_q(f)\ge 2D(\mu\Vert\nu)^2/(3qm)\). | proved | RESEARCH_CONTINUATION.md, Theorem 28; paper/second_attempt.tex; verification/verify_negative_replica_transport_obstruction.py | Direct differentiation gives the identity. The remainder follows from a one-coordinate inequality and entropy tensorization. Equality holds only for constant \(f\). The finite grids are checks, not the proof. |
| The fine negative moment factors exactly into the two graph marginals, the rectangular marginal, and a nonnegative conditional relative-switching alignment term. | proved | RESEARCH_CONTINUATION.md, Theorem 29; paper/second_attempt.tex; verification/verify_negative_replica_alignment.py | The code dimensions, quotient inclusion, state-sum normalization, and fiber average were independently rederived. The factor \(2^{N+1}\) is correct for the stated augmented state-sum convention. |
| The alignment term has zero-temperature slope \(F(n)+F(k)+B_\square(n,k)-F(n+k)\). | proved | same locations | This is a sequential finite-space \(q\to\infty\), then \(t\to\infty\) statement. It does not justify exchanging either limit with \(n\to\infty\). The balanced normalized liminf is at least \(2/\pi+\sqrt{2/\pi}-\sqrt2>0\), so alignment is leading order. |
| The first conditional alignment Hamiltonian is the displayed mixed trace \(\mathcal H_4\), with the exact Walsh decomposition and variance formula. | proved | RESEARCH_CONTINUATION.md, Theorem 30; verification/verify_negative_replica_alignment.py | Direct block multiplication and orthogonality of the relative-switching characters prove the formulas. Higher Eulerian terms are not uniformly negligible at \(u\asymp N^{-1/2}\), \(q\asymp N\). |
| Identical local graph and rectangular scalar partition curves determine conditional alignment. | flawed | RESEARCH_CONTINUATION.md, exact \(2+4\) collision; same verifier | Two exact rational examples agree in all three complete local absolute-energy histograms but have \(K_2=1.471325\ldots\) and \(1.135146\ldots\). Their \(\mathcal H_4\) laws already differ. Since \(\dim D_2=1\), the displayed relative gauges cover the true fiber twice; the verifier checks this redundancy and the normalized moments are unchanged. |
| The proposed alignment-transport inequality (AT) is a sharper independent settling lemma. | redundant | RESEARCH_CONTINUATION.md, after Theorem 30 | After defining \(a_j=\mathcal G_j(q_j,t_j)/q_j\), its left side is identically \(a_N-a_n-a_k\). AT is a valid sufficient statement only because it restates the desired almost-superadditivity; it is not yet a derived or weaker estimate. |
| For odd prime-power \(m\), the Paley conference matrix of order \(m^2+1\) has a Boolean eigenvector and maximum \(m(m^2+1)/2\). | proved | `paper/second_attempt.tex`, regular Paley theorem; `verification/verify_paley_subfield.py` | The additive-subfield coset proof is complete, but the result is known regular-conference theory, not new. It concerns that matrix, not the minimax \(F(m^2+1)\). |
| The order-11 normalized value is a record high for \(n\ge5\). | flawed | external report audited after revision `9946030` | \(F(11)/11^{3/2}=0.46597\ldots\), while \(F(7)/7^{3/2}=0.48595\ldots\). |
| \(F(13)=20\) and \(F(14)=21\) from a complete billion-record residual scan. | incomplete | external report received 2026-08-02; local artifact search; verification/check_conference_examples.py | The deduction is valid conditional on the reported complete order-12 optimum catalogue and extension minima. The shard logs, scanner, stream digests, survivor encodings, and extension/deletion receipts are absent from this checkout and accessible Git refs, so the lower bound cannot be replayed. Independently verified here: every order-13 principal submatrix of \(C_{14}\) has maximum 20 and \(M(C_{14})=21\), proving only \(F(13)\le20,\ F(14)\le21\). |
| Order 14 would be the first conference order proved exactly optimal. | flawed | same external report; verification/check_conference_examples.py | Already \(F(6)=5=M(C_6)\). Order 14 would be another exact conference order if the missing lower certificate replays. |
| Two optimal classes among \(1{,}018{,}997{,}864\) imply probability about \(2\times10^{-9}\) under random signings. | flawed | same external report | The denominator counts root-normalized unlabeled residual representatives, not uniformly labelled signings. A probability statement requires automorphism-weighted multiplicities. |

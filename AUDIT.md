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
| \(F(13)=20\) and \(F(14)=21\). | computational only | `verification/research_order13_certify.py`; `verification/order12_threshold_scan.c`; `verification/check_conference_examples.py` | Eight fail-closed nauty shards hash and scan all 1,018,997,864 residual graphs. Exactly two rooted records have order-12 maximum below 20; independent full-cube and incident-column enumeration gives \((M,E)=(18,24)\) for both. Every other predecessor has \(M\ge20\). The Bellman identity and a \(C_{14}\) principal witness give \(F(13)=20\); heredity, parity, and \(M(C_{14})=21\) give \(F(14)=21\). Completeness trusts nauty. |
| \(F(15)\in\{25,27\}\) from the near-minimal layer tower. | computational only | quarantined 2026-08-02 layer replay; external certificate bank | A fresh reconstruction matched the complete order-14 \(M=21\) and \(M=23\) layers and their extension values 27 and 29, giving the bracket. Integration remains pending: the supplied public verifier does not execute this tower, so the result is not yet included in the main exact-value table. Exact \(F(15)\) remains open because the order-13 predecessors with \((M,\delta_{\rm w})=(24,0)\) have not been classified. An external sweep is ongoing, but no partial count or eventual verdict is promoted without its banked v2 certificate. |
| Order 14 would be the first conference order proved exactly optimal. | flawed | same external report; verification/check_conference_examples.py | Already \(F(6)=5=M(C_6)\). Order 14 is another exact conference order; it is not the first. |
| Two optimal classes among \(1{,}018{,}997{,}864\) imply probability about \(2\times10^{-9}\) under random signings. | flawed | same external report | The denominator counts root-normalized unlabeled residual representatives, not uniformly labelled signings. A probability statement requires automorphism-weighted multiplicities. |
| \(M(C_{30})=75\), with 812 projective maximizers. | computational only | `verification/research_paley_alignment.c` | A source-built Gray-code scan exhausts all \(2^{29}\) projective Boolean states after verifying the Paley conference identity. The exact result was separately reproduced by a meet-in-the-middle evaluator during audit. It implies only \(F(30)\le75\). |
| The Paley alignment ratios for \(q=5,13,17,29\) increase and therefore converge to 1. | flawed | `verification/research_paley_alignment.c`; `RESEARCH_CONTINUATION.md`, Theorem 31 | The four exact ratios are increasing, but four selected prime-field cases prove neither monotonicity nor an asymptotic decay law. Square-\(q\) ratios are identically 1, so no single smooth fit covers all conference orders. |
| Paley spectral saturation is equivalent to asymptotic balance and vanishing minority-half Fourier leakage. | proved | `RESEARCH_CONTINUATION.md`, Theorem 31; `paper/second_attempt.tex`, Paley Fourier-leakage proposition | Additive Fourier diagonalization of the Paley core and exact optimization of the infinity sign give the conditions \(|S(f_q)|=o(q)\) and \(W(f_q)=o(q)\). This concerns \(M(C_{q+1})\), not the minimax \(F(q+1)\). On the proved multiplicatively dense sequence, the matching minimax lower bound would already be equivalent to the full limit \(1/2\), as shown below. |
| The relative-gauge state is the balanced pushforward convolution of the two graph Gibbs laws and the rectangular Gibbs law. | proved | `RESEARCH_CONTINUATION.md`, Theorem 32; `paper/second_attempt.tex`, microcanonical composition theorem | The state-space cardinalities, homomorphism, fiber size, and finite-temperature density normalization are explicit. Fourier coefficients factor into three local correlation tensors. |
| The \(2^{n+k-1}\)-st smallest product-triple deficit is a guaranteed block-composition gain. | proved | same theorem; `verification/verify_relative_profile_composition.py` | Exact max-plus identity plus one minimizing triple from each balanced fiber. The verifier exhausts all \(2+3\) block triples and includes a corrupted relative-orientation control. |
| The microcanonical order statistic captures the complete alignment gain. | flawed | exact \(2+4\) collision in the same verifier | Both examples have scalar order statistic zero but true gains 4 and 2. The order statistic is optimal only after discarding additive labels; the group geometry can improve it. |
| A half-interval sign function has Paley leakage at most \(2p/(\ell(p)-1)\), and asymptotically saturates the Paley spectral ceiling on a multiplicatively dense prime sequence. | proved | `RESEARCH_CONTINUATION.md`, Theorem 33; `paper/second_attempt.tex`, least-nonresidue theorem; `verification/verify_paley_least_nonresidue.py` | The Fourier formula and tail bound are elementary. Quadratic reciprocity and the prime number theorem in each fixed progression produce levels \(\ell(p)\to\infty\) with consecutive prime ratios tending to one. This remains a statement about \(M(C_{p+1})\), not a lower bound for \(F\). |
| A nonconstant Boolean function at prime order can have exact one-half Fourier support. | flawed | `paper/second_attempt.tex`, exact prime-order leakage obstruction | Vanishing at one nonzero frequency forces the coefficient polynomial to be \(\pm\Phi_p\), hence constant. The norm refinement gives \(W(f)\ge4(p-1)^2/p^3\). |
| The interval construction proves \(W=o(p)\) uniformly over all primes \(p\equiv1\pmod4\). | flawed | `RESEARCH_CONTINUATION.md`, after Theorem 33 | Along primes \(p\equiv5\pmod{12}\), the opposite-character frequency pairs \(\pm1\) and \(\pm3\) force \(\liminf W/p\ge8/(9\pi^2)\). The specially selected progressions are essential. |
| Heuristic witnesses at \(q=37,41,53,61\) disprove monotonicity of the exact Paley alignment ratio. | flawed | unbanked external search report audited 2026-08-02 | The displayed energies are certified lower bounds, not exact maxima. A lower best-found ratio at \(q=61\) cannot prove that the true maximum ratio falls. The witnesses may be useful finite evidence but no monotonicity claim follows. |
| The minimax problem is, up to an explicit additive \(O(n)\) term, the minimum two-sided cut deviation of a graph with \(\lfloor\binom n2/2\rfloor\) edges. | proved | proof below, “Fixed-half-density cut-deviation reduction”; `verification/verify_cut_discrepancy_equivalence.py` | With \(m=\binom n2\) and \(H_n\) as defined below, \(F(n)-1\le4H_n\le F(n)+\sqrt m+2\). Consequently \(F(n)/n^{3/2}\) converges iff \(4H_n/n^{3/2}\) converges, with the same limit. The script exhausts orders 2 through 6 and checks normalization and edge cases; the proof is analytic. |
| The Backurs--Bavarian cut-deviation theorem applies to the fixed-half-density reformulation and gives \(H_n=\Omega(n^{3/2})\). | proved | `LITERATURE.md`, “Cut deviation is an exact reformulation” | Their constant is unspecified and nonsharp. The theorem supplies neither the candidate constant \(1/8\) for \(H_n\) nor a relation between different orders, so it does not address existence of the normalized limit. |
| The Erdős--Goldberg--Pach--Spencer “bipartite discrepancy” is the same parameter as the all-cuts, two-sided cut deviation used here. | flawed | `LITERATURE.md`, same section | Their Section 3 parameter is a one-sided maximum over balanced complementary halves. Backurs--Bavarian revision 2, Definition 7.1, is the direct source for the absolute maximum over all cuts. |
| Projective \(|Q|\) graph states together with signed full-spin cross states satisfy the same exact max-plus identity as the augmented/projective relative-profile theorem. | flawed | `RESEARCH_CONTINUATION.md`, one-sided swapped-profile theorem; `verification/verify_swapped_profile_injection.py` | The alternate map has equal fibers and its fiber maximum dominates \(M(Y_g)\), which is enough for a valid order-statistic upper bound, but equality can fail. On the standard \(C_{14}\) split it gives \(\widetilde\Lambda=8\) and 596440 target states, versus \(\Lambda=10\) and 304908 for the exact theorem; the true gain is 22. |
| The one-sided swapped raw profile can prove balanced near-subadditivity independently of the eventual minimax constant. | flawed | `RESEARCH_CONTINUATION.md`, equations (41.5)--(41.7); `verification/verify_swapped_profile_floor.py` | Its certified energy is always at least \(F(r)+r\mu_r-O(r)\), so closure through this statistic requires \(\liminf F(r)/r^{3/2}\ge0.436377\ldots\). The obstruction is compatible with a limit \(1/2\). It transfers to the canonical exponential relaxation, but not to the raw exact augmented statistic. |
| The optimized canonical exponential certificates of the exact and one-sided profile conventions differ by only \(O(\sqrt N)\) at the mean-field scale. | proved | `RESEARCH_CONTINUATION.md`, equation (41.7); `paper/second_attempt.tex`, swapped-profile canonical comparison | Their product moment generating functions have pointwise ratio in \([1/2,4]\). Approximate minimizers and the entropy term \((N-1)\log2/t\) give losses \(L/(N-1)\) and \(2L/(N-1)\), respectively. This comparison does not identify the two raw order statistics. |
| Nontrivial Fourier mass in a labeled subthreshold shell can improve the scalar microcanonical gain. | proved | `RESEARCH_CONTINUATION.md`, Theorem 42; `verification/verify_labeled_shell_parseval.py` | If \(V_s\) is the nonconstant Fourier mass of the fiber occupancy, then \(\min b_s\le\mu_s-\sqrt{V_s/(2^{N-1}-1)}\). The exact \(2+4\) collision gains two units beyond the scalar order statistic. At the sampled \(C_{14}\) target the generic variance bound is too weak, so no asymptotic claim follows. |
| A finite integer-moment hierarchy can certify an empty labeled shell without explicitly minimizing its fiber occupancies. | proved | `RESEARCH_CONTINUATION.md`, Theorem 42A; `verification/verify_labeled_shell_moment_certificate.py` | Consecutive-root polynomials give lower bounds for the zero mass, and the equivalent localizing matrices are complete at every finite order. A degree-19 exact certificate detects an empty balanced \(C_{14}\) target fiber; direct reconstruction separately proves it is unique. This is finite closure only. |
| A polynomial-accuracy estimate of normalized localizing moments can robustly detect an isolated good gauge. | flawed | `RESEARCH_CONTINUATION.md`, equations (42.9)--(42.11); `verification/verify_labeled_vacancy_hierarchy.py` | With \(z\) holes every normalized negative localizing form has magnitude at most \(z/K\), and filling the holes by occupancy one changes every localizing matrix by exactly \((z/K)e_0e_0^{\mathsf T}\). The unique \(C_{14}\) hole has signal below \(1/8192\). An exact sign, exponential precision, or an abundance theorem is required. |
| A subexponential number of generic occupancy moments or mixed-cycle layers universally decides relative-gauge vacancy. | flawed | `RESEARCH_CONTINUATION.md`, Theorem 42A; `verification/verify_universal_moment_obstruction.py` | A restricted-coefficient polynomial collision produces a vacant and nonvacant \(K\)-point occupancy pair with matching moments through degree \(\Omega(\sqrt{K/\log K})\). For \(K=2^{N-1}\) this is exponential in \(N\). The pair is abstract and is not claimed realizable as signing shells. |
| A sparse Fourier--PSD minor can detect the unique empty fiber in the calibrated \(C_{14}\) shell. | flawed | same theorem; `verification/verify_labeled_vacancy_hierarchy.py` | If there are \(z\) holes and every positive occupancy is at least \(m\), all principal Fourier kernels on at most \((m-1)K/(mz)\) characters are positive semidefinite. At \(C_{14}\), any negative minor needs at least 6827 of 8192 characters. |
| Fixed occupancy inverse temperature or a bounded collision truncation is enough to convert the canonical shell law to vacancy. | flawed | same theorem and verifier | Uniformly, the soft minimum needs \(t\ge\log K=\Theta(N)\). At \(C_{14}\), the exact alternating factorial-moment expansion first certifies the hole at degree 87. This closes fixed-temperature and fixed-degree conversions, not exact all-order convolution. |
| Imposing exact half density can make the separately controlled rectangular cross discrepancy subleading. | flawed | `RESEARCH_CONTINUATION.md`, Proposition 43; `verification/verify_fixed_density_cross_floor.py` | Fixed total changes the optimum rectangular \(\infty\to1\) norm by at most \(O(\sqrt{nk})\), but every sign matrix still has norm at least \(\max(n\mu_k,k\mu_n)\). This rules out separate triangle-inequality control, not joint cancellation. |
| Complete/empty cloud amplification can be made asymptotically lossless by Seidel switching and \(O(N)\) fixed-density repair. | flawed | `RESEARCH_CONTINUATION.md`, Theorem 44; `verification/verify_equal_cloud_blowup.py` | Cloud-union cuts retain \(k^2d(G)\), while switching plus \(r\) repairs loses at most \((r+1)/2\). The normalized bound worsens by \(\sqrt k\). This theorem does not cover orthogonal cloud blocks. |
| The common symmetric-Hadamard four-fold lift of the order-five optimum reaches the lossless value 32 after a fixed-half diagonal completion. | computational only | `verification/verify_hadamard_cloud_lift.c` | Classification of all 64 symmetric order-four Hadamards, an explicit signed-five-cycle anti-isomorphism for global sign, and exhaustion of both representatives, all transformed diagonal completions, and every projective order-20 spin give sharp minimum 44. This is a no-go only for the common-\(H\), common-\(D\) family. |
| Orthogonal cloud blocks are locally incompatible with the lossless factor-four scale. | flawed | `RESEARCH_CONTINUATION.md`, order-16 construction; `verification/verify_cloud_dependent_hadamard_lift.py` | An explicit fixed-half order-16 signing with six zero-sum Hadamard cross blocks has maximum 32 and cross-only maximum 28. After reordering and signed-permutation gauges its edge signs form an order-four optimizer and all cross blocks share one nonsymmetric oriented Hadamard. Thus it is an exact one-step lift at \(4^{3/2}F(4)\), but its four internal blocks are tailored and no uniform or iterable operator is proved. |
| The common oriented-Hadamard order-16 frame can be completed to maximum 28 by choosing suitable internal signings. | flawed | `RESEARCH_CONTINUATION.md`, Theorem 45; `verification/verify_framed_hadamard_lift_30.py` | Six explicit cross-energy states force a three-cycle of inequalities among one block's framed response values. Maximum at most 28 would make the three values equal, but two differ by twice a sum of three signs. Thus the restricted minimum is 30. The explicit alternating completion has maximum 30, proving \(F(16)\le30\), but not equality. |
| \(F(16)\le30\), and 30 is the exact minimum within the pinned common oriented-Hadamard cross frame. | proved | `RESEARCH_CONTINUATION.md`, Theorem 45; `verification/verify_framed_hadamard_lift_30.py`; `verification/verify_framed_hadamard_lift_30.c` | The upper witness is an admissible sign matrix whose full histogram is independently reconstructed in Python and strict C. The restricted lower bound is the six-state pencil proof above. This does not provide a universal lower bound \(F(16)\ge30\). |
| For a symmetric conference matrix, the Boolean quadratic maximum is exactly controlled by the largest \(\ell_1\)-norm of a unit vector in either conference eigenspace. | proved | proof below, “Conference eigenspace formula” | If \(\alpha_\pm=n^{-1/2}\sup_{v\in E_\pm,\|v\|_2=1}\|v\|_1\), then \(2M(C)/(n\sqrt{n-1})=\max\{2\alpha_+^2-1,2\alpha_-^2-1\}\). Both signs are necessary because \(M\) contains an absolute value. |
| Dense-order Paley minimax rigidity is a strictly weaker intermediate lemma toward convergence. | flawed | proof below, “Dense Paley rigidity is full closure” | On the multiplicatively dense sequence where the Paley maximum is \((1/2+o(1))n^{3/2}\), the lower bound \(F(n)\ge(1-o(1))M(C_n)\), together with monotonicity, already implies \(\lim F(n)/n^{3/2}=1/2\). Conversely that limit implies the rigidity ratio. It is an equivalent full-strength target, not an intermediate lemma. |

## Fixed-half-density cut-deviation reduction

Put \(m=\binom n2\), \(e_0=\lfloor m/2\rfloor\), and

\[
 H_n=\min_{\substack{G\text{ on }[n]\\ e(G)=e_0}}
 \max_{S\subseteq[n]}
 \left|e_G(S,S^c)-\frac12|S||S^c|\right|.
\]

Given a signing \(A\), let \(G\) contain exactly the pairs with \(a_{ij}=-1\),
and write \(t=\sum_{i<j}a_{ij}=m-2e(G)\).  If \(x_S\) is \(+1\) on
\(S\) and \(-1\) on \(S^c\), and \(k=|S||S^c|\), direct expansion gives

\[
 Q_A(x_S)=t-2k+4e_G(S,S^c).                                      \tag{A.1}
\]

First restrict to signings with total sum
\(t_0=m-2e_0\in\{0,1\}\), and call their minimum maximum \(F_0(n)\).
Equation (A.1) gives

\[
 4\left(e_G(S,S^c)-\frac{k}{2}\right)=Q_A(x_S)-t_0.
\]

The vectors \(x_S\) exhaust the Boolean cube modulo global negation, so
the elementary inequality
\(\big|\max_x|u_x-c|-\max_x|u_x|\big|\le |c|\) yields

\[
 F_0(n)-1\le 4H_n\le F_0(n)+1.                                  \tag{A.2}
\]

It remains to compare \(F_0\) with \(F\).  Switching \(A\) by a sign vector
\(z\) replaces its total sum by \(Q_A(z)\).  For uniform \(z\), Walsh
orthogonality gives

\[
 \mathbb E Q_A(z)=0,
 \qquad
 \mathbb E Q_A(z)^2=m.
\]

Thus some switch has total sum \(t\) with \(|t|\le\sqrt m\).  Since
\(t\equiv t_0\pmod2\), changing
\(r=|t-t_0|/2\le(\sqrt m+1)/2\) edge signs reaches total sum \(t_0\).
If \(t>t_0\) flip \(r\) positive coefficients, and if \(t<t_0\) flip \(r\)
negative coefficients; the corresponding sign count is always at least
\(r\).
Each changed edge changes every value of \(Q_A\) by at most \(2\), hence
changes \(M(A)\) by at most \(2r\).  Starting from an optimizer for \(F(n)\)
therefore proves

\[
 F(n)\le F_0(n)\le F(n)+\sqrt m+1.                               \tag{A.3}
\]

Combining (A.2) and (A.3) gives the explicit one-sided finite bound

\[
 \boxed{F(n)-1\le4H_n\le F(n)+\sqrt{\binom n2}+2}.                \tag{A.4}
\]

The error is \(O(n)=o(n^{3/2})\), proving the normalized equivalence.

## Conference eigenspace formula

Let \(C=C^{\mathsf T}\) be a conference matrix of order \(n\), so
\(C^2=(n-1)I\), and put \(r=\sqrt{n-1}\).  Since \(\operatorname{tr}C=0\),
the eigenspaces \(E_\pm=\ker(C\mp rI)\) both have dimension \(n/2\), and

\[
 P_\pm=\frac12(I\pm C/r)
\]

are their orthogonal projections.  For any orthogonal projection \(P\) onto
\(E\),

\[
 \max_{x\in\{\pm1\}^n}x^{\mathsf T}Px
 =\left(\sup_{\substack{v\in E\\\|v\|_2=1}}\|v\|_1\right)^2.    \tag{A.5}
\]

Indeed, \((x^{\mathsf T}Px)^{1/2}=\|Px\|_2
=\sup_{v\in E,\|v\|_2=1}\langle v,x\rangle\); maximize over \(x\),
interchange the two suprema, and use
\(\max_{x\in\{\pm1\}^n}\langle v,x\rangle=\|v\|_1\).

Define

\[
 \alpha_\pm=\frac1{\sqrt n}
 \sup_{\substack{v\in E_\pm\\\|v\|_2=1}}\|v\|_1.
\]

The average of \(x^{\mathsf T}P_\pm x\) over the cube is
\(\operatorname{tr}P_\pm=n/2\), so \(\alpha_\pm^2\ge1/2\).  Therefore
no positive-part convention is hidden in the following calculation:

\[
 \max_x x^{\mathsf T}Cx=nr(2\alpha_+^2-1),
 \qquad
 \max_x(-x^{\mathsf T}Cx)=nr(2\alpha_-^2-1).
\]

Since \(Q_C(x)=x^{\mathsf T}Cx/2\),

\[
 \boxed{\frac{2M(C)}{n\sqrt{n-1}}
 =\max\{2\alpha_+^2-1,\,2\alpha_-^2-1\}}.                       \tag{A.6}
\]

Moreover, \(\alpha_\pm=1\) exactly when the corresponding eigenspace
contains a normalized sign vector.  Quantitatively, for a unit vector \(v\)
and \(s_i=\operatorname{sgn}(v_i)\), with either sign chosen when \(v_i=0\),

\[
 \left\|v-\frac{s}{\sqrt n}\right\|_2^2
 =2\left(1-\frac{\|v\|_1}{\sqrt n}\right).
\]

Thus conference spectral saturation is precisely asymptotic proximity of at
least one eigenspace to a cube direction.  If instead both
\(\alpha_\pm^2\le1-\eta\), then
\(M(C)\le(1/2-\eta)n\sqrt{n-1}\); this is the exact geometric obstruction
to a universal conference-alignment assertion.

## Dense Paley rigidity is full closure

Let \(n_j\) be the multiplicatively dense Paley orders supplied by the
least-nonresidue construction, so \(n_{j+1}/n_j\to1\) and
\(M(C_{n_j})=(1/2+o(1))n_j^{3/2}\).  If
\(F(n_j)\ge(1-o(1))M(C_{n_j})\), then for
\(n_j\le n<n_{j+1}\), monotonicity gives

\[
 \frac{F(n)}{n^{3/2}}
 \ge \frac{F(n_j)}{n_j^{3/2}}\left(\frac{n_j}{n}\right)^{3/2}
 =\frac12-o(1).
\]

Together with the audited limsup bound \(1/2\), this proves the full limit
and its value.  Conversely, if the full limit equals \(1/2\), then along this
sequence \(F(n_j)/M(C_{n_j})\to1\).  The dense-order rigidity statement is
therefore equivalent to full closure at constant \(1/2\), conditional only on
the already proved dense Paley alignment theorem.

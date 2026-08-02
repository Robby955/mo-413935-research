# Literature and concept map

## Scope and notation

This is a targeted map of results that bear directly on MathOverflow Question
413935.  It is not an exhaustive bibliography.  Every external reference below
is to a primary paper (or its official journal/arXiv page), and each entry
separates what the cited theorem supplies from what it does **not** supply.

Write

\[
 M(A)=\max_{x\in\{\pm1\}^n}\left|\sum_{i<j}a_{ij}x_i x_j\right|,
 \qquad
 F(n)=\min_{A}M(A),
\]

where \(A\) is symmetric, has zero diagonal, and has off-diagonal entries in
\(\{\pm1\}\).  The unresolved issue is convergence of
\(F(n)/n^{3/2}\), not merely its boundedness.

## Executive map

| Area | Closest usable fact | Exact gap to the MO problem |
|---|---|---|
| Mean-field spin glasses | Random-disorder free energies and ground states have thermodynamic limits, often with Parisi variational formulas. | The proofs average over independent disorder.  Here the disorder is chosen adversarially by an outer minimum. |
| Reverse hypercontractivity | Negative norms of positive Boolean functions compare sharply along the noise curve \(\rho^2=(1-p)/(1-q)\). | Applied to the disorder partition polynomial, the inequality has the correct parameter curve but the wrong direction for the needed block lower bound. |
| Spin-glass large deviations | A 2026 preprint proves an upper-tail LDP for the maximum of a Gaussian mixed \(p\)-spin Hamiltonian. | Minimizing over signings probes a lower tail of a Bernoulli two-spin ground state (and an absolute-value version), precisely the side not covered. |
| Discrepancy and Boolean polynomials | General results recover the \(n^{3/2}\) scale for degree-two sign polynomials or rectangular switching games. | They give constants or a different bilinear state space, not an \(o(n^{3/2})\) composition/rounding error. |
| Signed graphs and cut codes | Switching classes are cosets of the cut/cocycle code, and frustration is a coset-leader problem. | The code here is the cut code augmented by the all-one word; no cited result gives its complete-graph covering-radius deficit to second order. |
| Grothendieck/SDP | Boolean bilinear or quadratic optima have constant-factor semidefinite relaxations. | A fixed constant factor is leading-order at scale \(n^{3/2}\); convergence needs asymptotically lossless control. |
| Dense graph limits | Cut-metric convergence controls normalized dense energies/free energies at scale \(n^2\). | Near-optimal signings are fluctuation objects at scale \(n^{3/2}\); ordinary graphons erase that scale. |
| Seidel/conference matrices | Paley gives an infinite, dense family of symmetric conference matrices; Seidel switching and spectra are well developed. | A spectral ceiling need not be attained by a Boolean vector, and construction upper bounds cannot prove a matching lower bound for \(F(n)\). |
| Vector balancing | Powerful theorems choose signs so a vector sum lies in a prescribed convex body. | Applied to the exponentially many spin constraints, their generic normalization/constant losses do not yield a lossless amplification theorem. |

## 1. Mean-field spin glasses and thermodynamic limits

### Random-disorder interpolation

Guerra and Toninelli prove existence of the thermodynamic limit for the
Sherrington--Kirkpatrick and related mean-field models by smoothly
interpolating between one random system of size \(N\) and two **independent**
random subsystems; after taking the quenched expectation, the free energy is
subadditive.  They also obtain the analogous ground-state limit.

* Source: F. Guerra and F. L. Toninelli, [*The Thermodynamic Limit in Mean
  Field Spin Glass Models*](https://doi.org/10.1007/s00220-002-0699-y),
  *Communications in Mathematical Physics* 230 (2002), 71--79; open preprint
  [arXiv:cond-mat/0204280](https://arxiv.org/abs/cond-mat/0204280).
* Gives here: the correct model for what a successful finite-temperature
  interpolation would look like.
* Does not give: an inequality preserved by
  \(\min_{A\in\{\pm1\}^{\binom n2}}\).  Independence and disorder expectation
  are used in the interpolation; there is no quenched expectation after the
  adversary has selected the least energetic realization.

Auffinger and Chen identify the thermodynamic-limit ground-state energy of the
Gaussian mixed \(p\)-spin model by a zero-temperature Parisi variational
formula.

* Source: A. Auffinger and W.-K. Chen, [*Parisi formula for the ground state
  energy in the mixed p-spin model*](https://doi.org/10.1214/16-AOP1173),
  *Annals of Probability* 45 (2017), 4617--4631; [arXiv:1606.05335](https://arxiv.org/abs/1606.05335).
* Gives here: a rigorous finite-temperature-to-ground-state paradigm for a
  Gaussian Hamiltonian with a prescribed covariance as a function of spin
  overlap.
* Does not give: a variational formula after minimizing over deterministic
  coupling arrays.  Its Hamiltonian is a Gaussian process; the proof does not
  optimize its sample path.

Carmona and Hu prove universality of the limiting SK free energy with respect
to the independent environment distribution (under their moment assumptions).

* Source: P. Carmona and Y. Hu, [*Universality in Sherrington--Kirkpatrick's
  spin glass model*](https://doi.org/10.1016/j.anihpb.2005.04.001), *Annales
  de l'Institut Henri Poincare Probabilites et Statistiques* 42 (2006),
  215--222; [arXiv:math/0403359](https://arxiv.org/abs/math/0403359).
* Gives here: a warning that typical Gaussian and typical Bernoulli disorder
  can share a limiting free energy.
* Does not give: universality for an extreme disorder realization.  Replacing
  an expectation over i.i.d. signs by a minimum over all sign arrays is not a
  universality step.

### The relevant large-deviation side is still missing

Chen, Guionnet, Ko, Lacroix-A-Chez-Toine, and Mourrat consider a centered
Gaussian mixed \(p\)-spin field with covariance
\(N\xi(\sigma\cdot\tau/N)\), plus an external field.  Their Theorem 1.2 proves,
for the normalized maximum \(L_N\), an LDP for events \(L_N\ge r\) with
\(r\) at or above the typical ground-state value, at speed \(N\).  The paper
explicitly stresses that it treats **upper** deviations.  It says that, without
external field, lower deviations are expected to have speed \(N^2\), and cites
a proof only for spherical models.

* Source: H.-B. Chen, A. Guionnet, J. Ko, B. Lacroix-A-Chez-Toine, and
  J.-C. Mourrat, [*One-sided large deviations for the ground-state energy of
  spin glasses*](https://arxiv.org/abs/2603.06368), arXiv:2603.06368v1 (2026),
  especially Theorem 1.2 and the discussion immediately after Theorem 1.3.
* Gives here: a rigorous fractional-moment/convex-duality route for the upper
  tail of a random Gaussian ground state.
* Does not give: the lower-tail rate function, Bernoulli disorder, or the
  minimum over all \(2^{\binom n2}\) disorders.  The MO minimization is naturally
  an extreme lower-tail question, and the absolute value simultaneously
  constrains the maxima of \(Q_A\) and \(-Q_A\).  Thus this new preprint does not
  settle, or directly interpolate, the present problem.

This is the most relevant current spin-glass obstruction: a lower-tail LDP at
speed \(n^2\), with enough control to identify the endpoint reached among
exponentially many Bernoulli environments, would be genuinely new input.

### Negative moments and reverse hypercontractivity

Borell proved the sharp reverse-hypercontractive inequality for the two-point
space, which tensorizes to the Boolean cube.  In negative-exponent notation,
for a positive function and \(Q'>Q>0\),

\[
 \|T_\rho f\|_{-Q'}\ge\|f\|_{-Q}
 \quad\text{when}\quad
 \rho\le\sqrt{\frac{1+Q}{1+Q'}}.
\]

* Source: C. Borell, [*Positivity improving operators and
  hypercontractivity*](https://doi.org/10.1007/BF01318906), *Mathematische
  Zeitschrift* 180 (1982), 225--234.  A later product-space treatment is E.
  Mossel, K. Oleszkiewicz, and A. Sen, [*On reverse
  hypercontractivity*](https://arxiv.org/abs/1108.1210).
* Gives here: an exact comparison for the annealed-normalized disorder
  polynomial because \(f_{n,u}=T_{u/v}f_{n,v}\).  Its invariant curve is
  exactly \((1+q)u^2=\theta\), the curve generated by extensive block
  composition.
* Does not give: the lower comparison needed to combine the exact
  negative-replica supermultiplicativity across different sizes.  It proves
  only that the transport defect is nonnegative.  A power-saving upper bound
  on that defect is a new, strong conjecture, and finite data through order 9
  do not show saturation.

## 2. Discrepancy and Boolean quadratic forms

### Rectangular switching is close but not the same state space

Brown and Spencer study minimization for rectangular sign matrices under
independent row and column switches, the prototype of

\[
 \min_{C\in\{\pm1\}^{r\times s}}
 \max_{u\in\{\pm1\}^r,\,v\in\{\pm1\}^s}|u^{\mathsf T}Cv|.
\]

* Source: T. A. Brown and J. H. Spencer, [*Minimization of \(\pm1\) matrices
  under line shifts*](https://doi.org/10.4064/cm-23-1-165-171), *Colloquium
  Mathematicum* 23 (1971), 165--171 (with an erratum on p. 177).
* Gives here: the natural discrepancy scale for a complete bipartite
  cross-block, and a classical source for the Gale--Berlekamp connection.
* Does not give: the symmetric quadratic problem, where the two switching
  vectors are constrained to be the same.  In block composition the
  rectangular term is a new leading-order \(n^{3/2}\) cost, not a negligible
  remainder.

The exact cross-seed orbit in Proposition 23 is a weighted version of this
same switching game.  Its cross profile
\(c(xy^{\mathsf T})=|x^{\mathsf T}C_0y|\) is precisely the
Gale--Berlekamp orbit profile, but the internal energy
\(h(xy^{\mathsf T})\) turns the composition problem into the max-plus
convolution \(\min_g\max_R[h(R)+c(gR)]\).  Brown--Spencer controls the
cross-only term; it does not supply the weighted noncoverage theorem needed
here.

Modern work on Bennett's inequality gives constructive constant bounds for
sign bilinear forms and improves Gale--Berlekamp constants.

* Source: D. Pellegrino and A. Raposo Jr., [*Upper bounds for the constants of
  Bennett's inequality and the Gale--Berlekamp switching game*](https://doi.org/10.1112/mtk.12229),
  *Mathematika* 70 (2024); preprint [arXiv:2111.00445](https://arxiv.org/abs/2111.00445).
* Gives here: dimensionally sharp exponents for rectangular bilinear forms.
* Does not give: a limiting constant for symmetric zero-diagonal forms, nor an
  \(o((n+m)^{3/2})\) cross-block estimate.

### Fourier--Walsh inequalities explain the scale, not convergence

Defant, Mastylo, and Perez prove a dimension-free Boolean
Bohnenblust--Hille inequality: for a degree-\(d\) Boolean function bounded by
one, the \(\ell_{2d/(d+1)}\)-norm of its Fourier coefficients is bounded in
terms of \(d\) only.  For \(d=2\), applying this to a quadratic polynomial with
\(\Theta(n^2)\) unit coefficients forces a supremum of order \(n^{3/2}\).

* Source: A. Defant, M. Mastylo, and A. Perez, [*On the Fourier spectrum of
  functions on Boolean cubes*](https://doi.org/10.1007/s00208-018-1756-y),
  *Mathematische Annalen* 374 (2019), 653--680; [arXiv:1706.03670](https://arxiv.org/abs/1706.03670).
* Gives here: a general functional-analytic explanation of the correct
  exponent, independent of the Gaussian-sign proof in this repository.
* Does not give: a sharp degree-two constant for this restricted support, an
  optimizer classification, or any relation between different values of
  \(n\).  It therefore cannot by itself force convergence.

## 3. Cut codes, covering radii, and signed graphs

For an ordinary graph \(G\), the binary cocycle (cutset) code is the span of
its cuts.  Sole and Zaslavsky identify cosets of this code with switching
classes of signings and relate coset leaders to the line index of imbalance
(frustration).

* Source: P. Sole and T. Zaslavsky, [*A Coding Approach to Signed Graphs*](https://doi.org/10.1137/S0895480189174374),
  *SIAM Journal on Discrete Mathematics* 7 (1994), 544--553.
* Gives here: the exact signed-graph/coset language underlying the repository's
  covering-radius identity.
* Does not give: that exact identity without one modification.  Because the MO
  objective has an absolute value, the relevant code is the \(K_n\) cut code
  **augmented by the all-one word**, thereby identifying a signing with its
  global negation.  The 1994 paper treats the ordinary cocycle code and does
  not derive the order-\(n^{3/2}\) deficit in its covering radius.

The same authors determine covering radii of **cycle** codes through minimum
\(T\)-joins for several graph classes.

* Source: P. Sole and T. Zaslavsky, [*The covering radius of the cycle code of
  a graph*](https://doi.org/10.1016/0166-218X%2893%2990140-J), *Discrete Applied
  Mathematics* 45 (1993), 63--70.
* Gives here: evidence that graphical code structure can make a covering
  radius tractable.
* Does not give: the covering radius of the augmented cocycle code of \(K_n\).
  Passing to a dual code does not turn covering radius into minimum distance;
  the required second-order coset information remains.

Bowlin obtains an upper bound and extremal-family characterization for maximum
frustration of \(K_{l,r}\), including exact formulas for \(l\le7\).

* Source: G. S. Bowlin, [*Maximum Frustration in Bipartite Signed Graphs*](https://doi.org/10.37236/2204),
  *Electronic Journal of Combinatorics* 19(4) (2012), P10.
* Gives here: exact nonlinear-code/frustration information in the complete
  bipartite analogue.
* Does not give: the complete-graph quadratic problem.  Its independent
  switching variables on the two vertex classes are exactly the rectangular
  freedom absent from \(x_i x_j\) on \(K_n\).

MacWilliams/Krawtchouk or Delsarte calculations remain plausible, but a weight
enumerator alone does not determine covering radius.  What is needed is
uniform control of the heaviest coset leader, to accuracy \(o(n^{3/2})\), for a
family whose block decomposition introduces many cross-block cycles.  None of
the cited graphical-code results supplies that second-order theorem.

The post-audit one-vertex identity in `RESEARCH_CONTINUATION.md` introduces a
different exact object: the projective spin cube in which every configuration
has weight \((M(B)-|Q_B(x)|)/2\), and distance is Hamming distance modulo
antipodes. Its weighted covering radius reconstructs the optimum extension
value exactly and needs only the energy window of width
\(2\lfloor n/2\rfloor\) below the maximum. This is not the covering radius of
a fixed binary linear code, so the cited cycle- and cocycle-code theorems do
not directly apply. The order-7 counterexample further shows why an ordinary
covering code made only from exact maximizers loses necessary energy-layer
information. A relevant coding theorem would need uniform control of these
energy-weighted covering deficits, not merely a weight enumerator or the
unweighted ground-state code.

Indeed, two optimal order-9 classes in the exact catalogue have identical
absolute-energy histograms, and hence identical scalar partition functions at
every temperature, but different one-vertex extension values. This finite
collision explains why a spin-glass interpolation would need a geometric
order parameter rather than only the optimized scalar free energy. It does not
contradict the thermodynamic-limit results cited above, whose disorder is
sampled rather than adversarially optimized.

## 4. Grothendieck inequalities and semidefinite relaxations

Alon and Naor give a semidefinite algorithm that approximates the cut norm (or
equivalently a rectangular sign bilinear optimization, up to standard
transformations) within a universal constant using Grothendieck's inequality.

* Source: N. Alon and A. Naor, [*Approximating the Cut-Norm via Grothendieck's
  Inequality*](https://doi.org/10.1137/S0097539704441629), *SIAM Journal on
  Computing* 35 (2006), 787--803.
* Gives here: a tractable SDP certificate and constant-factor rounding for a
  closely related bilinear maximum.
* Does not give: a (1+o(1)) comparison, preserve the diagonal identification
  of the two Boolean vectors, or commute with the outer minimization over
  \(A\).  A universal constant loss is a leading-order error here.

Nesterov proves a \(\pi/2\)-type guarantee for SDP relaxation of a class of
nonconvex quadratic problems with diagonal quadratic constraints.

* Source: Y. Nesterov, [*Semidefinite relaxation and nonconvex quadratic
  optimization*](https://doi.org/10.1080/10556789808805690), *Optimization
  Methods and Software* 9 (1998), 141--160.
* Gives here: the standard benchmark for comparing a Boolean rank-one Gram
  matrix \(xx^{\mathsf T}\) with the elliptope
  \(\{R\succeq0:\operatorname{diag}R=1\}\).
* Does not give: an asymptotically vanishing integrality gap for indefinite
  zero-diagonal sign matrices or for the two-sided objective \(M(A)\).

Friedland and Lim formulate a genuinely symmetric Grothendieck inequality for
quadratic forms, including indefinite matrices, and place Nesterov's theorem
in a common conic framework.

* Source: S. Friedland and L.-H. Lim, [*Symmetric Grothendieck inequality*](https://arxiv.org/abs/2003.07345),
  arXiv:2003.07345 (2020).
* Gives here: the right language for a symmetric SDP comparison rather than
  forcing the problem into a rectangular form.
* Does not give: exactness for the adversarially best Seidel matrix.  Its
  universal comparison constants again do not shrink with \(n\).

The semidefinite route is still informative as a **separate relaxed model**.
But a thermodynamic limit of the SDP optimum cannot be transferred to \(F(n)\)
without a pointwise rounding theorem whose additive loss is
\(o(n^{3/2})\).

## 5. Graph limits and exchangeability

Lovasz and Szegedy show that convergent dense graph sequences have graphon
limits characterized by fixed subgraph densities.

* Source: L. Lovasz and B. Szegedy, [*Limits of dense graph sequences*](https://doi.org/10.1016/j.jctb.2006.05.002),
  *Journal of Combinatorial Theory, Series B* 96 (2006), 933--957.
* Gives here: compactness after quotienting by vertex relabeling at the usual
  \(n^2\) normalization.
* Does not give: a nontrivial limit for centered sign arrays whose cut-scale
  discrepancy is only \(n^{3/2}\).  Dividing by \(n^2\) sends precisely that
  information to zero.

Borgs, Chayes, Lovasz, Sos, and Vesztergombi prove equivalences involving
multiway cuts, dense graph ground-state energies, and free energies.

* Source: C. Borgs, J. T. Chayes, L. Lovasz, V. T. Sos, and K. Vesztergombi,
  [*Convergent sequences of dense graphs II: Multiway cuts and statistical
  physics*](https://doi.org/10.4007/annals.2012.176.1.2), *Annals of
  Mathematics* 176 (2012), 151--219.
* Gives here: a rigorous precedent for making graph energies continuous under
  a graph-limit topology.
* Does not give: continuity of a fluctuation-normalized energy.  Their dense
  energy/free-energy normalization is order \(n^2\) with a fixed finite spin
  model; rescaling a sign kernel by \(\sqrt n\) to retain the present signal
  destroys boundedness.

The \(L^p\) sparse-graphon theory permits unbounded graphons under uniform
upper-regularity/\(L^p\) hypotheses and relates metric, quotient, ground-state,
free-energy, and large-deviation convergence.

* Sources: C. Borgs, J. T. Chayes, H. Cohn, and Y. Zhao,
  [*An \(L^p\) theory of sparse graph convergence I*](https://arxiv.org/abs/1401.2906)
  and [*II*](https://arxiv.org/abs/1408.0744).
* Gives here: the most relevant established extension beyond bounded dense
  graphons.
* Does not give: compactness for the fluctuation rescaling here.  A
  \(\{\pm1\}\)-valued dense kernel multiplied by \(\sqrt n\) has diverging
  \(L^p\) norm for every finite \(p\ge1\), so the uniform hypotheses are not
  automatic.

Recent “higher-order graphon” theory studies fluctuations of finite
collections of motif densities in graphon-sampled random graphs.

* Source: A. Chatterjee, S. Dan, and B. B. Bhattacharya,
  [*Higher-Order Graphon Theory: Fluctuations, Degeneracies, and Inference*](https://arxiv.org/abs/2404.13822),
  arXiv:2404.13822 (2024).
* Gives here: proof that second- and higher-order graphon fluctuation theories
  exist for some observables.
* Does not give: control of a maximum over (2^n) correlated spin states or an
  adversarial choice of the underlying graph.  Finite motif CLTs do not imply
  continuity of this exponential-state ground energy.

Thus “take a graphon subsequence” is not a proof program until a topology is
specified that is compact for near-optimizers and continuous at the
\(n^{3/2}\) scale.

## 6. Seidel and conference matrices

A Seidel matrix is exactly a symmetric zero-diagonal \(\{\pm1\}\) matrix, so the
algebraic literature studies the same ambient objects, although usually through
switching classes, spectra, strongly regular graphs, or equiangular lines.

Paley's finite-field construction supplies symmetric conference matrices of
order \(q+1\) for prime powers \(q\equiv1\pmod4\); such a matrix satisfies
\(C^2=qI\).

* Original source: R. E. A. C. Paley, [*On Orthogonal Matrices*](https://doi.org/10.1002/sapm1933121311),
  *Journal of Mathematics and Physics* 12 (1933), 311--320.
* Gives here: the algebraic family behind the repository's asymptotic upper
  bound (together with taking principal submatrices and finding nearby eligible
  orders).
* Does not give: the Boolean maximum \(M(C)\).  The spectral inequality
  \(|x^{\mathsf T}Cx|\le\|C\|_{\mathrm{op}}\|x\|_2^2\) is only a ceiling; equality
  would require a Boolean eigenvector and is not automatic.

There is an important square-order exception.  When \(q=m^2\), with \(m\) an
odd prime power, switching according to a balanced choice of additive
\(\mathbb F_m\)-cosets makes the Paley conference matrix regular with row sum
\(m\).  Undoing the switch gives a Boolean eigenvector, so the full matrix has

\[
 M(C)=\frac{m(m^2+1)}2.
\]

* Sources: J. J. Seidel, [*A survey of
  two-graphs*](https://research.tue.nl/en/publications/a-survey-of-two-graphs/),
  Theorem 13.9, 1976; S. Goryainov, W. H. Haemers, V. V. Kabanov, and L.
  Shalaginov, [*Deza graphs with parameters \((n,k,k-1,a)\) and
  \(\beta=1\)*](https://doi.org/10.1002/jcd.21644), *Journal of Combinatorial
  Designs* 27 (2019), 188--202, Section 4.2.
* Gives here: exact spectral-ceiling attainment on an infinite subfamily of
  full Paley matrices.
* Does not give: \(F(m^2+1)=M(C)\), because \(F\) minimizes over all signings.
  The first case already separates them: the order-10 Paley matrix has maximum
  15 while \(F(10)=13\).  This result neither improves the minimax upper
  constant nor supplies a subsequential lower bound.

Szollosi and Ostergard enumerate Seidel matrices through order 13 up to the
relevant equivalences and determine spectral/algebraic data; related work
classifies small matrices with few eigenvalues.

* Source: F. Szollosi and P. R. J. Ostergard, [*Enumeration of Seidel matrices*](https://doi.org/10.1016/j.ejc.2017.10.009),
  *European Journal of Combinatorics* 69 (2018), 169--184;
  [arXiv:1703.02943](https://arxiv.org/abs/1703.02943).
* Gives here: established symmetry-reduction and invariant machinery for exact
  finite computation.
* Does not give: exact \(F(n)\) from its spectral tables, because \(M(A)\) is not
  a spectral invariant.  Nor can a finite enumeration prove an asymptotic
  limit.

Conference, Hadamard, Paley, or strongly regular constructions can establish
upper bounds along orders.  To prove nonexistence of the limit one would also
need lower bounds for \(F(n)\) on a competing subsequence; different energies
of two explicit constructions are not enough.

## 7. Vector balancing and rounding

The problem can be written as a vector-balancing instance: for every edge
\(e=(i,j)\), let its column be

\[
 v_e=(x_i x_j)_{x\in\{\pm1\}^n/\{x\sim-x\}},
\]

so \(F(n)=\min_{a_e\in\{\pm1\}}\|\sum_e a_ev_e\|_\infty\), with the second copy
of each row implicit in the absolute value.

Banaszczyk proves that sufficiently short Euclidean vectors can be signed so
their sum lies in any symmetric convex body having Gaussian measure at least
(1/2).

* Source: W. Banaszczyk, [*Balancing vectors and Gaussian measures of
  n-dimensional convex bodies*](https://doi.org/10.1002/%28SICI%291098-2418%28199807%2912%3A4%3C351%3A%3AAID-RSA3%3E3.0.CO%3B2-S),
  *Random Structures & Algorithms* 12 (1998), 351--360.
* Gives here: the strongest classical geometric template for rounding many
  coefficients simultaneously.
* Does not give: a lossless rounding theorem for this \(B\).  There are
  exponentially many highly dependent coordinates, the raw columns have
  exponentially large Euclidean norm, and the Gaussian-measure condition for
  the desired thin \(\ell_\infty\) body must be verified after normalization.
  The theorem has no mechanism that preserves a pre-existing near-optimal
  fractional coupling matrix while rounding it edge by edge.

Bansal gives an algorithmic partial-coloring/SDP implementation of entropy
method discrepancy bounds.

* Source: N. Bansal, [*Constructive Algorithms for Discrepancy Minimization*](https://doi.org/10.1109/FOCS.2010.7),
  FOCS 2010, 3--10; [arXiv:1002.2259](https://arxiv.org/abs/1002.2259).
* Gives here: an appropriate computational paradigm for testing whether a
  fractional signing can be progressively saturated while controlling all
  spin rows.
* Does not give: the needed quantitative conclusion automatically.  Any
  accumulated error of order \(n^{3/2}\) changes the leading constant; the MO
  problem requires \(o(n^{3/2})\), uniformly along an amplification or
  interpolation scheme.

This literature becomes actionable only after exploiting special algebraic
dependence among the rows \(x_i x_j\).  Treating them as an arbitrary
exponential constraint system discards the structure that a sharp theorem
would have to use.

## 8. What the literature most strongly suggests

The external results support four precise conclusions, none of which solves
the problem.

1. A thermodynamic-limit proof is plausible only if the interpolation survives
   adversarial optimization; the classical random-disorder interpolation does
   not.
2. A lower-tail disorder principle is more relevant than the usual Parisi
   formula or upper-tail LDP.  For the Ising model, the necessary \(n^2\)-speed
   lower-tail theory is not supplied by the cited work.
3. Convex, SDP, Grothendieck, or vector-balancing relaxations are useful only
   with an additive integrality/rounding loss \(o(n^{3/2})\).  Constant-factor
   approximation is insufficient.
4. Ordinary graphon compactness and ordinary code weight enumerators both
   forget the required second-order extremal information.  A successful limit
   object must retain fluctuations, while a successful coding argument must
   control worst cosets rather than only average weights.

The closest genuinely new target suggested by this map is therefore one of:

* an adversarial finite-temperature interpolation for the pure-sign model;
* an \(n^2\)-speed lower-deviation variational principle uniform enough to
  locate the best Bernoulli disorder;
* an \(o(n^{3/2})\) structured rounding theorem for the cut-character matrix;
* or a second-order covering-radius composition theorem for the augmented cut
  codes of complete graphs.

Any one of these would go materially beyond the cited literature.

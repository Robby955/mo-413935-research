# Second proof-search ledger

Status date: 2026-08-02
Status: **the limit question remains open in this attempt**

This file preserves the new derivations, counterexamples, and stopped routes.
The first-attempt ledger remains in `STATUS.md`; it has not been replaced.
Write

\[
Q_a(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
G_n(a)=\max_{x\in\{\pm1\}^n}|Q_a(x)|,
\qquad m=\binom n2.
\]

For sign couplings, \(F(n)=\min_{a\in\{\pm1\}^m}G_n(a)\).

## 1. Relaxed couplings: exact collapse

### 1.1 Ground state

**Theorem 1 (cube relaxation).**  \(G_n\) is a norm on
\(\mathbb R^m\). Consequently,

\[
\min_{a\in[-1,1]^m}G_n(a)=0,
\]

and the unique minimizer is \(a=0\).

**Proof.** Homogeneity and the triangle inequality are immediate. If
\(X\) is uniform on \(\{\pm1\}^n\), orthogonality of the degree-two Walsh
characters gives

\[
\mathbb E Q_a(X)^2=\sum_{i<j}a_{ij}^2.
\]

Thus \(G_n(a)=0\) implies \(a=0\). The origin belongs to the cube. \(\square\)

This answers two of the proposed relaxation questions sharply: the objective
is convex, and its minimum is not at an extreme point. In fact any convex
feasible set in the original coupling coordinates containing all sign vectors
contains their convex hull, the whole cube, and therefore contains the zero
optimizer when the objective remains \(G_n\).

The ground-state integrality gap is exactly \(F(n)\), hence

\[
\frac{n\sqrt{n-1}}\pi
\le
\min_{\sigma\in\{\pm1\}^m}G_n(\sigma)
-
\min_{a\in[-1,1]^m}G_n(a)
=F(n)
=O(n^{3/2}).
\]

It is therefore of leading order, not \(o(n^{3/2})\).

### 1.2 Finite temperature and its dual

Let \(\mathcal Z=\{s c_x:x\in\{\pm1\}^n, s\in\{\pm1\}\}\), retaining
the duplicated states, and define

\[
L_{n,t}(a)=\frac1t\log\sum_{x,s}
  \exp\bigl(t\langle a,s c_x\rangle\bigr).
\]

**Theorem 2 (finite-temperature collapse and dual).**  The function
\(L_{n,t}\) is strictly convex, its unique cube minimizer is \(0\), and

\[
\min_{a\in[-1,1]^m}L_{n,t}(a)=\frac{(n+1)\log2}{t}.
\]

Moreover its exact entropy dual is

\[
\min_{a\in[-1,1]^m}L_{n,t}(a)
=
\max_{p\in\Delta(\{(x,s)\})}
\left\{
\frac{H(p)}t-
\left\|\mathbb E_p[s c_x]\right\|_1
\right\}.
\]

**Proof.** The Hessian of log-sum-exp is the covariance of the vectors
\(s c_x\) under a full-support Gibbs measure. Their differences span
\(\mathbb R^m\), by the Walsh calculation above, so the Hessian is positive
definite. Symmetry makes the gradient at zero vanish. The value at zero is
the logarithm of \(2^{n+1}\) states. Finally use

\[
\frac1t\log\sum_z e^{t\langle a,z\rangle}
=\max_p\left(\langle a,\mathbb E_pz\rangle+\frac{H(p)}t\right)
\]

and Sion's theorem; minimizing a linear form over the cube gives
\(-\|\mathbb E_pz\|_1\). \(\square\)

At zero temperature the analogous Sion dual is

\[
\min_{a\in[-1,1]^m}\max_{z\in\mathcal Z}\langle a,z\rangle
=\max_p-\|\mathbb E_pz\|_1=0.
\]

Thus first-moment measure duality, the cut polytope, or the elliptope does not
repair the cube relaxation by itself. The uniform state measure already
annihilates every degree-two moment.

### 1.3 What randomized rounding can and cannot do

**Proposition 3 (variance-sensitive rounding).**  Let
\(a\in[-1,1]^m\), set \(V(a)=\sum_e(1-a_e^2)\), and round independently with
\(\mathbb E\sigma_e=a_e\). There is an outcome satisfying

\[
G_n(\sigma)
\le G_n(a)
+\sqrt{2V(a)(n+1)\log2}
+\frac43(n+1)\log2.
\]

**Proof.** For each distinct cut character, Bernstein's inequality applies
to \(\sum_e(\sigma_e-a_e)c_x(e)\), whose total variance is \(V(a)\) and
whose summands have absolute value at most \(2\). A union bound over the
antipodal state set with parameter \((n+1)\log2\) has probability strictly
less than one. Add the resulting uniform noise bound to \(G_n(a)\). \(\square\)

Hence \(V(a)=o(n^2)\) would be enough for an \(o(n^{3/2})\) rounding loss.
The actual cube optimizer has \(V(0)=m=\Theta(n^2)\); the theorem then gives
only a leading-order error, consistently with the exact integrality gap.

Two natural attempts to force saturation also fail:

- On the sphere \(\|a\|_2=\sqrt m\) without coordinate bounds,
  Parseval gives \(G_n(a)\ge\sqrt m\), and the one-edge vector
  \(a=\sqrt m e_{12}\) attains equality. The relaxed value is only \(O(n)\).
- Adding \(|a_e|\le1\) to that sphere forces every \(|a_e|=1\), so it is
  exactly the original discrete problem, not a relaxation.

**Stopped route.** A convex coupling relaxation plus \(o(n^{3/2})\) rounding
cannot start from the cube optimizer. A repair would need a nonconvex but
compact near-saturation constraint together with a new interpolation theorem.

## 2. Temperature scaling and interpolation

### 2.1 The temperature in the proposed functional is already zero temperature

For the functional in the question,

\[
\Phi_{n,\beta}
=\min_{A}\frac1{\beta n^{3/2}}
  \log\sum_x e^{\beta|Q_A(x)|},
\]

the elementary maximum-versus-log-sum-exp inequality gives, uniformly in \(A\),

\[
\frac{F(n)}{n^{3/2}}
\le \Phi_{n,\beta}
\le \frac{F(n)}{n^{3/2}}+\frac{\log2}{\beta\sqrt n}.
\]

Therefore, for every fixed \(\beta>0\), \(\lim_n\Phi_{n,\beta}\) exists if
and only if the original normalized limit exists, and the two limits are
equal. Fixed-temperature smoothing at this scaling is not an easier
intermediate problem.

The genuinely extensive mean-field normalization is

\[
\Psi_{n,\beta}
=\min_A\frac1{\beta n}
 \log\sum_{x,s}\exp\left(\frac{\beta sQ_A(x)}{\sqrt n}\right).
\]

It satisfies

\[
\frac{F(n)}{n^{3/2}}
\le\Psi_{n,\beta}
\le\frac{F(n)}{n^{3/2}}
 +\frac{(n+1)\log2}{\beta n}.
\]

Thus existence of \(\lim_n\Psi_{n,\beta}\) for arbitrarily large fixed
\(\beta\) would force convergence of \(F(n)/n^{3/2}\): its limsup-minus-liminf
would be at most \(\log2/\beta\) for every such \(\beta\).

### 2.2 Exact direct block estimate

Let

\[
h_n(t)=\min_A\log\sum_{x,s}e^{tsQ_A(x)}.
\]

Take optimizing within-block signings, choose the global sign of one block to
anti-align its two one-sided partition sums with those of the other block,
and average independent signs on the cross block. Some cross signing obeys

\[
h_{n+k}(t)
\le h_n(t)+h_k(t)-\log2+nk\log\cosh t.
\]

At \(t=\beta/\sqrt{n+k}\), with \(n/(n+k)\to\alpha\), the last term divided
by \(\beta(n+k)\) tends to

\[
\frac\beta2\alpha(1-\alpha),
\]

which is nonzero at every fixed \(\beta\) and grows in the zero-temperature
limit. This is the precise finite-temperature version of the old annealed
cross-block wall.

**Stopped route.** The direct annealed interpolation is not asymptotically
lossless. A successful interpolation must exploit quenched cancellation and
the outer adversarial minimum together; bounding the cross block separately
cannot work.

### 2.3 The absolute value cannot be dropped

Sign every edge incident with vertex \(1\) positively and every other edge
negatively. If \(T=\sum_{i>1}x_i\), then after fixing \(x_1=1\),

\[
Q_A(x)=\frac n2-\frac{(T-1)^2}{2}.
\]

Consequently \(\max_xQ_A(x)=\lfloor n/2\rfloor\), whereas
\(\max_x|Q_A(x)|=\binom n2\). Any one-sided interpolation can therefore
change the optimization by quadratic order.

## 3. Covariance and semidefinite formulations

Let \(\mathcal E_n=\{R\succeq0:\operatorname{diag}R=1\}\).

### 3.1 The full Gaussian covariance formula is exact but nonconvex

**Proposition 4.** For every signing \(A\),

\[
M(A)=\max_{R\in\mathcal E_n}
\left|\frac2\pi\sum_{i<j}a_{ij}\arcsin R_{ij}\right|.
\]

Gaussian sign rounding shows that every displayed value is at most \(M(A)\).
Conversely \(R=xx^{\mathsf T}\) is feasible and the formula equals \(Q_A(x)\).
Thus the arcsine covariance program is a reformulation, not a relaxation.

### 3.2 A relaxed model whose thermodynamic limit exists

Define

\[
S(A)=\max_{R\in\mathcal E_n}
\left|\frac12\operatorname{tr}(AR)\right|,
\qquad S_n=\min_A S(A).
\]

Rank-one correlation matrices show \(S(A)\ge M(A)\).

**Theorem 5 (adversarial elliptope limit).**

\[
\frac{n\sqrt{n-1}}2\le S_n,
\qquad
\lim_{n\to\infty}\frac{S_n}{n^{3/2}}=\frac12.
\]

At every symmetric conference order, equality holds in the finite lower
bound.

**Proof.** For an arbitrary sign matrix \(A\), put \(q=n-1\) and

\[
R^\pm=\frac12\left(I\pm\frac A{\sqrt q}\right)^2.
\]

Both matrices are positive semidefinite. Since every row of \(A\) has squared
norm \(q\), their diagonals are one, so \(R^\pm\in\mathcal E_n\). Also

\[
R^+-R^-=\frac{2A}{\sqrt q}.
\]

The two SDP objective values differ by

\[
\frac12\operatorname{tr}\bigl(A(R^+-R^-)\bigr)
=n\sqrt q.
\]

At least one has absolute value at least half of this. Conversely, for every
correlation matrix \(R\),
\(|\operatorname{tr}(AR)|\le\|A\|_{\rm op}\operatorname{tr}R
=n\|A\|_{\rm op}\). A conference matrix has operator norm \(\sqrt{n-1}\),
so it attains equality. For arbitrary \(n\), take a principal \(n\)-by-\(n\)
submatrix of a Paley conference matrix of order \(N=n(1+o(1))\). Its operator
norm is at most \(\sqrt{N-1}\), proving the matching asymptotic upper bound.
\(\square\)

The one-sided SDP dual is

\[
\min\left\{\sum_i z_i:\operatorname{Diag}(z)-A/2\succeq0\right\};
\]

the absolute objective takes the maximum of the programs for \(A\) and
\(-A\).

This is a genuine thermodynamic limit for a close relaxation. It does not
settle the Boolean problem: the missing statement is precisely an
asymptotically vanishing integrality gap. The finite gap is already strict at
order \(6\): \(F(6)=5\), whereas \(S_6=3\sqrt5\).

### 3.3 Exact barrier for the linearized Gaussian method

Define \(\Lambda(A)\) as the maximum of

\[
\sum_{i<j}a_{ij}(R^+_{ij}-R^-_{ij})
\]

over \(R^\pm\in\mathcal E_n\), subject to
\(a_{ij}(R^+_{ij}-R^-_{ij})\ge0\) for every edge.

**Theorem 6.**

\[
M(A)\ge\frac{\Lambda(A)}\pi,
\qquad
\Lambda(A)\ge n\sqrt{n-1},
\qquad
\lim_{n\to\infty}
\frac{\min_A\Lambda(A)}{n^{3/2}}=1.
\]

**Proof.** Gaussian sign rounding and \(\arcsin'(r)\ge1\), with the edgewise
orientation constraints, show that the two rounded expected energies differ
by at least \(2\Lambda(A)/\pi\). They both lie in \([-M(A),M(A)]\).
The square covariances \(R^\pm\) above are feasible and have objective
\(n\sqrt{n-1}\). On the other hand,

\[
\Lambda(A)=\frac12\operatorname{tr}(A(R^+-R^-))
\le n\|A\|_{\rm op}.
\]

Conference matrices attain both bounds; Paley principal submatrices fill all
orders asymptotically. \(\square\)

Thus the optimized two-covariance, edgewise-oriented linear certificate
\(\Lambda\) defined above has adversarial asymptotic lower-bound constant
exactly \(1/\pi\). It recovers the old lower bound but cannot improve its
leading constant without changing this certificate. More general covariance
constructions are not excluded.

## 4. The augmented cut code beyond the radius identity

Use bits, with sign \((-1)^b\), and let \(j\) be the all-one edge word. If
\(x_i=(-1)^{u_i}\), then the cut word is
\((\delta u)_{ij}=u_i+u_j\). Therefore

\[
D_n=\{\delta u+s j:u\in\mathbb F_2^n,\ s\in\mathbb F_2\}.
\]

**Theorem 7 (linear code and exact dual).** For \(n\ge3\), \(D_n\) is a
binary \([m,n]\) linear code and

\[
D_n^\perp=
\{H\subseteq E(K_n):\deg_H(v)\equiv0\pmod2\ \forall v,
\ |H|\equiv0\pmod2\}.
\]

For \(n\ge4\), the dual is spanned by 4-cycles, has dimension \(m-n\), and
minimum distance \(4\).

**Proof.** The cut space has dimension \(n-1\). The all-one word is not a
cut for \(n\ge3\), because cuts have even parity on each triangle and \(j\)
has odd parity. Orthogonality to all cuts is exactly even degree at every
vertex; orthogonality to \(j\) is even edge cardinality. The cycle space is
generated by triangles through one fixed vertex. Pairing such triangles
expresses every even-weight cycle-space word as a sum of 4-cycles. \(\square\)

The weight enumerator is consequently

\[
W_{D_n}(z)=\frac12\sum_{k=0}^n\binom nk
\left(z^{k(n-k)}+z^{m-k(n-k)}\right).
\]

### 4.1 Exact signed Eulerian partition formula

Let \(a_e=(-1)^{y_e}\) and

\[
Z_y(\beta)=\sum_{d\in D_n}
\exp\bigl(\beta[m-2\operatorname{wt}(y+d)]\bigr).
\]

Its largest exponent is \(M(A)\), and the coset MacWilliams identity gives

\[
Z_y(\beta)=2^n(\cosh\beta)^mP_y(\tanh\beta),
\]

where

\[
P_y(t)=\sum_{H\in D_n^\perp}(-1)^{|H\cap y|}t^{|H|}.
\]

Although the polynomial is alternating, \(P_y(t)>0\) for \(0\le t<1\), as
the normalized partition function shows. Uniformly in \(y\),

\[
M(A)\le\frac1\beta\log Z_y(\beta)
\le M(A)+\frac{n\log2}\beta.
\]

At each fixed \(\beta>0\), minimizing and dividing by \(n^{3/2}\) is therefore
equivalent to the original limit question. In MacWilliams coordinates the
unknown term is exactly

\[
\min_y\log P_y(\tanh\beta).
\]

It must cancel the explicit order-\(n^2\) term down to order \(n^{3/2}\).
Absolute estimates lose this signing-dependent cancellation.

Averaging \(P_y(t)\) over all \(y\) gives one, hence the best immediate
annealed bound is

\[
F(n)\le \inf_{\beta>0}
\left(\frac{n\log2}\beta+\frac{m\beta}{2}\right)
=n\sqrt{(n-1)\log2},
\]

which is weaker than the Paley constant \(1/2\).

### 4.2 Puncturing and mixed 4-cycles

For \(n\ge4\), projection away from the last vertex maps \(D_n\) onto
\(D_{n-1}\). Its kernel is the last-vertex star, so every old word has two
lifts whose new coordinates are complementary. Hence

\[
\rho(D_{n-1})\le\rho(D_n)
\le\rho(D_{n-1})+\left\lfloor\frac{n-1}{2}\right\rfloor,
\]

or equivalently

\[
F(n-1)+((n-1)\bmod2)\le F(n)\le F(n-1)+n-1.
\]

The \(n=3\) inequalities hold by direct inspection, although that projection
has four-word fibers. This recovers parity and the old one-vertex regularity,
not a second-order limit.

Across a vertex partition, the two internal codewords must share their
triangle-parity bit, and the two possible cross lifts are complementary
rank-one sign matrices. Dually, there are \(\binom n2\binom k2\) alternating
4-cycles supported across two blocks. At \(t=\lambda/\sqrt{n+k}\), their
total absolute fourth-order mass is already macroscopic. A product or cluster
expansion needs quenched cancellation; these terms are not a boundary error.

## 5. High-temperature expansion and its wall

For uniform \(X\in\{\pm1\}^n\), put \(Y=Q_A(X)\) and consider the symmetrized
partition \(\mathbb E\cosh(tY)\), equivalently the augmented state \(sY\).
Exact Walsh counting gives

\[
\mathbb EY^2=m,
\qquad
\mathbb EY^4=3m^2-2m
+24\sum_{C_4}\prod_{e\in C_4}a_e,
\]

and

\[
\kappa_4(Y)=3\operatorname{tr}(A^4)
-2n(n-1)(3n-4).
\]

Also

\[
\operatorname{tr}(A^4)=n(n-1)^2
+\sum_{i\ne j}(A^2_{ij})^2.
\]

The expansion
\[
\log\mathbb E\cosh(tY)
=\frac{m t^2}{2}+\frac{\kappa_4(Y)t^4}{24}+O(t^6)
\]
shows that, for each fixed \(n\), all sufficiently high-temperature
minimizers of the symmetrized partition minimize
\(\operatorname{tr}(A^4)\). The one-sided partition is different because a
signing-dependent cubic triangle term can occur. Spectral
Cauchy--Schwarz gives the lower bound, while for a Paley principal submatrix
inside a conference matrix of order \(N=n(1+o(1))\),
\[
\operatorname{tr}(A^4)
\le\|A\|_{\mathrm{op}}^2\operatorname{tr}(A^2)
\le(N-1)n(n-1).
\]
Consequently,

\[
\lim_n\frac{\min_A\operatorname{tr}(A^4)}{n^3}=1,
\qquad
\lim_n\frac{\min_A\kappa_4(Y)}{n^3}=-3.
\]

This does not determine the next term. Exhaustion of every root-normalized
order-seven signing gives
\[
\min_A\operatorname{tr}(A^4)=342,
\]
but those minimizers have three distinct exact sixth moments:
\[
50{,}781,\qquad53{,}661,\qquad59{,}421.
\]
Thus fourth-order data alone are insufficient.

**Stopped route.** No uniform remainder theorem was obtained that connects
this high-temperature expansion to the ground-state regime.

## 6. Ordinary graph and spectral limits erase the target scale

### 6.1 Cut-distance obstruction

**Theorem 8.** If \(M(A_n)=O(n^{3/2})\), the associated signed step graphons
converge to zero in cut norm. Nevertheless, there are two deterministic
correct-scale sequences with that same limit and separated normalized Boolean
maxima.

**Proof.** For sign vectors (x,y), split the vertices according to whether
\(x_i=y_i\). Symmetry cancels the cross terms and gives

\[
x^{\mathsf T}Ay=2Q_A(z^+)-2Q_A(z^-)
\]

for two vectors in \(\{-1,0,1\}^n\). A multilinear form attains its maximum
absolute value over the cube at a Boolean vertex, so
\(|x^{\mathsf T}Ay|\le4M(A)\). Expanding two indicator vectors in sign
vectors yields

\[
\|W_A\|_\square\le\frac{4M(A)}{n^2}=O(n^{-1/2}).
\]

Now take Paley conference matrices \(C_N\), for which
\(M(C_N)/N^{3/2}\le1/2+o(1)\). Choose
\(k=\lfloor2N^{3/4}\rfloor\) vertices and change every internal coefficient
to \(+1\), obtaining \(B_N\). Only \(O(N^{3/2})\) entries change, so both
graphons still tend to zero and \(M(B_N)=O(N^{3/2})\). Fix the planted spins
to \(+1\) and average the outside spins. Cross and outside terms have mean
zero, so some choice has energy at least \(\binom k2\). Thus

\[
\liminf_N\frac{M(B_N)}{N^{3/2}}\ge2.
\]

The normalized functional is not continuous at the common zero graphon.
\(\square\)

This is structural, not a missing compactness citation: every near-optimal
sequence already has the same ordinary graphon limit.

### 6.2 Limiting spectral distributions are also insufficient

The dense-scale empirical spectrum of \(A/n\) always tends to \(\delta_0\),
because its average squared eigenvalue is \((n-1)/n^2\).

At fluctuation scale, fix \(d>1/2\), start from a conference matrix \(C_N\),
and flip \(O(N^{3/2})\) negative edges until
\(Q_{B_N}(\mathbf1)\ge dN^{3/2}\). Then

\[
\|B_N-C_N\|_F^2=O(N^{3/2}),
\]

so Hoffman--Wielandt gives

\[
W_2\bigl(\mu_{B_N/\sqrt N},\mu_{C_N/\sqrt N}\bigr)\to0.
\]

Both empirical laws tend to
\(\tfrac12(\delta_{-1}+\delta_1)\), while their Boolean maxima differ by a
fixed normalized amount. This rules out limiting empirical spectral measures,
not spectral edges or eigenvector information.

### 6.3 Random signings are rigorously worse than the conference bound

For one infinite iid Rademacher edge array, choose spins greedily:

\[
x_1=1,\qquad
x_j=\operatorname{sgn}\sum_{i<j}a_{ij}x_i.
\]

Each new maximum-endpoint edge column is independent of the past, and
multiplication by the predictable earlier spins preserves its iid Rademacher
law. The resulting energy is therefore a sum of independent variables
distributed as \(|\varepsilon_1+\cdots+\varepsilon_k|\).
Bounded-difference concentration and Borel--Cantelli give almost surely

\[
\liminf_n\frac{M(A_n)}{n^{3/2}}
\ge\frac23\sqrt{\frac2\pi}=0.531923\ldots>\frac12.
\]

This does not lower-bound \(F(n)\); it says iid signings are not asymptotically
optimal. Their signed graphons still converge to zero in cut norm.

## 7. Exact small-order computation

The exhaustive search in `verification/research_exact_small_n.py` gives

| \(n\) | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| \(F(n)\) | 1 | 3 | 4 | 4 | 5 | 9 | 10 | 12 | 13 |
| \(F(n)/n^{3/2}\) | .3536 | .5774 | .5000 | .3578 | .3402 | .4860 | .4419 | .4444 | .4111 |

The reduction first switches every root edge positive and then enumerates one
unlabeled negative-edge graph on the remaining \(n-1\) vertices. Completeness
is checked against the known unlabeled graph counts and hashes of the full
`geng` streams. The energy formula is checked directly; all unreduced labeled
signings are independently enumerated through \(n=6\); optional NetworkX and
Z3 paths supply different decoders/solvers.

At \(n=10\), an exact Paley conference matrix over \(\mathbb F_9\) has Boolean
maximum \(15\), whereas \(F(10)=13\). Thus conference matrices need not even be
finite-order minimizers. The two optimal switching classes at \(n=10\) have
different characteristic polynomials, so the optimizer is not spectrally
unique.

These are reproducible finite computations, not evidence sufficient to infer
convergence, nonconvergence, or a limiting constant. The normalized values are
not monotone and show no certified persistent parity or conference-order gap.

## 8. Distinct live proof programs

### Program A: quenched extensive free energy

**Target lemma.** For every fixed \(\beta>0\), prove existence of
\(\lim_n\Psi_{n,\beta}\), with the theorem valid for arbitrarily large
\(\beta\).

**Why it settles the question.** The sandwich in Section 2 makes the
limsup--liminf gap at most \(\log2/\beta\); sending \(\beta\to\infty\) forces
it to zero.

**Attempt and wall.** Independent cross-coupling interpolation pays
\(\frac\beta2\alpha(1-\alpha)\) per spin. Standard spin-glass Gaussian
integration by parts averages fixed random disorder; it does not commute with
the adversarial minimum over a discrete disorder array.

**Possible repair.** A coupled variational principle in which the minimizing
disorder itself interpolates and the cross term is absorbed by the same order
parameter, rather than bounded separately.

**Assessment.** Structural for annealed/block methods; open for a genuinely
quenched minimax interpolation.

### Program B: near-saturated relaxation and rounding

**Target lemma.** Construct compact feasible sets \(K_n\supseteq\{\pm1\}^m\)
whose optimizers satisfy \(V(a)=o(n^2)\), admit an asymptotic composition law,
and have value no larger than \(F(n)+o(n^{3/2})\).

**Why it settles the question.** Proposition 3 rounds such optimizers at
negligible cost; a limit for the relaxed values transfers to \(F(n)\).

**Attempt and wall.** The cube optimizer is zero with leading integrality gap.
The Euclidean sphere either permits one-coordinate concentration or, with the
box imposed, is exactly discrete.

**Possible repair.** A spread/saturation constraint controlling both
\(\sum(1-a_e^2)\) and coordinate concentration, paired with dependent
rounding adapted to cut characters.

**Assessment.** The convex version is structurally impossible. A useful
nonconvex relaxation is conjectural and risks restating the sign constraint.

### Program C: signed Eulerian polynomial / covering radius

**Target lemma.** For one (equivalently every) fixed \(\beta>0\), prove that

\[
\frac1{\beta n^{3/2}}
\min_y\left[n\log2+m\log\cosh\beta
+\log P_y(\tanh\beta)\right]
\]

has a limit.

**Why it settles the question.** This expression differs from
\(F(n)/n^{3/2}\) by at most \(O(n^{-1/2})\).

**Attempt and wall.** Puncturing recovers only one-step regularity. Code
products are coupled by a common triangle-parity bit and rank-one cross lift.
Mixed minimum-weight dual 4-cycles contribute macroscopically; absolute
Delsarte or cluster estimates erase the necessary sign cancellation.

**Possible repair.** A quenched large-deviation theorem for the signed
even-Eulerian polynomial, or a second-order coset-weight composition inequality
that keeps the minimizing phase.

**Assessment.** Exact and sharply isolated. This is the cleanest algebraic
form of the remaining problem.

### Program D: second-order graph/spectral compactness

**Target lemma.** Build a compact, rate-sensitive topology for
\(A/\sqrt n\) in which \(M(A)/n^{3/2}\) is continuous and in which the
admissible sign arrays are closed under an asymptotically lossless size
composition.

**Why it settles the question.** A variational minimum over a compact limit
space would have an automatic asymptotic value if sampling/composition were
compatible.

**Attempt and wall.** Ordinary cut distance, dense spectra, and even the
fluctuation-scale limiting empirical spectral law provably lose the
functional. Rescaling a bounded graphon by \(\sqrt n\) destroys ordinary
uniform \(L^p\) compactness.

**Possible repair.** Retain spectral edges and eigenvectors, or a Gaussian
process/cut-process limit, not merely an empirical graphon or spectrum.

**Assessment.** Ordinary graph limits are structurally excluded; the required
second-order state space is not yet identified.

### Program E: nonexistence via algebraic subsequences

**Target lemma.** Prove separated lower bounds for \(F(n)/n^{3/2}\) on two
subsequences, not merely construct signings with separated upper bounds.

**Attempt and wall.** Conference, random, and planted signings have different
energies, but \(F(n)\) selects the best signing and ignores the inferior
constructions. Exact values through \(n=10\) fluctuate without a certifiable
subsequence mechanism.

**Assessment.** No current evidence supports a nonexistence claim. Algebraic
orders remain useful for upper bounds only until matching order-specific lower
bounds are found.

## 9. Sharpest remaining wall and best next lemma

No proof of convergence or nonconvergence was found. The strongest asymptotic
positive theorem is the exact thermodynamic limit
\(S_n/n^{3/2}\to1/2\) for the adversarial elliptope relaxation. The
post-audit continuation below adds exact optimizer non-heredity and the
weighted covering-radius Bellman identity. The strongest negative results are
that the natural convex relaxation has a leading integrality gap, ordinary
graphon and empirical-spectral limits erase the objective, and neither exact
minimizers nor their extremizing spin codes form a closed cavity state.

The single best next lemma is Program A's fixed-\(\beta\) statement:

> **Quenched minimax free-energy lemma.** For every \(\beta>0\), the limit
> \(\lim_{n\to\infty}\Psi_{n,\beta}\) exists.

Proving it for an unbounded set of \(\beta\)'s settles the MathOverflow
question immediately. In exact coding language it is the same as controlling
the signing-dependent cancellation in \(\min_y\log P_y\) at
\(t=\tanh(\beta/\sqrt n)\). Any successful proof must absorb, rather than
separately bound, the macroscopic mixed 4-cycle/cross-block contribution.

## 10. Post-audit optimizer-composition search

Exact cross-block optimization over every optimal switching-permutation class
through total order 10 disproved optimizer heredity. The sharp finite row is

\[
 \min_{M(A)=F(2),\ M(B)=F(8)}J(A,B)=15>F(10)=13.
\]

No optimal order-10 signing contains an optimal order-8 principal submatrix.
The obstruction is repaired by exactly two units of internal slack:

\[
 K_{2,8}(1,10)=15,\qquad K_{2,8}(1,12)=13.
\]

Among all 1,044 root-normalized order-8 representatives, 104 have maximum 12
and 68 of those admit an order-10 completion of maximum 13. The complete
proof, class distributions, witnesses, and solver boundary are preserved in
Section 7 of `RESEARCH_CONTINUATION.md`. This kills only composition of
arbitrary exact minimizers. It does not rule out a near-optimal family with
subleading internal slack.

## 11. Energy-weighted projective covering state

For an order-$n$ signing $B$, let $\mathcal P_n$ be the sign cube modulo
antipodes and define

\[
 w_B([x])=\frac{M(B)-|Q_B(x)|}{2},\qquad
 \rho_{\rm w}(B)=\max_{[b]}\min_{[x]}
 \bigl(d_\pm([b],[x])+w_B([x])\bigr).
\]

Pairing the two choices of a new vertex spin gives the exact theorem

\[
 E(B)=M(B)+n-2\rho_{\rm w}(B).
\]

Writing
$\delta_{\rm w}(B)=\lfloor n/2\rfloor-\rho_{\rm w}(B)$ therefore gives the
exact Bellman identity

\[
 F(n+1)=\min_B
 \left(M(B)+(n\bmod2)+2\delta_{\rm w}(B)\right).
\]

The complete optimizer catalogue through order 10 independently verifies the
identity against direct incident-sign enumeration. The ordinary covering
radius of the exact extremizers explains all order-9 and order-10 extension
class distinctions. It is not a sufficient state in general: all six optimal
order-7 classes have extremizer radius 3, but the next energy layer lowers the
weighted radius to 2 in exactly the two classes with extension value 12.

An order-9 collision is sharper. The optimal records `G?qmaw` and `GCpbaw`
have the same projective absolute-energy counts

\[
 \#\{|Q_B|=0,4,8,12\}=(60,111,60,25),
\]

so their scalar partition functions agree for every temperature, while their
extension values are 13 and 15. Thus even the full scalar free-energy curve is
not a closed cavity state; energy-layer geometry is essential.

Only the window

\[
 |Q_B(x)|\ge M(B)-2\lfloor n/2\rfloor
\]

can affect the weighted radius. Exhaustion of all root-normalized residual
graphs through order 8 gives the nontrivial full Pareto frontiers

\[
 \{(5,2),(7,1),(9,0)\}\quad(n=6),\qquad
 \{(10,1),(12,0)\}\quad(n=8)
\]

for $(M,\delta_{\rm w})$. Every displayed order-6 point gives $F(7)=9$ and
both order-8 points give $F(9)=12$, proving that internal energy and covering
deficit can compensate exactly in the Bellman objective.

**New wall.** One-vertex growth now has an exact state, but no asymptotic
theorem controls the Pareto profile
$(M(B),\delta_{\rm w}(B))$ over near-optimal signings. Any ground-state cavity
program that keeps only maximizers is structurally false; a viable compression
must retain quantitative energy-layer information.

## 12. Nonlinear Gaussian stability search

The earlier covariance proof discarded all curvature of the arcsine law.  For

\[
 R_s^\pm=\frac{(I\pm sA)^2}{1+s^2(n-1)},
\]

retaining the exact edgewise difference yields Theorem 17 of
`RESEARCH_CONTINUATION.md`.  At $s=1/\sqrt{n-1}$ it gives

\[
 M(A)\ge\frac{n(n-1)}\pi\arcsin\frac1{\sqrt{n-1}}
 +\frac{\|A^2-(n-1)I\|_F^2}
 {8\pi(n-1)(n-2)^{3/2}}.
\]

This has two genuine consequences.  It proves $F(21)\ge32$ after an exact
parity rounding, while the parity theorem for
$\Delta(A)=\|A^2-(n-1)I\|_F^2$ also proves $F(20)\ge30$.  It forces any
sequence attaining the lower constant $1/\pi$ to satisfy
$\|A^2-(n-1)I\|_F=o(n^2)$.

The earlier draft called this condition conference-like.  That wording was
too strong: a uniform random signing has
$\mathbb E\Delta=n(n-1)(n-2)=\Theta(n^3)$ and satisfies the stated
normalization in probability.  The proved conclusion is only vanishing
normalized Gram defect.

**Leading-order wall.**  The universal nonlinear gain is only
$\Theta(\sqrt n)$.  The trace-four excess can contribute at leading order for
matrices far from conference structure, but a minimizer may have excess
$o(n^4)$.  Therefore this refinement does not separate the liminf from
$1/\pi$ and does not prove convergence.

**Closed subroute.**  Optimizing the full parameter $s$ does not repair this.
If $G_n$ is the minimax value of the complete nonlinear right side after
optimizing $s$, then Theorem 22 proves
$G_n/n^{3/2}\to1/\pi$.  Principal blocks of Paley conference matrices of
order $n+o(n)$ give the matching upper barrier uniformly in $s$.  Any further
Gaussian route must combine laws in a way not reducible to this one-parameter
certificate or its ordinary convex combinations.

## 13. Multivertex weighted covering state

For fixed internal blocks $B,D$, Theorem 18 identifies the exact cross-block
optimization with the weighted covering radius of the projective rank-one
code.  Proposition 19 then gives

\[
 J(B,D)\le L(B,D)+\sqrt{2nk(\Xi(B,D)+\log4)},
\]

where $\Xi$ is the Gaussian-weighted effective size of the internal
near-ground rank-one words.

**What changed.**  This is an exact two-block Bellman state rather than an
upper bound that separates the internal and cross energies.  It identifies
precisely which energy layers matter and quantifies their effective entropy.

**Structural wall.**  The generic estimate
$\Xi\le(n+k-2)\log2$ reproduces the same leading
$\sqrt{nk(n+k)}$ cost as the old random cross-block argument.  Requiring
$\Xi=o(n+k)$ is also too crude: balanced repeated composition can tolerate a
linear entropy only if its constant matches the internal energy profile.

**Disproved target, retained for the ledger.**  The earlier proposal was to
find internal blocks $B,D$ with

\[
 L(B,D)+\sqrt{2nk(\Xi(B,D)+\log4)}
 \le(F(n)^{2/3}+F(k)^{2/3})^{3/2}
 +O((n+k)^{3/2-\varepsilon}).
\]

It would have given the power-saving near-subadditivity of $F^{2/3}$.
Proposition 23 proves it impossible: the left side is always at least
$\sqrt{2nk(n+k)\log2}$.  For balanced blocks this has constant
$2\sqrt{\log2}$, strictly above the required asymptotic ceiling $\sqrt2$.
Retaining the discarded cross term in the same all-state Hoeffding sum has
the same leading obstruction.

**Replacement target.**  Fix a low-operator-norm cross seed $C_0$ and retain
its row-and-column switching orbit.  With
$c(U)=|\langle C_0,U\rangle|$ and internal profile $h(R)$, the bad shifts at
target $K$ are exactly

\[
 \mathcal B_K=\bigcup_R R\{U:c(U)>K-h(R)\}.
\]

Prove $\mathcal B_K\ne\mathcal R_{n,k}$ at the power-saving composition
threshold.  This weighted noncoverage problem preserves the higher-order
dependence erased by the union bound.

The exact one-vertex increment additionally implies density-one control:
Bellman-optimal predecessors have both internal slack and covering deficit at
most $g(N)\sqrt N$ for all but $O(N/g(N))$ orders in $[N,2N]$.  This is not
the decisive scalar obstruction: Section 15 records a nonconvergent
countermodel with $O(\sqrt n)$ Bellman cost at every order.

## 14. Order-nine state-compression experiment

Complete enumeration gives $\mathcal B_9=\{(12,0)\}$ and a collision beyond
the earlier scalar-partition-function example.  The records `GHOgmo` and
`Gxd?Dc` agree in their complete energy histograms and in the ordered
pair-distance law of exact maximizers, but have extension values 15 and 17.

Thus the following state compressions are now finitely false:

1. the scalar maximum alone;
2. the exact-maximizer covering code alone, already false at order 7;
3. the complete scalar partition-function curve;
4. the scalar energy histogram plus the two-point law of exact maximizers.

The full energy-coloured two-point distribution separates all deficits among
the 12,346 rooted order-9 records.  This suggests a hierarchy of
energy-coloured overlap laws analogous to richer spin-glass order parameters.
It is evidence, not a theorem of sufficiency: higher-order collisions may
appear at larger orders, and even a complete finite overlap hierarchy still
needs an asymptotic composition law.

**Best next ground-state lemma.**  Prove the switching-orbit noncoverage
statement in Section 13 for a composable near-optimal family, first for
balanced doubling if necessary.  Any first-moment replacement must be checked
against the exact order-four diagnostic, where two good orbit shifts exist
despite a mean of $17/4$ violated constraints.

## 15. Exact scalar Bellman wall

Write

\[
 r_n=F(n+1)-F(n)-(n\bmod2)=\sigma_n+2\delta_n.
\]

The density-one result was initially interpreted as leaving sparse exceptional
orders.  A stronger countermodel disproves that diagnosis.  Starting from the
certified prefix through order ten and parity-rounding the increments of

\[
 t^{3/2}\left[\frac25+\frac1{25}\sin\log\log(t+e^2)\right]
\]

produces an integer sequence satisfying parity, monotonicity, the Gaussian
lower bound, the trivial upper bound, both scalar cavity consequences, and
$r_n=O(\sqrt n)$ at every order, while its normalized liminf and limsup are
$0.36$ and $0.44$.  This is an abstract scalar countermodel, not a family of
sign matrices.

The precise one-vertex diagnostic is

\[
 A_N=N^{-3/2}\sum_{n=N}^{2N-1}r_n.
\]

The normalized sequence converges if and only if $A_N$ converges; if
$A_N\to\Lambda$, its limit is $\Lambda/(2^{3/2}-1)$.  A useful sufficient
target is stabilization of the Cesaro mean

\[
 \frac1N\sum_{n\le N}\frac{r_n}{\sqrt n}.
\]

Actual matrix geometry must establish this stabilization; scalar magnitude
bounds cannot.

## 16. Order-ten temperature phases

The full rooted order-ten catalogue contains 274,668 records and 6,012
distinct absolute-energy histograms.  Exact polynomial comparison reduces
the minimax partition function to three phases.  With
$z=4\sinh^2t$, their difference polynomials factor as

\[
 P_0-P_1=8z^2(z-2)(z+1)(z+3)(z+4)^2,
\]

\[
 P_1-P_2=4z^2(z+4)^2(z^3+z^2-10z-8).
\]

Thus a conference signing with maximum 15 minimizes for
$0<t<0.658478948\ldots$, an intermediate maximum-15 signing minimizes until
$t=0.792460762\ldots$, and an order-ten ground-state signing with maximum 13
minimizes thereafter.  At $t=0$ every signing ties.  This exact finite result
shows that the outer minimizer changes with temperature; it does not prove an
asymptotic phase transition.

## 17. Negative-replica recursion and transport wall

The annealed-normalized negative moment

\[
 \mathcal G_n(q,t)=\log\mathbb E_A
 \left[\frac{Z_A(t)}{2^{n+1}(\cosh t)^{\binom n2}}\right]^{-q}
\]

is exactly superadditive in the order at fixed $(q,t)$.  The proof pairs the
two orientations of one internal block, uses Jensen twice, and loses no
normalizing factor.  This is the first exact block law found that is directly
compatible with the outer soft minimum.

It does not live at the extensive diagonal. Boolean reverse
hypercontractivity transports $\mathcal G_n/q$ on the exact invariant curve
$(1+q)\tanh^2t=\theta$, but its inequality points in the wrong direction.
The first proposed repair was a power-saving upper bound (PT) on that
transport defect. Section 20 now proves that PT is false with an
$\Omega(n^2)$ defect. The earlier finite data were adverse evidence, not the
disproof.

The negative-replica recursion itself survives. What failed was the scalar
transport that discarded the rectangular quotient and the conditional
relative-switching law. Sections 21--22 identify that missing state exactly.

## 18. Exact orders 11 and 12

The order-10 extension catalogue first gives the catalogue consequence
$F(11)\ge15$.  A complete exact nauty pass scans and hashes all 12,005,168
unlabeled residual graphs on ten vertices, filters the analytically eligible
classes, and finds no signing with maximum at most 15.  The mathematically reduced stream

```text
geng -q -D6 10 20:22
```

contains 2,153,606 records and gives the same result with the same evaluator;
sampled records at each eligible edge count are recomputed by a separate
adjacency formula.  With the
explicit witnesses,

\[
 F(11)=17,\qquad F(12)=18.
\]

These are computer-assisted exact values with a nauty completeness boundary.
They are not asymptotic progress by themselves.  The external claim that the
order-11 normalized value is a record for $n\ge5$ is false: the order-7 ratio
is larger.

### Exact orders 13 and 14

The formerly external billion-record claim has now been ingested and replayed
without using the supplied binaries.  A fail-closed verifier compiles a
fixed-order threshold scanner, relays all eight nauty shards through it while
hashing the exact input bytes, checks producer and consumer statuses, and
requires committed counts and SHA-256 digests.  The shard counts sum to
$1{,}018{,}997{,}864$.

Only `JCpVdXyxpz?` and `JCpdUg{[dM?` survive the exact order-12 threshold
$M\le18$.  Independent full-cube computation gives $M=18$ and direct
one-vertex extension minimum $24$ for both.  Every other predecessor has
$M\ge20$, so the Bellman identity gives $F(13)\ge20$.  The Paley
$C_{14}$ principal witness gives $F(13)\le20$; heredity and parity followed
by $M(C_{14})=21$ give

\[
 \boxed{F(13)=20,\qquad F(14)=21}.
\]

This remains a computer-assisted result with a nauty completeness boundary.
The scanner source is `verification/order12_threshold_scan.c`; the verifier
and pinned shard receipts are in
`verification/research_order13_certify.py`.

The accompanying layer bank was not accepted wholesale.  A separate fresh
reconstruction does verify the bracket $F(15)\in\{25,27\}$, but the banked
scripts do not reproduce that tower and the clean driver is not yet in this
repository.  Exact $F(15)$ is still blocked on the order-13
predecessors with $(M,\delta_{\rm w})=(24,0)$.
An external sweep of that branch is ongoing.  Its v2 artifacts are not in the
ingested bank, so neither a partial count nor any eventual order-15 verdict is
promoted here.

## 19. Square-order Paley audit

The reported Boolean ceiling attainment for Paley conference matrices of
order $m^2+1$ is correct but already known as the regular-conference
construction.  Additive $\mathbb F_m$-cosets give a sign eigenvector of
eigenvalue $m$, hence matrix maximum $m(m^2+1)/2$.

The original report overstated the consequence.  This is a theorem about the
specified full Paley matrix, not about $F(m^2+1)$.  At the first case,
$m=3$, the matrix maximum is 15 while $F(10)=13$.  The result neither improves
the asymptotic upper bound nor forces a subsequential value of the minimax
sequence.  It is retained as a structural correction and as a warning that
conference ceiling attainment depends on the order and switching class.

### Paley Fourier leakage and the separate rigidity wall

For every odd prime power $q\equiv1\pmod4$, diagonalizing the Paley core by
additive characters gives an exact identity.  For Boolean
$f:\mathbb F_q\to\{\pm1\}$, let $S(f)=\sum f$, and let $E_+(f),E_-(f)$ be
the unitary Fourier energies in the two nonzero frequency halves selected by
the signs of the quadratic Gauss-sum eigenvalues.  Put
$W(f)=\min(E_+(f),E_-(f))$.  Then

\[
 M(C_{q+1})=\frac{(q+1)\sqrt q}{2}
 -\sqrt q\min_f\left[
 W(f)+\frac{(|S(f)|-\sqrt q)^2}{2q}
 \right].
\]

The proof is Parseval followed by exact optimization of the infinity sign.
Therefore Paley Boolean alignment tends to the spectral ceiling if and only
if one can find $f_q$ with $|S(f_q)|=o(q)$ and minority-half leakage
$W(f_q)=o(q)$.  This is the precise character-sum target suggested by the
selected prime-field data.

A fresh source-built exhaustive scan verifies

\[
 M(C_6)=5,\quad M(C_{14})=21,\quad M(C_{18})=33,\quad
 M(C_{30})=75.
\]

The corresponding Boolean-to-spectral ratios are
$\sqrt5/3,3/\sqrt{13},11/(3\sqrt{17}),5/\sqrt{29}$ and increase on these
four selected orders.  This proves no monotonicity or fitted decay law.

Even asymptotic Paley alignment is only an upper-bound construction for $F$.
To force the original limit to be $1/2$, this route separately needs

\[
 F(q+1)\ge(1-o(1))M(C_{q+1})
\]

along a multiplicatively dense set of Paley prime orders.  That rigidity
statement is the hard minimax wall.  Exact optimality at orders 6 and 14 does
not approach a proof, and order 10 is a counterexample to unqualified
conference optimality.

## 20. Scalar negative-replica transport is closed

For fixed $\beta,\theta>0$, put

\[
 t_n=\frac\beta{\sqrt n},\qquad
 q_n=\frac\theta{\tanh^2t_n}-1,
\]

and let

\[
 \Delta_n=\frac{q_{2n}}{q_n}\mathcal G_n(q_n,t_n)
 -\mathcal G_n(q_{2n},t_{2n}).
\]

The augmented partition function is invariant under vertex switching and
global coefficient negation. For $n\ge3$ this free symmetry group has size
$2^n$, so the disorder cube splits into $2^{\binom n2-n}$ equal-size
orbits. Comparing the negative moment with its smallest value, then bounding
that value above and below by an $F(n)$-optimal state, gives the exact finite
inequality

\[
\Delta_n\ge q_{2n}\left[
 \binom n2(\log\cosh t_n-\log\cosh t_{2n})
 -(t_n-t_{2n})F(n)-n\log2
 -\frac{(\binom n2-n)\log2}{q_n}
\right].
\]

Using only $\limsup F(n)/n^{3/2}\le1/2$ therefore yields

\[
 \liminf_n\frac{\Delta_n}{n^2}
 \ge \frac{2\theta}{\beta^2}
 \left[\frac{\beta^2}{8}
 -\frac\beta2\left(1-\frac1{\sqrt2}\right)-\log2\right]-\log2.
\]

At $(\beta,\theta)=(4,8)$ the right side is

\[
 \sqrt2-2\log2=0.027919201253\ldots>0.
\]

Thus PT is false. More generally the bound is positive for sufficiently large
$\theta$ whenever $\beta>3.012373175\ldots$. Exact class-weighted data
through order 9 reproduce a much larger positive defect on the $(4,8)$ grid,
but the theorem does not depend on those computations.

There is also an exact mechanism for the gap. On a Boolean $m$-cube with
$m\ge1$, let $\theta>0$ and restrict to nonnegative times with
$q_s=\theta e^{2s}-1>0$. Put $f_s=T_{e^{-s}}f_0>0$ and

\[
 \mu_s=\frac{f_s^{-q_s}}{\mathbb E f_s^{-q_s}}\nu,
 \qquad D_s=D(\mu_s\Vert\nu),\qquad
 H_s=\frac1{q_s}\log\mathbb E f_s^{-q_s}.
\]

Then

\[
 \frac{dH_s}{ds}=-R_s,\qquad
 R_s=\mathbb E_{\mu_s}\frac{Lf_s}{f_s}
 -\frac{2(q_s+1)}{q_s^2}D_s,
\]

and the one-coordinate calculation plus entropy tensorization gives

\[
 \boxed{R_s\ge\frac{2D_s^2}{3q_sm}.}
\]

Equality holds only for constant $f_s$. Consequently, power-saving transport
would force a power-saving entropy loss in the escort disorder law. This
stability result diagnoses scalar non-saturation, but it does not recover the
joint block alignment discarded by scalar transport.

## 21. Exact conditional-alignment chain

For a split $N=n+k$, let $R_{n,k}$ be the rectangular row-and-column
switching code and set

\[
 C'_{n,k}=D_n\oplus D_k\oplus R_{n,k}.
\]

For $n,k\ge3$, $D_N\subset C'_{n,k}$ and
$|C'_{n,k}/D_N|=2^{N-1}$. The coarse quotient factors into the two graph
quotients and the rectangular quotient. If $h$ is a coarse coset, let $r_h$
be the conditional fine density on its fiber and
$K_q(h)=\mathbb E_{\rm fiber}r_h^{-q}$. Then

\[
 \boxed{
 \mathcal G_N=\mathcal G_n+\mathcal G_k
 +\mathcal G_{n,k}^{\square}+\mathcal A_{n,k}},
 \qquad
 \mathcal A_{n,k}=\log\mathbb E_{\widehat U_q}K_q(h)\ge0.
\]

Here $\widehat U_q$ is the factorized coarse law tilted by its negative
moment. In a relative gauge, the fiber is parameterized by
$\alpha\in\{\pm1\}^n/\{\pm\mathbf1\}$,
$\beta\in\{\pm1\}^k/\{\pm\mathbf1\}$, and $\tau\in\{\pm1\}$ through

\[
 Y_{\alpha,\beta,\tau}=
 \begin{pmatrix}
 D_\alpha AD_\alpha&C\\
 C^{\mathsf T}&\tau D_\beta BD_\beta
 \end{pmatrix}.
\]

With augmented graph state sums and the ordinary rectangular state sum,

\[
 r_{A,B,C}(\alpha,\beta,\tau)
 =\frac{2^{N+1}Z_Y}{Z_AZ_BZ_C^\square},
 \qquad \mathbb E_{\alpha,\beta,\tau}r=1.
\]

The normalization $2^{N+1}$ is tied to this state-sum convention; codeword
sums use different powers of two.

If
$\mathcal A^\infty(t)=\lim_{q\to\infty}\mathcal A(q,t)/q$, then sequential
finite-space limits give

\[
 \lim_{t\to\infty}\frac{\mathcal A^\infty(t)}t
 =F(n)+F(k)+B_\square(n,k)-F(N).
\]

Since $B_\square(n,n)\ge n\,\mathbb E|\varepsilon_1+\cdots+\varepsilon_n|$,

\[
 \liminf_n\frac{2F(n)+B_\square(n,n)-F(2n)}{n^{3/2}}
 \ge\frac2\pi+\sqrt{\frac2\pi}-\sqrt2
 =0.0202907708\ldots.
\]

Alignment is therefore a leading ground-state quantity. The rectangular
marginal alone cannot compensate the scalar transport gap in the large
parameter regime.

A proposed alignment-transport inequality was audited after this derivation.
If $a_j=\mathcal G_j(q_j,t_j)/q_j$, its purported new left side is identically
$a_N-a_n-a_k$. It is a correct sufficient target only because it restates the
desired almost-superadditivity. The actual open task is to derive such a bound
from a smaller, controlled geometric state.

## 22. First alignment Hamiltonian and exact scalar collision

The dual codes in the graph and rectangular quotients have minimum weight
four. Writing $P=D_\alpha AD_\alpha$ and $Q=D_\beta BD_\beta$, the coefficient
of $u^4$ in $\log r$ is

\[
 \mathcal H_4=\frac12\left[
 \operatorname{tr}(P^2CC^{\mathsf T})
 +\operatorname{tr}(Q^2C^{\mathsf T}C)
 +\tau\operatorname{tr}(PCQC^{\mathsf T})
 -nk(N-2)\right].
\]

Its exact Walsh form is

\[
\begin{aligned}
\mathcal H_4={}&
 \sum_{i<j}(A^2)_{ij}(CC^{\mathsf T})_{ij}\alpha_i\alpha_j
 +\sum_{a<b}(B^2)_{ab}(C^{\mathsf T}C)_{ab}\beta_a\beta_b\\
&+\tau\sum_{i<j,a<b}a_{ij}b_{ab}
 (C_{ia}C_{jb}+C_{ib}C_{ja})
 \alpha_i\alpha_j\beta_a\beta_b.
\end{aligned}
\]

The characters are orthogonal, so its variance is the sum of the squared
coefficients. In particular the mixed contribution is at least

\[
 4\binom n2\left\lfloor\frac{(k-1)^2}{4}\right\rfloor,
\]

and symmetrically with $n,k$ interchanged. Thus $\mathcal H_4$ is nonconstant
for $n,k\ge3$, and the conditional alignment is strictly nontrivial at
sufficiently high temperature. This does not permit truncating the expansion
at $u^4$ on the mean-field diagonal: higher Eulerian terms are not uniformly
small when $u\asymp N^{-1/2}$ and $q\asymp N$.

An exact split-$2+4$ collision makes the lost information concrete. At
$u=3/5$, two triples $(A,B,C)$ have identical complete local absolute-energy
histograms, hence identical graph and rectangular scalar partition curves at
every temperature, but

\[
 K_2^{(1)}=
 \frac{196585091273040100817}{133610891512185651200},
 \qquad
 K_2^{(2)}=\frac{6723290161}{5922841600}.
\]

Their $\mathcal H_4$ laws are respectively
$\{-12:4,-4:4,0:16,4:4,12:4\}$ and
$\{-4:8,0:16,4:8\}$. Therefore even all three local scalar free-energy curves
do not determine block alignment. Here $\dim D_2=1$, so the 32 displayed
relative gauges cover the true 16-point fiber twice; the normalized moments
are unchanged. The first viable non-tautological target is
a variational or composition theorem for the relative-switching Hamiltonian
hierarchy, beginning with these mixed traces and retaining enough higher
connected Eulerian clusters to be uniform on the mean-field scale.

## 23. Relative-gauge convolution and the microcanonical profile theorem

The conditional alignment law is not merely an unexplained correction term.
It is the pushforward of three independent local laws under the balanced
homomorphism

\[
 ((\sigma,[z]),(\eta,[w]),[xy^{\mathsf T}])
 \mapsto([zx],[wy],\sigma\eta).
\]

This gives an exact zero-temperature max-plus identity. If $L$ is the sum of
the two graph maxima and the rectangular maximum, then

\[
 L-M(Y_g)=\min_{\pi(a,b,r)=g}(d_A(a)+d_B(b)+d_C(r)).
\]

There are $2^{N-1}$ equal fibers. Therefore the $2^{N-1}$-st smallest local
triple deficit, denoted $\Lambda$, is always a realizable composition gain:

\[
 F(N)\le M(A)+M(B)+B_C-\Lambda.
\]

This differs structurally from the closed weighted union-bound route. It
counts actual near-maximal triples and does not randomize the cross block.
There is no universal leading floor in the argument itself. The exact
exponential relaxation is

\[
 \Lambda\ge\sup_{t>0}
 \frac{(N-1)\log2-\log D_A(t)-\log D_B(t)-\log D_C(t)}t.
\]

The scalar relaxation is also sharp at its stated information boundary. If
one knows only the unlabeled deficit multiset and equal fiber sizes, the
smallest $2^{N-1}$ values can be placed one per abstract fiber. Additive group
labels are the only possible source of a stronger general theorem. The known
$2+4$ collision has $\Lambda=0$ for both examples but true gains 4 and 2.

This creates a clean decision tree.

1. Prove a large-deviation bound showing that a composable near-optimal family
   has fewer than $2^{N-1}$ triples below the power-saving target deficit.
2. If that quantile is too small, retain extensive Fourier labels. The
   conditional Fourier coefficient factors into three local Gibbs
   correlations, and the mixed four-cycle trace is the first layer of this
   hierarchy.

The finite observation that exact maximizer counts remain small is relevant
but insufficient. A constant exact-maximizer count guarantees only the first
positive energy-lattice gap. The entire $\Theta(N^{3/2})$ near-maximal window
must have a controlled entropy profile.

## 24. Least-nonresidue interval construction closes Paley leakage

For prime $p\equiv1\pmod4$, take the Boolean half-interval with sum one. Its
nonzero Fourier magnitudes are exact Dirichlet-kernel values, and

\[
 |\widehat f(r)|^2\le p/\|r\|_p^2.
\]

If $\ell(p)$ is the least quadratic nonresidue, the minority-half leakage is
at most $2p/(\ell(p)-1)$. Quadratic reciprocity makes every integer through
$j$ a residue when

\[
 p\equiv1\pmod{8\prod_{\lambda\le j,\ \lambda\ {\rm odd\ prime}}\lambda}.
\]

The prime number theorem in each fixed progression permits a stage
construction with levels tending to infinity and consecutive prime ratios
tending to one. Thus Paley Boolean maxima reach their spectral ceilings on a
multiplicatively dense prime sequence.

This is a real theorem but not a minimax lower bound. It closes only the
Fourier half of the value-specific program. The remaining rigidity target is

\[
 F(p_j+1)\ge(1-o(1))M(C_{p_j+1})
\]

on that same sequence. If proved, monotonicity would force the original limit
to be $1/2$.

Two failed strengthenings are recorded rather than erased.

- Exact one-half Fourier support is impossible for a nonconstant Boolean
  function at prime order. Cyclotomic divisibility gives the qualitative
  obstruction, and an algebraic-norm argument gives
  $W(f)\ge4(p-1)^2/p^3$.
- The same interval is not uniform over all primes. Along
  $p\equiv5\pmod{12}$, the residue modes $\pm1$ and nonresidue modes $\pm3$
  force $\liminf W/p\ge8/(9\pi^2)$.

The arithmetic progression is therefore part of the construction, not a
technical convenience. Novelty is stated only as new to this repository;
specialist priority has not been established.

## 25. Fixed-half discrepancy, exact profile tails, and conference rigidity

This entry consolidates three routes that became precise only after the
relative-profile theorem. It records the exact reductions, retains the failed
concentration shortcuts with their proper scope, and separates the remaining
scalar and geometric lemmas.

### Fixed-half cut discrepancy is equivalent at the target scale

Put $m=\binom n2$ and

\[
 H(n)=\min_{\substack{G\text{ on }[n]\\
                       e(G)=\lfloor m/2\rfloor}}
 \max_{S\subseteq[n]}
 \left|e_G(S,S^{\mathsf c})-\frac12|S||S^{\mathsf c}|\right|.
\]

If the negative coefficients of a signing $A$ are the edges of $G$, write
$t=\sum_{i<j}a_{ij}=m-2e(G)$. For the spin which is positive on $S$ and
negative on $S^{\mathsf c}$,

\[
 Q_A(x_S)=t-2|S||S^{\mathsf c}|+4e_G(S,S^{\mathsf c}).       \tag{25.1}
\]

On the fixed-half layer, $t=t_0\in\{0,1\}$, so (25.1) identifies the two
maxima up to one unit. For an arbitrary optimizer, random switching replaces
$t$ by $Q_A(z)$, and

\[
 \mathbb E_zQ_A(z)=0,\qquad \mathbb E_zQ_A(z)^2=m.
\]

Some switch therefore has $|t|\le\sqrt m$. Flipping
$|t-t_0|/2$ coefficient signs reaches the fixed-half layer and changes every
quadratic energy by at most $|t-t_0|$. Consequently

\[
 \boxed{-1\le 4H(n)-F(n)\le\sqrt{\binom n2}+2.}             \tag{25.2}
\]

Thus $F(n)/n^{3/2}$ converges if and only if $4H(n)/n^{3/2}$ converges, and
the limits agree. Existing dense cut-discrepancy theorems recover the
$n^{3/2}$ scale but supply neither a sharp constant nor a relation between
orders. The reformulation is exact at the required scale, not a convergence
theorem.

### The scalar cumulant criterion is exact

For a split $N=n+k$, retain the local augmented graph energies $e_A,e_B$ and
the absolute rectangular energy $e_C$. If the three local states are uniform
and independent, put

\[
 X=e_A+e_B+e_C.
\]

The product space has $2^{2N-2}$ states and the relative-gauge map has
$2^{N-1}$ equal fibers. Strict exponential counting gives, for every $t>0$,

\[
 \boxed{
 \min_gM(Y_g)\le
 \frac{(N-1)\log2+\log\mathbb E e^{tX}}{t}.}              \tag{25.3}
\]

Equivalently, with

\[
 \kappa_N(\beta)=\frac1N\log\mathbb E
 \exp\left(\frac{\beta X}{\sqrt N}\right),
\]

a threshold $\theta_NN^{3/2}$ is certified whenever

\[
 \kappa_N(\beta)-\beta\theta_N
 \le-\left(1-\frac1N\right)\log2.                         \tag{25.4}
\]

This is exactly the canonical form of the microcanonical order-statistic
bound. Indeed,

\[
 \mathbb E e^{tX}
 =2^{-(2N-2)}e^{tL}D_A(t)D_B(t)D_C(t),
\]

so optimizing (25.3) recovers the earlier deficit-profile formula without an
additional relaxation. The complete cumulant function remains a live scalar
state.

Two standard proxies do not retain enough of it.

1. For balanced blocks $n=k=r$, any global quadratic majorant
   \[
    \log\mathbb E e^{t(X-\mathbb EX)}\le v_rt^2/2
   \]
   must have $v_r\ge r(r-1)$ by the second derivative at zero. Its optimized
   normalized certificate is therefore at least
   \[
    \sqrt{\frac{\log2}{2}}=0.588705011\ldots>\frac12.
   \]
   This closes Gaussian, subgaussian, or purely quadratic-cumulant
   certificates only. It is not a lower bound on the true profile quantile
   or on $F$; higher cumulants can change the far tail.

2. Ordinary degree-two hypercontractivity gives
   \[
    \|Q_A\|_p\le(p-1)\sqrt{\binom n2}.
   \]
   At threshold $cn^{3/2}$, optimizing in $p$ gives only
   $\exp[-\Theta_c(\sqrt n)]$, whereas (25.3) needs an
   $\exp[-\Theta(n)]$ tail. This rules out the standard degree-two moment
   certificate, not spectrum-sensitive inequalities, a full cumulant
   estimate, or the additive relative-gauge geometry.

These scoped walls are retained; neither is evidence that the exact scalar
criterion fails.

### Exact finite calibration of the scalar profile

The deterministic calibration exhausts every balanced vertex split of four
stored optimal witnesses: two rooted encodings of the same order-12 optimum,
one order-13 principal submatrix of $C_{14}$, and $C_{14}$. Let
$L=M(A)+M(B)+B_C$, let $\Lambda$ be the microcanonical gain, and let
$G=L-F(N)$ be the exact total relative-gauge gain. If $E_*$ denotes the
smallest energy with the parity of $\binom N2$ strictly above

\[
 [F(n)^{2/3}+F(k)^{2/3}]^{3/2},
\]

the exact results are

\[
\begin{array}{c|c|c|c|c|c|c}
N&\text{split}&\#\text{ splits}&\Lambda&L-\Lambda&G-\Lambda
 &E_*:\ \#\text{ target triples}/2^{N-1}\\ \hline
12&6+6&462\ \text{each}&0\ldots10&24\ldots26&6\ldots8
 &16:\ 213312\ldots269122/2048\\
13&6+7&1716&2\ldots10&28\ldots30&8\ldots10
 &20:\ 436864\ldots495632/4096\\
14&7+7&1716&6\ldots10&31\ldots33&10\ldots12
 &27:\ 231581\ldots305465/8192
\end{array}                                                \tag{25.5}
\]

The target count is the exact number of product triples with energy at least
$E_*$, equivalently strictly above the real target. In every tested split it
exceeds one fiber, so the unlabeled scalar pigeonhole theorem does not certify
the zero-error target on these witnesses. Additive labels recover a further
6--12 energy units.

This is finite evidence only. The best scalar excesses over the target are
$9.858\ldots$, $8.476\ldots$, and $5.544\ldots$; an $O(N)$ error would be a
valid power saving at scale $N^{3/2}$. The computation covers the listed
witnesses, not every optimal switching class at orders 13 and 14. Its true
gain also uses the separately certified values $F(12)=18$, $F(13)=20$, and
$F(14)=21$. Therefore (25.5) neither closes nor disproves the scalar route.

A convention audit separated two valid but different scalar theorems.  The
exact max-plus theorem uses augmented graph states and projective absolute
rectangular states.  If instead the graph factors are projective $|Q|$ states
and the rectangular factor consists of signed full-spin pairs, the map

\[
 ([z],[w],x,y)\longmapsto([zx],[wy],x_1)
\]

still has $2^{N-1}$ states per relative gauge.  Its fiber maximum only
dominates, rather than equals, $M(Y_g)$: a maximizing spin for $Y_g$ can be
injected into its fiber after choosing the two graph signs by absolute value.
Therefore its $2^{N-1}$-st deficit order statistic
$\widetilde\Lambda$ also satisfies

\[
 \min_gM(Y_g)\le L-\widetilde\Lambda.                  \tag{25.5a}
\]

On the standard $C_{14}$ split the one-sided model gives
$\widetilde\Lambda=8$ and 596440 target states.  The exact model gives
$\Lambda=10$ and 304908 target states, with true gain 22.  The first numbers
are not a corruption, but they do not verify the exact max-plus identity and
are weaker in this calibration.  Neither finite count disproves a scalar
power-saving theorem.

The one-sided swapped profile does have a leading raw floor. For balanced
$r+r$ blocks, if $\widetilde U=L-\widetilde\Lambda$ and
$\mu_r=\mathbb E|\sum_{i=1}^r\varepsilon_i|$, then

\[
 \widetilde U\ge\max\{M(A),M(B)\}+r\mu_r-3r.           \tag{25.5b}
\]

The proof counts at least $2^r$ signed cross pairs above
$r\mu_r-3r$ using a maximizing row spin and its radius-one Hamming ball.
Combining them with one maximizing graph state and all projective states of
the other graph reaches the order-statistic rank $2^{2r-1}$. Thus a balanced
near-subadditive theorem based on this alternate raw statistic requires

\[
 \liminf F(r)/r^{3/2}\ge
 \sqrt{2/\pi}/(2^{3/2}-1)=0.4363775564\ldots.           \tag{25.5c}
\]

This closes the swapped raw route only under a limiting constant below that
threshold. It remains compatible with the value-$1/2$ hypothesis and says
nothing adverse about the distinct exact max-plus statistic.

At fixed temperature, the normalized augmented and swapped product moment
generating functions have ratio in $[1/2,4]$. If
$U_*=\inf_{t>0}((N-1)\log2+\log P_*(t))/t$, positivity of the means and the
bound $U_*\le L$ give

\[
 U_{\rm aug}-L/(N-1)\le U_{\rm sw}
 \le U_{\rm aug}+2L/(N-1).                            \tag{25.5d}
\]

Thus the canonical exponential relaxations differ by only $O(\sqrt N)$ at
the mean-field scale. The threshold (25.5c) applies to both canonical
certificates, but not to the raw exact augmented order statistic.

### Conference alignment is an exact eigenspace problem

Let $C=C^{\mathsf T}$ be a conference matrix of order $n$, put
$r=\sqrt{n-1}$, and write

\[
 E_\pm=\ker(C\mp rI),\qquad
 \alpha_\pm=\frac1{\sqrt n}
 \sup_{\substack{v\in E_\pm\\\|v\|_2=1}}\|v\|_1.
\]

For an orthogonal projection $P$ onto a subspace $E$,

\[
 \max_{x\in\{\pm1\}^n}x^{\mathsf T}Px
 =\left(\sup_{\substack{v\in E\\\|v\|_2=1}}\|v\|_1\right)^2.
\]

Applying this to $P_\pm=(I\pm C/r)/2$ gives the exact formula

\[
 \boxed{
 \frac{2M(C)}{n\sqrt{n-1}}
 =\max\{2\alpha_+^2-1,\,2\alpha_-^2-1\}.}                 \tag{25.6}
\]

Both eigenspaces are necessary because $M$ contains an absolute value.
Spectral saturation is equivalent to at least one eigenspace containing unit
vectors with $\ell_1$ norm $(1-o(1))\sqrt n$. Thus a uniform gap of both
$\ell_1$ suprema below $\sqrt n$ is the exact possible obstruction for an
arbitrary conference family. The dense Paley construction proves alignment
for a special family; it does not prove anything adversarial about $F$.

On the multiplicatively dense aligned Paley sequence $n_j=p_j+1$,

\[
 F(n_j)\ge(1-o(1))M(C_{n_j})                              \tag{25.7}
\]

is equivalent to

\[
 \lim_{n\to\infty}\frac{F(n)}{n^{3/2}}=\frac12.          \tag{25.8}
\]

Indeed, (25.7), monotonicity, and $n_{j+1}/n_j\to1$ propagate the subsequence
lower bound to every order, while the known limsup gives the reverse bound.
Conversely, (25.8), Paley alignment, and $F(n_j)\le M(C_{n_j})$ force
$F(n_j)/M(C_{n_j})\to1$. Dense Paley minimax rigidity is therefore not an
easier intermediate statement: it is the entire value-$1/2$ closure in an
equivalent form.

### Sharp remaining lemmas

The scalar route now has one precise target. For proportional splits, find a
composable near-optimal family $A,B,C$ and $\delta>0$ such that fewer than
$2^{N-1}$ product triples exceed

\[
 [F(n)^{2/3}+F(k)^{2/3}]^{3/2}+O(N^{3/2-\delta}).         \tag{25.9}
\]

Equivalently, prove (25.4) at this threshold using the full local cumulant
function, not a quadratic proxy. This would give power-saving
near-subadditivity of $F^{2/3}$.

If the scalar quantile in (25.9) is too large, the sharp geometric target is
a power-saving composition or variational theorem for the labelled
relative-gauge convolution. It must control how the local near-maximal states
occupy the additive fibers, or equivalently enough of the factored Fourier
correlation hierarchy to prove

\[
 \min_gM(Y_g)
 \le [F(n)^{2/3}+F(k)^{2/3}]^{3/2}+O(N^{3/2-\delta}).      \tag{25.10}
\]

The mixed four-cycle Hamiltonian is only the first nonzero layer; higher
connected Eulerian terms cannot be discarded without a uniform estimate on
the mean-field scale. Equations (25.9) and (25.10) are respectively the sharp
scalar and geometric settling lemmas. The conference rigidity condition
(25.7) is a separate value-specific equivalent formulation of full closure,
not a shortcut around either composition problem.

## 26. Labeled-shell gain and the fixed-density cross floor

### A strict Fourier improvement over the scalar profile

For a threshold $s$, let $b_s(g)$ count local product triples in gauge fiber
$g$ whose total deficit is less than $s$, and put
$K=2^{N-1}$ and $\mu_s=K^{-1}\sum_g b_s(g)$. With normalized group Fourier
coefficients, Parseval gives

\[
 V_s=\sum_{\chi\ne1}|\widehat b_s(\chi)|^2
 =K^{-1}\sum_g(b_s(g)-\mu_s)^2.
\]

The exact inequality

\[
 \min_gb_s(g)\le\mu_s-\sqrt{V_s/(K-1)}                 \tag{26.1}
\]

follows from nonnegativity and the mean-zero constraint. If the right side is
less than one, integrality produces an empty fiber and hence
$\min_gM(Y_g)\le L-s$. This contains the scalar pigeonhole theorem and is
strict at total shell size $K$ whenever the occupancy is nonconstant.

Every coefficient factors into three labeled local shell transforms:

\[
 \widehat b_s(I,J,\epsilon)=K^{-1}
 \sum_{d_A+d_B+d_C<s}
 A_{d_A}(I,\epsilon)B_{d_B}(J,\epsilon)C_{d_C}(I,J).   \tag{26.2}
\]

For the exact $2+4$ collision at $s=2$, the scalar shell contains exactly
$K=32$ triples and hence gives no gain. The occupancy law is
$0^8,1^{16},2^8$, so $V=1/2$ and (26.1) certifies gain two. The actual gains
are four and two.

For the standard balanced $C_{14}$ split, the target shell has mean occupancy
$304908/8192=37.2202\ldots$, variance $151.4691\ldots$, 8159 nonzero
nontrivial coefficients, and occupancy range $0$ through $87$. Exactly one of
the 8192 fibers is empty, but (26.1) gives only $37.0842\ldots$. Thus variance is a real strict
improvement at the scalar boundary but is not the missing asymptotic theorem.
The next target is a low-tail or higher-moment inequality for the factored
labeled shell transform.

The higher-moment target is now exact. For any integer-valued occupancy $b$
and finite $A\subset\mathbb Z_{\ge1}$,

\[
 P_A(x)=(1-x)\prod_{a\in A}\frac{(x-a)(x-a-1)}{a(a+1)}
\]

equals one at zero and is nonpositive at every positive integer. Hence
$\Pr[b=0]\ge\mathbb E P_A(b)$. Its moments factor as zero-sum products of the
Fourier coefficients in (26.2). Equivalently, if
$H_r=(m_{i+j+1}-m_{i+j})_{i,j\le r}$, an everywhere-positive occupancy forces
$H_r\succeq0$.

At the balanced $C_{14}$ target, the nine adjacent-root locations

\[
 9,17,26,36,46,56,67,76,86
\]

give a positive exact degree-19 certificate. The corresponding order-nine
localizing quadratic form has negative numerator
$-584163517696745929254421003286532$. Thus higher labeled moments really do
recover the unique empty fiber missed by variance.

This does not yet scale automatically. The even and odd halves of the
$\operatorname{Bin}(D+1,1/2)$ multiplicity table have identical moments
through degree $D$, although the even half has one zero and the odd half none.
Consequently moments through degree $N-1$ cannot decide vacancy on an abstract
$2^{N-1}$-point fiber space. A signing-specific degree-$O(N)$ localizing
estimate, with controlled error after the shell factorization, is the new
settling target.

### Exact half-density does not remove the rectangular floor

For a bipartite cross graph encoded by $S\in\{\pm1\}^{n\times k}$ with total
$t$, its two-sided XOR cut deviation is exactly

\[
 D_\oplus=(|t|+\|S\|_{\infty\to1})/4.                 \tag{26.3}
\]

Moreover

\[
 \|S\|_{\infty\to1}\ge\max\{n\mu_k,k\mu_n\}.        \tag{26.4}
\]

If $R_{n,k}(t_0)$ minimizes this norm at fixed total $t_0$, row/column
switching followed by at most $|t-t_0|/2$ entry flips proves

\[
 R_{n,k}\le R_{n,k}(t_0)
 \le R_{n,k}+\sqrt{nk}+|t_0|.                          \tag{26.5}
\]

Thus enforcing the global half-density constraint costs only $O(\sqrt{nk})$,
but the leading rectangular floor survives. This is a no-go for triangle or
cut-norm arguments that control the cross term separately. It is not a no-go
for the labeled cancellation measured by (26.1)--(26.2).

## 27. Fixed-half blow-ups and the uniform Hadamard test

For equal clouds of size $k$, cloud-union cuts give

\[
 d(L)\ge\max_S\left|\sum_{i\in S,j\notin S}
 (e(V_i,V_j)-k^2/2)\right|.
\]

Thus a complete/empty blow-up of a fixed-half base graph $G$ has
$d(L)\ge k^2d(G)$. For $n\ge3$ its internal cloud signs can always be chosen
to make the full graph fixed-half. Seidel switching followed by $r$ edge
repairs preserves the lower bound up to $(r+1)/2$. Therefore $O(N)$ density
repair cannot remove the classical $\sqrt k$ normalized loss.

Orthogonal blocks evade the cloud-union argument, so no general Hadamard
no-go follows. The exact smallest uniform test still fails: for the negative
$C_5$ signing $A$, every symmetric Hadamard $H$ of order four and fixed-half
order-four signing $D$ satisfy

\[
 M(A\otimes H+I_5\otimes D)\ge44,
\]

whereas lossless four-fold scaling would give 32. Equality 44 is attained.
This exhausts only the common-$H$, common-$D$ family and leaves nonuniform
orthogonal block designs open.

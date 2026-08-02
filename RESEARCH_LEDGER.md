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

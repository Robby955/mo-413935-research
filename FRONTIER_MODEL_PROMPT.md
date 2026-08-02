# Frontier-model research prompt: MathOverflow 413935

You are being asked to make a serious, independent attempt at MathOverflow
Question 413935, “Min max of a quadratic form of plus-minus ones.”  The goal is
new mathematics, not a summary or a generic list of techniques.

If a repository checkout accompanies this prompt, read every tracked file
before beginning, especially `AUDIT.md`, `RESEARCH_LEDGER.md`,
`paper/second_attempt.tex`, `LITERATURE.md`, and every script under
`verification/`.  Treat even the audited statements below as hypotheses until
you have checked their proofs yourself.

## 1. Problem

For a symmetric zero-diagonal matrix $A=(a_{ij})$ with
$a_{ij}\in\{-1,1\}$, define

\[
 Q_A(x)=\sum_{1\le i<j\le n}a_{ij}x_ix_j
       =\frac12x^{\mathsf T}Ax,
 \qquad
 M(A)=\max_{x\in\{-1,1\}^n}|Q_A(x)|,
\]

and

\[
 F(n)=\min_A M(A).
\]

Determine whether

\[
 \lim_{n\to\infty}\frac{F(n)}{n^{3/2}}
\]

exists.  Its value is not required.  A proof of either convergence or
nonconvergence is acceptable.

## 2. Audited baseline

The following statements have complete proofs in the current research
package, but you must independently check any one you use:

\[
 \frac1\pi
 \le \liminf_n\frac{F(n)}{n^{3/2}}
 \le \limsup_n\frac{F(n)}{n^{3/2}}
 \le\frac12.
\]

More precisely,

\[
 F(n)\ge \frac{n\sqrt{n-1}}\pi,
 \qquad
 F(n)\le \frac12n\sqrt{N-1}
\]

whenever $A$ is chosen as an $n$-vertex principal submatrix of a
symmetric conference matrix of order $N$.  Paley orders with
$N=n(1+o(1))$ give the limsup bound.

There is one-vertex regularity

\[
 F(n)\le F(n+1)\le F(n)+n.
\]

In fact the augmented-code puncturing argument gives

\[
 F(n-1)+((n-1)\bmod 2)\le F(n)\le F(n-1)+n-1.
\]

Let $m=\binom n2$, let the cut word corresponding to $x$ be
$c_x=(x_ix_j)_{i<j}$, interpreted in binary coordinates, and augment by
global negation.  Then

\[
 F(n)=m-2\rho(D_n),
\]

where $D_n\subseteq\mathbb F_2^m$ is the augmented cut code.

## 3. Results already obtained in the second attempt

These results narrow the search.  Do not spend the attempt merely
rediscovering them.

### 3.1 The direct cube relaxation collapses

For

\[
 G_n(a)=\max_x|\langle a,c_x\rangle|,
 \qquad a\in\mathbb R^m,
\]

$G_n$ is a norm because

\[
 \mathbb E_x\langle a,c_x\rangle^2=\sum_ea_e^2.
\]

Consequently

\[
 \min_{a\in[-1,1]^m}G_n(a)=0
\]

with unique optimizer $a=0$.  The discrete-to-cube integrality gap is
exactly $F(n)=\Theta(n^{3/2})$.

At fixed unscaled positive temperature,

\[
 L_{n,t}(a)=\frac1t\log\sum_{x,s}
 e^{t\langle a,sc_x\rangle}
\]

is strictly convex with the same unique optimizer $0$.  Its entropy dual is

\[
 \min_{a\in[-1,1]^m}L_{n,t}(a)
 =\max_p\left\{\frac{H(p)}t-
 \left\|\mathbb E_p[sc_x]\right\|_1\right\}.
\]

Independent rounding with $\mathbb E\sigma_e=a_e$ yields some sign vector
$\sigma$ satisfying

\[
 G_n(\sigma)\le G_n(a)
 +\sqrt{2V(a)(n+1)\log2}
 +\frac43(n+1)\log2,
 \quad
 V(a)=\sum_e(1-a_e^2).
\]

This is $o(n^{3/2})$ only when $V(a)=o(n^2)$, whereas
$V(0)=\Theta(n^2)$.
A useful relaxation therefore needs a nonconvex near-saturation mechanism or
a fundamentally different lift.

### 3.2 Use the extensive temperature scaling

The initially proposed functional

\[
 \Phi_{n,\beta}=\min_A\frac1{\beta n^{3/2}}
 \log\sum_xe^{\beta|Q_A(x)|}
\]

is already uniformly within $\log2/(\beta\sqrt n)$ of
$F(n)/n^{3/2}$.  Existence of its fixed-$\beta$ limit is therefore
equivalent to the original problem, rather than an easier intermediate step.

The meaningful extensive free energy is

\[
 \Psi_{n,\beta}
 =\min_A\frac1{\beta n}
 \log\sum_{x\in\{\pm1\}^n,\ s\in\{\pm1\}}
 \exp\left(\frac{\beta sQ_A(x)}{\sqrt n}\right).
\]

Uniformly in $n$,

\[
 \frac{F(n)}{n^{3/2}}
 \le\Psi_{n,\beta}
 \le\frac{F(n)}{n^{3/2}}+\frac{(n+1)\log2}{\beta n}.
\]

Therefore the following is the cleanest settling lemma:

> **Quenched minimax free-energy lemma.**  For every fixed
> $\beta>0$, the limit $\lim_n\Psi_{n,\beta}$ exists.

It is enough to prove this for an unbounded set of fixed $\beta$'s.  The
limsup-minus-liminf gap of $F(n)/n^{3/2}$ would then be at most
$\log2/\beta$ for arbitrarily large $\beta$, hence zero.

The direct annealed block calculation does not suffice.  If

\[
 h_n(t)=\min_A\log\sum_{x,s}e^{tsQ_A(x)},
\]

then the best direct independent cross-block construction currently gives

\[
 h_{n+k}(t)\le h_n(t)+h_k(t)-\log2+nk\log\cosh t.
\]

At $t=\beta/\sqrt{n+k}$, the normalized cross term tends to

\[
 \frac\beta2\alpha(1-\alpha),
 \qquad \alpha=\lim\frac n{n+k},
\]

which is leading order.  A successful interpolation must absorb this cross
term into an optimized quenched variational structure, not bound it
separately.

### 3.3 A close SDP model does have a limit

Let

\[
 \mathcal E_n=\{R\succeq0:\operatorname{diag}R=1\},
 \quad
 S_n=\min_A\max_{R\in\mathcal E_n}
 \left|\frac12\operatorname{tr}(AR)\right|.
\]

Then

\[
 \frac{n\sqrt{n-1}}2\le S_n,
 \qquad
 \lim_n\frac{S_n}{n^{3/2}}=\frac12.
\]

The lower bound follows from the two feasible correlation matrices

\[
 R^\pm=\frac12\left(I\pm\frac A{\sqrt{n-1}}\right)^2.
\]

Conference matrices and nearby Paley principal submatrices give the matching
upper bound.  This does not solve the problem: a constant SDP-to-Boolean gap
is leading order.  Any use of this relaxation must produce an additive
$o(n^{3/2})$ comparison after the outer minimization.

### 3.4 Exact coding formulation

For $n\ge3$,

\[
 D_n=\{\delta u+s\mathbf1:u\in\mathbb F_2^n,
 s\in\mathbb F_2\}
\]

is a binary $[m,n]$ linear code, and

\[
 D_n^\perp=
 \{H\subseteq E(K_n):\deg_H(v)\equiv0\pmod2\ \forall v,
 |H|\equiv0\pmod2\}.
\]

For $n\ge4$, this dual is generated by 4-cycles.  If a signing is encoded by
$y\in\mathbb F_2^m$, define

\[
 P_y(t)=\sum_{H\in D_n^\perp}(-1)^{|H\cap y|}t^{|H|}.
\]

The exact coset MacWilliams identity is

\[
 Z_y(\theta)
 =\sum_{d\in D_n}e^{\theta[m-2\operatorname{wt}(y+d)]}
 =2^n(\cosh\theta)^mP_y(\tanh\theta).
\]

Although \(P_y\) is signed, \(P_y(t)>0\) for \(0\le t<1\).  At the extensive
scale, the open problem becomes control of

\[
 \min_y\log P_y\left(\tanh\frac\beta{\sqrt n}\right),
\]

on its natural order-\(n\) scale.  Absolute coefficient estimates destroy the
signing-dependent cancellation needed for an asymptotically lossless
composition law.  Mixed cross-block 4-cycles already contribute
macroscopically on that scale.

### 3.5 Ordinary compactness loses the signal

Every sequence with $M(A_n)=O(n^{3/2})$ converges to the zero signed graphon
in cut norm.  Deterministic planted-block perturbations can have the same
zero graphon limit while their normalized Boolean maxima are separated by a
fixed constant.  Likewise, two sequences can have the same limiting empirical
spectrum at scale $A/\sqrt n$ but separated Boolean maxima.  Ordinary
graphons and limiting empirical spectral measures are therefore insufficient;
any limit object must retain second-order process, spectral-edge, or
eigenvector information.

### 3.6 Exact finite data

Auditable exhaustive search gives

\[
\begin{array}{c|rrrrrrrrr}
n&2&3&4&5&6&7&8&9&10\\ \hline
F(n)&1&3&4&4&5&9&10&12&13.
\end{array}
\]

The normalized values are not monotone.  At order 10, a Paley conference
matrix has Boolean maximum 15, while the optimum is 13.  Finite computation is
guidance only; it cannot establish an asymptotic subsequence.

### 3.7 Results from the continued parallel attack

Define

\[
 r_n(t)=h_n(t)-\log2-m_n\log\cosh t,
 \qquad
 R_n(\beta)=r_n(\beta/\sqrt n).
\]

The block calculation has the exact annealed-normalized form

\[
 r_{n+k}(t)\le r_n(t)+r_k(t),
\]

and hence

\[
 R_{n+k}(\beta)
 \le R_n\left(\beta\sqrt{\frac n{n+k}}\right)
 +R_k\left(\beta\sqrt{\frac k{n+k}}\right).
\]

This is genuine subadditivity, but at contracted child temperatures.  An
explicit oscillating scalar countermodel satisfies this recursion, the exact
one-vertex cavity bounds, temperature monotonicity, and all current scalar
barriers.  These inequalities alone therefore do not force convergence.

The square covariances also give an entropy-aware finite-temperature bound.
For

\[
 R^\pm=\frac12\left(I\pm\frac A{\sqrt{n-1}}\right)^2,
 \qquad
 \Sigma^\pm_\lambda=(1-\lambda)I+\lambda R^\pm,
\]

Gaussian sign rounding, relative-entropy data processing, and the Gibbs
variational principle give

\[
 \begin{aligned}
 \mathcal F_{n,\beta}(A)\ge{}&
 \frac\lambda\pi\sqrt{\frac{n-1}{n}}
 +\frac{(n+1)\log2}{\beta n}\\
 &+\frac{\log\det\Sigma^+_\lambda+
 \log\det\Sigma^-_\lambda}{4\beta n}.
 \end{aligned}
\]

Uniformly in $A$, the determinant product is at least $d(\lambda)^n$, where

\[
 d(\lambda)=
 \begin{cases}
 (1-\lambda/2)^2,&\lambda\le2/3,\\
 2\lambda(1-\lambda),&\lambda\ge2/3.
 \end{cases}
\]

This supplies an explicit energy--entropy lower barrier, but it still tends
only to the ground-state constant $1/\pi$ as $\beta\to\infty$.

In coding variables, $p_n(u)=\min_yP_y(u)$ is exactly the least quotient
atom of iid edge noise, normalized by the quotient size.  Thus
$-\log p_n$ is a reverse max-divergence, $p_n$ is nonincreasing in $u$, and

\[
 (1-u^2)^{n/2}p_n(u)\le p_{n+1}(u)\le p_n(u).
\]

The normalized Eulerian polynomial is multi-affine, so it has the exact box
relaxation

\[
 \min_{a_e\in\{\pm1\}}\mathcal P_G(ua)
 =\min_{w\in[-u,u]^E}\mathcal P_G(w).
\]

Deletion--contraction is exact only after retaining exponentially many
boundary/parity sectors; scalar $P$ or $(P,R)$ does not close under vertex
extension.

Finally, if $A$ is uniform over all $2^{m_n}$ disorders, then for every
$q>0$,

\[
 \Psi_{n,\beta}
 \le-\frac1{q\beta n}\log\mathbb E_AZ_A^{-q}
 \le\Psi_{n,\beta}+\frac{m_n\log2}{q\beta n}.
\]

Taking $q=\lambda n$ leaves error $\log2/(2\lambda\beta)+o(1)$.  Existence
of this negative-replica free-energy limit for an unbounded set of
$\lambda$ is therefore a second precise settling lemma.

### 3.8 Exact hereditary cavity inequalities

For a vertex partition $S\sqcup T$, $|S|=s\ge1$, $|T|=k\ge1$, put
$C=A_{S,T}$. For every fixed $x\in\{\pm1\}^S$, exact Walsh orthogonality
gives

\[
 M(A)^2\ge Q_{A[S]}(x)^2+\|C^{\mathsf T}x\|_2^2+\binom k2.
\]

Choosing $x$ to maximize the induced signing and using parity yields

\[
 F(s+k)^2\ge F(s)^2+\binom k2+k(s\bmod2).
\]

Pairing the two completions $y$ and $-y$ gives a different exact inequality:

\[
 M(A)\ge |Q_{A[S]}(x)|
 +\mathbb E_y|\langle C^{\mathsf T}x,y\rangle|
 \ge |Q_{A[S]}(x)|+\frac1{\sqrt2}\|C^{\mathsf T}x\|_2.
\]

If $s$ is odd, coordinatewise monotonicity of the Rademacher $L^1$ norm
sharpens the scalar consequence to

\[
 F(s+k)\ge F(s)+\mu_k,
 \qquad
 \mu_k=\frac{k}{2^{k-1}}
 \binom{k-1}{\lfloor(k-1)/2\rfloor}.
\]

These are genuine matrix-level hereditary refinements, but their universal
terms are only $O(k^2)$ in $M^2$ or $O(\sqrt k)$ in $M$. They do not reach
the leading $n^{3/2}$ scale and, being lower bounds, do not construct an
upper composition.

They expose a second precise settling route. Define $H(n)=F(n)^{2/3}$.
Exact subadditivity fails because $H(4)>2H(2)$, but any uniform power-saving
estimate

\[
 H(n+m)\le H(n)+H(m)+O((n+m)^{1-\delta}),
 \qquad \delta>0,
\]

would force $H(n)/n$, and hence $F(n)/n^{3/2}$, to converge by balanced
binary composition. An unspecified $o(n+m)$ defect is not enough without a
summability condition across scales.

## 4. Known failed mechanisms that require a material change

Do not simply repeat any of the following:

1. Independent random cross edges in an ordinary block sum: their uniform
   discrepancy is itself order $(n+k)^{3/2}$.
2. Rank-one clone blow-ups: the main term scales as $k^2F(n)$, with an
   additional leading internal-block cost.
3. Direct annealed partition estimates: the $nk\log\cosh t$ cost is
   macroscopic at the extensive temperature.
4. Treating the conference spectral ceiling as attainable on Boolean vectors:
   explicit orders 6, 14, and 18 disprove this, and order 10 is not even
   optimal among signings.
5. Convexifying only the coupling cube and then rounding: the optimizer is
   zero and the integrality gap is leading order.
6. Ordinary graphon compactness or empirical spectral convergence: explicit
   deterministic counterexamples erase the target functional.
7. Dropping the absolute value: there are signings with
   $\max Q_A=O(n)$ but $\max|Q_A|=\Theta(n^2)$.
8. Edge-separable saturation penalties and near-saturated fractional cross
   blocks: both necessarily pay a leading $N^{3/2}$ cross cost.
9. Bounded-rank local lifts: a $k\times k$ sign block of rank $r$ has
   rectangular discrepancy at least $k^2/r^{3/2}$.
10. Canonical Seidel/Kronecker amplification: an optimal order-five signing
    has an exact order-25 counterexample, and the standard iterated tensor has
    exponentially worsening normalized energy.
11. Scalar use of the hereditary cavity inequalities above: their universal
    rewards are subleading and point in the lower-bound direction, while the
    missing convergence mechanism is an optimizer-compatible upper
    composition.

A route using one of these ideas is valid only if it supplies a genuinely new
mechanism that removes the displayed leading error.

## 5. Priority attack programs

Pursue several programs, but make a decision about which one is strongest.

### Program A: quenched minimax interpolation

Try to prove the fixed-$\beta$ lemma for $\Psi_{n,\beta}$.  Possible new
ingredients include:

- a Guerra--Toninelli-type interpolation in which the minimizing disorder is
  part of the order parameter;
- a variational formula over overlap distributions, exchangeable spin laws,
  or joint spin/disorder empirical measures;
- an approximate subadditivity theorem after optimizing the cross disorder
  conditionally on both blocks;
- a fractional-moment or lower-tail large-deviation formulation for the best
  Bernoulli disorder;
- a mixed-strategy minimax relaxation with a provably negligible purification
  gap.
- the negative-replica quantity with $q=\lambda n$, for which disorder is
  averaged but low-partition signings are exponentially favored.

For every proposed interpolation, calculate the cross derivative exactly and
show its sign or cancellation.  Do not exchange minimum, expectation,
logarithm, or limit without a theorem that applies to the actual compact
sets.

### Program B: signed Eulerian polynomial composition

Derive an exact vertex-extension or two-block identity for $P_y(t)$ that
retains the minimizing phase $y$.  Seek one of:

- a second-order subadditive inequality for $\min_y\log P_y$;
- a transfer operator whose normalized log spectral radius has a limit;
- a quenched large-deviation principle for signed even-Eulerian subgraphs;
- a deletion/contraction, cavity, or cluster expansion uniform at
  $t=\Theta(n^{-1/2})$;
- a Delsarte/Krawtchouk bound sensitive to worst cosets rather than only the
  weight enumerator.

Track mixed 4-cycles explicitly.  An error of $O(n^{3/2})$ in ground-state
units is not negligible.

### Program C: near-saturated relaxation and dependent rounding

Construct a compact nonconvex relaxation whose optimizers satisfy

\[
 \sum_e(1-a_e^2)=o(n^2),
\]

which composes across sizes and differs from the discrete value by
$o(n^{3/2})$.  Alternatively, prove that a broad natural class of such
relaxations cannot work.  Investigate dependent rounding that exploits the
degree-two Walsh structure rather than treating the $2^{n-1}$ constraints
as arbitrary.

### Program D: asymptotically lossless amplification

For every large base order $n$, try to construct signings at a dense set of
larger orders $N$ satisfying

\[
 \frac{M(A_N)}{N^{3/2}}
 \le \frac{F(n)}{n^{3/2}}+o_n(1).
\]

Test tensor products, randomized lifts conditioned on row sums, multi-block
designs, conference completions, and dependent cross blocks.  Prove uniform
Boolean bounds; operator-norm control alone is insufficient.

A particularly sharp target is to prove, for some $\delta>0$,

\[
 F(n+m)^{2/3}
 \le F(n)^{2/3}+F(m)^{2/3}+O((n+m)^{1-\delta}).
\]

A dyadically summable relative defect is essential for this argument; a power
saving is a clean sufficient condition. Calculate the error under repeated
balanced composition rather than calling a generic $o(n+m)$ term negligible.

### Program E: nonexistence

Do not assume convergence.  A nonexistence proof requires separated lower
bounds for $F(n)/n^{3/2}$ on distinct subsequences.  Different upper-bound
constructions are irrelevant unless accompanied by order-specific lower
bounds applying to every signing.

## 6. Proof and falsification standards

For each claimed result:

1. State the theorem precisely, including all quantifiers and scaling.
2. Give a complete proof or label the statement conjectural.
3. Check $n=2,3$, parity, and the absolute value.
4. Identify every use of compactness, convexity, minimax, concentration, or
   limit exchange and quote the required hypotheses.
5. Calculate whether every error is truly $o(n^{3/2})$.
6. Search for finite counterexamples before promoting a lemma.
7. Distinguish exact arithmetic, trusted-solver certificates, floating-point
   evidence, and proof.
8. If using a random construction, distinguish existence, high probability,
   almost sure statements, and expectation.
9. If using an SDP or spectral bound, do not assume the spectral optimizer is
   Boolean.
10. Preserve failed attempts with the exact line at which they fail.

Do not claim the MathOverflow problem is solved unless every step is complete
and independently auditable.

## 7. Requested output

Return a research report with these sections:

1. **Independent audit:** any error found in the baseline or second-attempt
   results.
2. **Strongest new theorem:** precise statement and complete proof.
3. **Main settling attempt:** full derivation toward convergence or
   nonconvergence.
4. **Counterexample search:** tempting statements tested and exact failures.
5. **Remaining wall:** the narrowest missing lemma, with exact scaling.
6. **Next experiment:** only if it can distinguish live mathematical
   hypotheses; include reproducible code or pseudocode and certificate
   boundaries.
7. **Status sentence:** one of “proved convergence,” “proved nonconvergence,”
   or an explicitly non-claiming partial result.

Continue until you obtain either a complete solution, a substantial theorem
that changes the frontier, or a sharply isolated settling lemma together with
failed counterexample searches.  “The problem is difficult” is not a valid
stopping point.

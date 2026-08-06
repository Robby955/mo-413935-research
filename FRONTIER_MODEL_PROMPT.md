# Broad frontier-model research prompt: MathOverflow 413935

> **Archive notice.** This prompt preserves the full exploration frontier.
> New research runs should use
> [`COMPOSITION_FRONTIER_PROMPT.md`](COMPOSITION_FRONTIER_PROMPT.md), which
> restricts work to the labeled relative-gauge composition problem and its
> power-saving cross-order target.

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

There is also an asymptotically exact formulation as a fixed-density graph
cut-discrepancy problem.  Put $m=\binom n2$, $e_0=\lfloor m/2\rfloor$, and

\[
 H_n=\min_{\substack{G\text{ on }[n]\\e(G)=e_0}}
 \max_{S\subseteq[n]}
 \left|e_G(S,S^c)-\frac12|S||S^c|\right|.
\]

Then the following explicit finite bound is proved:

\[
 \boxed{
 F(n)-1\le4H_n\le F(n)+\sqrt{\binom n2}+2.}
\]

For the reduction, encode the negative entries of a signing by a graph $G$,
write $t=\sum_{i<j}a_{ij}$ and $k=|S||S^c|$, and observe

\[
 Q_A(x_S)=t-2k+4e_G(S,S^c).
\]

Every switching class has a representative with $|t|\le\sqrt m$, because
$\mathbb E_zQ_A(z)^2=m$.  Editing at most $(\sqrt m+1)/2$ edge signs reaches
the layer $e(G)=e_0$, at a cost at most $\sqrt m+1$ in $M(A)$.  Thus
existence of the original limit is equivalent to existence of
$\lim 4H_n/n^{3/2}$, with the same value.

Backurs and Bavarian define the all-cuts absolute parameter

\[
 D_p(G)=\max_S|e_G(S,S^c)-p|S||S^c||
\]

and prove

\[
 D_{\rho_G}(G)
 =\Omega(\min(\rho_G,1-\rho_G)n^{3/2}).
\]

This applies to the layer defining $H_n$, since
$\rho_G=1/2+O(n^{-2})$ and $D_{\rho_G}=D_{1/2}+O(1)$.  It proves only the
correct order of magnitude with an unspecified nonsharp constant.  It gives
neither the candidate $1/8$ lower constant for $H_n$, nor convergence, nor a
cross-order inequality.  Use A. Backurs and M. Bavarian, “On the sum of L1
influences,” ECCC TR13-039, revision 2 (2014), Definition 7.1 and Theorem 7.2, as the
direct source.  Do not conflate this parameter with the one-sided maximum
over balanced complementary halves called bipartite discrepancy by Erdős,
Goldberg, Pach, and Spencer.

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
\begin{array}{c|rrrrrrrrrrrrr}
n&2&3&4&5&6&7&8&9&10&11&12&13&14\\ \hline
F(n)&1&3&4&4&5&9&10&12&13&17&18&20&21.
\end{array}
\]

The order-11 lower certificate scans and hashes all 12,005,168 unlabeled
residual graphs, exactly filters 2,153,606 eligible classes for cut
evaluation, and reproduces the result by generating that reduced stream with
the same evaluator.  Deterministic samples use a separate adjacency formula.  Its
completeness trusts nauty; order 12 then follows from puncturing, parity, and
an explicit witness.  At order 13, a fail-closed eight-shard scan hashes all
1,018,997,864 eleven-vertex residual graphs.  Exactly two rooted records have
order-12 maximum below 20, and direct incident-column enumeration gives
extension minimum 24 for both.  The Bellman identity and a principal block of
the order-14 Paley matrix give $F(13)=20$; heredity and parity then give
$F(14)=21$.  The normalized values are not monotone.  At order 10, a Paley conference
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

### 3.9 Exact finite optimizer-composition obstruction

For fixed signings $A,B$, define the optimized cross completion

\[
 J(A,B)=\min_C
 M\begin{pmatrix}A&C\\C^{\mathsf T}&B\end{pmatrix},
\]

and let $E(B)=J([0],B)$ be the best one-vertex extension. Complete exact
enumeration of optimal switching-permutation classes gives

\[
\begin{array}{c|c|c}
n&\#\text{ optimum classes}&\text{distribution of }E(B)\\ \hline
7&6&10^{\times4},\ 12^{\times2}\\
8&2&12^{\times2}\\
9&15&13^{\times4},\ 15^{\times11}\\
10&2&17^{\times1},\ 19^{\times1}.
\end{array}
\]

The later complete catalogue certificates give

\[
F(11)=17,\qquad F(12)=18,\qquad F(13)=20,\qquad F(14)=21.
\]

The last two values use the complete order-12 residual scan and the Paley
$C_{14}$ witnesses; they do not establish an asymptotic pattern.

Thus optimum signings with the same scalar value can have different future
extension costs. More sharply,

\[
 \boxed{
 \min_{M(A)=F(2),\ M(B)=F(8)}J(A,B)=15>F(10)=13.}
\]

No optimum order-10 signing contains an optimum order-8 principal submatrix.
The smallest internal sacrifice repairs the obstruction: if

\[
 K_{s,k}(u,v)=
 \min_{M(A)\le u,\ M(B)\le v}J(A,B),
\]

then exact enumeration gives

\[
 K_{2,8}(1,10)=15,\qquad K_{2,8}(1,12)=13.
\]

The critical $(2,8)$ result directly enumerates every cross signing for all
optimum order-8 representatives. A deterministic order-8 graph6 witness with
maximum 12 and a cross completion of maximum 13 proves sufficiency of the
two-unit slack. The full computation and its dependency boundary are in
`verification/research_cross_block_composition.py`.

This kills exact-optimizer heredity, not asymptotic composition with
subleading internal slack. A live ground-state target is to prove, for some
nonnegative $\eta_n=o(n^{3/2})$ and $\delta>0$,

\[
 K_{n,m}(F(n)+\eta_n,F(m)+\eta_m)^{2/3}
 \le F(n)^{2/3}+F(m)^{2/3}
 +O((n+m)^{1-\delta}).
\]

### 3.10 Exact weighted covering-radius Bellman state

Let

\[
 \mathcal P_n=\{\pm1\}^n/\{x\sim-x\},\qquad
 d_\pm([u],[v])=\min\{d_H(u,v),d_H(u,-v)\}.
\]

For an order-$n$ signing $B$, put $M=M(B)$ and

\[
 w_B([x])=\frac{M-|Q_B(x)|}{2},
\qquad
 \rho_{\rm w}(B)=\max_{[b]\in\mathcal P_n}
 \min_{[x]\in\mathcal P_n}
 \bigl(d_\pm([b],[x])+w_B([x])\bigr).
\]

Pairing the two values of a new vertex spin proves the exact identities

\[
 E(B)=M(B)+n-2\rho_{\rm w}(B)
\]

and, with
$\delta_{\rm w}(B)=\lfloor n/2\rfloor-\rho_{\rm w}(B)$,

\[
 F(n+1)=\min_B
 \left(M(B)+(n\bmod2)+2\delta_{\rm w}(B)\right).
\]

This is an exact one-vertex Bellman state. It is not enough to retain only
the projective code of configurations satisfying $|Q_B|=M(B)$. All six
optimal order-7 classes have extremizer-code covering radius 3, but lower
energy layers reduce the weighted radius to 2 in exactly the two classes with
extension value 12. At order 9, by contrast, radii 4 and 3 exactly distinguish
the extension values 13 and 15; the two order-10 radii 3 and 2 give extension
values 17 and 19.

More sharply, the optimal order-9 graph6 records `G?qmaw` and `GCpbaw` have
the identical projective absolute-energy histogram

\[
 \#\{|Q_B|=0,4,8,12\}=(60,111,60,25),
\]

so both scalar partition functions agree for every temperature, yet their
extension values are 13 and 15. A scalar free-energy curve is therefore not a
closed Bellman state. This does not disprove a nonlocal minimax interpolation;
it requires that interpolation to retain geometry or a richer order parameter.

The weighted radius depends only on configurations in the exact window

\[
 |Q_B(x)|\ge M(B)-2\lfloor n/2\rfloor;
\]

deeper energy layers can never attain its inner minimum. Complete enumeration
of all signings through residual order 8 gives the nontrivial Pareto frontiers

\[
 \mathcal B_6=\{(5,2),(7,1),(9,0)\},\qquad
 \mathcal B_8=\{(10,1),(12,0)\},
\]

for $(M,\delta_{\rm w})$. Every point of the first frontier gives $F(7)=9$,
and both points of the second give $F(9)=12$. Internal energy and covering
deficit can therefore compensate exactly.

A live route is to prove an asymptotic law or stability theorem for the
Pareto profile $(M(B),\delta_{\rm w}(B))$ over the near-optimal energy window.
Any proposed compression of this state must be tested against the order-7
counterexample.

### 3.11 Nonlinear Gaussian stability

The full square-covariance calculation is now known. For $s\ge0$, let

\[
 D_s=1+s^2(n-1),\quad d_s=\frac{2s}{D_s},\quad
 z_{ij,s}=\frac{s^2|(A^2)_{ij}|}{D_s}.
\]

Then

\[
 M(A)\ge\frac1\pi\sum_{i<j}
 [\arcsin(z_{ij,s}+d_s)-\arcsin(z_{ij,s}-d_s)].
\]

At $s=1/\sqrt{n-1}$, convexity of the arcsine difference gives

\[
 M(A)\ge\frac{n(n-1)}\pi\arcsin\frac1{\sqrt{n-1}}
 +\frac{\|A^2-(n-1)I\|_F^2}
 {8\pi(n-1)(n-2)^{3/2}}.
\]

The parity defect theorem gives

\[
 \|A^2-(n-1)I\|_F^2\ge
 \begin{cases}
 n(n-1),&n\text{ odd},\\
 2n(n-2),&n\equiv0\pmod4,\\
 0,&n\equiv2\pmod4,
 \end{cases}
\]

and consequently $F(20)\ge30$ and $F(21)\ge32$.  A sequence attaining the
asymptotic lower constant $1/\pi$ must satisfy
$\|A_n^2-(n-1)I\|_F=o(n^2)$, but this is only vanishing normalized Gram
defect; uniform random signings satisfy it in probability.

This entire one-parameter route has a sharp barrier. If $G_n$ is the minimax
value of the full displayed arcsine certificate after optimizing $s$, then

\[
 \lim_nG_n/n^{3/2}=1/\pi.
\]

Paley conference principal blocks of order $n+o(n)$ prove the upper bound
uniformly in $s$. Do not retry parameter optimization or ordinary convex
combinations of these certificates as a route to a larger universal constant.

### 3.12 Exact multivertex Bellman state

For fixed internal blocks $B,D$ of orders $n,k$, define

\[
 J(B,D)=\min_C M\begin{pmatrix}B&C\\C^{\mathsf T}&D\end{pmatrix},\qquad
 L(B,D)=\max_{x,y}|Q_B(x)+Q_D(y)|.
\]

Let $\mathcal R_{n,k}=\{[xy^{\mathsf T}]\}$ be the projective rank-one code,
with weight

\[
 w_{B,D}([xy^{\mathsf T}])
 =\frac{L(B,D)-|Q_B(x)+Q_D(y)|}{2}.
\]

Its weighted covering radius satisfies the exact identity

\[
 J(B,D)=L(B,D)+nk-2\rho_{\mathrm w}^{\square}(B,D).
\]

Equivalently,

\[
 F(n+k)=\min_{B,D}
 [L(B,D)+(nk\bmod2)+2\delta_{\mathrm w}^{\square}(B,D)].
\]

Only the window
$|Q_B(x)+Q_D(y)|\ge L(B,D)-2\lfloor nk/2\rfloor$ can affect the state.
This is an exact reformulation, not yet a composition estimate.

### 3.13 Weighted entropy bound and density-one control

Define

\[
 \Xi(B,D)=\log\sum_{R\in\mathcal R_{n,k}}
 \exp\left(-\frac{2w_{B,D}(R)^2}{nk}\right).
\]

A weighted Hamming-ball union bound proves

\[
 J(B,D)\le L(B,D)+\sqrt{2nk(\Xi(B,D)+\log4)}.
\]

The bound is correct, but the originally proposed use of it is impossible.
For every $B,D$,

\[
 L(B,D)+\sqrt{2nk(\Xi(B,D)+\log4)}
 \ge\sqrt{2nk(n+k)\log2}.
\]

For balanced blocks the left side has constant at least
$2\sqrt{\log2}$, while the required target has asymptotic ceiling $\sqrt2$.
The leading gap is $0.250895\ldots$. Retaining the discarded cross term in
the same all-state Hoeffding union sum has the same obstruction.

The live replacement retains dependence. Regard
$G=\mathcal R_{n,k}$ as a group, fix a cross seed $C_0$, and put
$c(U)=|\langle C_0,U\rangle|$ and
$h(R)=|Q_B(x)+Q_D(y)|$. At target $K$, the bad switching shifts are exactly

\[
 \mathcal B_K=\bigcup_{R\in G}R\{U:c(U)>K-h(R)\}.
\]

Prove $\mathcal B_K\ne G$ at the power-saving composition threshold. This is
a weighted sumset or max-plus convolution problem, not an entropy first
moment.

For one-vertex Bellman-optimal predecessors, exact telescoping also proves
that for every $g(N)\to\infty$, all but $O(N/g(N))$ orders in $[N,2N]$ have
internal slack at most $g(N)\sqrt N$ and weighted deficit at most
$\frac12g(N)\sqrt N$. An explicit abstract countermodel still oscillates
with $O(\sqrt n)$ scalar Bellman cost at every order, so sparse exceptions are
not the decisive wall. The exact diagnostic is stabilization of

\[
 N^{-3/2}\sum_{n=N}^{2N-1}
 [F(n+1)-F(n)-(n\bmod2)].
\]

### 3.14 Complete order-nine state test

Exhaustion of all 12,346 root-normalized order-nine signings gives

\[
 \mathcal B_9=\{(12,0)\}.
\]

The records `GHOgmo` and `Gxd?Dc` both have $M=14$, projective energy
histogram $(124,85,37,10)$ at energies $(2,6,10,14)$, and exact-maximizer
ordered pair-distance law $(10,16,14,20,40)$. Their profiles
$(\rho_{\rm ext},\rho_{\rm w},E)$ are nevertheless $(4,4,15)$ and
$(3,3,17)$. Thus the energy histogram plus the two-point geometry of exact
maximizers is still not a closed Bellman state.

The complete energy-coloured two-point distribution separates all weighted
deficits at order 9. This is finite evidence only. Test any proposed finite
overlap-state closure against higher orders before treating it as an
asymptotic theorem.

### 3.15 Complete order-ten temperature test

The complete rooted order-ten catalogue has 274,668 records and 6,012
absolute-energy histograms. The minimax partition function has exactly three
positive-temperature histogram phases. Their maxima are $15,15,13$, and the
transitions occur at

\[
 t=0.658478948\ldots,\qquad t=0.792460762\ldots.
\]

The high-temperature phase is conference, the middle phase is not a ground
state, and the low-temperature phase attains $F(10)=13$. At $t=0$ all
signings tie. Therefore the outer optimizer genuinely changes with
temperature. Any interpolation must allow optimizer phase changes rather than
track a single signing across all $\beta$.

### 3.16 Exact negative-replica recursion and its wall

For uniform order-$n$ disorder define

\[
 f_{n,t}(A)=\frac{Z_A(t)}{2^{n+1}(\cosh t)^{\binom n2}},\qquad
 \mathcal G_n(q,t)=\log\mathbb E_A f_{n,t}(A)^{-q}.
\]

Jensen over the random cross block, followed by pairing one internal signing
with its global negation, proves exactly

\[
 \mathcal G_{n+k}(q,t)\ge\mathcal G_n(q,t)+\mathcal G_k(q,t)
\]

for all $q>0$ and real $t$. No absolute value or normalization factor is
lost. With $u=\tanh t$, the signed Eulerian expansion also gives
$f_{n,u}=T_{u/v}f_{n,v}$.

Sharp Boolean reverse hypercontractivity implies

\[
 \frac{\mathcal G_n(Q',u')}{Q'}
 \le\frac{\mathcal G_n(Q,u)}Q
\]

when $(1+Q')u'^2=(1+Q)u^2$. This is the exact extensive block-parameter
curve, but the inequality points in the wrong direction. For

\[
 t_r=\beta/\sqrt r,\qquad
 q_r=\theta/\tanh^2t_r-1,
\]

the first proposed repair was the power-saving transport estimate

\[
0\le \frac{q_N}{q_n}\mathcal G_n(q_n,t_n)
-\mathcal G_n(q_N,t_N)\le C N^{2-\delta}
\]

on balanced child sizes. It is false. If

\[
 \Delta_n=\frac{q_{2n}}{q_n}\mathcal G_n(q_n,t_n)
 -\mathcal G_n(q_{2n},t_{2n}),
\]

then the orbit-counted finite-disorder Laplace bound gives

\[
 \liminf_n\frac{\Delta_n}{n^2}
 \ge\frac{2\theta}{\beta^2}
 \left[\frac{\beta^2}{8}
 -\frac\beta2\left(1-\frac1{\sqrt2}\right)-\log2\right]-\log2.
\]

At $(\beta,\theta)=(4,8)$ this is
$\sqrt2-2\log2=0.027919\ldots>0$. Thus scalar parameter transport is closed
at leading order.

The transport defect has an exact entropy-production representation. For a
positive function on an $m$-cube, along the invariant noise curve with
$m\ge1$, $q_s>0$, and nonnegative semigroup time,

\[
 \Delta=q_N\int R_s\,ds,
 \qquad
 R_s\ge\frac{2D(\mu_s\Vert\mathrm{Unif})^2}{3q_sm},
\]

where $\mu_s$ is the negative-moment escort law. Equality holds only for a
constant function. This quantifies non-saturation but does not restore the
block geometry discarded by scalar transport.

### 3.17 Exact conditional-alignment state

For a split $N=n+k$, set

\[
 C'_{n,k}=D_n\oplus D_k\oplus R_{n,k},
\]

where $R_{n,k}$ is the rectangular row-and-column switching code. For
$n,k\ge3$, $D_N\subset C'_{n,k}$ and every coarse coset has $2^{N-1}$ fine
cosets. If $r_h$ is the conditional fine density on a coarse fiber and
$K_q(h)=\mathbb E_{\rm fiber}r_h^{-q}$, then exactly

\[
 \mathcal G_N=\mathcal G_n+\mathcal G_k
 +\mathcal G_{n,k}^{\square}+\mathcal A_{n,k},
 \qquad
 \mathcal A_{n,k}=\log\mathbb E_{\widehat U_q}K_q(h)\ge0.
\]

A fiber has the concrete relative-switching states

\[
 Y_{\alpha,\beta,\tau}=
 \begin{pmatrix}
 D_\alpha AD_\alpha&C\\
 C^{\mathsf T}&\tau D_\beta BD_\beta
 \end{pmatrix},
\]

and, under the augmented state-sum convention,

\[
 r(\alpha,\beta,\tau)
 =\frac{2^{N+1}Z_Y}{Z_AZ_BZ_C^\square},
 \qquad\mathbb E_{\rm fiber}r=1.
\]

Sequentially taking $q\to\infty$ and then $t\to\infty$ shows that the
alignment term has slope

\[
 F(n)+F(k)+B_\square(n,k)-F(N).
\]

For balanced blocks its normalized liminf is at least

\[
 \frac2\pi+\sqrt{\frac2\pi}-\sqrt2
 =0.0202907708\ldots>0.
\]

Alignment is therefore leading order. The rectangular marginal alone cannot
repair scalar transport.

The first coefficient of $\log r$ is the mixed four-cycle Hamiltonian

\[
\mathcal H_4=\frac12\left[
 \operatorname{tr}(P^2CC^{\mathsf T})
 +\operatorname{tr}(Q^2C^{\mathsf T}C)
 +\tau\operatorname{tr}(PCQC^{\mathsf T})
 -nk(N-2)\right],
\]

where $P=D_\alpha AD_\alpha$ and $Q=D_\beta BD_\beta$. Exact split-$2+4$
examples have identical complete local graph and rectangular scalar
partition curves but different conditional moments

\[
 K_2=1.4713253466\ldots,\qquad K_2=1.1351460354\ldots.
\]

Here $\dim D_2=1$, so the relative-gauge parameterization covers the true
fiber twice; normalized moments are unchanged. Their $\mathcal H_4$ laws
already differ. Thus no collection of the three local scalar free-energy
curves determines alignment.

A proposed alignment-transport estimate is not yet progress: after defining
$a_j=\mathcal G_j(q_j,t_j)/q_j$, its left side is identically
$a_N-a_n-a_k$. The actual target is a variational limit or power-saving
composition theorem derived from a controlled relative-switching Hamiltonian
hierarchy. Higher connected Eulerian terms cannot be discarded without a
uniform mean-field cluster bound.

### 3.18 Microcanonical relative-switching profile theorem

There is a balanced homomorphism from the two augmented graph state spaces
and the projective rectangular rank-one state space onto the relative-gauge
group $G_{n,k}$ of size $2^{n+k-1}$. For local deficits $d_A,d_B,d_C$ and
$L=M(A)+M(B)+B_C$, it gives the exact max-plus identity

\[
 L-M(Y_g)=\min_{\pi(a,b,r)=g}[d_A(a)+d_B(b)+d_C(r)].
\]

The convention is exact: the graph factors are
$\{\pm1\}\times\mathbb P_n$ and $\{\pm1\}\times\mathbb P_k$, while the
rectangular factor is the projective absolute rank-one space identifying
$R$ with $-R$.  A different one-sided theorem uses projective $|Q|$ graph
states and signed full-spin cross pairs.  Its balanced map is

\[
 ([z],[w],x,y)\mapsto([zx],[wy],x_1),
\]

and its fiber maximum dominates, but need not equal, $M(Y_g)$.  Its
$2^{N-1}$-st deficit order statistic is therefore also a valid composition
gain.  On the standard balanced $C_{14}$ split it gives $\widetilde\Lambda=8$
and 596440 target states; the exact theorem gives $\Lambda=10$, true gain 22,
and 304908 target states.  Do not mix the two identities or call the alternate
numbers a verification of the exact max-plus theorem.

The alternate raw statistic has a universal balanced floor.  With
$\mu_r=\mathbb E|\sum_{i=1}^r\varepsilon_i|$,

\[
 L-\widetilde\Lambda
 \ge\max\{M(A),M(B)\}+r\mu_r-3r.
\]

It follows that this statistic can prove the balanced power-saving target only
if

\[
 \liminf F(r)/r^{3/2}\ge
 \sqrt{2/\pi}/(2^{3/2}-1)=0.4363775564\ldots.
\]

This does not kill the route under the value-$1/2$ hypothesis.  The optimized
canonical exponential bounds for the alternate and exact conventions differ
by at most $2L/(N-1)=O(\sqrt N)$ when $L=O(N^{3/2})$, so the same conditional
barrier applies to either canonical cgf.  It does not apply to the raw exact
augmented order statistic.

Let $\Lambda(A,B,C)$ be the $2^{n+k-1}$-st smallest product-triple deficit,
with multiplicity. Then

\[
 F(n+k)\le M(A)+M(B)+B_C-\Lambda(A,B,C).
\]

This is a deterministic composition theorem, not the failed weighted Hamming
union bound. It is optimal if one discards all additive labels and retains
only the three local deficit histograms. The exact $2+4$ collision has
$\Lambda=0$ in both cases but true gains 4 and 2, so group geometry can still
improve the scalar guarantee.

The scalar large-deviation criterion is exact and should be used in its full
form.  Put the uniform measure on the three local state spaces, write
$e_A=M(A)-d_A$, $e_B=M(B)-d_B$, $e_C=B_C-d_C$, and let

\[
 X=e_A+e_B+e_C
\]

for three independent uniform local states.  The product space has
$2^{2N-2}$ points and each relative-gauge fiber has $2^{N-1}$ points.
Exponential counting therefore proves, for every $t>0$,

\[
 \boxed{
 \min_gM(Y_g)
 \le
 \inf_{t>0}
 \frac{(N-1)\log2+\log\mathbb E e^{tX}}{t}.}
\]

Equivalently, with

\[
 \kappa_N(\beta)=\frac1N\log\mathbb E
 \exp\left(\frac{\beta X}{\sqrt N}\right),
\]

a target $\theta_NN^{3/2}$ is certified whenever

\[
 \kappa_N(\beta)-\beta\theta_N
 \le-\left(1-\frac1N\right)\log2.
\]

This is the canonical exponential relaxation of the exact microcanonical
order statistic, not the failed weighted Hamming-ball union bound. The
unknown object is the complete leading-scale cgf or a matching far-upper-tail
large-deviation theorem at probability $2^{-N+o(N)}$. Passing from a limiting
cgf to an exact rate requires the relevant uniform
Gärtner--Ellis/exponential-tightness hypotheses; do not call the two objects
equivalent without proving them.

Two scalar concentration shortcuts are closed, with limited scope:

1. For balanced blocks $n=k=r$, any global centered quadratic majorant

   \[
    \log\mathbb E e^{t(X-\mathbb EX)}\le v_rt^2/2
   \]

   must have $v_r\ge r(r-1)$.  Optimizing the cgf criterion with such a
   majorant has normalized certificate floor

   \[
    \sqrt{\frac{\log2}{2}}=0.588705011\ldots>\frac12.
   \]

   This rules out Gaussian, subgaussian, or variance-only cgf
   certificates.  It is not a lower bound on the true profile quantile:
   higher cumulants and the non-Gaussian far tail remain live.

2. Ordinary degree-two Bonami hypercontractivity gives

   \[
    \|Q_A\|_p\le(p-1)\sqrt{\binom n2}
   \]

   and only an $\exp[-\Theta(\sqrt n)]$ upper tail at threshold
   $cn^{3/2}$.  The fiber criterion needs $\exp[-\Theta(n)]$.  This rules
   out the standard bounded-moment certificate, not signing-sensitive
   higher-order inequalities or a full large-deviation theorem.

The sharp remaining unlabeled scalar statement can be written using exact
shell counts

\[
 H_A(d)=\#\{a:d_A(a)=d\},\quad
 H_B(d)=\#\{b:d_B(b)=d\},\quad
 H_C(d)=\#\{r:d_C(r)=d\}.
\]

Define

\[
 \mathcal S_N(s)=\frac1N
 \max_{d_A+d_B+d_C<sN^{3/2}}
 [\log H_A(d_A)+\log H_B(d_B)+\log H_C(d_C)].
\]

There are only polynomially many supported shell triples, so

\[
 \mathcal S_N(s)
 \le\frac1N\log\#\{(a,b,r):
 d_A(a)+d_B(b)+d_C(r)<sN^{3/2}\}
 \le\mathcal S_N(s)+\frac{6\log(N+1)}N,
\]

where the middle count is over local state triples, not merely deficit
values.  If

\[
 T_N=[F(n)^{2/3}+F(k)^{2/3}]^{3/2}+R_N,
 \qquad |R_N|\le C N^{3/2-\delta},
\]

with $C,\delta>0$ uniform over the proportional splits, and

\[
 s_N=\frac{M(A)+M(B)+B_C-T_N}{N^{3/2}},
\]

then the sufficient scalar shell lemma is

\[
 \boxed{
 \mathcal S_N(s_N)+\frac{6\log(N+1)}N
 <\left(1-\frac1N\right)\log2.}
\]

It must hold uniformly over the proportional splits and a composable
near-optimal choice of the three local signings.  Pointwise convergence of
shell entropies without a strict or power-saving margin is insufficient.

Exact balanced-split calibration has been carried out for every balanced
vertex split of the banked optimum witnesses at orders 12, 13, and 14.  The
unlabeled scalar profile bounds are respectively

\[
 24\ldots26,\qquad28\ldots30,\qquad31\ldots33,
\]

while the true full-signing values are $18,20,21$; the observed labeled
geometry bonuses are $6\ldots8$, $8\ldots10$, and $10\ldots12$.  At the
first admissible lattice energies above the near-subadditive targets, the
numbers of product triples are

\[
\begin{array}{c|c|c}
N&\text{target-triple range}&\text{fiber count}\\ \hline
12&213312\ldots269122&2048\\
13&436864\ldots495632&4096\\
14&231581\ldots305465&8192.
\end{array}
\]

Thus the unlabeled pigeonhole condition fails in all these finite
calibrations even though the labeled relative-gauge convolution has a large
gain.  This is evidence for retaining labels, not an asymptotic disproof of
the scalar shell lemma: an $O(N)$ finite geometry gap would still be
power-saving.

There is now an exact labeled improvement.  For a threshold $s$, let $b_s(g)$
count subthreshold triples in gauge fiber $g$, put $K=2^{N-1}$,
$\mu_s=K^{-1}\sum_g b_s(g)$, and let

\[
 V_s=\sum_{\chi\ne1}|\widehat b_s(\chi)|^2
 =K^{-1}\sum_g(b_s(g)-\mu_s)^2.
\]

Then

\[
 \boxed{\min_gb_s(g)\le\mu_s-\sqrt{V_s/(K-1)}.}
\]

If the right side is below one, an empty fiber gives
$\min_gM(Y_g)\le L-s$.  Each $\widehat b_s(I,J,\epsilon)$ is an exact
thresholded convolution of three local labeled shell transforms, so the
bound is computable without enumerating full signings.  It is strictly better
than the scalar profile in the existing $2+4$ collision: at $s=2$ the scalar
shell has exactly $K=32$ triples, but the occupancy law
$0^8,1^{16},2^8$ gives $V=1/2$ and certifies gain two.

The generic variance inequality is not closure.  On the standard balanced
$C_{14}$ split at target energy 27, $\mu=304908/8192$, there are 8159
nonzero nontrivial coefficients, $V=635307383/4194304$, and the true
occupancy range is 0 through 87, with exactly one empty fiber among 8192.
The bound gives $37.0842\ldots>1$ and does not detect it.  The sharp next target is a low-tail or
higher-moment inequality for $b_s$ derived from its factored Fourier
coefficients.

There is an exact moment hierarchy, so this target is no longer vague.  For
finite $A\subset\mathbb Z_{\ge1}$,

\[
 P_A(x)=(1-x)\prod_{a\in A}\frac{(x-a)(x-a-1)}{a(a+1)}
\]

satisfies $P_A(0)=1$ and $P_A(j)\le0$ for every positive integer $j$.
Therefore $\Pr[b_s=0]\ge\mathbb E P_A(b_s)$.  All required moments are
zero-sum products of the shell-factorized Fourier coefficients.  At the
balanced $C_{14}$ target, the nine roots

\[
 A=\{9,17,26,36,46,56,67,76,86\}
\]

give a positive exact degree-19 certificate.  Equivalently, the order-nine
localizing matrix $(m_{i+j+1}-m_{i+j})$ is not positive semidefinite.

The scalable wall is sharper than bounded order.  With $z$ empty fibers, a
localizing polynomial normalized by $q(0)=1$ has negative margin at most
$z/K$, and filling the holes by occupancy one changes every localizing matrix
by exactly $(z/K)e_0e_0^{\mathsf T}$.  The unique $C_{14}$ hole therefore has
an exponentially small normalized signal.  Its exact degree-nine witness has
magnitude only $7.6390\ldots\,10^{-6}$.

A box-principle collision gives a vacant and nonvacant $K$-point occupancy
pair with identical moments through degree
$\Omega(\sqrt{K/\log K})$.  For $K=2^{N-1}$ this is exponential in $N$.
These occupancies are abstract, not claimed realizable signing shells.  The
other generic conversions fail comparably at the calibration: the
Fourier--PSD hierarchy needs at least 6827 of 8192 characters; canonical
vacancy detection needs occupancy inverse temperature $\log K=\Theta(N)$;
and its alternating collision expansion first certifies the $C_{14}$ hole at
degree 87.

Do not ask merely for a degree-$O(N)$ estimate with polynomial absolute
error.  Seek an exact signing-specific character-sum sign, a theorem producing
polynomially many good gauges, or a shell-factorized estimate with exponential
precision justified by structure.

If the scalar shell inequality fails asymptotically, do not return to a
one-point union bound.  For each gauge define the labeled shell convolution

\[
 \mathcal N_g(s)=
 \#\{(a,b,r):\pi(a,b,r)=g,
 d_A(a)+d_B(b)+d_C(r)<sN^{3/2}\}.
\]

The actual goal is to prove $\mathcal N_g(s_N)=0$ for some $g$.  Its group
Fourier coefficients factor into the three local correlation tensors.
Control of this labeled convolution, or a justified hierarchy of its
low-weight characters beginning with the mixed four-cycle layer, is the
geometric fallback.  Any truncation must prove that the discarded characters
cannot refill every fiber.

### 3.19 Square-order Paley correction

For every odd prime power $m$, the full Paley conference matrix of order
$m^2+1$ has a Boolean eigenvector of eigenvalue $m$, obtained from additive
$\mathbb F_m$-cosets in $\mathbb F_{m^2}$. Hence that matrix has exact maximum
$m(m^2+1)/2$. This is the known regular-conference construction, not a new
theorem and not a lower bound on $F(m^2+1)$. At $m=3$, the matrix maximum is
15 while $F(10)=13$. Do not infer a subsequential value of $F$ from this
construction.

### 3.20 Paley Fourier-leakage identity

For an odd prime power $q\equiv1\pmod4$, let $C_{q+1}$ be the symmetric
Paley conference matrix. For Boolean $f:\mathbb F_q\to\{\pm1\}$, let
$S(f)=\sum f$ and partition its unitary additive Fourier energy into the two
nonzero frequency halves $E_+(f),E_-(f)$ selected by the signs of the
quadratic Gauss-sum eigenvalues. Put $W(f)=\min(E_+(f),E_-(f))$. Exact Fourier
diagonalization and optimization of the infinity sign give

\[
 M(C_{q+1})=\frac{(q+1)\sqrt q}{2}
 -\sqrt q\min_f\left[
 W(f)+\frac{(|S(f)|-\sqrt q)^2}{2q}
 \right].
\]

Hence the Paley Boolean maximum approaches its spectral ceiling exactly when
there are sign functions with

\[
 |S(f_q)|=o(q),\qquad W(f_q)=o(q).
\]

There is also an exact coordinate-free formulation for every symmetric
conference matrix, not only the Paley family. If

\[
 C^2=(n-1)I
\]

and $E_\pm$ are the eigenspaces for the eigenvalues
$\pm\sqrt{n-1}$, put

\[
 \alpha_\pm(C)=\frac1{\sqrt n}
 \sup_{\substack{v\in E_\pm\\\|v\|_2=1}}\|v\|_1.
\]

Projection onto the two eigenspaces and the identity

\[
 \max_{x\in\{\pm1\}^n}\|P_\pm x\|_2
 =\sup_{\substack{v\in E_\pm\\\|v\|_2=1}}\|v\|_1
\]

give the exact Boolean-eigenspace formula

\[
 \boxed{
 \frac{2M(C)}{n\sqrt{n-1}}
 =\max\{2\alpha_+(C)^2-1,\ 2\alpha_-(C)^2-1\}.}
\]

Both eigenspaces are required because the objective contains an absolute
value. In particular, a uniform gap
$\alpha_+(C)^2,\alpha_-(C)^2\le1-\eta$ gives the genuine Boolean
improvement

\[
 M(C)\le\left(\frac12-\eta\right)n\sqrt{n-1}.
\]

Conversely, spectral saturation is exactly the existence of an
$\ell_2$-unit vector in at least one eigenspace whose $\ell_1$ norm is
$(1-o(1))\sqrt n$. Merely controlling eigenvalues or empirical spectral
measures cannot decide this.

Exact source-built enumeration gives

\[
 M(C_6)=5,\quad M(C_{14})=21,\quad M(C_{18})=33,\quad
 M(C_{30})=75,
\]

with Boolean-to-spectral ratios
$\sqrt5/3,3/\sqrt{13},11/(3\sqrt{17}),5/\sqrt{29}$. They increase in these
four selected prime-field cases, but no monotonicity or asymptotic fit is
proved. Do not confuse this construction problem with the minimax lower bound.
Even spectral saturation of $C_{q+1}$ would still require the dense-order
minimax rigidity estimate

\[
 F(q+1)\ge(1-o(1))M(C_{q+1})
\]

to force the original limit to be $1/2$.

The Fourier half is no longer open on every sequence. For prime
$p\equiv1\pmod4$, the balanced half-interval function satisfies

\[
 S(f_p)=1,\qquad W(f_p)\le\frac{2p}{\ell(p)-1},
\]

where $\ell(p)$ is the least quadratic nonresidue. Quadratic reciprocity and
the prime number theorem in each fixed progression produce a multiplicatively
dense prime sequence with $\ell(p)\to\infty$. Hence

\[
 \frac{2M(C_{p+1})}{(p+1)\sqrt p}\to1
\]

on that sequence. Do not spend time reproving this leakage construction. The
remaining value-specific wall is the minimax rigidity estimate

\[
 F(p+1)\ge(1-o(1))M(C_{p+1})
\]

on the same sequence. Exact conference optimality at orders 6 and 14 does not
prove it, and order 10 is a finite counterexample to exact optimality.

Be precise about the logical strength of this wall. Let $n_j=p_j+1$ be the
multiplicatively dense aligned sequence just described, so that

\[
 \frac{n_{j+1}}{n_j}\to1,\qquad
 \frac{M(C_{n_j})}{n_j^{3/2}}\to\frac12.
\]

Then

\[
 \boxed{
 F(n_j)\ge(1-o(1))M(C_{n_j})
 \quad\Longleftrightarrow\quad
 \lim_{n\to\infty}\frac{F(n)}{n^{3/2}}=\frac12.}
\]

For the forward implication, monotonicity transfers the subsequence lower
bound to every $n_j\le n<n_{j+1}$, and the global conference upper bound
supplies the matching limsup. The reverse implication follows immediately
from the two limits along $n_j$. Thus dense Paley rigidity is a sharp
reformulation of a proposed value-$1/2$ solution, not an easier reduction.
Any proof must exploit a structural restriction applying to every signing;
additional numerical alignment of the Paley construction cannot establish
it.

## 4. Known failed mechanisms that require a material change

Do not simply repeat any of the following:

1. Independent random cross edges in an ordinary block sum: their uniform
   discrepancy is itself order $(n+k)^{3/2}$.
2. Rank-one clone blow-ups: the main term scales as $k^2F(n)$, with an
   additional leading internal-block cost.
3. Direct annealed partition estimates: the $nk\log\cosh t$ cost is
   macroscopic at the extensive temperature.
4. Treating the conference spectral ceiling as attainable on Boolean vectors:
   explicit orders 6, 14, 18, and 30 disprove this, and order 10 is not even
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
12. Composing arbitrary exact minimizers: optimum extension cost is
    class-sensitive, and the exact $(2,8)$ obstruction forces value 15 rather
    than $F(10)=13$. A viable composition must retain a boundary profile or
    permit controlled internal slack.
13. Retaining only exact maximizing spin configurations in a cavity state:
    two optimal order-7 classes have the same extremizer covering radius as
    the four extendible classes, but their next energy layer raises the exact
    extension value from 10 to 12.
14. Retaining the full scalar partition-function curve as the cavity state:
    two optimal order-9 classes have identical absolute-energy histograms at
    every level but different extension values.
15. Adding only the pair-distance law of exact maximizers to that histogram:
    the explicit order-9 records `GHOgmo` and `Gxd?Dc` still have different
    extension values.
16. Optimizing the parameter in the full nonlinear square-covariance
    certificate, or taking ordinary convex combinations of those
    certificates: its minimax asymptotic constant is exactly $1/\pi$.
17. The constant-matching weighted-entropy target formerly stated after the
    multivertex Bellman identity: a universal leading floor
    $2\sqrt{\log2}$ in balanced blocks exceeds the required $\sqrt2$.
18. Pointwise or density-one scalar Bellman-cost bounds by themselves: an
    explicit nonconvergent scalar model has $O(\sqrt n)$ cost at every order
    while satisfying all currently proved scalar cross-order inequalities.
19. Treating reverse hypercontractivity as an asymptotically sharp transport
    theorem: the proposed power-saving defect is rigorously false with an
    $n^2$-scale gap.
20. Adding the rectangular negative moment while discarding conditional
    alignment: the rectangular marginal cannot pay the scalar transport gap
    in the large-parameter regime, while alignment is provably leading.
21. Treating all three local scalar partition curves as a closed block state:
    the exact split-$2+4$ collision, with its twofold gauge redundancy
    accounted for, has different conditional alignment.
22. Calling the formal alignment-transport inequality a reduction: it is
    identically the desired almost-superadditivity in new notation.
23. Inferring a minimax lower bound from Boolean ceiling attainment by the
    square-order Paley matrices: already $F(10)=13<15$ for the first case.
24. Replacing the exact profile cgf by a global quadratic, Gaussian, or
    variance-only majorant: in balanced blocks its best possible certificate
    has normalized floor $\sqrt{\log2/2}>1/2$.
25. Applying ordinary degree-two Bonami hypercontractivity at the target
    threshold: it supplies only an $\exp[-\Theta(\sqrt n)]$ tail, whereas
    fiber noncoverage requires $\exp[-\Theta(n)]$. These two no-go statements
    do not rule out the true cgf, a signing-sensitive large-deviation theorem,
    or labeled cancellation.
26. Treating finite maximizer thinness, the order-$12$ through order-$14$
    balanced profile gaps, or a few conference alignment ratios as an
    asymptotic theorem. They are calibration data and counterexample tests,
    not a uniform shell LDP.
27. Invoking the Backurs--Bavarian fixed-density cut-discrepancy lower bound
    as a convergence theorem: it has a nonsharp constant and contains no
    cross-order comparison.
28. Calling dense Paley minimax rigidity a reduction of the problem: on the
    multiplicatively dense aligned Paley sequence already constructed, that
    rigidity statement is equivalent to the full limit being $1/2$.
29. Treating the one-sided swapped-profile theorem as the exact relative-gauge
    identity.  Projective $|Q|$ graph profiles plus signed full-spin cross
    pairs do have equal fibers and a valid domination bound, but their fiber
    maxima need not equal $M(Y_g)$.  Keep its order statistic distinct from
    the augmented-graph/projective-absolute-cross theorem.
30. Treating nonzero labeled Fourier variance as sufficient at dense target
    occupancy.  The Parseval bound is strict at the scalar boundary, but its
    right side is still $37.0842\ldots$ in the calibrated $C_{14}$ shell even
    though an empty fiber exists.  Higher moments or a genuine lower-tail
    theorem are required.
31. Assuming a subexponential occupancy-moment or small-character hierarchy is
    universally complete.  Abstract vacant/nonvacant laws can agree through
    $\Omega(\sqrt{K/\log K})=2^{\Omega(N)}$ moments, and the calibrated
    Fourier--PSD witness needs 6827 of 8192 characters.  An actual-shell proof
    needs exact realizable-shell structure.
32. Declaring the swapped raw profile either universally dead or universally
    sufficient.  Its exact floor gives the conditional threshold
    $0.436377\ldots$: it is obstructed below that constant and remains viable
    at the current value-$1/2$ hypothesis.
33. Returning to complete/empty equal-cloud amplification.  Fixed-half
    cloud-union cuts retain $k^2d(G)$, and Seidel switching plus $O(N)$ density
    repair changes this only by $O(N)$.  Orthogonal blocks evade that specific
    quotient bound, but the uniform order-four symmetric-Hadamard test already
    has sharp maximum 44 rather than the lossless value 32.
34. Inferring an iterable amplification theorem from either exact order-16
    Hadamard lift.  The first has $M=32=4^{3/2}F(4)$. A sharper fixed-frame
    completion has $M=30$ and proves $F(16)\le30$, but only at one order.
    Its alternating internal types are switching-permutation equivalent and
    both have maximum four, while repeating one literal type has minimum 38.
    Ordinary switching class is therefore not a closed substitution state.
    No uniform map, iteration law, or dense reachability statement is known.

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
- the exact conditional-alignment chain above. Derive a genuine variational
  or composition estimate from the relative-switching law; do not re-propose
  scalar PT, append only the rectangular marginal, or rename the desired
  almost-superadditivity as an alignment lemma.
- the exact microcanonical order statistic above. First test whether the
  three local deficit profiles cross the $2^{n+k-1}$ threshold at the needed
  leading deficit. If they do not, identify the additive Fourier labels that
  supply the missing gain rather than reverting to a one-point union bound.
- a controlled mixed-trace/connected-Eulerian hierarchy beginning with
  $\mathcal H_4$. Any truncation must be uniform when
  $u=\Theta(n^{-1/2})$ and $q=\Theta(n)$.

The first concrete attack is a full microcanonical shell theorem, not another
second-moment estimate. Fix $\varepsilon>0$ and proportional splits
$n/N\in[\varepsilon,1-\varepsilon]$. Permit composable internal signings with
$o(N^{3/2})$ sacrifice. Prove, uniformly over the chosen family, either

\[
 \mathcal S_N(s)=\Sigma_\alpha(s)+o(1),
 \qquad \alpha=\frac nN,
\]

locally uniformly at the required deficit, or directly the equivalent
speed-$N$ far-tail upper bound

\[
 \limsup_{N\to\infty}\frac1N
 \log\mathbb P\{X\ge\theta N^{3/2}\}
 \le-I_\alpha(\theta).
\]

The settling inequality is $I_\alpha(\theta)>\log2$ with a uniform strict
margin, equivalently $\Sigma_\alpha(s)<\log2$ at the near-subadditive
deficit. At a critical equality, supply an explicit power-saving or
dyadically summable second-order term; a bare pointwise LDP does not close
iteration. Determine the component shell rates first, then compute their
infimal convolution. The exact order-$12$ through order-$14$ balanced data
must be used as a regression test: a proposed scalar theorem that predicts
the observed true geometry bonus from only the three unlabeled histograms is
false.

If the scalar threshold fails, attack the labeled convolution itself. For
the relative-gauge group $G$, introduce the energy-refined Fourier
polynomials

\[
 P_A(\chi,z)=\sum_a\chi(\pi_A(a))z^{d_A(a)},\quad
 P_B(\chi,z)=\sum_b\chi(\pi_B(b))z^{d_B(b)},\quad
 P_C(\chi,z)=\sum_r\chi(\pi_C(r))z^{d_C(r)}.
\]

Then Fourier inversion gives the exact truncated coefficient formula

\[
 \mathcal N_g(s)=\frac1{|G|}\sum_{\chi\in\widehat G}\overline{\chi(g)}
 \sum_{d<sN^{3/2}}[z^d]\,
 P_A(\chi,z)P_B(\chi,z)P_C(\chi,z).
\]

Seek a support-noncoverage theorem proving that this is zero for at least one
$g$, even when its trivial-character average exceeds one. Viable mechanisms
include a full character large-deviation principle, an exact association-
scheme sign, an abundance theorem for empty fibers, or an additive-
combinatorial restriction on the three labeled low-deficit sets. A low-weight
character or low-moment expansion is admissible only with a uniform tail bound
at the exponentially small hole-density scale. The mixed four-cycle
Hamiltonian is the first mandatory term, not a complete state.

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

The smallest orthogonal calibration is now mixed rather than negative.  A
fixed-half order-16 signing made from four clouds of size four has six
zero-sum Hadamard cross blocks and exact maximum

\[
 32=4^{3/2}F(4).
\]

After reordering and signed-permutation gauges, the cross blocks share one
nonsymmetric oriented Hadamard and their edge signs form an $F(4)$ optimizer.
The internal block maxima are $6,6,4,4$, and no common-Hadamard gauge makes
all four internal blocks common even up to sign.  This proves a lossless-scale
one-step lift is locally possible, but also says the live state is multi-type:
try a Pareto/substitution recursion for the internal block types.  Do not
claim iteration until the type system closes with summable normalized error.
Even then, supply dense reachability or a compatible truncation/composition
argument; one fixed geometric subsequence does not force the full limit.

A sharper construction pins

\[
H=\begin{pmatrix}
1&1&1&1\\1&1&-1&-1\\1&-1&1&-1\\-1&1&1&-1
\end{pmatrix}
\]

on every upper cross block, with base signs $(+,+,+,-,+,-)$, and alternates
two order-four internal types $P,R$. It has exact maximum 30, so
$F(16)\le30$. A six-state proof shows 30 is the minimum over all four internal
completions in that fixed frame. Both $P,R$ have maximum four and are
switching-permutation equivalent, but their eight-entry framed response
vectors differ; the common-literal-internal subfamily has minimum 38. The
next legitimate finite-type question is therefore whether a bounded or
controlled family of framed response vectors closes under substitution with
a summable defect. A recursion using only $M(D)$, spectrum, or unframed
switching class is already falsified at this first step.

A particularly sharp scalar target is to prove, for some $\delta>0$,

\[
 F(n+m)^{2/3}
 \le F(n)^{2/3}+F(m)^{2/3}+O((n+m)^{1-\delta}).
\]

A dyadically summable relative defect is essential for this argument; a power
saving is a clean sufficient condition. Calculate the error under repeated
balanced composition rather than calling a generic $o(n+m)$ term negligible.

The finite non-heredity theorem shows that exact optimum blocks need not
realize this bound by direct completion. The more realistic constructive
target is the Pareto-profile inequality

\[
 K_{n,m}(F(n)+\eta_n,F(m)+\eta_m)^{2/3}
 \le F(n)^{2/3}+F(m)^{2/3}+O((n+m)^{1-\delta})
\]

with $\eta_n=o(n^{3/2})$. Determine what state, beyond the scalar maximum,
controls membership in the composable near-optimum family.

The exact multivertex identity supplies a sharper version of this program.
Do not attack the disproved weighted-entropy inequality from Section 3.13.
Instead prove the switching-orbit noncoverage statement there, first for
balanced doubling if necessary. Explore whether the energy-coloured two-point
law, a controlled hierarchy of coloured overlap laws, or additive-combinatorial
structure of the rank-one group can certify noncoverage with a power saving.

### Program E: fixed-half cut-discrepancy across orders

Use the exact-density parameter $H_n$ from the baseline as a genuinely
different language for the same asymptotic problem. The cleanest settling
target is, for some $\delta>0$,

\[
 \boxed{
 H_{n+k}^{2/3}
 \le H_n^{2/3}+H_k^{2/3}
 +O((n+k)^{1-\delta})}
\]

uniformly for proportional splits. Since $F(n)=4H_n+O(n)$, this has exactly
the same power-saving consequence as the target for $F^{2/3}$: conversion
through the $O(n)$ additive error costs only $O(n^{1/2})$ after taking the
$2/3$ power. An alternative is an asymptotically lossless amplification
theorem at edge density exactly $1/2+O(n^{-2})$, followed by the $O(N)$
edge-editing correction already proved.

Write

\[
 \Delta_G(U)=e_G(U,U^c)-\frac12|U||U^c|.
\]

For vertex blocks $V_1,V_2$ of sizes $n,k$, internal graphs $G_1,G_2$,
a bipartite cross graph $C$, and $P\subseteq V_1$, $R\subseteq V_2$, the
exact block identity is

\[
 \Delta_G(P\cup R)
 =\Delta_{G_1}(P)+\Delta_{G_2}(R)+\Delta_C(P,R),
\]

where

\[
 \Delta_C(P,R)
 =e_C(P,V_2\setminus R)+e_C(V_1\setminus P,R)
 -\frac12\left[|P|(k-|R|)+(n-|P|)|R|\right].
\]

The required cross-edge count is

\[
 e(C)=\left\lfloor\frac{\binom{n+k}{2}}2\right\rfloor
 -\left\lfloor\frac{\binom n2}2\right\rfloor
 -\left\lfloor\frac{\binom k2}2\right\rfloor
 =\frac{nk}{2}+O(1).
\]

The density condition itself is not the missing saving.  If
$S=2\mathbf1_C-1$, $t=\sum S_{ij}$, and
$D_\oplus(C)=\max_{P,R}|\Delta_C(P,R)|$, then exactly

\[
 D_\oplus(C)=\frac{|t|+\|S\|_{\infty\to1}}4,
 \qquad
 \|S\|_{\infty\to1}\ge\max\{n\mu_k,k\mu_n\}.
\]

Moreover, if $R_{n,k}(t_0)$ minimizes the rectangular norm at prescribed
feasible total $t_0$, row/column switching and entry editing prove

\[
 R_{n,k}\le R_{n,k}(t_0)
 \le R_{n,k}+\sqrt{nk}+|t_0|.
\]

Here the global fixed-half split has $|t_0|\le2$.  Thus exact density costs
only lower order, while the leading rectangular floor survives.  Any valid
composition must exploit cancellation with the two internal profiles; a
triangle or separate-norm estimate cannot close.

Do not choose $C$ independently. It must be optimized conditionally on the
two complete signed cut-deviation profiles so that the three terms cancel
simultaneously. Concrete attacks are:

- derive a dual or minimax description of the best such $C$ in terms of
  probability measures on pairs $(P,R)$, and determine whether its
  integrality gap is $o(N^{3/2})$;
- prove a puncturing, extension, or concatenation inequality that preserves
  the exact half-density layer up to $O(N)$ edits;
- identify a composable near-optimal Pareto state for the functions
  $P\mapsto\Delta_{G_1}(P)$ and $R\mapsto\Delta_{G_2}(R)$, rather than only
  their sup norms;
- test any proposed cross-order inequality against the exact small values,
  the order-$12$ through order-$14$ witnesses, and the exact $(2,8)$
  non-heredity obstruction.

Backurs--Bavarian supplies only a same-order $\Omega(n^{3/2})$ lower bound.
It does not address any item above. A claimed lower-bound route toward a
Paley subsequence must apply to every graph in the fixed half-density layer,
not merely to Paley, conference, quasirandom, or spectrally regular graphs.

### Program F: nonexistence

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
11. For a shell or cgf theorem, state the probability space, normalization,
    speed, proportional-split range, uniformity in the chosen signings, and
    the strict or power-saving margin. Verify every hypothesis used to pass
    between a cgf and a large-deviation rate.
12. For labeled convolution, state the relative-gauge group, quotient
    multiplicities, character normalization, and any low-dimensional gauge
    redundancy. A truncated Fourier hierarchy must include a proved tail
    bound.
13. For cut discrepancy, distinguish $D_{1/2}$ from $D_{\rho_G}$, all-cuts
    absolute discrepancy from one-sided balanced-half discrepancy, and exact
    half density from an asymptotically half-dense graph. Track the $O(n)$
    edge-editing loss explicitly.
14. For conference matrices, optimize over both eigenspaces, distinguish the
    maximum of a construction from the minimax value $F(n)$, and do not
    present dense Paley rigidity as progress unless a new universal
    restriction on every signing has actually been proved.
15. For finite computations, record class-enumeration completeness,
    switching/isomorphism multiplicities, exact witnesses, solver trust
    boundaries, independent replay, and corruption controls. Finite
    maximizer counts or alignment ratios are not asymptotic evidence without
    a stated statistical or structural theorem. The balanced profile
    calibration above is exact arithmetic over every balanced split of the
    banked witnesses, not an enumeration of every optimum class; replay
    verification/research_relative_profile_calibration.py and its corruption
    controls before extending the claim.

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

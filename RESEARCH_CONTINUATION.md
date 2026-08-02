# Continued attack after the second research note

Status date: 2026-08-01
Status: **new finite-temperature structure, no solution of the limit problem**

This ledger starts after `RESEARCH_LEDGER.md`.  It does not replace either
earlier ledger.  Write

\[
 Z_A(t)=\sum_{x\in\{\pm1\}^n,\ s\in\{\pm1\}}e^{tsQ_A(x)},
 \qquad
 h_n(t)=\min_A\log Z_A(t),
 \qquad m_n=\binom n2.
\]

## 1. Exact annealed-normalized block recursion

### Theorem 1 (block inequality and normalized subadditivity)

For every $n,k\ge1$ and $t\ge0$,

\[
 h_{n+k}(t)
 \le h_n(t)+h_k(t)-\log2+nk\log\cosh t.
\]

Consequently, if

\[
 r_n(t)=h_n(t)-\log2-m_n\log\cosh t,
\]

then

\[
 r_{n+k}(t)\le r_n(t)+r_k(t).
\]

At the extensive scaling, put

\[
 R_n(\beta)=r_n\left(\frac\beta{\sqrt n}\right).
\]

The exact rescaled inequality is

\[
 R_{n+k}(\beta)
 \le
 R_n\left(\beta\sqrt{\frac n{n+k}}\right)
 +R_k\left(\beta\sqrt{\frac k{n+k}}\right).
\]

#### Proof

Choose minimizers $A$ and $B$ for the two block sizes.  Write their
one-sided partition sums as

\[
 U_A^\sigma=\sum_xe^{t\sigma Q_A(x)},
 \qquad U_B^\sigma=\sum_ye^{t\sigma Q_B(y)},
 \qquad \sigma\in\{\pm1\}.
\]

Globally negating $B$ swaps $U_B^+$ and $U_B^-$.  For four nonnegative
numbers $a,b,c,d$,

\[
 \min\{ac+bd,ad+bc\}\le\frac{(a+b)(c+d)}2.
\]

Choose the global sign of $B$ using this inequality.  Now fill the $nk$
cross edges by independent uniform signs $C$.  For each fixed $x,y,\sigma$,

\[
 \mathbb E_C e^{t\sigma Q_C(x,y)}=(\cosh t)^{nk}.
\]

It follows that the average, over $C$, of the combined augmented partition
sum is at most

\[
 \frac12 Z_A(t)Z_B(t)(\cosh t)^{nk}.
\]

Some deterministic cross signing is no larger than this average.  Taking
logarithms proves the first inequality.  Since

\[
 m_{n+k}=m_n+m_k+nk,
\]

the second follows by exact cancellation.  In the last display, apply the
second inequality at $t=\beta/\sqrt{n+k}$ and rewrite that same $t$ as

\[
 \frac{\beta\sqrt{n/(n+k)}}{\sqrt n}
 \quad\hbox{and}\quad
 \frac{\beta\sqrt{k/(n+k)}}{\sqrt k}
\]

on the two subblocks.  \(\square\)

### Sharpening by one-sided imbalance

For a signing $A$, write

\[
 U_A^\sigma(t)=\sum_xe^{t\sigma Q_A(x)},
 \qquad
 \delta_A(t)=\frac{U_A^+(t)-U_A^-(t)}{Z_A(t)}.
\]

The orientation calculation in the proof is exact:

\[
 \begin{aligned}
 &\min\{U_A^+U_B^++U_A^-U_B^-,
 U_A^+U_B^-+U_A^-U_B^+\}\\
 &\hspace{25mm}=
 \frac{Z_AZ_B}{2}\left(1-|\delta_A\delta_B|\right).
 \end{aligned}
\]

Thus the block upper bound can be sharpened by

\[
 \log\left(1-|\delta_A(t)\delta_B(t)|\right).
\]

This correction is unavailable uniformly.  If the order is divisible by
four, label vertices cyclically and alternate edge signs along each orbit of
translation.  Every orbit has even length, including the diameter orbit of
length $n/2$, so translation sends $A$ to $-A$.  Hence it sends $Q_A$ to
$-Q_A$, giving $U_A^+=U_A^-$ and $\delta_A=0$ identically.  The original
$1/2$ orientation factor is therefore attained for arbitrarily large blocks.

### Proposition 2 (one-vertex cavity sandwich)

For every $n\ge1$ and $t\ge0$,

\[
 h_n(t)+\log2
 \le h_{n+1}(t)
 \le h_n(t)+\log2+n\log\cosh t.
\]

#### Proof

For an extension of $B$ by an incident sign vector
$b\in\{\pm1\}^n$, let $L_b(x)=\sum_i b_ix_i$.  Summing the new spin gives

\[
 Z_{B,b}(t)
 =2\sum_{x,s}e^{tsQ_B(x)}\cosh(tL_b(x)).
\]

The lower bound follows from $\cosh\ge1$, followed by
$Z_B(t)\ge e^{h_n(t)}$.  For the upper bound take a minimizing $B$ and
average over independent uniform $b_i$.  For every fixed $x$,

\[
 \mathbb E_b\cosh(tL_b(x))=(\cosh t)^n.
\]

Some deterministic $b$ is no larger than this average.  $\square$

### What this does and does not prove

For each fixed unscaled $t$, the annealed-normalized sequence $r_n(t)$ is
subadditive.  Fekete's lemma at fixed $t$ is not the desired theorem: here
$r_n(t)$ has a quadratic negative contribution, and the extensive problem
lies on the diagonal $t=\beta/\sqrt n$.  On that diagonal the exact recursion
changes each child temperature from $\beta$ to
$\beta\sqrt{n/(n+k)}$ or $\beta\sqrt{k/(n+k)}$.  Thus it is not an
ordinary subadditive inequality for $R_n(\beta)$ at a fixed $\beta$.

The missing repair is a comparison in the opposite temperature direction
whose total error is $o(n)$ in log-partition units.  A generic Lipschitz
comparison costs order $n$ when the block proportions stay bounded away
from zero, so the contraction cannot currently be discarded.

## 2. Entropy-aware Gaussian covariance lower bound

The earlier Gaussian-sign proof controls only a ground-state expectation.
The same covariances also carry a quantitative entropy certificate.

For a fixed signing $A$ of order $n\ge2$, let $q=n-1$ and

\[
 R^\pm=\frac12\left(I\pm\frac A{\sqrt q}\right)^2.
\]

For $0\le\lambda<1$, define positive-definite correlation matrices

\[
 \Sigma^\pm_\lambda=(1-\lambda)I+\lambda R^\pm.
\]

Let

\[
 \mathcal F_{n,\beta}(A)
 =\frac1{\beta n}\log
 \sum_{x,s}\exp\left(\frac{\beta sQ_A(x)}{\sqrt n}\right),
 \qquad
 \Psi_{n,\beta}=\min_A\mathcal F_{n,\beta}(A).
\]

### Theorem 3 (spectral-entropy certificate)

For every $A$, $\beta>0$, and $0\le\lambda<1$,

\[
 \begin{aligned}
 \mathcal F_{n,\beta}(A)
 \ge{}&
 \frac{\lambda}{\pi}\sqrt{\frac{n-1}{n}}
 +\frac{(n+1)\log2}{\beta n}\\
 &+\frac{
 \log\det\Sigma^+_\lambda+
 \log\det\Sigma^-_\lambda}{4\beta n}.
 \end{aligned}
\]

In particular, uniformly over all signings,

\[
 \boxed{
 \Psi_{n,\beta}
 \ge
 \frac{\lambda}{\pi}\sqrt{\frac{n-1}{n}}
 +\frac{(n+1)\log2}{\beta n}
 +\frac{\log d(\lambda)}{4\beta},}
\]

where

\[
 d(\lambda)=
 \begin{cases}
 (1-\lambda/2)^2,&0\le\lambda\le2/3,\\
 2\lambda(1-\lambda),&2/3\le\lambda<1.
 \end{cases}
\]

Since $d(\lambda)\ge(1-\lambda)^2$, the simpler bound with final term
$\log(1-\lambda)/(2\beta)$ also holds.

#### Proof

Both $R^\pm$ are correlation matrices: they are positive semidefinite, and
their diagonal entries equal one because \((A^2)_{ii}=q\).  Hence
$\Sigma^\pm_\lambda$ are positive-definite correlation matrices.

Let $X^\pm$ be the coordinatewise signs of centered Gaussian vectors with
covariances $\Sigma^\pm_\lambda$.  The Gaussian sign identity gives

\[
 \mathbb E[X_i^\pm X_j^\pm]
 =\frac2\pi\arcsin(\Sigma^\pm_{\lambda,ij}).
\]

For every edge,

\[
 a_{ij}(R^+_{ij}-R^-_{ij})=\frac2{\sqrt q}.
\]

Since $\arcsin'(u)\ge1$ on $(-1,1)$,

\[
 a_{ij}\left[
 \arcsin(\Sigma^+_{\lambda,ij})
 -\arcsin(\Sigma^-_{\lambda,ij})
 \right]
 \ge\frac{2\lambda}{\sqrt q}.
\]

After summing over the $m_n=nq/2$ edges,

\[
 \mathbb E Q_A(X^+)-\mathbb E Q_A(X^-)
 \ge\frac{2\lambda n\sqrt q}{\pi}.
\]

Construct a distribution on augmented states by choosing $S$ uniformly
from $\{\pm1\}$, then sampling $X^+$ conditional on $S=+1$ and
$X^-$ conditional on $S=-1$.  It satisfies

\[
 \mathbb E[SQ_A(X)]\ge\frac{\lambda n\sqrt q}{\pi}.
\]

It remains to bound its entropy.  Under the coordinatewise sign map, a
standard Gaussian becomes the uniform law $U_n$ on $\{\pm1\}^n$.  The
data-processing inequality for relative entropy gives

\[
 \begin{aligned}
 n\log2-H(X^\pm)
 &=D(\mathcal L(X^\pm)\|U_n)\\
 &\le D(N(0,\Sigma^\pm_\lambda)\|N(0,I))\\
 &=-\frac12\log\det\Sigma^\pm_\lambda,
 \end{aligned}
\]

where the last equality uses \(\operatorname{tr}\Sigma^\pm_\lambda=n\).
Therefore

\[
 H(S,X)
 \ge(n+1)\log2
 +\frac14\left(
 \log\det\Sigma^+_\lambda+
 \log\det\Sigma^-_\lambda\right).
\]

Apply the finite Gibbs variational principle

\[
 \mathcal F_{n,\beta}(A)
 =\max_p\left\{
 \frac{\mathbb E_p[SQ_A(X)]}{n\sqrt n}
 +\frac{H(p)}{\beta n}\right\}.
\]

This proves the pointwise bound.  For the uniform determinant estimate,
diagonalize $A/\sqrt q$ and denote one of its eigenvalues by $\mu$.  The
corresponding two eigenvalues of $\Sigma^\pm_\lambda$ are

\[
 1-\frac\lambda2+\frac\lambda2\mu^2\mathbin\pm\lambda\mu.
\]

Their product, as a function of $u=\mu^2\ge0$, is

\[
 f_\lambda(u)=
 \left(1-\frac\lambda2+\frac\lambda2u\right)^2-\lambda^2u.
\]

This convex quadratic is minimized on $u\ge0$ at $u=0$ when
$\lambda\le2/3$, and at $u=3-2/\lambda$ when $\lambda\ge2/3$.  Its minimum
is exactly $d(\lambda)$.  Multiplying over the $n$ eigenvalues gives

\[
 \det\Sigma^+_\lambda\det\Sigma^-_\lambda\ge d(\lambda)^n,
\]

which proves the boxed estimate.  The final simpler bound also follows
directly from $\Sigma^\pm_\lambda\succeq(1-\lambda)I$.  \(\square\)

### Corollary 4 (simple closed-form consequence)

Put

\[
 c_n=\frac1\pi\sqrt{\frac{n-1}{n}}.
\]

If $2\beta c_n\le1$, Theorem 3 with $\lambda=0$ gives

\[
 \Psi_{n,\beta}\ge\frac{(n+1)\log2}{\beta n}.
\]

If $2\beta c_n>1$, choosing

\[
 \lambda=1-\frac1{2\beta c_n}
\]

gives

\[
 \Psi_{n,\beta}
 \ge c_n+\frac{(n+1)\log2}{\beta n}
 -\frac{1+\log(2\beta c_n)}{2\beta}.
\]

Thus for fixed $\beta>\pi/2$,

\[
 \liminf_n\Psi_{n,\beta}
 \ge
 \frac1\pi+\frac{\log2}{\beta}
 -\frac{1+\log(2\beta/\pi)}{2\beta}.
\]

The right side tends to $1/\pi$ as $\beta\to\infty$, consistently with
the audited ground-state bound.

### Corollary 5 (sharper determinant consequence)

For every fixed $\beta>0$,

\[
 \liminf_n\Psi_{n,\beta}
 \ge \frac{\log2}{\beta}
 +\sup_{0\le\lambda<1}
 \left\{
 \frac\lambda\pi+\frac{\log d(\lambda)}{4\beta}
 \right\}.
\]

For example, if $\beta>3\pi/4$, the admissible choice
$\lambda=1-\pi/(4\beta)$ gives

\[
 \liminf_n\Psi_{n,\beta}
 \ge
 \frac1\pi+\frac{\log2}{\beta}-\frac1{4\beta}
 +\frac1{4\beta}\log\left[
 \frac\pi{2\beta}\left(1-\frac\pi{4\beta}\right)
 \right].
\]

This also tends to $1/\pi$ as $\beta\to\infty$, with a smaller entropy loss
than the determinant-free estimate in Corollary 4.

### Boundary of the result

This theorem proves that every adversarial signing has an explicit
energy--entropy tradeoff supported by two Gaussian-sign laws.  It is stronger
than using only one high-energy state, but it does not prove a thermodynamic
limit and does not improve the zero-temperature constant $1/\pi$.  The
determinant term retains the full spectrum of $A$; exploiting it beyond the
uniform bound would require a simultaneous lower bound relating spectral
flatness, Boolean free energy, and the adversarial minimum.

## 3. Exact code-noise and multiaffine structure

For a signing $A$ and $0\le u<1$, split its Eulerian polynomial into

\[
 \begin{aligned}
 P_A(u)&=\sum_{\substack{H\text{ Eulerian}\\|H|\text{ even}}}
 a_Hu^{|H|},\\
 R_A(u)&=\sum_{\substack{H\text{ Eulerian}\\|H|\text{ odd}}}
 a_Hu^{|H|}.
 \end{aligned}
\]

The one-sided high-temperature expansions are

\[
 U_A^\pm(t)=2^n(\cosh t)^{m_n}
 \left(P_A(\tanh t)\mathbin\pm R_A(\tanh t)\right).
\]

In particular, $\delta_A=R_A/P_A$.  If

\[
 p_n(u)=\min_A P_A(u),
\]

the block construction in Section 1 gives

\[
 p_{n+k}(u)\le p_n(u)p_k(u).
\]

### Theorem 6 (exact coset-noise law)

For $n\ge3$, let

\[
 \Gamma_n=\mathbb F_2^{m_n}/D_n,
 \qquad |\Gamma_n|=2^{m_n-n}.
\]

Let $X_u$ have independent edge bits with

\[
 \Pr(X_e=0)=\frac{1+u}{2},
 \qquad
 \Pr(X_e=1)=\frac{1-u}{2}.
\]

Then

\[
 \Pr(X_u+D_n=y+D_n)=2^{n-m_n}P_y(u),
\]

and consequently

\[
 p_n(u)=|\Gamma_n|\min_{\gamma\in\Gamma_n}
 \Pr(X_u+D_n=\gamma).
\]

Equivalently,

\[
 -\log p_n(u)=
 D_\infty\left(
 \operatorname{Unif}(\Gamma_n)
 \,\middle\|\,
 \mathcal L(X_u+D_n)
 \right).
\]

The function $p_n(u)$ is nonincreasing in $u$.  Moreover,

\[
 (1-u^2)^{n/2}p_n(u)\le p_{n+1}(u)\le p_n(u).
\]

#### Proof

Fourier inversion on the quotient gives

\[
 \begin{aligned}
 \Pr(X_u+D_n=y+D_n)
 &=\frac{|D_n|}{2^{m_n}}
 \sum_{H\in D_n^\perp}(-1)^{H\cdot y}u^{|H|}\\
 &=2^{n-m_n}P_y(u).
 \end{aligned}
\]

If $0\le u'\le u$, write $u'=us$.  The xor of independent edge-noise
vectors with biases $u$ and $s$ has bias $us$, so the quotient laws obey

\[
 \mu_{u'}=\mu_u*\mu_s.
\]

Convolution with a probability law cannot decrease the smallest atom.  This
proves monotonicity.  Finally use

\[
 h_n(t)=(n+1)\log2+m_n\log\cosh t+\log p_n(\tanh t)
\]

in Proposition 2 and the identity

\[
 \log\cosh t=-\frac12\log(1-\tanh^2t).
\]

\(\square\)

Thus, with $q_n=-\log p_n$,

\[
 0\le q_{n+1}(u)-q_n(u)
 \le-\frac n2\log(1-u^2).
\]

At $u=\tanh(\beta/\sqrt{n+1})$ this fixed-$u$ increment is
$O_\beta(1)$.  It does not control the simultaneous temperature shift from
order $n+1$ back to the same extensive $\beta$ at order $n$; monotonicity is
in the wrong direction for the required upper comparison.

### Theorem 7 (exact multiaffine box relaxation)

For a weighted graph or multigraph $G$, define

\[
 \mathcal P_G(w)=
 \sum_{\substack{\partial H=\varnothing\\|H|\text{ even}}}
 \prod_{e\in H}w_e.
\]

Equivalently,

\[
 \mathcal P_G(w)=2^{-(|V|+1)}
 \sum_{x,s}\prod_{e=ij}(1+s w_ex_ix_j).
\]

It is positive for $|w_e|<1$ and is affine in each edge coordinate.
Consequently, for every $0<u<1$,

\[
 \min_{a_e\in\{\pm1\}}\mathcal P_G(ua)
 =\min_{w\in[-u,u]^E}\mathcal P_G(w).
\]

#### Proof

Every factor in the averaged product is positive when $|w_e|<1$.  Expanding
the product and averaging $x,s$ retains exactly the even-cardinality Eulerian
subgraphs, proving the two representations.  At a box minimizer, hold all
coordinates but one fixed.  The objective is affine in the remaining
coordinate, so one of its two endpoints is no larger.  Push the coordinates
to endpoints one at a time.  \(\square\)

This has no integrality gap, but it does not contradict the collapsed cube
relaxation in the physical couplings: here $w_e=\tanh J_e$, and the edge
normalizers have already been divided out.

### Exact deletion--contraction and the state-space wall

Also put

\[
 \mathcal R_G(w)=
 \sum_{\substack{\partial H=\varnothing\\|H|\text{ odd}}}
 \prod_{e\in H}w_e.
\]

For a nonloop edge $e$, retaining loops and parallel edges after contraction,

\[
 \boxed{
 \begin{aligned}
 \mathcal P_G
 &=\mathcal P_{G\setminus e}
 +w_e\left(\mathcal R_{G/e}-\mathcal R_{G\setminus e}\right),\\
 \mathcal R_G
 &=\mathcal R_{G\setminus e}
 +w_e\left(\mathcal P_{G/e}-\mathcal P_{G\setminus e}\right).
 \end{aligned}}
\]

Indeed, Eulerian subgraphs of $G/e$ correspond to subgraphs of
$G\setminus e$ whose boundary is either empty or the two endpoints of $e$;
including $e$ swaps both the boundary requirement and edge parity.  Hence

\[
 \min_{w_e\in[-u,u]}\mathcal P_G
 =\mathcal P_{G\setminus e}
 -u\left|\mathcal R_{G/e}-\mathcal R_{G\setminus e}\right|.
\]

The recursion does not close on $p_n$: contraction creates signed
multigraphs, and odd and boundary sectors remain coupled.  The fully closed
state is

\[
 T_G(S,\varepsilon)=
 \sum_{\substack{\partial H=S\\|H|\equiv\varepsilon\pmod2}}w_H,
\]

which updates under insertion of $e=uv$ by

\[
 T_{G+e}(S,\varepsilon)=T_G(S,\varepsilon)
 +w_eT_G(S\triangle\{u,v\},\varepsilon\mathbin\oplus1).
\]

A new vertex requires every even boundary sector:

\[
 P_{n+1}=\sum_{\substack{S\subseteq[n]\\|S|\text{ even}}}
 \left(\prod_{i\in S}w_{i,n+1}\right)T_{K_n}(S,0).
\]

This exponential state is genuine.  An all-positive triangle has
$(P,R)=(1,u^3)$.  Extending it by star signs $(+,+,+)$ gives
$P=1+3u^4$, while signs $(-,+,+)$ give $P=1-u^4$.  Even the scalar pair
$(P,R)$ of the deleted graph cannot determine a vertex extension.

## 4. Amplification and near-saturation obstructions

### Theorem 8 (an exact nonconvex relaxation and its composition cost)

For every $n$,

\[
 F(n)=\min_{a\in[-1,1]^{m_n}}
 \left[G_n(a)+\sum_e(1-|a_e|)\right].
\]

#### Proof

Round each coordinate of $a$ to a nearest sign $\sigma_e$.  The triangle
inequality and $G_n(v)\le\|v\|_1$ give

\[
 F(n)\le G_n(\sigma)
 \le G_n(a)+G_n(\sigma-a)
 \le G_n(a)+\sum_e(1-|a_e|).
\]

An optimal sign vector gives equality in the other direction.  \(\square\)

More generally, let $p_n:[0,1]\to[0,\infty)$ satisfy $p_n(1)=0$, and suppose

\[
 \min_{a\in[-1,1]^{m_n}}
 \left[G_n(a)+\sum_ep_n(|a_e|)\right]=F(n).
\]

Testing the point $a=t\sigma$ for an optimal signing $\sigma$ gives

\[
 p_n(t)\ge\frac{(1-t)F(n)}{m_n},
 \qquad
 p_n(0)\ge\frac2{\pi\sqrt{n-1}}.
\]

Therefore a balanced two-block composition at total order $N$, with zero
cross couplings, pays cross penalty at least

\[
 \frac{2n_1n_2}{\pi\sqrt{N-1}}=\Theta(N^{3/2}).
\]

The same conclusion holds if the relaxed minimum is only
$F(n)-o(n^{3/2})$.  Thus every edge-separable exact or asymptotically exact
saturation penalty retains a leading composition wall.

### Theorem 9 (near-saturated cross discrepancy)

For $C\in[-1,1]^{n\times k}$, let

\[
 B(C)=\max_{x\in\{\pm1\}^n,\ y\in\{\pm1\}^k}|x^{\mathsf T}Cy|.
\]

Then

\[
 B(C)\ge
 \max\left\{
 \frac{\|C\|_F^2}{\sqrt{2n}},
 \frac{\|C\|_F^2}{\sqrt{2k}}
 \right\}.
\]

#### Proof

Randomize $y$ and then choose $x_i$ as the sign of row $i$ dotted with $y$.
The sharp elementary Khintchine lower bound gives

\[
 B(C)\ge\frac1{\sqrt2}\sum_i\|C_{i,*}\|_2.
\]

Since each squared row norm is at most $k$,

\[
 \sum_i\|C_{i,*}\|_2
 \ge\frac{\|C\|_F^2}{\sqrt k}.
\]

Apply the same argument to $C^{\mathsf T}$.  \(\square\)

If $\sum_{ij}(1-c_{ij}^2)=o(nk)$ and the aspect ratio stays bounded, this is
$\Theta((n+k)^{3/2})$.  In a two-block quadratic form, globally flipping one
block changes an energy $U+V$ to $U-V$, so

\[
 \max\{|U+V|,|U-V|\}=|U|+|V|\ge|V|.
\]

The full ground state is at least $B(C)$.  Near saturation therefore makes the
cross term leading before any rounding occurs.

### Theorem 10 (low-rank block-lift obstruction)

For a $k\times k$ sign matrix $H$, put

\[
 B(H)=\max_{u,v\in\{\pm1\}^k}|u^{\mathsf T}Hv|,
 \qquad \lambda=\|H\|_{\mathrm{op}}.
\]

Then

\[
 B(H)\ge\frac{\lambda^3}{k}.
\]

If $\operatorname{rank}H\le r$, then

\[
 B(H)\ge\frac{k^2}{r^{3/2}}.
\]

#### Proof

Let $p,q$ be unit top singular vectors.  Uniform row and column Euclidean
norms give

\[
 \|p\|_\infty,\|q\|_\infty\le\frac{\sqrt k}{\lambda}.
\]

Choose independent signs $v_j$ with
$\mathbb Ev_j=q_j/\|q\|_\infty$.  Some outcome satisfies

\[
 p^{\mathsf T}Hv\ge\frac\lambda{\|q\|_\infty}.
\]

For that outcome take $u=\operatorname{sign}(Hv)$.  Then

\[
 B(H)\ge\|Hv\|_1
 \ge\frac{|p^{\mathsf T}Hv|}{\|p\|_\infty}
 \ge\frac{\lambda^3}{k}.
\]

Finally $\|H\|_F=k$ implies $\lambda\ge k/\sqrt r$.  \(\square\)

For any cloud blow-up, the full quadratic maximum is at least $B(H_{ij})$ for
each cross block.  Fix cloud patterns, multiply every cloud by an independent
global sign, and observe that $u^{\mathsf T}H_{ij}v$ is a degree-two Walsh
coefficient of the resulting energy; every Walsh coefficient is bounded by
the sup norm.  Hence a correct-scale cross block $B(H)=O(k^{3/2})$ must have
rank $\Omega(k^{1/3})$.  Bounded-rank lifts necessarily remain on the old
$k^2$ scale.

### Exact Kronecker counterexample

Let $A$ be the order-five Seidel matrix with negative cycle edges and positive
chords.  Its Boolean energies are exactly $\{-4,0,4\}$, so Parseval and parity
give $M(A)=F(5)=4$.  For diagonal sign matrices $D,E$, define

\[
 C_{D,E}=(A+D)\otimes(A+E)-D\otimes E.
\]

This is a valid order-25 signing.  For product spins $u\otimes v$, if
$s=\operatorname{tr}D$ and $t=\operatorname{tr}E$, then

\[
 Q_{C_{D,E}}(u\otimes v)
 =\frac12[(s+2Q_A(u))(t+2Q_A(v))-st].
\]

Since $Q_A$ attains both signs of four,

\[
 M(C_{D,E})\ge4(8+|s|+|t|).
\]

This is at least $48>5^{3/2}M(A)$ unless $|s|=|t|=1$.  In the remaining
case, the dihedral symmetries of the cycle and the map $i\mapsto2i\pmod5$
reduce each diagonal completion to two representatives.  Exact spin
certificates give energies $88,-80,88$ for the three unordered representative
pairs.  Therefore every separated diagonal completion of this Kronecker
square has

\[
 M(C_{D,E})>5^{3/2}M(A).
\]

For the canonical strong tensor $P=A+I$, an exact order-25 witness has energy
100 for $P\otimes P-I$.  Iterating its product witness gives at order
$N=25^s$

\[
 \frac{M(P^{\otimes2s}-I_N)}{N^{3/2}}
 \ge\frac12(1.8^s-0.2^s).
\]

The matrices, all orbit checks, and corruption controls are reproduced in
`verification/verify_amplification_obstructions.py`.  This is an exact finite
counterexample to these tensor mechanisms, not a no-go theorem for every
nonlocal high-rank lift.

## 5. Quenched-state refinements and further walls

### Exact Bellman identity

For an order-$n$ signing $B$, let $\mu_{B,t}$ be the Gibbs law on $(X,S)$
proportional to $e^{tSQ_B(X)}$.  Extending by incident signs
$b\in\{\pm1\}^n$ gives

\[
 \frac{Z_{(B,b)}(t)}{Z_B(t)}
 =2\,\mathbb E_{\mu_{B,t}}\cosh(t b\cdot X).
\]

After removing the vertex factor and the $n$ new annealed edge factors, the
exact cavity reward is

\[
 \log\mathbb E_{\mu_{B,t}}\cosh(t b\cdot X)
 -n\log\cosh t.
\]

The minimizing $B$ in this Bellman step need not minimize the order-$n$
scalar free energy.  Thus a closed interpolation requires a state describing
the Gibbs law, not merely its optimized scalar value.

Covariance and the limiting pair-overlap law do not determine this reward.
For $n=2^{2k}$, let $\mu_H$ be uniform on the signed rows of a Sylvester
Hadamard matrix.  It has covariance $I$.  Parseval and convexity of
$y\mapsto\cosh(t\sqrt y)$ give

\[
 \min_b\mathbb E_{\mu_H}\cosh(t b\cdot X)
 \ge\cosh(t\sqrt n),
\]

with equality for a bent Boolean sign vector $b$.  The uniform law on the
whole cube also has covariance $I$, but its reward before normalization is
$(\cosh t)^n$.  At $t=\beta/\sqrt n$, their log costs tend respectively to
$\log\cosh\beta$ and $\beta^2/2$.  Both pair-overlap laws tend to
$\delta_0$.  Entropy distinguishes these examples, so this rules out only
covariance or pair overlap as a complete state, not a richer
energy--entropy order parameter.

### Exact negative-replica sandwich

Let $K_n=2^{m_n}$ be the number of signings and let $A$ be uniform over them.
For every $q>0$,

\[
 \Psi_{n,\beta}
 \le
 -\frac1{q\beta n}\log\mathbb E_A[Z_A^{-q}]
 \le
 \Psi_{n,\beta}+\frac{m_n\log2}{q\beta n}.
\]

#### Proof

If $Z_*=\min_AZ_A$, then

\[
 K_n^{-1}Z_*^{-q}\le\mathbb E_AZ_A^{-q}\le Z_*^{-q}.
\]

Take negative logarithms and divide by $q\beta n$.  \(\square\)

Fixed $q$ is too coarse: its error is $\Theta(n)$.  Taking $q=\lambda n$
gives asymptotic error $\log2/(2\lambda\beta)$.  Thus existence of the
negative-replica free-energy limit for unbounded fixed $\lambda$ would settle
the fixed-$\beta$ minimax limit.  Unlike an ordinary annealed moment, this
object emphasizes low-partition disorders and is compatible with the outer
minimum.

Typical iid disorder cannot simply replace that minimum.  The rigorous greedy
bound gives almost surely

\[
 \liminf_n\mathcal F_{n,\beta}(A_{\mathrm{iid}})
 \ge\frac23\sqrt{\frac2\pi}=0.531923\ldots,
\]

whereas the adversarial Paley upper bound gives

\[
 \limsup_n\Psi_{n,\beta}\le\frac12+\frac{\log2}{\beta}.
\]

These are separated for
$\beta>\log2/(\frac23\sqrt{2/\pi}-\frac12)=21.7131\ldots$.

### The scalar inequalities still do not force convergence

There is an explicit abstract countermodel satisfying the exact scalar block
and cavity inequalities, the annealed upper bound, temperature monotonicity,
and the uniform entropy barrier, while oscillating on the extensive diagonal.

Set $\varepsilon=c=0.1$, $b_1=0$, and for integers $n\ge2$ put

\[
 b_n=n^2(1+\varepsilon\sin\log n),
 \qquad
 r_n^*(t)=n\log2-cb_n\log\cosh t.
\]

The continuous function
$b(x)=x^2(1+\varepsilon\sin\log x)$ has

\[
 b''(x)=2+\varepsilon(\sin\log x+3\cos\log x)>0.
\]

It is increasing and convex, so $b_{n+k}\ge b_n+b_k$; the cases involving
$b_1=0$ follow directly from monotonicity.  Hence

\[
 r_{n+k}^*(t)\le r_n^*(t)+r_k^*(t).
\]

Also

\[
 b'(x)=x(2+2\varepsilon\sin\log x+\varepsilon\cos\log x),
\]

so $c(b_{n+1}-b_n)\le n$.  Therefore

\[
 \log2-n\log\cosh t
 \le r_{n+1}^*(t)-r_n^*(t)\le\log2,
\]

exactly the cavity interval.  The centered quantity
$r_n^*(t)-n\log2$ is nonincreasing in $t$, and $r_n^*(t)\le n\log2$.
Nevertheless,

\[
 \frac1nR_n^*(\beta)
 =\log2-\frac{c\beta^2}{2}
 (1+\varepsilon\sin\log n)+o(1),
\]

which has distinct subsequential limits.

This example can also be kept above Theorem 3's scalar lower barrier.  The
elementary inequality $d(\lambda)\le e^{-\lambda}$ implies

\[
 \sup_\lambda\left\{
 \frac{\beta\lambda}{\pi}+\frac14\log d(\lambda)
 \right\}\le\frac{\beta^2}{\pi^2}.
\]

Thus that barrier for $R_n/n$ is at most

\[
 \log2-\left(\frac14-\frac1{\pi^2}\right)\beta^2,
\]

whereas the countermodel is at least
$\log2-0.055\beta^2$.  The exact scalar inequalities proved so far therefore
cannot by themselves imply the desired limit.  New optimizer-specific state
information is logically necessary.  This abstract model is not an actual
signing partition function and does not encode every separate fact about the
problem, notably the Paley ground-state upper bound; its scope is only the
listed recursion, cavity, monotonicity, annealed, and entropy inputs.

## 6. Ground-state hereditary cavity inequalities

The following two inequalities were recovered from a separate frontier-model
run and then rederived and checked independently. They are exact statements
for every signing. They strengthen the hereditary information available at
zero temperature, but their unconditional terms are below the
$n^{3/2}$ scale.

Let the vertex set of a signing $A$ be partitioned as $S\sqcup T$, with
$|S|=s\ge1$ and $|T|=k\ge1$. Write $C=A_{S,T}$. For fixed
$x\in\{\pm1\}^S$, set

\[
 q=Q_{A[S]}(x),\qquad h=C^{\mathsf T}x,
 \qquad r(y)=Q_{A[T]}(y).
\]

Thus

\[
 Q_A(x,y)=q+h\mathbin\cdot y+r(y).
\]

### Theorem 11 (exact squared cavity inequality)

For every $x\in\{\pm1\}^S$,

\[
 \boxed{
 M(A)^2\ge q^2+\|C^{\mathsf T}x\|_2^2+\binom k2.}
\]

In particular, if $x$ attains $M(A[S])$, then

\[
 M(A)^2\ge M(A[S])^2+\binom k2+k(s\bmod2),
\]

and hence

\[
 \boxed{
 F(s+k)^2\ge F(s)^2+\binom k2+k(s\bmod2).}
\]

The companion inequality obtained by interchanging $s$ and $k$ also holds.

#### Proof

Under uniform $y\in\{\pm1\}^T$, the constant character, the degree-one
characters $y_j$, and the degree-two characters $y_i y_j$ are mutually
orthogonal. Parseval therefore gives the exact identity

\[
 \mathbb E_y Q_A(x,y)^2
 =q^2+\sum_{j\in T}h_j^2+\binom k2.
\]

The maximum of $Q_A(x,y)^2$ over $y$ is at least its average, and it is at
most $M(A)^2$. This proves the boxed matrix inequality. Each $h_j$ is a
sum of $s$ signs, so $h_j\equiv s\pmod2$ and consequently
$h_j^2\ge s\bmod2$. Choose $x$ attaining $M(A[S])$ and use
$M(A[S])\ge F(s)$. The resulting lower bound holds for every order-$(s+k)$
signing, so minimizing over $A$ proves the assertion for $F$. \(\square\)

For $k=1$, the parity lattice of the energies upgrades the last inequality
to the already known puncturing bound
$F(s+1)\ge F(s)+(s\bmod2)$. For proportional $s$ and $k$, however, the new
unconditional contribution is only $O(k^2)$ in $M(A)^2$, while
$F(s+k)^2$ and the squared maximum of a near-optimal signing are of order
$(s+k)^3$.

### Theorem 12 (exact block-pairing inequality)

For every $x\in\{\pm1\}^S$,

\[
 \boxed{
 M(A)\ge |Q_{A[S]}(x)|
 +\mathbb E_{y\in\{\pm1\}^T}
   |\langle C^{\mathsf T}x,y\rangle|.}
\]

Consequently, the sharp $p=1$ Khintchine inequality gives

\[
 M(A)\ge |Q_{A[S]}(x)|
 +\frac1{\sqrt2}\|C^{\mathsf T}x\|_2.
\]

#### Proof

Since $r(-y)=r(y)$, pairing $y$ and $-y$ gives

\[
 \begin{aligned}
 &\max\{|q+r(y)+h\mathbin\cdot y|,
          |q+r(y)-h\mathbin\cdot y|\}\\
 &\hspace{30mm}=|q+r(y)|+|h\mathbin\cdot y|.
 \end{aligned}
\]

Both energies on the left are bounded in absolute value by $M(A)$. Average
over $y$, then use Jensen's inequality and
$\mathbb E_y r(y)=0$:

\[
 M(A)\ge\mathbb E_y|q+r(y)|+\mathbb E_y|h\mathbin\cdot y|
 \ge |q|+\mathbb E_y|h\mathbin\cdot y|.
\]

This is the first assertion. The second is the sharp real Khintchine lower
bound
$\mathbb E_\varepsilon|\sum_jh_j\varepsilon_j|
\ge\|h\|_2/\sqrt2$. \(\square\)

There is also a sharper parity-only consequence. Put

\[
 \mu_k=\mathbb E\left|\varepsilon_1+\cdots+\varepsilon_k\right|
 =\frac{k}{2^{k-1}}
   \binom{k-1}{\lfloor(k-1)/2\rfloor}.
\]

If $s$ is odd, every coordinate of $h$ is a nonzero odd integer. The
Rademacher $L^1$ norm is coordinatewise nondecreasing in the absolute
coefficients: conditioning on all other signs reduces this to

\[
 \frac{|z+a|+|z-a|}{2}=\max\{|z|,|a|\}.
\]

Thus $\mathbb E|h\mathbin\cdot y|\ge\mu_k$, and choosing $x$ to attain
$M(A[S])$ yields

\[
 \boxed{F(s+k)\ge F(s)+(s\bmod2)\mu_k.}
\]

This is an exact additive refinement, but $\mu_k\sim\sqrt{2k/\pi}$, so its
unconditional gain is again below the leading scale. When $s$ is even, the
cross columns can all balance against the chosen $x$, so parity supplies no
positive universal reward.

### A second precise convergence target

Set

\[
 H(n)=F(n)^{2/3}.
\]

Exact subadditivity is already false at the first nontrivial test:
$H(4)=4^{2/3}>2H(2)=2$. A power-saving defect would nevertheless be enough.

### Theorem 13 (power-saving near-subadditivity would settle the problem)

Suppose there are constants $C<\infty$ and $\delta>0$ such that, for every
$n,m\ge1$,

\[
 H(n+m)\le H(n)+H(m)+C(n+m)^{1-\delta}.
\]

Then $H(n)/n$ converges, and therefore $F(n)/n^{3/2}$ converges.

#### Proof

Replacing $\delta$ by $\min\{\delta,1\}$ if necessary, assume first that
$0<\delta\le1$. The known upper bound gives $H(n)=O(n)$. Fix a block length
$k$. Combine
$q$ copies of a $k$-block along a balanced binary tree. At a level whose
subtrees contain at most $2^j$ leaves, there are at most
$\lceil q/2^j\rceil$ merges and each defect is
$O((2^jk)^{1-\delta})$. If $L=\lceil\log_2q\rceil$, the total defect is at
most a constant times

\[
 qk^{1-\delta}\sum_{j=1}^L2^{-j\delta}
 +k^{1-\delta}\sum_{j=1}^L2^{j(1-\delta)}
 =O(qk^{1-\delta}).
\]

For $\delta<1$, the second sum is
$O(k^{1-\delta}q^{1-\delta})$, and for $\delta=1$ it is
$O(\log q)$; both have the asserted bound. Hence, uniformly in $q$,

\[
 H(qk)\le qH(k)+O(qk^{1-\delta}).
\]

For arbitrary $N=qk+r$, $0\le r<k$, use the preceding estimate directly if
$r=0$. If $r>0$, one final application of the assumed inequality, together
with $H(r)=O(k)$, gives

\[
 \limsup_{N\to\infty}\frac{H(N)}N
 \le\frac{H(k)}k+O(k^{-\delta}).
\]

Choose a sequence $k\to\infty$ along which $H(k)/k$ tends to its liminf.
The displayed inequality makes the limsup no larger than that liminf.
Finally,
$F(n)/n^{3/2}=(H(n)/n)^{3/2}$. \(\square\)

Thus an alternative settling lemma is a sign-matrix composition theorem
with a power saving in the $H$ scale. Equivalently, using the known uniform
bounds $F(t)=\Theta(t^{3/2})$, it would suffice to prove a uniform estimate
of the form

\[
 F(n+m)
 \le
 \left(F(n)^{2/3}+F(m)^{2/3}\right)^{3/2}
 +O((n+m)^{3/2-\delta}).
\]

An unspecified $o(n+m)$ defect in the $H$ inequality is not by itself a
safe convergence criterion; the summability across dyadic scales is the
essential feature of the theorem above. The present cavity inequalities are
lower bounds and therefore do not supply this missing upper composition.

## 7. Current best continuation targets

The most concrete finite-temperature settling lemma is the negative-replica
limit. For
fixed $\beta,\lambda>0$, define

\[
 \Xi_{n,\beta,\lambda}
 =-\frac1{\lambda\beta n^2}
 \log\mathbb E_A\left[Z_A(\beta/\sqrt n)^{-\lambda n}\right],
\]

where $A$ is uniform over all signings.

> **Negative-replica limit lemma.**  For every fixed $\beta>0$ and for an
> unbounded set of $\lambda>0$, the limit
> $\lim_n\Xi_{n,\beta,\lambda}$ exists.

The exact sandwich shows that this lemma implies

\[
 \limsup_n\Psi_{n,\beta}-\liminf_n\Psi_{n,\beta}
 \le\frac{\log2}{2\lambda\beta}
\]

for arbitrarily large $\lambda$, and hence proves the fixed-$\beta$ minimax
limit.  Proving it for arbitrarily large $\beta$ then settles the original
MathOverflow question.

The ground-state alternative is Theorem 13's power-saving composition bound
for $H(n)=F(n)^{2/3}$. It would bypass temperature entirely. Current block
constructions fail at exactly this point: the internal block energies and the
cross field can cancel for the same spin assignment, while bounding the two
pieces separately incurs a leading $N^{3/2}$ cost. Theorems 11 and 12 give
new hereditary lower information, but no optimizer-compatible upper
construction.

The alternative structural target is an optimizer-specific stability theorem
for the Bellman state.  It would supplement

\[
 R_{n+k}(\beta)
 \le R_n(\beta\sqrt{n/(n+k)})
 +R_k(\beta\sqrt{k/(n+k)})
\]

with a reverse or same-temperature comparison of total cost $o(n+k)$.

The recursion alone is insufficient because it compares different
temperatures, and the abstract countermodel proves that scalar cavity bounds
cannot repair it.  Covariance and limiting pair overlap are also too coarse.
A successful state must retain enough energy--entropy or boundary-sector
information to determine the optimized cavity reward.  On the coding side,
the corresponding target is a strong data-processing or reverse-noise theorem
for the worst coset probability that remains effective when
$u=\tanh(\beta/\sqrt n)$.

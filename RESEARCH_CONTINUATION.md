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

## 7. Exact finite obstruction to optimizer heredity

The preceding walls show that the scalar value $F(n)$ does not close under
the available analytic recursions. Exact cross-block optimization now shows
that this loss of state already occurs among finite ground-state optimizers.

For signings $A$ and $B$ of orders $s$ and $k$, define

\[
 J(A,B)=\min_{C\in\{\pm1\}^{s\times k}}
 M\begin{pmatrix}A&C\\ C^{\mathsf T}&B\end{pmatrix}.
\]

For an order-$n$ signing $B$, its optimum one-vertex extension value is

\[
 E(B)=J([0],B).
\]

The equivalence relation in the following computation is generated by vertex
permutations and diagonal switching. Global negation is not identified when
it belongs to a different switching class.

### Theorem 14 (class-sensitive extension and non-heredity)

Exact exhaustive computation gives the following distribution of $E(B)$
over all switching-permutation classes satisfying $M(B)=F(n)$:

\[
\begin{array}{c|c|c|c}
n&\text{number of classes}&\text{distribution of }E(B)&\text{comparison}\\ \hline
2&1&3^{\times1}&3\\
3&2&4^{\times2}&4\\
4&1&4^{\times1}&4\\
5&1&5^{\times1}&5\\
6&1&9^{\times1}&9\\
7&6&10^{\times4},\ 12^{\times2}&10\\
8&2&12^{\times2}&12\\
9&15&13^{\times4},\ 15^{\times11}&F(10)=13\\
10&2&17^{\times1},\ 19^{\times1}&13\le F(11)\le17
\end{array}
\]

Thus an arbitrary optimum need not admit an optimum one-vertex extension.
At order 7, the two failing classes first become infeasible when the
configurations with $|Q_B(x)|=7$ are added to the maximizing configurations.
At order 9, the maximizing configurations with $|Q_B(x)|=12$ already make
extension at value 13 infeasible in all eleven failing classes.
At order 10 the two exact optimum classes have extension values 17 and 19;
the better class gives the rigorous new finite bound $F(11)\le17$. The
Gaussian lower bound and energy parity give $F(11)\ge13$. No exact value of
$F(11)$ is claimed. A direct upper witness uses order-10 graph6 record
`HCRbczQ` and incident-edge mask 440, with bit $j$ equal to one when the edge
to vertex $j$ is negative.

For two nontrivial blocks, the complete table through total order 10 is as
follows. Multiplicities count ordered pairs of switching classes after the
block orders have been fixed.

\[
\begin{array}{c|c|c}
(s,k)&\text{distribution of }J(A,B)
 &F(s+k)\\ \hline
(2,2)&4^{\times1}&4\\
(2,3)&4^{\times2}&4\\
(2,4)&5^{\times1}&5\\
(3,3)&5^{\times2},\ 7^{\times2}&5\\
(2,5)&9^{\times1}&9\\
(3,4)&9^{\times2}&9\\
(2,6)&10^{\times1}&10\\
(3,5)&10^{\times2}&10\\
(4,4)&10^{\times1}&10\\
(2,7)&12^{\times6}&12\\
(3,6)&12^{\times2}&12\\
(4,5)&12^{\times1}&12\\
(2,8)&15^{\times2}&13\\
(3,7)&13^{\times6},\ 15^{\times6}&13\\
(4,6)&13^{\times1}&13\\
(5,5)&13^{\times1}&13
\end{array}
\]

The $(2,8)$ row gives the sharp finite obstruction

\[
 \boxed{
 \min_{M(A)=F(2),\ M(B)=F(8)}J(A,B)=15>F(10)=13.}
\]

Consequently no optimal order-10 signing contains an optimal order-8
principal submatrix. In particular, there is no nested chain of principal
submatrices that is optimal at each of orders 8, 9, and 10.

The obstruction is repaired by the smallest possible internal slack. Among
the 1,044 root-normalized unlabeled order-8 representatives, 104 have maximum
12, and 68 of those can be joined to the order-2 optimum to give maximum 13.
The first deterministic witness has graph6 record `F?reg` and cross-mask
52010. The graph6 record encodes the negative residual edges after normalizing
the root edges positive. Cross-mask bit $8i+j$ is one exactly when the edge
from left vertex $i$ to right vertex $j$ is negative. Therefore, if

\[
 K_{s,k}(u,v)=
 \min_{\substack{M(A)\le u\\M(B)\le v}}J(A,B),
\]

then the exact Pareto-profile identities are

\[
 \boxed{K_{2,8}(1,10)=15,\qquad K_{2,8}(1,12)=13.}
\]

#### Proof and certificate boundary

Switching makes every edge incident to a fixed root positive. The remaining
negative edges form an ordinary graph on $n-1$ vertices. Nauty `geng`
enumerates one representative of every unlabeled residual graph; the script
asserts the complete graph counts and deterministic stream hashes. Direct
enumeration of the spin cube computes $M(B)$, and NetworkX groups the optimum
root-normalized records into switching-permutation classes.

For fixed $A$ and $B$, and a proposed bound $L$, the existence of a cross
signing is the finite system

\[
 -L\le Q_A(x)+Q_B(y)+
 \sum_{i,j}c_{ij}x_i y_j\le L
 \qquad(x,y\in\{\pm1\}).
\]

It is a binary linear feasibility problem in the $c_{ij}$. The verifier uses
a deterministic cutting-plane loop: solve the current subsystem, recompute
the candidate on every spin pair, and add every violated pair. An infeasible
subsystem certifies infeasibility of the full system; a feasible witness is
recomputed over the entire cube. Z3 is trusted for infeasibility when the
cross block has more than 16 entries.

Every one-vertex extension in the table has at most 10 cross entries and is
instead checked by direct enumeration. For every pair with at most 16 cross
entries, including the critical $(2,8)$ row, all cross signings are also
enumerated directly.
Negating the full cross block is induced by flipping every spin in one block,
so one cross coefficient may be fixed positive without loss.

Finally, the script enumerates all 1,044 residual graphs at order 8. None of
the four optimum root-normalized representatives admits a $(2,8)$ completion
below 15. Since order-8 maxima are even, the next possible internal value is
12; the displayed record and mask recompute to combined maximum 13. This
proves the two Pareto identities, subject only to the stated `geng`
completeness boundary. All calculations and corruption controls are in
`verification/research_cross_block_composition.py`. \(\square\)

### Theorem 15 (weighted covering-radius Bellman identity)

The entire one-vertex extension problem has an exact coding-theoretic state,
not just the finite class table above.  Write

\[
 \mathcal P_n=\{\pm1\}^n/\{x\sim-x\},\qquad
 d_\pm([u],[v])=\min\{d_H(u,v),d_H(u,-v)\}.
\]

For an order-$n$ signing $B$, put $M=M(B)$ and assign every projective spin
configuration the nonnegative integer weight

\[
 w_B([x])=\frac{M-|Q_B(x)|}{2}.
\]

Define its energy-weighted projective covering radius by

\[
 \rho_{\mathrm w}(B)
 =\max_{[b]\in\mathcal P_n}\min_{[x]\in\mathcal P_n}
 \bigl(d_\pm([b],[x])+w_B([x])\bigr).
\]

Then

\[
 \boxed{E(B)=M(B)+n-2\rho_{\mathrm w}(B).}
\]

Equivalently, with the nonnegative weighted covering deficit

\[
 \delta_{\mathrm w}(B)=\left\lfloor\frac n2\right\rfloor
 -\rho_{\mathrm w}(B),
\]

the original minimax sequence obeys the exact Bellman identity

\[
 \boxed{
 F(n+1)=\min_B\left(
 M(B)+(n\bmod2)+2\delta_{\mathrm w}(B)
 \right),}
\]

where the minimum is over all order-$n$ signings.

Let

\[
 X_B=\{[x]\in\mathcal P_n:|Q_B(x)|=M(B)\}
\]

be the projective extremizer code and let $\rho_{\rm ext}(B)$ be its ordinary
covering radius in $\mathcal P_n$.  If only the maximizing configurations are
retained in the extension constraints, their exact optimum is

\[
 E_{\rm ext}(B)=M(B)+n-2\rho_{\rm ext}(B).
\]

In particular,

\[
 \rho_{\mathrm w}(B)\le\rho_{\rm ext}(B),
 \qquad E(B)\ge E_{\rm ext}(B).
\]

#### Proof

Let $b\in\{\pm1\}^n$ be the incident-edge signs and let $t$ be the new
vertex spin.  For each fixed old spin vector $x$, pairing $t=1$ and $t=-1$
gives

\[
 \max_{t\in\{\pm1\}}
 \left|Q_B(x)+t\,b\mathbin\cdot x\right|
 =|Q_B(x)|+|b\mathbin\cdot x|.
\]

Moreover,

\[
 |b\mathbin\cdot x|=n-2d_\pm([b],[x]).
\]

All values of $Q_B$ have the parity of $\binom n2$, so $w_B([x])$ is an
integer.  Therefore, for fixed $b$,

\[
 \begin{aligned}
 \max_x\bigl(|Q_B(x)|+|b\mathbin\cdot x|\bigr)
 &=M+n-2\min_{[x]\in\mathcal P_n}
   \bigl(w_B([x])+d_\pm([b],[x])\bigr).
 \end{aligned}
\]

Minimizing the left side over $[b]$ changes the minimum on the right into a
maximum and proves the first boxed identity.  Since
$n-2\lfloor n/2\rfloor=n\bmod2$, the deficit form follows.  Minimizing
$E(B)$ over all $B$ and $b$ enumerates every signing of order $n+1$, proving
the Bellman identity.  Restricting the inner minimum to $w_B=0$ proves the
extremizer-code formula, and allowing all weighted configurations can only
decrease the covering radius. \(\square\)

There is an exact finite-width reduction inside this identity.  For every
center $[b]$, some extremizer has projective distance at most
$\lfloor n/2\rfloor$.  Hence a configuration with
$w_B([x])>\lfloor n/2\rfloor$ can never attain the inner minimum defining
$\rho_{\mathrm w}(B)$.  Therefore the weighted radius, and thus $E(B)$,
depends only on the near-ground window

\[
 \boxed{|Q_B(x)|\ge M(B)-2\left\lfloor\frac n2\right\rfloor.}
\]

The state is richer than the exact ground-state code but has energy width at
most $n$; configurations deeper in the landscape are rigorously irrelevant to
one-vertex extension.

The complete optimizer catalogue gives the following exact profiles.
Multiplicities again count switching-permutation classes, and each displayed
triple is $(\rho_{\rm ext},\rho_{\mathrm w},E)$:

\[
\begin{array}{c|c}
n&\text{profile distribution}\ \hline
2&(0,0,3)^{\times1}\\
3&(1,1,4)^{\times2}\\
4&(2,2,4)^{\times1}\\
5&(2,2,5)^{\times1}\\
6&(1,1,9)^{\times1}\\
7&(3,3,10)^{\times4},\ (3,2,12)^{\times2}\\
8&(3,3,12)^{\times2}\\
9&(4,4,13)^{\times4},\ (3,3,15)^{\times11}\\
10&(3,3,17)^{\times1},\ (2,2,19)^{\times1}
\end{array}
\]

Thus the extremizer code alone exactly distinguishes the successful and
failing order-9 classes and also explains the values 17 and 19 of the two
order-10 classes.  It is not a complete state in general: all six order-7
classes have extremizer radius 3, but the $|Q_B|=7$ layer lowers the weighted
radius from 3 to 2 in exactly the two classes whose extension value is 12.
This is a certified counterexample to the tempting claim that ground states
alone determine the cavity step.  The exact state is the energy-weighted
projective covering landscape.

### Proposition 16 (the scalar free-energy curve is not a Bellman state)

There are two optimal order-9 signings $B_+$ and $B_-$ having the identical
complete absolute-energy histogram

\[
 \#\{[x]\in\mathcal P_9:|Q_{B_\pm}(x)|=e\}
 =\begin{cases}
 60,&e=0,\\
 111,&e=4,\\
 60,&e=8,\\
 25,&e=12,
 \end{cases}
\]

but

\[
 \boxed{E(B_+)=13,\qquad E(B_-)=15.}
\]

Consequently, for every real $t$, both the absolute-value partition function
and the augmented partition function agree:

\[
 \sum_x e^{t|Q_{B_+}(x)|}
 =\sum_x e^{t|Q_{B_-}(x)|},
\]

\[
 \sum_{x,s}e^{tsQ_{B_+}(x)}
 =\sum_{x,s}e^{tsQ_{B_-}(x)}.
\]

#### Exact certificate

Take the root-normalized graph6 records `G?qmaw` and `GCpbaw`.  They lie in
different switching-permutation classes in the complete order-9 catalogue.
Direct spin enumeration gives the displayed histogram for each.  Their
extremizer covering radii are respectively 4 and 3, so Theorem 15 already
gives the two extension values; a separate enumeration of all incident sign
vectors recomputes 13 and 15.  Representative incident-edge masks are 142 and
106 respectively, with bit $i$ equal to one for a negative new edge.  The only
completeness dependency is the stated nauty catalogue boundary. \(\square\)

Thus even the full scalar partition-function curve, equivalently every
absolute-energy level count, loses the geometry needed for a cavity step.  It
does not refute the minimax free-energy limit program: a nonlocal interpolation
or a richer order parameter may still prove that limit.  It does rule out a
Bellman recursion whose state is only the optimized scalar free energy, even
if that scalar is retained at every temperature.

The verifier also exhausts every root-normalized residual graph, not only
optimal classes, through residual order 8.  Let

\[
 \mathcal B_n=\operatorname{ParetoMin}
 \{(M(B),\delta_{\mathrm w}(B)):B\text{ has order }n\},
\]

where both coordinates are minimized.  The complete exact frontiers are

\[
\begin{array}{c|c}
n&\mathcal B_n\\ \hline
2&\{(1,1)\}\\
3&\{(3,0)\}\\
4&\{(4,0)\}\\
5&\{(4,0)\}\\
6&\{(5,2),(7,1),(9,0)\}\\
7&\{(9,0)\}\\
8&\{(10,1),(12,0)\}.
\end{array}
\]

Every point in $\mathcal B_6$ minimizes the Bellman objective and gives
$F(7)=9$; both points in $\mathcal B_8$ give $F(9)=12$.  Thus a two-unit
increase in internal maximum can be exactly offset by a one-unit improvement
in weighted covering deficit.  This is positive finite evidence for using a
near-optimal Pareto family rather than exact minimizers, while making no
asymptotic claim.  Completeness again trusts the asserted `geng` counts and
stream hashes; every weighted-radius value is independently checked by direct
enumeration of all incident sign vectors.

### Theorem 17 (nonlinear Gaussian-sign lower bound)

Let $A$ be any order-$n$ signing.  For $s\ge0$ define

\[
 D_s=1+s^2(n-1),\qquad d_s=\frac{2s}{D_s},\qquad
 z_{ij,s}=\frac{s^2|(A^2)_{ij}|}{D_s}.
\]

Then the complete arcsine contribution in the square-covariance construction
gives

\[
 \boxed{
 M(A)\ge \frac1\pi\sum_{i<j}
 \left[
 \arcsin(z_{ij,s}+d_s)-\arcsin(z_{ij,s}-d_s)
 \right].}
 \tag{17.1}
\]

In particular, taking $s=1/\sqrt{n-1}$ gives

\[
 \boxed{
 M(A)\ge \frac{n(n-1)}\pi
 \arcsin\frac1{\sqrt{n-1}}.}
 \tag{17.2}
\]

For $n\ge3$ one has the stronger stability inequality

\[
 \boxed{
 M(A)\ge \frac{n(n-1)}\pi
 \arcsin\frac1{\sqrt{n-1}}
 +\frac{\operatorname{tr}(A^4)-n(n-1)^2}
 {8\pi(n-1)(n-2)^{3/2}}.}
 \tag{17.3}
\]

#### Proof

The matrices

\[
 R_s^\pm=\frac{(I\pm sA)^2}{D_s}
\]

are correlation matrices: they are positive semidefinite and their diagonal
entries are one.  Let $X^\pm$ be the coordinatewise signs of centered
Gaussian vectors with these correlations.  For $i\ne j$ put

\[
 u_{ij}=\frac{s^2a_{ij}(A^2)_{ij}}{D_s}.
\]

Multiplying the two off-diagonal correlations by $a_{ij}$ gives
$u_{ij}\pm d_s$.  The arcsine law and oddness of arcsine therefore yield

\[
 \mathbb E Q_A(X^+)-\mathbb E Q_A(X^-)
 =\frac2\pi\sum_{i<j}
 \bigl[\arcsin(u_{ij}+d_s)-\arcsin(u_{ij}-d_s)\bigr].
\]

The bracketed function is even in $u_{ij}$, so it equals the expression with
$z_{ij,s}=|u_{ij}|$.  Both expectations lie in $[-M(A),M(A)]$; division by
two proves (17.1).

For fixed $d\ge0$, set

\[
 f_d(z)=\arcsin(z+d)-\arcsin(z-d).
\]

On the correlation domain, $f_d$ is even and convex: for $z\ge0$,

\[
 f_d''(z)=h(z+d)-h(z-d)\ge0,
 \qquad h(t)=\frac{t}{(1-t^2)^{3/2}},
\]

because $h$ is increasing.  Hence $f_d(z)\ge f_d(0)=2\arcsin d$.
At $s=1/\sqrt{n-1}$, $D_s=2$ and $d_s=1/\sqrt{n-1}$; summing over the
$\binom n2$ edges proves (17.2).

For the quantitative term put $q=n-1$ and $d=1/\sqrt q$.  Since

\[
 h'(t)=\frac{1+2t^2}{(1-t^2)^{5/2}}
\]

is even and increasing in $|t|$, $f_d''(z)$ is increasing for $z\ge0$.
Taylor's theorem therefore gives

\[
 f_d(z)\ge f_d(0)+\frac12f_d''(0)z^2
 =2\arcsin d+\frac{q}{(q-1)^{3/2}}z^2.
\]

Here $z_{ij}=|(A^2)_{ij}|/(2q)$, and

\[
 \sum_{i<j}z_{ij}^2
 =\frac{\operatorname{tr}(A^4)-nq^2}{8q^2}.
\]

Substitution in (17.1) proves (17.3). $\square$

The excess has the exact sum-of-squares form

\[
 \operatorname{tr}(A^4)-n(n-1)^2
 =\|A^2-(n-1)I\|_F^2.
 \tag{17.4}
\]

Thus any sequence with
$M(A_n)=(1/\pi+o(1))n^{3/2}$ must obey
$\|A_n^2-(n-1)I\|_F=o(n^2)$.  This is a necessary conference-like stability
condition, not an existence theorem for such a sequence.  The universal gain
in (17.2) over the earlier linearized bound is only $\Theta(\sqrt n)$ and does
not change the leading constant.  It does have exact finite consequences.
For $n=21$, (17.2) is strictly greater than 30 and every quadratic energy is
even, so

\[
 \boxed{F(21)\ge32.}
\]

For a rational certificate of strictness, use $\pi/14<11/49$ and
$\sin t<t-t^3/6+t^5/120$ for $t>0$; direct rational squaring gives

\[
 \left(\frac{11}{49}-\frac{(11/49)^3}{6}
 +\frac{(11/49)^5}{120}\right)^2<\frac1{20}.
\]

Hence $\sin(\pi/14)<1/\sqrt{20}$ and
$\arcsin(1/\sqrt{20})>\pi/14$, which is exactly the required strict
inequality.  The old linearized bound is less than 30 and parity-rounds only
to 30.

### Theorem 18 (multivertex weighted Bellman identity)

Fix internal signings $B,D$ of orders $n,k$, and set

\[
 J(B,D)=\min_{C\in\{\pm1\}^{n\times k}}
 M\begin{pmatrix}B&C\\ C^{\mathsf T}&D\end{pmatrix},
 \qquad
 L(B,D)=\max_{x,y}|Q_B(x)+Q_D(y)|.
\]

Let $\mathcal P_{nk}=\{\pm1\}^{n\times k}/\{C\sim-C\}$ and define the
projective rank-one code

\[
 \mathcal R_{n,k}=\{[xy^{\mathsf T}]:x\in\{\pm1\}^n,
 y\in\{\pm1\}^k\}\subseteq\mathcal P_{nk}.
\]

Assign the well-defined integer weight

\[
 w_{B,D}([xy^{\mathsf T}])
 =\frac{L(B,D)-|Q_B(x)+Q_D(y)|}{2}
\]

and put

\[
 \rho_{\mathrm w}^{\square}(B,D)
 =\max_{[C]\in\mathcal P_{nk}}
 \min_{R\in\mathcal R_{n,k}}
 \bigl[d_\pm([C],R)+w_{B,D}(R)\bigr].
\]

Then

\[
 \boxed{J(B,D)=L(B,D)+nk-2\rho_{\mathrm w}^{\square}(B,D).}
 \tag{18.1}
\]

Equivalently, if

\[
 \delta_{\mathrm w}^{\square}(B,D)
 =\left\lfloor\frac{nk}{2}\right\rfloor
 -\rho_{\mathrm w}^{\square}(B,D),
\]

then

\[
 \boxed{
 F(n+k)=\min_{B,D}\left[
 L(B,D)+(nk\bmod2)+2\delta_{\mathrm w}^{\square}(B,D)
 \right].}
 \tag{18.2}
\]

#### Proof

For fixed $x,y$, pair $(x,y)$ with $(x,-y)$.  The internal energy is
unchanged and the cross energy changes sign, so

\[
 \max_{\eta=\pm1}
 |Q_B(x)+Q_D(y)+\eta x^{\mathsf T}Cy|
 =|Q_B(x)+Q_D(y)|+|x^{\mathsf T}Cy|.
\]

Moreover

\[
 |x^{\mathsf T}Cy|=nk-2d_\pm([C],[xy^{\mathsf T}]).
\]

Writing the internal absolute energy as $L-2w$ proves (18.1) after maximizing
over the rank-one code and minimizing over $[C]$.  The rank-one code contains
an exact internal maximizer of weight zero, so its distance from every center
is at most $\lfloor nk/2\rfloor$; hence the deficit is nonnegative.  Minimizing
over $B,D$ enumerates every order-$(n+k)$ signing and proves (18.2). $\square$

Exactly as in the one-vertex case, only the finite window

\[
 |Q_B(x)+Q_D(y)|
 \ge L(B,D)-2\left\lfloor\frac{nk}{2}\right\rfloor
\]

can affect the weighted radius.  Theorem 18 is the exact multivertex Bellman
state; by itself it does not estimate that state asymptotically.

### Proposition 19 (weighted entropy upper bound)

Define

\[
 \Xi(B,D)=\log\sum_{R\in\mathcal R_{n,k}}
 \exp\left(-\frac{2w_{B,D}(R)^2}{nk}\right).
\]

Then

\[
 \boxed{
 J(B,D)\le L(B,D)
 +\sqrt{2nk\bigl(\Xi(B,D)+\log4\bigr)}.}
 \tag{19.1}
\]

#### Proof

Choose the entries of $C$ independently and uniformly from $\{\pm1\}$.
For fixed rank-one $R$, Hoeffding's inequality gives, for $T\ge0$,

\[
 \Pr\{|\langle C,R\rangle|>T+2w(R)\}
 \le2\exp\left(-\frac{(T+2w(R))^2}{2nk}\right).
\]

After dropping the nonnegative cross term in the square and summing over
$R$, the union probability is at most

\[
 2\exp\left(-\frac{T^2}{2nk}+\Xi(B,D)\right).
\]

With $T=\sqrt{2nk(\Xi+\log4)}$ this is at most $1/2$.  Some $C$ therefore
satisfies $|\langle C,R\rangle|\le T+2w(R)$ for every $R$.  For that $C$,
$L-2w(R)+|\langle C,R\rangle|\le L+T$, and Theorem 18 proves (19.1).
$\square$

Since $|\mathcal R_{n,k}|=2^{n+k-2}$, the generic estimate
$\Xi\le(n+k-2)\log2$ recovers the known leading
$\sqrt{nk(n+k)}$ cross-block cost.  Thus the exact identity does not make the
ordinary block wall disappear.  A concrete sufficient profile theorem would
instead find $B,D$ for every large $n,k$ such that, for some $\varepsilon>0$,

\[
 L(B,D)+\sqrt{2nk(\Xi(B,D)+\log4)}
 \le\bigl(F(n)^{2/3}+F(k)^{2/3}\bigr)^{3/2}
 +O((n+k)^{3/2-\varepsilon}).
 \tag{19.2}
\]

It would imply power-saving near-subadditivity of $F^{2/3}$ and hence
convergence by Theorem 13.  Equation (19.2) is a settling target, not a proved
estimate.

### Proposition 20 (density-one control of Bellman-optimal predecessors)

For each $n$, choose $B_n$ attaining the minimum in Theorem 15 and put

\[
 \sigma_n=M(B_n)-F(n),\qquad
 \delta_n=\delta_{\mathrm w}(B_n).
\]

Then exactly

\[
 F(n+1)-F(n)=\sigma_n+(n\bmod2)+2\delta_n.
 \tag{20.1}
\]

Since $F(n)=O(n^{3/2})$, telescoping (20.1) gives

\[
 \sum_{n=N}^{2N}(\sigma_n+2\delta_n)=O(N^{3/2}).
 \tag{20.2}
\]

Consequently, for every function $g(N)\to\infty$, all but $O(N/g(N))$
indices $n\in[N,2N]$ admit a Bellman-optimal predecessor satisfying

\[
 \sigma_n\le g(N)\sqrt N,
 \qquad
 \delta_n\le\frac12g(N)\sqrt N.
 \tag{20.3}
\]

#### Proof

Equation (20.1) is Theorem 15 with $F(n)$ subtracted.  Both $\sigma_n$ and
$\delta_n$ are nonnegative.  Sum from $N$ to $2N$ and use the elementary
random-sign upper bound $F(t)=O(t^{3/2})$; the parity sum is only $O(N)$.
Markov's inequality applied to the nonnegative quantities
$\sigma_n+2\delta_n$ proves (20.3). $\square$

This is real asymptotic control of the exact Bellman state on a density-one
set.  It does not control sparse exceptional orders, which may still have
linear-size deficits, and therefore does not prove convergence.

### Complete order-9 weighted-geometry computation

A new exhaustive pass over all 12,346 root-normalized order-9 signings gives

\[
 \boxed{\mathcal B_9=\{(12,0)\}.}
\]

Among the 55 root-normalized records with $M=12$, 20 have
$\delta_{\mathrm w}=0$ and 35 have $\delta_{\mathrm w}=1$.  They form exactly
15 switching-permutation classes, split as four and eleven respectively.

There is also a stronger finite collision.  The graph6 records `GHOgmo` and
`Gxd?Dc` both have $M=14$, the same complete projective histogram

\[
 \#\{|Q|=2,6,10,14\}=(124,85,37,10),
\]

and the same ordered projective pair-distance enumerator of their ten exact
maximizers,

\[
 (N_0,N_1,N_2,N_3,N_4)=(10,16,14,20,40).
\]

Nevertheless their profiles $(\rho_{\rm ext},\rho_{\mathrm w},E)$ are
respectively $(4,4,15)$ and $(3,3,17)$.  Therefore even the complete energy
histogram together with the two-point overlap law of the exact maximizers
does not determine one-vertex extension.

Across the complete rooted catalogue, ten energy-histogram groups containing
874 records mix different weighted deficits.  Adding the exact-maximizer
pair-distance law leaves three mixed groups containing 112 records.  The
complete ordered energy-coloured two-point distance distribution

\[
 T_B(e,f,r)=\#\{([x],[y]):|Q_B(x)|=e, |Q_B(y)|=f,
 d_\pm([x],[y])=r\}
\]

separates every weighted deficit at order 9.  This last statement is an exact
finite catalogue result, not a claim that $T_B$ is sufficient at every order
or asymptotically.  The new verifier checks the catalogue count and nauty
stream digest, independently samples the graph6 decoder with NetworkX,
checks all weighted radii against direct extension enumeration, and includes
the two explicit collision records.

### Consequence for the convergence program

The finite obstruction does not refute asymptotic near-subadditivity: the
required internal sacrifice is only two units at order 8. It does refute any
proof that composes arbitrary exact minimizers or records only their scalar
maximum. The ground-state state space must retain a composability profile, or
allow subleading internal slack before optimizing the cross block.

A sharper settling target is therefore the following. Find
$r_n=o(n^{3/2})$, $C<\infty$, and $\delta>0$ such that

\[
 K_{n,m}\bigl(F(n)+r_n,F(m)+r_m\bigr)^{2/3}
 \le F(n)^{2/3}+F(m)^{2/3}
 +C(n+m)^{1-\delta}.
\]

Since $F(n+m)\le K_{n,m}(u,v)$ for every $u,v$, Theorem 13 would then prove
convergence. This target permits the exact kind of internal sacrifice forced
by the $(2,8)$ example while still requiring it to vanish on the normalized
scale.

Theorem 15 supplies a second, one-vertex state space.  The scalar $M(B)$ is
insufficient, while the Pareto data

\[
 \bigl(M(B),\delta_{\mathrm w}(B)\bigr)
\]

is exact for the next Bellman step.  A useful asymptotic theorem would control
this weighted covering deficit uniformly over the near-optimal energy window,
rather than only the number or geometry of exact maximizers.  The order-7
counterexample proves that an extremizer-only code theorem cannot provide
such control.

### Order-11 exploration boundary

The direct order-10 extension calculation is complete and proves
$F(11)\le17$. A separate cutting-plane exploration completed all fifteen
$(2,9)$ switching-class pairs, with values 17 in six cases and 19 in nine.
One $(3,8)$ pair certified value 17. The next $(3,8)$ lower-bound query did
not complete within the 120-second per-call policy and the exploratory run
was manually stopped before the remaining $(3,8)$, $(4,7)$, and $(5,6)$
pairs. No result from an uncompleted pair is used, and no exact value of
$F(11)$ is inferred. This is a solver boundary, not evidence that the missing
pairs are mathematically harder or have larger values.

## 8. Current best continuation targets

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

The ground-state alternative is now the Pareto-profile refinement following
Theorem 14. It seeks Theorem 13's power-saving composition bound after
allowing $o(n^{3/2})$ internal slack. Exact minimizers cannot be required:
$K_{2,8}(1,10)=15$ even though $F(10)=13$, whereas
$K_{2,8}(1,12)=13$. Current block constructions still fail because internal
block energy and the cross field can cancel for the same spin assignment,
while bounding the two pieces separately incurs a leading $N^{3/2}$ cost.
Theorems 11 and 12 give hereditary lower information but no such
optimizer-compatible upper construction.

Theorem 18 and Proposition 19 make this alternative more precise.  The exact
cross-block state is the weighted covering landscape of the projective
rank-one code, and the explicit constant-matching target is (19.2).  The
generic entropy estimate still has leading cost, so the new result isolates
the missing profile law rather than proving it.  Proposition 20 controls the
one-vertex state on a density-one set of orders, but a convergence argument
must still eliminate or bridge the sparse exceptional orders.

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

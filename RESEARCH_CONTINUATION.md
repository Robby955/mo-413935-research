# Continued attack after the second research note

Status date: 2026-08-02
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
10&2&17^{\times1},\ 19^{\times1}&F(11)=17
\end{array}
\]

Thus an arbitrary optimum need not admit an optimum one-vertex extension.
At order 7, the two failing classes first become infeasible when the
configurations with $|Q_B(x)|=7$ are added to the maximizing configurations.
At order 9, the maximizing configurations with $|Q_B(x)|=12$ already make
extension at value 13 infeasible in all eleven failing classes.
At order 10 the two exact optimum classes have extension values 17 and 19;
the better class gives the upper bound $F(11)\le17$. The later complete
order-11 certificate proves equality. A direct upper witness uses order-10 graph6 record
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
$\|A_n^2-(n-1)I\|_F=o(n^2)$.  This is vanishing normalized Gram defect, not a
conference-rigidity condition or an existence theorem.  The universal gain
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
convergence by Theorem 13.  This was initially recorded as a settling target.
Proposition 23 below proves that (19.2) is in fact impossible for every
$B,D$: the left side has a larger universal leading constant.  It is retained
here to make the failed route and its correction auditable.

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
set.  Theorem 24 below shows that sparse exceptional orders are not the
decisive scalar issue: even a uniform $O(\sqrt n)$ scalar Bellman cost is
compatible with nonconvergence.  A proof must stabilize the normalized
Bellman cost using matrix geometry not present in Proposition 20.

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

## 8. New corrections and structural results

The next results close two apparent routes from the preceding section and
replace them with narrower targets.  They do not settle the original limit.

### Theorem 21 (parity of the Gram defect)

For an order-$n$ signing put

\[
 \Delta(A)=\|A^2-(n-1)I\|_F^2
 =\sum_{i\ne j}(A^2)_{ij}^2.
\]

Then

\[
 \boxed{
 \Delta(A)\ge
 \begin{cases}
 n(n-1),&n\text{ odd},\\
 2n(n-2),&n\equiv0\pmod4,\\
 0,&n\equiv2\pmod4.
 \end{cases}}
 \tag{21.1}
\]

The three bounds are sharp on infinite Paley-derived subsequences: use a
symmetric conference matrix itself, or delete one or two vertices.

#### Proof

If $n$ is odd, every off-diagonal entry of $A^2$ is a sum of the odd number
$n-2$ of signs.  It is therefore a nonzero odd integer, and summing its
square over the $n(n-1)$ ordered pairs proves the first line.

Suppose $n\equiv0\pmod4$ and write $a_{ij}=(-1)^{e_{ij}}$.  Put
$d_i=\sum_{j\ne i}e_{ij}\pmod2$.  Direct reduction modulo four gives

\[
 \frac{(A^2)_{ij}}2\equiv1+d_i+d_j\pmod2.
 \tag{21.2}
\]

Thus $|(A^2)_{ij}|\ge2$ whenever $d_i=d_j$.  If $r$ of the $d_i$ are one,
the handshake lemma makes $r$ even, and

\[
 \Delta(A)\ge
 8\left[\binom r2+\binom{n-r}2\right]\ge2n(n-2).
\]

There is no positive universal bound when $n\equiv2\pmod4$, because a
symmetric conference matrix has $\Delta=0$.

For sharpness in the other two cases, take a conference matrix $C$ of order
$L$.  If one vertex is deleted, the retained principal block $A$ has every
off-diagonal entry of $A^2$ equal to a sign, so
$\Delta=n(n-1)$.  If two vertices are deleted, write the deleted columns as
an $n$ by $2$ matrix $U$.  The two columns of $U$ are orthogonal, so half of
its row products are $1$ and half are $-1$.  Hence precisely
$n(n-2)/2$ ordered row pairs have inner product of magnitude two, giving
$\Delta=2n(n-2)$.  Paley conference orders supply infinitely many examples.
$\square$

Combining (21.1) with (17.3) yields a second exact finite improvement:

\[
 \boxed{F(20)\ge30.}
 \tag{21.3}
\]

Indeed the right side of (17.3), using $\Delta\ge720$, is

\[
 \frac{380}{\pi}\arcsin\frac1{\sqrt{19}}
 {}+\frac{90}{19\pi\,18^{3/2}}>28.
\]

This strict inequality has a rational certificate.  Put $x=39/170$.
Then $x<1/\sqrt{19}$, $1/\sqrt{18}>70/297$,
$\arcsin x>x+x^3/6$, and $\pi<22/7$.  After multiplying by $\pi$, the
remaining margin is exactly

\[
 380\left(x+\frac{x^3}{6}\right)
 {}+\frac5{19}\frac{70}{297}-88
 =\frac{8798941}{2772405900}>0.
\]

Order-20 energies are even, so a strict lower bound above $28$ rounds to
$30$.  The earlier linear Gaussian bound rounds only to $28$.

### Theorem 22 (optimized nonlinear certificate and its sharp barrier)

Let $q=n-1$ and

\[
 \varepsilon(A)=\frac{\Delta(A)}{nq^3}.
\]

For every $0\le w<1$, define

\[
 a=\frac{1+w}{2}\sqrt{\varepsilon(A)},\qquad
 d=\frac{\sqrt{1-w^2}}{\sqrt q}.
\]

Then the full nonlinear bound (17.1) implies

\[
 \boxed{
 M(A)\ge\frac{nq}{2\pi}
 \left[\arcsin(a+d)-\arcsin(a-d)\right].}
 \tag{22.1}
\]

Consequently, if $\varepsilon(A_n)\to\varepsilon<1$, then

\[
 \boxed{
 \liminf_n\frac{M(A_n)}{n^{3/2}}
 \ge\frac1\pi\sqrt{\frac{2}{1+\sqrt{1-\varepsilon}}}.}
 \tag{22.2}
\]

#### Proof

Set $t=qs^2$ and $w=(t-1)/(t+1)$.  Negative $w$ is dominated by $-w$,
because it has the same $d$ and smaller nonnegative offsets.  For
$0\le w<1$,

\[
 d_s=d,\qquad z_{ij,s}=\frac{1+w}{2q}|(A^2)_{ij}|.
\]

For fixed $d$, let $f_d(z)=\arcsin(z+d)-\arcsin(z-d)$ and
$G(u)=f_d(\sqrt u)$.  The function $G$ is convex on the correlation domain.
Indeed, with $z=\sqrt u$,

\[
 G''(u)=\frac{zf_d''(z)-f_d'(z)}{4z^3}.
\]

If $h(v)=v(1-v^2)^{-3/2}$, then
$f_d''(z)=h(z+d)-h(z-d)$.  Moreover
$f_d'''(z)=h'(z+d)-h'(z-d)\ge0$, because $h'$ is even and increases with
absolute value.  Thus
$P(z)=zf_d''(z)-f_d'(z)$ satisfies $P(0)=0$ and
$P'(z)=zf_d'''(z)\ge0$.  Continuity handles $z=0$.

The mean of $z_{ij,s}^2$ over the $\binom n2$ edges is $a^2$, since
$2\sum_{i<j}(A^2)_{ij}^2=\Delta(A)$.  Thus $a$ is the root-mean-square of
the offsets and $a\le\max z_{ij,s}$, so $|a\pm d|\le1$.  Jensen's
inequality proves (22.1), with boundary cases obtained by continuity.
For fixed $w$ and $n\to\infty$, its normalized right side tends to

\[
 \frac1\pi\,
 \frac{\sqrt{1-w^2}}
 {\sqrt{1-\varepsilon(1+w)^2/4}}.
\]

Optimizing the ratio squared,

\[
 \frac{1-w^2}{1-\varepsilon(1+w)^2/4},
\]

gives $1+w=2/(1+\sqrt{1-\varepsilon})$ and yields (22.2).
$\square$

There is also a sharp limitation on this entire one-parameter mechanism.  Let
$B_s(A)$ denote the right side of (17.1), and put

\[
 G_n=\min_A\sup_{s\ge0}B_s(A).
\]

Then

\[
 \boxed{\lim_{n\to\infty}\frac{G_n}{n^{3/2}}=\frac1\pi.}
 \tag{22.4}
\]

The lower bound follows by taking $s=1/\sqrt q$.  For the upper bound, embed
order $n$ as a principal block of a Paley conference matrix of order
$n+r$, where $r=o(n)$.  If $U$ is the deleted block, then for $i\ne j$,
$|(A^2)_{ij}|=|(UU^{\mathsf T})_{ij}|\le r$.  Uniformly in $s$,

\[
 z_{ij,s}\le\frac rq=o(1),\qquad d_s\le\frac1{\sqrt q}.
\]

Writing the arcsine difference as an integral therefore gives

\[
 B_s(A)\le
 \frac{nq}{\pi\sqrt q}\,(1+o(1))
 =\left(\frac1\pi+o(1)\right)n^{3/2}
\]

uniformly in $s$.  This proves (22.4).  In particular, optimizing $s$,
retaining all the $A^2$ offsets, or taking ordinary convex combinations of
these certificates cannot improve the universal leading constant.

The earlier consequence
$\|A_n^2-(n-1)I\|_F=o(n^2)$ at constant $1/\pi$ should be described as
vanishing $n^{-2}$-normalized Gram defect, not as conference rigidity.
A uniform random signing has

\[
 \mathbb E\Delta(A)=n(n-1)(n-2)=\Theta(n^3)
\]

and therefore satisfies that weaker normalization in probability.

### Proposition 23 (the weighted-entropy target is impossible)

Let

\[
 E_{\rm ent}(B,D)
 =L(B,D)+\sqrt{2nk\bigl(\Xi(B,D)+\log4\bigr)}.
\]

For every pair of internal blocks,

\[
 \boxed{
 E_{\rm ent}(B,D)\ge\sqrt{2nk(n+k)\log2}.}
 \tag{23.1}
\]

#### Proof

There are $2^{n+k-2}$ projective rank-one words, and
$0\le w(R)\le L/2$.  Hence

\[
 \Xi\ge(n+k-2)\log2-\frac{L^2}{2nk}.
\]

Put $b=\sqrt{2nk(n+k)\log2}$.  If $L\le b$, then

\[
 E_{\rm ent}\ge L+\sqrt{b^2-L^2}\ge b;
\]

if $L>b$, the conclusion is immediate.  This proves (23.1). $\square$

For $n=k=N$, the unavoidable constant in (23.1) is
$2\sqrt{\log2}=1.665109\ldots$.  By the Paley upper bound, the right side of
the proposed target (19.2) is at most
$(\sqrt2+o(1))N^{3/2}$.  The gap

\[
 2\sqrt{\log2}-\sqrt2=0.250895\ldots
\]

is leading order.  Thus (19.2) is false for every choice of $B,D$, even with
an $O(N^{3/2-\eta})$ error.  Retaining the discarded cross term inside the
same Hoeffding union sum does not repair it: for a target $K\ge L$, that sum
is

\[
 U_K=2\sum_R
 \exp\left(-\frac{(K-h(R))^2}{2nk}\right),
 \qquad h(R)=|Q_B(x)+Q_D(y)|,
\]

and

\[
 U_K\ge2^{n+k-1}\exp\left(-\frac{K^2}{2nk}\right).
\]

Therefore this certificate cannot be below one unless
$K^2>2nk(n+k-1)\log2$.

The exact replacement keeps the dependence among the bad events.  Regard
$G=\mathcal R_{n,k}$ as a group under entrywise multiplication.  Fix a cross
seed $C_0$, put

\[
 c(U)=|\langle C_0,U\rangle|,\qquad
 h(R)=|Q_B(x)+Q_D(y)|,
\]

and let $C_g=C_0\odot g$ range over its row-and-column switching orbit.  Then

\[
 \boxed{
 J_{\rm orb}(C_0;B,D)
 =\min_{g\in G}\max_{R\in G}[h(R)+c(gR)].}
 \tag{23.2}
\]

For $K\ge L$, define $C_t=\{U:c(U)>t\}$.  The bad shifts are exactly

\[
 \boxed{\mathcal B_K=\bigcup_{R\in G}R\,C_{K-h(R)}.}
 \tag{23.3}
\]

Thus this orbit succeeds precisely when $\mathcal B_K\ne G$.  Equations
(23.2)--(23.3) turn the missing composition estimate into a weighted sumset,
or max-plus convolution, problem on $\mathbb F_2^{n+k-2}$.  Unlike (19.1),
it preserves all intersections among bad events.

The first nontrivial correction to a union bound is also exact.  Let $C$ have
independent uniform sign entries and let
$E_R=\{|\langle C,R\rangle|>u_R\}$.  If two rank-one words have projective
distance $d$, orient them to have ordinary distance $d$ and put $N_0=nk$.
Then

\[
\begin{aligned}
 p_{N_0,d}(u,v)
 =2^{-N_0}
 \sum_{a=0}^{N_0-d}\sum_{b=0}^{d}
 &\binom{N_0-d}{a}\binom db\\
 {}\times\mathbf1\{&
 |N_0-2(a+b)|>u,\\
 &|N_0-2d-2a+2b|>v\}.
\end{aligned}
 \tag{23.4}
\]

This is $\Pr(E_R\cap E_S)$: $a$ and $b$ count negative cross signs on the
coordinates where $R,S$ agree and disagree.  Consequently Hunter's
spanning-tree inequality gives, for every spanning tree $\mathcal T$ on the
rank-one words,

\[
 \Pr\left(\bigcup_RE_R\right)
 \le\sum_R\Pr(E_R)
 -\sum_{\{R,S\}\in\mathcal T}
 p_{N_0,d(R,S)}(u_R,u_S).
 \tag{23.5}
\]

Thus the full labelled energy-coloured distance graph supplies the edge
weights for a second-order inclusion--exclusion certificate.  The aggregate
two-point invariant from the order-nine experiment motivates this data but
need not determine the best spanning tree.  No asymptotic estimate strong
enough for composition has yet been proved from (23.5).

An exact $4+4$ diagnostic shows why higher dependence is genuinely relevant.
For

\[
 B=\begin{pmatrix}
 0&1&1&1\\
 1&0&1&-1\\
 1&1&0&1\\
 1&-1&1&0
 \end{pmatrix},
\]

the internal profile over the 64 projective rank-one states is

\[
 \#\{h=0,2,4,6,8\}=(14,24,16,8,2),
\]

and exhaustive cross-centre enumeration gives $J(B,B)=10$.  The exact iid
union sums at targets $10,12,14$ are respectively

\[
 \frac{17973}{4096},\qquad
 \frac{6073}{4096},\qquad
 \frac{1653}{4096};
\]

the first-moment argument certifies only $14$.  For an adaptive optimal seed,
the switching orbit has two good projective shifts at target $10$ even though
the mean number of violated constraints is $17/4$.  This finite example
rules out treating the first moment as a proxy for noncoverage.

### Theorem 24 (scalar Bellman control still does not force convergence)

Put

\[
 r_n=\sigma_n+2\delta_n
 =F(n+1)-F(n)-(n\bmod2).
 \tag{24.1}
\]

There is an explicit integer sequence $f(n)$ which agrees with the certified
values through order ten and satisfies all of the following:

1. $f(n+1)-f(n)\equiv n\pmod2$ and
   $0\le f(n+1)-f(n)\le n$;
2. $f(n)\ge n\sqrt{n-1}/\pi$ and $f(n)\le\binom n2$;
3. both scalar hereditary consequences of Theorems 11 and 12;
4. a numerically admissible scalar Bellman decomposition with $\sigma_n=0$ and
   $0\le r_n\le3\sqrt n$ at every $n\ge10$;
5. $\liminf f(n)/n^{3/2}=0.36$ and
   $\limsup f(n)/n^{3/2}=0.44$.

Hence even uniform $O(\sqrt n)$ control of the scalar Bellman cost is
insufficient.  This countermodel is not claimed to arise from sign matrices;
it proves that additional matrix geometry is logically necessary.

#### Construction and proof

Let $T=e^2$ and

\[
 g(t)=t^{3/2}\left[
 \frac25+\frac1{25}\sin\bigl(\log\log(t+T)\bigr)\right].
\]

Use the exact prefix

\[
 f(1),\ldots,f(10)=0,1,3,4,4,5,9,10,12,13.
\]

For $n\ge10$, let $d_n$ be the least positive integer congruent to
$n$ modulo two and at least $g(n+1)-g(n)$, and set
$f(n+1)=f(n)+d_n$.  For $t\ge10$,

\[
 0.52\sqrt t\le g'(t)\le0.68\sqrt t.
\]

Indeed the oscillatory derivative contributes at most
$\frac1{50}\sqrt t$, because
$t/((t+T)\log(t+T))\le1/2$.  It follows that

\[
 0\le d_n-[g(n+1)-g(n)]<2,\qquad d_n\le3\sqrt n\le n.
\]

Set $\sigma_n=0$ and
$\delta_n=[d_n-(n\bmod2)]/2$.  Parity makes $\delta_n$ a nonnegative
integer, while $d_n\le n$ gives
$\delta_n\le\lfloor n/2\rfloor$, the admissible weighted-deficit range.
This proves the parity, increment, and Bellman assertions.  Summation gives
$f(n)=13+g(n)-g(10)+O(n)$, proving the two normalized subsequential limits.
Indeed $\phi(t)=\log\log(t+e^2)$ is continuous, increasing, and unbounded,
while $\phi(n+1)-\phi(n)\to0$; integer subsequences therefore approach the
phases $\pi/2$ and $3\pi/2$ modulo $2\pi$.
The lower envelope $0.36n^{3/2}+13-0.44\,10^{3/2}$ exceeds
$n^{3/2}/\pi$ for $n\ge10$; the finite prefix is checked directly.
The upper bound follows inductively from $d_n\le n$.

For the squared hereditary inequality, note that $f(t)\ge t-1$.  If
$s\ge10$, every tail increment is positive, so
$f(s+k)-f(s)\ge k$.  If $s<10\le s+k$, then
$f(s+k)\ge s+k+3$ while the exact prefix has $f(s)\le s+3$, giving the
same inequality.  Cases with $s+k<10$ are direct.  Hence

\[
 f(s+k)^2-f(s)^2
 \ge k(2s+k-2)
 \ge\binom k2+k(s\bmod2).
\]

For Theorem 12, $f(s+k)-f(s)\ge k\ge\mu_k$ when $s$ is odd, again with the
finite prefix checked directly.  This completes the construction.

The exact one-vertex diagnostic is instead

\[
 A_N=\frac1{N^{3/2}}\sum_{n=N}^{2N-1}r_n.
 \tag{24.2}
\]

Then $F(n)/n^{3/2}$ converges if and only if $A_N$ converges.  If
$A_N\to\Lambda$, the limit is

\[
 \boxed{\frac{\Lambda}{2^{3/2}-1}.}
 \tag{24.3}
\]

Indeed, with $c_N=F(N)/N^{3/2}$ and $q=2^{-3/2}$, telescoping (24.1) gives

\[
 c_{2N}=q\,c_N+q\,A_N+O(N^{-1/2}).
\]

The implication from $c_N$ to $A_N$ is immediate.  Conversely, for arbitrary
$n$ put $m=\lfloor n/2\rfloor$.  The even recurrence above and
$F(2m+1)-F(2m)=O(m)$ give

\[
 c_n=q[c_m+\Lambda]+o(1).
\]

If $C=\Lambda/(2^{3/2}-1)$ and $e_n=c_n-C$, boundedness of $c_n$ yields
$\limsup|e_n|\le q\limsup|e_n|$, hence $\limsup|e_n|=0$.  This proves
all-order convergence.  A useful sufficient target is the Bellman--Cesaro law

\[
 \frac1N\sum_{n\le N}\frac{r_n}{\sqrt n}\longrightarrow a.
 \tag{24.4}
\]

Abel summation then gives
$\sum_{n<N}r_n=(2a/3)N^{3/2}+o(N^{3/2})$, and therefore
$F(N)/N^{3/2}\to2a/3$.

### Complete order-10 finite-temperature phase diagram

For an order-10 signing let $H=(h_0,\ldots,h_{22})$ be its projective
absolute-energy histogram, where $h_j$ counts states with
$|Q_A|=2j+1$.  Put

\[
 z=4\sinh^2t,\qquad R_0=1,\quad R_1=z+1,\quad
 R_{j+1}=(z+2)R_j-R_{j-1}.
\]

Then
$\cosh((2j+1)t)=\cosh(t)R_j(z)$, so the augmented partition function is the
common factor $4\cosh t$ times

\[
 P_H(z)=\sum_jh_jR_j(z).
\]

Exhaustion of all 274,668 root-normalized order-10 records gives 6,012
distinct histograms.  Every resulting polynomial coefficientwise dominates
at least one of the following three:

\[
\begin{array}{c|l|l}
&\text{nonzero histogram entries}&P(z)\text{ coefficients, low to high}\\
\hline
P_0&(h_0,h_3,h_4,h_7)=(180,180,140,12)
 &(512,2816,4512,3680,2120,792,156,12)\\
P_1&(108,88,96,84,60,48,24,4)
 &(512,2816,5280,4704,2232,576,76,4)\\
P_2&(h_1,h_2,h_5,h_6)=(200,192,80,40)
 &(512,2816,5792,5600,2520,520,40).
\end{array}
\]

Here the unlabelled $P_1$ tuple lists $h_0,\ldots,h_7$ in order.
The first phase contains one rooted conference record and has $M=15$; the
second has $M=15$; the third consists of ground-state records with $M=13$.
Exactly

\[
 P_0-P_1=8z^2(z-2)(z+1)(z+3)(z+4)^2,
\]

\[
 P_1-P_2=4z^2(z+4)^2(z^3+z^2-10z-8).
\]

The cubic has a unique positive root
$\zeta=3.083872\ldots$.  Therefore the minimizing histogram is $P_0$ for

\[
 0<t<\operatorname{arsinh}(1/\sqrt2),
\]

$P_1$ until

\[
 t=\operatorname{arsinh}(\sqrt\zeta/2)=0.792460\ldots,
\]

and $P_2$ thereafter, with adjacent ties at the two thresholds.  This is an
exact finite catalogue theorem, conditional only on the asserted nauty stream
completeness; it is not an asymptotic phase-transition claim.  It proves that
the outer minimizer genuinely changes with temperature and that neither the
high-temperature conference optimizer nor the ground-state optimizer alone
controls the full minimax free-energy curve.

### Consequence for the convergence program

The finite obstruction does not refute asymptotic near-subadditivity: the
required internal sacrifice is only two units at order 8. It does refute any
proof that composes arbitrary exact minimizers or records only their scalar
maximum. The ground-state state space must retain a composability profile, or
allow subleading internal slack before optimizing the cross block.

A sharper settling target is therefore the following.  Define

\[
 K_{n,m}(u,v)=
 \min_{\substack{M(B)\le u\\M(D)\le v}}J(B,D).
\]

Find nonnegative $\eta_n=o(n^{3/2})$, $C<\infty$, and $\delta>0$ such that

\[
 K_{n,m}\bigl(F(n)+\eta_n,F(m)+\eta_m\bigr)^{2/3}
 \le F(n)^{2/3}+F(m)^{2/3}
 +C(n+m)^{1-\delta}.
\]

For these feasible budgets,
$F(n+m)\le K_{n,m}(F(n)+\eta_n,F(m)+\eta_m)$, so Theorem 13 would then prove
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

### Exact computer-assisted values at orders 11 and 12

The order-10 extension catalogue already gives a solver-free intermediate
improvement.  If an order-11 signing had maximum 13, deleting any vertex
would give an order-10 signing of maximum at most 13, hence an optimum because
$F(10)=13$.  Every optimum order-10 class has minimum extension value 17 or
19, a contradiction.  Energy parity therefore gives the rigorous catalogue
consequence $F(11)\ge15$.

The full order-11 computation closes the remaining two units.  Switch all
edges incident with a root positive and let $G$ be the graph of negative
residual edges on the other ten vertices.  For $S\subseteq V(G)$ and
$e=|E(G)|$,

\[
 Q(S)=55-2e-2|S|(11-|S|)+4|\delta_G(S)|.                 \tag{25.1}
\]

If $M\le15$, the all-positive state gives $20\le e\le35$.  Global negation
followed by root normalization complements $G$, so one may take
$20\le e\le22$.  The singleton constraints then give maximum degree at most
5 for $e=20,21$ and at most 6 for $e=22$.

Nauty generates 12,005,168 unlabeled graphs on ten vertices.  The full stream
is hashed and exactly filtered to the classes allowed by the preceding
reductions; integer evaluation of (25.1) leaves no graph with $M\le15$.  A reduced
stream generated by

```text
geng -q -D6 10 20:22
```

contains 2,153,606 records and reproduces the zero-survivor result with the
same evaluator.  Deterministic samples at edge counts 20, 21, and 22 are
recomputed by a separate adjacency formula.  The full and reduced
newline-delimited graph6 stream hashes are respectively

```text
5650c7c979fdffd8c0f99a2f2ee8775938ec2a3dd69aa65be1207936824fc5b3
b62da4d7ebfaab4ccd801fd509f2fc85f6c2b815c8c1d2e969de7aa6a82c322d
```

The graph6 witness `ICRbczQMo` has maximum 17 over all 1,024 projective
spins.  Hence, with the explicit dependency boundary that nauty emits one
representative of every unlabeled graph,

\[
 \boxed{F(11)=17}.                                      \tag{25.2}
\]

Monotonicity gives $F(12)\ge17$, and every order-12 energy is even, so
$F(12)\ge18$.  The residual graph6 witness `JWUuDOR\K{?` has maximum 18 over
all 2,048 projective spins.  Therefore, at the same computer-assisted
standard,

\[
 \boxed{F(12)=18}.                                      \tag{25.3}
\]

Both witnesses are recomputed by the cut formula and by a separate direct
quadratic-form loop in `verification/research_order11_certify.py`.  No MILP,
SAT status, floating point, or timeout is used in (25.2)--(25.3).

### Exact orders 13 and 14

The next lower certificate is again a complete nauty enumeration, but now the
residual graph has eleven vertices.  The $1{,}018{,}997{,}864$ records are
split into eight `res/mod` shards.  The verifier hashes and counts the exact
bytes relayed from each generator process to a separately compiled C scanner,
checks both exit statuses, and compares each shard with committed record
counts, byte counts, SHA-256 digests, and survivor lists.

The scanner enumerates all $2^{11}$ projective spins of each root-normalized
order-12 signing.  It rejects a record as soon as it witnesses
$|Q|\ge 19$.  Every order-12 energy is even, so this is an exact
test for whether $M\le 18$.  Exactly two rooted residual records survive:

```text
JCpVdXyxpz?
JCpdUg{[dM?
```

A separate direct evaluator gives $M=18$ for each and enumerates all
$2^{11}$ projective incident columns, obtaining extension minimum $24$
for both.  All other order-12 predecessors have $M\ge20$, and one-vertex
extension never lowers $M$.  The exact Bellman identity therefore gives
$F(13)\ge20$.  Every one-point deletion of the order-14 Paley conference
matrix has maximum $20$, so


\[
 \boxed{F(13)=20}.                                      \tag{25.3a}
\]

Heredity gives $F(14)\ge20$, while all order-14 energies are odd.  Direct
evaluation gives $M(C_{14})=21$, hence

\[
 \boxed{F(14)=21}.                                      \tag{25.3b}
\]

The complete replay is in `verification/research_order13_certify.py`; its
fixed-size threshold kernel is `verification/order12_threshold_scan.c`.
This is a computer-assisted theorem with the explicit trust boundary that
nauty's eight shards cover the unlabeled residual graphs.  The proof does not
use the scanner's above-threshold energy histogram, a seeded incumbent, or a
solver status.

A separate quarantined reconstruction also matched the complete order-14
$M=21$ and $M=23$ layers and gives the computational bracket
$F(15)\in\{25,27\}$.  Indeed, if an order-15 signing had maximum at most
23, then each order-14 deletion would have odd maximum 21 or 23; the
reconstructed extension minima 27 and 29 rule this out.  Extending $C_{14}$
gives $F(15)\le27$, and order-15 parity leaves only 25 or 27.  That layer
driver has not yet been integrated into the repository, so the main
exact-value table stops at order 14.  The exact value $F(15)$ remains open;
the unclassified branch consists of order-13 predecessors with
$(M,\delta_{\rm w})=(24,0)$.
A later external update reports 470 such classes on one validation shard and
an ongoing 16.3-million-candidate sweep, but the corresponding v2 artifacts
are not in the ingested bank and no conclusion from that run is used here.

### Theorem 25 (exact negative-replica supermultiplicativity)

For uniform disorder $A$ of order $n$, define

\[
 f_{n,t}(A)=\frac{Z_A(t)}{2^{n+1}(\cosh t)^{m_n}},\qquad
 \mathcal G_n(q,t)=\log\mathbb E_A f_{n,t}(A)^{-q}.
\]

For every $n,k\ge1$, $q>0$, and real $t$,

\[
 \boxed{\mathcal G_{n+k}(q,t)\ge
 \mathcal G_n(q,t)+\mathcal G_k(q,t)}.                 \tag{25.4}
\]

#### Proof

Fix internal signings $A,B$ and let $C$ be a uniform $n\times k$ cross
signing.  Write

\[
 U_A^\sigma=\sum_xe^{t\sigma Q_A(x)},\qquad
 U_B^\sigma=\sum_ye^{t\sigma Q_B(y)}.
\]

Independent averaging of the cross edges gives

\[
 \mathbb E_C Z_{A,B,C}(t)=(\cosh t)^{nk}P,
 \quad P=U_A^+U_B^++U_A^-U_B^-.
\]

Replacing $B$ by $-B$ gives
$P'=U_A^+U_B^-+U_A^-U_B^+$, and
$P+P'=Z_AZ_B$.  Since $z\mapsto z^{-q}$ is convex, Jensen first over $C$
and then over the two orientations gives

\[
\begin{aligned}
&\frac12\left(\mathbb E_CZ_{A,B,C}^{-q}
 +\mathbb E_CZ_{A,-B,C}^{-q}\right)\\
&\quad\ge \frac{(\cosh t)^{-qnk}}2(P^{-q}+{P'}^{-q})\\
&\quad\ge 2^q(\cosh t)^{-qnk}(Z_AZ_B)^{-q}.
\end{aligned}
\]

Uniform disorder is invariant under $B\mapsto-B$.  Averaging $A,B$ proves
the unnormalized inequality.  Finally,
$m_{n+k}=m_n+m_k+nk$; the factors $2^q$ and
$(\cosh t)^{-qnk}$ cancel exactly against the definition of $f$.  This proves
(25.4).  At $t=0$ equality holds. \(\square\)

### Reverse-hypercontractive parameter mismatch

Put $u=\tanh t$ and regard $f_{n,u}=f_{n,\operatorname{arctanh}u}$ as a
function of the disorder signs.  Also abbreviate
$\mathcal G_n(q,u)=\mathcal G_n(q,\operatorname{arctanh}u)$.  The Fourier expansion of $f_{n,u}$ is the signed
even-Eulerian polynomial already derived above.  Therefore

\[
 f_{n,u}=T_{u/v}f_{n,v}\qquad(0<u<v<1),                \tag{25.5}
\]

where $T_\rho$ is the Boolean noise operator.  The sharp reverse
hypercontractive inequality for positive functions says, for $Q'>Q>0$,

\[
 \|T_\rho f\|_{-Q'}\ge\|f\|_{-Q}
 \quad\hbox{if}\quad
 \rho\le\sqrt{\frac{1+Q}{1+Q'}}.
\]

Since $\|f\|_{-Q}=\exp[-\mathcal G_n(Q,u)/Q]$ in this parametrization, it follows that

\[
 \frac{\mathcal G_n(Q',u')}{Q'}
 \le\frac{\mathcal G_n(Q,u)}Q                       \tag{25.6}
\]

on the invariant curve $(1+Q')u'^2=(1+Q)u^2$.  This is exactly the
parameter curve generated by block composition, but (25.6) points in the
wrong direction for closing the lower recursion (25.4).

Fix $\beta,\theta>0$ and, for all sufficiently large $r$, put

\[
 t_r=\frac\beta{\sqrt r},\qquad
 q_r=\frac\theta{\tanh^2t_r}-1
 =\frac\theta{\beta^2}r+\frac{2\theta}{3}-1+O(r^{-1}).
\]

A natural proposed repair was the following power-saving transport estimate:
for proportional child sizes,

\[
0\le \frac{q_N}{q_n}\mathcal G_n(q_n,t_n)
 -\mathcal G_n(q_N,t_N)\le C N^{2-\delta}.             \tag{PT}
\]

Only the first inequality follows from (25.6). If the upper bound held, then
for $a_n=\mathcal G_n(q_n,t_n)/q_n$, (25.4) would give

\[
 a_{n+k}\ge a_n+a_k-O((n+k)^{1-\delta})
\]

on balanced splits.  Also $0\le a_n\le m_n\log\cosh t_n=O_\beta(n)$.
and the following elementary blocking lemma would give convergence.

**Balanced almost-superadditivity lemma.**  If $b_n\ge0$, $b_n=O(n)$, and

\[
 b_{j+k}\ge b_j+b_k-K(j+k)^{1-\delta}
\]

whenever $j,k$ are between one third and two thirds of their sum, then
$b_n/n$ converges.

To prove it, fix a large $m$ and write $N=qm+r$, $0\le r<m$.  Partition $N$
into $q-1$ leaves of weight $m$ and one of weight $m+r$; when $r=0$, use $q$
leaves of weight $m$.  Recursively split the number of leaves as evenly as
possible, placing the single heavier leaf on the smaller side.  Every child
then has between one third and two thirds of its parent's weight.  Charge the
error $KW^{1-\delta}$ at a node of weight $W$ to its leaves proportionally to
their weights.  Along a root-to-leaf path the weights grow by at least $3/2$,
so the total charge per unit leaf weight is

\[
 K\sum_{h\ge0}((3/2)^hm)^{-\delta}=O(m^{-\delta}).
\]

The accumulated error is $O(Nm^{-\delta})$.  Dropping the nonnegative
contribution of the exceptional leaf gives

\[
 \liminf_{N\to\infty}\frac{b_N}{N}
 \ge \frac{b_m}{m}-O(m^{-\delta}).
\]

Letting $m$ tend to infinity along a limsup subsequence proves convergence.
Applying the lemma to $a_n$ would prove convergence of $a_n/n$ under (PT).
For

\[
 \mathcal S_{n,\beta,\theta}
 =-\frac1{q_n\beta n}\log\mathbb E_AZ_A(t_n)^{-q_n},
\]

the exact soft-minimum sandwich has limiting error at most
$\beta\log2/(2\theta)$. Thus the proposed estimate would have settled the
original problem for unbounded parameter sets. Theorem 27 below proves that
it is false.

For comparison, exact class-weighted enumeration gives the balanced defects

\[
 \Delta_n=\frac{q_{2n}}{q_n}\mathcal G_n(q_n,t_n)
 -\mathcal G_n(q_{2n},t_{2n}):
\]

\[
\begin{array}{c|rrrrrr}
n&4&5&6&7&8&9\\ \hline
\Delta_n/n^2,\ (\beta,\theta)=(1,1)
&.00492&.00878&.01166&.01387&.01560&.01700\\
\Delta_n/n^2,\ (\beta,\theta)=(2,4)
&.14850&.27713&.32759&.28609&.30762&.31597\\
\Delta_n/n^2,\ (\beta,\theta)=(4,8)
&.57748&1.04741&1.17936&.99417&1.10086&1.11063
\end{array}
\]

These floating-point evaluations use exact energy spectra and exact labeled
class multiplicities through order 9. They are finite regression evidence
only; the disproof is the analytic liminf estimate in Theorem 27.

### Theorem 26 (known square-order Paley ceiling attainment)

Let $m$ be an odd prime power.  A Paley conference matrix $C$ of order
$m^2+1$ has a Boolean eigenvector of eigenvalue $m$, and hence

\[
 \boxed{M(C)=\frac{m(m^2+1)}2}.                         \tag{26.1}
\]

This is the known regular Paley conference construction, included here to
correctly classify the external report rather than claim novelty.  Let
$K=\mathbb F_{m^2}$, $H=\mathbb F_m$, and let $\chi$ be the quadratic
character of $K$.  Index $C$ by $\{\infty\}\cup K$, with
$C_{\infty,a}=1$ and $C_{a,b}=\chi(a-b)$.  Choose signs
$\sigma:K/H\to\{\pm1\}$ with $(m+1)/2$ positive cosets, so
$\sum_c\sigma(c)=1$, and put

\[
 x_\infty=1,\qquad x_a=\sigma(a+H).
\]

Every element of $H^*$ is a square in $K$.  Thus the character sum over the
zero additive coset is $m-1$.  Multiplication by $H^*$ is transitive on the
nonzero cosets of the one-dimensional quotient $K/H$ and preserves $\chi$;
all their character sums are equal.  Since the total character sum on $K$ is
zero, each nonzero coset sum is $-1$.  Consequently, for $a$ in coset $c$,

\[
 (Cx)_a=1+(m-1)\sigma(c)-\sum_{d\ne c}\sigma(d)=m\sigma(c),
 \qquad (Cx)_\infty=m.
\]

Hence $Cx=mx$.  The conference identity $C^2=m^2I$ supplies the matching
spectral upper bound and proves (26.1).  Exact self-contained checks for
$m=3,5,7$ are in `verification/verify_paley_subfield.py`.

This theorem concerns those particular full Paley matrices, not the minimax
$F(m^2+1)$.  Already $m=3$ gives $M(C)=15$ at order 10, while $F(10)=13$.
It does not improve the $1/2$ upper bound and gives no lower bound on $F$.

### Theorem 27 (the proposed transport lemma is false)

Let

\[
 \Delta_n(\beta,\theta)
 =\frac{q_{2n}}{q_n}\mathcal G_n(q_n,t_n)
 -\mathcal G_n(q_{2n},t_{2n}),
 \qquad
 t_r=\frac\beta{\sqrt r},\quad
 q_r=\frac\theta{\tanh^2t_r}-1.
\]

If

\[
 U=\limsup_{n\to\infty}\frac{F(n)}{n^{3/2}},
\]

then

\[
 \boxed{
 \liminf_{n\to\infty}\frac{\Delta_n(\beta,\theta)}{n^2}
 \ge
 \frac{2\theta}{\beta^2}
 \left[\frac{\beta^2}{8}
 -\beta\left(1-\frac1{\sqrt2}\right)U-\log2\right]-\log2.}
 \tag{27.1}
\]

In particular, using $U\le1/2$ and $(\beta,\theta)=(4,8)$,

\[
 \boxed{
 \liminf_n\frac{\Delta_n(4,8)}{n^2}
 \ge\sqrt2-2\log2
 =0.027919201253\ldots>0.}                             \tag{27.2}
\]

Hence no bound $\Delta_n=O(n^{2-\delta})$ can hold uniformly in the proposed
parameter range.

#### Proof

Write $m=\binom n2$ and $f_n^*(t)=\min_A f_{n,t}(A)$. Vertex switching and
global coefficient negation preserve $f_{n,t}$. For $n\ge3$ these operations
form a free group of order $2^n$: global negation cannot be a switching,
because switching signs multiply to $+1$ around every triangle. Thus the
disorder cube splits into $2^{m-n}$ equal-size orbits, and

\[
 -\log f_n^*(t)-\frac{(m-n)\log2}{q}
 \le\frac{\mathcal G_n(q,t)}q
 \le-\log f_n^*(t).                                   \tag{27.3}
\]

For an $F(n)$-optimal signing,
$Z_A(t)\le2^{n+1}e^{tF(n)}$, so

\[
 \log f_n^*(t)\le tF(n)-m\log\cosh t.                 \tag{27.4}
\]

Conversely every signing has a maximizing projective pair $x,-x$. Choosing
the favorable augmented sign for those two states gives
$Z_A(t')\ge2e^{t'F(n)}$, hence

\[
 \log f_n^*(t')\ge t'F(n)-m\log\cosh t'-n\log2.       \tag{27.5}
\]

Apply (27.3) at $(q_n,t_n)$ and $(q_{2n},t_{2n})$, then use
(27.4)--(27.5). This gives the finite inequality

\[
\begin{aligned}
 \Delta_n\ge q_{2n}\bigg[&
 m(\log\cosh t_n-\log\cosh t_{2n})
 -(t_n-t_{2n})F(n)-n\log2\\
 &-\frac{(m-n)\log2}{q_n}\bigg].                     \tag{27.6}
\end{aligned}
\]

Now $q_n/n\to\theta/\beta^2$, $q_{2n}/n\to2\theta/\beta^2$, and

\[
 \frac mn(\log\cosh t_n-\log\cosh t_{2n})\to\frac{\beta^2}{8}.
\]

Taking the liminf in (27.6) proves (27.1), and substitution proves (27.2).
$\square$

More generally, the bracket before the final entropy penalty is positive for

\[
 \beta>3.012373175204\ldots,
\]

and then sufficiently large $\theta$ gives a leading gap. At $\beta=4$ the
threshold is $\theta>7.690245425859\ldots$.

### Theorem 28 (quantitative entropy production)

Let $m\ge1$, let $\nu$ be uniform on $\{\pm1\}^m$, and let $\theta>0$.
On any interval of nonnegative times for which
$q_s=\theta e^{2s}-1>0$, let

\[
 L=\frac12\sum_{i=1}^m(\tau_i-I),\qquad
 f_s=T_{e^{-s}}f_0>0,\qquad q_s+1=\theta e^{2s},
\]

and define

\[
 \mu_s=\frac{f_s^{-q_s}}{\mathbb E_\nu f_s^{-q_s}}\nu,
 \quad D_s=D(\mu_s\Vert\nu),
 \quad H_s=\frac1{q_s}\log\mathbb E_\nu f_s^{-q_s}.
\]

Then

\[
 \boxed{\frac{dH_s}{ds}=-R_s},\qquad
 R_s=\mathbb E_{\mu_s}\frac{Lf_s}{f_s}
 -\frac{2(q_s+1)}{q_s^2}D_s,                          \tag{28.1}
\]

and

\[
 \boxed{R_s\ge\frac{2D_s^2}{3q_sm}.}                 \tag{28.2}
\]

Equality in (28.2) holds only when $f_s$ is constant.

#### Proof

Differentiating $H_s$, using $\partial_sf_s=Lf_s$ and
$q_s'=2(q_s+1)$, gives (28.1) because

\[
 D_s=-q_s\mathbb E_{\mu_s}\log f_s
 -\log\mathbb E_\nu f_s^{-q_s}.
\]

For one coordinate, condition on all other coordinates and call the two
function values $a,b$. Put

\[
 z=\frac q2\left|\log\frac ba\right|,
 \quad d(z)=z\tanh z-\log\cosh z.
\]

Here $d(z)$ is the conditional KL divergence. The conditional generator
contribution is

\[
 e_q(z)=\frac12\left(
 \frac{\cosh((1+2/q)z)}{\cosh z}-1\right).
\]

For
$r_q(z)=e_q(z)-2(q+1)d(z)/q^2$, put $c=2/q$. Then

\[
 q e_q(z)=\frac{\cosh(cz)-1+\tanh z\sinh(cz)}c
 \ge z\tanh z+c\,d(z),
\]

where $d(z)\le z^2/2$ was used. Therefore

\[
 qr_q(z)\ge2\log\cosh z-z\tanh z=:h(z).
\]

Finally $h(z)\ge2d(z)^2/3$. Indeed the derivative of
$h-2d^2/3$ is nonnegative after using $d(z)\le z^2/2$ and
$\sinh z\cosh z\ge z+2z^3/3$. Thus

\[
 r_q(z)\ge\frac{2d(z)^2}{3q}.                         \tag{28.3}
\]

Average (28.3) over all coordinate conditionals. If
$\bar d_i=\mathbb E d_i$, entropy tensorization and Jensen--Cauchy give

\[
 D_s\le\sum_i\bar d_i,\qquad
 \sum_i\mathbb E d_i^2\ge\sum_i\bar d_i^2
 \ge\frac{(\sum_i\bar d_i)^2}{m}\ge\frac{D_s^2}{m}.
\]

More exactly,

\[
 R_s=\sum_i\mathbb E r_q(z_i)
 +\frac{2(q+1)}{q^2}\left(\sum_i\bar d_i-D_s\right),
\]

whose second term is the nonnegative dual total correlation. This proves
(28.2). The scalar inequality is strict for $z>0$, so equality forces every
conditional $z_i=0$. Hence $f_s$ is invariant under every coordinate flip
and is constant. $\square$

For transport endpoints $s_n<s_N$ on the invariant curve,

\[
 \Delta_{n,N}=q_N\int_{s_n}^{s_N}R_s\,ds.             \tag{28.4}
\]

For proportional endpoints with fixed $\beta,\theta$, one has
$q_s\asymp n$, $m\asymp n^2$, and a fixed-length integration interval.
Thus a power-saving defect would force
$\int D_s^2ds=O(n^{4-\delta})$. An escort law with $D_s\ge cn^2$ on a
fixed-length subinterval necessarily gives a leading $\Omega(n^2)$ defect.

### Theorem 29 (exact conditional-alignment chain)

For a split $N=n+k$ with $n,k\ge3$, let $R_{n,k}$ be the binary rectangular
row-and-column switching code and set

\[
 C'_{n,k}=D_n\oplus D_k\oplus R_{n,k}.
\]

Then $D_N\subset C'_{n,k}$ and

\[
 |C'_{n,k}/D_N|=2^{N-1}.                              \tag{29.1}
\]

Let $\mathcal G_{n,k}^{\square}(q,t)$ be the negative-moment logarithm for
the annealed-normalized rectangular partition function. On a coarse coset
$h$ let $r_h$ be the conditional fine density relative to the uniform fiber
law, and put

\[
 K_q(h)=\mathbb E_{\rm fiber}r_h^{-q}.
\]

Then

\[
 \boxed{
 \mathcal G_N(q,t)=\mathcal G_n(q,t)+\mathcal G_k(q,t)
 +\mathcal G_{n,k}^{\square}(q,t)+\mathcal A_{n,k}(q,t),}
 \tag{29.2}
\]

where

\[
 \mathcal A_{n,k}(q,t)=\log\mathbb E_{\widehat U_q}K_q(h)\ge0. \tag{29.3}
\]

Here $\widehat U_q$ is the product coarse quotient law tilted by its negative
moment.

#### Proof

The dimensions are

\[
 \dim D_n=n,\quad \dim D_k=k,\quad
 \dim R_{n,k}=N-1,
\]

so $\dim C'_{n,k}=2N-1$, while $\dim D_N=N$. The restriction of every global
cut and the global all-one word lies in $C'_{n,k}$, proving (29.1). The coarse
quotient factors because the three coordinate sets are disjoint and the raw
edge-noise law is a product.

For a coarse representative $(A,B,C)$, write

\[
 Y_{\alpha,\beta,\tau}=
 \begin{pmatrix}
 D_\alpha AD_\alpha&C\\
 C^{\mathsf T}&\tau D_\beta BD_\beta
 \end{pmatrix},
\]

where $\alpha,\beta$ are projective sign vectors and $\tau\in\{\pm1\}$.
These are the $2^{N-1}$ fine states. With

\[
 Z_C^\square(t)=\sum_{x,y}e^{t x^{\mathsf T}Cy},
\]

direct changes of spin variables give

\[
 \mathbb E_{\alpha,\beta,\tau}Z_Y(t)
 =\frac{Z_A(t)Z_B(t)Z_C^\square(t)}{2^{N+1}}.
\]

Consequently

\[
 r_{A,B,C}(\alpha,\beta,\tau)
 =\frac{2^{N+1}Z_Y(t)}{Z_A(t)Z_B(t)Z_C^\square(t)},
 \qquad \mathbb E_{\rm fiber}r=1.                    \tag{29.4}
\]

Grouping the uniform fine negative moment first by coarse fibers gives
(29.2)--(29.3); convexity of $z\mapsto z^{-q}$ gives $K_q(h)\ge1$.
$\square$

Define

\[
 \mathcal A_{n,k}^{\infty}(t)
 =\lim_{q\to\infty}\frac{\mathcal A_{n,k}(q,t)}q.
\]

For

\[
 h_n(t)=\min_A\log Z_A(t),\qquad
 h_{n,k}^{\square}(t)=\min_C\log Z_C^\square(t),
\]

finite-space Laplace limits in (29.2) give

\[
 \mathcal A_{n,k}^{\infty}(t)
 =h_n(t)+h_k(t)+h_{n,k}^{\square}(t)
 -(N+1)\log2-h_N(t).
\]

Taking $t\to\infty$ sequentially gives

\[
 \boxed{
 \lim_{t\to\infty}\frac{\mathcal A_{n,k}^{\infty}(t)}t
 =F(n)+F(k)+B_\square(n,k)-F(N).}                     \tag{29.5}
\]

This does not exchange $q,t,$ and $N$ limits. Since
$B_\square(n,n)\ge n\,\mathbb E|\varepsilon_1+\cdots+\varepsilon_n|$,
(29.5) and the audited asymptotic bounds give

\[
 \liminf_n\frac{2F(n)+B_\square(n,n)-F(2n)}{n^{3/2}}
 \ge\frac2\pi+\sqrt{\frac2\pi}-\sqrt2
 =0.0202907708\ldots>0.                               \tag{29.6}
\]

Thus alignment is a leading ground-state quantity.

The rectangular marginal alone cannot repair transport. Put
$\theta=\lambda\beta^2$ and

\[
 T_n=\frac{\mathcal G_n(q_n,t_n)}{q_n}
 -\frac{\mathcal G_n(q_{2n},t_{2n})}{q_{2n}}.
\]

The rectangular lower bound
$B_\square(n,n)\ge n\mathbb E|\varepsilon_1+\cdots+\varepsilon_n|$ and the
finite-state max-divergence estimate give

\[
 \limsup_n\frac1n\frac{\mathcal G_{n,n}^{\square}(q_{2n},t_{2n})}{q_{2n}}
 \le\frac{\beta^2}{4}-\frac\beta{\sqrt\pi}+2\log2.
\]

Inserting the audited $U\le1/2$ in (27.1) and dividing by
$q_{2n}\sim2\lambda n$ gives the corresponding lower bound on $T_n/n$.
Subtraction yields

\[
\begin{aligned}
\limsup_n\frac1n\left[
 \frac{\mathcal G_{n,n}^{\square}(q_{2n},t_{2n})}{q_{2n}}-2T_n\right]
\le{}&4\log2+\frac{\log2}{\lambda}\\
&-\beta\left[\frac1{\sqrt\pi}
 -\left(1-\frac1{\sqrt2}\right)\right].              \tag{29.7}
\end{aligned}
\]

The coefficient of $\beta$ is $0.2712963647\ldots>0$. First take
$n\to\infty$ for each fixed $(\beta,\lambda)$. Along a subsequent
zero-temperature sequence with $\beta\to\infty$ and
$\lambda\beta\to\infty$, the right side tends to $-\infty$. Thus any
cancellation-based settling argument must obtain a leading positive
contribution from conditional alignment.

### Theorem 30 (mixed four-cycle alignment Hamiltonian)

Put $P=D_\alpha AD_\alpha$, $Q=D_\beta BD_\beta$. The coefficient of $u^4$
in the conditional log-density $\log r$ is

\[
\boxed{
\mathcal H_4=\frac12\left[
 \operatorname{tr}(P^2CC^{\mathsf T})
 +\operatorname{tr}(Q^2C^{\mathsf T}C)
 +\tau\operatorname{tr}(PCQC^{\mathsf T})
 -nk(N-2)\right].}                                   \tag{30.1}
\]

Equivalently,

\[
\begin{aligned}
\mathcal H_4={}&
 \sum_{i<j}(A^2)_{ij}(CC^{\mathsf T})_{ij}\alpha_i\alpha_j
 +\sum_{a<b}(B^2)_{ab}(C^{\mathsf T}C)_{ab}\beta_a\beta_b\\
&+\tau\sum_{i<j,a<b}a_{ij}b_{ab}
 (C_{ia}C_{jb}+C_{ib}C_{ja})
 \alpha_i\alpha_j\beta_a\beta_b.                   \tag{30.2}
\end{aligned}
\]

#### Proof

The graph and rectangular dual codes have minimum weight four. Their signed
weight-four sums are

\[
 c_4(S)=\frac18[\operatorname{tr}S^4-v(v-1)(2v-3)],
\]

and

\[
 r_4(C)=\frac14[\operatorname{tr}((CC^{\mathsf T})^2)-nk(N-1)].
\]

Direct block multiplication gives

\[
\begin{aligned}
\operatorname{tr}Y^4={}&\operatorname{tr}P^4+\operatorname{tr}Q^4
 +2\operatorname{tr}((CC^{\mathsf T})^2)\\
&+4\operatorname{tr}(P^2CC^{\mathsf T})
 +4\operatorname{tr}(Q^2C^{\mathsf T}C)
 +4\tau\operatorname{tr}(PCQC^{\mathsf T}).
\end{aligned}
\]

Subtracting the three marginal weight-four coefficients proves (30.1).
Expanding the switched entries and cancelling the diagonal baseline proves
(30.2). $\square$

The Walsh characters in (30.2) are orthogonal. Hence

\[
\begin{aligned}
\mathbb E\mathcal H_4^2={}&
 \sum_{i<j}[(A^2)_{ij}(CC^{\mathsf T})_{ij}]^2
 +\sum_{a<b}[(B^2)_{ab}(C^{\mathsf T}C)_{ab}]^2\\
&+\sum_{i<j,a<b}(C_{ia}C_{jb}+C_{ib}C_{ja})^2.        \tag{30.3}
\end{aligned}
\]

For each row pair the last sum counts four times the pairs of equal signs in
the row product, so

\[
 \mathbb E\mathcal H_4^2
 \ge4\binom n2\left\lfloor\frac{(k-1)^2}{4}\right\rfloor. \tag{30.4}
\]

Thus the first alignment Hamiltonian is nonconstant for $n,k\ge3$. Higher
Eulerian terms cannot be discarded uniformly at
$u\asymp N^{-1/2},q\asymp N$.

There is also an exact scalar-state collision at split $2+4$. Use right-edge
order $(01,02,03,12,13,23)$ and

\[
 A=(1),\qquad B=(-1,-1,-1,-1,-1,1),
\]

with

\[
 C_1=\begin{pmatrix}-1&-1&-1&-1\\-1&-1&1&1\end{pmatrix},
 \qquad
 C_2=\begin{pmatrix}-1&-1&-1&-1\\-1&1&-1&1\end{pmatrix}.
\]

The two triples have identical complete local absolute-energy histograms,
and therefore identical three scalar partition curves. At
$u=3/5$, so $t=\log2$, exact rational evaluation gives

\[
 K_2(C_1)=
 \frac{196585091273040100817}{133610891512185651200},
 \qquad
 K_2(C_2)=\frac{6723290161}{5922841600}.               \tag{30.5}
\]

Their $\mathcal H_4$ laws are respectively

\[
 \{-12:4,-4:4,0:16,4:4,12:4\},\qquad
 \{-4:8,0:16,4:8\}.
\]

Here $\dim D_2=1$, so the displayed 32 relative-gauge tuples cover the true
16-point fiber twice. Every density and $\mathcal H_4$ value is repeated
twice, and the normalized fiber moments in (30.5) are unchanged. Thus all
three local scalar free-energy curves fail to determine alignment.

Fix $\beta,\theta>0$, put
$t_j=\beta/\sqrt j$, $q_j=\theta/\tanh^2(t_j)-1$, and define

\[
 a_j=\frac{\mathcal G_j(q_j,t_j)}{q_j},\qquad
 T_{j\leftarrow N}
 =a_j-\frac{\mathcal G_j(q_N,t_N)}{q_N}.
\]

The tempting alignment-transport bound

\[
 \frac{\mathcal G_{n,k}^{\square}+\mathcal A_{n,k}}{q_N}
 \ge T_{n\leftarrow N}+T_{k\leftarrow N}-O(N^{1-\delta})
\]

is sufficient but not yet a reduction. The exact identity is

\[
 \frac{\mathcal G_{n,k}^{\square}+\mathcal A_{n,k}}{q_N}
 -T_{n\leftarrow N}-T_{k\leftarrow N}
 =a_N-a_n-a_k.
\]

Thus the displayed inequality is the desired almost-superadditivity in new
notation. The non-tautological target is a variational or power-saving
composition theorem for the relative-switching Hamiltonian hierarchy itself.

### Theorem 31 (exact Paley Fourier-leakage identity)

Let $q\equiv1\pmod4$ be an odd prime power and let $C$ be the Paley
conference matrix of order $q+1$. Fix a nontrivial additive character of
$\mathbb F_q$. For a Boolean function $f:\mathbb F_q\to\{\pm1\}$, let
$S(f)=\sum_a f(a)$ and use the unitary additive Fourier transform. Partition
the nonzero frequencies into $H_+$ and $H_-$ according to the sign of the
quadratic Gauss-sum eigenvalue, and set

\[
 E_\pm(f)=\sum_{r\in H_\pm}|\widehat f(r)|^2,
 \qquad W(f)=\min(E_+(f),E_-(f)).
\]

Then

\[
 \boxed{
 M(C)=\frac{(q+1)\sqrt q}{2}
 -\sqrt q\min_f\left[
 W(f)+\frac{(|S(f)|-\sqrt q)^2}{2q}
 \right].}                                             \tag{31.1}
\]

Indeed, the Paley core has Fourier multipliers $0,+\sqrt q,-\sqrt q$.
Parseval gives

\[
 f^{\mathsf T}Kf=\sqrt q(E_+-E_-),
 \qquad E_++E_-=q-\frac{S(f)^2}{q}.
\]

Writing a Boolean vector as $(s,f)$ and optimizing its infinity sign gives

\[
 \max_{s=\pm1}|Q_C(s,f)|
 =|S(f)|+\frac{\sqrt q}{2}
 \left(q-\frac{S(f)^2}{q}-2W(f)\right).
\]

Subtracting from the spectral ceiling and completing the square proves
(31.1). Consequently,

\[
 \frac{2M(C)}{(q+1)\sqrt q}\longrightarrow1            \tag{31.2}
\]

along a sequence of Paley orders if and only if there are Boolean $f_q$ with

\[
 |S(f_q)|=o(q),\qquad W(f_q)=o(q).                      \tag{31.3}
\]

Thus Paley spectral alignment is exactly a one-sided additive-Fourier
concentration problem. The infinity coordinate is optimized exactly; the
substantive issue is minority-half Fourier leakage.

Fresh source-built exhaustive scans give

\[
\begin{array}{c|rrrr}
q&5&13&17&29\\ \hline
M(C_{q+1})&5&21&33&75\\
2M/((q+1)\sqrt q)&\sqrt5/3&3/\sqrt{13}&
11/(3\sqrt{17})&5/\sqrt{29}.
\end{array}
\]

The four selected prime-field ratios increase, but this is finite evidence
only. It proves neither monotonicity nor a power-law deficit. The exact checks,
including $2^{29}$ projective states at order 30, are in
`verification/research_paley_alignment.c`.

Even (31.2) would not settle the minimax problem: it is an upper-bound fact
about one signing. A sufficient independent rigidity statement is

\[
 F(q+1)\ge(1-o(1))M(C_{q+1})                            \tag{31.4}
\]

along primes $q\equiv1\pmod4$. Equations (31.2)--(31.4), the prime number
theorem in that progression, monotonicity of $F$, and the existing limsup
bound would force the original limit to be $1/2$. Exact optimality at orders
6 and 14 is far too little evidence for (31.4), and order 10 is a finite
counterexample to unconditional conference optimality.

## 9. Current best continuation targets

The finite-temperature route now has an exact scalar no-go theorem and an
exact nonscalar state. Theorem 27 closes PT with a leading transport gap, and
Theorem 28 quantifies non-saturation through escort-law entropy production.
Theorem 29 shows what scalar transport discarded: a conditional law on
relative block switchings whose zero-temperature slope is exactly the joint
optimizer-composition gain. Theorem 30 identifies its first mixed Hamiltonian
and proves that even all three local scalar partition curves do not determine
it.

Theorem 31 adds a separate, value-specific route. Its tractable half is the
Paley leakage target (31.3), a Boolean Fourier/character-sum problem. Its hard
half is the dense-order minimax rigidity estimate (31.4). The increasing
alignment ratios through order 30 do not provide a lower bound for $F$ and do
not displace the relative-switching free-energy route as the main
value-agnostic program.

The next finite-temperature problem is therefore not another scalar reverse
hypercontractive estimate. It is a variational limit or power-saving
composition theorem for the relative-switching alignment free energy,
uniform on the mean-field diagonal. A claim phrased only as the displayed
alignment-transport inequality is tautological; it must be derived from a
smaller state such as a controlled hierarchy of connected mixed Eulerian
clusters.

The ground-state alternative is now the Pareto-profile refinement following
Theorem 14. It seeks Theorem 13's power-saving composition bound after
allowing $o(n^{3/2})$ internal slack. Exact minimizers cannot be required:
$K_{2,8}(1,10)=15$ even though $F(10)=13$, whereas
$K_{2,8}(1,12)=13$. Current block constructions still fail because internal
block energy and the cross field can cancel for the same spin assignment,
while bounding the two pieces separately incurs a leading $N^{3/2}$ cost.
Theorems 11 and 12 give hereditary lower information but no such
optimizer-compatible upper construction.

Theorem 18 gives an exact cross-block state, but Proposition 23 proves that
the weighted-entropy target (19.2) is impossible by a leading constant gap.
The live multivertex target is instead the switching-orbit noncoverage
statement (23.3), or another estimate that exploits the intersections among
bad rank-one constraints.  Proposition 20 controls the one-vertex state on a
density-one set, while Theorem 24 shows that even uniform scalar
$O(\sqrt n)$ control would not suffice.  The one-vertex target is
stabilization of (24.2), for example through the Bellman--Cesaro law (24.4).

The ground-state form of the same alignment problem is an optimizer-specific
stability theorem for the Bellman state. It would supplement

\[
 R_{n+k}(\beta)
 \le R_n(\beta\sqrt{n/(n+k)})
 +R_k(\beta\sqrt{k/(n+k)})
\]

with a reverse or same-temperature comparison of total cost $o(n+k)$.

The recursion alone is insufficient because it compares different
temperatures, and the abstract countermodel proves that scalar cavity bounds
cannot repair it. Covariance and limiting pair overlap are also too coarse.
A successful state must retain enough relative eigenvector and mixed-cycle
geometry to determine the optimized cavity reward. On the coding side, the
corresponding target is a conditional strong data-processing or reverse-noise
theorem for fine cosets inside the graph--rectangular coarse quotient, still
effective when $u=\tanh(\beta/\sqrt n)$.

## 10. Relative-profile composition and Paley leakage closure

### Theorem 32 (microcanonical relative-switching composition)

Let

\[
 \mathbb P_r=\{\pm1\}^r/\{\pm\mathbf1\},\qquad
 G_{n,k}=\mathbb P_n\times\mathbb P_k\times\{\pm1\},
 \qquad |G_{n,k}|=2^{N-1},\quad N=n+k.
\]

For fixed graph signings $A,B$ and rectangular signing $C$, take the local
augmented graph state spaces

\[
 \Omega_A=\{\pm1\}\times\mathbb P_n,\qquad
 e_A(\sigma,[x])=\sigma Q_A(x),
\]

and analogously for $B$. Project signed rank-one matrices by $R\sim-R$ and
put $e_C([R])=|\langle C,R\rangle|$. The map

\[
 \pi((\sigma,[z]),(\eta,[w]),[xy^{\mathsf T}])
 =([zx],[wy],\sigma\eta)
\]

is onto $G_{n,k}$ and every fiber has $2^{N-1}$ elements.

Write

\[
 L=M(A)+M(B)+B_C,\qquad B_C=\max_{x,y}|x^{\mathsf T}Cy|,
\]

and let $d_A,d_B,d_C$ be the corresponding local energy deficits. For
$g=(\alpha,\beta,\tau)$ define

\[
 Y_g=\begin{pmatrix}
 D_\alpha AD_\alpha&C\\
 C^{\mathsf T}&\tau D_\beta BD_\beta
 \end{pmatrix}.
\]

Then exactly

\[
 L-M(Y_g)=\min_{\pi(a,b,r)=g}[d_A(a)+d_B(b)+d_C(r)].       \tag{32.1}
\]

Let $\Lambda(A,B,C)$ be the $2^{N-1}$-st smallest product-triple deficit,
with multiplicity. Pigeonholing one minimizing triple in each fiber gives

\[
 F(N)\le M(A)+M(B)+B_C-\Lambda(A,B,C).                    \tag{32.2}
\]

Equivalently, if fewer than $2^{N-1}$ triples have deficit below $s$, some
relative gauge has maximum at most $L-s$.

This is not the failed weighted Hamming union bound. It is a deterministic
quantile theorem for the actual three local near-maximal profiles. If

\[
 D_A(t)=\sum_a e^{-td_A(a)}
\]

and similarly for $B,C$, exponential counting gives

\[
 \Lambda\ge\sup_{t>0}\left[
 \frac{(N-1)\log2-\log D_A(t)-\log D_B(t)-\log D_C(t)}t
 \right]_+.                                                \tag{32.3}
\]

The bound is optimal after discarding the additive labels: distribute the
$2^{N-1}$ smallest deficits one per abstract equal-sized fiber. The actual
group geometry can do better. In the exact $2+4$ scalar collision, both
unlabeled profiles have $\Lambda=0$, whereas the true relative-gauge gains
are 4 and 2.

At finite temperature, $\pi$ pushes the product of the two graph Gibbs laws
and the rectangular $2\cosh(t\langle C,R\rangle)$ law to the conditional
alignment law. Its Fourier coefficients factor into three local correlation
tensors:

\[
 \widehat r(I,J,\epsilon)
 =\mathbb E_{\mu_A}[\sigma^\epsilon z_I]
  \mathbb E_{\mu_B}[\eta^\epsilon w_J]
 \mathbb E_{\mu_C}[x_Iy_J].                               \tag{32.4}
\]

This supplies a monotone projected hierarchy between the scalar quantile and
the complete alignment state. The mixed four-cycle Hamiltonian is its first
nonzero layer.

The exact scalar settling target is now a large-deviation statement. For a
composable near-optimal family, fewer than $2^{N-1}$ triples must lie below
the deficit needed to reduce

\[
 M(A)+M(B)+B_C
\]

to

\[
 [F(n)^{2/3}+F(k)^{2/3}]^{3/2}+O(N^{3/2-\delta}).
\]

Then (32.2) yields power-saving near-subadditivity. A bounded number of exact
maximizers is not enough; one needs the whole leading-scale near-maximal
profile.

### Theorem 33 (multiplicatively dense prime Paley alignment)

For a prime $p\equiv1\pmod4$, let $\ell(p)$ be its least positive quadratic
nonresidue. Put $m=(p-1)/2$ and take the balanced interval function

\[
 f_p(a)=-1\quad(0\le a<m),\qquad f_p(a)=1\quad(m\le a<p).
\]

Then $S(f_p)=1$ and, for $r\ne0$,

\[
 |\widehat f_p(r)|^2
 =\frac4p\frac{\sin^2(\pi rm/p)}{\sin^2(\pi r/p)}.
\]

The inequality $\sin(\pi\|r\|_p/p)\ge2\|r\|_p/p$ gives

\[
 W(f_p)\le E_-(f_p)
 \le2p\sum_{d=\ell(p)}^\infty d^{-2}
 \le\frac{2p}{\ell(p)-1}.                                 \tag{33.1}
\]

For

\[
 L_j=8\prod_{\substack{\lambda\le j\\\lambda\ {\rm odd\ prime}}}\lambda,
\]

quadratic reciprocity shows that $p\equiv1\pmod{L_j}$ implies
$\ell(p)>j$. The prime number theorem in each fixed progression lets one
concatenate increasingly deep progression tails so that consecutive prime
ratios tend to one while the levels $j$ tend to infinity. Therefore there
is a multiplicatively dense prime sequence $p_i$ with

\[
 |S(f_{p_i})|=1=o(p_i),\qquad W(f_{p_i})=o(p_i).
\]

By Theorem 31,

\[
 \frac{2M(C_{p_i+1})}{(p_i+1)\sqrt{p_i}}\longrightarrow1.  \tag{33.2}
\]

This closes the Fourier-leakage half of the value-specific Paley route. It
does not prove a lower bound for $F$. The remaining statement is the dense
sequence rigidity estimate

\[
 F(p_i+1)\ge(1-o(1))M(C_{p_i+1}).                          \tag{33.3}
\]

If (33.3) held, monotonicity and multiplicative density would force the
original limit to equal $1/2$.

There are two exact limitations. First, nonconstant Boolean functions at
prime order never have $W=0$; a cyclotomic norm argument gives

\[
 W(f)\ge\frac{4(p-1)^2}{p^3}.
\]

Second, one interval does not work uniformly. Along primes
$p\equiv5\pmod{12}$, the residue pair $\pm1$ and nonresidue pair $\pm3$
force

\[
 \liminf W(f_p)/p\ge8/(9\pi^2).
\]

The nested arithmetic progressions are therefore essential to this
construction.

The deterministic checks are
`verification/verify_relative_profile_composition.py` and
`verification/verify_paley_least_nonresidue.py`. They verify the finite
identities, exact witnesses, order statistics, and corruption controls. The
prime-sequence conclusion rests on the analytic proof, not finite sampling.

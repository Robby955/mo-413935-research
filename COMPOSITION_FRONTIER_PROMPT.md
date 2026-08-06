# Self-contained frontier prompt: relative-gauge composition

You are being asked to make a serious, auditable attempt at the remaining
composition problem behind MathOverflow Question 413935:

https://mathoverflow.net/questions/413935/min-max-of-a-quadratic-form-of-plus-minus-ones

The public research repository is:

https://github.com/Robby955/mo-413935-research

Read the entire repository before beginning. Record the exact commit you read.
The broad exploration frontier was frozen at tag
`research-frontier-2026-08-05`. On the current default branch, prioritize:

- `ACTIVE_RESEARCH.md`
- `paper/composition_framework.tex`
- `paper/finite_results.tex`
- `paper/second_attempt.tex`
- `AUDIT.md`
- `RESEARCH_LEDGER.md`
- `RESEARCH_CONTINUATION.md`
- every script and receipt under `verification/`

Do not rely only on this prompt. Independently check every theorem you use.

## Problem

For a symmetric zero-diagonal sign matrix A, define

```text
Q_A(x) = sum_{1 <= i < j <= n} a_ij x_i x_j,
M(A)   = max_{x in {+-1}^n} |Q_A(x)|,
F(n)   = min_A M(A).
```

Determine whether

```text
F(n) / n^(3/2)
```

has a limit. Existence or nonexistence is enough; its value is not required.
The question remains open.

## Audited baseline

The repository proves

```text
1/pi <= liminf F(n)/n^(3/2)
     <= limsup F(n)/n^(3/2) <= 1/2,

F(n) <= F(n+1) <= F(n) + n.
```

It also proves the exact augmented cut-code identity

```text
F(n) = choose(n,2) - 2 rho(D_n),
```

and the exact values

```text
F(2),...,F(14) = 1,3,4,4,5,9,10,12,13,17,18,20,21.
```

The order-11 and order-13 lower bounds are computer-assisted and have an
explicit nauty completeness boundary. The pinned oriented-Hadamard
construction proves only

```text
F(16) <= 30.
```

Do not promote this to equality. Do not import unbanked claims about F(15) or
higher orders.

## The convergence criterion

Put

```text
H(n) = F(n)^(2/3).
```

The repository proves that the following estimate would settle the problem:

```text
H(n+k) <= H(n) + H(k) + O((n+k)^(1-delta))
```

for some fixed `delta > 0`, uniformly over all sizes. Equivalently at the
known scale, it is enough to prove

```text
F(n+k)
 <= (F(n)^(2/3) + F(k)^(2/3))^(3/2)
    + O((n+k)^(3/2-delta)).
```

An unspecified little-o error is not automatically safe under repeated
composition. Require a dyadically summable modulus, such as a power saving.

## The exact surviving state

Fix internal signings A and B of orders n and k, a rectangular signing C,
and N=n+k. Let

```text
P_r      = {+-1}^r / global sign,
G_{n,k}  = P_n x P_k x {+-1},
|G|      = 2^(N-1).
```

Use augmented graph states

```text
Omega_A = {+-1} x P_n,
e_A(sigma,[x]) = sigma Q_A(x),
```

and analogously for B. Let the projective rank-one cross states be

```text
Rbar_{n,k} = {[xy^T]},
e_C([R])   = |<C,R>|.
```

There is a balanced homomorphism

```text
pi((sigma,[z]),(eta,[w]),[xy^T])
  = ([zx],[wy],sigma eta)
```

onto `G_{n,k}`. Its product domain has `2^(2N-2)` elements and every fiber
has exactly `2^(N-1)` elements.

Define

```text
M_A = M(A),
M_B = M(B),
B_C = max_{x,y} |x^T C y|,
L   = M_A + M_B + B_C,

d_A = M_A - e_A,
d_B = M_B - e_B,
d_C = B_C - e_C.
```

For `g=(alpha,beta,tau)`, define the full signing

```text
Y_g = [[D_alpha A D_alpha, C],
       [C^T, tau D_beta B D_beta]].
```

The exact max-plus identity is

```text
L - M(Y_g)
 = min_{pi(a,b,r)=g} [d_A(a)+d_B(b)+d_C(r)].
```

This is the central state. Do not replace it by separate maxima or by an
unlabeled product distribution unless you explicitly account for the loss.

## Labeled occupancy target

For a deficit threshold s, define

```text
b_s(g) = number of product triples in fiber g
         with d_A + d_B + d_C < s.
```

Then

```text
b_s(g)=0  if and only if  M(Y_g) <= L-s.
```

At the desired composition ceiling

```text
T = (F(n)^(2/3)+F(k)^(2/3))^(3/2)
    + O(N^(3/2-delta)),
s = L-T,
```

the sharp remaining task is to prove that a suitable choice of A, B, and C
has at least one empty fiber. A stronger and more stable target is to prove
that exponentially many fibers are empty or have the required gain.

Characters of `G_{n,k}` are indexed by even subsets I and J and a bit epsilon.
For each local deficit shell define its signed character sums

```text
A_d(I,epsilon), B_d(J,epsilon), C_d(I,J).
```

The exact Fourier factorization is

```text
hat b_s(I,J,epsilon)
 = 1/2^(N-1) * sum_{d_A+d_B+d_C<s}
   A_dA(I,epsilon) B_dB(J,epsilon) C_dC(I,J).
```

This formula, not another scalar relaxation, is the preferred starting point.

## Finite-temperature version

The annealed-normalized negative moment has the exact chain

```text
Gamma_N
 = Gamma_n + Gamma_k + Gamma_square_{n,k} + Alignment_{n,k},
Alignment_{n,k} >= 0.
```

The alignment term is the reverse moment of the conditional density on the
same relative-switching fibers. Its sequential zero-temperature slope is

```text
F(n) + F(k) + B_square(n,k) - F(n+k).
```

It is provably leading order. The first nonconstant Eulerian term in its log
density is

```text
H4(alpha,beta,tau) = 1/2 [
  tr(P^2 C C^T)
  + tr(Q^2 C^T C)
  + tau tr(P C Q C^T)
  - nk(N-2)
],
```

where `P=D_alpha A D_alpha` and `Q=D_beta B D_beta`.

You may pursue a variational limit for this full state, but merely rewriting
the desired almost-superadditivity in terms of `Alignment` is tautological.
You must derive a quantitative inequality from a smaller controlled object.

## What has been closed

Do not spend the run rediscovering any of the following unless your mechanism
materially changes the state retained:

1. The relaxed coupling cube has unique optimizer zero and a leading-order
   integrality gap.
2. The adversarial elliptope relaxation has normalized limit 1/2 but does not
   control the discrete minimax.
3. Ordinary graphons and empirical spectral laws erase the n^(3/2) scale.
4. Ordinary block composition and cross-block triangle bounds lose the full
   leading order.
5. Scalar negative-replica parameter transport is false with a power-saving
   defect; reverse hypercontractivity has the wrong-sided comparison.
6. Adding only the rectangular marginal cannot compensate for that defect.
7. The weighted Hamming-ball union-bound certificate has a universal
   leading-order floor above 1/2.
8. Complete local scalar partition curves do not determine alignment.
9. The sharp unlabeled microcanonical order statistic can miss a good labeled
   gauge. In the balanced C14 calibration, subthreshold triples outnumber the
   fiber budget while the exact labeled occupancy still has an empty fiber.
10. Generic low-degree moment, polynomial vacancy, and small Fourier-support
    certificates cannot reliably detect an isolated empty fiber.
11. Isolated Hadamard constructions without an iterable framed invariant do
    not imply convergence.
12. Paley spectral alignment is only an upper construction. The missing
    minimax rigidity lower bound on a dense sequence is equivalent in
    difficulty to proving the full limit is 1/2.

## Active tasks, ranked

### Task 1: good-fiber abundance

Try to prove a theorem of the following form for proportional splits:

```text
There exist near-optimal A,B and a cross seed C such that at the target
deficit, at least 2^(cN) relative gauges g satisfy b_s(g)=0,
```

or a quantitatively equivalent low-occupancy statement with a summable energy
loss. The abundance requirement is not mandatory for logical sufficiency, but
it avoids the precision wall created by a single exceptional gauge.

### Task 2: exact character-sum sign

Use the labeled shell factorization to force `min_g b_s(g)=0` directly. Look
for algebraic signs, parity, congruences, or Fourier support constraints that
cannot occur in an everywhere-positive integer occupancy.

### Task 3: closed relative-switching variational state

Find a finite or compact state with:

- an exact or controlled composition rule;
- enough information to recover optimizer-compatible cancellation;
- a power-saving error under proportional composition;
- a proof that higher Eulerian layers are either included or uniformly
  controlled.

The framed order-16 Hadamard theorem shows that unframed switching class is
not sufficient: orientation relative to the cross frame changes the maximum
by eight units. A candidate state must pass this finite test.

## Standards

For every proposed theorem:

1. State it with exact quantifiers and normalization.
2. Prove every analytic and combinatorial step.
3. Track the absolute value and projective multiplicities exactly.
4. Check all powers of two in the relative-gauge fibers.
5. Distinguish a sufficient restatement from a genuinely easier lemma.
6. Show that every error is `O(N^(3/2-delta))` in F scale, or otherwise
   dyadically summable in H scale.
7. Test small cases with independent implementations and corruption controls.
8. Do not treat timeout or heuristic search failure as infeasibility.
9. Do not claim the original question is solved unless the cross-order theorem
   is uniform and complete.

## Requested output

Produce:

1. A short audit of the exact relative-gauge statements you use.
2. At least three serious derivations within the active labeled-state route,
   not three unrelated branches of mathematics.
3. A complete proof of any new theorem, including edge cases.
4. Counterexamples to false intermediate claims.
5. Reproducible scripts and expected output for every finite check.
6. A patch against the exact repository commit you read, without deleting the
   existing failed-route ledger.
7. A final status sentence using one of:
   - convergence proved;
   - nonconvergence proved;
   - substantial new cross-order theorem;
   - active target remains open, with the sharpest new obstruction stated.

Do not broaden the project with category theory, Langlands theory, another
generic relaxation, or another exact value unless you can exhibit a concrete
map from that work to the labeled fiber inequality above.

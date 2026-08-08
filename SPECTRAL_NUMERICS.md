# Why the upper bound 1/2 is hard to improve

This file records the measurements behind the current active route. It is
evidence and orientation, not proof: nothing here changes the audited bounds

```text
1/pi <= liminf F(n)/n^(3/2) <= limsup F(n)/n^(3/2) <= 1/2.
```

## The ceiling, and the one quantity that controls it

For any signing `A`, Cauchy-Schwarz gives `M(A) <= n*||A||/2`, and
`||A|| >= sqrt(n-1)` always, with equality exactly for conference matrices.
So the whole upper bound is

```text
F(n) <= n*sqrt(n-1)/2,
```

and no spectral argument can do better. Write the alignment fraction

```text
m(A) = M(A) / (n*||A||/2),
```

the share of the ceiling that sign vectors actually reach. Two bounds hold for
every `A`: `m <= 1` trivially, and `m >= 2/pi = 0.6366` because Nesterov
rounding of the semidefinite relaxation is within `2/pi`, and for a conference
matrix the relaxation value equals the ceiling exactly. The gap between `1/pi`
and `1/2` is precisely the gap between those two facts.

## A first-moment cap for flat spectra

Suppose `A` has a flat spectrum, all eigenvalues `+-sqrt(n)`, with `P+` the
projection onto the top eigenspace of dimension `n/2`. Then

```text
x'Ax = sqrt(n) * (2|P+x|^2 - n),
```

so the question is only how well a sign vector aligns with `P+`. For a Haar
random subspace and fixed `x` with `|x|^2 = n`, the quantity `|P+x|^2/n` is
Beta(n/4, n/4), whose large-deviation rate at `t` is
`(1/4)*log(1/(4t(1-t)))` per coordinate. Setting `t = (1+m)/2` gives
`4t(1-t) = 1-m^2`, and balancing against the `2^n` sign vectors,

```text
(1/4)*log(1/(1-m^2)) = log 2   =>   1-m^2 = 1/16   =>   m <= sqrt(15)/4 = 0.9682.
```

So a flat-spectrum matrix with no arithmetic structure cannot reach more than
about `0.968` of its ceiling, which would give a constant near `0.484`, below
`1/2`.

This is a first-moment count for Haar-random flat matrices. It is not a
theorem about `+-1` matrices, and the next section shows it is false for
Paley.

## Measurements

Haar-random flat-spectrum matrices, diagonal zeroed to match the constraint,
maximum found by multi-start local search, so these are lower bounds on `m`
(`verification/generic_flat.py`):

```text
n      64      128     256     512     1024
m      0.9554  0.9423  0.9344  0.9277  0.9250
```

They sit under the cap and drift down slowly.

Paley conference matrices, exact for `n <= 26` and search lower bounds beyond
(`verification/paley_scan.py`):

```text
n           6      14      18      30      38      42      54      74     102     114
M/ceiling   0.745  0.832   0.889   0.929   0.943   0.959   0.962   0.940  0.978   0.949
M/n^(3/2)   0.340  0.401   0.432   0.456   0.465   0.474   0.476   0.467  0.486   0.472
```

At `n = 102` the Paley alignment is at least `0.978`, above the `0.9682` cap
that Haar-random flat matrices obey. Since these are lower bounds on `M`, that
comparison is safe in the direction it is used.

## The obstruction, stated plainly

A `+-1` matrix with a perfectly flat spectrum **is** a conference matrix, and
the constructible conference matrices are arithmetic. Paley exceeds the
generic cap rather than obeying it, in agreement with the least-nonresidue
alignment theorem recorded in `RESEARCH_CONTINUATION.md`. So the quasirandom
heuristic that predicts a constant below `1/2` fails on exactly the family one
can write down, and the Paley route cannot improve the upper bound.

Improving `limsup` therefore needs a `+-1` family that provably behaves like a
generic flat matrix, or a different mechanism altogether.

## The optimizers are not the flattest matrices

Writing `rho = ||A||/sqrt(n-1)`, so `rho = 1` means conference
(`verification/opt_spectrum.py`):

```text
optimum        M     M/n^(3/2)   rho      m
n = 10         13    0.411       1.311    0.661
n = 13         20    0.427       1.041    0.853
n = 14         21    0.401       1.000    0.832
```

Order 10 is the informative one. Its optimum inflates the spectral norm by 31
percent and buys an alignment of `0.661`, close to the Nesterov floor
`0.6366`, and beats the conference matrix 13 to 15 by refusing to be flat. At
orders 13 and 14 the optimum is flat again, and at 14 the conference matrix is
optimal. The extremal object changes character with `n`, which is one reason
the small table does not extrapolate.

## Scope

Every large-`n` maximum above is a local-search value and therefore a lower
bound on the true maximum, not an exact one. The `sqrt(15)/4` cap is a
first-moment estimate for a random model, stated here to explain what the
measurements are being compared against. Zeroing the diagonal of a Haar-flat
matrix moves the spectral norm by `O(1)` against a scale of `sqrt(n)`, which
is why it is ignored. None of this bounds `F(n)` for any finite `n`.

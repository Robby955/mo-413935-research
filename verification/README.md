# Verification guide

All pass/fail decisions in the new identity checker and exhaustive search use
integer or rational arithmetic. Floating point is used only for displayed
normalized values and for the original Gaussian sampling check.

Use Python 3.10 or newer. The original two scripts need only the standard
library:

```bash
python3 verification/verify_attempt.py
python3 verification/check_conference_examples.py
```

The new exact identity checker also needs only the standard library:

```bash
python3 verification/verify_new_results.py
python3 verification/verify_continuation.py
python3 verification/verify_coding_continuation.py
python3 verification/verify_amplification_obstructions.py
python3 verification/verify_cavity_hereditary.py
```

It checks the linear augmented code and its even-Eulerian dual through
\(n=6\), the coset MacWilliams identity in exact rational arithmetic, Walsh
and fourth-cumulant identities, the SDP square-covariance construction,
higher-moment counterexamples, and the absolute-value counterexample. It has
deterministic corruption controls.

Expected output:

    code_dual_orthogonality_checks=33864
    parseval_cumulant_signings_checked=1258
    order_7_trace_minimizers_checked=3024
    order_7_sixth_moment_classes=3
    macwilliams_cosets_checked=56
    covariance_sdp_matrices_checked=9
    absolute_value_orders_checked=13
    deterministic_seed=413935
    corruption_controls=PASSED

The continuation checker independently enumerates the new cavity and block
partition identities, tests the Gaussian covariance and determinant algebra,
verifies the finite negative-replica sandwich, checks the Hadamard cavity
counterexample, and exercises the abstract oscillating countermodel. Expected
output:

    cavity_extensions_checked=3294
    partition_inequalities_checked=24
    covariance_edges_checked=50280
    determinant_products_checked=4904
    negative_replica_checks=16
    hadamard_cavity_checks=2
    abstract_countermodel_checks=2199
    false_conference_determinant_bound_detected=TRUE
    deterministic_seed=413935
    corruption_controls=PASSED

Its Gaussian calculations use floating point only to check finite algebraic
instances. The entropy theorem itself rests on the analytic data-processing
inequality and Gibbs variational principle; this script is not a proof of
those results.

The amplification-obstruction checker uses exact integer arithmetic for the
order-five Kronecker counterexample and its orbit reduction. Expected output:

    order_5_energy_set=-4,0,4
    diagonal_pairs_checked=1024
    balanced_completion_orbits=10,10
    certificate_energies=NN:88,NA:-80,AA:88,uniform:100
    tensor_iterations_checked=2
    corruption_controls=PASSED

The coding-continuation checker uses exact rational arithmetic for the
coset-noise law, multiaffine box check, transfer sectors, and
deletion--contraction formulas. Expected output:

    deletion_contraction_checks=9
    transfer_sector_checks=192
    noise_cosets_checked=144
    multiaffine_grid_points=729
    k3_extension_values=28/27,80/81
    strict_noise_monotonicity_detected=TRUE
    omitted_deletion_term_detected=TRUE
    corruption_controls=PASSED

The hereditary-cavity checker exhausts every signing, every vertex subset,
every fixed spin on the retained block, and every paired completion through
order 5. It verifies both exact cavity inequalities, the closed formula for
$\mu_k$, all induced consequences for the independently recomputed small
$F(n)$ values, and the exact failure of subadditivity for
$H(n)=F(n)^{2/3}$ at $2+2=4$. Expected output:

    orders_exhausted=1..5
    signings_checked=1099
    partitions_checked=33866
    fixed_x_checks=254253
    paired_y_checks=1065508
    mu_formula_checks=6
    f_consequence_checks=15
    small_F_values=F(1)=0,F(2)=1,F(3)=3,F(4)=4,F(5)=4
    mu_values_k_0_to_5=0,1,1,3/2,3/2,15/8
    h_exact_subadditivity_counterexample=2+2->4
    corruption_controls=walsh_edge_double_count,quadratic_flip_parity,mu_strengthening,h_exact_subadditivity
    cavity_hereditary_verification=PASSED

## Exact cross-block composition

`research_cross_block_composition.py` requires nauty `geng`, NetworkX, and
the `z3-solver` package:

```bash
python3 verification/research_cross_block_composition.py
```

The script reconstructs every optimum switching-permutation class through
order 9. For each pair of nontrivial optimum blocks with total order at most
10, it solves the cross-edge feasibility problem by a deterministic Z3
cutting-plane loop and recomputes every returned witness on the full spin
cube. Every pair with at most 16 cross edges is independently checked by
enumerating all cross signings.

One-vertex extensions use direct enumeration only. The critical `2+8`
non-heredity theorem also uses direct enumeration for every optimum order-8
root representative. A second exhaustive pass through all 1,044 residual
graphs at order 8 proves that internal maximum 12 is the least slack that can
participate in an order-10 optimum.

For each optimal extension class, the script also computes the projective
covering radius `rho_ext` of the exact extremizers and the energy-weighted
radius `rho_weighted` of the full energy landscape. It independently asserts

```text
E(B) = M(B) + n - 2*rho_weighted
extremizer-only value = M(B) + n - 2*rho_ext
```

against direct enumeration of incident-edge signings.

It then enumerates every root-normalized residual graph, not just optimal
classes, through order 8 and computes the complete Pareto frontier of
`(M(B), delta_weighted)`. Every radius-derived extension value is checked by a
second direct enumeration of all incident sign vectors.

The expected critical records are:

```text
pair 2+8: class distribution {15: 2}, while F(10)=13
pair 3+7: class distribution {13: 6, 15: 6}
order-7 extensions: {10: 4, 12: 2} across 6 classes
order-9 extensions: {13: 4, 15: 11} across 15 classes
order-10 extensions: {17: 1, 19: 1} across 2 classes; value-17 witness HCRbczQ, mask 440
order-7 covering profiles: {(3,3,10): 4, (3,2,12): 2}
order-9 covering profiles: {(4,4,13): 4, (3,3,15): 11}
order-10 covering profiles: {(3,3,17): 1, (2,2,19): 1}
partition-function collision: G?qmaw and GCpbaw have histogram {0:60,4:111,8:60,12:25}, extensions 13 and 15
full Bellman frontier at order 6: {(5,2), (7,1), (9,0)}; all give F(7)=9
full Bellman frontier at order 8: {(10,1), (12,0)}; both give F(9)=12
extension obstruction layers: order_7_abs_7=2, order_9_abs_12=11
minimum order-8 block maximum for an order-10 optimum: 12
order-8 value-12 root representatives composing to 13: 68 of 104
slack witness: graph6 F?reg, cross coefficient mask 52010
direct_crosschecks=23
status=PASSED
```

For the slack witness, the graph6 word records negative residual edges after
the root edges have been switched positive. Cross-mask bit `8*i+j` equals one
when the edge from left vertex `i` to right vertex `j` is negative.

Certificate boundary: the optimizer catalogue trusts nauty `geng`, asserted
graph counts, and committed stream hashes. Switching-class counts trust
NetworkX isomorphism. Z3 infeasibility is trusted only for pair rows with more
than 16 cross entries. The principal `2+8` obstruction, every one-vertex
extension, and the internal-slack witness do not rely on Z3. Completeness of
the stored order-10 optimum-class list is independently reproduced by

```bash
python3 verification/research_exact_small_n.py \
  --min-n 10 \
  --max-n 10 \
  --classify-switching-optima \
  --strict-stream-digests \
  --labeled-crosscheck-max-n 0
```

## Exhaustive values through n = 10

`research_exact_small_n.py` requires Brendan McKay's nauty `geng` executable.
On Homebrew systems it is normally installed by `brew install nauty`; the
script also accepts `--geng /absolute/path/to/geng`. NetworkX is optional and
is used only for an independent graph6 decoder and switching-class grouping.

The normal audit command is:

```bash
python3 verification/research_exact_small_n.py \
  --max-n 10 \
  --networkx-crosscheck \
  --classify-switching-optima \
  --strict-stream-digests
```

The expected exact values are

```text
n:     2  3  4  5  6  7   8   9  10
F(n):  1  3  4  4  5  9  10  12  13
```

The search switches all root edges positive and enumerates one unlabeled
negative-edge graph on the remaining vertices. It asserts the A000088 graph
counts, hashes every full graph6 stream, checks the reduced energy formula
directly, enumerates every unreduced labeled signing through \(n=6\), tests
switching invariance, and reconstructs an optimizer spectrum. The full
\(n=10\) run checks 274,668 residual graphs. Exact output is JSON Lines so it
can be archived or diffed.

Certificate boundary: completeness trusts nauty `geng` plus the asserted
counts and committed stream hashes. This is an auditable exhaustive
computation, not a formal proof-assistant certificate.

## Independent Z3 route

Install the `z3-solver` Python package and run:

```bash
python3 verification/research_z3_certify.py
```

The default run independently proves UNSAT/SAT pairs for \(n=7,8,9\):

```text
n=7: UNSAT at 7,  SAT at 9
n=8: UNSAT at 8,  SAT at 10
n=9: UNSAT at 10, SAT at 12
```

Energy parity makes each pair exact. The returned SAT witness is recomputed
over every spin vector by plain Python and hashed. The solver seed is fixed at
413935. Z3 is a trusted solver in this workflow; the script does not export a
DRAT or LFSC proof. The \(n=10\), bound-11 UNSAT query timed out in the audited
120-second trial, so the nauty enumeration—not Z3—is the lower certificate for
\(F(10)=13\).

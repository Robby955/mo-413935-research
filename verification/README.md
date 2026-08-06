# Verification guide

## Focused manuscript builds

The focused finite-results and composition-framework manuscripts are built
independently of the broad research archive:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -jobname=mo-413935-finite-results \
  -output-directory=paper paper/finite_results.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -jobname=mo-413935-composition-framework \
  -output-directory=paper paper/composition_framework.tex
```

Expected artifacts:

```text
paper/mo-413935-finite-results.pdf
paper/mo-413935-composition-framework.pdf
```

The broad archival manuscript remains `paper/second_attempt.tex` and is not
generated from either focused source.

All pass/fail decisions in the new identity checker and exhaustive search use
integer or rational arithmetic. Floating point is used only for displayed
normalized values and for the original Gaussian sampling check.

Use Python 3.10 or newer. The original two scripts need only the standard
library:

```bash
python3 verification/verify_attempt.py
python3 verification/check_conference_examples.py
```

The conference checker also verifies that all fourteen order-13 principal
submatrices of the order-14 Paley matrix have maximum 20. The complete
matching lower certificate is described below.

`research_paley_alignment.c` independently constructs the prime-field Paley
conference matrices at orders 6, 14, 18, and 30, verifies their conference
identities, and scans every projective Boolean state by Gray code. Build it
from source; no banked executable is used:

```bash
cc -std=c11 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror \
  verification/research_paley_alignment.c -lm \
  -o /tmp/research_paley_alignment
/tmp/research_paley_alignment
```

Expected output:

```text
q=5 order=6 conference=PASS M=5 projective_maximizers=12 ratio=sqrt(5)/3 ratio_squared=5/9 decimal=0.745355992499930
q=13 order=14 conference=PASS M=21 projective_maximizers=156 ratio=3/sqrt(13) ratio_squared=9/13 decimal=0.832050294337844
q=17 order=18 conference=PASS M=33 projective_maximizers=204 ratio=11/(3*sqrt(17)) ratio_squared=121/153 decimal=0.889297291799888
q=29 order=30 conference=PASS M=75 projective_maximizers=812 ratio=5/sqrt(29) ratio_squared=25/29 decimal=0.928476690885259
corruption_control=symmetric_edge_flip_detected
all_checks=PASS
```

All maxima, maximizer counts, and squared ratios are asserted in exact integer
arithmetic. Only the displayed decimal ratios use floating point. The four
selected ratios do not constitute a monotonicity or asymptotic theorem.

The least-nonresidue checker verifies the arithmetic inputs and exact finite
witnesses for the interval construction behind the new asymptotic Paley
alignment theorem:

```bash
python3 verification/verify_paley_least_nonresidue.py
```

Expected output:

```text
p=73 level=3 least_nonresidue=5 S=1 Q_core=264 Q_full=265 leakage_bound=36.500000000000 witness_ratio=0.838267675864
p=241 level=5 least_nonresidue=7 S=1 Q_core=1704 Q_full=1705 leakage_bound=80.333333333333 witness_ratio=0.907675246292
p=2521 level=7 least_nonresidue=11 S=1 Q_core=60624 Q_full=60625 leakage_bound=504.200000000000 witness_ratio=0.957525265511
p=9241 level=11 least_nonresidue=13 S=1 Q_core=430728 Q_full=430729 leakage_bound=1540.166666666667 witness_ratio=0.969635838899
fourier_case=p73 E_residue=67.392013544710 E_nonresidue=5.594287825153 parseval=PASS
corruption_controls=composite_modulus,interval_endpoint
paley_least_nonresidue_verification=PASSED
```

The prime-sequence conclusion itself uses the prime number theorem in each
fixed arithmetic progression and is proved in the note. The script does not
pretend to verify that analytic theorem by finite sampling.

The relative-profile checker exhausts all `2 + 3` block triples, verifies the
balanced relative-gauge map and exact max-plus identity, and checks the
microcanonical order statistic and exponential-profile bound:

```bash
python3 verification/verify_relative_profile_composition.py
```

Expected output:

```text
block_cases_checked=1024
balanced_fibers_checked=16384
product_triples_checked=262144
exponential_profile_checks=672
strict_geometry_gains=768
scalar_collision=lambda:0,0 true_gains:4,2
corruption_control=relative_orientation_omission_detected
relative_profile_composition_verification=PASSED
```

The `2 + 4` collision check independently confirms that identical unlabeled
scalar profiles can have zero order-statistic guarantee but different true
relative-gauge gains. All decisions in this checker use integer or rational
arithmetic.

## Balanced relative-profile calibration at orders 12--14

`research_relative_profile_calibration.py` computes the complete augmented
graph-deficit and absolute rectangular-deficit histograms for every balanced
vertex split of four banked optimal witnesses: the two rooted encodings of the
same order-12 optimum, a principal order-13 submatrix of `C14`, and `C14`
itself. It then convolves the three local histograms, evaluates the exact
microcanonical order statistic `Lambda`, and compares the resulting scalar
bound with the true gain:

```bash
python3 verification/research_relative_profile_calibration.py
```

Expected output:

```text
case=JCpVdXyxpz? order=12 partitions=462 profile_types=12 lambda=0..10 scalar_bound=24..26 geometry_gap=6..8 child_optimal_splits=7 best_bound=24 target=14.142135623731 excess=9.857864376269 target_lattice_energy=16 target_triples=213312..269122 fiber_count=2048 sha256=2c29479f48ff9644bfb860640ba35450af837502371b17f59153653b708fbc01
gain_pairs=0/6:6,2/8:1,2/10:70,4/12:135,8/16:120,10/18:130
case=JCpdUg{[dM? order=12 partitions=462 profile_types=12 lambda=0..10 scalar_bound=24..26 geometry_gap=6..8 child_optimal_splits=7 best_bound=24 target=14.142135623731 excess=9.857864376269 target_lattice_energy=16 target_triples=213312..269122 fiber_count=2048 sha256=5ff5291913cd0cf4db0e3cb2648934e99607fc59f184ff734abb9a765c31d444
gain_pairs=0/6:6,2/8:1,2/10:70,4/12:135,8/16:120,10/18:130
case=C14-minus-infinity order=13 partitions=1716 profile_types=9 lambda=2..10 scalar_bound=28..30 geometry_gap=8..10 child_optimal_splits=52 best_bound=28 target=19.524318098857 excess=8.475681901143 target_lattice_energy=20 target_triples=436864..495632 fiber_count=4096 sha256=b9ce312b08962293639003979cac33a0c6a852dd40209a00bd23653eb437d8aa
gain_pairs=2/10:52,6/14:156,6/16:338,8/16:156,8/18:546,10/20:468
case=C14 order=14 partitions=1716 profile_types=3 lambda=6..10 scalar_bound=31..33 geometry_gap=10..12 child_optimal_splits=624 best_bound=31 target=25.455844122716 excess=5.544155877284 target_lattice_energy=27 target_triples=231581..305465 fiber_count=8192 sha256=e96d8070d53fcdd2a47ce95eddb865d666c74dd02557df53ae8ec993d3fbd4b7
gain_pairs=6/16:364,6/18:260,10/22:1092
alternate_swapped_profile=lambda:8,target_count:596440
corruption_controls=balanced_double_count,cross_absolute_value
relative_profile_calibration=PASSED
```

The real target for a split of orders `n + k` is

```text
(F(n)^(2/3) + F(k)^(2/3))^(3/2).
```

Every local product-triple energy has the parity of `binom(n+k,2)`. Exact
rational inequalities place the three displayed targets strictly in the
intervals `(14,16)`, `(18,20)`, and `(25,27)`. Thus the first admissible
energy strictly above the target is exactly 16, 20, and 27, respectively.
`target_triples` counts product triples at or above this lattice energy, so it
is exactly the number whose energy is strictly above the real target; floating
point is used only to display the target and excess. The verifier asserts the
complete `(independent ceiling, Lambda, target count)` distribution embedded
in its fixture, not only the printed range.

The full product profile has `2^(2N-2)` entries and the relative-gauge map has
`2^(N-1)` fibers. In all of these finite cases the target-triple count exceeds
the fiber count, so the unlabeled scalar pigeonhole theorem alone does not
reach the zero-error target. The observed scalar-to-true gaps of `6..8`,
`8..10`, and `10..12` are finite evidence only. In particular, they do not
rule out an asymptotically acceptable `O(N)` or other power-saving defect.

The exact true-gain calculation uses the separately certified values
`F(12)=18`, `F(13)=20`, and `F(14)=21`: the identity relative gauge is one of
the scanned full signings, while every other gauge is another signing of the
same order. This script scans all balanced splits of the listed witnesses; it
does not certify that the displayed order-13 and order-14 witnesses exhaust
all optimal switching classes. Here `profile_types` means distinct
eight-scalar summary records, not distinct complete histogram triples.

The corruption controls detect omission of the unordered-balanced-split guard
and replacement of the absolute rectangular energy by a signed one. The script
also evaluates a separate one-sided profile that replaces the augmented graph
states by projective `|Q|` states while expanding the cross factor to signed
full-spin pairs. Its balanced fiber maximum dominates rather than equals the
full signing maximum. On the standard `C14` split this alternate theorem gives
`Lambda=8` and target count `596440`; the exact max-plus theorem gives
`Lambda=10`, true gain `22`, and target count `304908`. Exact partition
counts, subset-labelled stream digests,
gain-pair distributions,
target-count distributions, and the equality of the two order-12 aggregate
summary distributions provide additional fail-closed checks.

## One-sided swapped-profile injection

`verify_swapped_profile_injection.py` checks the separate scalar theorem with
projective absolute graph states and signed full-spin cross states. It verifies
equal relative-gauge fibers, constructs the gauge-maximizer injection, checks
the resulting order-statistic bound, and confirms by strict examples that the
fiber relation is domination rather than an exact max-plus identity:

```bash
python3 verification/verify_swapped_profile_injection.py
```

Expected output:

```text
split=2+2 normalized_cases=2 represented_labelled=64 gauges=16 states=128 strict_fiber_dominations=4 lambda_gain_pairs=0/0:1,0/2:1
split=2+3 normalized_cases=8 represented_labelled=1024 gauges=128 states=2048 strict_fiber_dominations=64 lambda_gain_pairs=2/4:8
split=2+4 normalized_cases=64 represented_labelled=32768 gauges=2048 states=65536 strict_fiber_dominations=688 lambda_gain_pairs=0/2:12,0/4:12,2/4:24,4/4:14,4/6:2
split=3+3 normalized_cases=64 represented_labelled=32768 gauges=2048 states=65536 strict_fiber_dominations=628 lambda_gain_pairs=2/6:24,4/6:40
normalized_cases_checked=138
represented_labelled_cases=66624
balanced_gauges_checked=4240
swapped_states_checked=133248
strict_fiber_dominations=1384
stream_sha256=638daefed306506cac5f7a724a64b4717d902601a80b2a5b7e73a9da768fbec9
corruption_controls=relative_sign_omission,max_plus_strictness
swapped_profile_injection=PASSED
```

All arithmetic is integer-exact. The exhaustive cases are switching-normalized;
the reported labelled-case count is recomputed from orbit multiplicities.

`verify_swapped_profile_floor.py` checks the separate all-order lower bound on
the swapped raw order statistic. It exhausts every rectangular signing through
`4 x 4`, then checks every switching-normalized balanced block triple through
order `3 + 3`:

```bash
python3 verification/verify_swapped_profile_floor.py
```

Expected output:

```text
rectangular_matrices_checked=66066
rectangular_pairs_checked=16810248
raw_balanced_cases_checked=66
raw_product_states_checked=65664
stream_sha256=9552e95262c5cbc683e64dfea427da2952cc3b686442af3146385f67c7567b32
corruption_controls=radius_one_multiplicity,incompatible_double_maximum
swapped_profile_floor=PASSED
```

The analytic theorem counts a radius-one Hamming ball around each maximizing
row spin. It implies the conditional asymptotic threshold
`sqrt(2/pi)/(2^(3/2)-1) = 0.436377...`; the finite exhaustion checks the rank,
multiplicities, and normalization, not the all-order proof.

## Labeled-shell Parseval refinement

`verify_labeled_shell_parseval.py` verifies the exact Fourier refinement of
the scalar profile. It checks the strict `2 + 4` gain, reconstructs the full
fiber occupancy for one balanced `C14` split by a shellwise Walsh transform,
and verifies Parseval in exact rational arithmetic:

```bash
python3 verification/verify_labeled_shell_parseval.py
```

Expected output:

```text
small_high=lambda:0 occupancy:0x8,1x16,2x8 fourier:2:16,12:16 V:1/2 certified_gain:2 true_gain:4
small_low=lambda:0 occupancy:0x8,1x16,2x8 fourier:21:-16,27:16 V:1/2 certified_gain:2 true_gain:2
c14_split=maxima:11,11,21 B:304908 K:8192 nontrivial_fourier:8159 V:635307383/4194304 occupancy:0..87 zero_fibers:1 l2_certifies:false
corruption_controls=relative_orientation,target_cutoff,walsh_inverse
arithmetic=integer,fraction
labeled_shell_parseval_verification=PASSED
```

The `2 + 4` scalar shell has exactly one fiber's worth of triples and hence
certifies no gain without labels. Its nonconstant Fourier mass certifies gain
two. For `C14`, empty fibers exist, but the generic variance inequality is too
weak because the mean occupancy is about `37.22`; exact inversion finds only
one empty fiber among 8192. This is a scoped failure of the `L2` bound, not of
the full labeled convolution.

## Labeled-shell moment certificate

`verify_labeled_shell_moment_certificate.py` independently reconstructs the
same `C14` occupancy by direct sparse XOR convolution, rather than Walsh
inversion. It verifies a degree-19 consecutive-root certificate and the
equivalent order-nine localizing-matrix witness in exact arithmetic:

```bash
python3 verification/verify_labeled_shell_moment_certificate.py
```

Expected output:

```text
c14_direct=maxima:11,11,21 B:304908 K:8192 occupancy:0x1,min_positive:6,max:87
degree19=adjacent_roots:9,17,26,36,46,56,67,76,86 expectation:1707454816960049615/99244391564512637853696 certified_empty_fibers:1
localizing_order:9 positive_one_minus_b_numerator:584163517696745929254421003286532
chebyshev_generic_order:14 range:6..87
corruption_controls=leading_sign,nonconsecutive_factor,target_cutoff
arithmetic=integer,fraction
labeled_shell_moment_certificate=PASSED
```

The positive expectation rigorously forces at least one empty fiber. Direct
reconstruction separately shows there is exactly one. The analytic moment
hierarchy is complete for finite occupancies.

`verify_labeled_vacancy_hierarchy.py` checks the sharp generic limitations of
that hierarchy on the same independently reconstructed occupancy:

```bash
python3 verification/verify_labeled_vacancy_hierarchy.py
```

Expected output:

```text
c14_vacancy=K:8192 B:304908 zeros:1 positive_range:6..87
soft_inverse_temperature=log(K) first_collision_degree=87
fourier_psd_safe_dimension=6826 maximum_nontrivial_coefficient=1223/512
normalized_localizing_margin=-36795384082687448302747606657/4816759830492505652837357886720000
corruption_controls=filled_vacancy,unshifted_fourier_kernel,fill_one_localizer
arithmetic=integer,fraction
labeled_vacancy_hierarchy=PASSED
```

The exact checks distinguish four statements: the unique-hole margin is at
most `1/K`; the canonical soft test needs inverse temperature `log(K)`; the
alternating collision expansion first succeeds at degree 87; and every
Fourier principal minor on at most 6826 characters is necessarily positive.

`verify_universal_moment_obstruction.py` checks the explicit six-point
vacant/nonvacant pair and the derivative-vector counts behind the general
box-principle obstruction:

```bash
python3 verification/verify_universal_moment_obstruction.py
```

Expected output:

```text
explicit_pair=K:6 matched_moments:0..2 cubic:837/873
derivative_count_cases=3631 simple_collision_cases=449 exact_collision_cases=469 exact_only_cases=20
corruption_controls=zero_removal,unequal_padding,zeroth_range
arithmetic=integer
universal_moment_obstruction=PASSED
```

The all-order theorem is analytic. The finite script checks its exact support
polynomial, padding step, derivative-range product, simpler sufficient bound,
and deterministic corruption controls.

## Fixed-half cut-discrepancy check

The fixed-density equivalence is an analytic theorem. Its finite checker
independently exhausts all signings through order six, recomputes `F(n)` and
`4H(n)`, checks the pointwise fixed-layer identity and switching mean-square
step, and deliberately corrupts the cut centering:

```bash
python3 verification/verify_cut_discrepancy_equivalence.py
```

Expected output:

```text
orders=2:F=1:4H=2,3:F=3:4H=4,4:F=4:4H=4,5:F=4:4H=4,6:F=5:4H=8
fixed_layer_pointwise_checks=6711
switching_mean_square_checks=33866
corruption_control=wrong_cut_centering_detected
cut_discrepancy_equivalence_verification=PASSED
```

All decisions use integer arithmetic. The finite exhaustion checks the proof's
normalization and edge cases; it is not the proof of the all-order theorem.

## Fixed-density rectangular cross floor

`verify_fixed_density_cross_floor.py` exhausts every rectangular signing
through `4 x 4`. It checks the exact XOR-cut identity, both Khintchine floors,
the fixed-total optimum, and the constructive switching-and-editing proof:

```bash
python3 verification/verify_fixed_density_cross_floor.py
```

Expected output:

```text
rectangles=1x1..4x4
xor_identity_checks=74954
khintchine_floor_checks=149908
fixed_total_minimum_checks=116
constructive_switch_edit_checks=116
minimum_profile_sha256=2f0e4e234d225e8b5ad4513212900132817ad323e91f3e6e0eb0874847a3c0f1
corruption_controls=wrong_half_centering_detected,wrong_total_parity_detected
fixed_density_cross_floor_verification=PASSED
```

The finite checker confirms normalization and edge cases. The all-order
rectangular floor and balancing inequality are proved analytically in the
research note.

## Fixed-half cloud amplification checks

`verify_equal_cloud_blowup.py` exhausts small fixed-half complete/empty cloud
blow-ups, every cloud-union cut, and switch-plus-rebalance instances. It also
checks the exceptional failure of fixed-half completion at base order two:

```bash
python3 verification/verify_equal_cloud_blowup.py
```

Expected output:

```text
orders_and_clouds=3x2,3x3,4x2
fixed_half_blowups_checked=26
cloud_union_cuts_checked=368
switch_rebalance_checks=6848
n2_infeasible_cloud_sizes=2..6
switch_profile_sha256=3e0d47a7b780432841530c2231881d0bddc7169c0abc9f5238bab1e9347ed535
corruption_controls=cross_edge,n2_feasibility
equal_cloud_blowup_verification=PASSED
```

The all-order `k^2` cloud-union bound and the `(r+t_N)/2` repair loss are
proved analytically; the script checks normalization, parity, and small edge
cases.

`verify_hadamard_cloud_lift.c` exhausts the uniform four-fold symmetric-
Hadamard lift of the order-five optimum:

```bash
cc -std=c11 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror \
  verification/verify_hadamard_cloud_lift.c \
  -o /tmp/verify_hadamard_cloud_lift
/tmp/verify_hadamard_cloud_lift
```

Expected output:

```text
symmetric_hadamards=64
trace_counts=-4:8,0:48,4:8
base_global_sign_reduction=0,2,4,1,3
representative_minima=trace0:44x4,trace4:48x8
maxima_table_fnv64=273ea01435d2c8a5
fixed_half_total_sign=0
fixed_half_minimum_attained=44
product_spin_maximum=32
fine_spin_witness=44
corruption_controls=hadamard_entry,base_anti_isomorphism,fixed_half_total,fine_spin_witness
hadamard_cloud_lift_verification=PASSED
```

This is a complete finite result for the common symmetric `H` and common
diagonal completion `D`.  The checked permutation `(0,2,4,1,3)` sends the
signed five-cycle base to its negative, which justifies quotienting the 64
Hadamards by global sign before the maximum table is exhausted.  It is not a
theorem about general orthogonal or cloud-dependent lifts.

`verify_cloud_dependent_hadamard_lift.py` checks a complementary exact
order-16 construction:

```bash
python3 verification/verify_cloud_dependent_hadamard_lift.py
```

Expected output:

```text
clouds=4x4
quotient_free_hadamard_blocks=6
fixed_half_internal_blocks=4
projective_spins_checked=32768
base_F4=4
cross_only_maximum=28
lift_maximum=32
projective_maximizers=14
energy_histogram_sha256=b4950bbfbac3a6223c6260572ce698fd54168777da04adfa317c95456baa38ca
order4_hadamards_checked=768
common_H_nonzero_cloud_orders=0132:2;0231:2;1023:2;1320:2;2013:2;2310:2;3102:2;3201:2
common_H_total_gauge_solutions=16
canonical_common_H_cloud_order=0132
canonical_common_H_is_symmetric=FALSE
canonical_base_signs=+,+,+,-,+,-
canonical_base_maximum=4
canonical_internal_maxima=6,6,4,4
common_H_common_D_up_to_sign_gauge_solutions=0
stored_orientation_fixed_template_profile=M4_to_38:8,M4_to_40:40,M6_to_32:8,M6_to_40:8
corruption_controls=hadamard_entry,lower_block_transpose,synthetic_common_H
cloud_dependent_hadamard_lift_verification=PASSED
```

Direct and blockwise evaluations agree on every projective spin. The
transpose-aware scan over all 24 cloud orders is essential: the construction
has a common nonsymmetric oriented Hadamard representation, and its intrinsic
edge signs are an `F(4)` optimizer. It is an exact one-step lift with four
tailored internal blocks, two of them deliberately nonoptimal. It is not a
uniform or iterable amplification theorem. The stored-orientation template
profile is only a coordinate-level diagnostic.

`verify_framed_hadamard_lift_30.py` checks the sharper completion of that
common oriented frame:

```bash
python3 verification/verify_framed_hadamard_lift_30.py
```

Expected output:

```text
order=16
cross_frame=common_oriented_H4
base_signs=+,+,+,-,+,-
internal_pattern=P,R,P,R
internal_maxima=4,4,4,4
P_R_switching_permutation_equivalent=TRUE
projective_spins_checked=32768
cross_only_maximum=28
lift_maximum=30
projective_maximizers=38
full_maximizers=76
matrix_sha256=352392c57458568ddbf2920d4cb487f67d21fb491ab8c75c6862e9c7fc6a9181
energy_histogram_sha256=875e4f931630501ec7730abad70df7ab029602f7a48b4b10d56df9de1e319388
six_state_obstruction_cycle=0->3->2->0
fixed_H_arbitrary_internal_minimum=30
common_literal_internal_minimum=38
common_literal_internal_profile=M38:6,M40:18,M42:32,M46:8
corruption_controls=hadamard_entry,transpose,six_state_pair
framed_hadamard_lift_30_verification=PASSED
```

The upper bound is an admissible order-16 signing evaluated independently by
direct and blockwise formulas. The fixed-frame lower bound is a separate
six-state certificate: maximum at most 28 would force a directed cycle of
three framed response inequalities, hence equality of response values whose
difference is twice a sum of three signs. The script also recomputes `F(4)`,
checks the switching-permutation equivalence of the two internal blocks,
exhausts the 64 literal common-internal completions, pins the full energy
histogram, and exercises deterministic corruption controls. This proves
`F(16) <= 30` and the fixed-frame minimum 30; it does not prove `F(16) = 30`.

An independent strict-C reconstruction checks the same matrix, complete
histogram, six-state obstruction, and common-internal profile:

```bash
cc -std=c11 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror \
  verification/verify_framed_hadamard_lift_30.c \
  -o /tmp/verify_framed_hadamard_lift_30
/tmp/verify_framed_hadamard_lift_30
```

Expected output:

```text
projective_spins_checked=32768
direct_block_histogram_match=TRUE
lift_maximum=30
projective_maximizers=38
cross_only_maximum=28
six_state_obstruction=PASSED
common_literal_internal_minimum=38
strict_c_framed_hadamard_lift_30_verification=PASSED
```

The new exact identity checker also needs only the standard library:

```bash
python3 verification/verify_new_results.py
python3 verification/verify_continuation.py
python3 verification/verify_coding_continuation.py
python3 verification/verify_amplification_obstructions.py
python3 verification/verify_cavity_hereditary.py
python3 verification/verify_nonlinear_bellman.py
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

The nonlinear-Bellman checker exhausts every signing through order 5 at four
covariance parameters, verifies the trace-four sum-of-squares identity and
the quantitative and optimized-Jensen arcsine bounds, and exhausts the parity
Gram-defect theorem through order 6. It compares the multivertex
weighted-radius identity against direct block enumeration in 23 small cases
and includes exact rational certificates for `F(20) >= 30` and
`F(21) >= 32`. Expected output:

    nonlinear_arcsine_checks=4392
    optimized_jensen_checks=4384
    parity_defect_signings_checked=33866
    multivertex_bellman_checks=23
    order_20_baseline=27.998982684199
    order_20_strengthened=28.018726463062
    order_20_parity_rounded=30
    order_21_old_bound=29.894026823521
    order_21_new_bound=30.148921557028
    order_21_parity_rounded=32
    ordinary_distance_corruption_detected=TRUE
    parity_defect_corruption_detected=TRUE
    deterministic_seed=413935
    nonlinear_bellman_verification=PASSED

Only the displayed arcsine values use floating point. Both strict finite
comparisons are separately certified using rational arithmetic.

`verify_frontier_walls.py` uses exact arithmetic to check the universal
weighted-entropy floor, exhausts all 32,768 projective cross centres for the
order-four block diagnostic, and checks 200,000 orders of the abstract
Bellman countermodel. Floating point is used only for the displayed logarithm
constant and the smooth countermodel phase. Expected output:

    entropy_floor_exact_checks=1401
    balanced_entropy_leading_gap=0.250895659942
    b4_cross_centers_checked=32768
    b4_orbit_shifts_checked=128
    b4_height_profile=0:14,2:24,4:16,6:8,8:2
    b4_iid_union_sums=K10:17973/4096,K12:6073/4096,K14:1653/4096
    b4_cross_optimum=10
    b4_optimal_projective_centers=92
    countermodel_orders_checked=200000
    countermodel_hereditary_checks=251500
    dyadic_identities_checked=100000
    maximum_increment_over_sqrt_n=1.264911064067
    deterministic_seed=413935
    corruption_controls=entropy_sign,adaptive_seed,union_threshold,increment_parity,dyadic_parity
    frontier_walls_verification=PASSED

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

## Complete order-9 weighted geometry

`research_order9_weighted_geometry.py` requires nauty `geng`, NetworkX, and
NumPy:

```bash
python3 verification/research_order9_weighted_geometry.py
```

It evaluates all 12,346 root-normalized order-9 signings. Every energy,
projective distance, weighted radius, and direct extension value is computed
in exact integer arrays. The script asserts the catalogue count and committed
nauty stream digest, samples an independent NetworkX graph6 decoder, and
groups the 55 optimum root representatives into switching-permutation
classes. It also verifies the two explicit non-optimal collision records.

Expected output:

    order_9_root_records=12346
    order_9_pareto={(12,0)}
    order_9_M12_root_pairs={(12,0):20,(12,1):35}
    order_9_M12_switching_classes={delta0:4,delta1:11}
    histogram_mixed_groups=10,records=874
    histogram_plus_maximizer_distance_mixed_groups=3,records=112
    energy_colored_two_point_mixed_groups=0
    collision=GHOgmo:(4,4,15),Gxd?Dc:(3,3,17)
    collision_histogram={2:124,6:85,10:37,14:10}
    collision_maximizer_pair_distance=(10,16,14,20,40)
    geng_stream_sha256=6b740e1c1ec4f6c7d5539e2e236da0f1ad6aa3120d534590b0ea1f09ddc0b345
    deterministic_seed=413935
    order_9_weighted_geometry=PASSED

Certificate boundary: completeness trusts nauty `geng`, with the exact count
and stream digest asserted. Switching-class counts trust NetworkX graph
isomorphism. The weighted-radius values and direct extension optima do not use
a SAT, MILP, or floating-point solver. Separation by the energy-coloured
two-point invariant is an order-9 finite fact only.

## Complete order-10 temperature phases

`research_order10_temperature.py` requires nauty `geng` and NumPy. It
reuses the independently audited graph6 parser from
`research_exact_small_n.py` and fail-closes on the complete stream count and
digest:

```bash
python3 verification/research_order10_temperature.py
```

It evaluates all 274,668 root-normalized records, obtains 6,012 distinct
absolute-energy histograms, and proves the three-phase envelope using exact
integer polynomial arithmetic. Expected output:

    order_10_root_records=274668
    order_10_absolute_energy_histograms=6012
    phase_records=1,4,4
    phase_maxima=15,15,13
    phase_transition_z=2,root(z^3+z^2-10z-8)
    positive_cubic_root=3.083872359436
    temperature_thresholds=0.658478948462,0.792460761565
    extensive_beta_thresholds=2.082293268414,2.505980962859
    geng_stream_sha256=ce9c5d4d27c8e55de5f0c6348ec781a650382e16bdff26b6c3418fa00a9cfcf9
    independent_recomputation=PASSED
    corruption_controls=histogram,coefficientwise_direction,phase_factorization
    order_10_temperature_phases=PASSED

Catalogue completeness trusts nauty plus the asserted count and digest. The
phase envelope and thresholds use no floating-point comparison except for
the displayed decimal approximations.

## Negative-replica transport experiment

`research_negative_replica_transport.py` requires nauty `geng` and `countg`.
It directly checks the exact supermultiplicative inequality in small labeled
orders, then weights every unlabeled residual graph by
`(n-1)!/|Aut(G)|` to reconstruct the uniform labeled-disorder moment through
order 9:

```bash
python3 verification/research_negative_replica_transport.py
```

Expected output:

    negative_replica_supermultiplicativity_checks=24
    transport_beta_1_theta_1_delta_over_n2=4:0.0049185348,5:0.0087837537,6:0.0116569500,7:0.0138659979,8:0.0156035779,9:0.0170007335
    transport_beta_2_theta_4_delta_over_n2=4:0.1484980817,5:0.2771315062,6:0.3275853627,7:0.2860893527,8:0.3076235046,9:0.3159700451
    transport_beta_4_theta_8_delta_over_n2=4:0.5774823306,5:1.0474135703,6:1.1793621066,7:0.9941740615,8:1.1008552230,9:1.1106325514
    transport_orders=4..9
    deterministic_seed=413935
    corruption_controls=unlabeled_multiplicity,t0_normalization
    negative_replica_transport_verification=PASSED

Energies, automorphism orders, and labeled multiplicities are exact. The
moment logarithms and displayed transport defects use floating point and are
finite evidence only. Completeness trusts nauty's streams and automorphism
orders. The analytic leading-gap theorem does not depend on this script.

## Negative-replica transport obstruction and alignment

`verify_negative_replica_transport_obstruction.py` checks the constants in
the analytic disproof of power-saving transport, the exact Boolean-cube
entropy-production decomposition, and its one-coordinate stability
remainder:

```bash
python3 verification/verify_negative_replica_transport_obstruction.py
```

Expected output:

    scalar_stability_checks=3600
    cube_entropy_production_checks=240
    maximum_derivative_error=1.468e-09
    minimum_scalar_grid_slack=1.998e-11
    PT_liminf_beta4_theta8=0.027919201253204
    PT_beta_threshold=3.012373175204
    PT_theta_threshold_at_beta4=7.690245425859
    deterministic_seed=413935
    corruption_controls=scalar_sign,orbit_entropy
    negative_replica_transport_obstruction_verification=PASSED

The scalar and cube grids are regression checks, not proofs. The proof of the
stability inequality is the elementary one-coordinate argument in the note.

`verify_negative_replica_alignment.py` uses only integer and rational
arithmetic. It checks the quotient dimensions, conditional-density
normalization, the exact `2 + 4` scalar-state collision, and the mixed
four-cycle trace/Walsh identity:

```bash
python3 verification/verify_negative_replica_alignment.py
```

Expected output:

    code_chain_checks=16
    universal_mixed_variance_checks=74240
    K2_high=196585091273040100817/133610891512185651200
    K2_low=6723290161/5922841600
    H4_high={-12: 4, -4: 4, 0: 16, 4: 4, 12: 4}
    H4_low={-4: 8, 0: 16, 4: 8}
    collision_true_fiber_states=16
    arithmetic=integer,fraction
    corruption_controls=edge_order,density_mean,K2_separation,H4_trace
    negative_replica_alignment_verification=PASSED

The two witnesses have identical complete local graph and rectangular
absolute-energy histograms, so their three local scalar partition functions
agree at every temperature. Their different rational `K2` values are an exact
conditional-alignment distinction. Because `dim(D_2)=1`, the verifier's 32
relative-gauge tuples cover the true 16-point fiber twice. It checks this
twofold redundancy explicitly; normalized fiber moments are unchanged.

## Square-order Paley Boolean eigenvectors

`verify_paley_subfield.py` is a standard-library exact checker for prime
subfields `m=3,5,7`:

```bash
python3 verification/verify_paley_subfield.py
```

Expected output:

    paley_subfield_witnesses=m=3:N=10:M=15,m=5:N=26:M=65,m=7:N=50:M=175
    corruption_controls=coset_constancy
    paley_subfield_verification=PASSED

It constructs `GF(m^2)` directly, verifies the conference identity and the
Boolean eigenvector equation, and evaluates the witnessed quadratic energy.
The paper's proof covers every odd prime power; the script deliberately avoids
a finite-field dependency and checks prime `m` only.

## Exact values at orders 11 and 12

`research_order11_certify.py` requires nauty `geng`. Its default reduced pass
uses the exact all-positive, complement, and singleton reductions:

```bash
python3 verification/research_order11_certify.py
```

Expected output:

    geng_mode=filtered
    geng_records=2153606
    eligible_order11_records=2153606
    order11_subset_checks=48173339
    order11_M_le_15_survivors=0
    independently_crosschecked_edge_counts=20,21,22
    geng_stream_sha256=b62da4d7ebfaab4ccd801fd509f2fc85f6c2b815c8c1d2e969de7aa6a82c322d
    order11_witness=ICRbczQMo maximum=17 maximizers=5
    order12_witness=JWUuDOR\K{? maximum=18 maximizers=20
    certified_values=F(11)=17,F(12)=18
    deterministic_seed=413935
    corruption_controls=graph6_ordering,absolute_value,false_bound,energy_parity
    order11_order12_verification=PASSED

The slower full-stream mode is:

```bash
python3 verification/research_order11_certify.py --full-stream
```

It scans all 12,005,168 unlabeled graphs on ten residual vertices, internally
selects the same 2,153,606 eligible records for cut evaluation, and requires full stream SHA-256
`5650c7c979fdffd8c0f99a2f2ee8775938ec2a3dd69aa65be1207936824fc5b3`.
Every cut decision and both witness replays use integers. Deterministic samples
at each eligible edge count are recomputed with a separate adjacency formula.
Completeness trusts
nauty to emit one representative of every isomorphism class; no MILP, SAT,
floating-point, or timeout result is used.

## Exact values at orders 13 and 14

`research_order13_certify.py` uses only the Python standard library plus the
local conference checker. It compiles `order12_threshold_scan.c` from source
with warnings as errors; supplied executables are never used. The quick pass
checks the two rooted survivors, every projective spin, every projective
incident column, all vertex deletions, the Paley witnesses, and three
corruption controls:

```bash
python3 verification/research_order13_certify.py
```

Expected output:

```text
order12_survivor=JCpVdXyxpz? M=18 maximizers=20 extension_minimum=24 optimal_centers=772 deletion_maxima=17x12
order12_survivor=JCpdUg{[dM? M=18 maximizers=20 extension_minimum=24 optimal_centers=772 deletion_maxima=17x12
paley_C14_M=21 principal_order13_maxima=20x14
corruption_controls=graph6_padding,empty_residual,edge_flip
order13_order14_quick_verification=PASSED
deterministic_seed=413935
```

The full certificate requires nauty `geng`. It generates all eleven-vertex
residual graphs in eight disjoint shards, relays every byte to the scanner,
and checks the producer and consumer statuses, byte counts, record counts,
SHA-256 digests, and survivor sets:

```bash
python3 verification/research_order13_certify.py \
  --full-stream --jobs 8 --geng /absolute/path/to/geng
```

The asserted shard receipts are:

```text
shard  records      bytes       sha256                                                            survivors
0      119431209    1433174508  0ae80a506fc9ef5aca212fb41eea05a453dd8326dc078fb1e2228b0c96cc2d5c  -
1      128496882    1541962584  54cd49e34942d1905213b802bff7adc565bd2c3251ac2268414ea179857559c5  -
2      128472053    1541664636  92ad1076355bd829ed0596dc70b9203dcffd88b5a4bc549de84aa37e6ceea318  JCpVdXyxpz?
3      121592284    1459107408  946fdda70e727892256a9e10c6e8d9823f2c1366f867499baa9a451b5d46e2e8  -
4      119556409    1434676908  b3dbbe8d186b4be67dbb119d41f34cb1637444bfa0207710549b9d04557a2798  JCpdUg{[dM?
5      134239743    1610876916  6acefb47e065175ccb85bb157d4cc50d6b1c8519d782146f2c449bc2582c2520  -
6      143004566    1716054792  6e2d97e700e524f244c14ddd1f0a832404b1661366bbebf4d5324640efc3d275  -
7      124204718    1490456616  73551d55bb414836abc29cf2700aee4552c9d476210cfba156434a9530dcc522  -
```

Their record counts sum to 1,018,997,864. Exactly two rooted records have
order-12 maximum at most 18; both have extension minimum 24. Every other
predecessor has maximum at least 20. Together with the verified Paley
witnesses this certifies `F(13)=20` and `F(14)=21`. The explicit remaining
trust boundary is nauty's completeness for unlabeled graphs. The full local
receipt is `order13_full_receipt.txt`; it is documentary only and is not read
by the verifier.

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

---
id: HYP-3245
title: LRC14 equioscillation / autocorrelation atlas
status: EVIDENCE / exact motif scout plus synthesis; not an LRC14 proof
source: codex-2026-06-28
tangent: T1309
technique: LTI-309
tournament_technique: LTT-209
script: 04-computation/lrc14_equioscillation_autocorrelation_atlas_codex_20260628.py
result: 05-knowledge/results/lrc14_equioscillation_autocorrelation_atlas_codex_20260628.out
reflection: 07-reflections/lrc14-equioscillation-autocorrelation-atlas-codex-20260628.md
related:
  - HYP-3244
  - HYP-3243
  - HYP-3242
  - HYP-3241
  - HYP-3240
  - HYP-3239
  - HYP-3238
  - HYP-3237
  - HYP-3236
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3230
  - HYP-3231
  - HYP-3229
  - HYP-3228
  - HYP-3227
  - HYP-3226
  - HYP-3225
  - HYP-3224
  - HYP-3223
  - HYP-3222
  - HYP-3219
  - HYP-3218
  - HYP-3217
  - HYP-3214
  - HYP-3213
  - HYP-3212
  - HYP-3205
  - HYP-3204
  - HYP-3203
  - HYP-3202
  - HYP-3163
  - HYP-3132
  - OPEN-Q-108
---

# HYP-3245: LRC14 Equioscillation / Autocorrelation Atlas

## Claim

HYP-3214 gives the exact bridge:

```text
Fejer F_7
  = (de Moivre cubic)^2
  = Chebyshev/V_7 equioscillation kernel
  = positive-definite Delsarte kernel
  = triangular autocorrelation of the AP interval
```

HYP-3245 extends this into a reusable principle:

```text
equioscillation is the primal/root-side shadow;
autocorrelation is the Fourier/Gram-side certificate.
```

For LRC14, the new concrete out-correlation signal is that every non-AP
HYP-3202 trap moves ordinary speed-support autocorrelation mass outward:
short-lag mass decreases and long-lag mass increases by the same amount.  The
coarse out-correlation transport is universal across traps, while HYP-3225's
Green/Rayleigh/Plucker sidecars explain the finer ripple type.

Post-rebase integration: incoming HYP-3227 is the conductance-graph companion.
It proves the same finite trap boundary remains connected to legal certificate
coordinates, even after deleting Toeplitz and even through Green-only
coordinates.  HYP-3245 should therefore be read as the autocorrelation
projection of that trap-discharge graph: outward pair-mass transport is the
lag-space shadow of conductance to normal-fan exits.

Second post-rebase integration: incoming HYP-3228 is the shell-magic
companion.  HYP-3228 packages the k=8 `L_y` dual as the finite shell vector
`[10,0,0,1,0,0,10]`; HYP-3245 supplies the lag/Gram transport language that
could turn its shell slack into AP support, Toeplitz, conductance, or
ordered-tail exits.  The split is useful because positive Fejer/PSD sector
mass, shell contact, and ordinary speed-support autocorrelation are nearby but
not identical clocks.

Third post-rebase integration: incoming HYP-3229, HYP-3230, and HYP-3231
sharpen what to measure next.  HYP-3229 keeps Gamma0(7), divisor-fiber
coefficients, and subshift autocorrelation as labelled sidecars for the
Fejer/Gram route.  HYP-3230 supplies the three-gap/Farey recursion of the
cap-kernel clock, while HYP-3231 records the scale-normal recursion ledger for
legal proof states.  Therefore HYP-3245 should add three transport fields in
the next scout: a modular-coefficient sidecar for each lag profile, a
three-gap/Farey breakpoint label for each outward-lag surplus, and a
scale-normal survival bit saying which quotient coordinates remain legal after
rescaling.

Fourth post-rebase integration: incoming HYP-3232 and HYP-3217 explain why
the new `scale_survival_bit` and `contact_word` fields should be first-class.
HYP-3232 locates the cap-kernel difficulty at the Eisenstein fold where
modulus-covariance breaks at the apex half: speeds below the apex carry clean
three-gap scale behavior, while the antipode half produces the binding
deviation.  HYP-3217 says the Mobius/Eisenstein/Legendre modes sit inside the
`Q(zeta_7)` subfield tower, with the cubic de Moivre mode as the missing
degree-3 contact language and Fejer as the sextic square.  Therefore HYP-3245
should tag each lag-transport residual by `apex_fold_side` and
`cyclotomic_mode` before using it as proof evidence.

Fifth post-rebase integration: incoming HYP-3233 and HYP-3234 sharpen the same
labels in two compatible ways.  HYP-3234's signed address chart-change sheaf
says the Mobius/Eisenstein/Legendre recursions are local charts whose letter
names cannot be globally reused without a chart-change map.  HYP-3233's
cyclotomic-factor packet says each mode has a factor signature
`(x-1)^depth * Phi_d`, with the Phi_7/de Moivre factor carrying the apex
hardness.  Therefore HYP-3245 should also tag outward-lag residuals by
`chart_change_class` and `cyclotomic_factor_signature`; a lag transport that
looks invariant in one chart may be hiding a cancelled-slot or Phi_7 debt in
another.

Sixth post-rebase integration: incoming HYP-3235 and HYP-3218 turn the
autocorrelation transport into a field/conductor test.  HYP-3235 says the cap
and binding dip live in `Q(cos(2pi/7))`, with conductor `7^2` appearing exactly
in the binding rows, and that `F_7=(de Moivre cubic)^2` is a totally-positive
cyclotomic square.  HYP-3218 says the equidistribution finish is the explicit
Fejer/Vaaler certificate: Lee-Yang margin `phi(7)`, Gauss-sum norm `7`, and AP
self-duality are one critical surface.  Therefore HYP-3245 should add
`cap_field_conductor`, `fejer_square_status`, `gauss_sum_margin`, and
`ap_self_dual_fixed_point_status` to the next lag-transport scout.

Seventh post-rebase integration: upstream HYP-3236, HYP-3237, and HYP-3219
make the same out-correlation packet three-dimensional.  HYP-3236 says the
finite trap boundary has an all-bank Green conductance projection: AP and
doubled AP maximize `lambda2` and minimize Green resistance, while every
non-AP trap has resistance excess.  HYP-3237 says the level-7
equioscillation is the core side of the Vitali wall: measure handles the
bulk, while primitive `Phi_14` witnesses and the Brouwer max-min saddle handle
the measure-zero core.  HYP-3219 says the odd obstruction should split into a
Brouwer trace sign times an SOS magnitude.  Therefore HYP-3245 should add
`green_resistance_slack`, `lambda2_conductance_rank`,
`negative_covariance_leakage`, `thomson_current_profile`,
`fiedler_bottleneck_id`, `vitali_wall_side`, `brouwer_saddle_sign`,
`phi14_core_witness`, `core_bulk_transport_status`, `brouwer_trace_sign`,
`degree_sos_factorization`, and `even_odd_bonferroni_node_slack` before using
out-correlation as a terminal scalar.

Eighth post-rebase integration: upstream HYP-3238 and HYP-3239 show that the
out-correlation packet must also remember the symmetry type of the lost
coordinate.  HYP-3238 recasts the frontier as an even-positive /
odd-negative compression audit: Fejer squares, covariance, and positive Green
conductance are legal only when the odd/negative payload is zero,
reconstructible, dual-annihilated, or retained.  HYP-3239 refines the
topological side from Brouwer to Borsuk-Ulam at `n=14`: `14=|D_7|`, the de
Moivre cubic is the three two-dimensional `D_7` irreps, and the complement is
an orientation-reversing anti-automorphism because `7=3 mod 4`.  Therefore
HYP-3245 should add `positive_negative_duality_status`,
`odd_negative_payload_reconstruction`, `dihedral_irrep_label`,
`complement_antiautomorphism_sign`, `borsuk_ulam_index`,
`imaginary_gauss_sum_sign`, and `phi4_bimodal_extremizer_rank`.  The new
guess is that outward autocorrelation transport is the lag-space shadow of a
sign-isotypic defect unless these labels show the odd/negative coordinate has
already been discharged.

Ninth post-rebase integration: upstream HYP-3240 and HYP-3241 turn the core
side from a qualitative label into witness arithmetic.  HYP-3241 says the
`Phi_14` core witnesses form three antipodal pairs, and the count is the
equioscillation saddle index `(p-1)/2`; AP and the Goddyn-Wong tight row share
that same primitive witness set.  HYP-3240 says covering-tight dilations
`d*{1..13}` are constructed by promoted witnesses `t=1/(14d)` in
`Phi_{14d}`, while the binding dip is not a single `Q(sqrt(-7))` norm scalar.
Therefore HYP-3245 should add `equioscillation_saddle_index`,
`phi14_core_universality_status`, `dilation_witness_grid`,
`core_witness_break_reason`, and `imaginary_norm_route_status` before using
outward autocorrelation transport as a core-side proof witness.

Tenth post-rebase integration: upstream HYP-3242 says the cap is the measured
Euler characteristic of the danger-cover nerve, the lonely point is the cover
hole, and the Borsuk-Ulam antipodal witness pair is the topological face of
the same sign obstruction.  Therefore HYP-3245 should add
`danger_nerve_euler_characteristic`, `lonely_hole_status`,
`cech_betti_sidecar`, `topological_shadow_class`, and
`cover_hole_witness_pair`.  The new guess is that outward autocorrelation
transport is a one-dimensional shadow of hole creation/destruction in the
danger nerve: AP keeps the hole aligned with Fejer/equioscillation contact,
while traps move lag mass outward and must prove that the cover hole survives
or name the topological shadow that failed.

Eleventh post-rebase integration: upstream HYP-3243 packages the visual proof
routes as typed carriers: circle endpoint arrangements, oriented topes and
cocircuits, Cech safe-component nerves, `D_7` index packets,
`Phi_14/Phi_{14d}` witness strata, Green graphs, Toeplitz/Fejer normal-fan
faces, Lee-Yang root motion, ear payload graphs, finite chamber atlases, and
state-lift obligations.  Its finite theorem schema is the missing legality
gate for this atlas:

```text
primitive row
  -> open safe tope
     or AP/GW Phi_14 equality
     or dilation Phi_{14d} equality
     or finite Toeplitz/Green/root-motion chamber discharge
     or state-lift H=7 contradiction
     or named residual debt.
```

Therefore HYP-3245 should add `oriented_matroid_tope_status`,
`circle_endpoint_arrangement_cell`, `cech_safe_component_rank`,
`finite_chamber_schema_status`, `state_lift_H7_obstruction`, and
`proof_carrier_tournament_rank`.  The new guardrail is that
out-correlation transport is a theorem-facing coordinate only after this
carrier schema says which lost endpoint, sign, witness, chamber, or state-lift
payload has been retained or legally discharged.

Twelfth post-rebase integration: upstream HYP-3226 is now executed rather
than merely reserved.  Its motif ledger ranks the comb-overlap Gram kernel
M073, AP self-dual Fejer equidistribution M093, danger-cover nerve hole M102,
and finite chamber carrier atlas M103 near the proof-facing top.  This makes
HYP-3245's motif list testable: Fejer/AP autocorrelation belongs to the
M073/M093 payload cluster, while out-correlation transport is a candidate
finite-boundary column only if it states its HYP-3226 `payload_preserved`,
`payload_destroyed`, `repair_sidecar`, and `terminal_risk_label`.  The new
guardrail is the HYP-3226 no-free-slider rule: Barker, Golay, Zadoff-Chu,
Welch, difference-set, Remez, and tournament motifs are not proof objects
until they enter as typed payload atoms rather than famous-pattern names.

Thirteenth post-rebase integration: upstream HYP-3244 supplies the
tiling/half-tiling controlled-forgetting span.  That changes the status of
out-correlation residuals: a lag or motif quotient is proof-grade only after
it declares whether the witness was built in the fixed-Hamiltonian-path
tiling cover, whether it descends through the half-tiling/orbit quotient, and
which sidecar carries the descent.  The sidecars are path-presentation fiber
`H(T)/|Aut(T)|`, parent-automorphism incident-word orbit,
rectangle/hourglass coboundary residue, tail/tip deletion signature, and the
n=4 canary/filler `S`-fiber.  The new guardrail is controlled forgetting:
scalar lag transport may be a true AP-facing signal, or it may be a failed
descent certificate whose lost witness lives in the tiling cover.

## Exact Fejer Bridge

The script rechecks HYP-3214 numerically and then pins the autocorrelation
interpretation:

```text
max |(de_moivre_cubic(2cos t))^2 - F_7(t)| = 7.923e-12
max F_7(2pi*j/7), j=1..6 = 2.868e-30
max central slope at those zeros = 1.351e-10
F_7 coefficients = (7-|lag|)_+
interval_7 autocorrelation equals F_7 coefficients: True
interval_8 autocorrelation coefficients = (8-|lag|)_+
```

Thus the double-zero/equioscillation condition is the same object as a
triangular autocorrelation law.  This is the operational content behind the
phrase "AP autocorrelation is Fejer."

HYP-3214 also adds the necessary guardrail: the naive pairing of the orbit
Fejer kernel with the sector Fejer kernel does not equal the cap.  The
7-sector Fejer function is the coverage-side magic function, while the cap
side lives on the 14-clock Johnson/pair-Pascal scheme.

## Multitude Atlas

The script records several analogous motifs:

```text
Fejer_AP_interval:
  double zeros / Chebyshev sharpness <-> triangular autocorrelation

Johnson_pair_Pascal:
  14-clock cap face <-> Johnson-scheme pair distribution

Christoffel_Darboux_OPUC:
  reproducing-kernel peaks <-> Toeplitz moment / Verblunsky data

Welch_ETF_simplex:
  equal Gram entries <-> Welch frame-potential equality

difference_set_two_level:
  all nonzero cyclic lags equal <-> two-level periodic autocorrelation

CAZAC_Zadoff_Chu:
  constant magnitude spectrum <-> zero nonzero periodic autocorrelation

Golay_complementary_pair:
  paired ripple cancellation <-> summed aperiodic sidelobes vanish

Barker_low_sidelobe:
  small side-lobe ripple <-> bounded aperiodic autocorrelation

equiripple_FIR_Remez:
  alternating error extrema <-> magnitude-square filter autocorrelation

tournament_out_correlation:
  score/outdegree balance <-> out-neighborhood overlap A A^T
```

The exact computed examples include a Fano `(7,3,1)` difference set with
periodic autocorrelation `[3,1,1,1,1,1,1]`, a Zadoff-Chu length-7 sequence
with nonzero periodic autocorrelation at numerical zero, a Golay length-8
pair with summed sidelobes all zero, a Barker length-7 sequence with
side-lobes in `{0,-1}`, and the regular simplex meeting the Welch equality.

## Abstract Signal Dictionary

The shared pattern is stronger than "flatness."  Equioscillation is a
contact-side statement: a signed residual hits alternating or repeated tight
nodes.  Autocorrelation is the same payload after the map
`h -> h * reverse(h)`, hence a Gram/Bochner statement with nonnegative
Fourier mass.  Out-correlation is the defect when total pair mass survives
but the short-lag/AP part is exported to lags that need named sidecars.

This suggests new proof-facing fields:

```text
contact_word:
  ordered signs/multiplicities at equality or root-contact nodes

lag_barycenter:
  center of mass of R_E(d)_+ minus R_E(d)_-

transport_cost:
  earthmover distance from AP triangular autocorrelation to E autocorrelation

Fejer_annihilator_projection:
  component of the lag residual killed by the F_7 double-zero kernel

shell_lag_commutator:
  difference between HYP-3228 shell slack and HYP-3245 lag-transport slack

sidecar_entropy:
  number of Green/Rayleigh/Plucker/modular/scale labels needed to make the
  autocorrelation quotient legal

small_pattern_motif_id:
  HYP-3226 motif id carrying the lag signal as a typed payload atom

payload_preserved:
  LRC coordinate retained by the motif or autocorrelation projection

payload_destroyed:
  coordinate lost by scalarizing the motif or lag residual

repair_sidecar:
  sidecar that restores the destroyed coordinate or names terminal debt

terminal_risk_label:
  HYP-3226 direct / sidecar / analogy / raw risk tag for the motif

scale_survival_bit:
  whether the signal survives the HYP-3231 scale-normal quotient or must be
  routed to named destroyed-coordinate debt

apex_fold_side:
  below-apex clean three-gap sector or antipode-half Eisenstein-break sector

cyclotomic_mode:
  Mobius / Eisenstein / Legendre / cubic de Moivre / sextic Fejer payload

chart_change_class:
  local signed-address chart and cancelled-slot debt needed to compare modes

cyclotomic_factor_signature:
  `(x-1)^depth * Phi_d` factor carried by the residual projection

cap_field_conductor:
  totally-real field/conductor tag for the binding row

fejer_square_status:
  whether the residual respects the totally-positive Fejer square certificate

gauss_sum_margin:
  HYP-3218 Fejer/Vaaler margin seen by the residual

ap_self_dual_fixed_point_status:
  whether the AP self-dual fixed point survives the quotient

green_resistance_slack:
  HYP-3236 effective-resistance excess relative to AP

lambda2_conductance_rank:
  algebraic-connectivity rank of the row in the conductance projection

negative_covariance_leakage:
  covariance sign data destroyed by positive-part graph compression

thomson_current_profile:
  unit-current demand profile through the trap bottleneck

fiedler_bottleneck_id:
  Fiedler cut/channel that first detects the out-correlation trap

vitali_wall_side:
  bulk-measure side or cyclotomic-core side of the witness packet

brouwer_saddle_sign:
  topological orientation of the max-min equioscillation saddle

phi14_core_witness:
  primitive 14th-root witness status for the core row

core_bulk_transport_status:
  whether the information lives in the measure bulk or the cyclotomic core

brouwer_trace_sign:
  HYP-3219 sign of the odd obstruction

degree_sos_factorization:
  Brouwer degree/trace times SOS magnitude factorization status

even_odd_bonferroni_node_slack:
  local slack at the even/odd Bonferroni split

positive_negative_duality_status:
  HYP-3238 even-positive / odd-negative packet split

odd_negative_payload_reconstruction:
  whether clipped odd/negative data is zero, restored, or dual-annihilated

dihedral_irrep_label:
  HYP-3239 `D_7` irrep carrying the residual

complement_antiautomorphism_sign:
  whether complement acts as orientation-reversing anti-automorphism

borsuk_ulam_index:
  free `Z_2` odd-degree certificate tag

imaginary_gauss_sum_sign:
  `i*sqrt(7)` obstruction / sign-irrep orientation tag

phi4_bimodal_extremizer_rank:
  rank of the row under the bimodal/phi4 extremizer package

equioscillation_saddle_index:
  HYP-3241 antipodal `Phi_14` witness-pair count, equal to `(p-1)/2`

phi14_core_universality_status:
  whether AP and Goddyn-Wong share the same primitive core witnesses

dilation_witness_grid:
  HYP-3240 promoted `Phi_{14d}` witness grid for dilation rows

core_witness_break_reason:
  why the base `Phi_14` witness fails or survives for the row

imaginary_norm_route_status:
  whether the `Q(sqrt(-7))` norm scalar route is legal or refuted

danger_nerve_euler_characteristic:
  HYP-3242 measured Euler characteristic of the danger-cover nerve

lonely_hole_status:
  whether the witness is carried by a cover hole rather than measure bulk

cech_betti_sidecar:
  Cech/Betti label for the danger-cover obstruction

topological_shadow_class:
  whether the active shadow is combinatorial, homological, representation,
  algebraic, or Diophantine

cover_hole_witness_pair:
  Borsuk-Ulam antipodal witness pair attached to the cover hole

oriented_matroid_tope_status:
  HYP-3243 open-safe tope or boundary cocircuit state

circle_endpoint_arrangement_cell:
  endpoint/wall-crossing cell carrying the lag residual

cech_safe_component_rank:
  safe-component nerve rank after open/closed boundary separation

finite_chamber_schema_status:
  whether the row is open safe, known witness equality, finite chamber
  discharge, state-lift contradiction, or named residual debt

state_lift_H7_obstruction:
  whether a residual atom lifts to the forbidden `H=7` endpoint

proof_carrier_tournament_rank:
  HYP-3243 carrier priority before scalar lag transport is trusted

tiling_witness_lift_status:
  whether the lag/motif witness was built in the HYP-3244 fixed-path tiling
  cover before quotienting

half_tiling_descent_certificate:
  automorphism/fiber/coboundary proof that the half-tiling quotient preserves
  the target LRC predicate

path_presentation_fiber_weight:
  `H(T)/|Aut(T)|` or named fiber debt attached to the carrier

parent_aut_word_orbit_id:
  incident-word orbit under parent automorphisms

rectangle_hourglass_residue:
  GF(2) coboundary residue detecting illegal quotient compression

tail_tip_deletion_signature:
  tip/tail deletion payload retained through the tiling-to-half-tiling span

controlled_forgetting_span_status:
  lift/compress/descent/fail-debt status for the motif or residual
```

The bold hypothesis is that AP extremality is the unique zero-commutator
point: its contact word, Fejer annihilator, triangular lag law, shell magic,
and scale-normal cap recursion all align.  Non-AP traps are not random
failures; they are controlled commutators between these projections.

## LRC Trap Out-Correlation

For a speed support `E`, define ordinary difference autocorrelation:

```text
A_E(d) = |{x in E : x+d in E}|
```

For the AP row `E={0,1,...,7}`, this is the triangular Fejer coefficient row:

```text
A_AP(d) = 8-d for 0 <= d <= 7, and 0 for d > 7.
```

For every non-AP HYP-3202 trap, the residual

```text
R_E(d) = A_E(d) - A_AP(d)
```

has:

```text
sum_{d=1}^7 R_E(d) < 0
sum_{d=8}^{14} R_E(d) > 0
sum_{d=1}^{14} R_E(d) = 0
```

So the finite trap rows do not merely have "less AP autocorrelation."  They
move correlation outward.  This is a clean reading of the user's
out-correlation prompt: the obstruction is a transport of pair mass from
short lags, where AP saturates the Fejer/interval law, to outward lags, where
local exchange and endpoint sidecars become visible.

Representative rows:

```text
green split row (0,1,2,3,11,12,13,14):
  low_lag_deficit=-16, out_lag_surplus=+16, sign_changes=1
  residual=[0,-1,-2,-3,-4,-3,-2,-1,1,2,3,4,3,2,1]

rank-2 Plucker row (0,1,3,5,7,9,11,13):
  low_lag_deficit=-9, out_lag_surplus=+9, sign_changes=3
  residual=[0,-6,0,-4,1,-2,2,0,3,1,2,1,1,1,0]

mixed sidecar row (0,8,9,10,11,12,13,14):
  low_lag_deficit=-7, out_lag_surplus=+7, sign_changes=1
  residual=[0,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1]
```

Class averages:

```text
rank2_pair_plucker_bottleneck: count=6, avg_L1=20.000, avg_sign_changes=3.333
green_low_connectivity_bottleneck: count=2, avg_L1=28.000, avg_sign_changes=2.000
AFM_frustrated_high_rayleigh_debt: count=2, avg_L1=23.000, avg_sign_changes=2.000
mixed_green_lorentzian_sidecar: count=1, avg_L1=14.000, avg_sign_changes=1.000
```

This is not a replacement for HYP-3225.  It is a new coarse projection of the
same finite boundary.  The out-correlation residual detects the AP interval
autocorrelation failure; Green/Rayleigh/Plucker sidecars explain the legal
proof route.

## Proof Route

HYP-3245 suggests a three-level certificate:

```text
Level 1: Fejer / equioscillation
  AP interval has triangular autocorrelation and double-zero Fejer sharpness.

Level 2: out-correlation transport
  Every non-AP trap shifts finite speed-support autocorrelation mass outward.

Level 3: sidecar discharge
  Toeplitz/Green/Rayleigh/Plucker/HB and HYP-3244 tiling-descent sidecars
  route the outward mass into the HYP-3224 normal-fan certificate rather than
  treating it as a raw scalar.
```

The formal target is not "autocorrelation distance from AP is positive."
That is too weak.  The target is a signed transport lemma:

```text
low-lag deficit + outward-lag surplus
  -> AP support gap / Toeplitz slack / core witness index /
     dilation witness sidecar / cover-hole survival /
     finite proof-carrier chamber status / state-lift obstruction /
     tiling lift and half-tiling descent status /
     ordered-tail q0+q6 loss
```

If such a lemma exists, HYP-3214's Fejer kernel is the equality model,
HYP-3245's out-correlation table is the finite trap projection, and HYP-3225
supplies the boundary sidecar types.

## Tournament Analysis

Vertices are certificate/motif families, not runners or raw arcs.  The
pairwise observable is which motif retains both the equioscillation payload and
the autocorrelation/positive-definite payload while preserving the LRC
predicate.  The measured tournament is transitive:

```text
score_hist={100:1,95:1,91:1,86:1,80:1,76:1,72:1,65:1,60:1,45:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
priority_path=
  Fejer_AP_interval ->
  Johnson_pair_Pascal ->
  Christoffel_Darboux_OPUC ->
  Welch_ETF_simplex ->
  difference_set_two_level ->
  CAZAC_Zadoff_Chu ->
  Golay_complementary_pair ->
  Barker_low_sidelobe ->
  equiripple_FIR_Remez ->
  tournament_out_correlation
```

## Assumption Challenge

Alternate vertices considered:

```text
runners
sectors
Fourier modes
autocorrelation lags
zeros/equioscillation nodes
Toeplitz moments
trap sidecars
Johnson pairs
difference-set blocks
codewords
filter bands
proof obligations
oriented-matroid topes/cocircuits
circle endpoint cells
Cech safe-component nerves
finite chamber atlases
state-lift obligations
flip bits
deletion roots
rectangle/hourglass defects
parent-automorphism word orbits
tiling-cover fibers
```

Chosen vertices: certificate/motif families plus LRC trap autocorrelation
residuals.

Preserved predicate: HYP-3214 Fejer magic-function identity, HYP-3224
normal-fan AP extremality, and HYP-3225 finite trap sidecar discharge.

Destroyed information: exact LRC safe-set geometry, PGF root radii, endpoint
owners, and analytic equidistribution obligations unless routed to named
sidecars; after HYP-3243 this also includes finite chamber address, open
versus closed boundary status, and state-lift payload; after HYP-3244 it also
includes Hamiltonian-path presentation, canary/filler fiber mass, parent
word-orbit identity, coboundary residue, and tail/tip deletion payload.

Challenged assumption: equioscillation is only a root-location story.  The
atlas treats it as the primal shadow of autocorrelation and
positive-definiteness.

## Next Tests

1. Try to prove the short-lag to outward-lag transport lemma for the bounded
   k=8 bank, not only the `12` HYP-3202 traps.
2. Test whether low-lag deficit controls AP support gap or Toeplitz slack
   after adding HYP-3225 sidecar classes.
3. Add Johnson/pair-Pascal autocorrelation fields so the 7-sector Fejer
   coverage side and 14-clock cap side can be compared without the false
   naive Fejer pairing.
4. Bind the motif rows to executed HYP-3226 ids with explicit
   `small_pattern_motif_id`, `payload_preserved`, `payload_destroyed`,
   `repair_sidecar`, and `terminal_risk_label` fields.
5. Search for a Fejer-Riesz or Christoffel-Darboux factor whose coefficients
   turn outward lag surplus into HYP-3204's `q0+q6` loss pricing central `q3`.
6. Add `contact_word`, `lag_barycenter`, `transport_cost`,
   `Fejer_annihilator_projection`, `shell_lag_commutator`,
   `sidecar_entropy`, `scale_survival_bit`, `apex_fold_side`, and
   `cyclotomic_mode` to the next bounded-bank scout.
7. Add HYP-3234's `chart_change_class`, HYP-3233's
   `cyclotomic_factor_signature`, and HYP-3235/HYP-3218 conductor/margin
   labels so the same residual can be compared across signed recursion charts,
   cyclotomic factors, and Fejer field certificates without reusing global
   letter names illegally.
8. Add HYP-3236 Green conductance, HYP-3237 Vitali/Brouwer core-wall, and
   HYP-3219 Brouwer trace-sign/SOS labels so outward lag transport can be
   charged simultaneously to resistance slack, core witness topology, and
   even/odd obstruction factorization.
9. Add HYP-3238 even-positive/odd-negative compression labels and HYP-3239
   dihedral/Borsuk-Ulam labels so the lag transport cannot hide a clipped sign
   irrep, anti-automorphism, or imaginary Gauss-sum obstruction.
10. Add HYP-3241 saddle-index/core-universality labels and HYP-3240
    dilation-witness/non-norm guardrails so core-side transport distinguishes
    AP/Goddyn-Wong primitive witnesses from promoted `Phi_{14d}` dilation
    witnesses and rejected imaginary-quadratic norm scalars.
11. Add HYP-3242 danger-nerve/Euler-characteristic labels so lag transport
    can be compared with cover-hole survival, Cech/Betti sidecars, and the
    Borsuk-Ulam antipodal witness pair rather than only scalar slack.
12. Add HYP-3243 finite chamber/topology-graph carrier labels so every
    accepted lag residual names its open tope, boundary cocircuit, Cech
    component, Toeplitz/Green/root-motion chamber discharge, state-lift
    contradiction, or residual-debt slot.
13. Add HYP-3244 tiling/half-tiling descent labels so every accepted
    out-correlation residual names its tiling witness lift,
    half-tiling descent certificate, path-presentation fiber weight,
    parent-automorphism word orbit, rectangle/hourglass residue, and
    tail/tip deletion signature before scalar quotienting.

-> HYP-3245, HYP-3244, HYP-3243, HYP-3242, HYP-3241, HYP-3240, HYP-3239, HYP-3238, HYP-3237, HYP-3236, HYP-3235, HYP-3234, HYP-3233, HYP-3232, HYP-3230, HYP-3231, HYP-3229, HYP-3228, HYP-3227, HYP-3219, HYP-3218, HYP-3217, HYP-3214, HYP-3226, HYP-3225, HYP-3224, HYP-3223, HYP-3222,
HYP-3213, HYP-3212, HYP-3205, HYP-3204, HYP-3203, HYP-3202, HYP-3163,
HYP-3132, T1309, LTI-309, LTT-209, OPEN-Q-108.

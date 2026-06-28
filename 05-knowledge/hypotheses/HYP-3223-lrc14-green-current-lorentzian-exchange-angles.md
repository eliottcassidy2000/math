---
id: HYP-3223
title: Green-current certificates and Lorentzian exchange circuits are two new LRC14 k=8 proof angles
status: SYNTHESIS / proof-angle proposal; not an LRC14 proof
source: codex-2026-06-28
tangent: T1323
technique: LTI-323
tournament_technique: LTT-223
reflection: 07-reflections/lrc14-green-current-lorentzian-exchange-angles-codex-20260628.md
related:
  - HYP-3236
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3231
  - HYP-3230
  - HYP-3229
  - HYP-3228
  - HYP-3227
  - HYP-3226
  - HYP-3225
  - HYP-3224
  - HYP-3222
  - HYP-3221
  - HYP-3218
  - HYP-3217
  - HYP-3214
  - HYP-3216
  - HYP-3213
  - HYP-3212
  - HYP-3211
  - HYP-3210
  - HYP-3204
  - HYP-3203
  - HYP-3202
  - HYP-3201
  - HYP-3200
  - HYP-3163
  - HYP-3162
  - HYP-3161
  - HYP-3160
  - HYP-3154
  - HYP-3153
  - HYP-3139
  - HYP-3138
  - HYP-3132
  - OPEN-Q-108
---

# HYP-3223: Green-Current And Lorentzian Exchange Angles

## Claim

Two proof-certificate angles should be tried against the remaining LRC14 k=8
targets:

```text
Angle A: Green-current / effective-resistance certificate
Angle B: Lorentzian / valuated-exchange circuit certificate
```

Both are downstream of the exact evidence already on mainline.  HYP-3200 and
HYP-3202 say consecutive speeds are the bounded-bank maximizer for total
empty-sector covariance and for each cyclic-distance covariance layer.  HYP-3202
also says local exchange moves almost work but leave a finite trap manifold.
HYP-3203 says raw cyclotomic closeness and naive compression are the wrong
quotients; the useful root-locus data is AP-polarized support plus orbit
sidecars.  HYP-3210, the Toeplitz/Szego thread, and the Perron/Joukowski
thread say the same object has a spectral face.

Push-time rebase added three relevant guardrails.  HYP-3211 separates the
apex-7 multiplicative tournament face from the additive cyclotomic LRC face.
HYP-3212 makes the cap a Chebyshev equioscillation / Delsarte magic-function
target.  HYP-3221 says config-blind algebraic certificates hit the apex-7
Lee-Yang obstruction and the final rescue must be config-specific analytic
equidistribution.

Rebase also added HYP-3204 as the ordered-tail exchange scout: the current
best coefficient route prices central `q3` mass against bimodality
`q0+q6`.  HYP-3223 is complementary: it asks whether the covariance/trap side
has a Green-current or valuated-exchange certificate that can feed the same
ordered-tail exchange packet.

Namespace note: this packet began as a local HYP-3204 draft, moved to
HYP-3222 after HYP-3204 landed on mainline, and then moved again to HYP-3223
after mainline claimed HYP-3222 for the Joukowski/Hermite-Biehler/Perron exact
certificate packet.  The ideas here are downstream certificate languages, not
the exact Joukowski/HB/Perron certificate itself.

The later HYP-3205 spectral-dictionary packet gives the compatibility target:
AP/consec and doubled AP are the only simultaneous maximizers across the
current `q0`, `L_y`, covariance-layer, AP-support, and Toeplitz coordinates.
HYP-3223 should therefore be tested as a new dictionary-coordinate proposal,
not as another free scalar.  Its Green-current and Lorentzian rows should
predict the same first-failed-coordinate trap discharge that HYP-3205 asks for.

HYP-3223 therefore asks for certificates that could consume the
even/covariance and finite-trap evidence without requiring a globally monotone
speed path.  It is not proposed as a replacement for the HYP-3221 analytic
equidistribution obligation; it is a way to classify and discharge the bounded
k=8 certificate debt before that final analytic step.

## Angle A: Green-Current Certificate

Treat the empty-sector covariance matrix as a Green kernel rather than only as
a PSD or Perron matrix.  The tentative dictionary is:

```text
sector covariance matrix C_E       = Green kernel / response matrix
Sigma kappa_2 = 1^T C_E 1          = all-ones current energy
distance layers D1,D2,D3           = boundary conductance channels
exchange move E -> E'              = local star-mesh / Schur-complement edit
HYP-3202 traps                     = finite bottleneck networks
AP/consec                          = extremal capacity / minimum resistance row
```

The proof target becomes a Thomson/Rayleigh principle:

```text
the AP network minimizes the effective resistance, or maximizes the discharged
all-ones current capacity, among primitive legal k=8 packets; every non-AP
packet either has a Rayleigh-improving network edit or a named bottleneck debt.
```

This is sharper than "covariance is positive" because it records where the
current flows.  Positive-pair rows that beat plain FKG should fail by one of:

```text
boundary_leakage_conductance too large
distance-layer current misaligned
Schur complement has negative trap debt
all-ones current is not the Perron/ground current
```

The immediate scout is not a new bank census; it is a certificate emitter over
the known HYP-3202 rows:

```text
green_kernel_id
effective_resistance_profile
thomson_energy_slack
rayleigh_monotonicity_status
schur_complement_exit
star_mesh_reduction_word
boundary_leakage_conductance
trap_network_bottleneck_id
all_ones_current_alignment
dilation_orbit_gauge_status
terminal_network_discharge_or_debt
```

Wild but testable guess: the `11` arbitrary-swap traps from HYP-3202 are not
random local maxima.  They may be exactly the rows whose Schur-complement
reduction preserves all three distance layers until a boundary leakage term is
exposed.  In that reading, the trap discharge is a finite electrical
network-reduction lemma rather than a finite list of rational inequalities.

## Post-HYP-3236 Execution Note

HYP-3236 executed this Green-current route on the full anchored bounded k=8
bank.  The positive-part covariance conductance graph is AP-tight:
consecutive and doubled AP are the only all-bank maximizers of algebraic
connectivity `lambda2` and total positive conductance, and the only minimizers
of Kirchhoff, mean, max, and distance-layer effective resistance.  The `11`
non-AP arbitrary-exchange traps all have Green resistance excess, split
primarily into `3` Kirchhoff-excess traps and `8` max-resistance bottleneck
traps.

The execution also sharpens the guardrail.  The map
`C_E -> max(C_E,0)` is a lossy compression of the signed covariance kernel.
The actual proof packet is therefore the conductance graph plus negative
covariance leakage, Schur boundary terminals, and odd Worpitzky/Hermite-Biehler
sidecars.  Raw positive association or raw `lambda2` alone is not the proof
object.

Concurrent HYP-3225 executed the trap-local half of the same proposal: the
`11` non-AP exchange traps all discharge first through Toeplitz `lambda_min`,
then split into Green, Rayleigh, and rank-2 Plucker sidecar classes.  Thus
HYP-3223's "next scout" has become two complementary packets: HYP-3236 for
all-bank conductance/connectivity extremality, and HYP-3225 for typed finite
trap fingerprints.

Incoming HYP-3214 adds the harmonic dual face: the 7-sector magic function is
the positive-definite Fejer kernel `F_7=(de Moivre cubic)^2`, equal to AP
autocorrelation, while the cap is the separate 14-clock Johnson pair-Pascal
functional.  This makes the Green-current proposal sharper: try to derive the
conductance Dirichlet form as a finite shadow of the Fejer/autocorrelation
certificate, but do not identify the Green scalar, Fejer kernel, and cap until
the lost coordinates have sidecars.

## Angle B: Lorentzian / Valuated Exchange Certificate

Treat co-emptiness probabilities as a set function:

```text
F_E(A) = P(all sectors in A are empty)
```

Then pair covariance, triple associator debt, AP support direction, and
compression traps are all finite differences of the same object:

```text
Rayleigh difference:      F(ij) - F(i)F(j)
third cumulant/cocycle:   Delta_i Delta_j Delta_k log F
AP support inequality:    exposed normal to the miss-vector cone
exchange trap:            local failure of a basis-exchange inequality
```

The proposed proof target is not "the miss polynomial is Lorentzian" as a
blind slogan.  The target is more local:

```text
after keeping dilation/mirror/two-block sidecars, the relevant k=8 covariance
functional lies in a Lorentzian or valuated-matroid exchange chamber whose
normal cone exposes the AP row.
```

If true, HYP-3202's exchange traps are tropical Plucker circuit events: they
are rows where a naive exchange inequality is tied or points the wrong way
because the proof has forgotten the sidecar that selects the chamber.

The next scout should emit:

```text
rayleigh_difference_matrix
lorentzian_hessian_signature
valuated_exchange_slack
tropical_plucker_defect
matroid_circuit_sidecar_id
basis_exchange_trap_class
strong_rayleigh_status
m_concavity_status
ap_chamber_normal_vector
sidecar_restored_exchange_status
terminal_exchange_discharge_or_debt
```

This angle is compatible with HYP-3203.  The AP-polarized support inequality
would become a normal-fan statement; orbit-aware compression would become the
sidecar set needed to stay inside the correct chamber.

## Combined Proof Sketch

The two angles may be the same certificate in dual languages:

```text
Green-current side:
  minimizing resistance / maximizing all-ones current capacity

Lorentzian side:
  normal-fan exchange chamber / valuated basis exchange

bridge:
  Schur complements of response matrices are valuated exchange moves
```

The useful remaining proof target is therefore not one scalar and not one
local monotone path.  It is a pair of certificate ledgers:

```text
for every primitive k=8 row E:
  either electrical edit improves toward AP,
  or valuated exchange emits the missing sidecar,
  or the row is one of a finite trap networks discharged directly.
```

Rebase guardrail: if this ledger becomes config-blind, HYP-3221 predicts it
will fail at the same algebraic obstruction as Bonferroni/Asano/SOS.  The
certificate must keep the row's arithmetic sidecars or explicitly hand them to
the equidistribution route.

## Tournament Analysis

Vertices are proof carriers, not runners, sectors, or tournament arcs:

```text
green_current_certificate
effective_resistance_profile
rayleigh_schur_trap_discharge
lorentzian_hessian_certificate
valuated_exchange_certificate
tropical_plucker_defect
finite_trap_sidecar
polarized_cyclotomic_support
raw_covariance_scalar
plain_speed_exchange
plain_positive_association
```

Pairwise observable: which carrier preserves the LRC14 `Sigma kappa_2` /
coverage-cap predicate while naming destroyed coordinates and terminal debt.

Switch/gauge: orient from a carrier that needs an unproved monotone path toward
one that provides a certificate, sidecar, or finite discharge.  Ties are broken
by whether the carrier also preserves HYP-3203's AP support direction.

Fingerprint:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1}
directed_3cycles = 0
hamiltonian_path_count = 1
selected path:
green_current_certificate
-> effective_resistance_profile
-> rayleigh_schur_trap_discharge
-> lorentzian_hessian_certificate
-> valuated_exchange_certificate
-> tropical_plucker_defect
-> finite_trap_sidecar
-> polarized_cyclotomic_support
-> raw_covariance_scalar
-> plain_speed_exchange
-> plain_positive_association
```

## Assumption Challenge

Considered vertex sets: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
Toeplitz moments, OPUC coefficients, exchange moves, trap rows, network nodes,
matroid circuits, and proof obligations.

For this session the chosen vertices are proof carriers.  Runner or arc
tournaments destroy the certificate payload; raw row statistics destroy where
current leaks and which exchange sidecar is missing.  The preserved predicate
is primitive k=8 covariance/coverage extremality plus the AP support direction.
The destroyed information that must be named is dilation orbit, mirror orbit,
two-block trap status, boundary leakage, and odd Worpitzky/third-cumulant debt.

## Next Tests

1. For the `12` arbitrary-swap local maxima from HYP-3202, compute
   effective-resistance profiles and Schur-complement reduction words.  Check
   whether the `11` decoys share a small set of bottleneck types.
2. For the same rows, compute Rayleigh differences, Hessian signatures of the
   relevant homogenized co-emptiness polynomial, and tropical Plucker defects.
   Check whether sidecar-restored exchange removes the local maxima.
3. Test whether the HYP-3203 AP support normal is also a normal vector of the
   valuated exchange chamber, whether the Green-current certificate produces
   the same exposed face, and whether HYP-3205's dictionary coordinate that
   first fails on a trap predicts the same discharge.
4. Keep the odd Worpitzky/associator channel as named debt.  Do not force it
   into the electrical or Lorentzian certificate until the third-cumulant
   sidecar has its own discharge.

## Execution Addendum: HYP-3227

HYP-3225 executes the local trap-fingerprint audit promised above, and
HYP-3227 executes the full-bank conductance-graph extension on the exact
HYP-3205/HYP-3224 bounded k=8 bank.  HYP-3227 confirms that several
Green-current coordinates are AP-tight sidecars:
all-ones Green energy, positive-covariance sector graph lambda2/Kirchhoff,
precision graph lambda2/Kirchhoff, and precision killing have no primitive
beaters beyond the doubled AP equality row.  It also names the failure mode:
precision M-matrix defect has `181` primitive beaters, so inverse-Green
conductance is not a terminal scalar.

The main upgrade is the trap-discharge graph.  The `11` non-AP exchange traps
remain connected to certificate coordinates after Toeplitz is deleted
(`lambda2=2.537866286`) and even using only Green-current coordinates
(`lambda2=1.208613477`).  Therefore the next electrical target is a finite
Schur-complement/star-mesh lemma for the trap graph, not another scalar
conductance ranking.  Angle B should now be tested against the same trap
graph: does the Lorentzian/valuated-exchange defect see the same
Fiedler-positive M-matrix-defect island?

## Links

HYP-3223 is a proof-angle synthesis over HYP-3222, HYP-3221, HYP-3213, HYP-3212,
HYP-3211, HYP-3210, HYP-3205, HYP-3204, HYP-3203, HYP-3202, HYP-3201, HYP-3200,
HYP-3163, HYP-3162, HYP-3161, HYP-3160, HYP-3154, HYP-3153, HYP-3139,
HYP-3138, HYP-3132, HYP-3227, HYP-3226, HYP-3225, THM-577, LTI-323, LTI-325, LTT-223, LTT-225,
T1323, T1325, and OPEN-Q-108.

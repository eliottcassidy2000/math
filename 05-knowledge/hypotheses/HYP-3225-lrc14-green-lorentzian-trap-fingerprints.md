---
id: HYP-3225
title: LRC14 Green-current / Lorentzian trap fingerprints
status: EVIDENCE / exact finite trap-neighborhood scout; not an LRC14 proof
source: codex-2026-06-28
tangent: T1308
technique: LTI-308
tournament_technique: LTT-208
script: 04-computation/lrc14_green_lorentzian_trap_fingerprints_codex_20260628.py
result: 05-knowledge/results/lrc14_green_lorentzian_trap_fingerprints_codex_20260628.out
reflection: 07-reflections/lrc14-green-lorentzian-trap-fingerprints-codex-20260628.md
related:
  - HYP-3227
  - HYP-3224
  - HYP-3223
  - HYP-3222
  - HYP-3221
  - HYP-3213
  - HYP-3212
  - HYP-3211
  - HYP-3210
  - HYP-3205
  - HYP-3204
  - HYP-3203
  - HYP-3202
  - HYP-3201
  - HYP-3200
  - HYP-3163
  - HYP-3132
  - HYP-3143
  - HYP-3148
  - HYP-3149
  - HYP-3199
  - HYP-3226
  - OPEN-Q-108
---

# HYP-3225: LRC14 Green-Current / Lorentzian Trap Fingerprints

## Claim

HYP-3225 executes the HYP-3223 next hook on the HYP-3202 arbitrary-swap
local maxima.  The main conclusion is two-layered:

```text
Toeplitz lambda-min is still the universal first trap discharge.
The residual trap mechanisms split into Green/Rayleigh/Plucker sidecar types.
```

Thus HYP-3224's moment-cone chart is stronger, not weaker, after the
Green-current and Lorentzian audit.  Toeplitz is not just one scalar
explanation.  It behaves like a boundary chart that absorbs several
nonisomorphic local obstruction mechanisms.

## Exact Scout

The scout imports the `12` HYP-3202 arbitrary-swap local maxima, including
consecutive, and evaluates each one together with all one-swap neighbors:

```text
bank = anchored bounded k=8, E={0} union A, A subset [1,14], |A|=7
rows_evaluated = 577
local_maxima = 12
target = (0,1,2,3,4,5,6,7)
```

The consecutive baseline is:

```text
toeplitz_lambda_min = +0.042304730706
Sigma_kappa2 = 6237419/8643600
positive_laplacian_lambda2 = +0.192033074001
effective_resistance_by_distance = (6.710086, 7.345148, 7.675710)
positive_cut_min = +0.176504234347
hessian_signature = (pos, zero, neg) = (6, 0, 0)
conditional_rayleigh_negatives = 56
pair_tropical_plucker_max_gap = +0.458959814453
```

The AP row already has conditional Rayleigh negatives and a positive pair
Plucker gap.  Therefore those raw signs are not themselves terminal proof
obstructions.  They become useful only as relative trap fingerprints after the
Toeplitz discharge chart and orbit sidecars have been retained.

## Classification

All `12` rows select `Toeplitz_lambda_min` as the first dictionary discharge
coordinate in the HYP-3224/HYP-3205 sense.  The `11` non-AP traps then split:

```text
AP_terminal: 1
rank2_pair_plucker_bottleneck: 6
green_low_connectivity_bottleneck: 2
mixed_green_lorentzian_sidecar: 1
AFM_frustrated_high_rayleigh_debt: 2
```

The scout's correlations also show why a one-signal replacement is unlikely:

```text
corr(lambda2_ratio, toeplitz_deficit) = -0.636805
corr(plucker_gap, toeplitz_deficit) = +0.285960
corr(rayleigh_negative_count, toeplitz_deficit) = -0.247417
```

The moment-cone deficit is not identical to any one measured electrical or
Lorentzian scalar.  Its value is that it sees all of the trap classes at once.

## Map A: Certificate Geometry

The current proof geometry should be read as a multi-chart atlas:

```text
exchange/covariance chart:
  proves improvement away from the finite trap manifold

Toeplitz moment-cone chart:
  discharges every measured trap boundary

Green-current chart:
  explains low-connectivity and high-resistance bottlenecks

Lorentzian / valuated-exchange chart:
  explains rank-2 pair Plucker bottlenecks and tied exchange circuits

Rayleigh / AFM chart:
  identifies frustrated conditional-dependence debt

Hermite-Biehler / Joukowski chart:
  glues the even covariance side to the odd Worpitzky side
```

This strengthens the HYP-3224 route.  A proof can aim for:

```text
bulk lemma:
  outside the trap manifold, AP support / covariance layers / ordered tails
  improve by a legal exchange or compression.

finite boundary theorem:
  on the trap manifold, Toeplitz curvature is strictly deficient unless the
  row is AP or the known nonprimitive doubled AP equality artifact.

sidecar theorem:
  each Toeplitz-deficient trap has a typed Green, Rayleigh, or pair-Plucker
  certificate explaining why local exchange alone stalled.
```

The ambitious version is that the sidecars are the visible finite sections of
one Chebyshev/Cohn-Elkies/Delsarte magic-function dual over the
`Q(cos(2pi/7))` face from HYP-3212/HYP-3213.

## Map B: Proof-Circuit / Filler-Core Reading

The Erdos-870-inspired tournament/filler work suggests a second map.  A hard
finite theorem can be made tractable by isolating:

```text
live core       = HYP-3202 exchange traps
deterministic filler = Toeplitz moment-cone boundary chart
canary fields   = Green bottleneck, Rayleigh debt, Plucker circuit
terminal gate   = AP/doubled-AP equality decision
```

In circuit-complexity language, HYP-3225 proposes a small proof dispatcher
rather than a monolithic scalar classifier:

```text
input row
  -> exchange-gradient gate
  -> dictionary-discharge gate
  -> Toeplitz boundary gate
  -> Green/Rayleigh/Plucker sidecar gate
  -> HB/Joukowski odd-side gluing gate
  -> analytic/equidistribution handoff
```

The proof circuit stays below the bounded-core degree-four wall noted in the
HYP-3132/HYP-3150/HYP-3151 thread: no generic quintic scalarization is needed
if each gate carries its destroyed coordinate as a named sidecar.

Post-pull integration: HYP-3226 reserves the next small-pattern atlas.  Its
motif table should treat the five HYP-3225 trap classes as typed payload
atoms, not as free-standing numerology.  In particular, `11/12` trap counts
are useful only when attached to Toeplitz discharge plus Green/Rayleigh/
Plucker repair sidecars.

Post-rebase integration: HYP-3227 extends this local trap-neighborhood audit
to the full anchored bounded k=8 bank and builds trap/certificate conductance
graphs.  The two packets should be read together: HYP-3225 names the residual
trap types, while HYP-3227 shows the trap/certificate graph remains connected
without Toeplitz and even with Green-only coordinates.

## Tournament Analysis

The tournament vertices are proof certificates and trap sidecars, not runners
or arcs.  The pairwise observable is: which certificate explains trap
discharge while retaining the most normal-fan, moment, and trap-local payload?
The gauge is `A -> B` when `A` keeps more proof-relevant information than `B`
or makes `B` a shadow.

The measured tournament is transitive:

```text
score_hist = {99:1, 91:1, 84:1, 78:1, 72:1, 66:1, 42:1, 30:1, 9:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
priority_path =
  toeplitz_trap_discharge ->
  green_effective_resistance_profile ->
  conditional_rayleigh_surplus_census ->
  rank2_tropical_plucker_gap ->
  schur_payload_exit ->
  covariance_hessian_signature ->
  raw_positive_covariance_graph ->
  plain_exchange_local_maximum ->
  runner_or_arc_tournament
```

This is exactly the assumption-challenge result: runner or arc tournaments are
too lossy for this target.  Certificate tournaments preserve the LRC predicate
and the sidecar obligations.

## Assumption Challenge

Alternate vertex sets considered:

```text
runners
arcs
sectors
covariance pairs
positive-conductance graph cuts
conditional Rayleigh events
pair Plucker quadruples
exchange moves
proof-obligation charts
```

Chosen vertices: certificate types and trap-local sidecars.

Preserved predicate: primitive k=8 AP/consecutive covariance and coverage
extremality, together with HYP-3224 normal-fan trap discharge.

Destroyed information: raw speed order, full PGF root curve, and the
Hermite-Biehler odd side unless named sidecars are kept.

Challenged assumption: "Toeplitz discharge should be one scalar explanation."
The scout instead treats Toeplitz as the chart switch that absorbs several
Green/Rayleigh/Plucker trap mechanisms.

## Next Tests

1. Prove the finite trap classification symbolically: derive exact
   inequalities defining the five trap classes above.
2. Test whether the same sidecar families persist for k=9, k=10, and larger
   bounded neighborhoods of the HYP-3202 trap rows.
3. Express `rank2_pair_plucker_bottleneck` as an explicit valuated matroid or
   tropical Plucker relation on co-emptiness weights.
4. Express `green_low_connectivity_bottleneck` as a Schur complement,
   Verblunsky, or Fejer-Riesz slack.
5. Connect the dispatcher to HYP-3204's central exchange-rate inequality:
   determine whether the Toeplitz slack projects linearly to `q0+q6` loss
   pricing `q3` gain.
6. Use the HYP-3222 HB/Perron packet to attach the odd Worpitzky sidecar to
   each finite trap type.

-> HYP-3225, HYP-3224, HYP-3223, HYP-3222, HYP-3221, HYP-3213, HYP-3212,
HYP-3211, HYP-3210, HYP-3205, HYP-3204, HYP-3203, HYP-3202, HYP-3201,
HYP-3200, HYP-3163, HYP-3132, HYP-3143, HYP-3148, HYP-3149, HYP-3199,
HYP-3226, T1308, LTI-308, LTT-208, OPEN-Q-108.

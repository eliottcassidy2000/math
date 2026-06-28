---
id: HYP-3260
title: LRC14 unit equioscillation nullspace
status: EVIDENCE / exact finite scout; not an LRC14 proof
source: codex-2026-06-28
tangent: T1348
technique: LTI-348
tournament_technique: LTT-248
script: 04-computation/lrc14_unit_equioscillation_nullspace_codex_20260628.py
result: 05-knowledge/results/lrc14_unit_equioscillation_nullspace_codex_20260628.out
reflection: 07-reflections/lrc14-unit-equioscillation-nullspace-codex-20260628.md
related:
  - HYP-3300
  - HYP-3265
  - HYP-3259
  - HYP-3258
  - HYP-3257
  - HYP-3256
  - HYP-3255
  - HYP-3254
  - HYP-3253
  - HYP-3252
  - HYP-3251
  - HYP-3250
  - HYP-3249
  - HYP-3248
  - HYP-3247
  - HYP-3246
  - HYP-3245
  - HYP-3243
  - HYP-3242
  - HYP-3241
  - HYP-3132
  - HYP-2909
  - THM-523
  - OPEN-Q-108
---

# HYP-3260: LRC14 Unit Equioscillation Nullspace

## Claim

The HYP-3246/HYP-3247 Chebyshev frame is real but incomplete as a terminal
proof object.  The six unit active gradients at `a/14`, `a in (Z/14)*`, have
rank exactly `3` over the residue coordinates `1..13`:

```text
a= 1:  +e_1  - e_13
a=13:  -e_1  + e_13
a= 3:  +e_5  - e_9
a=11:  -e_5  + e_9
a= 5:  +e_3  - e_11
a= 9:  -e_3  + e_11
```

Thus the local unit-index coordinate is the three antipodal binding-pair
coordinate.  Its zero columns are

```text
2, 4, 6, 7, 8, 10, 12
```

and the covering residue `0 mod 14` is also invisible until it kills all unit
points.  Therefore a proof using the index/equioscillation frame must retain a
blind residue/height sidecar plus safe-component topology.

## Exact Readout

The script audits named rows and one-swap AP collars up to added speed `84`.

Key named contrast:

```text
AP                 mass=0,      unit projection unchanged
GW 12->24          mass=0,      unit projection unchanged
near 12->36        mass=1/1260, unit projection unchanged
petal 10->20       mass=1/980,  unit projection unchanged
blind decoy 2->16  mass=11/364, unit projection unchanged and mod-14 delta 0
```

So the unit projection cannot distinguish AP from GW, the first positive
near-miss `12->36`, or the same-residue height move `2->16`.  The mod-14
nonunit ledger distinguishes `12->24` from `12->36`, but it still misses
`2->16`; height/scale data are not optional.

One-swap collar scan:

```text
rows=923
unit-blind rows=317
unit-blind boundary rows=1          # 12->24
unit-blind positive rows=316
smallest unit-blind positive row=12->36, mass=1/1260
```

The covering branch is separately visible because a multiple of `14` changes
the unit status from `E6 K0 S0` to `E0 K6 S0`: all six unit points are killed.

## Proof Consequence

HYP-3246/HYP-3247 should be used as an index-coordinate theorem, not as a
complete classifier.  The exact gluing target is:

```text
unit rank-3 index
  + blind residue/height ledger
  + strict-safe component atlas
  + covering 14-multiple kill switch
```

This connects the new Chebyshev/index frame back to the older sector and
Q108 machinery: the local unit object proves the 14-free witness fragment
and names the tight core, while the blind directions are exactly where
Goddyn-Wong, near-miss positive rows, height moves, and covering-floor rows
live.

Fetch-time integration: incoming HYP-3250, HYP-3251/HYP-3252, HYP-3253/S81,
HYP-3254's Q(sqrt-7) floor-SPEC reorganization, HYP-3255/S82, HYP-3256,
HYP-3258, HYP-3259, and HYP-3265 sharpen this from both sides. HYP-3250's S80 tests
say the tight locus behaves like AP/GW dilations plus a uniform-margin
complement, with witness peaks `=(p-1)*d`.  HYP-3251/HYP-3252 say the index is
ambient/descriptive rather than S-dependent proof content.  HYP-3253/S81 makes
the construction/unit witnesses rigorous and the bounded single-swap margin
rigorous, while isolating resonant survivor-positivity as the analytic core.
HYP-3254 says the Q(sqrt-7) split partially organizes the floor SPEC but does
not replace the S-dependent floor proof.  HYP-3255 reframes resonant danger as
an on-grid core/off-grid bulk split, and HYP-3256 says the apex-7 residue layer
organizes but does not decide the magnitude-level tight-vs-loose census.
HYP-3258 gives the census version of the same split: unit runners form the
binding skeleton, nonunit runners form the covering/magnitude layer, and only
`12 -> 24` remains flexible.
Incoming HYP-3257/S83 adds the C3 orbit view of the three binding pairs; this
is compatible with the rank-3 readout because the rank-3 unit coordinate is
the three-point C3 witness orbit, still separated from the blind covering layer.
HYP-3259 gives the real-manifold version: unit binding speeds are
infinitesimally rigid, while covering speeds supply the flex directions.
HYP-3265 gives the contact-graph case-split version: the unit contact
matching is the equality-branch classifier, and this nullspace audit is the
warning that the classifier must retain blind residue/height and covering
sidecars.
HYP-3300 then promotes these sidecars into observability columns and Morse
boundary data for a finite-chamber descent proof attempt.
That is exactly the rank-3 nullspace warning: the unit index describes the
saddle target, while the blind residue/height ledger and safe-component/floor
branches carry the S-dependent work.

## Tournament Analysis

Vertices are proof carriers rather than runners:

```text
strict_safe_component_atlas
unit_plus_blind_residue_height_ledger
covering_14_multiple_kill_switch
unit_active_gradient_rank
full_unit_residue_signature
raw_unit_values
```

Pairwise observable: which carrier decides the threshold branch while
retaining the unit index and blind sidecars.

Switch/gauge: `A -> B` when `A` wins a majority of branch-decision,
unit-index, blind-residue/height, covering, endpoint, and formalization
criteria.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1]
hamiltonian_path_count=1
priority_path=strict_safe_component_atlas
  -> unit_plus_blind_residue_height_ledger
  -> covering_14_multiple_kill_switch
  -> unit_active_gradient_rank
  -> full_unit_residue_signature
  -> raw_unit_values
```

## Assumption Challenge

Alternate vertices considered: runners, unit points, complement-pair binders,
unit-gradient rows, nonunit residues, height/scale moves, cover arcs, safe
components, and proof obligations.

Chosen vertices: proof carriers plus the unit-gradient/nullspace split.

Preserved predicate: the HYP-3246/HYP-3247 unit threshold certificate and the
exact strict-safe/open-boundary distinction at `delta=1/14`.

Destroyed information: nonunit residue placement, height/scale placement,
higher safe peaks, endpoint-owner topology, and covering-floor behavior unless
the blind sidecar is retained.

Challenged assumption: equioscillation at the six unit points is already a
complete local proof object.  It is a rank-3 index coordinate and must be
glued to blind residue/height ledgers plus safe-component topology.

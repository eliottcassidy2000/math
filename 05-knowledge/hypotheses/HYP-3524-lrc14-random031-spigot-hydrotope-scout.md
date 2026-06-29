---
id: HYP-3524
title: LRC14 random031 spigot/hydrotope residual scout
status: EVIDENCE / spigot-emitter and sliced-box chamber scout; not an LRC14 proof
source: codex-2026-06-29 inspired by the spigot algorithm page and arXiv:2606.28280, integrating HYP-3522, HYP-3521, HYP-3520, HYP-3513, HYP-3512, HYP-3494, HYP-3511, HYP-3510, HYP-3493, HYP-3490, and HYP-3486
tangent: T1524
technique: LTI-524
tournament_technique: LTT-424
script: 04-computation/lrc14_random031_spigot_hydrotope_scout_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_spigot_hydrotope_scout_codex_20260629.out
reflection: 07-reflections/lrc14-random031-spigot-hydrotope-scout-codex-20260629.md
related:
  - HYP-3522
  - HYP-3521
  - HYP-3520
  - HYP-3513
  - HYP-3512
  - HYP-3494
  - HYP-3511
  - HYP-3510
  - HYP-3493
  - HYP-3490
  - HYP-3486
  - HYP-3485
  - HYP-3484
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - OPEN-Q-108
---

# HYP-3524: LRC14 Random031 Spigot/Hydrotope Residual Scout

## Claim

The spigot-algorithm analogy is useful if it is read as an online proof
emitter, not as a numerical metaphor.  In `random_covering_031`, owner labels
should be emitted only when a local certificate makes them stable, and the
remaining tail should never be allowed to re-enter through a scalar quotient.

HYP-3522 already gives the exact finite data:

```text
seam owners                    = (23,45,93,113,147,169,173)
transport owners               = (23,93,113)
branch-boundary owners         = (23,93,147,169)
bracket lift                   = (147,169)
transport plus branch boundary = (23,93,113,147,169)
residual tail                  = (45,173)
```

HYP-3524 packages this as a spigot schedule:

```text
S0 forbidden seam input:
  tail = (23,45,93,113,147,169,173)

S1 transport emitter:
  emit (23,93,113)
  tail = (45,147,169,173)

S2 branch-boundary lift emitter:
  emit (147,169)
  tail = (45,173)

S3 residual tail:
  prove the two-owner puncture/apex boundary lemma
```

The computed safety checks are:

```text
cumulative_sizes = (0,3,5,5)
tail_sizes       = (7,4,2,2)
monotone_cumulative = True
monotone_tail       = True
no_duplicate_emit   = True
unemitted_tail      = (45,173)
```

Thus the residual target from HYP-3522 can be stated as a no-hidden-tail
theorem: once transport and bracket owners have been emitted, no later legal
quotient may smuggle back the discarded owners or confuse the residual pair
with bracket/transport shadows.

## External Prompt Integration

The spigot source supplies the proof discipline:

```text
head is processed and discarded only modulo a certificate;
tail remains bounded and explicit;
digits are emitted left-to-right and not reused.
```

The hydrotope paper supplies a second diagnostic: chamber sign patterns for a
box sliced by a hyperplane.  In HYP-3524, the sliced box is not the theorem.
It is a canary for whether a threshold/chamber quotient has kept enough owner
information.

The arXiv paper is:

```text
Surface Water Wave Scattering and the Hydrotope
arXiv:2606.28280
```

The transferable part is the organization of sign-pattern chambers by a sliced
box / hyperplane volume, not the water-wave physics.

## Exact Readout

Executable scout:

```text
04-computation/lrc14_random031_spigot_hydrotope_scout_codex_20260629.py
05-knowledge/results/lrc14_random031_spigot_hydrotope_scout_codex_20260629.out
```

Branch-boundary and mirror cross-check:

```text
branch=0 bypass_u=(679,680,681,682,683,684)
  phases=(7,8,9,10,11,12)
  left owners=(93,147,169)
  right owners=(23,169)
  intersection=(169,)

branch=1 bypass_u=(527,528,529,530,531,532)
  phases=(2,3,4,5,6,7)
  left owners=(23,169)
  right owners=(93,147,169)
  intersection=(169,)

mirror_pair_owner_words={(23,93,113):6}
```

This is exactly the HYP-3522 bracket lift: owner `169` is the common bracket
hinge on both sheets, while `147` appears in the ordinary boundary union and
not in the pure bypass transport word.

## Hydrotope Chamber Audit

For each weight system, the script forms all `2^7` owner subsets and compares
their chamber signs against three thresholds:

```text
transport
transport_plus_boundary
residual
```

The surprise is that the natural residue shadows do not separate the proof
targets:

```text
residue_mod14:
  residual bucket size = 3
  residual examples = (45,173), (113,147), (147,169)

centered_residue:
  residual bucket size = 5
  residual examples include (23,45), (45,93), (45,173), (113,147), (147,169)

filtration_layer:
  residual bucket size = 15
```

But the owner-support-cell weight system is singleton for the live targets:

```text
owner_support_cells weights =
  {23:82,45:60,93:112,113:96,147:122,169:102,173:114}

transport bucket size               = 1
transport_plus_boundary bucket size = 1
residual bucket size                = 1
```

This is not a proof by volume.  It says the chamber sign quotient becomes
legal only after it has secretly retained enough owner-support geometry to
separate the residual pair.  Residue-only and layer-only chamber quotients are
canaries for illegal forgetting.

## Sliced-Box Volumes

The scout also computes normalized inclusion-exclusion volumes for

```text
{0 <= x_i <= w_i, sum x_i <= h}
```

at the transport, boundary, and residual thresholds.  These are useful
diagnostics but not legal terminal quotients.  They deliberately forget owner
identity.

Example readout:

```text
owner_support_cells:
  V_transport = 36194029738549/147743976360960
  V_boundary  = 1459908749489263/1477439763609600
  V_residual  = 17531014120337/1477439763609600
  tail_after_boundary = 17531014120337/1477439763609600
```

The equality of `V_residual` and `tail_after_boundary` here is a chamber
diagnostic in this chosen coordinate system, not an allowed replacement for
the residual owner lemma.

## Quotient Safety

Clean or conditionally clean:

```text
full_filtration_word:
  safe=True
  reconstructs transport, bracket_lift, residual

transport_plus_boundary_word:
  safe=True only with the HYP-3522 reconstruction split

hydrotope_signature_owner_support_cells:
  safe=True as a labelled chamber diagnostic

route_sidecar_R:
  safe=True and required by HYP-3513 unless route reconstruction is proved
```

Unsafe:

```text
bypass_owner_word_only:
  cannot see bracket lift or residual

hydrotope_signature_residue_mod14:
  mixed buckets residual=3, boundary=3

hydrotope_signature_centered_residue:
  mixed buckets residual=5, boundary=5

hydrotope_signature_filtration_layer:
  mixed buckets residual=15, boundary=15

sliced_box_volume_*:
  scalar volume forgets owner identities

raw_counts_7_5_2:
  sees the tail sizes but forgets labels
```

## Proof Pull

The next proof-facing theorem should be an online emitter theorem:

```text
emit transport owners (23,93,113)
then emit bracket-lift owners (147,169)
then prove no hidden tail remains except residual pair (45,173)
with route sidecar R retained
```

This suggests a Lean interface:

```text
EmitterState :=
  seam_word
  emitted_word
  tail_word
  certificate
  route_sidecar

emit_transport : State seven -> State four
emit_bracket   : State four  -> State two
close_tail     : State two   -> terminal_discharge
```

HYP-3524 does not close the residual lemma.  It makes the residual lemma more
formal: prove that `(45,173)` is the only legal terminal tail after the two
emitters, and that no quotient used afterward identifies `(45,173)` with
`(147,169)`, `(113,147)`, or any other chamber-mate seen in residue shadows.

## Tournament Analysis

Vertices are proof emitters and quotient sidecars, not runners or raw arcs.

Pairwise observable:

```text
tail shrinkage
owner-label reconstruction
route-sidecar legality
scalar-forgetting penalty
```

Switch/gauge: higher retained proof payload; ties use emitter order.

Fingerprint:

```text
score_hist={100:1,94:1,89:1,83:1,78:1,39:1,31:1,13:1}
directed_3cycles=0
hamiltonian_path=
  full_filtration_spigot_packet
  -> residual_pair_tail_lemma
  -> transport_plus_boundary_emitter
  -> hydrotope_chamber_audit_with_owner_labels
  -> route_sidecar_R_guard
  -> sliced_box_volume_shadow
  -> raw_threshold_sign_shadow
  -> raw_owner_count_shadow
```

## Assumption Challenge

Candidate vertices considered: owners, emitted digits, branch-boundary events,
hydrotope chambers, sliced-box thresholds, `u`-fibers, mirror pairs, route
sidecars, and proof obligations.

Chosen vertices are proof emitters and quotient sidecars.  This preserves the
random031 residual predicate only when emitted owner digits and route sidecar
`R` can be reconstructed.  Scalar chamber volumes deliberately destroy that
information and remain diagnostics.

---
id: HYP-3484
title: LRC14 random031 forbidden-seam flow geometry
status: EVIDENCE / geometric-topological reframing; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3455, HYP-3460, HYP-3477, HYP-3478, HYP-3479, HYP-3480, HYP-3481, HYP-3482, and HYP-3483
tangent: T1444
technique: LTI-444
tournament_technique: LTT-344
script: 04-computation/lrc14_random031_forbidden_seam_flow_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_forbidden_seam_flow_codex_20260629.out
reflection: 07-reflections/lrc14-random031-forbidden-seam-flow-codex-20260629.md
related:
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3480
  - HYP-3479
  - HYP-3478
  - HYP-3477
  - HYP-3476
  - HYP-3475
  - HYP-3472
  - HYP-3460
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3450
  - HYP-3438
  - HYP-2241
  - THM-523
  - OPEN-Q-108
---

# HYP-3484: LRC14 Random031 Forbidden-Seam Flow Geometry

## Claim

HYP-3484 is the seam-surgery refinement of the HYP-3482 punctured-cylinder
atlas and the HYP-3483 recursion-flow comparator.  The hard pair in
`random_covering_031` should be treated as a forbidden seam, not a wall that
the phase flow must cross.  The row is a mirror-punctured cylinder:

- four isolated dead components form two mirror-paired punctures;
- the max-delta mirror gate pair on components `43` and `54` carries the
  seven-owner gluing charge `(23,45,93,113,147,169,173)`;
- all q=`14V` phase witnesses avoid that max-delta seam;
- the same components are still touched exactly `12` times through a
  lower-delta mirror bypass.

This reframes the named random031 clause as:

```text
mirror-punctured cylinder
+ seven-owner forbidden seam
+ lower-delta phase bypass
+ C=27 carry / doubling split
=> random031 terminal gluing packet
```

## Exact Readout

Executable scout:

```text
04-computation/lrc14_random031_forbidden_seam_flow_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_random031_forbidden_seam_flow_codex_20260629.out
```

Puncture/cylinder data:

```text
components=98
component_class_hist={'both_alive':50,'branch0_only':22,'branch1_only':22,'dead_both':4}
dead_islands=4
dead_mirror_pairs=((12,85),(13,84))
dead_projection_edges=0
projection_components=4
```

Forbidden seam:

```text
max_delta=7
seam_gate_count=2
seam_components=(43,54)
seam_owner_union=(23,45,93,113,147,169,173)
seam_endpoint_speeds=(23,113)
seam_endpoint_residues_mod14=((23,9),(113,1))
c27_duplicate_residues={12:2}
```

Phase-flow complement:

```text
V=173
q=14V=2422
actual_witnesses=282
any_gate_hits=242
no_gate_hits=40
hard_seam_phase_hits=0
lower_delta_same_component_bypass_hits=12
bypass_component_branch_delta_hist={(43,'branch1',2):6,(54,'branch0',2):6}
```

The two bypass gates are lower-delta owner-current gates on the same
components, with owner union `(23,93,113)`, not the seven-owner seam.

HYP-3480 now supplies the zero-edge singleton-current contrast:

```text
random031_terminal_class=random031_hard_currentless_control
control_components_with_complete_branch_unit_touch=0/4
control_mirror_pairs_with_branch_unit_mirror_gate=0/2
```

So the six small-touch rows can plausibly close by mirror-unit singleton
current, but `random_covering_031` cannot.  It remains a separate seam-flow
terminal packet.

## Surgery Tests

The scout performs two local gate surgeries.

Deleting the two max-delta seam gates changes nothing in phase routing:

```text
remove_seam:
  any_gate_hits=242
  no_gate_hits=40
  delta_any=0
  delta_no=0
```

Deleting the two lower-delta bypass gates moves exactly the same-component
traffic from gate-hit to no-gate:

```text
remove_bypass:
  any_gate_hits=230
  no_gate_hits=52
  lost_any=12
  gained_no=12
```

So the hard pair is not a phase wall.  It is a forbidden seam whose complement
carries the phase flow; the lower-delta bypass is the actual phase channel on
those components.

## Recursion Connection

This packet makes the HYP-3483 controlled-span reading surgical: `n+2`
describes the owner-boundary seam debt, while `n*2` describes the two-adic
phase-flow bypass.

For `n=14`, the carry modulus from the triangular bridge is:

```text
C=2n-1=27
8*C(n,2)+1=729=(2n-1)^2
```

The additive `n+2` face sees the punctures and the `C=27` carry seam.  The
multiplicative `n*2` face sees phase flow under `u=2t mod 1`.  Random031 is an
interface defect: the carry face sees a seven-owner max-delta seam, while the
doubling face routes around it and pays only a 12-hit lower-delta bypass toll.

## Experiment Design

Four next tests are now precise:

1.  Seam-deletion invariance: deleting the max-delta seam should remain
    invisible to q=`14V` phase routing.
2.  Bypass-channel sensitivity: deleting lower-delta same-component bypass
    gates should remove exactly the component-local phase traffic unless a
    wider gate basis supplies a secondary route.
3.  Carry-lift stress: perturb or lift seam owners by `C=27` classes and test
    whether preserving endpoint residues `(1,9) mod 14` predicts seam/bypass
    stability better than preserving raw owner size or raw delta.
4.  Puncture filling: compare HYP-3478/HYP-3480 singleton pockets with
    random031's four punctures and test whether isolated dead islands become
    dangerous only once a seven-owner seam connects their surrounding flow
    components.

## Tournament Analysis

Vertices are geometric/topological proof carriers, not runners or arcs.

```text
pairwise_observable =
  puncture isolation + forbidden-seam avoidance + lower-delta bypass
  + C27 carry/doubling split + scalar firewall
score_hist={5:1,58:1,63:1,64:2,67:1,69:1}
directed_3cycles=0
hamiltonian_path =
  F00_forbidden_seam_complement_flow
  -> F01_lower_delta_mirror_bypass_channel
  -> F03_c27_carry_vs_phase_doubling_split
  -> F04_seven_owner_gluing_clause
  -> F02_mirror_punctured_cylinder
  -> F05_small_touch_singleton_pocket_shadow
  -> F06_raw_wall_or_delta_scalar
```

Assumption challenge: considered runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid
circuits, survivor gates, dead islands, phase-flow hits, carry residues, and
proof obligations.  The chosen carrier preserves the LRC branch-relocation
predicate plus the q=`14V` CRT witness predicate, while raw wall or delta
scalars forget the seam/complement distinction.

## Proof Pull

Use HYP-3484 as the surgery addendum to the HYP-3482/HYP-3483 random031 named
clause downstream of HYP-3479:

```text
hard-orbit discharge
<= separating-current transfer
 + random031 forbidden-seam flow clause.
```

The random031 clause should not ask phase flow to cross the hard pair.  It
should prove that a seven-owner forbidden seam in a mirror-punctured cylinder
cannot terminally glue once the lower-delta bypass, C=27 carry coordinate, and
component-cover escapes are retained.  HYP-3480's `0/4` branch-unit control
readout says this clause should not be replaced by the mirror-unit singleton
lemma for the six other zero-edge rows.

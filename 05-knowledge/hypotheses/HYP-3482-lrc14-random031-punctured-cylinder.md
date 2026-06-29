---
id: HYP-3482
title: LRC14 random031 punctured-cylinder forbidden-seam atlas
status: EVIDENCE
source: codex-2026-06-29
tangent: T1442
technique: LTI-442
tournament: LTT-342
script: 04-computation/lrc14_random031_punctured_cylinder_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_punctured_cylinder_codex_20260629.out
reflection: 07-reflections/lrc14-random031-punctured-cylinder-codex-20260629.md
related:
  - HYP-3481
  - HYP-3480
  - HYP-3479
  - HYP-3478
  - HYP-3477
  - HYP-3476
  - HYP-3460
  - HYP-3455
  - HYP-3451
  - HYP-3438
  - THM-523
---

# HYP-3482: LRC14 random031 punctured-cylinder forbidden-seam atlas

## Claim

`random_covering_031` should be treated as a mirror-punctured cylinder with a
forbidden seam, not as a failed wall-current row.

The row is the unique currentless hard-orbit overlap in HYP-3479.  Its
max-delta pair has no compatible q=`14V` phase-grid hits, while the same two
hard components are hit twelve times by lower-delta opposite-branch gates.  So
the hard pair behaves like a seam in the boundary-gluing model: the complement
carries the phase flow, and the seam itself is a forbidden crossing clause.

## Exact Atlas

Computed by
`04-computation/lrc14_random031_punctured_cylinder_codex_20260629.py` and
stored in
`05-knowledge/results/lrc14_random031_punctured_cylinder_codex_20260629.out`.

Row data:

```text
row=random_covering_031
V=173
P=(12,)
E=(0,4,26,53,60,80,89,103,115,118,128,150)
components=98
dead_components=4
low_rank_escape_components=94
component_class_hist={'both_alive':50,'branch0_only':22,'branch1_only':22,'dead_both':4}
```

Cell-model reading:

```text
cylinder_chi=0
dead_island_count=4
punctured_model_chi=-4
```

The four isolated dead islands form two mirror pairs.  Their owner grammar is:

```text
island 0: interval [113/840,55/406], mirror 3, owners (23,45), span 22
island 1: interval [57/406,69/490], mirror 2, owners (93,113), span 20
island 2: interval [421/490,349/406], mirror 1, owners (93,113), span 20
island 3: interval [351/406,727/840], mirror 0, owners (23,45), span 22
```

Thus the dead-island owner union is `(23,45,93,113)` and the span word is
`(22,20,20,22)`.

The hard pair is the seam:

```text
hard_orbit_components=(43,54)
hard_orbit_delta=7
typed_pair=('B0:9|B0:1','B1:1|B1:9')
structural_pair=(('B0|B0','branch0','two_sided',3,4),
                 ('B1|B1','branch1','two_sided',4,3))
hard_seam_owner_union=(23,45,93,113,147,169,173)
rescue_subset=(23,45,93,113,147,169)
rescue_graph_connected=True
```

The phase-flow pullback from HYP-3460 is:

```text
phase_witnesses=282
phase_actual_count=282
hard_gate_hits=0
same_component_lower_delta_hits=12
hard_component_hit_counter={(54,'branch0',2):6,(43,'branch1',2):6}
gate_delta_counter={2:118,3:40,4:62,5:16,6:6}
no_gate_hits=40
```

This is the main structural punchline: the q=`14V` flow never crosses the
max-delta seam, but it touches the same two components twelve times through a
lower-delta mirror bypass.

## Bold Reframes

### 1. Forbidden seam theorem target

Replace "the hard pair is a wall" by:

> A max-delta mirror seam with zero q=`14V` phase hits and mirror-balanced
> lower-delta same-component hits is a boundary-gluing clause.  Its complement,
> not the seam, carries the phase flow.

The proof target is to show that seam saturation cannot coexist with the
connected seven-owner rescue graph and the mirror-balanced lower-delta bypass
unless a low-rank escape remains.

### 2. Mirror-punctured cylinder theorem target

Replace "random031 is a strange exception" by:

> It is a cylinder with four mirror-paired punctures.  The four dead islands
> are isolated holes; the hard pair is the forbidden seam joining the two
> boundary branches; the 282 phase witnesses live on the seam complement.

In this model the Euler-characteristic shadow is not the terminal invariant.
The retained predicate is: every phase-flow component on the complement either
reaches a low-rank E/branch escape or must cross the forbidden seam.  Since the
stored flow crosses it zero times and still hits the two seam components through
lower-delta gates, the next theorem should classify that bypass as the terminal
discharge.

### 3. Additive rim versus doubling fold

The owner rim looks additive:

```text
(23,45), (93,113), (147,169) plus apex 173
span word=(22,20,22)
```

The phase pullback is multiplicative/two-adic:

```text
u = 2t mod 1
```

So the old `n+2` versus `n*2` recursion split becomes geometric here.  The
owner labels walk around the rim by additive jumps, while the branch coordinate
is produced by the two-adic fold.  The hard seam is where the additive owner
rim and multiplicative branch fold disagree most strongly.  That disagreement
is not raw obstruction; it is the place where the bypass should be proved.

## Experiment Design

Build the seam-complement graph:

1. Remove the two max-delta hard seam gates from the survivor-gate graph.
2. Route all `282` exact q=`14V` phase witnesses through the remaining
   branch-compatible gates.
3. Mark whether each connected phase-flow component reaches one of the `94`
   low-rank escape components before it would have to cross the seam.
4. Repeat for the other seven hard mirror orbits from HYP-3477 to distinguish
   zero-hit seams from genuine phase walls.

Expected theorem shape:

```text
max_delta_zero_hit_seam
+ connected seven-owner rescue graph
+ mirror-balanced lower-delta same-component bypass
=> terminal low-rank escape or forbidden seam contradiction
```

## Tournament Analysis

Vertex set: proof carriers, not runners.  Vertices are forbidden seam,
punctured-cylinder topology, seven-owner layered seam word, phase-branch bypass
worldlines, additive-rim/doubling-fold recursion lens, dead-island owner pairs,
and raw-count shadows.

Pairwise observable: retained terminal LRC payload, exact component geometry,
phase-flow compatibility, owner-layer detail, recursion interpretability, and
scalar-forgetting penalty.

Switch/gauge: orient toward the carrier that preserves more of the terminal
proof predicate; ties follow the Hamiltonian path

```text
C00_forbidden_seam_complement_flow
-> C01_mirror_punctured_cylinder_model
-> C02_seven_owner_layered_seam_word
-> C03_phase_branch_bypass_worldlines
-> C04_additive_rim_vs_doubling_fold_lens
-> C05_dead_island_owner_pairs
-> C06_raw_counts_shadow
```

Fingerprint:

```text
score_hist={8:1,55:1,56:1,60:1,64:2,67:1}
```

## Assumption Challenge

Do not take the tournament vertices to be runners, arcs, or raw gates.  Those
choices forget exactly the geometry now visible in `random031`.

Alternate vertex sets considered:

- dead islands / punctures
- seam components
- fixed circle sections
- section boundaries
- wall-crossing events
- residues and owner layers
- cover arcs
- phase-flow worldlines
- proof obligations

Preserved predicate: terminal discharge of the unique hard/currentless
`random031` overlap in the LRC14 proof stack.

Destroyed by bad quotients: interval addresses, mirror mate, branch
orientation, seven-owner seam layering, lower-delta bypass, phase-hit count,
and the additive-rim versus two-adic-fold distinction.

## Status

Evidence only, not a proof.  HYP-3482 turns `random031` from an exception name
into a finite topological proof packet.  The next hard work is the
seam-complement flow theorem.

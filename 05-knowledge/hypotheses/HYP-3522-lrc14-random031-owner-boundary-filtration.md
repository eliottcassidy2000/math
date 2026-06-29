---
id: HYP-3522
title: LRC14 random031 owner-boundary bracket filtration
status: EVIDENCE / exact owner-filtration certificate; not an LRC14 proof
source: codex-2026-06-29 refinement of HYP-3520 and HYP-3521, now guarded by HYP-3513, and continuing HYP-3512, HYP-3494, HYP-3511, HYP-3510, HYP-3493, HYP-3490, HYP-3486, HYP-3485, HYP-3484, HYP-3483, HYP-3482, HYP-3481, HYP-3477, HYP-3460, and HYP-3455
tangent: T1522
technique: LTI-522
tournament_technique: LTT-422
script: 04-computation/lrc14_random031_owner_boundary_filtration_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_owner_boundary_filtration_codex_20260629.out
reflection: 07-reflections/lrc14-random031-owner-boundary-filtration-codex-20260629.md
related:
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
  - HYP-3477
  - HYP-3460
  - HYP-3455
  - HYP-3402
  - HYP-3034
  - THM-523
  - OPEN-Q-108
---

# HYP-3522: LRC14 Random031 Owner-Boundary Bracket Filtration

## Claim

HYP-3520 makes the four-owner owner-current boundary debt
`(45,147,169,173)` exact, and HYP-3521 places the pure bypass inside the
terminal `230+40+12` certificate ledger.  HYP-3522 refines the owner side one
step further: two of those four seam-only owners are branch-boundary bracket
owners, leaving a smaller residual pair.

The pure-bypass owner problem in `random_covering_031` should not be stated as
"transport all seven seam owners through the lower-delta bypass."  That is the
wrong quotient.  The exact finite object is an owner filtration:

```text
seven-owner seam
= transport word
+ branch-boundary bracket lift
+ residual puncture/apex debt
```

Computed exactly:

```text
hard_seam_owners=(23,45,93,113,147,169,173)
pure_bypass_transport_owners=(23,93,113)
branch_boundary_bracket_owners=(23,93,147,169)
bracket_lift_owners=(147,169)
transport_only_owners=(113,)
residual_after_branch_boundary=(45,173)
```

Thus the terminal random031 owner target can shrink from HYP-3520's
four-owner cobordism word to a two-owner residual boundary lemma for
`(45,173)`, after transport constancy and branch-boundary bracketing are
proved.

## Exact Readout

Computed by:

```text
04-computation/lrc14_random031_owner_boundary_filtration_codex_20260629.py
05-knowledge/results/lrc14_random031_owner_boundary_filtration_codex_20260629.out
```

Gate-level persistence:

```text
hard_seam_gate_count=2
hard_gate_phase_hits=0
lower_bypass_gate_count=2
lower_bypass_gate_owner_union=(23,93,113)
lower_bypass_owners_equal_transport=True
seam_owners_avoid_transport_debt=True
```

Pure bypass stalk:

```text
component_size=12
branch_hist={0:6,1:6}
hit_components=(43,54)
endpoint_ranks=(2,)
bypass_owner_word_hist={(23,93,113):12}
mirror_pair_count=6
mirror_owner_word_persistent=True
```

Branch-order brackets:

```text
branch=0 bypass_u=(679,680,681,682,683,684)
         phases=(7,8,9,10,11,12)
         left owners=(93,147,169)
         right owners=(23,169)

branch=1 bypass_u=(527,528,529,530,531,532)
         phases=(2,3,4,5,6,7)
         left owners=(23,169)
         right owners=(93,147,169)
```

In both branch sheets the adjacent ordinary cells have owner intersection
`(169,)`.  Their union supplies `(23,93,147,169)`: the original transport
owners `23,93` plus the bracket lift `147,169`.  Owner `113` is transported
but not supplied by the branch bracket; owners `45,173` remain residual.

Owner-layer table:

```text
owner  layer
23     transport+branch_boundary
45     residual_dead_island_boundary
93     transport+branch_boundary
113    transport_only
147    branch_boundary_only
169    branch_boundary_only
173    residual_apex_boundary
```

## Integration Guardrails

Incoming HYP-3513 changes how the HYP-3490 firewall should be attached:
private-firewall status is already pure for existing colored axes `C/F/N/T`,
but the full HYP-3490 route is not reconstructed without route sidecar `R` or
a separate route theorem.  Therefore the HYP-3522 residual lemma should not
only say "compatible with HYP-3490"; it should carry the random031 route sidecar
or prove that the route is reconstructed from the owner filtration.

The extended HYP-3520 seam-sheaf canary gives the compression rule for this
filtration.  Safe component quotients are `flow_class`, `allowed_exit`,
`owner_union`, and `sheet_pgf_bucket`.  Unsafe quotients include owner-count,
endpoint-rank, branch histogram, component size, and mirror-closure alone.
HYP-3522 uses the safe side: transport, branch-boundary, residual, and route
sidecars are kept before any scalar compression is allowed.

## Proof Pull

The sharpened random031 terminal interface is now four lemmas:

1. Transport-word constancy: every pure-bypass cell and mirror mate carries
   exactly owners `(23,93,113)`.
2. Branch-boundary bracket lift: the ordinary cells adjacent to the bypass add
   exactly the missing bracket owners `(147,169)` while preserving the
   transport word.
3. HYP-3510/HYP-3511 separation: connected phase transport and free-hole
   bracketing discharge phase/no-gate cells without forcing `(45,173)` through
   the bypass transport word.
4. Residual two-owner boundary lemma: the remaining puncture/apex pair
   `(45,173)` cannot support a counterexample after the first three lemmas and
   the HYP-3490/HYP-3513 private-label firewall route sidecar is attached.

This gives a proof-facing target that is both smaller and more legal than the
old seven-owner gluing clause.

## Why This Matters

HYP-3493/HYP-3494 named owner-boundary persistence as the live obstruction, but
still left the target too large.  HYP-3522 says the bypass is not supposed to
carry every seam owner.  The bypass carries the stable transport word; the
ordinary branch boundary carries a bracket correction; the only live debt is
the two-owner residual.

This reframes the random031 topology as a short exact sequence of owner
charges:

```text
0 -> transport (23,93,113)
  -> transport+boundary (23,93,113,147,169)
  -> residual (45,173)
```

The sequence is not asserted as algebraic exactness yet; it is the finite
certificate shape the next proof should formalize.

## Tournament Analysis

Vertices are owner-boundary proof carriers, not runners or raw gate arcs.

Pairwise observable: retained owner charge, transport/boundary split, mirror
persistence, and quotient safety.

Switch/gauge: orient toward the carrier that preserves the most terminal proof
payload; ties use the filtration order.

Fingerprint:

```text
score_hist={18:1,70:1,76:1,82:1,87:1,94:1,98:1}
directed_3cycles=0
hamiltonian_path=
  transport_owner_word_constant
  -> mirror_pair_owner_word_persistence
  -> branch_boundary_bracket_lift
  -> residual_puncture_apex_debt
  -> free_hole_bracket_separation
  -> coarse_connected_phase_carrier
  -> raw_seven_owner_shadow
```

## Assumption Challenge

Candidate vertices considered: runners, hard gates, lower-delta bypass gates,
owner labels, dead islands, branch-boundary ordinary cells, mirror pairs,
free-hole packets, coarse phase components, and proof obligations.

Chosen vertices are owner-boundary carriers.  This preserves the random031
terminal predicate while deliberately forgetting raw runner order only after
transport, bracket, residual, and mirror sidecars have been attached.

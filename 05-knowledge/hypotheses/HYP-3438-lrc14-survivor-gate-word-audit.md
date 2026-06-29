---
id: HYP-3438
title: LRC14 survivor-gate word audit
status: EVIDENCE / exact survivor-gate classification; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3436 minimal bad-core cover extractor
tangent: T1399
technique: LTI-399
tournament_technique: LTT-299
script: 04-computation/lrc14_survivor_gate_word_audit_codex_20260629.py
result: 05-knowledge/results/lrc14_survivor_gate_word_audit_codex_20260629.out
reflection: 07-reflections/lrc14-survivor-gate-word-audit-codex-20260629.md
related:
  - HYP-3437
  - HYP-3436
  - HYP-3435
  - HYP-3434
  - HYP-3431
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3422
  - HYP-3417
  - HYP-3129
  - THM-523
---

# HYP-3438: LRC14 Survivor-Gate Word Audit

HYP-3438 executes the immediate HYP-3436 local-to-global follow-up.  Instead of
collapsing the HYP-3436 survivor mass to a scalar, it decomposes every mixed
`E_safe` component into an exact gate word whose letters are adjacent bad-core
blocks and survivor gaps.  Each survivor gate retains:

```text
survivor interval
+ left/right bad-core block if present
+ minimal B0/B1 odd-owner cover on each bad block
+ endpoint wall labels
+ parent even wall labels
+ legal sidecar route
```

## Exact Readout

On the `135` primitive covering rows inherited from HYP-3436:

```text
total_even_safe_components=17164
fully_bad_components=3770
mixed_components=6228
clean_components=7166
survivor_gates=8702
rows_with_mixed_owner_residual=79
route_hist={
  edge_singleton_parent_gate: 5950,
  edge_survivor_residual: 1080,
  mixed_owner_residual: 878,
  owner_current_small_delta: 794
}
survivor_branch_mask_hist={both:1064, branch0:3819, branch1:3819}
gate_adjacency_hist={left_bad_edge:3515, right_bad_edge:3515, two_sided:1672}
```

The dominant legal object is now an edge boundary gate: a survivor gap adjacent
to one `(1,1)` bad-core block and touching the parent even wall.  This is a
meaningful correction to the earlier corridor intuition.  The canonical
`{1,...,11,13,84}` row has exactly four mixed gates, all
`edge_singleton_parent_gate`.

The branch-mask ledger from the merged audit is also useful: branch-0 and
branch-1 survivor gates occur symmetrically (`3819` each), while `1064` gates
are genuinely both-branch survivors.  The hardest residuals are therefore not
branch-asymmetry artifacts; they are adjacent-cover and endpoint-gluing
obligations.

## Canonical Mod-35 Law

The canonical family

```text
S_m = {1,2,3,4,5,6,7,8,9,10,11,13,84m}
```

was audited through a full residue period `m=1..35`.  The observed law has no
failures:

```text
outer_13_7 edge gates exist exactly when 7 does not divide m
inner_11_5 edge gates exist exactly when 5 does not divide m
if 35 divides m, no mixed survivor gates remain and the row survives
through clean E_safe components
```

Bucket count across one period:

```text
outer_plus_inner: 24 residue classes
outer_only:        6 residue classes
inner_only:        4 residue classes
clean_only:        1 residue class
```

This suggests a much sharper canonical theorem target than the previous
singleton-tail statement:

```text
Prove the mod-35 canonical gate law for all m.
```

The proof should split the outer `(13,7)` gates from the inner `(11,5)` gates
and then prove that the simultaneous loss at `35 | m` is not a counterexample
but a clean-component survivor case.

## General Proof Target

Every mixed `E_safe` component emits a finite survivor gate.  A global bad-core
gluing proof must show that every gate is one of:

```text
edge singleton-parent gate
interior endpoint-even-wall gate
small owner-current delta
named residual routed to HYP-3437/HYP-3439/HYP-3440
two-adic descent / exact-period / state-lift / signed-SPEC debt
```

The remaining hard mass is no longer the existence of survivors but the
classification of `878` mixed-owner residual gates and `1080` edge-survivor
residual gates.  Those should be fed into the overlap-tax Menger route
(HYP-3437), rescue-core bridge (HYP-3439), endpoint-cut vocabulary (HYP-3440),
or owner-current/two-adic/state-lift exits.

## Tournament Analysis

Vertices are survivor-gate proof carriers, not runners or arcs.  Alternate
vertices considered: runners, gaps, fixed sections, section boundaries,
wall-crossing events, residues, cover arcs, endpoint walls, branch-bad owners,
cover-pair deltas, mixed even components, survivor gaps, and proof obligations.

Pairwise observable:

```text
predicate retention + exact gate word + adjacent cover delta
+ endpoint wall payload + legal sidecar route
```

Tie Hamiltonian path:

```text
G00_gate_word_exactness
-> G01_adjacent_cover_delta
-> G02_endpoint_wall_alternation
-> G03_corridor_singleton_recognizer
-> G04_endpoint_spine_sidecar
-> G05_owner_current_route
-> G06_overlap_tax_menger_route
-> G07_two_adic_loss_route
-> G08_signed_SPEC_route
-> G09_raw_survivor_measure
```

Fingerprint: score histogram `{6:1,98:1,104:1,109:1,112:1,114:1,124:1,126:1,129:1,140:1}`,
no directed `3`-cycles, and `9` edge flips against the illegal
raw-survivor-first gauge.

Quotient guardrail: this quotient preserves the local branch-relocation
predicate inside mixed `E_safe` components.  It destroys raw runner order and
scalar survivor mass, while retaining endpoint and owner payload as gate
labels.  A future quotient may forget those labels only after reconstructing
or legally routing them.

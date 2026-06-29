# HYP-3476: LRC14 Dead-Cover Pair-Current Exception Audit

**Status:** EVIDENCE / two-gate exception audit; not a proof.
**Created:** codex-2026-06-29
**Owner:** codex
**Artifacts:**
- `04-computation/lrc14_dead_cover_pair_current_codex_20260629.py`
- `05-knowledge/results/lrc14_dead_cover_pair_current_codex_20260629.out`
- `07-reflections/lrc14-dead-cover-pair-current-codex-20260629.md`

## Question

HYP-3472 proves on the current `135`-row bank that every dead-component row
has a low-rank E/branch gate touching the dead-cover projection, but a single
gate gives a projection-edge cut only on `123/130` dead rows and a separating
current only on `121/130`.

This hypothesis asks whether the named HYP-3472 exceptions are an artifact of
using one gate at a time.  For the edge-cut exception rows

```text
random_covering_001
random_covering_031
random_covering_039
random_covering_062
random_covering_074
random_covering_086
random_covering_101
```

and the two additional separating-current exceptions

```text
covering_AP_with_84
ap_omit_12_tail_84x01
```

test whether a pair of low-rank E/branch gates, or a mirror-orbit-aware pair,
removes enough adjacent blocker labels from the HYP-3451 dead-cover projection
to create a projection-edge cut or a separating branch current.

## Result

The targeted audit ran on the nine HYP-3472 exception rows:

```text
bank_rows=135
bank_dead_rows=130
audit_scope=hyp3472_exception_targets
audited_dead_rows=9
audited_single_e_branch_cuts=330
audited_pair_e_branch_cuts=9113
exact_mirror_pairs=165
```

Two-gate E/branch currents do not close the seven random projection-edge
exceptions:

```text
edge_exception_rows_closed_by_pair_edge_cut=0/7
edge_exception_rows_closed_by_pair_separating_current=0/7
```

They do close the two AP separating-only exceptions:

```text
separating_exception_rows_closed_by_pair_separating_current=2/9
['covering_AP_with_84', 'ap_omit_12_tail_84x01']
```

The sharper diagnostic is the edge-support gap.  The seven random exceptions
have edgeless dead-cover projections:

```text
random rows: projection edges = 0
edge_support_label_count = 0
gate_edge_support_label_count = 0
```

Thus no union of E/branch gate-adjacent labels can create a projection-edge
cut there; the projection has no edge to cut.  The AP rows are different:
their dead-cover projection is connected with `90` edges, `14` edge-support
labels, and the four AP gates touch `8` of those labels.  A non-mirror pair
with labels

```text
B0:11, B0:13, B1:5, B1:7
```

removes `36` projection edges and gives `largest_drop=2`,
`component_gain=2` on both `covering_AP_with_84` and
`ap_omit_12_tail_84x01`.  The exact mirror pair removes `40` edges but does
not separate, so mirror symmetry is an edge-cut carrier but not the AP
separating-current carrier.

## Method

Reuse the HYP-3472 projection and gate-adjacency code.  For each target row,
enumerate unordered pairs of low-rank E/branch survivor gates.  For each pair,
remove the union of adjacent blocker labels from the dead-cover projection and
measure:

- remaining edge count and component count;
- largest-component drop;
- whether the pair is a projection-edge cut;
- whether the resulting components are branch-separated;
- whether the pair is the exact mirror partner of HYP-3475 or belongs to the
  same mirror-orbit family.

## Assumption Challenge

Tournament vertices will be proof obligations and finite current carriers, not
runners or arcs.  Candidate vertex sets considered before fixing this quotient:
runners, gaps, fixed circle sections, section boundaries, wall-crossing
events, residues, cover arcs, Fourier modes, matroid circuits, survivor gates,
mirror orbits, gate pairs, blocker labels, graph cuts, and formal proof
obligations.

The quotient preserves the predicate:

```text
dead-cover obstruction on a target row has a <=2-gate boundary-current carrier.
```

It destroys individual interval geometry, ordered branch orientation inside a
mirror pair, and full row identity unless the script retains typed sidecars.
HYP-3474 is the guardrail: count or color profiles are telemetry unless they
reconstruct the projection/cut predicate.

## Proof Pull

HYP-3476 splits the HYP-3472 exceptions into two proof types.

1.  AP separating-only rows are local pair-current packets: single gates cut
    edges but do not separate; selected two-gate pairs separate.
2.  The seven random rows are not pair-current failures.  They are zero-edge
    singleton dead projections, so the next theorem must use a different
    invariant: owner-current imbalance across isolated dead components,
    HYP-3455/HYP-3451 gluing/conductance, HYP-3460 phase-branch bypass,
    two-adic descent, signed-SPEC/Rprime, or state-lift debt.

The HYP-3473 terminal interface should therefore not promise a universal
`<=2`-gate projection cut.  A better interface is:

```text
dead row =>
  touching E/branch gate
  and either edge-support pair-current packet
  or zero-edge singleton-current terminal debt.
```

Incoming connected work after this audit further subdivides the seven
zero-edge random rows.  HYP-3477's hard mirror-orbit reservation intersects
this set at `random_covering_031`.  The colored gate unit-delta split
(`HYP-3472-lrc14-colored-gate-unit-delta-split`) intersects it at
`random_covering_039` and `random_covering_074`, both in the 19-row
cover-delta sidecar packet.  The remaining zero-edge rows
`random_covering_001`, `random_covering_062`, `random_covering_086`, and
`random_covering_101` should be tested first as unit-delta singleton-current
rows before adding new gate tuple machinery.

## Links

HYP-3476 depends on HYP-3475, HYP-3474, HYP-3473, HYP-3472, HYP-3471,
HYP-3462, HYP-3470, HYP-3461, HYP-3460, HYP-3459, HYP-3458, HYP-3455,
HYP-3453, HYP-3451, HYP-3450, HYP-3438, HYP-3436, HYP-2595, THM-523,
T1436, LTI-436, LTT-336, and OPEN-Q-108.

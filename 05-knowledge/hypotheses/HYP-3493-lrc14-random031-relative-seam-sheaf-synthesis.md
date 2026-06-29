# HYP-3493: LRC14 Random031 Relative Seam-Sheaf Synthesis

**Status:** EVIDENCE / synthesis plus executable seam-sheaf scaffold; not a proof.

**Created:** codex-2026-06-29

## Claim

`random_covering_031` should be promoted from a forbidden-seam graph picture to
a relative seam sheaf over the mirror-punctured cylinder.

The hard max-delta pair is not a wall.  It is a deleted seam whose complement
carries the phase flow.  The proof object should keep three stalks at once:

1. **Fiber stalk.**  HYP-3486 gives the exact two-adic cylinder graph after
   seam deletion: `282` q=`14V` witness cells on `258` occupied base fibers,
   split into `242` endpoint-rank-2 routed cells, `40` free-hole cells, and
   one pure `12`-cell lower-delta bypass component on hard components
   `(43,54)`.
2. **Label-firewall stalk.**  HYP-3490 says adjacent-label current deletion is
   structurally blocked on random031: every E/branch-touched blocker label is
   private, and random031 is the only row that is both private-label and hard.
3. **Relative-topology stalk.**  HYP-3485/HYP-3484/HYP-3482/HYP-3481 say the
   correct base is `cylinder - four mirror-paired punctures - forbidden
   max-delta seam`, with the seven owner labels as seam boundary charge and
   the lower-delta bypass as phase flow on the complement.

The predicted proof interface is a local-to-global gluing theorem for these
stalks, not another scalar wall count.

```text
relative_seam_sheaf_cell =
  (u_fiber,
   branch_sheet,
   mirror_partner,
   flow_class in {rank2_routed, free_hole, pure_bypass},
   endpoint_rank,
   owner_boundary_word,
   private_label_status,
   relative_cycle_id,
   sheet_pgf_bucket,
   allowed_exit)
```

Legal exits:

```text
rank2_routed -> endpoint-rank-2 discharge
free_hole -> mirror free-hole lemma
pure_bypass -> 12-cell bypass lemma plus owner-boundary debt
private_label_firewall -> no projection-current route; keep seam sheaf
residual -> named owner-current / two-adic / signed-SPEC / state-lift debt
```

## Why this is new after the concurrent work

HYP-3485 supplied the connection atlas, HYP-3486 supplied the exact finite
fiber graph, and HYP-3490 supplied the current-deletion firewall.  Read
together, they say the remaining random031 clause is not one obstruction but a
three-layer compatibility problem:

```text
fiber trichotomy + private labels + relative seam complement
```

The useful theorem should therefore resemble a relative sheaf exactness
statement: local witness cells are discharged in their fiber class unless the
owner boundary word emits a named global debt.  Projection-current deletion is
not an available simplification on this row.

## Imports from past and random-topic searches

The repo search connects this packet to older work in a specific, guarded way.

- HYP-3025/HYP-3034 suggest using closed arc-Cech and path-lift persistence,
  but only after owner deletion and endpoint labels are retained.  Raw homology
  of the punctured cylinder is too coarse.
- HYP-3140 suggests a local fiber-PGF or sheet-count moment, but HYP-3486 gives
  the guardrail: vertical half-turn gluing is not class-preserving and mixes
  ordinary with free-hole cells.  Sheet moments must live on the legal
  horizontal+mirror graph.
- HYP-3422/HYP-3428/HYP-3483 say `u=2t` is an address coordinate for the
  `n*2` phase pullback.  It must be paired with the `n+2` owner-boundary seam
  word; neither recursion wins alone.
- HYP-3311/HYP-3408 style `Q(sqrt(-7))`, residue, and two-adic sidecars are
  legal only as owner-lift stability data.  They are not scalar substitutes
  for the seam sheaf.
- The older non-archimedean / local-to-global motif is useful as a gluing
  metaphor: prove local stalk discharge and then check owner-boundary
  compatibility.  Do not turn it into a p-adic congruence proof without the
  owner labels.
- The van-der-Corput / discrepancy motif should be tested on the ordered
  six-hit bypass blocks from HYP-3483, not on random equidistribution.  The
  relevant discrepancy object is deterministic phase ordering on the bypass
  component.

## Executed scaffold

Companion script:

```text
04-computation/lrc14_random031_relative_seam_sheaf_codex_20260629.py
05-knowledge/results/lrc14_random031_relative_seam_sheaf_codex_20260629.out
```

The scout imports HYP-3486 and annotates each legal horizontal+mirror component
as a relative seam-sheaf stalk.  Exact readout:

```text
component_count=79
flow_class_hist={'free_hole':14,'pure_bypass':1,'rank2_routed':64}
allowed_exit_hist={
  'endpoint_rank2_discharge':64,
  'mirror_free_hole_packet':14,
  'pure_bypass_plus_owner_boundary_debt':1
}
mirror_closed_components=79/79
owner_union_size_hist={0:14,2:26,3:15,4:18,5:3,6:3}
endpoint_rank_hist={2:65}
sheet_pgf_bucket_count=10
```

There are no `mixed_or_debt` components in the current finite carrier.  Every
legal component is mirror-closed, so the next Cech/path-lift test should target
owner-boundary persistence rather than mirror closure.  The pure bypass stalk
is

```text
R01 class=pure_bypass size=12
u=[527,684] branches={0:6,1:6}
hits=(43,54) owners=(23,93,113)
seam_debt=(45,147,169,173)
exit=pure_bypass_plus_owner_boundary_debt
```

This exactly matches the HYP-3483 span: the bypass carries the `n*2` phase
block while the missing owner word is the `n+2` seam debt.

## Assumption challenge

Do not default the tournament vertices to runners or arcs.  Candidate vertices
considered here are runners, gaps, fixed circle sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits,
witness cells, fiber components, owner-boundary words, and proof obligations.

Chosen vertex set: relative seam-sheaf proof packets.  This quotient preserves
the terminal random031 predicate: every q=`14V` phase-flow cell must either
reach an allowed low-rank/free-hole/bypass exit or emit named debt before a
forbidden seam crossing is used.  It destroys raw interval geometry, exact
cell positions inside a fiber, and full blocker-label incidence unless those
are restored by the listed stalk fields.

Challenged assumption: the hard pair is not a geometric wall separating phase
flow.  HYP-3484/HYP-3486 show the phase flow avoids it; the seam is forbidden
boundary data, while the complement carries the actual dynamics.

## Tournament Analysis

Vertices are proof carriers:

```text
S00_relative_seam_sheaf_packet
S01_fiber_graph_trichotomy
S02_private_label_firewall
S03_ordered_twoadic_bypass_blocks
S04_cech_owner_persistence
S05_sheet_pgf_moment
S06_nonarchimedean_owner_lift_gate
S07_raw_scalar_shadow
```

Pairwise observable: retained terminal random031 predicate; legality of the
quotient; retention of fiber class, branch sheet, mirror mate, owner boundary,
private-label status, Cech/path-lift persistence, and sheet-count payload.

Switch/gauge: orient toward the carrier with more retained proof payload and
fewer illegal quotient losses.  Ties follow the displayed Hamiltonian path.

Designed synthesis fingerprint: transitive carrier tournament with score
histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}`, singleton SCCs, no directed
3-cycles, and Hamiltonian path

```text
S00_relative_seam_sheaf_packet
 -> S01_fiber_graph_trichotomy
 -> S02_private_label_firewall
 -> S03_ordered_twoadic_bypass_blocks
 -> S04_cech_owner_persistence
 -> S05_sheet_pgf_moment
 -> S06_nonarchimedean_owner_lift_gate
 -> S07_raw_scalar_shadow
```

This fingerprint is a route-order declaration, not a new row-bank computation.

## Next experiments

1. Add a Cech/path-lift persistence test: deleting the forbidden seam should
   not merge components unless the owner-boundary sidecar records the debt.
2. Add a sheet-PGF moment on legal components only, then check whether the
   rank-2/free-hole/bypass trichotomy is separated by low-degree moments.
3. Run the same sheaf table on the other HYP-3477 hard mirror orbits as
   negative controls.

## Pointers

HYP-3493, HYP-3490, HYP-3486, HYP-3485, HYP-3484, HYP-3483, HYP-3482,
HYP-3481, HYP-3480, HYP-3479, HYP-3477, HYP-3460, HYP-3455, HYP-3451,
HYP-3438, HYP-3428, HYP-3422, HYP-3140, HYP-3034, HYP-3025, HYP-3023,
HYP-3311, HYP-3408, THM-523, T1453, LTI-453, LTT-353, OPEN-Q-108.

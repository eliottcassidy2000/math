---
id: HYP-3480
title: LRC14 zero-edge singleton-current audit
status: EVIDENCE / finite singleton-current certificate audit plus Lean ledger; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3479 hard-orbit current join, HYP-3478 small-touch geometry, HYP-3476 pair-current/router, and HYP-3477 hard-orbit discharge
tangent: T1440
technique: LTI-440
tournament_technique: LTT-340
script: 04-computation/lrc14_zero_edge_singleton_current_codex_20260629.py
result: 05-knowledge/results/lrc14_zero_edge_singleton_current_codex_20260629.out
lean_module: 04-computation/lean/TournamentH7/TournamentH7/LRCSingletonCurrentLedger.lean
lean_result: 05-knowledge/results/lrc14_singleton_current_lean_codex_20260629.out
reflection: 07-reflections/lrc14-zero-edge-singleton-current-codex-20260629.md
related:
  - HYP-3486
  - HYP-3485
  - HYP-3490
  - HYP-3484
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3479
  - HYP-3478
  - HYP-3477
  - HYP-3476
  - HYP-3475
  - HYP-3472
  - HYP-3471
  - HYP-3460
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3438
  - THM-523
  - OPEN-Q-108
---

# HYP-3480: LRC14 Zero-Edge Singleton-Current Audit

HYP-3476 shows that the random boundary-current exceptions are not missing a
larger E/branch pair cut: their dead-cover projections are edgeless singleton
packets.  HYP-3478 shows that the six non-hard rows are mirror-paired
singleton owner packets, while HYP-3479 isolates `random_covering_031` as the
unique hard/currentless overlap.  This audit checks the component-level
singleton-current carrier that remains after those splits.

## Scope

Primary rows:

```text
random_covering_001
random_covering_039
random_covering_062
random_covering_074
random_covering_086
random_covering_101
```

Control row:

```text
random_covering_031
```

The script joins HYP-3472 dead-cover projections, HYP-3476 route labels,
HYP-3478 singleton mirror geometry, HYP-3479 hard-overlap flags, and
unit-delta/two-adic owner sidecars.  It tests whether each zero-edge singleton
component has a complete branch-unit E/branch touch and whether mirror partner
components admit mirror-compatible unit gate pairs.

Existing HYP-3478 companions already split the six primary rows into three
clean best-touch rows (`062`, `086`, `101`), asymmetric branch-unit row `001`,
and cover-delta sidecar rows `039`, `074`.  HYP-3480 preserves that
row-minimum/best-touch split but strengthens it: component-level complete
branch-unit mirror pairs exist on all six primary rows.  HYP-3490 now explains
why all seven random projection-edge exceptions resist adjacent-label
pair-current cuts: every touched blocker label is private.  This audit supplies
the six private/no-hard singleton-current terminal packet, while
HYP-3486/HYP-3485/HYP-3484/HYP-3483/HYP-3482/HYP-3481 supply the
`random_covering_031` control as a seam-complement fiber graph, relative
topology bridge, forbidden-seam flow, controlled recursion-flow span,
punctured-cylinder seam, and mirror-punctured topology packet.

## Result

Exact aggregate readout:

```text
audited_rows=7 target_rows=6 control_rows=1
route_hist={'random031_overlap_hard_and_currentless': 1, 'small_touch_no_hard_current_exception': 6}
terminal_class_hist={'mirror_unit_singleton_packet': 4, 'mirror_unit_singleton_packet_cover_delta_min_shadow': 2, 'random031_hard_currentless_control': 1}
target_projection_edge_hist={0: 6}
target_min_gate_kind_hist={'branch_unit_delta': 4, 'delta_sidecar_packet': 2}
target_dead_components=14
target_components_with_complete_branch_unit_touch=14/14
target_mirror_pairs_with_branch_unit_mirror_gate=7/7
cover_delta_min_shadow_rows_with_unit_certificate=('random_covering_039', 'random_covering_074')
control_components_with_complete_branch_unit_touch=0/4
control_mirror_pairs_with_branch_unit_mirror_gate=0/2
hard_currentless_control_rows=('random_covering_031',)
```

Lean ledger:

```text
04-computation/lean/TournamentH7/TournamentH7/LRCSingletonCurrentLedger.lean
05-knowledge/results/lrc14_singleton_current_lean_codex_20260629.out
```

The Lean module records the seven audited rows, the `4+2+1=7` terminal split,
the `14/14` component-touch and `7/7` mirror-pair readouts, and the HYP-3480
tournament carrier scores.  The stored build output reports no axioms for:

```text
hyp3480_target_row_count
hyp3480_random031_is_control
hyp3480_audited_row_partition
hyp3480_all_target_components_have_branch_unit_touch
hyp3480_all_target_mirror_pairs_have_unit_gate
hyp3480_dispatch_complete
hyp3480_dispatch_matches_counts
```

The carrier-count and score-extremal checks use native-decision axioms, in the
same style as the neighboring HYP-3479 finite ledger.

Thus the six small-touch/no-hard rows share a sharper terminal carrier than
HYP-3478's row-minimum split suggested: even `random_covering_039` and
`random_covering_074`, whose absolute minimum E/branch gates are
cover-delta-sidecar gates, still have complete branch-unit mirror pairs at the
component level.  The component-level packet has seven swapped-owner mirror
pairs:

```text
001: (165,179), (81,153)
039: (63,129)
062: (9,81)
074: (15,99)
086: (169,133)
101: (7,175)
```

The control row `random_covering_031` has four singleton dead components and
zero projection edges, but only two components have any complete E/branch
touch and none have a complete branch-unit touch.  It remains exactly the
HYP-3455/HYP-3460/HYP-3479 hard/currentless gluing clause, now refined by
HYP-3486/HYP-3485/HYP-3484/HYP-3483/HYP-3482/HYP-3481 into a random031
seam-complement/fiber/forbidden-seam/recursion packet.

## Proof Pull

A plausible next finite theorem is:

```text
mirror-swapped singleton B0/B1 owner pair
+ mirror-compatible complete branch-unit E/branch gate pair
+ route sidecar R
+ owner residue/two-adic word
=> legal singleton-current terminal packet.
```

If that lemma is proved, the non-AP currentless random frontier outside the
already named hard control collapses to the single mirror-unit packet.  The
remaining non-AP hard/currentless debt is `random_covering_031`, now routed to
the HYP-3486 fiber-graph trichotomy, HYP-3485 relative topology bridge,
HYP-3484 forbidden-seam flow lemma, HYP-3483 controlled-span lemma, HYP-3482
seam lemma, and HYP-3481 topology lemma.  HYP-3490 supplies the adjacent-label
firewall that separates these two lanes before either terminal packet fires.

## Tournament Analysis

Vertices are terminal singleton-current proof carriers, not runners, raw row
names, or scalar gate counts.  Candidate vertices include owner-label pairs,
dead singleton components, unit-delta gate words, route labels, hard-overlap
flags, two-adic residue payloads, signed-SPEC exits, and formal proof
obligations.  The quotient preserves whether a zero-edge row has a mirror-unit
singleton-current certificate or is the named random031 hard/gluing overlap
whose topology/recursion sidecar is
HYP-3486/HYP-3485/HYP-3484/HYP-3483/HYP-3482/HYP-3481.  The HYP-3490
private-label firewall explains why neither lane should keep searching for a
larger adjacent-label projection-edge current.

Tournament fingerprint from the script:

```text
score_hist={5: 1, 8: 1, 45: 1, 52: 1, 59: 2, 62: 1, 67: 1}
directed_3cycles=0
hamiltonian_path=Z00_mirror_unit_singleton_current_packet -> Z01_component_complete_touch_certificate -> Z02_random031_hard_control_clause -> Z03_route_sidecar_R_join -> Z04_cover_delta_min_gate_shadow -> Z05_owner_residue_two_adic_shadow -> Z06_raw_zero_edge_projection_count -> Z07_raw_row_name_list
```

---
id: HYP-3044
title: LRC14 residual topology-exception teeth
status: EVIDENCE / finite exception audit and theorem target; not a proof
source: codex-2026-06-26-S207
tangent: T1125
script: 04-computation/lrc14_residual_exception_teeth_codex_s207.py
result: 05-knowledge/results/lrc14_residual_exception_teeth_codex_s207.out
related:
  - HYP-3041
  - HYP-3040
  - HYP-3039
  - HYP-3038
  - HYP-3037
  - HYP-3036
  - HYP-3035
  - HYP-3034
  - HYP-3033
  - HYP-3032
  - HYP-3031
  - HYP-3030
  - HYP-3029
  - HYP-3028
  - HYP-3027
  - HYP-3024
  - HYP-2963
  - THM-572
  - LTI-192
  - LTT-090
  - OPEN-Q-108
---

# HYP-3044: Residual Topology-Exception Teeth

## Claim

The two HYP-3035 residual fibers not separated by compact arc topology are
not new residual route families.  They are two owner-strip collars that can
be stated as a finite exception lemma.

Read through HYP-3040's hidden statement ledger and HYP-3039's
hidden-coordinate rule, compact arc topology is the visible coordinate and the
owner stalk / primitive deck is the hidden coordinate that must be retained
when topology forgets too much.  HYP-3038's q=23 drop/add square gives the
same warning in a different local grammar: raw scalar coincidence is harmless
only after the endpoint-owner coordinate is reattached.

HYP-3041's AP-tail puncture atlas gives the closest family-scale analogue:
a mod-14 owner strip can be correct and still forget the exact-period
puncture clock.  The topology exceptions below appear to be the same kind of
controlled-forgetting event, with compact topology remembering the collar
shape while the owner stalk and primitive deck remember the route-splitting
period data.

On the S194/S199 residual ledger, compact arc topology fails on exactly two
same-topology buckets:

```text
single swap 9->99    COVERING-MOMENT
single swap 9->155   Q-WITNESS

single swap 11->121  COVERING-MOMENT
single swap 11->163  Q-WITNESS
```

All four rows are strict-open.  They are all single-swap collars, at drops
`9` and `11`.  In both buckets, the covering row has zero primitive safe mass
for `2 <= q <= 13`, while the Q-witness row has first primitive safe
denominator equal to the dropped speed and to its `q_threshold`.  The coarse
largest-safe-component stalk splits both buckets.

Thus a topology failure after the status gate should be read as:

```text
compact arc topology fails
  => same-topology single-swap collar
  => owner-stalk split
  => primitive q<=13 deck split
```

The candidate theorem is a local exception lemma: if a post-status residual
coarse ET+unit fiber is not split by compact arc topology, then it is an
owner-stalk collar whose primitive `q<=13` deck also separates the
Q-witness/covering route.

## Computation

Script:

```text
04-computation/lrc14_residual_exception_teeth_codex_s207.py
```

Stored output:

```text
05-knowledge/results/lrc14_residual_exception_teeth_codex_s207.out
```

The script imports the S199 residual tooth atlas and the S200 primitive-period
scheduler.  It rebuilds only the `38` stored residual packets and filters the
topology-mixed buckets.

Key census:

```text
topology_exception_fibers=2
topology_exception_buckets=2
topology_exception_rows=4
exception_drop_speeds=[9, 11]
all_exception_rows_single_swap=True
exception_status_counts={'open': 4}
covering_zero_primitive_deck_2_13=True
q_witness_first13_equals_drop_and_q0=True
coarse_stalk_splits_all_exceptions=True
primitive_deck_2_13_splits_all_exceptions=True
```

The split audit on the four exception rows:

```text
arc_topology_compact           buckets=2  mixed_route=2
coarse_safe_stalk              buckets=4  mixed_route=0
exact_safe_stalk               buckets=4  mixed_route=0
magnitude_cocycle              buckets=4  mixed_route=0
q_or_covering_certificate      buckets=3  mixed_route=0
first_primitive_safe_q_2_13    buckets=3  mixed_route=0
primitive_deck_2_13            buckets=3  mixed_route=0
route_label_sink               buckets=2  mixed_route=0
```

## Tournament Analysis

Vertices are exception repair carriers, not runners.  The pairwise observable
records exception coverage, route split, non-route legality, local owner
signal, primitive-period alignment, compression, and proof cost.  The switch
is majority comparison of these carrier vectors.

The tournament is transitive:

```text
score_hist={0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1}
directed_3cycles=0
hamiltonian_path_count=1
score_order=topology_then_owner_stalk_rule
  > primitive_deck_2_13
  > coarse_safe_stalk
  > exact_safe_stalk
  > arc_topology_compact
  > route_label_sink
```

The ranking is not a proof by itself.  Its role is to keep the proof-object
honest: the preferred carrier is not "route label" and not "topology alone",
but the rule that topology is first, and the only topology failures are then
owner-stalk exceptions with an independent primitive-period split.

## Next Pull

Add packet sidecar fields:

```text
residual_topology_exception_id
topology_exception_drop
topology_exception_stalk_key
topology_exception_first_primitive_q
topology_then_owner_stalk_rule
```

Then prove, familywise if possible, that no other S194/S199-style residual
topology failure exists beyond the two owner-stalk collars above.

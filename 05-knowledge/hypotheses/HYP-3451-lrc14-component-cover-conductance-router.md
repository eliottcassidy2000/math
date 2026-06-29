---
id: HYP-3451
title: LRC14 component-cover conductance router
status: EVIDENCE / graph-conductance proof router; not an LRC14 proof
source: codex-2026-06-28 continuation of HYP-3450 component-cover obstruction extractor
tangent: T1411
technique: LTI-411
tournament_technique: LTT-311
script: 04-computation/lrc14_component_cover_conductance_router_codex_20260628.py
result: 05-knowledge/results/lrc14_component_cover_conductance_router_codex_20260628.out
reflection: 07-reflections/lrc14-component-cover-conductance-router-codex-20260628.md
related:
  - HYP-3450
  - HYP-3435
  - HYP-3434
  - HYP-3432
  - HYP-3429
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3422
  - HYP-3417
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3451: LRC14 Component-Cover Conductance Router

## Claim

HYP-3450 identifies the local obstruction object: each dead component of
`E_safe` carries a branch-0 minimal odd-bad cover and a branch-1 minimal
odd-bad cover, while each audited row still has an endpoint-rank `<= 2`
survivor.

HYP-3451 turns those paired covers into a graph ledger.  Dead components are
vertices.  Two dead components are adjacent when their minimal paired covers
share a branch-coloured blocker such as `B0:5` or `B1:13`.  The resulting
component/blocker projection is the finite object on which a Menger cut,
Green-current certificate, or algebraic-connectivity obstruction should live.

The strongest new lesson is that "largest paired-cover rank" is not the same
as "closest to counterexample."  The unique rank-`6` row has many escape
components.  The dangerous rows are the AP-with-`84m` tails: they have very
high dead-component fraction, connected dead-cover projection, and only four
low-rank escapes.

## Exact Readout

Script:

```text
04-computation/lrc14_component_cover_conductance_router_codex_20260628.py
```

Stored result:

```text
05-knowledge/results/lrc14_component_cover_conductance_router_codex_20260628.out
```

Aggregate audit on the HYP-3450 bank:

```text
rows_audited=135
rows_with_low_rank_escape=135/135
max_dead_pair_rank=6 at random_covering_082
max_dead_fraction=0.962963 at ap_omit_12_tail_84x05
blocker_entropy_range=[2.000000,3.872947]
effective_blocker_range=[4.000000,14.651206]
```

Highest danger rows by the current finite router:

```text
ap_omit_12_tail_84x02:
  components=48
  dead=44
  alive=4
  low_rank_escape=4
  max_pair_rank=3
  projection_components=1
  lambda2_largest_projection=4.386932
  cheeger_proxy=0.104451

covering_AP_with_84 / ap_omit_12_tail_84x01:
  components=26
  dead=22
  alive=4
  low_rank_escape=4
  max_pair_rank=3
  projection_components=1
  lambda2_largest_projection=3.046340
  cheeger_proxy=0.138470
```

The largest paired-cover rank row is not the tight proof target:

```text
random_covering_082:
  components=138
  dead=32
  alive=106
  low_rank_escape=106
  max_pair_rank=6
  danger_score=0.013126
  lambda2_largest_projection=2.461547
```

## Proof Pull

The conductance formulation suggests the following finite lemma:

```text
No primitive covering row can make the two-colour blocker projection saturate
all E_safe components without leaving an endpoint-rank <= 2 escape.
```

Equivalently, a counterexample would be a graph with no escape vertices: every
`E_safe` component would be `dead_both`, and the branch-coloured blocker
projection would cover the whole component set.  HYP-3451 says the current
audited bank always has a low-rank escape, and the tight AP-with-`84m` rows are
the right place to prove the saturation obstruction first.

The proof should try one of three routes:

```text
Menger cut: find a bounded blocker cut whose removal exposes an escape gap.
Green current: assign signed current to branch blockers and prove nonzero
  boundary at the four AP-tail escape components.
Algebraic connectivity: show any connected dead-cover projection with the AP
  tail blocker alphabet still has a forced endpoint-rank <= 2 boundary.
```

## Compression Guardrail

The blocker entropy range is a warning.  Compressing the dead-cover graph to
the number of dead components or total blocker mass loses which branch-coloured
blockers tie the components together.  The retained payload is:

```text
component id + branch-coloured minimal blockers + projection adjacency +
escape endpoint rank.
```

This is the graph-theoretic version of the HYP-3434 overlap-tax warning.

## Tournament Analysis

Tournament vertices are graph proof carriers, not runners or raw interval
counts.

Pairwise observable:

```text
predicate retention, cut payload, conductance route, entropy firewall
```

Fingerprint:

```text
score_hist={19:1, 25:1, 50:1, 51:1, 55:1, 56:1, 57:1}
directed_3cycles=0
hamiltonian_path =
  component_blocker_projection_graph
  -> green_current_escape_certificate
  -> menger_cut_saturation_obstruction
  -> algebraic_connectivity_router
  -> blocker_entropy_firewall
  -> raw_dead_component_count
  -> raw_dead_mass_scalar
```

## Next Hook

Prove the AP-with-`84` base case as a graph theorem: its `22` dead components
form one blocker-projection component, yet the row has four endpoint-rank `2`
escapes.  Then lift to the AP-with-`84m` tail rows and finally to arbitrary
primitive rows by treating `random_covering_082` as the high-rank but
low-danger negative control.

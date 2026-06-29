---
id: HYP-3450
title: LRC14 component-cover obstruction extractor
status: EVIDENCE / exact component obstruction audit; not an LRC14 proof
source: codex-2026-06-28 continuation of HYP-3435/HYP-3434/HYP-3429 after branch-cover and overlap-tax prompts
tangent: T1410
technique: LTI-410
tournament_technique: LTT-310
script: 04-computation/lrc14_component_cover_obstruction_extractor_codex_20260628.py
result: 05-knowledge/results/lrc14_component_cover_obstruction_extractor_codex_20260628.out
reflection: 07-reflections/lrc14-component-cover-obstruction-extractor-codex-20260628.md
related:
  - HYP-3435
  - HYP-3434
  - HYP-3433
  - HYP-3432
  - HYP-3431
  - HYP-3430
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3417
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3450: LRC14 Component-Cover Obstruction Extractor

## Claim

HYP-3434 identified the exact correction term lost by scalarizing the
one-branch cover:

```text
branch0_mass = naive_slack + overlap_tax.
```

HYP-3435 then kept both branch cells and active endpoint gates.  HYP-3450
pushes those two observations to the component level.  For each component `J`
of `E_safe`, it records whether branch `0` survives, branch `1` survives, both
survive, or both die.  When a branch dies on `J`, the script extracts the
smallest odd-bad interval subcover of `J`.  When a branch survives, it records
the survivor endpoint rank and labels.

This changes the theorem target.  A failed LRC14 branch-cover row would have to
put every `E_safe` component into `dead_both`, and every such component would
carry two exact minimal odd-bad covers plus the even endpoint gates that define
the component.  The proof route should classify those paired covers, not chase
a raw scalar average.

## Exact Readout

Script:

```text
04-computation/lrc14_component_cover_obstruction_extractor_codex_20260628.py
```

Stored result:

```text
05-knowledge/results/lrc14_component_cover_obstruction_extractor_codex_20260628.out
```

Aggregate audit on the HYP-3435 bank:

```text
rows_audited=135
components_total=17164
rows_with_branch_survivor=135/135
rows_with_endpoint_rank_le_2_survivor=135/135
rows_with_best_rank_gt_2=0
component_class_hist={
  both_alive: 6492,
  branch0_only: 3451,
  branch1_only: 3451,
  dead_both: 3770
}
dead_components_with_two_min_covers=3770/3770
dead_pair_cover_rank_hist={2:3492, 3:222, 4:44, 5:10, 6:2}
```

The tight canonical row remains proof-facing:

```text
covering_AP_with_84:
  components=26
  dead_both=22
  both_alive=4
  best_rank=2
  best_interval=[33/196,6/35]
  best_labels=L[E:84] R[B1:5]
  largest_dead_cover=[43/588,55/588]
  b0_cover=(1,(1,))
  b1_cover=(2,(11,13))
  pair_rank=3
```

The largest paired dead-cover rank in the audited bank is `6`, attained by
`random_covering_082`, where one dead component needs three branch-0 blockers
and three branch-1 blockers:

```text
b0_cover=(3,(77,89,105))
b1_cover=(3,(63,69,125))
```

## Proof Pull

The new finite lemma target is:

```text
Every primitive LRC14 covering row has an E_safe component whose branch-union
survivor has endpoint rank <= 2.
```

The contrapositive is more useful for proof search:

```text
If every E_safe component is dead_both, then the paired odd-bad minimal covers
and even endpoint gates form a forbidden finite obstruction.
```

This is a natural place for Menger cuts, Green currents, and conductance:
components are vertices, odd bad intervals are blockers, even endpoints are
gates, and a counterexample would be a paired cover flow saturating every
component.  HYP-3450 says the audited rows always leave a rank-`<=2` escape
component before that saturation can happen.

## Compression Guardrail

The information-theoretic lesson is explicit.  Replacing a union of odd-bad
intervals by total mass forgets not only overlaps, as in HYP-3434, but also the
component on which the cover is attempted.  The missing payload is:

```text
component id + two branch-cover statuses + minimal blocker covers + endpoint rank.
```

That payload is not commutative scalar data.  It is a local finite-ruler
certificate.

## Tournament Analysis

Tournament vertices are component proof obligations, not runners or arcs.

Pairwise observable:

```text
predicate retention, cover minimality, endpoint-gate rank,
two-adic branch payload, conductance-router value, scalar-risk control
```

Fingerprint:

```text
score_hist={12:1, 23:1, 50:1, 51:1, 55:2, 58:1}
directed_3cycles=0
hamiltonian_path =
  minimal_component_obstruction_extractor
  -> endpoint_gate_rank_lemma
  -> paired_odd_bad_cover_certificate
  -> branch_cover_sensitivity_ledger
  -> green_current_conductance_router
  -> raw_branch_union_measure
  -> scalar_harmonic_tail_slogan
```

## Next Hook

Prove or refute:

```text
paired dead-cover rank cannot saturate all E_safe components.
```

Start from `random_covering_082`, the current largest dead-pair rank example,
and from the canonical row `{1..11,13,84}`.  Try to build a bipartite
component/blocker incidence graph, then look for a bounded Menger cut or a
dual Green-current certificate that forces at least one rank-`<=2` survivor.

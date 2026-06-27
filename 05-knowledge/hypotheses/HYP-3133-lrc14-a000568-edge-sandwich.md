---
id: HYP-3133
title: LRC14 A000568 edge-witness sandwich
status: EVIDENCE / executable quotient scout; not a proof
source: codex-2026-06-27-A000568
tangent: T1200
technique: LTI-261
tournament_technique: LTT-159
script: 04-computation/lrc14_a000568_edge_sandwich_codex_20260627.py
result: 05-knowledge/results/lrc14_a000568_edge_sandwich_codex_20260627.out
related:
  - HYP-3131
  - HYP-3129
  - HYP-3128
  - HYP-3127
  - HYP-3125
  - HYP-3124
  - HYP-3123
  - HYP-3121
  - HYP-2968
  - OPEN-Q-108
---

# HYP-3133: LRC14 A000568 Edge-Witness Sandwich

## Claim

The user's observation

```text
10 < 12 < 16
20 < 56 < 80
```

is not only sequence numerology.  It is evidence for a useful middle quotient
between HYP-3124's four-sector edge word and its paired tail/tip child deck.

For an edge witness on `m` tournament vertices:

```text
sector_word_count(m) < A000568(m+1) < paired_child_signature_count(m)
```

holds in the first nontrivial tested cases `m=4,5,6`:

```text
m=4: 10 < 12 < 16
m=5: 20 < 56 < 80
m=6: 35 < 456 < 632
```

Here A000568 counts unlabeled tournaments one vertex higher.  The proposed
interpretation is:

```text
edge_extension_sandwich_certificate =
  four_sector_tetrahedral_word
  + A000568_unlabeled_extension_shadow
  + paired_tail_tip_child_deck
  + SPEC_resonance_lattice_sidecar_or_named_debt
```

The sector word is the equinumerosity layer.  The A000568 class is the
equidistribution/free-extension layer.  The paired child deck is the
equidecomposability layer.

## Evidence

The script
`04-computation/lrc14_a000568_edge_sandwich_codex_20260627.py` computes exact
small edge-witness counts through `m=6`, using the small A000568 values for the
unlabeled tournament middle quotient.  Output is stored at
`05-knowledge/results/lrc14_a000568_edge_sandwich_codex_20260627.out`.

The exact readout:

| edge `m` | sector words | A000568(`m+1`) | paired child signatures | status |
|---:|---:|---:|---:|---|
| 2 | 1 | 2 | 1 | small boundary |
| 3 | 4 | 4 | 4 | equality |
| 4 | 10 | 12 | 16 | proper sandwich |
| 5 | 20 | 56 | 80 | proper sandwich |
| 6 | 35 | 456 | 632 | proper sandwich |

The `m=6` computation touches all labelled tournaments and directed edge
instances at that size.  The largest `m=6` sector refinement, the balanced
sector `(1,1,1,1)`, has `47` paired child signatures inside one four-sector
word, so the child-pair layer is not a minor decoration: it is a real
equidecomposition refinement.

## LRC Use

HYP-3129 has already shifted the multi-far floor away from EH/BV and into an
elementary resonance-lattice certificate: exact finite low frequencies plus a
Parseval-tail bound.  HYP-3133 gives a finite quotient stratifier for the
remaining constant chase.  Instead of sorting hard rows by a scalar `Rprime` or
raw sector word, sort them by:

```text
sector word -> A000568 extension shadow -> paired endpoint children.
```

This should be tested on the HYP-3129 worst finite SPEC rows.  A bad row should
either improve under one endpoint deletion, land in a named resonance-lattice
debt class, or acquire an Asano/Lee-Yang/phi4/Cech/H7 terminal sidecar.

HYP-3128 is the guardrail: naive Asano over the crowded `R` tail fails.  The
A000568 middle shadow is therefore not a zero-free proof engine.  It is a
controlled-forgetting diagnostic between sector equinumerosity and full
tail/tip equidecomposition.

## Post-Rebase Connection With HYP-3131

After this scout was rebased, incoming HYP-3131/S69 independently recorded an
A000568 warm-up window: `C(n,3) <= A000568(n) <= 2(n-1)!/3` holds for
`n=4..7` and breaks at `n=8`.  That sharpens the status of HYP-3133.  The
sandwich should not be promoted to an all-`n` recurrence.  It is a finite
apex-7/LRC14 tameness diagnostic, exactly the range where HYP-3125/HYP-3129
need row stratifiers for the `r=2..6` few-apex floor.

## Tournament Analysis

Tournament vertices are proof carriers and quotient shadows, not runners or
scalar sequences.  The pairwise observable compares retained proof payload
over LRC predicate retention, edge-sandwich explanation power, tail/tip
recursion, equidecomposition power, equinumerosity guard, SPEC resonance
bridge, finite constant chase, formalization readiness, and scalar guardrails.

The scout tournament has:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

The selected path is:

```text
edge_floor_two_ended_packet
-> SPEC_resonance_lattice_certificate
-> A000568_unlabeled_extension_shadow
-> paired_tail_tip_child_deck
-> closed_form_constant_chase
-> Asano_tip_tail_legality_guard
-> even_graph_cycle_space_parity
-> four_sector_tetrahedral_word
-> raw_sequence_numerology
```

## Next Test

Attach `a000568_extension_shadow` to HYP-3125 edge-floor packets and then rerun
the HYP-3129 finite low-frequency row ledger.  The target is a monotonicity
dichotomy:

1. the middle shadow already predicts a positive SPEC floor;
2. one of the paired endpoint deletions improves the floor;
3. or the row is a finite named resonance-lattice debt.

That is the nearest proof-facing use of the 12/56 observation.

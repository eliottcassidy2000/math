---
id: HYP-3134
title: A000568 sits inside the edge-witness envelope as a global-consistency quotient
status: EVIDENCE / executable edge-envelope scout; not a proof
source: codex-2026-06-27-S272
tangent: T1201
technique: LTI-262
tournament_technique: LTT-160
script: 04-computation/a000568_edge_envelope_lrc14_codex_s272.py
result: 05-knowledge/results/a000568_edge_envelope_lrc14_codex_s272.out
related:
  - HYP-3136
  - HYP-3135
  - HYP-3133
  - HYP-3132
  - HYP-3131
  - HYP-3130
  - HYP-3129
  - HYP-3128
  - HYP-3127
  - HYP-3125
  - HYP-3124
  - HYP-3123
  - HYP-3121
  - HYP-3106
  - HYP-3054
  - HYP-3050
  - HYP-3049
  - HYP-3047
  - OPEN-Q-108
---

# HYP-3134: A000568 As Edge-Envelope Global Consistency

## Claim

The prompt observation is structurally real:

```text
10 < 12 < 16
20 < 56 < 80
```

In the HYP-3124 edge-witness census, `1,4,10,20,...` is the lower envelope of
raw four-sector size decks around a directed edge, while `1,4,16,80,...` is
the upper envelope after keeping both endpoint-deletion children.  A000568
appears one vertex later as a middle global quotient:

```text
sector_count(n) <= A000568(n+1) <= sector_plus_child_pair_count(n).
```

S272 extends the pattern one row farther:

```text
n=4: 10 < A000568(5)=12  < 16
n=5: 20 < A000568(6)=56  < 80
n=6: 35 < A000568(7)=456 < 632
```

The interpretation is controlled forgetting.  Four-sector decks are too
coarse: they count local equinumerosity/equidistribution around an edge but
forget the endpoint-deletion children.  Paired child signatures are safe but
overcomplete: they remember enough local data to avoid illegal quotienting.
A000568 sits between them because global tournament isomorphism classes are a
consistency quotient of local edge-witness data.

## Evidence From S272

The executable scout
`04-computation/a000568_edge_envelope_lrc14_codex_s272.py` computes the
edge-envelope counts directly through `n=6`, using the standard small A000568
values for `U(n+1)`.

| edge size n | sector deck count | sector+tail child | sector+tip child | sector+paired child | A000568(n+1) | wedge position |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | 1 | 2 | NA |
| 3 | 4 | 4 | 4 | 4 | 4 | NA |
| 4 | 10 | 14 | 14 | 16 | 12 | `1/3` |
| 5 | 20 | 52 | 52 | 80 | 56 | `3/5` |
| 6 | 35 | 274 | 274 | 632 | 456 | `0.705193` |

The lower count is exactly `C(n+1,3)`, the number of four-sector size
compositions for the `n-2` outside vertices.  The upper count is not a simple
composition count; it is the observed count after retaining the two endpoint
deletion children.  The wedge says the local HYP-3124 packet is intentionally
over-resolved before quotienting.

## LRC14 Transfer

The remaining LRC14 proof has the same shape.  HYP-3129 supplies an elementary
SPEC floor, but a scalar Fourier certificate is too coarse unless it is
attached to the packet whose predicate it preserves.  HYP-3125 supplies an
edge-floor packet, but the full local edge packet may be overcomplete unless
we name the global gluing rule.  The A000568 wedge suggests the missing
middle column:

```text
edge_floor_packet +=
  envelope_position
  + global_consistency_class
  + edge_child_gluing_status
  + resonance_lattice_class
  + SPEC_bound_status
```

In words: start with the safe HYP-3124 two-ended edge packet, then quotient
only by a named global-consistency rule.  For LRC14, that rule should be the
combined HYP-3129 resonance-lattice/SPEC certificate, HYP-3125 edge-floor
packet, and HYP-3123 normal-fan/Cech or observer-gluing legality exit.

Post-fetch integration: HYP-3133 supplies the direct A000568 edge-sandwich
scout, HYP-3136 assembles the multi-far floor into a finite constant chase,
HYP-3130 supplies the Gaussian/Beurling-Selberg minorant side of the same
tool stack, and HYP-3131 says far elements push miss-PGF zeros outward so the
multi-far part is not the main obstruction.  HYP-3132 names the bounded-core
k=8 De Moivre/phi4 tail node, and HYP-3135 keeps the corresponding middle
resolvent/SPEC payload.  HYP-3134 is compatible with all of these: it is not
a new analytic floor, but the quotient discipline for deciding when those
certificates may forget the HYP-3124 paired child payload.

## Tournament Analysis

Tournament vertices are proof carriers and quotient operators, not runners or
raw counts.  Pairwise observable compares A000568 wedge explanation, LRC
predicate retention, controlled forgetting, tail/tip payload, global gluing,
SPEC-floor transfer, finite checkability, and failure guardrails.

S272 fingerprint:

```text
score_hist={0:1,1:1,3:3,5:1,6:1,7:1,8:1}
directed_3cycles=1
scc_sizes=[3,1,1,1,1,1,1]
hamiltonian_path_count=3
selected_path =
  global_consistency_quotient
  -> edge_floor_packet_schema
  -> paired_tail_tip_child_envelope
  -> observer_extension_cut_payload
  -> A000568_middle_orbit_count
  -> resonance_lattice_SPEC_certificate
  -> asano_lee_yang_dichotomy
  -> raw_four_sector_composition
  -> raw_count_numerology
```

The selected path is useful because it puts the quotient rule before both
local over-refinement and raw counting.  Counts alone are not a proof; they
become proof-relevant only when they specify which local edge children glue to
one global packet without changing the LRC predicate.

## Next Theorem Target

For the `r=2..6` covering-floor rows, build a joined packet ledger:

```text
R-safe -> Q-safe edge packet
  + HYP-3124 paired child signatures
  + HYP-3129 resonance lattice/SPEC certificate
  + global_consistency_class
  + edge_child_gluing_status
  + terminal_exit_or_named_debt
```

Then prove that quotienting by `global_consistency_class` preserves the
positive floor certified by the SPEC bound.  This is the LRC analogue of
A000568 sitting between `20` and `80`: the final certificate should be neither
the raw scalar sector count nor the fully over-refined local child packet, but
the legal global quotient between them.

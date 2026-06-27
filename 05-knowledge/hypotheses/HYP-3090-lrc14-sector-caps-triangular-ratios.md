---
id: HYP-3090
title: LRC14 sector caps are triangular pair-avoidance ratios
status: VERIFIED EXACT cap skeleton / remaining deviation debt; not a proof of LRC14
source: kind-pasteur-2026-06-27-S31ag + codex-2026-06-27-postrebase integration
theorem: THM-576
script: 04-computation/lrc14_cap_is_pair_avoidance_kps.py
result: 05-knowledge/results/lrc14_cap_is_pair_avoidance_kps.out
related:
  - THM-576
  - HYP-3083
  - HYP-3085
  - HYP-3088
  - HYP-3089
  - HYP-3091
  - HYP-3092
  - HYP-3093
  - HYP-3094
  - OPEN-Q-108
---

# HYP-3090: LRC14 Sector Caps Are Triangular Pair-Avoidance Ratios

This hypothesis records the incoming S31ag cap result as a first-class
proof-interface object for the LRC14 finite-address spine.

The covering caps have the clean triangular form

```text
cap_k = C(k+1, 2) / C(14, 2) = C(k+1, 2) / 91
```

exactly for `k=10,11,12,13`:

```text
cap_10 = 55/91
cap_11 = 66/91
cap_12 = 78/91 = 6/7
cap_13 = 91/91 = 1
```

The two binding rows dip below the triangular skeleton:

```text
cap_9 = 1979/4004 = 45/91 - 1/4004
cap_8 = 2243/5880 = 36/91 - 1081/76440
```

THM-576 recasts the cap side as pairwise avoidance, verifies/proves the small
sector cases through `j<=3`, and exposes the next visible obstruction as an
order-3 break rather than another pairwise term.

Numbering note: THM-575 is the separate raw-time Conjecture 7.1 refutation.
The pairwise-avoidance cap theorem is THM-576.

The concrete small-part version is:

```text
min_{|P|=j} meas(lonely(P)) = C(14-j, 2) / C(14, 2)
```

for `j=1,2,3`, with minimizers `{1}`, `{1,13}`, and `{1,12,13}`.
The binding-row minimizers recorded by the exact scout are `{1,11,12,13}`
for `j=4` and `{1,5,7,8,9}` for `j=5`.

## Proof-Spine Meaning

This does not close LRC14.  It changes the shape of the remaining cap/moment
debt:

```text
clean pairwise cap skeleton
  + two deviation constants (k=8,9)
  + first higher-order break packets
  -> normalized arc floor / O2 discharge / named residual debt
```

So future finite-address ledgers should not store a generic "cap bound"
coordinate.  They should store:

```text
cap_ratio_or_deviation_status
pairwise_avoidance_sector_count
triangular_cap_value
deviation_constant
higher_order_break_status
handoff_to_normalized_arc_floor
```

HYP-3091, HYP-3092, and HYP-3093 are downstream of this carrier in different
ways.  HYP-3092 refines the cap value itself as pair-normalized Pascal mass.
HYP-3091 says the cap packet should retain its equinum/equidecomp/equidist
fiber, especially the `D=41` bounded core and `1/lmax=V*` apex invariant.
HYP-3094 then decides whether a cap/deviation packet can feed
nested-refinement covering discharge or must be rerouted to
cross-handoff/K33/THM-572 debt.  The key guardrail remains the same
controlled-forgetting rule: a cap scalar is proof-usable only if the pairwise
skeleton, deviation, higher-order break address, three-sameness payload, and
HYP-3093 forgetting-cost tuple that produced it are retained.

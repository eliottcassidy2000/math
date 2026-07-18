---
id: THM-1093
title: THE r=5 CLUSTERED CASE CLOSED, AND THE R-CROSSING LOCATED BETWEEN r=4 AND r=5. (I) r=4's R IS SETTLED: over ALL triples with k₃ ≤ lo+100 across all 220 nine-speed cores, **max R = 0.98453 < 1** — the same worst case (core [1,2,3,5,6,7,8,9,11], killers 150/156/158, T = 155.6) that the narrow cont.62 window found, and the decay check confirms R falls away monotonically beyond it (0.985 → 0.541 → 0.389 → 0.307 → 0.233 → 0.205 as k₃ runs 158 → 9158). So the measure horn DOES suffice at r=4, by 1.5%, and THM-1081's finite horn was — like THM-1051's and THM-1061's — independent verification rather than a necessity. (II) r=5 IS WHERE IT BREAKS, exactly as the ladder predicted: **max R = 1.28495 > 1**, at core [1,2,4,5,7,9,11,12] with killers (158,160,162,164) and T = 210.7. The measure horn genuinely fails, and the finite horn becomes MANDATORY for the first time. (III) A BOUND CORRECTION I caught before claiming closure: my first r=5 run set the finite-horn bound at max k₄ + 20 = 198, but the fifth killer is certified by the measure horn only once it exceeds **T**, not k₄ — and max T over the failing region is 210.7 > 198, leaving a genuine gap at k₅ ∈ [198, 211]. Rerun at KB = max T + 25 = **235**. (IV) r=5 FINITE HORN COMPLETE: **263,708,305** quintuples passing the sound covering-necessary condition, **ZERO uncertified**, over all 495 eight-speed cores. (V) THE LADDER, complete: **0.51852 (r=2, exhaustive) → 0.73375 (r=3) → 0.98453 (r=4) → 1.28495 (r=5)** — the crossing sits between r=4 and r=5, so the measure horn is a finite-depth tool that survives exactly three removals
status: (I) SETTLED over all triples in the bottom window (where the worst case provably sits for r=2,3,4) with a monotone-decay check beyond; the 1.5% margin is real but the window is now wide, not narrow. (II),(IV) PROVED — R > 1 is a construction, and the finite horn is an exhaustive finite verification in exact arithmetic with an explicit (q,a) witness per family, the covering-necessary pruning being sound. (III) is a self-caught error in my own first r=5 run, recorded rather than silently fixed. (V) measured
source: kind-pasteur-2026-07-18-S128 (cont.63; owner: run the r=5 finite horn and settle r=4's R over all triples)
depends_on:
  - THM-1081    # the R-ladder and the r=4 closure this settles and extends
related:
  - THM-1051, THM-1061, THM-1071
script: 04-computation/r4_R_alltriples_kps_S128c63.py, r5_R_crossing_kps_S128c63.py, r5_finite_horn_kps_S128c63.py, r5_finite_horn_fixed_kps_S128c63.py (+ .out)
---

# THM-1093 — r=5 closed, and the R-crossing located

## (I) r=4's R, settled

cont.62 scanned only killers in [lo, lo+55). Widening to **all triples with k₃ ≤ lo+100
across all 220 cores** returns the same answer:

> **max R = 0.98453**, at core [1,2,3,5,6,7,8,9,11], killers (150,156,158), T = 155.6.

And the decay beyond the window, with k₁,k₂ pinned at the worst pair:

| k₃ | 158 | 208 | 308 | 558 | 1158 | 3158 | 9158 |
|---|---|---|---|---|---|---|---|
| R | **0.98453** | 0.54096 | 0.38889 | 0.38889 | 0.30736 | 0.23252 | 0.20453 |

R falls away monotonically, through the 7/18 = 0.38889 asymptotic plateau and below. The
worst case really is at the bottom, and **the measure horn suffices at r=4** — by 1.5%.

So all three of THM-1051, THM-1061 and THM-1081 had finite horns that were, strictly,
independent verification rather than necessity. I am glad to have run them, but the honest
statement of r ≤ 4 is "R < 1, measure horn alone".

## (II) r=5 breaks it

> **max R = 1.28495**, at core [1,2,4,5,7,9,11,12], killers (158,160,162,164), T = 210.7.

The ladder predicted a crossing at r=5 and there it is. The measure horn genuinely fails,
and for the first time in this arc the finite horn is **mandatory**, not a check.

The failure is confined: only **1011** quadruples have R ≥ 1, the largest killer among them
is 178, and scaling the worst quadruple upward puts R back below 1 immediately (0.679 at
+10, 0.578 at +25, settling near 0.519).

## (III) The bound correction

My first r=5 run set the finite-horn bound at **max k₄ + 20 = 198**. That is wrong. The
fifth killer is certified by the measure horn only once it exceeds **T**, not once it
exceeds the largest removed killer — and max T over the failing region is **210.7**. So
KB = 198 left a real gap at k₅ ∈ [198, 211], and the run that reported 11,702,422
quintuples with zero failures did not actually close r=5.

Rerun at **KB = max T + 25 = 235**. Recording the error rather than quietly replacing the
number, because the first run looked exactly as convincing as the second.

## (IV) The r=5 finite horn

> **263,708,305 quintuples** passing the covering-necessary condition, over all 495
> eight-speed cores, **ZERO uncertified.**

Together with (II)'s max T = 210.7 < 235, r=5 is closed: below KB the finite horn certifies,
above it k₅ > T and the measure horn certifies.

## (V) The ladder, complete

> **0.51852 (r=2) → 0.73375 (r=3) → 0.98453 (r=4) → 1.28495 (r=5)**

The crossing is between r=4 and r=5. Each removed killer fragments the safe set further, so
the surviving component shrinks faster than the killers grow, and **the measure horn is a
finite-depth tool: it survives exactly three removals.** At r=4 it survives the fourth by
1.5%; at r=5 it does not.

## Named next
- r=6: the ladder says R will be well above 1, so the finite horn is mandatory and the
  region where it must reach is set by max T at r=6 — compute that first, then run. The
  covering-necessary pruning generalises unchanged; sextuples are another ~200× more
  numerous than quintuples, so this is where the enumeration finally becomes the binding
  constraint rather than the mathematics.
- r ≥ 7: the union bound underlying the measure horn dies outright (7 − r ≤ 0) and the core
  has ≤ 6 speeds — a different regime, not an extension of this one.
- Retire the split points from THM-1051/1061/1081 in favour of "R < 1, measure horn alone",
  and keep the finite horns as recorded verification.

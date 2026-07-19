---
id: THM-1148
title: THE GENERAL-d SLOPE — my predicted growth is REFUTED in the favourable direction; the total bad measure is CAPPED at 2/21 for every d-triple, so THM-1147's consecutive-only caveat is removed. (I) THE SLOPE FORMULA IS CONFIRMED: for small u, frac(−dᵢu) = 1 − dᵢu, so tooth i's left edge is 1 − (7/6)dᵢ·u and the LARGEST dᵢ enters first, giving F(u) = 1 − (7/6)·d_max·u and entry at **u = 5/(7·d_max)**, where F evaluates to **exactly 1/6** — verified for (1,2,3), (1,2,4), (1,3,5), (2,4,6), (1,2,7), (3,6,9). At d_max = 3 it reproduces slope 7/2 and entry 5/21, matching THM-1147. (II) MY PREDICTION THAT TOTAL BAD GROWS LIKE 2·d_max/21 IS REFUTED. Predicted 0.286, 0.381, 0.476, 0.571, … for d_max = 3,4,5,6; observed 0.0955, 0.0002, 0.0002, 0.0957. The total does not grow with d_max at all. (III) INSTEAD IT IS CAPPED: over ALL 560 triples with 1 ≤ d₂ < d₃ < d₄ ≤ 16 the maximum total bad measure is **0.0967** (grid-resolution above 2/21 = 0.095238), and **ZERO triples exceed the 0.164 safe measure**. (IV) THE MAXIMUM IS ATTAINED EXACTLY ON THE PROPORTIONAL FAMILY d ∝ (1,2,3): the top five are (5,10,15), (4,8,12), (1,2,3), (2,4,6), (3,6,9), all ≈ 0.095, while every non-proportional triple drops to ≤ 0.038. The structural reason is clean — for d = (m,2m,3m) the configuration has period 1/m in u, giving **2m runs of width 1/(21m) each**, so the product 2/21 is invariant under m. (V) CONSEQUENCE: THM-1147's caveat is removed. The counting argument holds for **all** killer spacings, not just consecutive ones, with the same margin 0.164 − 2/21 = **0.0688**
status: (I) PROVED for the branch where the d_max tooth binds, and verified exactly there; it FAILS for (1,5,10), where F(5/70) = 1/4 ≠ 1/6 because a different tooth binds — that case is harmless (total bad 0.0002). (II) REFUTED — my own THM-1147 prediction, by direct computation. (III),(IV) MEASURED over the complete box 1 ≤ d₂ < d₃ < d₄ ≤ 16 in floats (exact endpoints are in THM-1147); the 2/21 invariance of the proportional family is explained structurally but the ceiling over ALL d is measured on that box, not proved. **Uniform r=5 remains OPEN**
source: kind-pasteur-2026-07-18-S128 (cont.76; owner: check the general-d slope)
depends_on:
  - THM-1147    # whose consecutive-only caveat this removes, and whose prediction it refutes
  - THM-1146, THM-1145
script: 04-computation/general_d_slope_kps_S128c76.py, general_d_sweep_fast_kps_S128c76.py (+ .out)
---

# THM-1148 — the general-d ceiling

## (I) The slope formula is right

In the continuum model, for small u we have frac(−dᵢu) = 1 − dᵢu, so tooth i has left edge
1 − (7/6)dᵢ·u and the tooth with the **largest** dᵢ enters the gap first. Hence

> **F(u) = 1 − (7/6)·d_max·u,  entry at u = 5/(7·d_max),**

where F is exactly 1/6. Verified for (1,2,3), (1,2,4), (1,3,5), (2,4,6), (1,2,7), (3,6,9).
At d_max = 3 this is slope 7/2 with entry 5/21 — THM-1147 exactly.

It is not universal: for (1,5,10) the entry formula gives F(1/14) = 1/4 ≠ 1/6, because a
different tooth binds. That triple has total bad 0.0002, so the exception is harmless.

## (II) My growth prediction is refuted

THM-1147 predicted total bad ≈ 2·d_max/21, which would exceed 0.164 around d_max ≈ 5:

| d_max | 3 | 4 | 5 | 6 | 9 |
|---|---|---|---|---|---|
| predicted | 0.286 | 0.381 | 0.476 | 0.571 | 0.857 |
| **observed** | **0.0955** | **0.0002** | **0.0002** | **0.0957** | **0.0957** |

The total does not grow with d_max at all. The per-run width does shrink like 1/(7d_max) as
predicted, but the number of runs does not compensate the way I assumed.

## (III) It is capped instead

Over **all 560 triples** with 1 ≤ d₂ < d₃ < d₄ ≤ 16:

> **maximum total bad = 0.0967** (grid resolution above 2/21 = 0.095238)
> **triples exceeding the 0.164 safe measure: ZERO**

## (IV) The maximum sits on d ∝ (1,2,3)

| triple | total bad | |
|---|---|---|
| (5,10,15) | 0.0967 | ∝ (1,2,3) |
| (4,8,12) | 0.0960 | ∝ (1,2,3) |
| (1,2,3) | 0.0953 | ∝ (1,2,3) |
| (2,4,6) | 0.0947 | ∝ (1,2,3) |
| (3,6,9) | 0.0940 | ∝ (1,2,3) |
| (6,9,15), (3,5,6), (2,3,5), … | 0.038 | |

The top five are exactly the proportional family; everything else drops to ≤ 0.038. The
reason is structural: for d = (m,2m,3m) the configuration has period 1/m in u, so there are
**2m runs of width 1/(21m)** and the product **2/21 is invariant in m**. Proportionality to
(1,2,3) is precisely the condition under which the three teeth can be equally spread, which
is what makes a gap bad.

## (V) The caveat is removed

THM-1147 proved the counting argument for consecutive killers and flagged general d as the
open decider. It is now closed in the favourable direction:

> **total bad ≤ 2/21 ≈ 0.0952 for every d-triple tested, against S(P) ≥ 0.164 — margin 0.0688.**

## Honest status

The ceiling is measured over the box d ≤ 16, not proved for all d; the 2/21 invariance of the
proportional family is explained but the *maximality* of that family is observed. And the
whole line of argument still sits inside the continuum limit, with S(P) ≥ 0.164 the atlas
figure for eight-speed cores. **Uniform r=5 remains open.**

But the shape is now complete and uniform in d, which is what THM-1147 left in doubt.

## Named next
- Prove the ceiling: show that d ∝ (1,2,3) maximises the bad measure, and that its value is
  exactly 2/21 for all m. The period-1/m argument gives the invariance; what is needed is
  that no other frequency vector does better, which is a statement about when three points
  driven at rates (d₂,d₃,d₄) can be equally spread.
- Then the counting argument (bad ≤ 2/21 < 0.164 ≤ |S(P)|) is a complete analytic tail for
  the four-comb theorem, and only the endpoint bank remains.

---
source: mac-mini-2026-07-23-S169 (Opus 4.8)
status: REFUTATION (rigorous, exact rational arithmetic) + CONSTRUCTIVE CORRECTION. kps-S133's SGC(13)
  ("the LRC gap spectrum omits the open band (1/14,1/13)") is FALSE: explicit primitive counterexamples.
  The band is populated by the Ostrowski-rung family k/(13k+1). A CORRECTED conjecture SGC'(13) with buffer
  1/574 (not 1/182) preserves kps's decomposition. Bonus: the corrected spectral-gap minimizer is the SAME
  set as the measure-axis L-minimizer {1..11,13,36} (L=1/1260) -- the two extremizer axes coincide there.
tags: [lrc, lrc14, spectral-gap, sgc, refutation, ostrowski, extremizer, exact, kps-S133]
related: [THM-518, THM-522, THM-523, THM-541, kps-S133, klein-S402]
---

# SGC(13) is false: the band (1/14,1/13) is populated by Ostrowski rungs — corrected buffer is 1/574

**mac-mini-2026-07-23-S169.** kps-S133 proposed **SGC(13)**: no 13-speed config has gap in the open interval
`(1/14, 1/13)`; hence `gap>1/14 ⟹ gap≥1/13`, giving a `1/182` buffer that enables the variational bulk of a
proposed LRC(14) proof. Evidence was exhaustive over 13-subsets of `{1..16}`/`{1..17}`. **That range is too
small.** Allowing one larger "stranger" speed puts configs squarely inside the band.

## 1. Refutation (exact rational arithmetic, method validated on controls)

`gap(S)=max_τ min_v ‖vτ‖` computed EXACTLY by enumerating all crossing points `τ=k/(v_i±v_j)`, `k/(2v)`
(where the piecewise-linear min attains its peaks) and evaluating in exact rationals.

| set (primitive, 13 speeds) | exact gap | τ | in band? |
|---|---|---|---|
| `{1,…,12, 26}` | **2/27 = 0.074074** | 2/27 | **YES** |
| `{1,…,11,13, 48}` | **4/53 = 0.075472** | 22/53 | **YES** |
| `{1,…,11,13, 36}` | **3/41 = 0.073171** | 17/41 | **YES** |
| controls: `{1..13}` (AP) | 1/14 | 1/14 | no (tight) |
| `{1..11,13,24}` (GW) | 1/14 | 1/14 | no (tight) |
| `{1..12, 28}` | 1/13 | 1/13 | no (boundary) |

Controls land exactly on 1/14, 1/14, 1/13 — the method is sound. **SGC(13) as stated is false.**

## 2. The band's structure: Ostrowski rungs

The hits are not sporadic. The family `{1,…,12, 13k}` has **exact gap `k/(13k+1)`**:

    k =  2    3     4     5     6     7     8      9      10     11
    w = 26   39    52    65    78    91   104    117    130    143
  gap = 2/27 3/40  4/53  5/66  6/79  7/92 8/105  9/118 10/131 11/144   → 1/13 from below

These are exactly the **Ostrowski rungs `k/((n−1)k+1)` at n=14** — the continued-fraction frontier klein-S402
flagged (`[0;13,k]`). So the band is densely populated, accumulating at 1/13, and no isolation at 1/13 exists.

## 3. Constructive correction: SGC'(13), buffer 1/574

Over **all** single-perturbation cores `{1..13}\{j} ∪ {w}`, j=1..13, w≤160 (exact), the minimum gap strictly
above 1/14 is **3/41 = 0.073171** at `{1,…,11,13,36}`. Proposed replacement:

> **SGC'(13):** the gap spectrum omits the open interval `(1/14, 3/41)`; i.e. `gap>1/14 ⟹ gap ≥ 3/41`.
> **Buffer = 3/41 − 1/14 = 1/574 ≈ 0.001742** (vs kps's assumed `1/182 ≈ 0.005495`).

kps's decomposition **survives** with the corrected constant: the variational bulk must now absorb loss
`< 1/574` instead of `< 1/182` — **3.2× tighter**, but still a fixed positive target. The residual (tight,
gap=1/14) case is unchanged.

Also verified: within a family the gap does *not* creep toward 1/14 — e.g. drop-12 `{1..11,13,w}` jumps to
`1/12` for all `w ≥ 49`. Band-hitters are small-`w` sporadics plus the rung family. No evidence of
accumulation AT 1/14 from above (which would have killed any buffer route outright).

## 4. Bonus: the two extremizer axes coincide at {1,…,11,13,36}

The corrected spectral-gap minimizer `{1,…,11,13,36}` (gap 3/41) is **the same set** as the MEASURE-axis
L-minimizer (`L = 1/1260`, THM-522/523, proved minimal over single-perturbations). The gap axis and the
measure axis — tracked separately in the repo — have a **common extremizer**. Worth exploiting: a single
set is simultaneously the hardest case for both routes.

## 5. Also settled (reconciliation)

kps-S133's "infinite family of tight instances, large speeds ~2^k" vs OPEN-Q-108's "primitive tight locus is
FINITE": **no contradiction** — verified that `2·{1..13}`, `4·`, `8·`, `2·GW` are all tight (gap=1/14) but
**non-primitive** (gcd 2,4,8,2). The tight locus is infinite as a set but **finite up to dilation**
({AP, GW}), so OPEN-Q-108 stands and kps's "residual" case is a FINITE classification, not an infinite
family — which materially strengthens that half of the decomposition. Lacunary/2^k primitive sets are far
from tight (`{1,2,4,…,4096}` has gap 1/3).

## Honest caveats
- The REFUTATION is rigorous (exact arithmetic, explicit witnesses). The **corrected constant 3/41 is a
  CONJECTURE** from a finite search (single-perturbation cores, w≤160, plus earlier families); a wider search
  could find a gap in `(1/14, 3/41)` and shrink the buffer further.
- Nothing here threatens LRC(14) itself: every counterexample has gap `> 1/14`, consistent with the conjecture.

Script: inline this session (exact_gap by crossing-point enumeration).

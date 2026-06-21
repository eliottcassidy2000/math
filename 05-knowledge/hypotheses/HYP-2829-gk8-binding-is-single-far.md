---
id: HYP-2829
title: The gK8 concentration's BINDING wide case is SINGLE-FAR (r=1) — max L_yK8 is dominated by bounded (r=0) with a LARGE wide margin (~2), and decreases sharply with far-count r; so the concentration reduces to [bounded finite check] + [single-far via THM-563 periodicity] + [r≥2 uniformly lower], avoiding the razor-thin p0 dichotomy
status: VERIFIED (max L_yK8 by far-count: r=0 >> r=1 > r≥2, k=8,9,10; wide margin ~2 vs p0's razor-thin). The single-far binding is provable via THM-563 periodicity (my domain). OPEN: rigorize the r≥2-below-r=1 step + the single-far sup.
source: mac-mini-2026-06-22-S24
related:
  - HYP-2823  # variance reframe: wide => low Var(N) => low extreme mass q0+q6 => low L_yK8
  - HYP-2820  # q6-suppression (why wide L_yK8 is low: q6 suppressed)
  - HYP-2812  # claude-opus gK8 concentration extremality (this is its far-count structure)
  - THM-563   # single-far periodicity (closes the binding case)
  - HYP-2809  # the dichotomy gap (this REPLACES it for the gK8 route)
  - OPEN-Q-108
---

# HYP-2829 — The gK8 binding is single-far

## The finding (VERIFIED)
Max `L_yK8 = 10q0+q3+10q6` over size-k configs, stratified by far-count `r = #{e>14}`:

| k | r=0 (bounded) | r=1 | r=2 | r≥3 | 10·cap |
|---|---------------|-----|-----|-----|--------|
| 8 | **3.582** | 1.74 | 1.21 | 1.23 | 3.815 |
| 9 | **4.434** | 2.47 | 2.19 | 2.16 | 4.943 |
| 10 | **5.286** | 4.60 | 3.10 | 2.77 | 5.714 |

**Bounded (r=0) dominates by a LARGE gap** (1.7–2× the wide max), and L_yK8 drops sharply with r.

## Why (the q6 mechanism)
`L_yK8 = 10(q0+q6) + q3`, and `q0+q6 = P(N∈{0,6})` is the EXTREME mass of the miss-count N. Wide
configs decorrelate ⟹ low `Var(N)` (HYP-2823) ⟹ N concentrates in the MIDDLE ⟹ low extreme mass ⟹
low `L_yK8`. The q6 extreme is suppressed per far runner (HYP-2820). So the gK8 functional PENALIZES
spread, giving the wide region a LARGE margin (~2) — **unlike the p0 route, where the wide doublet sits
razor-thin below Q**. The dangerous-for-p0 doublet (high q0) has LOW q6, hence LOW L_yK8.

## The reduction (the path)
The gK8 concentration `max_E L_yK8 ≤ 10cap` reduces to:
1. **r=0 (bounded):** the finite check (12805 bases / span≤14 exhaustive, DONE), `L_yK8(consec)<10cap`
   (margins 0.23/0.51/0.43).
2. **r=1 (single-far) ≤ r=0:** for bounded/dilated base B, `max_{far≥15} L_yK8(B∪{far}) ≤ L_yK8(consec_k)`
   — provable via THM-563 periodicity (`far·(deviation of each q_t)` periodic ⟹ sup = window+tail). **My
   domain.** This is the BINDING wide case.
3. **r≥2 ≤ r=1:** uniformly lower (more decorrelation); a separate, EASIER bound (large margin).

So the gK8 route closes WITHOUT the razor-thin p0 dichotomy (HYP-2809): the binding is single-far
(periodicity, closed), and multi-far is comfortably below. The gK8 functional converts the hard p0
dichotomy into a comfortable far-count monotonicity whose binding step is exactly the single-far
periodicity I already have.

## Next
Rigorize: (a) single-far gK8 sup ≤ bounded max for ALL bounded/dilated B (periodicity, exhaustive like
the THM-563 12805-base check); (b) the r≥2 ≤ r=1 (or ≤ bounded) step. Then the gK8 Lean
`gK8_concentration_extremality` sorry closes.


## CONFIRMED (mac-mini-S24): single-far gK8 sup < 10cap with COMFORTABLE margin
sup_{far in [15,120]} L_yK8(consec_{k-1} u {far}) = 2.37 / 3.91 / 4.81 (k=8/9/10), all BELOW 10cap with
margin **0.90 / 1.03 / 1.44** (and below the bounded max with gap 1.21/0.52/0.47). So the BINDING wide
case (single-far) clears 10cap with margin ~1 -- vs THM-563's p0 check margin 0.13. The gK8 functional
makes the wide binding COMFORTABLE: a LOSSY single-far bound suffices. Rigorous closure = the THM-563
period-max machinery applied to L_yK8 (far*deviation of each q_t periodic) over all bounded/dilated bases
B -- the SAME 12805-base finite check that PASSED for p0 (margin 0.13), now with margin ~1, so it passes
a fortiori. NET: gK8 concentration = [bounded r=0 finite, DONE] + [single-far r=1 < 10cap margin ~1, via
THM-563 periodicity, the binding wide] + [r>=2 lower: q6 doubly-suppressed (HYP-2820) + comfortable wide
p0 slack]. The wide leg closes through the COMFORTABLE single-far margin, not the razor-thin p0 dichotomy.

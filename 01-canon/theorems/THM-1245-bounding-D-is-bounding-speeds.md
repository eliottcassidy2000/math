---
id: THM-1245
title: BOUNDING D REDUCES EXACTLY TO BOUNDING THE SPEEDS — since D = M·s, and 8.5M FAMILIES OF EVIDENCE — eliminating t* at the straddling pair gives g·(vᵢ+vⱼ) = D, i.e. **D = M·s exactly**; with M ∈ (1/14, 3/41) this makes D ≈ s/14, so bounding D is EQUIVALENT to bounding the active-pair sum s ≤ 2·max(V). Hence "bound D" is precisely "for PRIMITIVE families, does M near the floor force bounded speeds?" — and primitivity is essential, since dilating {1,…,13} by k holds M = 1/14 while sending max speed to 13k. Evidence: a pure-integer early-exit scan (gap > num/den ⟺ min(vp mod q, q−vp mod q)·den > num·q) covered ~8.5 MILLION random primitive families at speed caps 20–110; the only family found with M ≤ 3/41 was {1,…,13} itself, and at caps ≥ 30 there were ZERO. The three known near-floor families have max speeds 13, 24, 36 — all small, consistent with a bound but not proving one
status: D = M·s is exact algebra and the reduction to bounding speeds is exact. The 8.5M-family scan is SEARCH, and weak in a specific way I record: random sampling of 13 speeds from [1,N] essentially never produces the structured near-floor families (it missed both {1,…,11,13,24} and {1,…,11,13,36}), so its value is scale, not coverage. The bound on D is NOT proved and (1/14, 3/41) remains unsettled
source: opus-2026-07-19-S398 (owner: work the bound on D)
depends_on: [THM-1240 (which named bounding D as the obstruction), THM-1235 (slack, the gap edge), THM-1205 (D = M·s in the (D,s) coordinates), THM-1050 (dilation — why primitivity is required)]
scripts: 04-computation/bound_D_opus_S398.py -> 05-knowledge/results/bound_D_opus_S398.out
---

# THM-1245 — what bounding D actually is

## The reduction

At the straddling active pair, vᵢt* = aᵢ + M and vⱼt* = aⱼ − M; eliminating
t* gives M(vᵢ+vⱼ) = vᵢaⱼ − vⱼaᵢ, i.e.

> **D = M·s**, exactly.

Verified on all three near-floor families. With M ∈ (1/14, 3/41) this forces
D ≈ s/14, and since s = vᵢ + vⱼ ≤ 2·max(V):

> **bounding D ⟺ bounding the active-pair sum ⟸ bounding max(V).**

So THM-1240's named obstruction becomes a sharper question:

> **For PRIMITIVE families, does M near the floor force bounded speeds?**

Primitivity is not a technicality — it is the whole content. Dilating
{1,…,13} by k holds M = 1/14 exactly while sending max(V) to 13k, so without
restricting to gcd = 1 no bound of any kind can exist (THM-1050).

## The evidence

The test is purely integral — gap at p/q exceeds num/den iff
min(vp mod q, q − vp mod q)·den > num·q — so the hot loop needs no rational
arithmetic. That plus early exit on the first failing speed made the scan
roughly a hundred times faster than the Fraction version:

| speed cap N | primitive families | with M ≤ 3/41 |
|---|---|---|
| 20 | 1,643,851 | 13 — all {1,…,13} |
| 30 | 1,727,608 | **0** |
| 45 | 1,754,179 | **0** |
| 70 | 1,710,444 | **0** |
| 110 | 1,639,810 | **0** |

About **8.5 million** primitive families; the only one reaching M ≤ 3/41 is
the classical extremal, and nothing at all at caps ≥ 30. Sanity checks pass:
the three known near-floor families do not beat 3/41, and {1,…,12,14}
(M = 1/13) does.

## The weakness I want on record

Random sampling of 13 speeds from [1,N] essentially never produces the
structured near-floor families — this scan missed **both** {1,…,11,13,24}
and {1,…,11,13,36}, which are known to sit at 1/14 and 3/41. So its
contribution is **scale, not coverage**, and the structured searches of
S396/S397 remain the stronger evidence. Eight million families sounds
decisive and is not.

What the scan does support is the *shape* of a bound: the three known
near-floor families have max speeds 13, 24 and 36, all small, and raising
the cap moves random families away from the floor rather than toward it.

## Status

The reduction is exact and, I think, the useful output: "bound D" and "M near
the floor forces bounded primitive speeds" are the same statement, and the
latter is the one to attack. The bound itself is **not proved**, and
(1/14, 3/41) remains unsettled.

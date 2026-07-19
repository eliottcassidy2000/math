---
id: THM-1210
title: D = 1 IS EXACTLY THE CLASSICAL SIEVE — AND THE HARD STRATUM IS SAVED BY HIGH-D PAIRS, NOT BY SIEVING — the determinant D = vᵢaⱼ − vⱼaᵢ = 1 forces gcd(vᵢ,vⱼ) = 1, and the canonical instance (vᵢ,vⱼ) = (1, s−1) puts the maximizer at t* = 1/s with both runners at distance exactly 1/s and every other runner ≥ 1/s away **iff s divides no speed** — precisely the classical sieve at modulus s, delivering gap exactly 1/s. So D = 1 with s = 14 IS the tight family (both known ones have active pair (1,13), D = 1, s = 14 = q₀), and D = 1 with s ≥ 15 would be a sieve witness too weak to reach 1/14. The sharp consequence, verified: hard-stratum families (q₀ > 14, where no small-modulus sieve exists) must have D ≥ 2 — measured **0 of 8 with D = 1**, all at D ≥ 6. LRC(14) therefore runs on TWO mechanisms: sieving (D = 1) for q₀ ≤ 14, and high-determinant pairs for the rest
status: the identification D = 1 ⟺ classical sieve is exact algebra; the tight-family values are exact; the hard-stratum prediction is verified on 8 hard-stratum families (0 with D = 1) alongside 6 easy ones, and every family tested satisfies s ≤ 14D. The two-mechanism decomposition is a structural reading, not a proof of either half
source: opus-2026-07-19-S391 (owner: work the D=1 active pair question)
depends_on: [THM-1205 (the active-pair ratio, which posed this question), boxeph-S120 / HYP-7745 (the located-maximizer theorem), THM-1105 (the hard stratum q₀ > 14), THM-1035/1040 (the classical sieve)]
scripts: 04-computation/dequals1_opus_S391.py -> 05-knowledge/results/dequals1_opus_S391.out
---

# THM-1210 — what D = 1 actually is

## The identification

D = vᵢaⱼ − vⱼaᵢ = 1 forces gcd(vᵢ,vⱼ) = 1. Take the canonical instance
vᵢ = 1, vⱼ = s − 1, so s = vᵢ + vⱼ. Then t* = 1/s, and

- both active runners sit at distance **exactly 1/s**;
- any other runner v satisfies ‖v/s‖ ≥ 1/s **iff s ∤ v**.

That is precisely the **classical sieve at modulus s**, and it delivers gap
exactly 1/s. So:

| D | s | meaning |
|---|---|---|
| 1 | **14** | the sieve at q₀ = 14 — **the tight families** |
| 1 | ≥ 15 | a sieve witness giving 1/s < 1/14 — too weak |
| ≥ 2 | — | a genuinely different mechanism |

Both known tight families confirm it exactly: {1,…,13} and {1,…,11,13,24}
have t* = 1/14, g = 1/14, active pair **(1,13)**, **D = 1**, s = 14 = q₀.

## The sharp consequence, verified

If a family lies in the **hard stratum** (q₀ > 14 — every q ≤ 14 divides some
speed, so no small-modulus sieve exists), then a D = 1 active pair would be a
sieve at modulus s = q₀ ≥ 15, giving gap 1/s < 1/14. Since LRC(14) holds on
that stratum, the maximizer must come from elsewhere. Prediction: **hard-stratum
families have D ≥ 2.**

Measured: **0 of 8 hard-stratum families have D = 1**; all sit at D ≥ 6
(D = 6, 7, 12, 12, 12, 15, 15, 17). Every family tested — easy and hard —
satisfies s ≤ 14D.

## The two-mechanism reading

LRC(14) is carried by two distinct mechanisms, and the D-value tells them
apart:

1. **D = 1, q₀ ≤ 14** — the classical sieve. Covers ~91% of families
   (THM-1105), and contains the extremal cases at s = 14 exactly.
2. **D ≥ 2, q₀ > 14** — the hard ~9%, where no sieve is available and the
   maximizer is a high-determinant pair instead.

This is the cleanest decomposition of the problem the programme has
produced: the sieve is not merely *a* tool that happens to work on most
families, it is *exactly* the D = 1 case of the located-maximizer formula,
and the residual stratum is precisely where D must exceed 1.

## Answering the posed question

**"Is there a 13-family whose active pair has D = 1 and sum ≥ 15?"** Such a
family would be a counterexample to LRC(14), and structurally it would need
the classical sieve at a modulus ≥ 15 to be the *best available* point. The
hard-stratum families — the only ones where small-modulus sieving fails —
do not do this; they find D ≥ 6 maximizers with comfortable margins. So the
D = 1 branch of THM-1205 is, empirically, not where a counterexample lives.

**Not proved.** This locates and explains the branch rather than closing it;
8 hard-stratum families is evidence, not a theorem, and the D ≥ 2 branch
(needing s > 14D) is untouched here.

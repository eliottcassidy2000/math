---
id: THM-774
title: The safe-measure floor over 12-cores — exact bounded-height census (min |G'| = 7/858 at {1..13}∖{6}, the detuning extremal again), the unconditional Lipschitz tail, and the ρ bridge that bounds regime 2's normalized band domain
status: PROVED-EXACT (the maxP ≤ 18 census; the Lipschitz tail; the ρ bridge) + VERIFIED (adversarial scale rays, tooth insertions, hill-descent to height 2500) + CONJECTURE (the asymptotic floor, named precisely)
source: opus-2026-07-14-S301 (owner directive: make progress or prove route no-gos; separate signal from noise)
depends_on:
  - LRC(<=13)  # M(P) >= 1/13 for 12-cores (the Lipschitz tail)
  - THM-755    # v* = r_P/(pi |G'_P|) (the bridge's target)
related: [THM-757, THM-761, THM-767, HYP-6055, HYP-6830, HYP-6835, HYP-6840, MISTAKE-145]
verification: 04-computation/lrc14_gprime_floor_decision_opus_S301.py
  (+ 05-knowledge/results/lrc14_gprime_floor_decision_opus_S301.out)
---

# THM-774 — the safe-measure floor and the ρ bridge

**Why this object.** Regime 2 of the ≥4-far endgame (HYP-6830) needs the band
(maxP, v*] bounded in normalized coordinates: ρ(P) = v*(P)/maxP ≤ K. After
MISTAKE-145 (raw fragmentation r_P is NOT bounded by divisor structure), the only
surviving route through r_P ≤ Σ_P is a POSITIVE FLOOR on the safe measure
|G'_P| = |{t : min_{p∈P} ||pt|| ≥ 1/14}| over 12-cores. This file decides what is
decidable today, exactly, and names the remaining conjecture sharply.

## (1) The ρ bridge (PROVED, one line)

> For every 12-core P:  **ρ(P) = v*(P)/maxP = r_P/(π·|G'_P|·maxP) ≤ 12/(π·|G'_P|).**

*Proof.* The 1/14-safe set is the complement of Σ_{p∈P} p ≤ 12·maxP open arcs, so
its component count satisfies r_P ≤ Σ_P ≤ 12·maxP. Divide. ∎

## (2) The unconditional Lipschitz tail (PROVED, elementary)

> For every 12-core P:  **|G'_P| ≥ 2·(M(P) − 1/14)/maxP ≥ 1/(91·maxP).**

*Proof.* min_p ||p t|| is maxP-Lipschitz; around a maximin point the clearance
stays ≥ 1/14 on an interval of half-width (M(P) − 1/14)/maxP, and M(P) ≥ 1/13 by
settled LRC(13). ∎

(The tail decays with height — the scale-quotient wall — which is exactly why the
floor question below is about SHAPES: |G'| is dilation-invariant, so the infimum
over all 12-cores equals the infimum over primitive shapes.)

## (3) The exact bounded-height census (PROVED-EXACT)

Over ALL primitive 12-element subsets of {1..H}, exhaustively (1,820 / 6,188 /
18,564 shapes for H = 16 / 17 / 18), exact rational interval arithmetic:

> **min |G'_P| = 7/858 ≈ 0.008159, attained uniquely at
> P* = {1,2,3,4,5,7,8,9,10,11,12,13} = {1..13} ∖ {6}, for every H ∈ {16,17,18}.**

| shape | exact `\|G'\|` | note |
|---|---|---|
| {1..13} ∖ {6} | **7/858** | the census minimizer — kps's DETUNING EXTREMAL, one level down |
| {1..11, 13} (GW-core) | 426/35035 ≈ 0.01216 | runner-up |
| {1..11, 26} | 811/35035 ≈ 0.02315 | |
| {1..12} | 6617/194040 ≈ 0.03410 | the deep-well core |

**The double-extremal signal.** {1..14} ∖ {6} is the exact minimizer of M over
covering 13-sets (kps HYP-6055: M = 2/23, "the least detuned divisor-complete
family"); its 12-core twin {1..13} ∖ {6} is now the exact minimizer of the safe
MEASURE over bounded 12-cores. Least-detuned near-AP shapes minimize both the
margin and the room: signal that the same extremal family governs the margin
functional and the measure functional. (Also 858 = 2·3·11·13 — pair-sum ruler
arithmetic; 7/858 = 7·(1/858).)

## (4) The adversarial no-decay verification (VERIFIED, battery discipline of MISTAKE-140/145)

- **Scale rays flat:** {c, 2c, …, 11c, 12c±1} has |G'| → ≈ 0.0483 (constant in c
  through c = 40); GW-dilate one-slot perturbations ≈ 0.048; no decay.
- **Tooth insertions flat:** {1..11, N} ≈ 0.048 through N = 2003; {1..10,12,N},
  {2..12,N} ≈ 0.09; the MISTAKE-145 falsifier mechanism does not depress |G'|.
- **Hill-descent (heights ≤ 2500, replacement + local moves, 5 seeds incl. spread
  and compressed):** every run converges back to the bounded shapes; global best
  found = the census minimizer 7/858 itself.

## (5) The asymptotic floor (CONJECTURE, stated sharply)

> **inf over ALL primitive 12-cores of |G'_P(1/14)| = 7/858, attained only at
> {1..13} ∖ {6}.**

If true, with (1): **ρ(P) ≤ 12·858/(7π) = 10296/(7π) < 469 uniformly** — the
regime-2 band domain is uniformly bounded in normalized (peel-relative)
coordinates, and the ≥4-far endgame becomes [sheets above scale 43 (THM-761)] ∪
[deck pierces (THM-767/771)] ∪ [a normalized band atlas with ρ < 469]. The
measured ρ maximum is 9.335 (HYP-6830) — the bridge is lossy by ~50× but the
LOGIC only needs finiteness. Falsifier shape (for future adversaries): a primitive
shape family with |G'| < 7/858 must beat the least-detuned near-AP at its own
game — cover ≥ 719/726 of the circle with twelve combs of total budget 12/7;
the census says nothing below height 18 does, and every structured ray tested
converges to ≥ 0.048.

## (6) What this closes and what it does not

- CLOSES: the regime-2 boundedness question on every stratum where the floor is
  proved (maxP ≤ 18 exactly; every adversarially tested family), via (1)+(3).
- LEAVES OPEN: the asymptotic floor (5) — a shape-space compactness statement.
  The honest route to it: |G'| < ε forces near-perfect comb covering, which
  forces additive structure (DMNR-flavored rigidity of near-exact arc-comb
  covers); the census minimizer's structure (near-AP, one gap) is where that
  argument should land.
- SEPARATES SIGNAL FROM NOISE: the floor route is SIGNAL (no decay found under
  the full battery discipline; the extremal is a named, already-load-bearing
  shape); the raw-r_P route was NOISE (MISTAKE-145); the Lipschitz tail (2) is
  proved but is NOT the mechanism — measure does not decay along real families,
  heights do not enter.

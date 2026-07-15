---
id: THM-777
title: The safe-measure floor over 12-cores — exact bounded-height census (min |G'| = 7/858 at {1..13}∖{6}, the detuning extremal again), the unconditional Lipschitz tail, and the ρ bridge that bounds regime 2's normalized band domain
status: PROVED (a global explicit positive floor via THM-780; the Lipschitz tail; the rho bridge) + PROVED-EXACT (maxP <= 18 census) + VERIFIED (adversarial battery) + CONJECTURE (the sharp global value 7/858 and uniqueness)
source: opus-2026-07-14-S301 (owner directive: make progress or prove route no-gos; separate signal from noise)
depends_on:
  - LRC(<=13)  # M(P) >= 1/13 for 12-cores (the Lipschitz tail)
  - THM-755    # v* = r_P/(pi |G'_P|) (the bridge's target)
  - THM-780    # explicit global floor 182^(-12)
related: [THM-757, THM-761, THM-767, HYP-6055, HYP-6830, HYP-6835, HYP-6840, MISTAKE-145]
verification: 04-computation/lrc14_gprime_floor_decision_opus_S301.py
  (+ 05-knowledge/results/lrc14_gprime_floor_decision_opus_S301.out)
---

# THM-777 — the safe-measure floor and the ρ bridge

**Why this object.** Regime 2 of the ≥4-far endgame (HYP-6830) needs the band
(maxP, v*] bounded in normalized coordinates: ρ(P) = v*(P)/maxP ≤ K. After
MISTAKE-145 (raw fragmentation r_P is NOT bounded by divisor structure), the only
surviving route through r_P ≤ Σ_P is a POSITIVE FLOOR on the safe measure
|G'_P| = |{t : min_{p∈P} ||pt|| ≥ 1/14}| over 12-cores.  THM-780 now proves
that positive floor explicitly; this file also identifies the much sharper
candidate value.

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

## (5) The global floor (PROVED positive; sharp value CONJECTURED)

THM-780 applies the settled `1/13` margin for every twelve-core at the smaller
threshold `1/14`.  Its phase-cell pigeonhole gives, uniformly and without a
height or primitivity hypothesis,

> **`|G'_P| >= 182^(-12) =
> 1/1320859596446125189798629376`.**

Together with (1), this proves the explicit normalized bound

> **`rho(P) <= 12*182^12/pi < 5.046*10^27`.**

This constant is intentionally crude, but it closes the existence of a finite
normalized regime-2 atlas.  The sharp statement remains:

> **inf over ALL primitive 12-cores of |G'_P(1/14)| = 7/858, attained only at
> {1..13} ∖ {6}.**

If the sharp statement is true, (1) improves the bound to
**rho(P) <= 12*858/(7pi) = 10296/(7pi) < 469**.  The >=4-far endgame then has
the practical decomposition [sheets above scale 43 (THM-761)] union
[deck pierces (THM-767/771)] union [a normalized band atlas with rho < 469]. The
measured ρ maximum is 9.335 (HYP-6830) — the bridge is lossy by ~50× but the
logic only needs finiteness, now supplied unconditionally by THM-780.  Any
falsifier of the sharp value—a primitive shape family with `|G'|<7/858`—must
beat the least-detuned near-AP at its own game: cover ≥ 719/726 of the circle
with twelve combs of total budget 12/7;
the census says nothing below height 18 does, and every structured ray tested
converges to ≥ 0.048.

## (6) What this closes and what it does not

- CLOSES: the global regime-2 boundedness question, with an explicit but huge
  cutoff, by THM-780 plus (1).
- LEAVES OPEN: the sharp asymptotic value and uniqueness in (5), not positivity.
  The honest route to it: |G'| < ε forces near-perfect comb covering, which
  forces additive structure (DMNR-flavored rigidity of near-exact arc-comb
  covers); the census minimizer's structure (near-AP, one gap) is where that
  argument should land.
- SEPARATES SIGNAL FROM NOISE: the floor route is SIGNAL (no decay found under
  the full battery discipline; the extremal is a named, already-load-bearing
  shape); the raw-r_P route was NOISE (MISTAKE-145); the Lipschitz tail (2) is
  proved but is NOT the mechanism — measure does not decay along real families,
  heights do not enter.

---
id: THM-777
title: The sharp safe-measure floor over 12-cores — exact bounded-height census, unconditional ρ bridge, and the 7/858 conjecture
status: PROVED-EXACT (the maxP ≤ 18 census; the Lipschitz tail; the ρ bridge) + VERIFIED (the explicitly listed finite adversarial probes) + CONJECTURE (the sharp 7/858 floor and uniqueness only; THM-780 later proves a much smaller uniform positive floor)
source: opus-2026-07-14-S301 (owner directive: make progress or prove route no-gos; separate signal from noise)
depends_on:
  - LRC(<=13)  # M(P) >= 1/13 for 12-cores (the Lipschitz tail)
  - THM-755    # v* = r_P/(pi |G'_P|) (the bridge's target)
related: [THM-757, THM-761, THM-767, THM-780, THM-783, THM-784, THM-792, THM-793, HYP-6830, HYP-6835, HYP-6840, MISTAKE-141, MISTAKE-145]
verification: 04-computation/lrc14_gprime_floor_decision_opus_S301.py
  (+ 05-knowledge/results/lrc14_gprime_floor_decision_opus_S301.out)
---

# THM-777 — the safe-measure floor and the ρ bridge

**Why this object.** Regime 2 of the ≥4-far endgame (HYP-6830) needs the band
(maxP, v*] bounded in normalized coordinates: ρ(P) = v*(P)/maxP ≤ K. After
MISTAKE-145 (raw fragmentation r_P is NOT bounded by divisor structure), the only
surviving route through r_P ≤ Σ_P is a POSITIVE FLOOR on the safe measure
|G'_P| = |{t : min_{p∈P} ||pt|| ≥ 1/14}| over 12-cores. This file decides what is
decidable today, exactly, and names the remaining conjecture sharply.

**Later closure of the qualitative floor.** THM-780 uses the settled `1/13`
margin and a heavy joint-phase cell to prove the unconditional height-free
bound `|G'_P|>=182^(-12)`. Thus the existence of some positive floor and some
finite normalized `rho` bound is no longer conjectural. What remains open here
is the sharp value `7/858`, its uniqueness, and a practically sized cutoff.

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
| {1..13} ∖ {6} | **7/858** | the unique census minimizer through height 18 |
| {1..11, 13} (GW-core) | 426/35035 ≈ 0.01216 | runner-up |
| {1..11, 26} | 811/35035 ≈ 0.02315 | |
| {1..12} | 6617/194040 ≈ 0.03410 | the deep-well core |

**A bounded detuning signal, not a double extremum.** The deletion of speed 6
recurs in the repository's near-AP detuning studies, but the former claim that
`{1..14} ∖ {6}` globally minimizes `M` over covering 13-sets was withdrawn in
MISTAKE-141. No cross-functional extremality follows here. The proved statement
is only that `{1..13} ∖ {6}` uniquely minimizes safe measure in the complete
primitive census with `maxP<=18`. (Also `858=2*3*11*13`, compatible with the
pair-sum endpoint arithmetic; this factorization is diagnostic, not a proof.)

## (4) Fixed-base exclusion and finite adversarial probes

- **Fixed-base tooth insertion (PROVED, THM-793):** for every fixed 11-speed
  base `B` with safe mass `mu>0` and `r_B` components,
  `|G'_{B union {N}}|>=6mu/7-2r_B/(7N)`. Hence
  `liminf_{N->infinity}|G'_{B union {N}}|>=6mu/7>0`; one new high frequency over
  a fixed base cannot be a safe-mass-decay counterfamily.
- **Scale samples (VERIFIED rows only):** the two perturbations
  `{c,2c,...,11c,12c±1}` were evaluated at
  `c=2,3,5,8,13,21,34,40`, and GW perturbations at `c=2,5,13,29`. Their values
  vary but all listed values lie between `0.04573` and `0.05025`.
- **Tooth samples (VERIFIED rows only):** `{1..11,N}` was evaluated at
  `N=101,503,1009,2003`; the two other named bases were evaluated at
  `N=101,1009`. These samples lie above `0.04800` and `0.09016`, respectively.
- **Hill descent (HEURISTIC):** five 500-step seeded searches with proposal
  height at most 2500 found nothing below `7/858`. One seed was the census
  minimizer itself, and the spread seed stopped at height 727 with mass
  `0.129333`; this is a probe, not convergence or global-optimality evidence.

## (5) The sharp asymptotic floor (CONJECTURE, stated sharply)

> **inf over ALL primitive 12-cores of |G'_P(1/14)| = 7/858, attained only at
> {1..13} ∖ {6}.**

If true, with (1): **ρ(P) ≤ 12·858/(7π) = 10296/(7π) < 469 uniformly** — the
regime-2 peel-ratio coordinate is uniformly bounded, and the ≥4-far endgame
splits into [sheets above scale 43 (THM-761)] ∪ [deck pierces (THM-767/771)] ∪
[the remaining shape/residue atlas inside ρ < 469]. This is not by itself a
finite classification of that atlas. The
measured ρ maximum is 9.335 (HYP-6830) — the bridge is lossy by ~50× but the
LOGIC only needs finiteness. Falsifier shape (for future adversaries): a primitive
shape family with |G'| < 7/858 must cover more than `719/726` of the circle with
twelve combs of total budget `12/7`. The complete census says no primitive
shape through height 18 does; the finite samples above do not decide an
unbounded structured ray.

## (6) What this closes and what it does not

- CLOSES: the practical regime-2 ratio bound on the exact stratum `maxP<=18`, via (1)+(3),
  and rules out a single high-frequency insertion over any fixed base as a
  safe-mass-decay mechanism, via THM-793.
- THM-780 LATER CLOSES: the qualitative global floor and hence the existence of
  a bounded normalized ratio coordinate, with the deliberately crude bound
  `rho(P)<=12*182^12/pi`.
- LEAVES OPEN: the sharp asymptotic identity and uniqueness in (5), plus a
  structurally usable rather than astronomical cutoff and atlas.
  The honest route to it: |G'| < ε forces near-perfect comb covering, which
  forces additive structure (DMNR-flavored rigidity of near-exact arc-comb
  covers); the census minimizer's structure (near-AP, one gap) is where that
  argument should land.
- SEPARATES SIGNAL FROM NOISE: the sharp-floor route remains viable but unproved; the
  raw-`r_P` route is refuted (MISTAKE-145). THM-793 rigorously removes fixed-base
  tooth insertion from the decay mechanisms. The unresolved possibilities are
  varying-base degeneration and iterated or multiscale insertion are no longer
  candidates for making safe mass tend to zero, by THM-780, but they can still
  defeat practical transport or endpoint alignment. A global floor bounds
  `rho`; it does not by itself classify endpoint-owner and
  residue transport.

---
id: HYP-+2897
title: LRC(14) proof<->disproof dialectic -- the q-witness splits 13-sets into NON-covering (witnessed t=1/q) and COVERING (over-determined => VERY lonely, M~0.095>>1/14); the measure inequality is insufficient for barely-safe seeds (equidistribution is the mechanism); M({1..11,13})=1/12 is the adversarial seed
status: dialectic (proof<->disproof), GROUNDED computation; one of my own ideas REFUTED (measure inequality too weak)
source: mac-mini-2026-06-22-S43 (owner: M({1..11,13})=1/12 & zeta(-1)=-1/12; prove AND disprove, go back and forth; in/finite families)
related:
  - THM-523    # kps q-witness: omit a mult of q<=14 => lonely at t=1/q (= my S42 t=1/n witness, generalized)
  - THM-566    # kps adversarial {1..11,13, 84*lcm}: no bounded-denom witness, equidistribution saves it
  - THM-560    # structured + GW sporadic tight tilers (the non-covering boundary)
  - HYP-+2878  # my covering-system / over-determination route
  - HYP-2856   # 3/pi^2 = 1/(2 zeta(2)) Farey floor
---

# HYP-+2897: the proof<->disproof dialectic -- covering sets are forced LONELY

## The premise as a SEED (M({1..11,13})=1/12)
M({1,...,11,13})=1/12 because the set OMITS a multiple of 12, so t=1/12 is a witness (||s/12||>=1/12
for all s, none ==0 mod 12). This is the q=12 case of kps THM-523 (= my S42 t=1/n witness generalized
to all q<=14): **omit a multiple of any q<=14 => lonely at t=1/q.** So {1..11,13} is the SEED of kps's
adversarial THM-566 family: add a large multiple of 12 to make it a covering set.

## The dialectic (back and forth)
**Disproof 1** (multiples of 14 only): too weak -- THM-523 needs a multiple of EVERY q in {2..14}.
**Proof move:** a counterexample must be a COVERING set (contain a multiple of every q<=14). With 13
runners and 13 constraints (q=2..14), this is over-determined.
**Disproof 2** (search covering sets): random covering 13-sets give **min M ~ 0.095 >> 1/14=0.071**
(float lower bound on M, so M confirmed >= 0.095). The over-determination FORCES covering sets to be
spread => very lonely. No counterexample. (argmin {2,3,4,8,10,11,12,14,18,20,24,32,39}.)
**Disproof 3** (adversarial seed + large N, {1..11,13,N}): I first claimed the MEASURE INEQUALITY
[meas(Safe(B)) > |L|*(2/14) => B∪L safe, union bound] defeats it. **REFUTED (discipline):**
meas(Safe({1..11,13}) @ 1/14) = 426/35035 = 0.0122 < 2/14 = 0.143, so the union bound does NOT save it.
The real defense is EQUIDISTRIBUTION (the large N's arcs cannot align with the seed's tiny 0.0122 safe
set without resonance) -- kps THM-566 (no bounded-denom witness, but an equidistributed one exists).

## Synthesis (the emerging proof shape)
LRC(14) <=> every covering set has M >= 1/14. The q-witness (THM-523) handles all NON-covering sets
(M >= 1/q >= 1/14); the tight instances M=1/14 (consec, GW; THM-560) are the NON-covering boundary
(they omit q=14, witnessed by t=1/14). The COVERING sets are forced over-determined => lonely
(M ~ 0.095 >> 1/14 in search). The residual rigor: bounded covering sets (compactness, THM-527) +
unbounded (equidistribution, THM-566) -- the measure inequality is too weak alone (it needs
meas(Safe(B)) > |L|*2/14, which fails for barely-safe seeds; equidistribution does the real work).

## The zeta(-1) = -1/12 connection (thematic, honest)
The Farey floor is 3/pi^2 = 1/(2 zeta(2)) (HYP-2856); zeta(-1) = -1/12 is its FUNCTIONAL-EQUATION
PARTNER (s=2 <-> 1-s=-1: zeta(-1) = 2^{-1} pi^{-2} sin(-pi/2) Gamma(2) zeta(2) = -1/12). The seed's
M=1/12 and zeta(-1)=-1/12 share the number 12; the "1+2+3+...=-1/12" regularization is the large-N ->
infinity limit of the adversarial family (the infinite added speed, equidistributing). INSPIRATION /
thematic, not a rigorous tool -- logged honestly. The rigorous floor uses zeta(2), not zeta(-1).
-> THM-523, THM-566, THM-560, HYP-2878, HYP-2856.

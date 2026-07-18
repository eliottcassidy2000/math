---
id: THM-1100
title: THE EXTENDED RATIONAL-POINT SIEVE — a new route built from the surviving tools: "t = p/q is lonely" is the band condition min(vp mod q, q − vp mod q)·14 ≥ q, which depends ONLY on the speeds mod q, so bounded-denominator loneliness is a FINITE arithmetic condition rather than an analytic one; extending the classical seven moduli from q ≤ 14 to q ≤ 30 raises the kill rate from 86.5% to 100% on 400 random families, INCLUDING all 54 that block every classical modulus. Dilation obstructs (2^m·{1,…,13} needs denominator 14·2^m) so the route runs on PRIMITIVE families, where the minimal denominator was 25 under random sampling, 32 under light adversarial search and 39 under harder search — an increasing sequence that forbids claiming a bound. The natural "largest gap" reformulation is REFUTED as a control: gaps are ~0.001 (predicting denominator ~1000) while actual minima are ~39, so the gaps are not large, they are SITUATED at low-denominator rationals
status: kill rates measured on 400 random families and the 54-family hard stratum; the dilation obstruction verified exactly; the adversarial witness at q = 39 verified exhaustively (no lonely rational with q < 39, gcd = 1). The bounded-denominator statement is a CONJECTURE — three rounds of search gave 25 → 32 → 39, which is the signature of a supremum the search is not reaching
source: opus-2026-07-17-S373 (owner: work the surviving tools toward a new route)
depends_on: [THM-1035/1040 (the classical seven-moduli sieve, the q ≤ 14 case), THM-1050 (dilation invariance, which forces the primitive restriction), THM-1095 (which retired the ledger route this replaces), MISTAKE-154/156 (the sampling lesson applied here)]
scripts: 04-computation/extended_sieve_opus_S373.py, sieve_adversarial_opus_S373.py, largest_gap_opus_S373.py -> 05-knowledge/results/
---

# THM-1100 — a new route, and what it would take

> **CENTRAL CONJECTURE REFUTED (opus-S374), see THM-1105.** The bounded-denominator conjecture below is FALSE. Blocking a modulus q needs only ONE speed with q | v (that runner sits at the origin for every p), so a single speed divisible by lcm(1..Q) blocks every q <= Q at once: V = {lcm(1..Q)} u {12 coprime speeds} is primitive with no lonely rational of denominator <= Q, for arbitrary Q. This is why the searches recorded here climbed 25 -> 32 -> 39 -- the supremum is infinite. WHAT SURVIVES: the position law min-denominator = q0(V), the first modulus dividing no speed (96.8%), with 0/454 failures for q0 <= 14 and 34.8% for q0 > 14.

THM-1095 retired the uniform Bonferroni ledger. This file builds a route
from what survived, centred on the tool with the best kill rate.

## The observation the classical sieve under-uses

t = p/q is lonely iff, for every speed v,

> **min(vp mod q, q − (vp mod q)) · 14 ≥ q**

This depends **only on v mod q**. So "does a lonely rational of denominator
≤ Q exist" is a condition on the speeds mod lcm(1..Q) — **finite and
arithmetic**, not analytic. That is a different kind of problem from the
uniformity that blocked the ledger, which is the point of the pivot.

The classical seven-moduli sieve is exactly the case q ≤ 14, where the
band degenerates to "q divides no speed" (for q ≤ 14 every nonzero residue
passes). For q > 14 the band is a genuine but WIDE constraint — a random
speed fails with probability only ~1/7 — and there are φ(q) numerators to
try per q.

## The kill rate

| Q | random 13-families | the hard stratum (blocks all seven) |
|---|---|---|
| 14 | 86.5% | 0/54 |
| 20 | 99.5% | 52/54 |
| **30** | **100%** | **54/54** |

Every named adversarial family falls, and cheaply: the THM-1055 primitive
failure at 7/15, the THM-1060 L=31 family at 3/23, the AP d=8 at 1/2, the
sum-free {1,3,…,25} at 1/2, and {1,…,13} at 1/14 — the last with equality,
which is exactly right for the extremal family.

## The dilation obstruction, stated correctly

k·{1,…,13} needs denominator 14, 28, 14, 14, 98 for k = 1,2,3,5,7. So
dilation inflates the denominator **only along the 2,7-part of k** —
dilating by any k coprime to 14 leaves it at 14. But 2^m·{1,…,13} needs
14·2^m, which is unbounded, so **no uniform Q exists over all families**
and the route must run on PRIMITIVE families, with THM-1050 reducing the
rest.

## What is NOT established, and why I will not claim it

The minimal denominator over primitive families measured:

| method | max found |
|---|---|
| 600 random families | 25 |
| light adversarial hill-climb | 32 |
| harder hill-climb (45 restarts, 700 steps) | **39** |

Each increase in search effort raised the maximum. That is precisely the
signature of a supremum the search is not reaching, and after MISTAKE-154,
THM-1055 and MISTAKE-156 I will not read a sampled maximum as a bound. The
hill-climb optima do cluster tightly (31–39), which is mildly reassuring,
but clustering is not a proof. The q = 39 witness

> V = {31, 32, 36, 74, 145, 210, 231, 260, 304, 459, 500, 552, 616}

was verified exhaustively: no lonely rational with q < 39, gcd(V) = 1.

## A reformulation I tried and must retract

The natural control is "the largest uncovered gap is bounded below by
L₀ > 0", since an interval of length L contains a rational of denominator
~1/L. **This is refuted as a control.** Measured largest gaps are ~0.001
(the q=39 witness: 0.001018, predicting denominator ~983) while the actual
minimal denominator is 39 — a factor of 25 too pessimistic. Over 150
primitive families the smallest largest-gap was 0.000912, predicting ~1096
against observed minima under 40.

**The gaps are not large — they are SITUATED.** Uncovered gaps cluster at
low-denominator rationals, because that is where the combs align. So the
right object is the *arithmetic position* of the gaps, not their size, and
that is what a proof of the route must exploit.

## The route, stated

> **Bounded-denominator conjecture.** There is an absolute Q₀ such that
> every primitive 13-speed family admits a lonely rational p/q with q ≤ Q₀.

If true with explicit Q₀, LRC(14) for primitive families reduces to a
finite check on residue classes mod lcm(1..Q₀) — infinitely many families
collapsing to finitely many classes, which is exactly the uniformity that
every previous route failed to supply. Evidence: Q₀ ≥ 39; nothing tested
exceeds it; the classical sieve is the Q₀ = 14 fragment covering 86.5%.

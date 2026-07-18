---
id: THM-1105
title: THE ARITHMETIC POSITION LAW — the uncovered gaps sit at q₀(V), the FIRST modulus dividing no speed: min-denominator = q₀ in 96.8% of primitive families, and the regime split is exactly the classical boundary — 0 failures in 454 families with q₀ ≤ 14 (the classical sieve lemma, confirmed), against 34.8% failures in the 46 with q₀ > 14, which is precisely the hard stratum. This also REFUTES THM-1100's bounded-denominator conjecture by an explicit one-line construction: a single speed divisible by lcm(1..Q) blocks every q ≤ Q at once (that runner sits at the origin for every p), so V = {lcm(1..Q)} ∪ {12 coprime speeds} is primitive with no lonely rational of denominator ≤ Q, and Q is arbitrary. The excess min-den − q₀ is likewise unbounded by my searches (5 random → 19 adversarial)
status: the refutation is EXPLICIT and exact (verified: for Q=30 no modulus ≤ 30 fails to divide lcm(1..30); min-denominator 32). The position law is measured over 500 primitive families; the q₀ ≤ 14 regime showing 0/454 failures is the classical lemma and is PROVED elsewhere. The excess bound is NOT established — searches gave 5 then 19, the same escalating signature as before
source: opus-2026-07-17-S374 (owner: work the arithmetic position of the gaps)
depends_on: [THM-1100 (whose central conjecture this refutes), THM-1035/1040 (the classical sieve lemma, which is the q₀ ≤ 14 half of the law), THM-1050 (dilation), MISTAKE-154/156 (the escalating-search lesson, applied again)]
scripts: 04-computation/denominator_unbounded_opus_S374.py, divisibility_law_opus_S374.py -> 05-knowledge/results/
---

# THM-1105 — where the gaps actually sit

S373 established that the uncovered gaps are not large but *situated*, and
named their arithmetic position as the thing to study. This file finds the
position, and in doing so kills S373's own conjecture.

## The refutation, which I should have run before proposing the conjecture

Blocking a modulus q needs only ONE speed with q | v: at t = p/q that
runner sits exactly at the origin, ‖v·p/q‖ = 0 < 1/14, **for every p**. So
a single speed divisible by lcm(1..Q) blocks every q ≤ Q simultaneously.

> V = {lcm(1..Q)} ∪ {12 speeds coprime to it}

is primitive, has 13 speeds, and admits no lonely rational of denominator
≤ Q. Q is arbitrary, so **the minimal denominator is unbounded and
THM-1100's bounded-denominator conjecture is FALSE.** Verified: for Q = 30,
no modulus ≤ 30 fails to divide lcm(1..30), and the family's minimal
denominator is 32.

This is also why the S373 searches kept climbing 25 → 32 → 39: the
supremum is infinite and the hill-climbs were crawling toward it. See
MISTAKE-157.

## The law that survives

Define **q₀(V) = the smallest modulus dividing no speed**. Over 500
primitive families:

| min-den − q₀ | count | share |
|---|---|---|
| **0** | **484** | **96.8%** |
| 1 | 8 | 1.6% |
| 2–5 | 8 | 1.6% |

**The gaps sit at q₀.** And the failures split exactly along the classical
boundary:

| regime | families | law fails |
|---|---|---|
| q₀ ≤ 14 | 454 | **0 (0.0%)** |
| q₀ > 14 | 46 | 16 (34.8%) |

The q₀ ≤ 14 half is the classical sieve lemma — q divides no speed ⟹ 1/q
is lonely — which needs q ≤ 14 because ‖v/q‖ ≥ 1/q must beat 1/14. The
computation confirms it with zero exceptions. **All open content is the
q₀ > 14 regime**, and note that q₀ > 14 says every q ≤ 14 divides some
speed, which is exactly the ~11% hard stratum the classical sieve leaves.
The two independent characterisations of "hard" coincide.

## What is NOT established

The salvageable conjecture is that the *excess* E(V) = min-den − q₀ is
bounded. It is not established: random families gave max 5, adversarial
hill-climbing gave **19** (q₀ = 15, min-den = 34, V = {29, 146, 179, 182,
191, 209, 216, 264, 299, 307, 361, 391, 400}). That is the same escalating
signature that has now misled me three times, so I record it as open and
lean against it.

Curiously the lcm construction has excess **0** for Q = 20, 30, 40 — the
law is exact on the families that refute the absolute bound. The
adversarial families for the excess are a different population: small q₀
(15) with a large minimal denominator.

## Where this leaves the problem

LRC(14) localises cleanly:

1. **q₀ ≤ 14** (≈91% of families): closed by the classical sieve lemma.
2. **q₀ > 14** (≈9%, the hard stratum): needs the extended statement "q
   divides no speed ⟹ some p/q is lonely" for q > 14. True in ~65% of
   these; when false the excess is small but not provably bounded.

That is a sharper localisation than anything the measure-theoretic route
produced in ten sessions, and it is stated in purely arithmetic terms —
divisibility of the speed set — with no analytic content at all.

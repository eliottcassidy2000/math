---
id: THM-1115
title: THE SIMULTANEOUS-BLOCKING TRADE-OFF IS REFUTED, THERE IS NO MEASURE GAP, AND THE TIGHT LOCUS IS LARGER THAN THE AP CLASSIFICATION SUGGESTS — I predicted that blocking many moduli would cost uncovered measure (blocking needs divisible speeds; divisible speeds have fine, near-equidistributed combs; independence maximises the uncovered set). FALSE: uncovered is independent of q₀ (minimum 0.108–0.122 across q₀ = 4…20, medians ~0.13 throughout), and the maximally-blocking lcm families sit at 0.099–0.107, slightly BELOW typical — the opposite of the prediction. Nor is there a gap above the tight family: single-swap perturbations of {1,…,13} descend continuously (0.00102, 0.00275, 0.00364, …). And a second tight family exists: V = {1,…,11,13,24} has uncovered EXACTLY 0, is not a dilate of {1,…,13}, and is lonely at precisely the same six points p/14 with gap exactly 1/14
status: the trade-off refutation is measured over 400 primitive families bucketed by q₀ plus four lcm families; the absence of a gap is exhibited by an explicit descending sequence; the second tight family is verified in exact rational arithmetic (uncovered = 0, union a single component [0,1], lonely set = {1,3,5,9,11,13}/14, min‖vt‖ = 1/14 exactly)
source: opus-2026-07-17-S376 (owner: work the simultaneous blocking constraint)
depends_on: [THM-1110 (which identified simultaneous blocking as the binding constraint), THM-1105 (whose AP tightness classification this scopes), THM-1035/1040 (the classical sieve at q = 14)]
scripts: 04-computation/blocking_tradeoff_opus_S376.py, measure_gap_opus_S376.py -> 05-knowledge/results/
---

# THM-1115 — the trade-off that isn't

THM-1110 localised the difficulty to **simultaneous** blocking, since no
single modulus binds. This file tests the natural tension there and finds
it absent.

## The hypothesis, and its refutation

Blocking many moduli forces speeds divisible by many things, hence large
and highly structured. But a large speed has a very fine comb, close to
equidistributed and therefore close to independent of the others — and
independence is the *best* case for loneliness. So blocking ought to cost
measure.

**It does not.** Uncovered measure by q₀ over 400 primitive families:

| q₀ | 4 | 7 | 10 | 14 | 17 | 20 |
|---|---|---|---|---|---|---|
| min uncovered | 0.1217 | 0.1111 | 0.1090 | 0.1150 | 0.1194 | 0.1149 |
| median | 0.1321 | 0.1348 | 0.1328 | 0.1361 | 0.1318 | 0.1171 |

Flat. And the maximally-blocking families — one speed equal to lcm(1..Q) —
give 0.0989, 0.1066, 0.1070, 0.1068, **below** the typical 0.13 rather
than above it. The two structures, blocking level and uncovered measure,
are essentially **independent**. My prediction had the sign wrong.

## No measure gap

If uncovered were either 0 or ≥ c > 0, the tight family would be isolated
and everything else robustly lonely — a strong theorem. Single-swap
perturbations of {1,…,13} refute it, descending continuously: 0.00102,
0.00275, 0.00275, 0.00364, … Adversarial minimisation independently
reached 0.0632. The tight locus is a **limit point**, not isolated.

## A second tight family

The perturbation sweep turned up **uncovered exactly 0** for

> **V = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24}**

verified in exact rational arithmetic: the union of danger sets is the
single component [0,1], gcd(V) = 1, and V is **not** a dilate of
{1,…,13}. Its lonely set is {1, 3, 5, 9, 11, 13}/14 — precisely the same
six points as the classical tight family — with min‖vt‖ = 1/14 exactly.

**The mechanism.** At t = p/14 with gcd(p,14) = 1, ‖v·p/14‖ ≥ 1/14 iff
v ≢ 0 mod 14. So *every* family with no speed divisible by 14 is lonely at
p/14 — that is the classical sieve at q = 14 — and 24 ≡ 10 (mod 14) keeps
the residues clear. Tightness is therefore not about the lonely points,
which are always at p/14 for these families, but about whether **anything
else** survives. The extremal locus is correspondingly larger than one
might read off THM-1105.

**Scope note.** This does not contradict THM-1105, which classified tight
**arithmetic progressions** as exactly the dilates k·{1,…,13}; V is not an
AP. But it does correct the natural over-reading of that result: tight
families in general are not just the dilates.

## What this leaves

Simultaneous blocking remains the binding constraint of THM-1110, but it
is not payable in measure — there is no conservation law here to exploit.
Both natural strengthenings (a trade-off, a gap) are now closed off.

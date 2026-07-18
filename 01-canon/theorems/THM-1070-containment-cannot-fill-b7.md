---
id: THM-1070
title: THE CONTAINMENT/FRAGMENTATION TECHNIQUE CANNOT FILL THE B7 LEDGER SLOTS — both new slots admit a VALID bound and neither is usable: the S6 containment floor is tight only at k=2 and loses a factor of ~5 per additional speed (exact/floor 3.5, 24.5, 114, 200, 2101 at k=2..6), and the iterated-fragmentation S7 upper bound is ~1190x too loose (target 0.00208, bound ~2.5); this REFUTES the THM-1065 proposal that the same technique is the natural candidate for each remaining slot — it filled the k=2 slot precisely because k=2 is where it is sharp
status: measured exactly (rational arithmetic; fragmentation step 0/400 violations, iterated bound 0/120 violations, containment floor 0 violations at k=4,5,6 — all bounds VALID; the degradation factors are medians over 12 tuples per k, a measured trend)
source: opus-2026-07-17-S367 (owner: work the S6 lower bound and S7 upper bound)
depends_on: [THM-1065 (whose proposed route this refutes), THM-1026 (the five-slot ledger), THM-1025/1012 (the sharp k=2 slot), THM-932 (the fragmentation machinery)]
scripts: 04-computation/b7_slots_opus_S367.py, 04-computation/b7_degradation_opus_S367.py -> 05-knowledge/results/b7_slots_opus_S367.out
---

# THM-1070 — why the ledger does not extend by the same means

> **ROUTE FOUND (opus-S368), see THM-1075.** The k-fold folded identity called for below is the RESONANCE-LATTICE expansion: mu = sum over Lam(a) = ker(Z^k -> Z) of prod hhat(n_i), rank k-1, in cumulant form mu = sum_S (2lam)^(k-|S|) delta(S). This also explains the degradation measured here: k=2 is sharp because rank 1 is CYCLIC (one Bernoulli closes it = THM-965), and rank >= 2 admits no such closure.

THM-1065 located the LRC(14) obstruction as the need for a UNIFORM
level-7 certificate, and proposed that the remaining slots be filled the
way S2 was: "the same sawtooth/containment technique is the natural
candidate for each — S4 already has its separated-regime floor." This
file tests that proposal on the two new B7 slots and refutes it.

## The S7 upper bound (a new direction — every prior slot was a floor)

Pairwise bounds are hopeless here: mu(7-fold) <= rho ~ 1/49, and summing
over C(13,7) = 1716 subsets gives 35 against a target of 0.00208. The
k-fold intersection is genuinely far smaller than any pair overlap, so
the upper bound needs the FRAGMENTATION side. The step

> mu(A n D_x) <= 2*lam*mu(A) + kappa(A)*(2*lam/x)

is **valid** (0/400 violations). Iterated seven times it is also valid
(0/120), and **~1190x too loose** (median; min 433). The boundary terms,
carrying the component count kappa which grows like the sum of the
speeds, swamp the (2*lam)^7 main term.

## The S6 lower bound (iterated containment, the S360 method one level up)

Also **valid** (0 violations at k=4,5,6). Its accuracy:

| k | exact/floor | vs previous |
|---|---|---|
| 2 | 3.50 | — (proved sharp: THM-1012/1025) |
| 3 | 24.50 | x7.0 |
| 4 | 114.33 | x4.7 |
| 5 | 200.08 | x1.8 |
| 6 | 2100.88 | x10.5 |

Roughly a factor of ~5 lost per additional speed; three orders of
magnitude by k=6. It also carries a **regime cost**: positivity needs
a_{i+1} > 7 a_i at every step, so certifying 13 speeds requires a span of
8^12 = 6.9e10 — a regime the seven-moduli sieve already disposes of. The
floor is therefore doubly inadequate: too loose, and valid only where it
is not needed.

## The diagnosis

Each containment step assumes a worst-case alignment of the new comb
against the current set, and each fragmentation step assumes worst-case
boundary placement. Those per-step worst cases are independent
assumptions that **compound multiplicatively while being unable to
co-occur** — a comb cannot be adversarially aligned against all its
predecessors at once. At k=2 there is only one alignment to assume, which
is exactly why the technique is sharp there and nowhere else.

## What this changes

The B7 ledger needs slots accurate to O(1) relative error. Both available
techniques are off by ~10^3. So the route is not "more slots by the same
method"; it is a bound that uses the JOINT alignment structure of the k
combs rather than iterating a per-step worst case. THM-965's folded
identity is the model to imitate — it is exact for k=2 because it tracks
the actual relative position (a+b and b-a mod 14) instead of bounding it
away. The honest target is a **four-variable folded identity** (the S4
slot, THM-1035) and its k-fold generalisation, not another iteration.

**Standing correction.** THM-1065's closing sentence — "every slot filled
so far was filled by containment and counting rather than by new analytic
machinery" — was true and misleading: one slot has been filled, at the
one value of k where the method is sharp. See MISTAKE-155.

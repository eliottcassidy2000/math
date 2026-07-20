---
id: THM-1291
title: THE CF ACTIVE-LEG LAW, PROVED (HYP-7970 promoted) — at any maximizer t* = a/q of any finite family (0 < M < 1/2): (0) the active speeds occupy a SINGLE ± residue class mod q, every straddling pair sums to a multiple kq of q, and D = k·(qM) is quantized; (A) the first positive integer u* whose distance at t* drops to ≤ M is ALWAYS a convergent denominator q_m of a/q, with UNCONDITIONAL proof; (B) under the scope hypothesis H ("u* is itself a speed" ⟺ no non-speed integer below the smallest active speed beats M), u* IS the smallest active speed, it is a convergent denominator, every straddling pair through it has s = kq and D = k·δ_m where δ_m = |q_m·a − p_m·q| — the determinant is the continued-fraction remainder of the active convergent. H is necessary (390 random violations found, (A) intact in all; the S403 controls are exactly H-failures) and holds across the entire known near-floor corpus
status: PROVED (elementary; unimodular-basis lemma + best-approximation argument + the pinch slope analysis; full proofs below). Referee: (R1) 0 failures / ~200k random (a,q,m̃) for (A); (R2) 0/400 for Prop 0; (R3) corpus verdicts exact incl. the two-scale D-set of {1..12,14}; (R4) 390 H-violations, zero (A)-failures. Lean-shaped: Lemma 1 + (A) are Euclid + case analysis; Prop 0 is three lines of modular arithmetic
source: opus-2026-07-19-S404 (owner: prove HYP-7970, full session); law observed opus-S403
depends_on: [HYP-7970 (the observed law), THM-1205/1210 (the (D,s) pinch this refines), classical best-approximation theory (proof self-contained below)]
scripts: 04-computation/lrc_cf_law_proof_referee_opus_S404.py -> 05-knowledge/results/lrc_cf_law_proof_referee_opus_S404.out; discovery data 04-computation/lrc_cf_active_leg_law_opus_S403.py
---

# THM-1291 — the CF active-leg law

**Setup.** V a finite set of nonzero integers (WLOG positive; ‖−vt‖ = ‖vt‖). M = M(V) =
max_t min_{v∈V} ‖vt‖, attained (piecewise-linear, nonzero slopes ±v — hence also NO
plateaus). Assume 0 < M < 1/2, and let t* = a/q be a maximizer in lowest terms. Write
δ(u) := q·‖u·a/q‖ = min_p |ua − pq| ∈ ℤ for u ∈ ℤ. For every active speed v (‖vt*‖ = M):
q ∤ v (else M = 0) and δ(v) = qM =: m̃ ∈ ℤ≥1. Call v *above-integer* if frac(vt*) = M,
*below-integer* if frac(vt*) = 1 − M.

## Prop 0 (active class; straddle quantization) — unconditional

All active speeds are ≡ ±v₀ (mod q) for any one active v₀: two same-side actives v, w
have (v − w)t* ∈ ℤ ⟹ q | v − w; two opposite-side actives have (v + w)t* ∈ ℤ ⟹
q | v + w. Consequently every straddling pair has sum s = kq (k ≥ 1) and represented
determinant D = M·s = k·m̃ ∈ ℤ. ∎

## Pinch Lemma — unconditional

At any maximizer both types are present. Right-derivative of ‖vt‖ at t* is +v for
above-integer actives, −v for below-integer ones (M < 1/2 keeps t* off kinks of inactive
constraints' relevance: the min's one-sided slopes are determined by actives). If no
below-integer active existed, the min would be strictly increasing to the right of t*
(all active slopes +v > 0, inactives slack), contradicting maximality; symmetrically on
the left. ∎ (Hence "pinched" is automatic; S403's operational term.)

## Lemma 1 (unimodular bound) — the classical engine, self-contained

Let p_m/q_m be the convergents of a/q (standard algorithm; when a > q/2 the two
denominators q₀ = q₁ = 1 repeat — take m maximal below). Signed remainders
r_m = q_m·a − p_m·q alternate in sign with |r_m| = δ_m strictly decreasing to
δ_last = 0 (at p/q = a/q itself). **Claim: for 1 ≤ u < q_{m+1}: δ(u) ≥ δ_m, with
equality only if u = q_m** (valid whenever δ_m ≥ 1, i.e. q_{m+1} ≤ q).
*Proof.* (q_m, p_m), (q_{m+1}, p_{m+1}) is a ℤ²-basis (determinant ±1). Write
(u, p) = x(q_m, p_m) + y(q_{m+1}, p_{m+1}) with p a minimizer of |ua − pq|; then
ua − pq = x·r_m + y·r_{m+1}, and r_m, r_{m+1} have opposite signs. Cases: y = 0 ⟹
u = x·q_m with x ≥ 1, δ(u) = x·δ_m ≥ δ_m, equality iff x = 1 iff u = q_m. x = 0 ⟹
u = y·q_{m+1} ≥ q_{m+1}, excluded. x, y same sign ⟹ u ≥ q_m + q_{m+1} > q_{m+1},
excluded. x, y opposite signs ⟹ |x·r_m + y·r_{m+1}| = |x|δ_m + |y|δ_{m+1} > δ_m. ∎

## (A) The first-beating integer is a convergent denominator — unconditional

Let u* := min{u ≥ 1 : δ(u) ≤ m̃} (well-defined: δ(q) = 0). **Then u* = q_m for some m.**
*Proof.* If u* = q it is the last convergent denominator. Otherwise take m maximal with
q_m ≤ u*; then u* < q_{m+1}, so Lemma 1 gives δ(u*) ≥ δ_m = δ(q_m). Thus q_m also
satisfies δ(q_m) ≤ m̃, and minimality forces u* ≤ q_m, i.e. u* = q_m. ∎
(Referee R1: 0/200k failures.)

## (B) The law — under the scope hypothesis H

**H:** u* ∈ V. (Equivalently: every integer u below the smallest active speed has
‖u·t*‖ > M — the two forms are interchangeable by the argument below.)

**Then:** δ(u*) ≥ m̃ (u* is a speed, M is the min) and ≤ m̃ (definition), so u* is
ACTIVE with δ(u*) = m̃ = δ_m exactly; every speed w < u* has δ(w) > m̃ (minimality), so
**u* is the smallest active speed and is a convergent denominator q_m of t***. By the
Pinch Lemma an opposite-side active w exists; by Prop 0 every such pair has
s = u* + w = kq and

> **D = M·s = k·δ_m = k·|q_m·a − p_m·q| — the determinant is the m-th CF remainder of
> the maximizer, times the scale k.**

At k = 1 the partner is w = q − q_m — the q-mirror of the convergent (δ(q − u) = δ(u)
always, opposite side always), which is why the corpus pairs read (q_m, q − q_m):
(1,13) at 1/14, (5,36) at 17/41, (7,120) at 55/127, (5,54) at 23/59, (7,240) at 70/247. ∎

## Corollary 1 (ladder unification)

Families sharing the CF prefix of their maximizers share the active convergent: the
ladder {1,…,11,13,12m} has maximizers with convergent chain 1, 2, 5, 12, … and active
leg 5 throughout; D = δ₂ grows linearly along the ladder (3, 4, 5, 6 at m = 3..6). The
reduced/represented distinction (MISTAKE-173) is structural: at m = 5 the value reduces
(M = 1/13) while the maximizer stays at q = 65 because the speed 13 kills every a/13;
the law lives at the represented level (k = 1 at q = 65).

## Corollary 2 (the 4/55-hunt table)

Every fraction a/s ∈ (1/14, 3/41) has CF [0; 13, 1, k, …] with k ≥ 2 (THM-1268's window
arithmetic), hence convergent denominators 1, 13, 14, 14k+13, … with remainder table
δ(1) = a, δ(13) = s − 13a, δ(14) = 14a − s, δ(14k+13) = …. **Under H, any 4/55-type
realizer's smallest active speed and its numerator are pinned to this finite table**
(e.g. a D = 4, s = 55 realizer with smallest active speed 13 forces a = 13·, δ(13) =
55 − 13a = 4 ⟹ a = 51/13 ∉ ℤ — impossible; leg 14 forces 14a − 55 = 4·(s_pair/55)·…):
the (u*, a) compatibility table is finite and explicit per rung. Feed to the mod-p CRT
stack (boxeph) and the gate frame (death-star): the three constraint systems now compose
on the SAME finite parameter.

## Scope: H is necessary, and exactly delineates S403's boundary

Random families: 390 H-violations found, (A) intact in all (u* still a convergent
denominator — it just isn't a speed, so the smallest ACTIVE speed needn't be one).
The S403 "failures" are precisely H-failures: primes13 (u* = 1 ∉ V; t* = 1/4) and the
mixed control (u* = 2 ∉ V; t* = 7/15). The near-floor corpus satisfies H throughout —
all its families contain {1, …, 11}-type prefixes, and near-floor cores of Wall A type
(covering, near-AP) contain the small integers that make H automatic whenever the
active convergent is small. When H fails, Prop 0 and (A) still hold; only the
identification "active leg = convergent" is lost.

## Handoffs

Lean: Lemma 1 + (A) are kernel-friendly (ℤ²-basis + Euclid; no reals needed — work with
δ on ℤ/q). The (D,s) consumers: THM-1268/1269's "bound D" now reads "bound k·δ_m over
admissible maximizers under H" — Ostrowski territory. mac-mini-S123's rung-lattice
(dominance-crossing moduli) and this law (CF remainders) should be composed — both are
faces of the (D,s) rung lattice. HYP-7970 status: PROVED as this theorem, scope H.

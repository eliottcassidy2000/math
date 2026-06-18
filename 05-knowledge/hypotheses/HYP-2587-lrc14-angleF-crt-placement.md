---
id: HYP-2587
title: ANGLE F — the multi-band CRT placement route for LRC(14)-S3 (HYP-2581d route b). The CRT/covering-system reformulation of the witness, the fixed-modulus refutation, and the q=14·Vmax constructive witness equivalent to ρ*>0.
status: MIXED. PROVED: the residue reduction (witness a/q ⟺ {v·a mod q}⊆SAFE_q ⟺ the 13 speed-forbidden residue classes A_bad(v) do NOT cover Z/q); the equivalence (i)=(iii) (full witness ⟺ discretized lonely density ρ_q>0) at q=14·Vmax; the structural fact every q-covering 13-set contains a multiple of 14. VERIFIED: no fixed modulus q works for all S (each of 91,98,168,182,210,252 fails on a positive fraction of large-Vmax sets — the earlier "q=91 covers all 40" was a small-sample artifact); q=14·Vmax gives a constructive witness for every in-scope set tested (0/400 + 0 near-tight fails); q_min(S)≤~2·Vmax, unrelated to the optimal denominator D. CONCLUSION: CRT does NOT bypass the ρ*≥c0>0 crux — it RE-EXPRESSES it as a uniform covering-system question. LRC(14) NOT proved.
source: mac-mini-2026-06-18-S2 (ANGLE F)
depends_on:
  - THM-527    # the lonely-density reformulation; ρ* = lim ρ_q
  - HYP-2581d  # the unified tight residual, route (b) = multi-band CRT placement
related:
  - OPEN-Q-108 # the uniform fattening lemma = ρ*≥c0>0 = the covering-system non-cover
  - HYP-2584   # bounded-spread compact reduction
external: Lonely Runner Conjecture; covering systems of congruences; Chinese Remainder Theorem.
---

# HYP-2587 — ANGLE F: the CRT/covering-system placement route

## The route (HYP-2581d route b)
Construct the LRC witness τ for an S3 covering 13-set `S = P ∪ L` by Chinese Remainder:
take τ = a/q rational, and choose `a mod q` so every runner lands in its safe band.

## [T1] The residue reduction (PROVED, elementary)
For fixed modulus `q`, a rational τ = a/q is a level-1/14 witness for speed-set `S` iff
for every `v ∈ S`,
  `(v·a mod q) ∈ SAFE_q := { r : 14·min(r, q−r) ≥ q }`
(proof: `‖v(a/q)‖ = ‖(v mod q)(a/q)‖` by 1-periodicity). The danger band
`DANGER_q = {r : min(r,q−r) ≤ ⌊(q−1)/14⌋}` is a single interval of `2⌈q/14⌉−1` residues
around 0. Equivalently a witness a/q exists **iff the 13 speed-forbidden residue classes**
`A_bad(v) = { a : v·a mod q ∈ DANGER_q }` **do NOT cover Z/q** — a covering-system question.
Depends only on `{v mod q}`, so any "fixed q works for all S" claim is a FINITE statement.

## [T2] No fixed modulus works (VERIFIED — refutes the easy hope)
Over a broad search (12 000+ primitive q-covering S3 sets, Vmax up to ~600), EVERY tested
fixed modulus `q ∈ {91, 98, 168, 182, 210, 252}` fails on a positive fraction (5–15%) of
sets once Vmax is large. (The earlier "q=91 covers all 40" was a small-sample artifact —
corrected.) In particular `S* = {1,2,3,5,7,8,9,10,11,12,13,38,42}` has NO witness at q=91.
**So the CRT modulus must SCALE with the set.**

## [T3] The natural CRT modulus q = 14·Vmax (VERIFIED, constructive)
At `q = 14·Vmax` a witness `a/q` exists for **every** in-scope set tested (0 fails over
400 broad sets with Vmax to 3000, plus 0 fails over a near-tight stress sweep). The true
minimal modulus `q_min(S)` is even smaller — typically `≤ 2·Vmax`, often `O(Vmax/10)` — and
is **UNRELATED to the optimal-τ denominator D** (which grows with the binding pair but is a
red herring: constructing a witness ≥1/14 needs far less resolution than the max-min).
Example: `S*` has optimum `M=1/14` at τ=47/56, but a witness also exists at q_min=56 with
τ=9/56, and at q=14·Vmax=588 (even though 56 ∤ 588).

## [T4] Exact equivalence to the measure floor (PROVED reduction)
At `q = 14·Vmax`, the discretized lonely density `ρ_q(S) := #{good a}/q` equals the witness
density, and `ρ_q → ρ*(P,E)` (THM-527). The exact chain (proved up to the THM-527
finite-Vmax `O(1/Vmax)` boundary error, which only makes the ruler heuristic an OVERCOUNT,
`(i) ⊆ (ii)`):
  CRT witness exists  ⟺  ρ_q > 0  ⟺  the 13 classes A_bad(v) miss Z/q.
**So CRT does NOT bypass the ρ*>0 crux** — it re-expresses it. The precise solvability
obstruction to a UNIFORM (all-S) closure is exactly
  `Union_{v∈S} A_bad(v) ≠ Z/{14·Vmax}` uniformly over S  ⟺  ρ* ≥ c0 > 0  (= OPEN-Q-108).

## [T5] What CRT buys (the honest gain)
1. A **constructive, finite, per-set witness with an explicit modulus 14·Vmax** — no search
   over irrational τ, and the optimum's denominator D is irrelevant.
2. A clean **combinatorial restatement of the open crux**: "the 13 speed-forbidden residue
   classes never form a covering system of Z/{14·Vmax}". Each class has ≈ 2·Vmax residues;
   their union must miss Z/q. Heavy overlap (`Σ|A_bad| ≈ 26·Vmax ≫ 14·Vmax`) is exactly why
   free residues (= witnesses) survive — the covering-system view makes the survival margin
   explicit and is a fresh handle on OPEN-Q-108.

## [T6] CRT decoupling and a structural fact
- **STRUCTURAL FACT (PROVED):** every q-covering 13-set contains a speed divisible by 14
  (covering needs a multiple of 14; no speed ≤13 qualifies). So `Vmax` coprime to 14 forces
  the multiple-of-14 to be a NON-max cluster member.
- For `q = 14·m`, `gcd(m,14)=1`, CRT gives `a ↔ (a mod 14, a mod m)`. The level band is a
  single interval mod q, NOT a product of a mod-14 and a mod-m condition, so the decoupling
  is **approximate**: `a mod 14` governs the 2|v / 7|v speeds to leading order, `a mod m`
  governs the cluster, but the boundary couples them at scale 1/14. (Witnesses found cleanly
  for m=169,197,199,211,223,9999.)

## Net
ANGLE F yields a **constructive finite witness** (q=14·Vmax) and a **clean covering-system
reformulation** of the residual, but the uniform-over-S solvability is provably the SAME
inequality as the measure floor ρ*≥c0>0 (OPEN-Q-108). CRT is a constructive packaging and a
combinatorial handle, **not** a bypass. LRC(14) remains open.

Files: `04-computation/lrc14_angleF_crt_placement_mac-mini-S2.py`,
`lrc14_angleF_crt_decouple_mac-mini-S2.py`, `lrc14_angleF_fixedq_witness_mac-mini-S2.py`,
`lrc14_angleF_obstruction_mac-mini-S2.py`, `lrc14_angleF_equivalence_mac-mini-S2.py`,
`lrc14_angleF_match_debug_mac-mini-S2.py`, `lrc14_angleF_final_mac-mini-S2.py`;
results in `05-knowledge/results/lrc14_angleF_*_mac-mini-S2.out`.
→ THM-527, HYP-2581d, OPEN-Q-108, HYP-2584.

---
id: THM-1110
title: THE EXTENDED SIEVE LEMMA FOR q > 14 IS FALSE — explicit counterexample (V = {11,70,77,137,144,156,175,213,226,232,246,262,281}, primitive, 15 divides no speed, yet all eight numerators p ∈ (ℤ/15)* are blocked); the sharp reason is that the forbidden window W_q = {r : min(r,q−r)·14 < q} has |W_q| = 2⌊(q−1)/14⌋+1, so |W_q| = 1 EXACTLY when q ≤ 14 — the classical hypothesis "q divides no speed" is precisely "no speed lands in W_q" and is sufficient only there, which makes the classical threshold sharp and explains THM-1105's regime split structurally. No single-q counting proof can exist at 13 speeds (bad numerators ≤ 13(|W_q|−1) ≈ 1.857q > φ(q), the same 13/7 > 1 obstruction as everywhere). But the counting theorem DOES fire below 13: if s·k_q < φ(q) with k_q = #{w ∈ W_q\{0} : gcd(w,q)=1}, a lonely p/q exists — valid for s ≤ 11 at q = 90, closing the analogue for up to 12 runners and missing LRC(14) by exactly 2 speeds
status: the counterexample is explicit and exhaustively verified (every p coprime to 15 displayed with its blocking speed); the window formula and the |W_q| = 1 ⟺ q ≤ 14 equivalence are elementary and proved; the counting theorem is proved (a union bound) and its reach s ≤ 11 computed over q < 200
source: opus-2026-07-17-S375 (owner: prove the extended sieve lemma for q > 14)
depends_on: [THM-1105 (which posed this as the open target and measured 34.8% failure), THM-1035/1040 (the classical lemma, whose threshold this shows sharp)]
scripts: 04-computation/extended_lemma_false_opus_S375.py, blocking_cost_opus_S375.py -> 05-knowledge/results/
---

# THM-1110 — the lemma is false, and exactly why

THM-1105 named this as the open target: "q divides no speed ⟹ some p/q is
lonely", for q > 14. It had already measured 34.8% failure in that regime.
This file settles it — the statement is **false** — and extracts the sharp
structure.

## The forbidden window, and why 14 is the boundary

p/q is lonely iff, for every speed v, v·p mod q ∉ W_q where

> **W_q = {r : min(r, q−r)·14 < q} = {0, ±1, …, ±⌊(q−1)/14⌋}**,
> **|W_q| = 2⌊(q−1)/14⌋ + 1**

so **|W_q| = 1 ⟺ q ≤ 14**. For q ≤ 14 the only forbidden residue is 0, so
"q divides no speed" *is* the full hypothesis — that is the classical
lemma, and the threshold is sharp rather than an artifact of its proof.
For q ≥ 15 the window widens (W₁₅ = {0,1,14}) and the hypothesis stops
being sufficient.

## The counterexample

> V = {11, 70, 77, 137, 144, 156, 175, 213, 226, 232, 246, 262, 281},
> primitive, and 15 divides no speed.

| p | blocked by | p | blocked by |
|---|---|---|---|
| 1 | v=226 (≡1) | 8 | v=77 (≡1) |
| 2 | v=232 (≡14) | 11 | v=11 (≡1) |
| 4 | v=11 (≡14) | 13 | v=232 (≡1) |
| 7 | v=77 (≡14) | 14 | v=226 (≡14) |

All eight numerators blocked, by just four speeds each killing exactly
|W₁₅|−1 = 2. The covering is tight.

## No counting proof exists at 13 speeds

A speed v coprime to q kills the numerators p = v⁻¹w for w ∈ W_q\{0}, so

> #bad p ≤ 13·(|W_q| − 1) = 26⌊(q−1)/14⌋ ≈ 1.857 q > q > φ(q)

and the union bound never fires — for q = 15, 23, 29, 41, 71 the bad
count (26, 26, 52, 52, 130) exceeds φ(q) (8, 22, 28, 40, 70) every time.
This is the **same 13/7 > 1 obstruction** that killed S₁, the Bonferroni
ledger, and every union bound in this program, now appearing in the
purely arithmetic setting.

## What IS provable — and by exactly how much it misses

Refining: only w coprime to q gives a p coprime to q, so v kills exactly
k_q = #{w ∈ W_q\{0} : gcd(w,q) = 1} numerators. Hence:

> **THEOREM.** If s speeds, none divisible by q, satisfy **s·k_q < φ(q)**,
> then some p/q is lonely.

Maximising φ(q)/k_q over q < 200 gives q = 90 (φ = 24, k = 2, ratio 12),
so the theorem fires for **s ≤ 11**. Runners-up: q = 23 (ratio 11), q = 25,
41, 66, 150 (ratio 10).

So this elementary argument closes the analogue for **up to 12 runners**
and falls **exactly 2 speeds short** of LRC(14). That margin is the
precise, quantitative form of "why 14 is hard": not a vague difficulty,
but a union bound that runs out at 11 speeds.

## The blocking cost

Without a multiple of q, blocking q costs at least ⌈φ(q)/k_q⌉ speeds. Over
q ∈ [15,60] the costs run 3–11, all within the 13-speed budget, so **every
individual q > 14 is blockable**. The binding constraint is therefore
never a single modulus but **simultaneous** blocking across many q — and
the cheap route for that is divisibility (q | v kills every p at once),
which is exactly what the THM-1105 lcm construction exploits.

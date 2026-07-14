---
id: THM-779
title: The r=8 token-walk blocking criterion at the prime-7 lens — full deck blocking is an integer-decidable walk condition (piece surjectivity + wall rainbow + no simultaneous walls); the walk exits within K0 = 5 walls under adversarial tuning, so components with more walls are pierced
status: PROVED (the criterion and its exact decision procedure; referee 4000-point token check + S301 cross-validation) + VERIFIED (K0 = 5 adversarial census; the survivor's exact factorization) + OPEN (the unconditional surjectivity-exit bound, named)
source: opus-2026-07-14-S302, built directly on boxeph/codex-S6's THM-773 token algebra (their suggested pull: absent-eighth-owner transport)
depends_on:
  - THM-773   # the token k_a = -w_a^{-1} round(w_a x), the X^7 - X factorization criterion
  - THM-767   # zero-variance at the prime lens; chamber locking
  - THM-771   # the exact seven-owner defect frame
related: [HYP-6835, HYP-6840, THM-777]
verification: 04-computation/lrc14_r8_token_walk_criterion_opus_S302.py
  (+ 05-knowledge/results/lrc14_r8_token_walk_criterion_opus_S302.out)
---

# THM-779 — the r=8 token-walk blocking criterion

**Frame.** Lens c = 7, core P (|P| = 5), exceptions W = {w_1..w_8}, 7 ∤ w_a. By
THM-773's token algebra (refereed here exactly on 4,000 random rational points):
off its walls, owner a blocks exactly the sheet

`k_a(x) = -w_a^{-1} · round(w_a x)  (mod 7)`,

its walls sit at `x ∈ (1/2 + Z)/w_a` (mesh 1/w_a), crossing a wall steps the token
by `-w_a^{-1}`, and AT its own wall the owner blocks nothing (clearance exactly
1/14). The full-period structure: `x → x+1` carries every token by −1; the walk is
periodic in x with period 7.

## (1) The blocking criterion (PROVED; integer-decidable)

> The deck is fully blocked over an interval J (every sheet strictly blocked at
> every x ∈ J) **iff**
> 1. **(pieces)** on every wall-free piece of J, the eight tokens cover F_7
>    (eight values on seven letters — exactly one collision pair);
> 2. **(walls)** at every wall in J, the seven non-walling tokens are EXACTLY
>    F_7 — equivalently, the walling owner is a member of the collision pair
>    just before its wall;
> 3. **(no simultaneity)** no two owners wall at the same x in J (a double wall
>    leaves six tokens, which cannot cover F_7 — pierced instantly).

*Proof.* Zero-variance at the prime lens (THM-767(2)) gives each owner exactly
one blocked sheet off walls and none at them; coverage of Z_7 by the eight
tokens is (1); at a wall of a the other seven are constant and must cover all
seven sheets, which for seven values means exactly F_7 — and since the full
multiset on the adjacent pieces is F_7 plus one duplicate, removing a leaves F_7
iff a carried the duplicated value (2); (3) is counting. Conversely these
conditions clearly give coverage at every x. ∎

The whole check is an integer walk: sort the walls of the eight meshes, maintain
tokens, test surjectivity — no interval arithmetic. This replaces HYP-6840's
chamber/rainbow search (1,164 pieces of exact Fractions) with O(#walls) integer
steps, and reduces the r=8 escape-hatch question (Q1) to a symbolic-dynamics
question: **can the deterministic token walk remain inside the surjectivity
region SURJ ⊂ (F_7)^8 (density 28·7!/7^8 ≈ 2.45%) across many walls?**

## (2) The chain structure (why runs die)

After a's wall the collision pair is (a, γ) with γ = the unique owner holding
`token_a − w_a^{-1}`; condition (2) at the NEXT wall (owner b) forces
**b ∈ {a, γ}** — the wall-owner schedule (fixed by the eight meshes) must agree
with the deterministic hop-target chain (fixed by the token algebra). Every
owner switch in the schedule is a ~1/7 algebraic coincidence; runs of covered
walls are exponentially expensive.

## (3) The census (VERIFIED, adversarial)

- 120 random 8-tuples (w ≤ 500), quarter-period windows: **median maximal
  blocking run = 1 wall; 90th percentile 2; maximum 4.**
- Annealed to MAXIMIZE the run (400 steps): **K0 = 5 consecutive covered walls**
  at {8, 13, 19, 23, 92, 359, 372, 438} — the adversarial ceiling found.
- HYP-6840's one-moment survivor (P = {5,7,8,13,14}, W = {108,169,143,213,206,
  197,30,162}, x = 19/216) factors exactly as THM-773 predicted: tokens
  [None, 6, 5, 3, 1, 4, 2, 0] — absent owner 108 at its wall, the other seven a
  perfect rainbow — and its blocking run is EXACTLY 1 wall (piece-level failure
  immediately outside). The survivor was the algebra's minimal case, not the
  seed of a blocking family.

## (4) The pierce consequence

> **Any core-safe component containing more than K0 walls cannot be fully
> blocked** — a free sheet (a full 1/14-witness moment) exists inside it. Since
> a component of length λ contains ≈ λ·Σw_a walls, r = 8 full blocking at the
> prime-7 lens fails for every family whose longest core-safe component sees
> more than K0 walls — e.g. Σw ≳ 5/λ. Empirically (S301 + this census) that is
> every family tested; HYP-6840's Q1 is resolved-negative in practice and Q2
> (the positive r=8 pierce) holds with K0 = 5 as the working constant.

## (5) What remains (honest)

- **The unconditional exit lemma (OPEN, named):** prove K0 is finite/small for
  all 8-tuples — e.g. via (2): a blocking run of length K forces K−1 successive
  schedule-vs-algebra coincidences; the schedule's owner-switch pattern is
  three-distance-rigid while the hop targets move on a large orbit unless the
  w's satisfy explicit mod-7 degeneracies; classify the degeneracies (the rigid
  packet class) or bound their runs. This is the r=8 analogue of the
  chamber-locking proof, and DMNR-style tools apply to the hop-target orbit.
- Small components (≤ K0 walls) are outside the census argument; they fall to
  the direct piece check (criterion (1)) — finite per family.
- Lens generalization (7g | c strata, r = p+1 at other primes p ≤ 13 after
  rescaling δ) is mechanical but unwritten.

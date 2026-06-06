---
id: THM-418
title: The prime-shell window lemma — a uniform M ≥ 2/p LRC bound from cyclotomic multipliers
status: PROVED (elementary); the optimal 2/(2n-1) is the ramified refinement (HYP-2240)
source: claudebox-2026-06-03-S625
depends_on:
  - THM-412   # the twisted-shell dodge (single-shell version)
  - THM-415   # quantitative C'(n): optimal min M = 2/(2n-1)
related:
  - HYP-2240  # the cyclotomic-tower refinement to 2/(2n-1) (the ramified frontier)
  - HYP-2230  # the unit-distance / grid-disproof bridge (this is the window-lemma port)
  - THM-407   # doubling shell-transitive iff 2n-1 prime (unramified) — the same dichotomy
---

# THM-418 — the prime-shell window lemma

## Statement

Convention (THM-398): `n` runners, `S = {v_1,…,v_{n-1}}` distinct positive integers, gap `1/n`,
`M(S) = max_t min_i ‖v_i t‖`.

> **THM-418.** Let `p` be a prime with `p ∤ v_i` for all `i`. Among the `p-1` multipliers
> `a ∈ (ℤ/p)^*`, at most `2(n-1)` are "banded"; so the number with `min_i ‖v_i·(a/p)‖ ≥ 2/p` is at
> least `(p-1) - 2(n-1)`. Hence **if `p ≥ 2n` then `M(S) ≥ 2/p`.**
>
> Taking `p₀` = the least prime `≥ 2n` with `p₀ ∤ v_i` (= least prime `≥ 2n` when all `v_i < 2n`),
> `M(S) ≥ 2/p₀`, and by Bertrand `p₀ < 4n`, so **`M(S) > 1/(2n)`** uniformly in `n`.

## Proof

At the witness `t = a/p`, `‖v_i a/p‖ = (1/p)·dist(a v_i mod p, 0)`. A runner is *banded*
(`dist ≤ 1`, i.e. `a v_i mod p ∈ {0, 1, p-1}`) iff `a ∈ {0, v_i^{-1}, -v_i^{-1}}`; since `a` is a
unit and `v_i` a unit (as `p ∤ v_i`), `a v_i ≠ 0`, so exactly the **2** values `a = ±v_i^{-1}` band
runner `i`. Thus at most `2(n-1)` multipliers are banded, and the rest have every runner at
`dist ≥ 2`, giving `min_i ‖v_i a/p‖ ≥ 2/p`. A good multiplier exists once `p-1 > 2(n-1)`, i.e.
`p ≥ 2n`. ∎

## The cyclotomic reading (the window lemma)

The multipliers `(ℤ/p)^* = Gal(ℚ(ζ_p)/ℚ)`; the involution `a ↦ -a` is complex conjugation, so the
banded multipliers come in `±`-pairs (one pair per runner) and the available perspectives are the
`φ(p)/2` pairs of the real subfield `ℚ(ζ_p)^+`. The count `p-1 > 2(n-1)` is exactly the **finite
window lemma**: *there are more modulus-1 rotations (`p-1`) than the width-2 band can block
(`2` per runner)* — the elementary, single-shell analogue of the grid-disproof's
geometry-of-numbers window (cf. [[HYP-2230]]).

## Honest scope: constant-factor, not optimal

`p₀ ≥ 2n` forces `2/p₀ < 2/(2n-1)`, the optimal value (THM-415). The loss is a constant factor
(`ratio (2n-1)/p₀ ∈ [0.64, 0.97]`, `lrc_prime_shell_dodge_s625.out`); `M > 1/(2n)` is
**trivial-strength** — half the conjectured `1/n`. The factor 2 is intrinsic to the counting: at the
*exact* extremal shell `m = 2n-1` the count is tight (`p-1 = 2(n-1)`, guarantee `= 0`) even though a
good multiplier does exist for the extremal configs (THM-415 is attained). Closing the gap needs the
shell `2n-1` itself, where:

- **`2n-1` prime (unramified):** `φ = 2n-2 = 2(n-1)`; the window is critically full. A good
  multiplier exists for the extremal (THM-415) but is not handed over by counting.
- **`2n-1` composite (ramified), e.g. `n=14`, `2n-1 = 27 = 3³`:** the non-unit runners (multiples of
  3) are *automatically* unbanded (a non-unit never equals `±1`), shrinking the obstruction, but the
  shell now carries the 3-adic tower structure — only `2` good multipliers survive for the
  AP-minus-2 extremal (`{13,14} mod 27`), found, not guaranteed.

This unramified/ramified split is exactly THM-407 (doubling is shell-transitive iff `2n-1` is prime
with primitive root 2). The optimal bound `2/(2n-1)` is the **cyclotomic-tower refinement**
([[HYP-2240]]) — the LRC echo of the grid-disproof needing a class field *tower* to pass from
constant-factor to the true exponent.

## Verification
`04-computation/lrc_prime_shell_dodge_s625.py` → `05-knowledge/results/lrc_prime_shell_dodge_s625.out`
(counting bound: 0 violations / 3000 configs; bound `M ≥ 2/p₀` holds exhaustively n=4..7).

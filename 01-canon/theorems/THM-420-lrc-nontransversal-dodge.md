---
id: THM-420
title: The non-transversal dodge — for 2n-1 prime, every non-±-transversal multiple-of-n config is loose (M ≥ 2/(2n-1))
status: PROVED (elementary); reduces C′(n) [2n-1 prime] to the ±-transversal core (HYP-2341)
source: claudebox-2026-06-06-S642
depends_on:
  - THM-398   # LRC(n) <= C′(n): multiple-of-n ⟹ loose
  - THM-412   # the twisted-shell dodge
  - THM-415   # quantitative C′: min M = 2/(2n-1)
related:
  - HYP-2341  # the ±-transversal residual = the quasi-random/Paley core
  - HYP-2321  # Paley = QR / the ±-pair structure
---

# THM-420 — the non-transversal dodge

Convention (THM-398): `n` runners, speed set `S` = `n-1` distinct positive integers, gap `1/n`,
`M(S)=max_t min_i ‖v_i t‖`; `S` loose if `M(S)>1/n`. `C′(n)`: a multiple-of-n config is loose
(⟹ LRC(n)). Here `2n-1 = p` is **prime** (the unramified family: n=15,19,21,22,24,…, all `> 13` open).

## Statement

> Let `p = 2n-1` be prime and `S` a speed set with `p ∤ v_i` for all `i`. If the residues
> `{v_i mod p}` are **not a ±-transversal** of `(ℤ/p)^*` — i.e. some two residues lie in the same
> ±-pair `{u,-u}`, or fewer than `(p-1)/2 = n-1` distinct ±-pairs are hit — then there is a multiplier
> `a ∈ (ℤ/p)^*` with `a v_i ∉ {0, ±1} (mod p)` for all `i`, so the witness `t = a/p` gives
> ```
> M(S) ≥ 2/p = 2/(2n-1) > 1/n   —   S is loose.
> ```

## Proof

At shell `p`, a multiplier `a` is *bad* if some runner is banded, `a v_i ∈ {0,1,p-1}` (distance `≤1`).
Since `p ∤ v_i`, each `v_i` is a unit, so `a v_i = 0` is impossible and `a v_i = ±1 ⟺ a = ±v_i^{-1}`.
Thus the bad set is `B = ⋃_i {v_i^{-1}, -v_i^{-1}}`, of size `≤ 2(n-1) = p-1` (all units), with
**equality iff the `2(n-1)` values are distinct iff the `{v_i^{-1}}` (equivalently the `{v_i}`) hit
`n-1` distinct ±-pairs — a ±-transversal.** If `S` is not a transversal then `|B| < p-1`, so a good
`a ∉ B` exists; for it `min_i dist(a v_i mod p, 0) ≥ 2`, whence `M(S) ≥ 2/p`. Finally `2/(2n-1) > 1/n`
since `2n > 2n-1`. ∎

## Consequence: the residual is the ±-transversals (a major reduction)

A ±-transversal requires the `n-1` residues mod `p` to tile the `(p-1)/2 = n-1` ±-pairs with **no
collision** — vanishingly rare among configs (verified: **0** transversals in 6000 random
multiple-of-15 and multiple-of-19 configs; the `4277/4277` and `4147/4147` non-transversals are all
loose by this dodge, `lrc_transversal_reduction_s642.out`). So:

> **For `2n-1` prime, C′(n) holds for all multiple-of-n configs except possibly the rare
> ±-transversals.** The generic case is now an *elementary one-line dodge*; the entire residual of the
> unramified LRC family is the ±-transversal core.

The remaining transversal configs are verified loose by other means (1-clocks / lower shells / the
measure dodge B′ of THM-398): e.g. the AP-with-top-bumped `{1,…,n-2,n}` is a transversal with
`M = 1/(n-1)` (S630); sampled transversal residuals have `M = 1/4, 1/6, …` ≫ `1/n`
(`lrc_transversal_M_s642.out`) — none a counterexample. See [[HYP-2341]] for the residual structure.

## Why the transversals are the hard core (the reframe)

A ±-transversal mod `p` hits **every** ±-class once — it is the *maximally spread / spectrally flat*
residue pattern, the **Paley/quasi-random** structure (the QR set is one such transversal; the random
tournament's residues are transversal-like, S641). So the elementary dodge clears every config *except*
the quasi-random ones, and the LRC frontier (for unramified `n`) **is exactly the quasi-random core** —
the same place the character-ratio/Gauss-sum (S638) and the Rado-tournament limit (S641) live. The
divisible case (`p | v_i`) is excluded above: such `v_i ≥ 2n-1` and, if also `≡0 mod n`, is
`≥ n(2n-1)` (dominant ⟹ B′).

---
id: THM-402
title: Round (locally-transitive) tournaments are 2-dichromatic; χ is constant (=2) on the LRC-tight set
status: PROVED
source: opus-2026-06-03-S592
depends_on:
  - HYP-2133   # LRC accessible tournaments are round (S591)
related:
  - THM-401    # the additive (interval) modulus
---

# THM-402 — round tournaments are 2-dichromatic; χ ≡ 2 on the LRC-tight set

## Statement

Let `T` be a **round** (equivalently *locally transitive*: every out- and in-
neighbourhood induces a transitive subtournament) tournament on `m` vertices. Then its
**dichromatic number** (the least number of transitive parts covering the vertices)
satisfies
```
χ(T) ≤ 2,   with χ(T) = 1  iff  T is transitive,  and  χ(T) = 2  iff  T has a 3-cycle.
```
Consequently `χ` is **constant** on every set of round non-transitive tournaments — in
particular on the **LRC-tight set** (all of whose members are round and contain a
3-cycle): `χ ≡ 2`.

## Proof

A round tournament is realizable on the circle: there exist distinct angles
`θ_1,…,θ_m ∈ ℝ/ℤ` with `i → j` iff the clockwise arc from `θ_i` to `θ_j` is `< 1/2`
(Moon; this is the standard characterization of locally transitive tournaments).

Choose a diameter of the circle through no point `θ_i` (possible since the points are
finite). It splits `{θ_i}` into two open semicircles `A` and `B`. Take any two vertices
`u, w ∈ A` with `θ_u` before `θ_w` in clockwise order *within* `A`. Since both lie in a
common semicircle, the clockwise arc `θ_u → θ_w` is `< 1/2`, so `u → w`. Hence the
clockwise order on `A` is a linear order realizing `T[A]` as **transitive**; likewise
`T[B]`. So `{A, B}` is a partition into two transitive parts and `χ(T) ≤ 2`.

`χ(T)=1` iff `T` is transitive (no 3-cycle). If `T` has a 3-cycle it is not transitive,
so `χ(T) ≥ 2`, whence `χ(T)=2`. ∎

## LRC corollary

By HYP-2133/S591 the LRC comparator (the half-turn / circular-position rule) produces
**only round tournaments**, and every LRC-tight configuration's runner tournament is
round and non-transitive (it contains 3-cycles — e.g. the regular rotational `R_m`). By
THM-402, `χ = 2` for **all** of them — the regular `R_m` *and* the non-vertex-transitive
sporadics alike. Therefore:

> **`χ` is constant (`= 2`) on the entire LRC-tight set** — closing the qualifier left
> open in S591/HYP-2133. (Paley `m≥7` is **not** round, has `χ=3`, and is LRC-inaccessible
> — the foil confirming LRC is the additive/interval, not the QR/multiplicative,
> phenomenon, THM-401.)

## Verification

`lrc_round_2dichromatic_s592.py`: over all round tournaments on `m=5` (544 of them: 120
transitive `χ=1`, 424 with `χ=2`, **0 with `χ≥3`**) and a sample on `m=7` (147 round, all
`χ≤2`) — exactly matching the theorem.

**Artifacts:** `04-computation/lrc_round_2dichromatic_s592.py` (+`.out`). Builds on
HYP-2133 (LRC ⟹ round), THM-401 (additive interval). New: **HYP-2134**.

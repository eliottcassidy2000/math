# THM-635 — Translate-block bound: {m,…,m+11} has reach ≥ 2/25 for m ≥ 2

**Status:** VERIFIED + FORMALIZED sorry-free & kernel-pure (`LRCCoveringReach.lean`)
**Author:** mac-mini-2026-07-06-S34
**Relevance:** closes opus-S127's branch-3 uniform-k escape of the finite-covering crux (C)

## Statement

The consecutive 12-block `B(m) = {m, m+1, …, m+11}` (m a positive integer) has, for
every `m ≥ 2`,

> **M(B(m)) ≥ m/(2m+11) ≥ 2/25.**

Only `m = 1` (the AP `{1,…,12}`) is tight (`M = 1/13`). Formally
(`reach = sSup (margin (Blk m) '' [0,1])`):
`translate_block_reach : 2 ≤ m → 2/25 ≤ reach(Blk m)`.

## Proof

Witness `t = 1/(2m+11)`. Every speed `v ∈ {m,…,m+11}` satisfies `m ≤ v < 2m+11`, so
`v mod (2m+11) = v ∈ [m, m+11] = [m, (2m+11) − m]` — off the forbidden band `{0,…,
±(m−1)}`. Hence `‖v·t‖ ≥ m/(2m+11)`, so `M ≥ m/(2m+11)`. Finally `m/(2m+11) ≥ 2/25 ⟺
25m ≥ 4m+22 ⟺ m ≥ 2`. ∎

Exact values (`M = m/(2m+11)`, verified): `1/13, 2/15, 3/17, 4/19, 5/21, 6/23, 7/25,
…` — the AP at the bottom, then a jump to `2/15 ≫ 2/25`.

## Role in (C)

This is a **direct instance of the covering-reach atom** `reach_ge_of_covering`
(THM/S34) at `q = 2m+11, c = 1, μ = m`. In opus-S127's four-branch skeleton of (C):
branch 3 is the L-lift escape `V ≡ AP mod L`; its **uniform-k** case is exactly a
translate of the AP, `{1+Lk,…,12+Lk} = B(1+Lk)` with `1+Lk ≥ 2`, so THM-635 makes it
loose. (The mixed-k case is closed separately by decorrelation / the scale gap.)

## Formalization

`04-computation/lean/.../LRCCoveringReach.lean`, axioms `[propext, Classical.choice,
Quot.sound]`: `Blk`, `translate_block_reach`.

→ opus-S127 (L-lift escape is loose), kps-S43/S46 (finite covering, no height bound);
THM-633 (d=1 ladder, same reach style); THM/S34 `reach_ge_of_covering` (the atom).

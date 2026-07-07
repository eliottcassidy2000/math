# THM-636 — The decorrelation atom and the escape families' looseness (r ≤ 11)

**Status:** VERIFIED + FORMALIZED sorry-free & kernel-pure (`LRCDecorrelation.lean`)
**Author:** mac-mini-2026-07-06-S38 (owner's height-descent inspiration)
**Relevance:** the first rigorous handling of the S36 covering-escape obstruction (the compressed `≡ AP mod L` families that broke the finite covering)

## The decorrelation atom

For `V, K, b : Fin 12 → ℤ` with `vᵢ = bᵢ + L·kᵢ` (a bounded base + an `L`-scaled
lift), `|bᵢ| ≤ B`, `0 < L`:

> **`reach(K) − B/L ≤ reach(V)`**   (`reach X = sSup (margin X '' [0,1])`).

Proof: `distZ` is 1-Lipschitz, so at the witness `t = t_K/L` (`t_K` a max-margin
point of `K`), for each `i`
`‖vᵢ·t_K/L‖ = ‖kᵢ t_K + bᵢ t_K/L‖ ≥ ‖kᵢ t_K‖ − |bᵢ t_K/L| ≥ margin(K,t_K) − B/L`.
So `margin(V, t_K/L) ≥ reach(K) − B/L`, hence `reach(V) ≥ reach(K) − B/L`. ∎
(`margin_decorr` per-time; `reach_decorr` the reach form.)

## The escape families are loose (r ≤ 11)

The S36 covering-escape families are `vᵢ = bᵢ + L·kᵢ` with `bᵢ ∈ {1,…,12}` (`B =
12`) the AP residues and `kᵢ ≥ 1` the lifts. Let `r = #{distinct kᵢ}`.

> **`r ≤ 11` (a repeated lift = a close pair):** the lift family `K` has ≤ 11
> distinct nonzero speeds, so LRC(≤12) gives `reach(K) ≥ 1/12`, hence
> **`reach(V) ≥ 1/12 − 12/L > 2/25` for `L > 3600`** — LOOSE, not in the gap.
> (`1/12 − 12/3600 = 2/25` exactly.)

`escape_loose_of_lift_floor` (kernel-pure): given `reach(K) ≥ 1/12` (the LRC(≤12)
citation), `L > 3600`, `vᵢ = bᵢ + L kᵢ`, `|bᵢ| ≤ 12`, concludes `reach(V) > 2/25`.
This handles the originally-found escape families (`k ∈ {1,2}`, `r = 2`).

## The `r = 12` case (descent)

All `kᵢ` distinct ⟹ `K` is a smaller 12-family (height ÷ L). `reach_decorr` gives
`reach(V) ≥ reach(K) − 12/L`; by induction on height (base = the finite covering,
which works at bounded height — no escape there), `reach(K) ≥ 2/25`, so
`reach(V) > 2/25` **when `reach(K) > 2/25` strictly**. The boundary case
`reach(K) = 2/25` (e.g. `K` = the block-lift) gives only `2/25 − 12/L` from
decorrelation, but the true `reach(V) = 1/12 > 2/25` (verified) — loose, needing a
sharper bound (the base structure). So `r = 12` closes modulo the `reach(K) = 2/25`
refinement.

## Significance

This is the correct technique (Tao-style height descent) for the escape families —
the exact class that broke the finite covering (S36/S37). The `r ≤ 11` case is fully
formalized. Note: per opus-S130 (MISTAKE-117), `(C)` does not prove LRC(14) (the J-K
top link is invalid); but the decorrelation-descent is a real advance on `(C)` (a
true, hard statement) and is the right tool for a **direct** LRC(14) descent
(Route 1 / Tao) as well.

## Formalization

`LRCDecorrelation.lean`, axioms `[propext, Classical.choice, Quot.sound]`:
`distZ_lipschitz`, `margin_decorr`, `reach_decorr`, `escape_loose_of_lift_floor`.

→ owner's inspiration (M(V) → M(K), r < 12 = close pair loose, r = 12 descend);
S36/S37 (the escape families); THM-635 (translate = the uniform-lift case),
THM-633 (d=1 ladder), the reach atom S34; LRC(≤13) citation.

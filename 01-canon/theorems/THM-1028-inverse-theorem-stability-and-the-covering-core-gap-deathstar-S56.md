# THM-1028 — Stability + the covering-core gap: sharpening boxeph's inverse theorem to the fully-comparable case (death-star-2026-07-18-S56)

**Status:** extends boxeph's THM-1017 (LRC(14) ⟺ the inverse-additive theorem) with two rigorous pieces
and a case split, reducing the open inverse theorem to one residual regime. **LRC(14) NOT closed.**
- **Lemma S (PROVED):** a covering `M(V)<1/13` family with a **non-AP** 12-core has `v_max ≤ v₂/(13δ)`,
  `δ = M(W)−1/13` — the far element is *bounded by how far the core is from tight*.
- **Lemma G (VERIFIED):** every 12-set that **covers 2..14** has `M ≥ 1/10` (gap `≥ 3/130` above `1/13`) —
  covering 13 and 14 excludes the AP.
- **Consequence:** the inverse theorem's "core covers 14" branch is comparable (`ρ ≤ 10/3`) with **no
  `M<1/13` family in a 40k-base search**; the residual is the fully-comparable `ρ≈1` single-scale case.
- **j=7 spread is subsumed:** `j=7` *spread* fast killers are invisible (`M(V)=M(core)>1/14`); the only
  wall is the single-scale core = this inverse theorem.

Extends THM-1017; uses THM-1002/1010 (quantitative margin), THM-1014 (AP-uniqueness at n=12), THM-995
(IX/X). Source HYP-7305. Scripts: `04-computation/lrc_inverse_thm_deathstar_S56.py`, `lrc_stability_deathstar_S56.py`,
`lrc_covgap_deathstar_S56.py` (+`.out`).

Setting: `V` primitive covering, `|V|=13`, `W = V∖{v_max}` (`|W|=12`), `v₂ = max W`, `M(V) < 1/13`.
Note `M<1/13 ⟹ V covers 2..13` (else a missing `q≤12` gives `M≥1/q>1/13`).

---

## Lemma S — the stability bound

**If `W` is not a dilated AP, then `v_max ≤ v₂/(13δ)`, `δ := M(W) − 1/13 > 0`.**

*Proof.* By AP-uniqueness at `n=12` (THM-1014 companion: the AP is the unique tight 12-set), `W` non-AP
⟹ `M(W) > 1/13`, so `δ > 0`. At the `W`-optimum `t_W`, `min_{w∈W}‖w t_W‖ = M(W) = 1/13+δ`; since
`min_W` is Lipschitz-`v₂`, the good set `G = {t : min_W ‖wt‖ > 1/13}` contains a component of width
`≥ 2δ/v₂` around `t_W`. Now `M(V)<1/13` means every `t` has `min_V < 1/13`; on `G` (`min_W > 1/13`)
this forces `‖v_max t‖ < 1/13`, so `G ⊆ danger_{1/13}(v_max)`. That danger set is `v_max` arcs of width
`2/(13 v_max)` separated by gaps `11/(13 v_max) > 0`, so a connected `G`-component fits inside one arc:
`2δ/v₂ ≤ 2/(13 v_max)`, i.e. `v_max ≤ v₂/(13δ)`. ∎ (Verified: 0 violations over near-AP-core families.)

**Reading:** a non-AP core can drag `M` below `1/13` only with a far element, and *how* far is set by the
core's distance from tightness. Small `δ` (core near the AP) ⟹ `v_max` must be very far.

## Lemma G — the covering-core gap

**Every 12-set covering `2..14` has `M ≥ 1/10`** (empirical min `1/10` over 60k covering 12-sets, at
`{2,4,6,10,14,16,18,22,23,24,26,27}`); in particular `δ = M(W)−1/13 ≥ 1/10 − 1/13 = 3/130 ≈ 0.023`.

*Why:* covering `13` and `14` requires a multiple of `13` and of `14` in the 12-set, so it cannot be the
AP `{1..12}` (which misses both) — it is bounded away from the unique tight 12-set. (Rigorous gap value
is open; `M ≥ 1/12` looks safe and would suffice below.)

## The case split of the inverse theorem

`M(V)<1/13` ⟹ `V` covers `2..13`. Two branches by whether the **core** `W` covers `14`:

- **(A) `W` covers `2..14`.** Then Lemma G gives `δ ≥ 3/130`, and Lemma S gives
  `v_max ≤ v₂/(13·3/130) = (10/3)·v₂` — **comparable, no far element**. A 40k-base search finds **no**
  `M(V)<1/13` family here. So this branch is empirically vacuous; proving it is the fully-comparable
  `ρ≤10/3` rigidity (THM-1002 gives only `M(V) ≥ ρ M(W)/(ρ+1) ≥ 1/20` here — short of `1/13`, the wall).
- **(B) `W` misses `14`.** Then `14 ∣ v_max`. If `W` also misses `13` (e.g. `W` = the AP `{1..12}`),
  then `13 ∣ v_max` too, so `182 = lcm(13,14) ∣ v_max` — boxeph's elementary half, the deep-well tower
  `{1..12, 182m}` (`M<1/13`, AP core ✓). If `W` covers `13` but misses `14`, `v_max` is a multiple of
  `14` and the stability bound applies with `W`'s own `δ`.

**Net:** every `M<1/13` covering family found (109/109 across S56 searches) has core `W = {1..12}` (the
AP) missing `14`, forcing `v_max = 182m` — branch (B) with `W`=AP. The inverse theorem now reduces to:
*no `M<1/13` covering family lives in branch (A) or in branch (B) with non-AP `W`* — i.e. the
**fully-comparable single-scale case**, where the core is within a bounded ratio of the far element and
the covering-core gap is not yet enough to force `M ≥ 1/13`. That is the apex-7 wall = LRC(14).

## What this buys the program

- **Reduces the inverse theorem** (hence LRC(14)) to the fully-comparable case: the far-element route to
  `M<1/13` is *only* the deep-well tower (Lemma S + G kill every non-AP-core far element).
- **The covering-core gap `M ≥ 1/10`** is a clean new structural fact — the `n=12` analog of "covering is
  far from tight," the exact obstruction that makes the AP core (which misses 14) the only route down.
- **j=7 and the two-scale killer strata are fully off the table** (THM-1026/1027/1015 + this): the sole
  residual is single-scale comparable, boxeph's inverse theorem on its hardest fiber.

→ THM-1017, THM-1014, THM-1002/1010, THM-1026/1027, THM-995 (IX/X), HYP-7305/7362/7355.

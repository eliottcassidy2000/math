# THM-1028 — Stability + the covering-core gap: sharpening boxeph's inverse theorem to the fully-comparable case (death-star-2026-07-18-S56)

**Status:** CONDITIONAL + EMPIRICAL NOTE; it does **not** rigorously reduce the open inverse theorem.
**LRC(14) NOT closed.**
- **Lemma S (CONDITIONAL):** assuming the still-open global n=12 AP-uniqueness statement
  (HYP-7310), a covering `M(V)<1/13` family with a **non-AP** 12-core has
  `v_max ≤ v₂/(13δ)`, `δ = M(W)−1/13`. THM-1014 proves only the dilated-AP
  one-killer family; it does not discharge this assumption for arbitrary 12-cores.
- **Lemma G (EMPIRICAL):** a 60k-family census found every 12-set that **covers 2..14**
  had `M ≥ 1/10`. No uniform proof of this numerical gap is supplied here.
- **Conditional consequence:** Lemma S + Lemma G would make the "core covers 14" branch
  comparable (`ρ ≤ 10/3`). The reported absence of `M<1/13` families in a 40k-base
  search is evidence, not a proof that the residual is the only regime.
- **j=7 interpretation (conditional):** the same evidence suggests spread fast killers are
  invisible and the comparable core is the wall; this file does not prove that global dispatch.

Explores THM-1017; uses THM-1002/1010 (quantitative margin), assumes HYP-7310
(global AP-uniqueness at n=12; THM-1014 covers only one structured family), THM-995
(IX/X). Source HYP-7305. Scripts: `04-computation/lrc_inverse_thm_deathstar_S56.py`, `lrc_stability_deathstar_S56.py`,
`lrc_covgap_deathstar_S56.py` (+`.out`).

Setting: `V` primitive covering, `|V|=13`, `W = V∖{v_max}` (`|W|=12`), `v₂ = max W`, `M(V) < 1/13`.
Note `M<1/13 ⟹ V covers 2..13` (else a missing `q≤12` gives `M≥1/q>1/13`).

---

## Lemma S — the conditional stability bound

**Assuming global n=12 AP uniqueness: if `W` is not a dilated AP, then
`v_max ≤ v₂/(13δ)`, `δ := M(W) − 1/13 > 0`.**

*Conditional proof.* By the open AP-uniqueness conjecture at `n=12` (HYP-7310), `W` non-AP
⟹ `M(W) > 1/13`, so `δ > 0`. At the `W`-optimum `t_W`, `min_{w∈W}‖w t_W‖ = M(W) = 1/13+δ`; since
`min_W` is Lipschitz-`v₂`, the good set `G = {t : min_W ‖wt‖ > 1/13}` contains a component of width
`≥ 2δ/v₂` around `t_W`. Now `M(V)<1/13` means every `t` has `min_V < 1/13`; on `G` (`min_W > 1/13`)
this forces `‖v_max t‖ < 1/13`, so `G ⊆ danger_{1/13}(v_max)`. That danger set is `v_max` arcs of width
`2/(13 v_max)` separated by gaps `11/(13 v_max) > 0`, so a connected `G`-component fits inside one arc:
`2δ/v₂ ≤ 2/(13 v_max)`, i.e. `v_max ≤ v₂/(13δ)`. ∎ (Verified: 0 violations over near-AP-core families.)

**Conditional reading:** once global AP uniqueness supplies `δ>0`, how far a non-AP core can drag
`M` below `1/13` is controlled by its distance from tightness. Small `δ` permits a farther element.

## Lemma G — the empirical covering-core gap

**Empirical claim:** every tested 12-set covering `2..14` has `M ≥ 1/10` (minimum `1/10` over 60k
covering 12-sets, at
`{2,4,6,10,14,16,18,22,23,24,26,27}`); in particular `δ = M(W)−1/13 ≥ 1/10 − 1/13 = 3/130 ≈ 0.023`.

*Why:* covering `13` and `14` requires a multiple of `13` and of `14` in the 12-set, so it cannot be the
AP `{1..12}` (which misses both) — it is bounded away from the unique tight 12-set. (Rigorous gap value
is open; `M ≥ 1/12` looks safe and would suffice below.)

## The conditional case split of the inverse theorem

`M(V)<1/13` ⟹ `V` covers `2..13`. Two branches by whether the **core** `W` covers `14`:

- **(A) `W` covers `2..14`.** Then Lemma G gives `δ ≥ 3/130`, and Lemma S gives
  `v_max ≤ v₂/(13·3/130) = (10/3)·v₂` — **comparable, no far element**. A 40k-base search finds **no**
  `M(V)<1/13` family here. So this branch is empirically vacuous; proving it is the fully-comparable
  `ρ≤10/3` rigidity (THM-1002 gives only `M(V) ≥ ρ M(W)/(ρ+1) ≥ 1/20` here — short of `1/13`, the wall).
- **(B) `W` misses `14`.** Then `14 ∣ v_max`. If `W` also misses `13` (e.g. `W` = the AP `{1..12}`),
  then `13 ∣ v_max` too, so `182 = lcm(13,14) ∣ v_max` — boxeph's elementary half, the deep-well tower
  `{1..12, 182m}` (`M<1/13`, AP core ✓). If `W` covers `13` but misses `14`, `v_max` is a multiple of
  `14` and the stability bound applies with `W`'s own `δ`.

**Empirical net:** every `M<1/13` covering family found (109/109 across S56 searches) has core
`W = {1..12}` (the AP) missing `14`, forcing `v_max = 182m` — branch (B) with `W`=AP. If the
global AP-uniqueness assumption in Lemma S and the uniform numerical gap in Lemma G were proved,
this route would reduce the remaining search to the **fully-comparable single-scale case**. As
written, the census identifies that regime as a bottleneck; it does not prove the reduction.

## What this buys the program

- **Conditional route:** global AP uniqueness plus a proved covering-core gap would reduce the
  inverse-theorem search to a comparable regime; this file establishes neither global input.
- **Candidate invariant:** the observed covering-core gap `M ≥ 1/10` is a useful structural
  conjecture, not a theorem.
- **Scope warning:** these searches do not put `j=7` or all two-scale killer strata rigorously off
  the table. They identify the comparable single-scale regime as an empirical bottleneck.

→ THM-1017, THM-1014, THM-1002/1010, THM-1026/1027, THM-995 (IX/X), HYP-7305/7362/7355.

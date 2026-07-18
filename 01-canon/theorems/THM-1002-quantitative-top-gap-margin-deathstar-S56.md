# THM-1002 — The quantitative top-gap margin (death-star-2026-07-17-S56)

**Status: PROVED.** A sharp, quantitative refinement of THM-1000. For any family of positive speeds,
```
M(V) ≥ Vmax · M(V∖{Vmax}) / (Vmax + v₂),
```
`v₂` = second-largest speed. With `M(V∖{Vmax}) ≥ 1/(n−1)` (LRC(n−1), settled for n≤14):
`M(V) ≥ Vmax/((n−1)(Vmax + v₂))`. Verified on AP, GW, deep well, far-element, covering-min, 2·AP,
comparable families — holds everywhere, and is **tight in the far-element regime** (`{1..12,500}`:
bound `0.0751` vs actual `M=1/13=0.0769`). Source HYP-7305; extends THM-1000. Script
`04-computation/lrc_quant_margin_deathstar_S56.py` (+`.out`). WLOG positive speeds (`‖vt‖=‖|v|t‖`).

---

## Statement and proof

**Theorem.** Let `V` be a finite family of positive integers, `w = Vmax = max V`, `V' = V∖{w}`,
`v₂ = max V'`. Then `M(V) ≥ w·M(V')/(w + v₂)`.

**Proof.** Fix any `λ < w·M(V')/(w+v₂)`; we show `M(V) ≥ λ`.
1. **A wide good-component of `V'`.** Let `t₀` attain `M(V')`: `min_{v∈V'} ‖v t₀‖ = M(V')`. The map
   `f(t) = min_{v∈V'} ‖vt‖` is Lipschitz with constant `v₂` (each `‖vt‖` has slope `≤ v ≤ v₂`), so
   `f(t) > λ` for `|t − t₀| < (M(V') − λ)/v₂`. Hence the component `I₀` of `{f > λ}` containing `t₀`
   has length `|I₀| ≥ 2(M(V') − λ)/v₂`.
2. **`w`'s arcs can't cover `I₀`.** The set `{‖wt‖ ≤ λ}` is a union of arcs of width `2λ/w` centered
   at `k/w`, separated by gaps of width `1/w − 2λ/w = (1−2λ)/w > 0`. From
   `λ < w M(V')/(w+v₂)` one gets `2(M(V')−λ)/v₂ > 2λ/w`, so `|I₀| > 2λ/w = ` one arc width. A
   connected interval wider than one arc must contain a point of some gap — i.e. a `t* ∈ I₀` with
   `‖w t*‖ > λ`.
3. **`V` is λ-lonely at `t*`.** There `min_{v∈V'} ‖v t*‖ > λ` (as `t* ∈ I₀`) and `‖w t*‖ > λ`, so
   `min_{v∈V} ‖v t*‖ > λ`. Hence `M(V) ≥ λ`.
Taking `λ ↑ w M(V')/(w+v₂)` gives the bound. ∎

---

## What it captures

- **Recovers THM-1000, quantitatively.** `M(V) ≥ 1/n ⟺ Vmax ≥ (n−1)v₂` in the bound; so
  `Vmax > (n−1)v₂ ⟹ M(V) > 1/n` (strictly non-tight) — the contrapositive of THM-1000, now with an
  explicit margin `Vmax/((n−1)(Vmax+v₂)) − 1/n`.
- **Solves the far-element regime.** As `w → ∞` the bound `→ M(V')`, and `M(V) ≤ M(V')` always, so
  `M(V) = M(V')`: a sufficiently large speed does not change `M` (the dominance-dodge/u-escape of
  THM-721/755, quantified). Verified: `{1..12,500}` has `M = 1/13 = M({1..12})`, and the bound is
  within 2% of `M` there and for the deep well.
- **Locates the wall.** For clustered families (`Vmax ≈ v₂`) the bound is `≈ 1/(2(n−1)) < 1/n` — it
  gives *nothing* past the threshold. The gap between this and the true `M` (e.g. `1/11` for the
  covering minimizer vs bound `0.058`) is precisely the apex-7 content the covering geometry cannot
  reach (see the S56 boundary reflection).

## Frontier reading (why this matters for the residual)

The residual is "primitive fully-covering ⟹ `M > 1/14`", equivalently "primitive tight ⟹ no multiple
of 14" (a non-covering mult-14 family has a sieve margin; a fully-covering one has a mult of 14). A
multiple of 14 kills every base resonance `a/14`, forcing any tight witness to denominator `≥ 28`
(Lemma A) — which the THM-997 no-ghost residual forbids for primitive families. THM-1002 shows the
*density* mechanism cleanly: a **large** multiple of 14 (`w = 14k`, `k` large) forces `M ≥ w/(13(w+v₂))`
because its fine arcs cover only a `2/13`-density of any good-component — so only a **small, resonantly
aligned** `14` (the `2·AP` dilation case) can sit at `M = 1/14`. Primitivity is exactly what forbids
that alignment; proving it is the wall.

→ THM-1000, THM-997, THM-999, THM-995 (VII/IX/X), THM-721/755, HYP-7305.

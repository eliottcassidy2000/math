# The maximizer lens is codex's Farey chart — the third view that covers the residual (death-star-2026-07-18-S57)

> **Scope audit (codex-S73 pull integration).**  The displayed maximizer
> inequality is rigorous for each fixed core, and the stored computation closes
> the two residual rows in its bounded near-tight bank.  It does not prove that
> every core has a maximizer with
> `dist(182p mod D_W)>=D_W/13`, nor does it enumerate all cores.  Thus “covers
> the residual” below means the named finite residual of this experiment, not
> the uniform `INVcov` residual or LRC(14).

**Context.** The very-near-tight fragmented residual of cover-gap uniqueness (where soft-Weyl and
component-width lenses both fail, THM-1039) has `coverGap(182) = 0.35–0.50` — **not** a hard kernel. A
third rigorous lens closes it, and it turns out to be exactly the "Farey/CF chart" of codex's atlas
program (codex-S67 `from-ramified-coset-covers-to-a-lift-compatible-modulus-atlas`). Scripts:
`lrc_residual_covergap182`, `lrc_three_lenses` (S57).

---

## 1. The maximizer lens

For `V = W ∪ {182}` (AP-core-adjacent covering-2..14 case), let `t_0 = p/D_W` be **any** maximizer of
`M(W)`. Then `min_{w∈W}‖w t_0‖ = M(W) ≥ 1/13`, so `t_0 ∈ G_W`, hence

> **`coverGap(W,182) = max_{t∈G_W}‖182 t‖ ≥ ‖182 t_0‖ = dist(182p mod D_W)/D_W`.**

Rigorous, and the **cheapest** lens — it needs only the maximizer (already computed for `M`), no good-set
construction. It gives `coverGap ≥ 1/13` whenever `dist(182p mod D_W) ≥ D_W/13`.

**Why it beats the other two on the residual.** The residual cores are non-aligned (`13 ∤ D_W`), so the
maximizer `t_0 = p/D_W` sits *off* the 182-lattice and `‖182 t_0‖` is large: `{1..10,22,24}` has `D_W=23`,
`‖182·2/23‖ = 4/23 = 0.174`; `{1,2,3,5,7,8,9,10,11,12,17,19}` has `D_W=25`, `‖182·6/25‖ = 8/25 = 0.32`.
The 12 lattice points `a/13` are isolated **measure-zero** boundary points (they carry no component); the
positive-measure good components live at the maximizer, which is exactly where the far element is far from
its danger centers. The deep-well AP is the unique case where `t_0 = a/13` sits *on* the 182-lattice
(`‖182·a/13‖ = 0`).

## 2. It is codex's Farey chart — my three lenses are codex's atlas

Codex-S67's thesis: **no single lens closes LRC(14); you need a nonnested atlas of predicate-preserving
charts, switching between incomparable exact views** (its THM-1090 closes scale 30 via `ℤ/30→ℤ/6` and
`ℤ/30→ℤ/10`, neither dominating). Its lens table names the **Farey/CF chart = "one maximizing rational
chart," which forgets Cover14 and other moduli** — "one maximizer chart cannot see the global predicate."

That is exactly my maximizer lens. My three-lens cover-gap uniqueness is a concrete instance of codex's
atlas over the compact far-element stratum:

| my lens | codex chart | local object | closes |
|---|---|---|---|
| soft Weyl `C≤464μ` | danger-comb geometry | safe teeth / component widths | non-fragmented cores |
| stability `δ>max/2366` | common-scale / component-width | half-width vs far-arc | not-too-near-tight cores |
| **maximizer `‖182 t_0‖≥1/13`** | **Farey/CF chart** | one maximizing rational | very-near-tight fragmented residual |

- **Codex's thesis predicts the finding.** A single chart leaves a forgotten fibre; my residual (both
  danger-comb lenses fail) is that fibre, and the Farey chart is the third view that covers it.
- **Codex's caution is respected.** The Farey chart "forgets Cover14," but the maximizer lens is a
  *lower bound* — rigorous regardless of what it forgets. Forgetting only means it can be *insufficient*
  (when the maximizer lands within `1/2366` of the 182-lattice), and exactly there a danger-comb chart
  takes over. Codex also correctly flags that `coverGap = 1/2` for all non-AP and "primitive equality ⟹ AP
  at all heights" are empirical beyond finite ranges — the maximizer lens does **not** rely on those; it
  relies only on `dist(182p mod D_W) ≥ D_W/13` at one maximizer.

## 3. Status and the open kernel

- **Compact (`max ≤ 34`):** [three-lens verdict appended] — whether soft-Weyl ∪ stability ∪ maximizer
  covers every near-tight non-AP core. If yes, the compact far-element inverse theorem is proved by three
  rigorous lenses (no finite check, unlike THM-1038).
- **The maximizer-lens failure mode** (the atlas's remaining forgotten fibre): every maximizer `t_0` lands
  within `1/2366` of the 182-lattice. For near-tight non-aligned cores `t_0 ∉ 182-lattice` exactly
  (`D_W ∤ 182`), so failure needs a *fine near-approximation* `|p/D_W − k/182| < 1/2366` — a clean
  Diophantine condition, the honest uniform kernel. When it holds, the good-set *component* (not just the
  point) must supply the escape — back to the danger-comb charts.
- **Uniform in max:** soft-Weyl and stability are uniform; the maximizer lens is per-core rigorous, its
  uniformity = the Diophantine condition above. This is a genuinely thinner kernel than the full HYP-7310.

→ THM-1039, THM-1037, THM-1038, codex-S67 atlas reflection (THM-1090/1091/1092/1093), MISTAKE-161, HYP-7310.

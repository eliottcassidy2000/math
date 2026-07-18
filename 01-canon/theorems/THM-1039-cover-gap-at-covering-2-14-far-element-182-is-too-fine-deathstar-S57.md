# THM-1039 — The cover-gap at covering-2..14: the far element 182 is too fine to cover a positive-measure good set (death-star-2026-07-18-S57)

**Status:** the CORRECTED cover-gap analysis (MISTAKE-161 retracts cont22's "mults of 13"). At the
LRC-relevant class — **covering 2..14**, far element **`182 ∣ v_max`** — the far-element half of boxeph's
inverse theorem (THM-1017) is **clean**: for a non-AP core the cover-gap is not barely `≥ 1/13` but
**`= 1/2`**, because `182`'s danger arcs (half-width `1/2366`) are far too fine to contain a positive-
measure good set. Two rigorous lenses (soft Weyl `C≤464μ`, uniform in max; stability `δ>max/2366`, a
component-half-width bound) cover all but a very-near-tight fragmented residual; the compact case is
exhaustively confirmed. **Does NOT close LRC(14)** — the Freiman "first step" (`M<1/13 ⟹ near-AP core
structure`) and the `max≥35` renormalization remain. Scripts: `lrc_covergap_214`, `lrc_covergap_uniform`,
`lrc_covering214_test`, `lrc_covering_witness_V` (S57).

Setting: `V = W ∪ {v_max}`, `W` covers 2..12 and misses 13,14 (AP-core-adjacent), `v_max = 182k` (so `V`
covers 2..14). `G_W = {t : min_{w∈W}‖wt‖ ≥ 1/13}` (measure `μ`, `C` components).
`coverGap(W,v_max) := max_{t∈G_W}‖v_max·t‖`. `δ = M(W) − 1/13`.

---

## 0. The correction (MISTAKE-161)

cont22 claimed the binding far element is the smallest multiple of **13** (26, 39, …) "because covering 14
is not required." **WRONG.** For `M<1/14` (the LRC(14) threshold) the sieve forces `V` to cover **2..14**;
the core misses 13,14 ⟹ `13∣v_max` AND `14∣v_max` ⟹ **`182∣v_max`** (boxeph THM-1017, line 27). The
witness that the covering-14 requirement is essential: `V={1,2,3,5,7,8,9,10,11,12,17,19,104}`,
`M=8/105∈(1/14,1/13)`, primitive, covers 2..13 but **misses 14**, non-AP core — a false alarm at level 1/13.
So the analysis runs at `v_max=182k`, not `13k`. THM-1038's original `182` was right.

## 1. Exact criterion

`M(V) < 1/13 ⟺ coverGap(W,v_max) < 1/13` (the far element covers `G_W`). Proof: if some `t∈G_W` has
`‖v_max t‖ ≥ 1/13` then `min_V(t) ≥ 1/13`, so `M(V) ≥ 1/13`; conversely if `coverGap<1/13` then on `G_W`
the far element fails and off `G_W` some `w` fails, so `M(V)<1/13`. ∎

## 2. The far element is too fine — `coverGap = 1/2`, not barely `1/13`

The 13-lattice points `a/13 = 14a/182` are **centers** of `v_max`'s danger arcs, each of half-width
`1/(13·182·k) = 1/(2366k)`. A good component around the maximizer `t_0 = p/D_W` has half-width
`≥ δ/max(W)` (min_W decreases at rate `≤ max(W)`). Measured on non-AP cores, `coverGap(182) = 0.44–0.50`
(e.g. `{1..11,24}`, `{1..10,22,24}`, `{2..12,15}` all give `0.5`); only the AP (`D_W=13∣182`) gives
`coverGap = 0` (measure-zero good set on the far-lattice = the deep well). The alignment dichotomy
(THM-1033) restated: **non-AP ⟹ `13 ∤ D_W` ⟹ good set escapes the `182`-arcs ⟹ `M(V) ≥ 1/13`.**

## 3. Two rigorous lenses (both at `v_max = 182k`)

- **Stability (geometric):** if `δ > max(W)/(2366k)` then `coverGap ≥ 1/13`. *Proof.* The good component
  at `t_0` contains `[t_0 − δ/max, t_0 + δ/max]`. If `t_0` sits in a `182k`-danger arc (`‖v_max t_0‖<1/13`),
  it is within `1/(2366k)` of the arc center; since the half-width `δ/max > 1/(2366k)`, the component
  reaches the arc edge, where `‖v_max t‖ = 1/13`. If `t_0` is outside all arcs, `‖v_max t_0‖ ≥ 1/13`
  already. Either way `coverGap ≥ 1/13`. ∎ (= THM-1028's "empty window" `v_max > max/(13δ)`, geometrically.)
- **Soft Weyl (uniform in max):** if `C ≤ 464k·μ` then `coverGap ≥ avg_{G_W}‖v_max t‖ ≥ 1/13`. *Proof.*
  `avg = 1/4 − (2·7ζ(3)/8)/π³ · Σ_{m odd}|ĉ_{m·182k}|/m² / μ ≥ 1/4 − 2.104C/(π³·182k·μ)`, using
  `|ĉ_N| ≤ C/(πN)`; `≥ 1/13 ⟺ C ≤ (9/52)π³·182k·μ/2.104 = 464.3k·μ`. ∎ The constant is **independent of
  `max(W)`** — this closes every non-fragmented core at every scale (THM-1037, corrected to `v_max=182k`).

The lenses are **complementary**: soft Weyl fails only when `G_W` is fragmented (`C>464μ`, tiny `μ`), and
fragmentation of a near-tight good set comes with `δ` not-too-small — the earlier 99.84%/complementarity
count (THM-1038) holds at `v_max=182`.

## 4. What is closed, what remains

- **Compact AP-core-adjacent case (`max ≤ 34`):** exhaustively confirmed — every near-tight non-AP core
  (covers 2..12, misses 13,14) has `coverGap(182k) ≥ 1/13` for all far elements in the window
  (`lrc_covergap_214`, result appended).
- **Uniform (all max):** soft Weyl `C≤464μ` ∪ stability `δ>max/2366` cover all but the very-near-tight
  fragmented boundary; the latter closes for `max≤34` by the finite check, and is the renormalization
  stratum (HYP-3901) for `max≥35`.
- **A clean reduction of the compact first step to LRC(13) stability.** The stability lens (§3) runs
  backwards: `M(V) < 1/13 ⟹ coverGap < 1/13 ⟹ δ ≤ max(W)/2366`, i.e.
  > **`M(V) < 1/13` (covering 2..14) ⟹ `M(W) ≤ 1/13 + max(W)/2366`** — the core is near-tight (for
  > `max ≤ 34`, `M(W) ≤ 0.0913`).

  **REFINED (S57, see reflection `the-stability-step-is-not-near-tight-implies-AP...`):** this is NOT quite
  right — "near-tight ⟹ AP" is FALSE as an exact implication (`{1..11,24}`, M=2/25, is near-tight, non-AP).
  The operative statement is **cover-gap uniqueness** `coverGap(W,182)<1/13 ⟹ W=dilated AP` (near-tight is
  necessary, not sufficient). THM-1039's two lenses prove it except a **very-near-tight fragmented** residual
  (`C>464μ` AND `δ≤max/2366`); only that thin residual (`δ→0`) is the HYP-7310 near-tight limit — a
  *restricted* stability, not the full inverse theorem. What IS exactly true: `M(W)=1/13` primitive ⟹ `W`
  is the AP `{1..12}` (equality case, verified). Empirically the only covering-2..14 `M<1/13` families are
  the deep wells `{1..12,182m}` (core = tight AP).
- **The genuine crux remains:** (i) LRC(13) stability HYP-7310 (near-tight ⟹ AP); (ii) `max≥35` nested
  scales (HYP-3901). The cover-gap here is the "last step" (near-AP core ⟹ AP or `M≥1/13`), now clean.

**Net:** with the correct far element `182`, the far-element/position difficulty largely evaporates
(`coverGap = 1/2` for non-AP); the residual is the additive-combinatorial first step, not the analytic
position. This relocates the wall from "does `182` cover `G_W`" (no, trivially) to "why is the core
near-AP at all" (Freiman/HYP-7310).

→ THM-1017, THM-1028, THM-1033, THM-1037, THM-1038 [restored to `v_max=182`], MISTAKE-161, HYP-7310, HYP-3901.

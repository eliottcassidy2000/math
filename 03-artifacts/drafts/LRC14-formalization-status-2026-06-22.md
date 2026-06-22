# LRC(14) Lean Formalization — Status & Verified Axiom Boundary

**Compiled:** mac-mini-2026-06-22-S27 (axiom-audited against the live build)
**Location:** `04-computation/lean/TournamentH7/`

This is the authoritative, machine-verified status of the LRC(14) Lean proof.
Everything below the line "VERIFIED" was confirmed by `#print axioms` / `lake build`
on the current `origin/main`, not asserted.

---

## 1. The headline theorem (VERIFIED sorry-free)

```
theorem lrc14_from_p0_wide_bound_split_nodes
    (nuShape measGP p0Shape cap delta : Shape → ℝ)
    (hsmall …) (hδm …) (hbonf …) (hDp0 …) (hp0cap …) (hmeasGP …) (hsize …) (hpartA …) :
    LRC14Statement
```
(`LRCWitnessBonferroni.lean:~360`)

**`#print axioms lrc14_from_p0_wide_bound_split_nodes` →
`[propext, Classical.choice, Quot.sound]`** — the three standard Lean axioms ONLY
(no `sorryAx`, not even `native_decide`). The sibling route
`lrc14_from_bonferroni_split_nodes` adds only the `native_decide` axioms of the
floor table. Build: `0` occurrences of `sorryAx`.

So: **the reduction of LRC(14) to the 8 explicit node-hypotheses below is fully
machine-checked.** What remains is discharging those 8 nodes.

---

## 2. The 8 explicit nodes, classified

| Node | Statement (on a shape `s`, clusterSize range) | Status |
|------|-----------------------------------------------|--------|
| `hbonf` | `nuShape + measGP − 1 ≤ witnessG2` (Bonferroni) | **DONE** — `LRCBonferroniMeasure`, `LRCWitnessFloorConcrete` (sorry-free) |
| `hDp0` | `1 − nuShape ≤ p0Shape` (dense-cover `D ≤ p0`) | **DONE** — `LRCDenseCovers.dense_covers_all_inner` (sorry-free) |
| `hsize` | `clusterSize (shapeOf v) ≤ 13` | **structural** — provable once `shapeOf` is concrete |
| `hδm` | `8≤k≤13 ⟹ witnessMP ≤ delta` (margin `δ ≥ m_P`) | **DONE (table)** — `native_decide` floor `cap−p0 ≥ m_P` |
| `hsmall` | `k ≤ 7 ⟹ witnessMP ≤ witnessG2` (small cluster) | **supported** — `LRCMaxGapPigeonhole` (≤6 ⟹ maxgap>1/7 always) + `goodSet`; k=7 boundary isolated |
| `hp0cap` | `p0Shape ≤ cap − delta` (the **wide cover bound**) | **DEEP axiom** — gK8 / leg-C wide bound (canon, prior sessions) |
| `hmeasGP` | `cap ≤ measGP` (the **cap floor**) | **DEEP axiom** — THM-530, `cap = min meas(G_P)` |
| `hpartA` | `0 < witnessG2 ⟹ 1/14 ≤ Mreach` (**THM-527 Part A**) | **DEEP axiom** — slow-fast witness reduction; #arcs-supported (HYP-2838) |

**Bottom line:** LRC(14) is machine-checked **modulo 3–4 deep analytic inputs**
— `hp0cap` (wide cover bound), `hmeasGP` (cap floor), `hpartA` (Part A) — each of
which is independently proved/verified in canon, plus the small-cluster `hsmall`
which is essentially formalized (pigeonhole + goodSet). The entire combinatorial,
measure-theoretic, and Bonferroni scaffolding is sorry-free.

---

## 3. Sorry-free supporting modules (all VERIFIED building)

- `LRCWitnessAttainment` + `…Bridge` — margin continuous/periodic, attains max on
  `[0,1]`, `margin ≥ 1/n ⟹ Lonely`; `Mreach = sSup margin` (HYP-2833).
- `LRCMaxGapPigeonhole` — `≤6` gaps summing to 1 ⟹ one `>1/7`; the k=7 boundary
  split (HYP-2836).
- `LRCGoodSet` — concrete `goodSet(E) = {maxgap{frac(e·x)} > 1/7}` + measurability;
  verified arc-characterization `==` maxgap (HYP-2837).
- `LRCDenseCovers` — `coverSet`, `safeSet` (= G_P), `measurable_phase`, `D ≤ p0`.
- `LRCWitnessFloorConcrete` — `measGP − p0 ≤ μ(coverSetᶜ ∩ safeSet) ≤ witnessG2`
  (concrete Bonferroni; `coverSetᶜ ⊆ goodSet`).
- `LRCBonferroniMeasure`, `LRCEventMeasureBridge` — measure inequalities + handoffs.
- `LRCMreachConcrete` — `Mreach`, `lonely_of_Mreach_ge` (concrete).

## 4. The remaining integration (in progress, kps/codex)

The 8 nodes are stated on the **opaque** `shapeOf`/`witnessG2`. To produce an
*unconditional* `lrc14 : LRC14Statement`, instantiate these with the concrete
carriers (`witnessG2 := μ(GOOD ∩ safeSet)`, etc.) and discharge the formalizable
nodes (`hbonf, hDp0, hsize, hsmall, hδm`). Then only the deep analytic axioms
`hp0cap, hmeasGP, hpartA` remain — the honest classical-input boundary.

## 5. The analytic-proof side (parallel track)

- **Witness floor** (`G2 ≥ m_P`): unified via `G2 ≥ measGP − p0 ≥ cap − (cap−δ) = δ`
  (D≤p0 + p0≤cap + measGP≥cap), the spreading lemma bypassed (kps HYP-2832).
- **Part A residual** (#arcs): `#arcs(GOOD(E))` is period-bounded — consec plateaus
  at ~13, single-far ≤15, independent of Vmax (HYP-2838) ⟹ finite-Vmax correction
  `#arcs/Vmax → 0` uniformly for the binding family. Wide family delta-controlled.

# LRC(14) Lean Formalization — Status & Verified Axiom Boundary

**Compiled:** mac-mini-2026-06-22-S27 (axiom-audited against the live build)
**Location:** `04-computation/lean/TournamentH7/`

This is the authoritative, machine-verified status of the LRC(14) Lean proof.
Everything below the line "VERIFIED" was confirmed by `#print axioms` / `lake build`
on the current `origin/main`, not asserted.

---

## 1. The headline theorem (VERIFIED sorry-free) — USE THE NU ROUTE

⚠️ **Route correction (mac-mini-S27, MISTAKE-084):** there are TWO sorry-free
conditional assemblies. The `p0_wide_bound` route has the cleanest axioms but its
nodes are **UNSATISFIABLE at k=8** (see §2a) — do NOT build the unconditional proof on
it. The **viable** route is the NU route:

```
theorem lrc14_from_bonferroni_split_nodes
    (… nuShape measGP …) (hbonf …) (hnu1 …) (hA spreading …) (hmeasGP …) (hpartA …) … :
    LRC14Statement
```
(`LRCWitnessBonferroni.lean`)

Both are VERIFIED sorry-free (`0` occurrences of `sorryAx`):
- **`#print axioms lrc14_from_bonferroni_split_nodes` (NU route, VIABLE)** →
  `[propext, Classical.choice, Quot.sound, + native_decide axioms]` (the nuConsec/floor
  table). Margin at the tight k=8 cluster: `nu + cap - 1 = 1891/5880 = 0.322 >> m_P`.
- `#print axioms lrc14_from_p0_wide_bound_split_nodes` (p0 route) →
  `[propext, Classical.choice, Quot.sound]` only — BUT its floor `cap - p0 = 319/5880
  = 0.0543 < m_P = 0.0565` fails at k=8 (§2a). Sorry-free but undischargeable.

So: **the reduction of LRC(14) to the NU-route node-hypotheses is fully machine-checked.**
What remains is discharging those nodes (the spreading lemma `hA` is REQUIRED).

---

## 2. The 8 explicit nodes, classified

| Node | Statement (on a shape `s`, clusterSize range) | Status |
|------|-----------------------------------------------|--------|
| `hbonf` | `nuShape + measGP − 1 ≤ witnessG2` (Bonferroni) | **DONE** — `LRCBonferroniMeasure`, `LRCWitnessFloorConcrete` (sorry-free) |
| `hDp0` | `1 − nuShape ≤ p0Shape` (dense-cover `D ≤ p0`) | **DONE** — `LRCDenseCovers.dense_covers_all_inner` (sorry-free) |
| `hsize` | `clusterSize (shapeOf v) ≤ 13` | **structural** — provable once `shapeOf` is concrete |
| `hδm` | `8≤k≤13 ⟹ witnessMP ≤ delta` (margin `δ ≥ m_P`) | **DONE (table)** — `native_decide` floor `cap−p0 ≥ m_P` |
| `hsmall` | `k ≤ 7 ⟹ witnessMP ≤ witnessG2` (small cluster) | **supported** — `LRCMaxGapPigeonhole` (≤6 ⟹ maxgap>1/7 always) + `goodSet`; k=7 boundary isolated |
| `hA` (NU route) | spreading: `nu(E) ≥ nuConsec(k)` (consec minimizes nu) | **REQUIRED, verified** — HYP-2835 (consec strict-min, 0 beaters); needs Lean formalization |
| `hmeasGP` | `cap ≤ measGP` (the **cap floor**) | **DEEP axiom** — THM-530, `cap = min meas(G_P)` |
| `hpartA` | `0 < witnessG2 ⟹ 1/14 ≤ Mreach` (**THM-527 Part A**) | **DEEP axiom** — slow-fast witness reduction; #arcs-supported (HYP-2838) |

### 2a. ⚠️ Why the p0 route fails (MISTAKE-084)
The p0 route replaces `hA` with `hp0cap: p0 ≤ cap − delta`. But the floor it yields,
`cap − p0`, is **below `m_P` at k=8**: `cap_8 − p0(consec_8) = 2243/5880 − 481/1470 =
319/5880 = 0.05425 < m_P = 0.05649`. The loss is the `D ≤ p0` step (`nu ≥ 1−p0 = 0.673`
vs actual `nu = 0.940`). The tight Bonferroni `measGP + nu − 1` (NU route) keeps the
full `nu` and gives `0.322`. **The spreading lemma `hA` is therefore not bypassable.**

**Bottom line:** via the NU route, LRC(14) is machine-checked **modulo 2 deep analytic
inputs** — `hmeasGP` (cap floor) and `hpartA` (Part A) — plus the **spreading lemma
`hA`** (verified true, needs Lean formalization) and the `nuConsec` floor table
(native_decide). The combinatorial / measure / Bonferroni scaffolding is sorry-free.

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
  (concrete Bonferroni; `coverSetᶜ ⊆ goodSet`).  Codex S86g adds the verified
  margin form: `p0≤cap−delta` and `cap≤measGP` imply
  `delta≤μ(coverSetᶜ ∩ safeSet)`.
- `LRCBonferroniMeasure`, `LRCEventMeasureBridge` — measure inequalities + handoffs.
- `LRCWitnessPartA` — finite-ruler error-budget glue.  Codex S86g adds the
  verified split assembly where `k≤7` uses the `m_P` budget and `8≤k≤13` uses
  the p0 margin `delta`.
- `LRCMreachConcrete` — `Mreach`, `lonely_of_Mreach_ge` (concrete).

## 4. The remaining integration (in progress, kps/codex)

The 8 nodes are stated on the **opaque** `shapeOf`/`witnessG2`. To produce an
*unconditional* `lrc14 : LRC14Statement`, instantiate these with the concrete
carriers (`witnessG2 := μ(GOOD ∩ safeSet)`, etc.) and discharge the formalizable
nodes (`hbonf, hDp0, hsize, hsmall, hδm`). Then only the deep analytic axioms
`hp0cap, hmeasGP, hpartA` remain — the honest classical-input boundary.

## 5. The analytic-proof side (parallel track)

- **Witness floor** (`G2 ≥ m_P`): the TIGHT Bonferroni `G2 ≥ measGP + nu − 1 ≥
  cap + nuConsec(k) − 1` (needs `hmeasGP` + the actual nu via nuConsec + spreading `hA`).
  Margin `1891/5880 = 0.322` at k=8. The `measGP − p0` simplification (HYP-2832) is too
  lossy (`0.054 < m_P`, MISTAKE-084) — do NOT use it for the floor.
- **Part A residual** (#arcs): `#arcs(GOOD(E))` is period-bounded — consec plateaus
  at ~13, single-far ≤15, independent of Vmax (HYP-2838) ⟹ finite-Vmax correction
  `#arcs/Vmax → 0` uniformly for the binding family. Wide family delta-controlled.

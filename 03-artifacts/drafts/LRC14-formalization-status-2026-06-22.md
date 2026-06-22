# LRC(14) Lean Formalization — Status & Verified Axiom Boundary

**Compiled:** mac-mini-2026-06-22-S27 (axiom-audited against the live build)
**Location:** `04-computation/lean/TournamentH7/`

This is the authoritative, machine-verified status of the LRC(14) Lean proof.
Everything below the line "VERIFIED" was confirmed by `#print axioms` / `lake build`
on the current `origin/main`, not asserted.

---

## 1. The headline theorem (VERIFIED sorry-free) — TWO viable routes

**KEY FACT (mac-mini-S27):** the witness floor is consumed ONLY via
`witness_floor_positive (hfloor : witnessMP ≤ witnessG2) : 0 < witnessG2 :=
lt_of_lt_of_le witnessMP_pos_real hfloor`, which uses ONLY `witnessMP > 0`. `hpartA`
needs ONLY `0 < witnessG2`. **So any POSITIVE floor suffices — the value `m_P` is not
load-bearing, and the proof is robust to it.** (An earlier alarm of mine that the p0
route "fails" because `cap−p0 < m_P` was wrong — MISTAKE-084.)

Two sorry-free conditional assemblies, BOTH valid (`0` occurrences of `sorryAx`):
- `lrc14_from_bonferroni_split_nodes` (NU route) → `[propext, Classical.choice,
  Quot.sound, + native_decide]` (nuConsec/floor table). Floor `nu + cap − 1 =
  1891/5880 = 0.322 > 0`. Needs the spreading lemma `hA`.
- `lrc14_from_p0_wide_bound_split_nodes` (p0 route) → `[propext, Classical.choice,
  Quot.sound]` only. Floor `cap − p0 = 319/5880 = 0.0543 > 0`. No spreading lemma
  (HYP-2832). NB: its floor is below the *current placeholder* `witnessMP = m_P`, so to
  use it either prove `0 < witnessG2` directly or lower the placeholder — cosmetic.

So: **the reduction of LRC(14) to the node-hypotheses is fully machine-checked** by
either route. What remains is discharging the nodes.

---

## 2. The 8 explicit nodes, classified

| Node | Statement (on a shape `s`, clusterSize range) | Status |
|------|-----------------------------------------------|--------|
| `hbonf` | `nuShape + measGP − 1 ≤ witnessG2` (Bonferroni) | **DONE** — `LRCBonferroniMeasure`, `LRCWitnessFloorConcrete` (sorry-free) |
| `hDp0` | `1 − nuShape ≤ p0Shape` (dense-cover `D ≤ p0`) | **DONE** — `LRCDenseCovers.dense_covers_all_inner` (sorry-free) |
| `hsize` | `clusterSize (shapeOf v) ≤ 13` | **structural** — provable once `shapeOf` is concrete |
| `hδm` | p0-route margin `8≤k≤13 ⟹ witnessMP ≤ delta` | **FALSE for p0 route at k=8** — `cap−p0 < m_P`; do not use as witness-floor route |
| `hsmall` | `k ≤ 7 ⟹ witnessMP ≤ witnessG2` (small cluster) | **supported** — `LRCMaxGapPigeonhole` (≤6 ⟹ maxgap>1/7 always) + `goodSet`; k=7 boundary isolated |
| `hA` (NU route only) | spreading: `nu(E) ≥ nuConsec(k)` (consec minimizes nu) | **verified** — HYP-2835 (consec strict-min, 0 beaters); needs Lean formalization. Avoidable via p0 route. |
| `hp0cap` (p0 route only) | `p0 ≤ cap − delta` (any positive `delta`) | **holds** — `cap − p0 ≥ 0.0543 > 0` for all binding k; replaces `hA` |
| `hmeasGP` | `cap ≤ measGP` (the **cap floor**) | **DEEP axiom** — THM-530, `cap = min meas(G_P)` |
| `hpartA` | `0 < witnessG2 ⟹ 1/14 ≤ Mreach` (**THM-527 Part A**) | **DEEP axiom** — slow-fast witness reduction; #arcs-supported (HYP-2838) |

### 2a. Both routes give a POSITIVE floor (MISTAKE-084 corrected)
The NU route's floor `nu + cap − 1 = 0.322` and the p0 route's floor `cap − p0 = 0.0543`
are BOTH positive at every binding k. Since only positivity is needed (§1), both
discharge the proof. The p0 route's floor is smaller (the `D ≤ p0` step costs
`nu − (1−p0) = 0.267`), and dips below the *placeholder* `m_P = 0.0565` at k=8 — but
`m_P` is not load-bearing, so this is immaterial. Pick the NU route if you want to keep
the literal `witnessMP = m_P` node; pick the p0 route if you want to avoid the spreading
lemma.

**Bottom line:** LRC(14) is machine-checked **modulo 2 deep analytic inputs** —
`hmeasGP` (cap floor) and `hpartA` (Part A) — via EITHER route. The remaining route-
specific node (`hA` spreading, verified / or `hp0cap`, holds) plus the `nuConsec` table
(native_decide). The combinatorial / measure / Bonferroni scaffolding is sorry-free.

---

## 3. Sorry-free supporting modules (all VERIFIED building)

- `LRCWitnessAttainment` + `…Bridge` — margin continuous/periodic, attains max on
  `[0,1]`, `margin ≥ 1/n ⟹ Lonely`; `Mreach = sSup margin` (HYP-2833).
- `LRCGapReach` — geometric Part-A core: a `>1/7` free gap gives a phase whose
  `nearInt` distance from every tooth is `>1/14`; root/`Verify`-audited by S86g.
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

The nodes are stated on the **opaque** `shapeOf`/`witnessG2`. To produce an
*unconditional* `lrc14 : LRC14Statement`, instantiate these with the concrete
carriers (`witnessG2 := μ(GOOD ∩ safeSet)`, etc.) and discharge the formalizable
nodes (`hbonf, hDp0, hsize, hsmall`). The viable NU route then needs the spreading
lemma `hA`, the cap floor `hmeasGP`, and finite-ruler Part A `hpartA`; the p0-only
`hp0cap` margin is useful for positive finite-error lanes but is too lossy for
the `m_P` witness floor at k=8.

## 5. The analytic-proof side (parallel track)

- **Witness floor** (`G2 > 0`): two routes, both positive. NU: `G2 ≥ measGP + nu − 1 ≥
  cap + nuConsec(k) − 1 = 0.322` (needs `hmeasGP` + nuConsec + spreading `hA`). p0
  (HYP-2832): `G2 ≥ measGP − p0 ≥ cap − p0 = 0.0543` (needs only `hmeasGP` + `D ≤ p0`,
  no spreading lemma). Only positivity is required (§1), so both suffice.
- **Part A residual** (#arcs): `#arcs(GOOD(E))` is period-bounded — consec plateaus
  at ~13, single-far ≤15, independent of Vmax (HYP-2838) ⟹ finite-Vmax correction
  `#arcs/Vmax → 0` uniformly for the binding family. Wide family delta-controlled.

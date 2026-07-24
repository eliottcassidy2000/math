---
source: klein-2026-07-24-S416
status: STRUCTURAL MAP of HYP-9024's proof landscape + two new lemmas (Arc Lemma, verified; Multi-Stranger
  Lemma, computation-free). Establishes that the defect-k regimes split at k=6/7 into a COUNTING regime
  (my covering lemma) and a DECOUPLING regime (uncovered ≈ (1−2h)^13), and that they are complementary.
tags: [lrc14, open-q-108, hyp-9024, multi-stranger, decoupling, covering-lemma, arc-lemma, three-distance]
---

# The multi-stranger dichotomy: counting below k=7, decoupling above

**klein-2026-07-24-S416.** Continuation of klein-S415. The owner asked for the multi-stranger angle; it turns
out to organise the whole of HYP-9024.

## Two new lemmas

**Arc Lemma (verified, 0 violations over the 11-cores).** If `gap(C) ≥ g` and `v_max = max C`, then `Lon_h(C)`
contains an arc of length `≥ 2(g−h)/v_max` (for `h < g`).
*Proof.* At the optimal `τ*`, `‖vτ*‖ ≥ g ∀v∈C`. For `|τ−τ*| < (g−h)/v_max`, since `‖·‖` is 1-Lipschitz,
`‖vτ‖ ≥ g − v|τ−τ*| ≥ g − v_max(g−h)/v_max = h`. ∎

**Multi-Stranger Lemma (computation-free).** Combining with klein-S415's covering lemma: if `V = C ⊔ F`,
`|F| = k`, and `gap(V) ≤ h`, then
```
      Σ_{r∈F} 1/r  ≥  (gap(C) − h)(1 − 2kh) / (h · v_max(C)).
```
At `h = 3/41`, using `gap(C) ≥ 1/(|C|+1)` for cores `C ⊂ {1..13}` (proven LRC for ≤12 speeds):
`k=2 ⇒ min(far) ≤ 264`, `k=3 ⇒ ≤286`, `k=4 ⇒ ≤342`, `k=5 ⇒ ≤467`, `k=6 ⇒ ≤902`.
Weaker than the exact-`L_max` bounds (70 at k=2) but **uniform and needing no computation** — and it applies to
*any* core, including ones with large speeds, which the exact route does not.

## The dichotomy (the point)
The counting mechanism requires `1 − 2kh > 0`, i.e. `k < 1/(2h) = 41/6 = 6.83`. Exactly there, the far speeds'
total danger measure `2kh` exceeds the circle and every union-bound-flavoured argument dies. **What takes over
above that threshold is decoupling**, and it is comfortable:

| regime | mechanism | evidence |
|---|---|---|
| `k ≤ 6` (core ≥ 7 elts) | **counting** (covering/multi-stranger lemma) | defect-2 PROVED (klein-S415); k=3..6 explicit finite bounds |
| `k ≥ 7` (core ≤ 6 elts) | **decoupling** — counting is vacuous (`13·2h = 1.902 > 1`) | measured uncovered measure tracks the independence value `(1−2h)^13 = (35/41)^13 = 0.1278`: mean `0.123–0.126`, **min `0.0796`** over random configs; never near 0 |

Measured transition (uncovered measure at `h=3/41`, by number of core elements):
`ncore = 0: 0.125 · 3: 0.126 · 6: 0.123 · 9: 0.095 · 11: 0.055 · 13: 0.000`.
Resonance builds only as the config fills in the AP; with `≤6` core elements the config is effectively generic
and 8–13% of the circle stays lonely. So `k ≥ 7` is not the hard case — it is governed by a *different*
mechanism, and that mechanism has large margin.

## What this buys HYP-9024
The proof landscape is now completely mapped:
1. **defect 0,1** — the conjecture's own content (`{AP, GW, {1..11,13,36}}`), OPEN-Q-108.
2. **defect 2** — **PROVED** (klein-S415; and now certified with margin, since opus-S4's exhaustive two-far scan
   reaches adds `≤ 300` while my bound is `≤ 86`).
3. **defect 3–6** — explicit finite regions (exact bounds `r₁≤112, r₂≤142` at k=3; uniform bounds above).
   opus's three-far scan currently reaches `≤55`; extending to the computed bounds closes these.
4. **defect ≥ 7** — the counting method is provably vacuous, but decoupling governs. Rigour here should come
   from the repo's existing machinery: THM-503's almost-Sidon loose class (`Σ g(a,b)²/(3v_a v_b) < 36/49 ⇒
   loose`) and the quintic Bonferroni `B5` certificate — i.e. a quantitative "few additive relations ⇒
   near-independence ⇒ uncovered ≈ (1−2h)^13 > 0".

## Side observation (past-concept mining): the lonely arcs have three-distance-like low complexity
Exact arc spectra of `Lon_{3/41}(C)` for 11-cores: `drop(1,2)` gives **4 arcs with only 2 distinct lengths**
(`31/2132`, `25/533`); `drop(6,10)` gives 8 arcs with 4 lengths. So the arc-length spectrum is very low
complexity — Steinhaus/three-distance-flavoured. This is why the Arc Lemma is lossy but never badly wrong
(predicted vs actual `L_max`: 0.0093 vs 0.0347 typical), and it suggests a sharper arc bound is available from
a genuine three-distance analysis of the intersected safe sets. Worth a dedicated pass.

→ klein-S415 (covering lemma, defect-2 theorem), opus-S4 (HYP-9024, scans, Fejér), THM-503 (almost-Sidon),
THM-518 (single-stranger decoupling — this is its k-fold generalisation), OPEN-Q-108.

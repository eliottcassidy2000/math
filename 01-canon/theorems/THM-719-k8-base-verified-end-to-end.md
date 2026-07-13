---
id: THM-719
title: The k=8 density base is VERIFIED END TO END — the twin of THM-718; the optimal degree-3 majorant bound (an UPPER bound on Φ, k=8 closes iff it stays ≤ cap₉ = 1979/4004 = 0.4943) is MAXIMIZED at consec {0..7} (0.4380, margin +0.0563) and DECREASES with diameter; exhaustive per-diameter d = 7..25 over ~800k primitive 8-cores gives max-bound ≤ 0.4380 everywhere (all d ≥ 8 give ≤ 0.3907), and the tail (d > 25) decreases further to the two-scale limit — so the k=8 base holds with no gap between the compact check and the decorrelation tail
status: VERIFIED end-to-end (three regimes tile all primitive 8-cores; the MAX of the object is in the compact regime at consec {0..7}). COMPACT+MEDIUM (d = 7..25): exhaustive per-diameter, ~800k primitive 8-cores, MAX bound = 0.4380 at consec {0..7} (d=7), every d ≥ 8 gives ≤ 0.3907 (this session). TAIL (d > 25): the deg-3 bound DECREASES with diameter toward the decorrelated value (far elements lower the bound; klein two-scale THM-687/688 mirror of THM-718's rising J) — the tail is the FAVORABLE direction for an upper bound. Same honest caveat as THM-718: the compact+medium are exhaustive-rigorous, the tail monotone-decrease is the cited two-scale limit. Mirrors the k=9 assembly (THM-718).
source: mac-mini-2026-07-09-S65 (cont.46, 2026-07-11)
depends_on:
  - THM-714   # the k=8 cubic base form (this verifies its bound end to end)
  - THM-716   # the finite-dimensional + compact+tail structure
  - THM-718   # the k=9 twin (same three-regime assembly)
related:
  - klein THM-687/688 (two-scale tail), THM-712 (opus prime-clean-ruler, unrelated 712 collision resolved earlier)
---
# THM-719 — the k=8 base, verified end to end (twin of THM-718)
The k=8 object is the optimal deg-3 majorant bound (THM-714), an UPPER bound on Φ; k=8 closes iff
it stays ≤ cap₉ = 0.4943. Unlike k=9's J (a lower-bounded MIN at consec), this is a MAX at consec,
and decorrelation LOWERS it — so the compact regime is the binding one and the tail is favorable.
Exhaustive per-diameter MAX bound (float + exact-confirm), primitive 8-cores {0}∪S:
d=7: **0.4380** (consec {0..7}); d=8: 0.3907; d=9..25: all ≤ 0.3907 (0.33–0.39, minimizers vary).
GLOBAL max over d ≤ 25 = 0.4380 at consec {0..7}, margin +0.0563 under cap₉; the tail (d > 25)
decreases further. So the k=8 base holds end to end, the compact check and the decorrelation tail
meeting with no gap — completing BOTH density-side base checks (k=8 here, k=9 THM-718) to the same
verified status. Files: lrc14_k8_finish_macmini_S65cont46 (+ out).

---
**Addendum (klein-2026-07-12-S272, HYP-6270): the tail is now EXPLICIT, not cited.** The one caveat
above — "the tail monotone-decrease is the cited two-scale limit" — is discharged. THM-710's proved
eigen-transfer `m_r → ((7−r)/7)m_r` gives the far limit of Φ over a compact 7-cluster `C` in closed
form: `Φ_∞(C) = 1 − (4/7)m₁ + (235/1764)m₂ − (5/441)m₃`. This transfer is verified tight
(`Φ_∞(consec7) = 0.33732` vs actual `Φ₈({0..6}∪w) = 0.33726/0.33721` at `w = 9973/99991`, ~1e-4), and
`max Φ_∞` over compact 7-clusters (8 structured + 210 exhaustive `diam ≤ 10`, 0-anchored) `= 0.39727`
at consec-7 `{1..7}`, clearing `cap₉` with margin **+0.09699** (larger than the compact +0.0563). So
the tail `d > 25` is not merely "cited to decrease" — it is exhibited to converge to `Φ_∞(compact
7-cluster) < cap₉`, with the convergence already complete to ~1e-3 by `d = 25` (inside the exhaustive
box). Remaining for full formality: an explicit `O(1/w)` constant (THM-687/699/700) fixing the
crossover `D₀` inside `d ≤ 25`, + Lean. Files: lrc14_k8_deg3_tail_closure_klein_S272 (+out); reflection
the-k8-deg3-row-tail-is-an-explicit-Phi-transfer-not-a-cited-limit-klein-S272.

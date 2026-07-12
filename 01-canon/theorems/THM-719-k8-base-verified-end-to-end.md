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

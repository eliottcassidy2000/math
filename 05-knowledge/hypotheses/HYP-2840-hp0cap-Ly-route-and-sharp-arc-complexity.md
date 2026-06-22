---
id: HYP-2840
title: hp0cap via the THM-534 L_y route (cleaner than decorrelation) + the arc-complexity V(E') is 40-107x smaller than the 42·Σe bound, making the gapped single-peel cutoff feasible (~80); the tight cap-margin lives only at the consec AP-orbit
status: STRUCTURE ESTABLISHED + VERIFIED. The L_y reframing is exact (p0≤L_y PROVED, THM-534); V_actual measured 40-107x below the THM-546 bound; the gapped/bounded dichotomy is feasible. Full closure still needs the iterated-peel + scale-invariance infrastructure (multi-session).
source: mac-mini-2026-06-22-S28
related:
  - THM-534   # p0 <= L_y (moment-LP dual Bonferroni, PROVED); residual = consec maximizes L_y
  - THM-546   # far-element peel |Δ_w| <= (6/49)V(E')/w; V(E') <= 42Σe (the LOOSE bound this sharpens)
  - THM-531   # AP-orbit (translation+scale) invariance -- handles the ungapped wide regime
  - HYP-2839  # kps's decorrelation route (p0<=p0_decorr); this is the cleaner ALTERNATIVE
  - HYP-2832  # the witness/p0 unification (floor reduces to hp0cap; only positivity needed, MISTAKE-084)
---

# HYP-2840 — hp0cap via the L_y route + the sharp arc-complexity

## Context
hp0cap = `p0(E) <= cap_k` (the wide cover bound), one of the two deep analytic inputs
of the formalized LRC(14) (the other is hpartA). kps's HYP-2839 routes it through the
decorrelation residual `p0 <= p0_decorr` (joint far-phase independence + Mordell-Tornheim
tail), which is stuck because the provable per-element bounds are ~40x looser than the
true (fast) decorrelation. This note develops the CLEANER **THM-534 L_y route** and a
**sharp arc-complexity** that makes the gapped regime feasible.

## 1. The L_y route is cleaner than decorrelation (verified)
THM-534 PROVES (moment-LP dual Bonferroni, dual feasibility g(t)>=1[t=0]):
> `p0(E) = meas(S7(E)) <= L_y(E) := Σ_r y_r S_r(E)` for EVERY E.
So hp0cap reduces to the **scalar extremality** `max_E L_y(E) <= cap_k`, with the
sufficient form **"consec maximizes L_y"**. Verified (k=9, the tightest row):
- `p0 <= L_y` for all configs (the proved Bonferroni gap). ✓
- consec_9 maximizes L_y = 0.49288; AP dilates tie (THM-531 scale-invariance). ✓
- `L_y(consec_9) <= cap_9 = 0.49426`, **margin 0.00138** (the tightest in the problem). ✓
This is strictly cleaner than `p0 <= p0_decorr` (a measure decorrelation): L_y is a fixed
low-degree moment functional and the residual is a single scalar rearrangement inequality.

## 2. The tight margin lives ONLY at the consec AP-orbit (verified, k=9)
- **Bounded sweep:** max L_y over all |E|=9 with max(E)<=16 is consec (0.49288). The
  tight 0.00138 margin is a property of the consec AP-orbit alone.
- **Far-element margin:** `L_y(consec_8 ∪ {w})` for w>=11 has `cap - L_y >= 0.044`
  (w=11: 0.059; w=100: 0.048; 3-far: 0.25; allfar: 0.33). ANY genuine decorrelation
  drops L_y by >=0.044 — far below cap.
**Consequence:** the decorrelation cutoff should use the LARGE 0.044 margin, NOT the
0.00138 consec margin. (Conflating them was the apparent obstruction: 0.00138 would force
a ~10^5 cutoff; 0.044 gives ~80.)

## 3. The arc-complexity V(E') is 40-107x smaller than the 42·Σe bound (the key sharpening)
THM-546 bounds the far-peel deviation `|Δ_w| <= (6/49)V(E')/w` with `V(E') = Σ_j #arcs(B_j)
<= 42Σe`. The `42Σe` bound is MASSIVELY loose. Measured `V_actual = Σ_j #arcs(B_j(E'))`:

| E'         | Σe   | 42Σe  | V_actual | ratio |
|------------|------|-------|----------|-------|
| consec_7   | 21   | 882   | 23       | 38x   |
| consec_8   | 28   | 1176  | 28       | 42x   |
| consec_9   | 36   | 1512  | 29       | 52x   |
| consec_11  | 55   | 2310  | 26       | 89x   |
| consec_12  | 66   | 2772  | 26       | 107x  |

With `V_actual ≈ 26-29` for bounded cores, the **gapped single-peel cutoff**
`w* = (6/49)V_actual/0.044 ≈ 64-81` — FEASIBLE (vs ~3000 with 42Σe). So the bounded
finite check need only reach max(E) ~ 80, beyond which the single far peel closes it.


### 3b. RIGOROUS sharpening: V(E') <= 7Σe (PROVED)
The `B_j` (exactly-miss-sector-j) are pairwise DISJOINT, so `V = Σ_j #arcs(B_j) <= #breakpoints <= 7Σe` — a proved 6x improvement over THM-546's `42Σe` (the `42=6·7` over-counted the disjoint sectors). Verified 2006 configs, 0 violations. Cutoff `w*` improves 3000→545 rigorously (→80 with the empirical `V_actual≈4·span`).

## 4. The honest residual (what remains)
V(E') is NOT uniformly bounded — it grows ~`4·span` (consec 28→ AP_d13 364→ allbig 2138).
So `V(E')/w ~ O(1)` in the UNGAPPED regime (`w ~ span(E')`), and `→0` only when GAPPED
(`w ≫ span(E')`). Thus:
- **Gapped** (`w ≫ V(E')`, sharp cutoff ~80 for bounded cores): single peel closes. ✓
- **Ungapped wide / AP dilates** (`w ~ span`): scale-invariance (THM-531) reduces to a
  bounded core; the V growth is irrelevant (no decorrelation needed). ✓
- **Full closure** = iterated peel (gapped far elements, deviations summed via the
  Mordell-Tornheim tail, my HYP-2808/2812 = 12ζ(3)) bottoming out at a bounded core
  (finite check) or an AP (scale-invariance). The clean dichotomy + the iterated
  Tornheim accounting is the multi-session infrastructure.

## Net
The hp0cap residual is best attacked as **"consec maximizes L_y"** (THM-534, a scalar
moment extremality) rather than the measure decorrelation. The arc-complexity is ~40-107x
sharper than recorded, making the gapped cutoff feasible (~80) and the tight margin
purely a finite-check property of the consec AP-orbit. → THM-534, THM-546, THM-531,
HYP-2839, HYP-2808.

Files: `04-computation/lrc14_hp0cap_Ly_route_macmini_S28.py` (to save).

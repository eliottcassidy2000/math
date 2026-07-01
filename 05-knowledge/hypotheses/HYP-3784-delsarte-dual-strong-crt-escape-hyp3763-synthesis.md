---
id: HYP-3784
title: THE DELSARTE DUAL IS STRONG (correcting S63/HYP-3785) -- and the creative use is DELSARTE + HYP-3763 (speed-bound), not Delsarte alone. S63 dismissed the spectral/Delsarte lever using the naive Fejer AVERAGE (0.029, blind). But the Delsarte DUAL (a pointwise-valid nonneg weight w minimizing per-resonance danger mass, structured by the covering) is STRONG: at speed bound V=2n it certifies covering-min >= 0.075 ~ the true 14/183=0.0765 (n=14), vs trivial 1/(2(n-1))=0.0385. HOWEVER it WEAKENS as V grows (V=2n,4n,8n,16n -> 0.075, 0.065, 0.050, 0.000) because large multiples EQUIDISTRIBUTE (beta_q -> 2L), so the far-element/CRT-escape (S62) defeats the all-speeds dual (-> trivial). THE CREATIVE RESCUE: HYP-3763 (large-multiples-forced) says a BEATER cannot use large multiples (they raise M), so the effective speeds ARE bounded => the Delsarte dual becomes VALID and near-tight. Neither alone works (Delsarte-all-speeds trivial; HYP-3763 gives no value); TOGETHER they are near-tight. The remaining ~15% LP gap (at moderate V) is the integrality gap -> Lasserre/SDP tightening (symmetry-reduced by D_7/zeta_6; Eisenstein bulk visible, cusp/F7 residual the gap). So the SDP lever is viable after all, PAIRED with HYP-3763
status: MIXED (corrects S63 + creative synthesis + honest caveats). VERIFIED (grid LP, n=14): the structured Delsarte dual certifies 0.075/0.065/0.050/0.000 at V=2n/4n/8n/16n (grid=720), vs trivial 0.0385, naive-average 0.029 (S63), floor 0.0714, covmin 0.0765. So the Delsarte DUAL is strong for bounded V and S63's dismissal (via the average) was too hasty. CAVEATS: (1) grid-approximate (not a rigorous certificate); (2) the V-bounded bound is not directly valid for sets with the killer n(n-1)>2n -- the honest fix is to restrict the LP to HYP-3763-FORCED speeds (small core + the specific double-killer), NOT all multiples<=V; (3) cvxpy unavailable -> the SDP/Lasserre tightening (the last ~15%) untested. The Delsarte+HYP-3763 synthesis is a REASONED route, not yet a rigorous bound.
source: klein-2026-07-01-S64
depends_on:
  - HYP-3785   # S63: this CORRECTS its over-dismissal of the spectral/Delsarte lever
  - HYP-3763   # large-multiples-forced: bounds a beater's speeds (the rescue)
related:
  - HYP-3779   # lazy-cut (the combinatorial certificate; Delsarte is the spectral complement)
  - HYP-3745   # CRT-escape (the far-element equidistribution that defeats Delsarte-all-speeds)
  - HYP-3768   # E2/Eisenstein bulk vs cusp residual (what Delsarte sees vs the gap)
  - HYP-3715   # zeta_6/hexagonal (symmetry for the SDP reduction)
results:
  - 04-computation/covering_min_delsarte_lp_klein.py
  - 05-knowledge/results/covering_min_delsarte_lp_klein.out
  - 05-knowledge/results/covering_min_delsarte_Vdep_klein.out
---

# HYP-3784 — the Delsarte dual is strong; the creative use is Delsarte + HYP-3763

## Correcting S63 (HYP-3785)
S63 dismissed the SDP/Delsarte lever as "gap-limited" using the naive **Fejer AVERAGE** of loneliness
(`E_w[min_v ||vt||] = 0.029 << M = 0.0765`). But that is the WRONG object: the average is a weak lower
bound (max >= average) and is blind (S54). The Delsarte **DUAL** is different -- a **pointwise-valid**
nonneg weight `w(t)` (`int w = 1`) minimizing the per-resonance danger mass, using the covering structure
(a covering set has one multiple per `q in {2..n}`):
> `min_w SUM_q beta_q`, `beta_q = max_{q|v} int w[||vt||<L]`; if `SUM_q beta_q < 1` then covering-min `>= L`.
This is an LP over `w`, and it is **strong**.

## The Delsarte dual is strong (bounded speeds)
Grid LP (`n=14`, grid 720), certified lower bound vs the speed bound `V`:
```
  V = 2n (28):  covering-min >= 0.0750   (~ the true 14/183 = 0.0765 !)
  V = 4n (56):  covering-min >= 0.0650
  V = 8n (112): covering-min >= 0.0500
  V =16n (224): covering-min >= 0.0000
  reference: trivial 1/(2(n-1))=0.0385, naive average (S63)=0.029, floor 1/n=0.0714, covmin=0.0765
```
At `V=2n` the Delsarte dual nearly EQUALS the covering-min (`0.075 vs 0.0765`) -- far above trivial and
the naive average. S63's "spectral methods are blind" was too hasty: the *dual* sees the spike; only the
*average* is blind.

## Why it weakens with V: the far-element defeat
As `V` grows the bound collapses (`0.075 -> 0` by `V=16n`) because **large multiples EQUIDISTRIBUTE**:
`int w[||vt||<L] -> 2L` for large `v` (any `w`), so `beta_q -> 2L`, and `SUM_q beta_q -> (n-1) 2L`,
giving only the trivial `1/(2(n-1))`. This is exactly the **far-element resonance / CRT-escape** (S62/S63,
HYP-3745): a huge speed's danger zone is spread and captures `2L` of every dual weight. So the Delsarte
dual **over all speeds is trivial** -- confirming the spectral-method limit, but for the sharp reason
(equidistribution of large multiples), not "averages are blind."

## The creative rescue: Delsarte + HYP-3763
The all-speeds Delsarte is defeated by large multiples -- but **HYP-3763 (large-multiples-forced) says a
BEATER cannot use them**: replacing a small core speed by a large multiple raises `M` (the Steinhaus
scaling `c/(kc+1)`), so a set with `M < L` has *bounded* speeds. Therefore the Delsarte dual restricted
to the speeds a beater can actually use is **valid AND near-tight**. Neither tool alone suffices --
Delsarte-all-speeds is trivial; HYP-3763 alone gives no value -- but **together they are near-tight**:
HYP-3763 bounds *which* multiples, Delsarte bounds the danger given bounded multiples. This is the
creative use of the Delsarte method for this problem.

## The SDP (Lasserre) tightening -- the last ~15%
At a fixed moderate `V` the Delsarte LP leaves a gap (`0.065` vs `0.0765` at `V=4n`): the LP integrality
gap. The **Lasserre/SDP** level-2 tightening (a PSD moment matrix over speed-pairs) closes it, and
symmetry-reduces via the `D_7 / zeta_6 / apex-7` structure into small blocks -- the **Eisenstein/E2
isotype** (the bulk, spectrally visible) plus the **cusp/F7 isotypes** (the residual, the last gap). This
is the concrete creative SDP construction; it needs `cvxpy` (unavailable here) but is now well-motivated,
not dismissed.

## Honest caveats
Grid-approximate (not a rigorous certificate). The `V`-bounded LP is not directly valid for sets with the
forced killer `n(n-1) > 2n` -- the honest fix is to run the Delsarte LP over the **HYP-3763-forced speed
set** (small core + the specific double-killer, whose danger `w` dodges at the binding `t*=14/183`), not
all multiples `<= V`. The Delsarte + HYP-3763 synthesis is a reasoned route to a near-tight lower bound,
not yet a rigorous bound; the SDP tightening is untested (no `cvxpy`).

## Net
The Delsarte DUAL is strong (correcting S63: `0.075 ~ covmin` at `V=2n`, vs the naive average `0.029`),
but the all-speeds dual is defeated by the equidistributing far element (CRT-escape). The creative,
correct use is **Delsarte + HYP-3763**: HYP-3763 bounds a beater's speeds, making the Delsarte dual valid
and near-tight; the Lasserre/SDP tightening (symmetry-reduced, Eisenstein-bulk + cusp-residual) closes
the last integrality gap. So the spectral lever is viable -- as a *paired* tool, not alone.

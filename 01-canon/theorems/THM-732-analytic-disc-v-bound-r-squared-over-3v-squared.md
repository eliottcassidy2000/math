---
id: THM-732
title: The analytic disc_v bound — disc_v ≤ r²/(3v²) (r = #arcs of the leave-one-out good set G'_{~v}), PROVED rigorously from the trivial endpoint bound |U(ℓ)|≤2r. Fed into THM-731 it gives the fully explicit RIGOROUS certificate L ≥ (1/7)(6|G'_{~v}| − √2·r/v), so L>0 whenever r < 3√2·v|G'_{~v}| — a COMBINATORIAL (arc-count) condition, no longer harmonic analysis. This certificate CERTIFIES L>0 for the covering-min extremals (deep well ratio 0.46, min-L residue ratio 0.92), the binding families every other method failed on. The crude constant is too lossy only for non-extremal families with a small far element (e.g. {1,3..14}, large true L=0.030) — which need the sharper shared endpoint cancellation (the density Q_s object) or a compactness argument
status: PROVED (the bound) + PARTIAL (universal closure). disc_v ≤ r²/(3v²) is RIGOROUS (|U(ℓ)|≤2r endpoints, Σ_{m≠0}1/m²=π²/3). Hence L ≥ (1/7)(6|G'_{~v}|−√2 r/v) is a RIGOROUS explicit lower bound for every core v. VERIFIED to certify L>0 on the covering-min extremals (deep well, near-AP residue) and a census of large-far-element covering families. NOT universal: for covering families with a small far element and moderate |G'_{~v}| (e.g. {1,3,4,…,14}, far element 14, true L=0.030) the crude constant exceeds the threshold at every peel, though the TRUE disc_v (THM-731) still certifies (+0.018) — so those families need the sharper endpoint-sum cancellation (shared with the density Q_s, THM-729) or a separate large-L/compactness argument. The analytic bound requested is PROVED; it discharges the extremals and reduces the remainder to a combinatorial arc-count inequality plus the shared cancellation for small-far-element sets.
source: klein-2026-07-13-S288
depends_on:
  - THM-731   # the certificate this feeds: L=(6/7)|G'_{~v}|−ε_v, |ε_v|²≤(6/49)disc_v
  - THM-729   # the endpoint-sum cancellation the small-far-element families still need (shared)
related:
  - HYP-6495  # klein-S288 (this bound)
  - HYP-6485  # klein-S287 (THM-731, the certificate)
  - THM-724   # deep well = M-extremal (certified here)
  - THM-726   # deep well unique global covering-min
---

# THM-732 — the analytic disc_v bound: disc_v ≤ r²/(3v²)

The step that upgrades THM-731's certificate from VERIFIED to a **rigorous explicit inequality**, and
converts the remaining difficulty from harmonic analysis to arc-counting.

## The bound (RIGOROUS)

Write the leave-one-out good set `G'_{~v}=∪_{i=1}^r[a_i,b_i]` as `r` arcs, with `2r` endpoints carrying
signs `ε_p=±1`. Its Fourier coefficients are `ĉ_ℓ = U(ℓ)/(2πiℓ)` for `ℓ≠0`, where
`U(ℓ)=Σ_p ε_p e(−ℓp)` is the endpoint sum. The trivial triangle bound `|U(ℓ)| ≤ Σ_p|ε_p| = 2r` gives
$$disc_v=\sum_{m\ne0}|\hat c_{mv}|^2=\sum_{m\ne0}\frac{|U(mv)|^2}{(2\pi mv)^2}
\le \sum_{m\ne0}\frac{(2r)^2}{(2\pi mv)^2}=\frac{4r^2}{4\pi^2v^2}\cdot\frac{\pi^2}{3}
=\boxed{\dfrac{r^2}{3v^2}}.$$

Nothing beyond `|U|≤2r` and `Σ_{m≠0}m^{-2}=π²/3` is used — fully rigorous, universal.

## The explicit certificate (RIGOROUS)

Feeding `disc_v ≤ r²/(3v²)` into THM-731 (`L=(6/7)|G'_{~v}|−ε_v`, `|ε_v|≤√((6/49)disc_v)`):
$$L \ge \tfrac67|G'_{~v}| - \sqrt{\tfrac{6}{49}\cdot\tfrac{r^2}{3v^2}}
= \tfrac67|G'_{~v}| - \tfrac{\sqrt2}{7}\,\tfrac{r}{v}
= \frac1{7}\Big(6\,|G'_{~v}| - \sqrt2\,\tfrac{r}{v}\Big).$$
So **`L>0` whenever `r < 3\sqrt2\, v\,|G'_{~v}|`** — a condition on the arc count `r`, the peel element
`v`, and the measure `|G'_{~v}|`. The harmonic analysis is discharged; what remains is combinatorial.

## What it certifies — and where the crude constant is too lossy

Peeling the far element (`disc` shrinks like `1/v²`), the explicit bound gives (NG=2²¹, verified):

| family | `v` | `|G'_{~v}|` | `r` | `L ≥ (1/7)(6|G'|−√2 r/v)` | ratio `√2 r/(6v|G'|)` |
|---|---|---|---|---|---|
| deep well `{1..12,182}` (M-min) | 182 | 0.03410 | 12 | **+0.0159** | 0.46 |
| near-AP residue `{1..11,13,84}` (L-min) | 84 | 0.01216 | 4 | **+0.0008** | 0.92 |
| `{1..10,12,13,154}` | 154 | 0.06192 | 10 | +0.0400 | 0.25 |
| `{2..14}` | 14 | 0.07143 | 2 | +0.0324 | 0.47 |

**The covering-min extremals certify** — the deep well (the proven global `M`-min, THM-724/726) and the
near-AP residue (the min-`L` `|core|=1` body, kps cont.70), the families every structural and elementary
method failed on. The arc count `r` is **small** (≤12, not the `~78` worst case) because the good sets are
heavily-overlapped small-measure unions — that smallness is what makes the crude bound suffice here.

**Where it is not enough.** For covering families with a *small* far element and moderate `|G'_{~v}|`
— e.g. `{1,3,4,…,14}` (far element 14, true `L=0.030`) — the crude `r²/(3v²)` exceeds the threshold at
*every* peel, so this rigorous certificate does not fire, even though the family is *easy* (large `L`) and
the **true** `disc_v` (THM-731) certifies it comfortably (`+0.018` at `v=8`). The gap is the endpoint-sum
cancellation `|U(mv)|≪2r` that the crude bound discards — the **same** cancellation the density route needs
for `Q_s` (THM-729). So these families reduce to the shared endpoint-cancellation estimate, or to a
separate large-`L`/compactness argument for small-far-element sets.

## Status of the covering case after THM-732

- The analytic `disc_v` bound is **PROVED** (`disc_v ≤ r²/(3v²)`).
- The covering-min **extremals are certified** `L>0` rigorously and explicitly.
- Universal closure reduces to: **(i)** the combinatorial arc-count inequality `r < 3√2 v|G'_{~v}|` for
  large-far-element covering sets (arc-counting, tractable), and **(ii)** the shared endpoint cancellation
  `|U(mv)|≪2r` (= density `Q_s`, THM-729) for the residual small-far-element sets.

*Files: `04-computation/lrc14_disc_v_bound_klein_S288.py`, `lrc14_disc_v_census_klein_S288.py`,
`lrc14_disc_v_failfamily_klein_S288.out` (+.out). HYP-6495. Upgrades THM-731; the residual (ii) is
THM-729 shared with density.*

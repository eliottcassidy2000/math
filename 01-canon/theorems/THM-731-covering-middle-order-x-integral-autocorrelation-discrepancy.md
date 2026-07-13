---
id: THM-731
title: The covering middle-order x-integral — the per-core correlation ε_v = Cov(1_{D_v},1_{G'_{~v}}) satisfies |ε_v|² ≤ (6/49)·[v-grid discrepancy of the good-set autocorrelation A_{~v}], a POSITIVE-DEFINITE spatial object (the covering analogue of THM-729). Combined with the peeling identity L=(6/7)|G'_{~v}|−ε_v it gives a rigorous lower-bound CERTIFICATE L ≥ (6/7)|G'_{~v}|−√((6/49)·disc_v) that avoids the S266 Fourier divergence (keeps 1_{G'} intact). VERIFIED to CERTIFY L>0 on the two covering-min extremals (deep well, near-AP residue) + stress cases, tight to 7–21%, with the certificate ordering matching the L ordering monotonically — the METRIC surrogate that passes mac-mini-S83's acid test (structural deficits fail it)
status: PARTIALLY PROVED. The inequality (*) |ε_v|² ≤ (6/49)·disc_v and the peeling identity L=(6/7)|G'_{~v}|−ε_v are RIGOROUS (Cauchy–Schwarz + Wiener–Khinchin + the Poisson v-grid sampling formula; ε_v-spectrum supported on vℤ because 1_{D_v}(t)=h(vt)). Hence L ≥ L_cert(v):=(6/7)|G'_{~v}|−√((6/49)disc_v) is a RIGOROUS lower bound for every core v. What is VERIFIED-not-proved: that L_cert>0 (i.e. the certificate is positive and near-tight) — confirmed numerically on 4 covering families (NG=2²¹), best peel = the far element, gap 7–21% (near-exact for small-far {2..14}). To make L>0 a THEOREM one needs an ANALYTIC upper bound on disc_v (the v-grid discrepancy of the good-set autocorrelation) — now a POSITIVE, GEOMETRIC quantity (set-overlap), governed by the good-set spectral decay (edge structure), NOT a divergent signed sum.
source: klein-2026-07-13-S287
depends_on:
  - THM-729   # the density analogue this mirrors (autocorrelation Riemann-discrepancy, w-grid)
  - THM-724   # covering-min rigidity: deep well is the M-extremal
  - THM-726   # deep well unique global covering-min
related:
  - HYP-6485  # klein-S287 (this construction)
  - HYP-6455  # klein-S285 (one lattice, three cosets — this is the covering-coset metric form)
  - HYP-6475  # klein-S286 (relation-order stratification — this realizes the "metric x-integral" recommendation)
  - HYP-6465  # opus-S269 (ε_v multi-linear-dominated; this routes AROUND it by keeping A_{~v} intact)
  - HYP-6480  # mac-mini-S83 (structural surrogates all severed; this is the metric one that tracks L)
---

# THM-731 — the covering middle-order x-integral is a positive-definite autocorrelation discrepancy

The faithful covering analogue of THM-729, built on the owner's S286 recommendation: *do not Fourier-expand
`1_{G'}` (the S266 divergence); build the positive-definite `x`-integral of the middle-order sum directly.*

## Setup

Good set `G' = {t∈[0,1) : ‖wt‖ ≥ 1/14 ∀w∈S}` (the `1/14`-lonely times), `L=|G'|`. Bad arc
`D_w={‖wt‖<1/14}`, `|D_w|=1/7`. Leave-one-out good set `G'_{~v}=∩_{w≠v}(complement of D_w)`. opus's
per-core correlation is `ε_v = Cov(1_{D_v},1_{G'_{~v}}) = |D_v∩G'_{~v}| − (1/7)|G'_{~v}|`.

**Peeling identity (RIGOROUS).** `1_{G'}=1_{G'_{~v}}(1−1_{D_v})`, so
$$L=(6/7)\,|G'_{~v}| - \varepsilon_v. \tag{peel}$$
`ε_v>0` shrinks the good set; `ε_v<0` grows it. A lower bound on `L` needs an upper bound on `ε_v`.

**This targets `L=p₀` directly (consistent with mac-mini-S84).** mac-mini-S84 showed the inclusion-exclusion
*moment* expansion `L=Σ_k(−1)^k E_k = p₀` is a trivial binomial identity — the "±20 middle *moment*-orders
`|T|=6,7`" are a binomial artifact of that expansion, not analytic structure. THM-731 does **not** go through
the moment expansion: the peeling identity reaches `L=p₀` (the safe measure) directly, and the genuine
*resonance*-order content (opus-S269's noncore-pair resonances, the Fourier-side `|a|₁`) is absorbed **intact**
into the autocorrelation `A_{~v}`. So "middle-order x-integral" here means the metric object mac-mini-S84
identified as the only live tool — not the deflated moment order.

## The inequality (*) — RIGOROUS

`1_{D_v}(t)=h(vt)` with `h=1_{‖·‖<1/14}` the measure-`1/7` arc; its Fourier coefficients live on the
`v`-grid: `\hat{1_{D_v}}(ℓ)=\hat h(ℓ/v)` if `v∣ℓ`, else `0`. Writing `ĉ_ℓ=\hat{1_{G'_{~v}}}(ℓ)`,
$$\varepsilon_v=\sum_{m\ne0}\hat h(m)\,\overline{\hat c_{mv}}.$$
Cauchy–Schwarz in `m`, then two evaluations:
$$|\varepsilon_v|^2\le\Big(\sum_{m\ne0}|\hat h(m)|^2\Big)\Big(\sum_{m\ne0}|\hat c_{mv}|^2\Big)
=\frac{6}{49}\cdot\Big[\frac1v\sum_{j=0}^{v-1}A_{~v}(j/v)-|G'_{~v}|^2\Big]. \tag{*}$$
- `Σ_{m≠0}|\hat h(m)|² = ∫h² − |\hat h(0)|² = 1/7 − 1/49 = 6/49` (`h` an indicator).
- `Σ_{m≠0}|ĉ_{mv}|² = Σ_m|ĉ_{mv}|² − |ĉ_0|²`. By the Poisson sampling identity `Σ_m\hat A(mv)=(1/v)Σ_{j}A(j/v)`
  applied to the **autocorrelation** `A_{~v}(τ)=|G'_{~v}∩(G'_{~v}−τ)|` (whose transform is `|ĉ_ℓ|²≥0`), and
  `ĉ_0=|G'_{~v}|`, this equals `(1/v)Σ_{j=0}^{v-1}A_{~v}(j/v)−|G'_{~v}|²`.

The bracket `disc_v := (1/v)Σ_{j=0}^{v-1}A_{~v}(j/v)−|G'_{~v}|² = Σ_{m≠0}|ĉ_{mv}|² ≥ 0` is the
**`v`-grid discrepancy of the good-set autocorrelation**: how far the `v`-grid average of the set-overlap
function `A_{~v}` sits above its true mean `|G'_{~v}|²`. It is a positive-definite **spatial** (x-)integral,
manifestly `≥0`, and involves **no Fourier expansion of the product** `∏(1−1_{D_w})` — the divergence of
S266 is entirely avoided; the multi-linear middle-order content (opus-S269) is absorbed, intact, into
`A_{~v}`.

## The certificate

Combining `(*)` with `(peel)` (`ε_v ≤ |ε_v|`):
$$\boxed{\,L \ge L_{\mathrm{cert}}(v) := \tfrac67|G'_{~v}| - \sqrt{\tfrac{6}{49}\,disc_v}\,}\qquad\text{for every core }v.$$
This is a **rigorous lower bound**; `L_cert(v)>0` certifies a `1/14`-lonely time exists. One is free to
peel the **best** core; empirically that is the **far element** (large `v` ⟹ fine `v`-grid ⟹ small
`disc_v`), which **inverts** opus-S269's difficulty: the large core `v≥17` that is hardest for the
cluster/Fourier route is the **easiest** for the metric route.

## What is verified (NG=2²¹, `lrc14_covering_autocorr_leaveoneout_klein_S287.py`)

`(*)` holds for every core on all families. Best-peel certificate vs true `L`:

| family | true `L` | best `L_cert` | gap | ρ (tightness) |
|---|---|---|---|---|
| deep well `{1..12,182}` (M-min) | 0.02390 | **0.02212** (peel 182) | 7% | 0.56 |
| near-AP residue `{1..11,13,84}` (L-min) | 0.00536 | **0.00424** (peel 84) | 21% | 0.67 |
| small-far `{2..14}` | 0.06122 | **0.06118** (peel 14) | ~0% | — |
| variant `{1,3..13,182}` | 0.02630 | **0.02494** (peel 182) | 5% | — |

**All four certify `L>0`.** The certificate ordering `0.0042<0.0221<0.0249<0.0612` is a **perfect monotone
match** to the `L` ordering `0.0054<0.0239<0.0263<0.0612` — the certificate correctly flags the near-AP
residue as the **binding/most-stuck** family. This is exactly the ordering every *structural* deficit gets
**wrong** (mac-mini-S83: `{1..11,13,84}` dominates the deep well on Schur/Q₂ deficit yet has smaller `L`).
So the autocorrelation discrepancy is a **faithful metric surrogate** — the object mac-mini-S83 proved no
structural invariant can be.

## The remaining analytic step

**DONE (THM-732, klein-S288):** `disc_v ≤ r²/(3v²)` (`r`=#arcs of `G'_{~v}`), proved from the trivial
endpoint bound `|U(ℓ)|≤2r`. This gives the fully explicit rigorous certificate `L ≥ (1/7)(6|G'_{~v}|−√2 r/v)`,
which **certifies `L>0` on the covering-min extremals** (deep well, min-`L` residue) and reduces universal
closure to the **combinatorial** `r < 3√2 v|G'_{~v}|`. The crude constant is too lossy only for
small-far-element easy families (`{1,3..14}`), which reduce to the shared endpoint cancellation
`|U(mv)|≪2r` (= the density `Q_s`, THM-729). So the analytic step here is discharged for the hard cases;
see THM-732.

*Files: `04-computation/lrc14_covering_autocorr_leaveoneout_klein_S287.py` (+.out),
`lrc14_covering_autocorr_xintegral_klein_S287.py` (+.out). HYP-6485. Mirrors THM-729; realizes the
klein-S286 (HYP-6475) "build the metric x-integral" recommendation for the covering side.*

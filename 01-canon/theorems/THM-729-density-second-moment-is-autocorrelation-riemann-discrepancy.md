---
id: THM-729
title: The density √-cancellation — the endpoint-sum 2nd moment Q_s=Σ_ℓ|U_s(ℓw)|²/ℓ² equals (2πw)²·[Riemann-discrepancy of the autocorrelation A of 1_{R_s} at the w-grid], a 1-D discrepancy (NOT multi-linear); Q_s=O(diam) (verified) ⟹ |S|=O(√diam) by Cauchy–Schwarz ⟹ Error→0, closing the density row on a finite box. The density route needs only this 1-D autocorrelation discrepancy, where the covering route needs the multi-linear Gowers cancellation — the tractability asymmetry made concrete
status: PARTIALLY PROVED. The exact identity (Q_s = autocorrelation Riemann-discrepancy) and the arc-wise DIAGONAL bound (diag = Σ_i 2π²{w·w_i}(1−{w·w_i}) = O(M), matching the full Q_s) are RIGOROUS. Q_s=O(diam) is EMPIRICAL (7 clusters, Q_s/diam∈[1,1.7]); |S|≤C√diam and Error=C√diam/d→0 follow. The full rigorous upper bound = a 1-D discrepancy of the piecewise-linear A at the w-grid (the crude min-separation large sieve fails: endpoint δ→0 under clustering; the arc-differencing self-cancels tiny arcs, which the autocorrelation identity captures but a clean bound must exhibit).
source: klein-2026-07-13-S280
depends_on:
  - THM-727   # the Fourier reduction |S| ≤ (1/2π²)Σ_s Σ_ℓ |U_s(ℓw)||sin(πℓ/7)|/ℓ²
  - THM-728   # U_s = the endpoint sum; the derivative structure
related:
  - HYP-6415  # klein-S280
  - HYP-6410  # klein-S279 — the multi-linear crux (covering needs sharp; density needs only this √)
---

# THM-729 — the density 2nd moment is the autocorrelation's Riemann-discrepancy

## The √-reduction

From THM-727, `|S|=|w·Error| ≤ (1/2π²)Σ_s Σ_{ℓ≠0}|U_s(ℓw)||\sin(πℓ/7)|/ℓ²`. Cauchy–Schwarz in `ℓ`:

$$\sum_{\ell}\frac{|U_s(\ell w)||\sin(\pi\ell/7)|}{\ell^2}\le \sqrt{Q_s}\cdot\sqrt{\textstyle\sum_\ell \sin^2(\pi\ell/7)/\ell^2},
\qquad Q_s:=\sum_{\ell\ne0}\frac{|U_s(\ell w)|^2}{\ell^2}.$$

So `|S| = O(√(max_s Q_s))`. **If `Q_s=O(diam)` then `|S|=O(√diam)`**, hence on the peel `w=d≥diam`,
`Error=|S|/w=O(√diam)/d=O(1/√d)→0`. For `d>D_0` (finite) `Error<cap₉−0.397`; the band `26≤d≤D_0` is
the structured check (S275, max `Φ=0.347`), and resonant peels (`w=lcm≫Σe'`) are S275-harmless. So
`Q_s=O(diam)` closes the density row.

## The exact identity (RIGOROUS)

`|U_s(N)|=2\pi|N|\,|\hat f(N)|`, `f=1_{R_s}`. With the Dirac-comb identity
`Σ_{ℓ∈ℤ}|\hat f(ℓw)|²=(1/w)Σ_{m=0}^{w-1}A(m/w)`, `A(t)=|R_s∩(R_s−t)|` the autocorrelation (`A(0)=|R_s|`,
`∫_0^1A=|R_s|²`):

$$\boxed{\ Q_s=(2\pi w)^2\Big[\tfrac1w\sum_{m=0}^{w-1}A(m/w)-\int_0^1 A\Big]\ }$$

— `Q_s` is exactly `(2πw)²` times the **left-Riemann-discrepancy of the autocorrelation `A` at the
`w`-grid**. `A` is piecewise-linear with `O(M)` breakpoints (`M=#\partial R_s=O(diam)`). This is a
**1-D discrepancy**, not a multi-linear cancellation.

## What is bounded

- **Empirical:** `Q_s=O(diam)` (`Q_s/diam∈[1.0,1.7]` over 7 clusters, `diam` to 199) ⟹ the √-cancellation
  and the row-closure above.
- **Diagonal (rigorous):** arc-wise, `Q_s=\sum_{i}2\pi^2\{w w_i\}(1-\{w w_i\})+(\text{off-diag})`,
  `w_i=`arc widths; the diagonal is `≤(π^2/2)\cdot(M/2)=O(M)` and **empirically ≈ the full `Q_s`** (the
  off-diagonal cancels). So the `O(M)` scaling has a rigorous, dominant backbone.
- **The crude large sieve fails:** `Q_s≤M(2+\tfrac43δ^{-1})`, `δ=\min\|w(p-p')\|`, is useless because
  `δ→0` under endpoint clustering (measured `δ=0.012→0.0000`). The tiny-arc contributions self-cancel
  (`e(-Na)-e(-Nb)\to0`), which the autocorrelation identity captures but the min-separation bound does not.

## The asymmetry made concrete (vs covering)

The density route needs only `Q_s=O(diam)` — a **1-D discrepancy of `A` at the `w`-grid** (`A` piecewise
linear, `O(M)` breakpoints). The covering route needs the **multi-linear (Gowers) cancellation**
(opus-S262/S263, mac-mini-S76: averaging provably insufficient, 3rd-order Schur required). So the density
√-route is genuinely lower-order and more tractable, exactly the S279 asymmetry — now pinned to a concrete
1-D discrepancy object.

## Files
`04-computation/lrc14_second_moment_klein_S280.py`, `lrc14_delta_sep_klein_S280.py` (+ outs).

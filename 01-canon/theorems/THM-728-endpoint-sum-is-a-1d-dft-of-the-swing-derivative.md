---
id: THM-728
title: The coupled endpoint sum is a 1-D DFT of a swing-DERIVATIVE sequence — U_s^{e'}(N)=e'·χ̂(N mod e') EXACTLY (no k-dim Weyl), and ch_j=cond_s(leave_j)−cond_s(enter_j) is a discrete derivative giving U_s^{e'}(N)=e'Σ_{n≡κ mod e'}ĝ(n)e(ns/7e')(e(n/7e')−1) with the sin(πn/7e') derivative-gain; the far frequency and e' DROP OUT and |U_s^{e'}| is bounded (O_k(1), verified e'→400) ⟹ |S|=O(k). Remaining rigor = the conditionally-convergent 1-D sum needs Beurling–Selberg mollification (SHARED with the covering side)
status: PARTIALLY PROVED. The 1-D DFT reduction and the swing-derivative identity are RIGOROUS. The boundedness |U_s^{e'}(N)|≤C_k (independent of e' and of the far frequency N=ℓw) is EMPIRICAL (verified to e'=400: ≤4 for 6 others, ≤13 for 7 others; net imbalance χ̂(0)≤4). The rigorous constant needs the derivative-gain sum Σ_{n≡κ}ĝ(n)(e(n/7e')−1)=O(1), a conditionally-convergent 1-D Fourier sum requiring Beurling–Selberg mollification of the cond_s indicator.
source: klein-2026-07-13-S278
depends_on:
  - THM-727   # the exact Fourier reduction S=Σ_sΣ_ℓ(−1/2πiℓ)U_s(ℓw)ĝ_s(ℓ)
related:
  - HYP-6400  # klein-S278 — this session
  - opus-S261 # covering-side Beurling–Selberg mollification — the SAME analytic tool
---

# THM-728 — the endpoint sum is a 1-D DFT of the swing-derivative

## The reduction is 1-dimensional (RIGOROUS)

Grouping `U_s(N)=Σ_{e'}U_s^{e'}(N)` by owning offset, the endpoints of `R_s` owned by `e'` sit at
`(j+σ/7)/e'`, and — the decisive simplification — the **coupling indicator `χ(j)` is periodic mod `e'`
in `j`** (it depends on the other offsets' phases `frac(e''(j+σ/7)/e')=frac((e''/e')j+…)`, all periodic
mod `e'`). Hence the inner sum is a **1-D DFT** over `ℤ/e'`, with **no `k`-dimensional Weyl sum**:

$$U_s^{e'}(N)=\sum_\sigma e(-N\sigma/7e')\cdot e'\,\hat\chi_{s,\sigma}(N\bmod e'),\qquad
\hat\chi(\kappa)=\tfrac1{e'}\sum_{j=0}^{e'-1}\chi(j)e(-\kappa j/e').$$

The far frequency enters only through `κ = N mod e' = ℓw mod e'`.

## The coupling is a discrete DERIVATIVE (RIGOROUS)

The `R_s`-endpoint indicator is `ch_j = χ_{s}(j) = 1_{\text{cond}_s}(\text{leave}_j)−1_{\text{cond}_s}(\text{enter}_j)`,
where `cond_s(x)=[`the other offsets occupy exactly `\{0..6\}\setminus\{s\}]`, and `leave_j-enter_j=δ=1/(7e')`.
So `ch` is a **discrete derivative** of the sampled `cond_s` indicator `g=1_{cond_s}`. With `\hat g(n)`
its Fourier coefficients,

$$U_s^{e'}(N)=e'\!\!\sum_{n\equiv\kappa\,(e')}\!\!\hat g(n)\,e(ns/7e')\big(e(n/7e')-1\big),\qquad
|e(n/7e')-1|=2|\sin(\pi n/7e')|.$$

The `\sin(\pi n/7e')` factor is the **derivative gain**: it is `≈ πn/(7e')` for low `n`, cancelling one
power against `\hat g(n)\sim \rho_{cond}/|n|` and killing the diagonal `n=0` (the net-imbalance term
`\hat\chi(0)` is thereby `O(1)`, not `O(e')` — verified `≤4`).

## Boundedness ⟹ O(k) (EMPIRICAL, decisive)

`|U_s^{e'}(N)| = e'|\hat\chi(N\bmod e')|` is **bounded independent of `e'` and of the far frequency**:
verified to `e'=400`, `≤4` (6 other offsets) and `≤13` (7 others); the swing-count `T` is constant in
`e'`. Therefore, via THM-727,

$$|S|\le\frac1{2\pi^2}\Big(\sum_{\ell\ne0}\frac{|\sin(\pi\ell/7)|}{\ell^2}\Big)\sum_s\sum_{e'}|U_s^{e'}(\ell w)|
= O_k(1)=O(k),$$

for bounded `k` (the row needs `k≤8`) — an absolute constant. The estimate **closes the density-row
tail**, confirming the S276 total bound `|S|≤0.61R` with a mechanism, not a fit.

## Remaining rigor (one clean 1-D sum) — and the unification

The bound `|U_s^{e'}(N)|≤C_k` reduces to `Σ_{n≡κ (e')}\hat g(n)(e(n/7e')-1)=O(1)`. The `cond_s`
indicator `g` has non-summable Fourier coefficients (`\hat g(n)\sim1/n`), so this conditionally-convergent
sum needs **Beurling–Selberg mollification** of `g` (smooth minorant/majorant, rapid Fourier decay, at a
boundary cost). This is **exactly** the tool the covering side reduced to (opus-S261: mollified
discrepancy of the coprime core). Both LRC(14) routes now terminate on the same Beurling–Selberg
mollification estimate.

## Files
`04-computation/lrc14_chi_dft_klein_S278.py`, `lrc14_chi_largeeprime_klein_S278.py` (+ outs).

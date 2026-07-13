---
id: THM-727
title: The two-scale error times w reduces EXACTLY to the R_s-endpoint exponential sums (w cancels) — |S| ≤ (1/2π²)Σ_s Σ_ℓ |U_s(ℓw)|·|sin(πℓ/7)|/ℓ², with U_s(N)=Σ_{endpoints}ε_p e(−Np); the UNCOUPLED diagonal is exactly 7e'·[e'∣ℓw] per offset (⟹ O(gcd(w,e')²/e'), Denjoy–Koksma), and the full O(k) bound reduces to a k-dimensional Weyl estimate on the COUPLED endpoint sums Σ_{j∈J}ε_j e(−Nj/e'). CORRECTS S276: per-offset is NOT ≤1 (reaches ~1.6); the robust bound is the TOTAL |S| ≤ 0.61·R
status: PARTIALLY PROVED. The Fourier reduction (exact identity + triangle inequality) and the uncoupled diagonal bound are RIGOROUS. The coupled endpoint sum bound (the remaining step) is a standard-but-uncarried k-dimensional Weyl / Erdős–Turán estimate; the total |S| ≤ 0.61·R and per-offset ≤ ~1.6 are EMPIRICAL (30-case sweep, prime grid).
source: klein-2026-07-13-S277
depends_on:
  - THM-725   # the per-interval reduction (Error = Σ_I ∫_I g_{s_I}(wx))
  - THM-700   # the cover-measure two-scale error (this is its exact Fourier form)
related:
  - HYP-6350  # klein-S276 — the O(k) / resonance-sum bound (this reduces it + corrects the ≤1)
  - HYP-6380  # klein-S277 — this session
---

# THM-727 — the two-scale error's exact reduction to endpoint exponential sums

## The reduction (RIGOROUS)

`Error_cover = Σ_s ∫_0^1 1_{R_s}(x) g_s(wx) dx`, `R_s = {E'` misses exactly sector `s}` (THM-700/725),
`g_s(y)=1{y∈[s/7,(s+1)/7)}−1/7`. Fourier-expand both factors; orthogonality collapses to the diagonal:

$$S := w\cdot\text{Error}_{\text{cover}} = \sum_s \sum_{\ell\ne 0} \frac{-1}{2\pi i\,\ell}\,U_s(\ell w)\,\hat g_s(\ell),
\qquad U_s(N)=\sum_{p\ \text{endpt of }R_s}\varepsilon_p\,e(-Np),$$

the factor `w` **cancelling** against `\hat{1_{R_s}}(-\ell w)=U_s(-\ell w)/(2\pi i(-\ell w))`. With
`|\hat g_s(\ell)|=|\sin(\pi\ell/7)|/(\pi|\ell|)`:

$$|S|\ \le\ \frac{1}{2\pi^2}\sum_{s=0}^{6}\sum_{\ell\ne 0}\frac{|U_s(\ell w)|\,|\sin(\pi\ell/7)|}{\ell^2}.$$

Everything reduces to the **endpoint exponential sums** `U_s(\ell w)` — a clean, exact reformulation.
(The trivial `|U_s(\ell w)| ≤ 2\rho_s` recovers THM-700's `O(Σe')`; the content is doing better.)

## The uncoupled diagonal (RIGOROUS)

Group `U_s(N)` by the offset `e'` owning each endpoint (`p=(j+\sigma/7)/e'`, evenly spaced at
`i/(7e')`, boundary `i\bmod 7`). If *every* crossing were an endpoint of one type (no coupling),
`Σ_{j=0}^{e'-1} e(-Nj/e') = e'·[e'\mid N]` **exactly** (integer `N=\ell w`, full period). So the
uncoupled `U_s^{e'}(\ell w)=7e'·[e'\mid \ell w]` — nonzero only when `q\mid\ell` (`q=e'/\gcd(w,e')`),
and

$$\sum_{\ell}\frac{|U_s^{e'}(\ell w)|\,|\sin(\pi\ell/7)|}{\ell^2}=\frac{7e'}{q^2}\sum_r\frac{|\sin(\pi rq/7)|}{r^2}=O\!\Big(\frac{\gcd(w,e')^2}{e'}\Big).$$

`O(1/e')` at clean `w` (harmless), `O(e')` at resonance — the Denjoy–Koksma per-offset law.

## The remaining step (the coupling)

The genuine `U_s^{e'}(N)=Σ_\sigma e(-N\sigma/(7e'))Σ_{j\in J_{s,\sigma,e'}}\varepsilon_j e(-Nj/e')`,
where `J` = the `j` giving actual `R_s`-endpoints and `\varepsilon_j\in\{\pm1\}` — both set by the
*other* offsets (each an arithmetic progression in `j` of step `e''/e'`). So `\varepsilon_j=\chi(\alpha j)`
for a fixed-complexity torus indicator `\chi` (`k−1` frequencies), and the inner sum is a
**`k`-dimensional Weyl sum** `Σ_j \chi(\alpha j)e(-Nj/e')`. Bounding it by `O_k(\min(e',1/\|N/e'\|))`
(Erdős–Turán + `\chi`'s Fourier decay, `k` fixed) yields the coupled per-offset bound and hence the
full `O(k)`. This is the one uncarried step — a standard multi-dimensional Weyl estimate.

## Correction to S276 (`per-offset ≤1` was wrong)

Direct decomposition `S=Σ_{e'}S_{e'}` (30-case sweep, prime grid): **`\max_{e'}|S_{e'}| = 1.588`**
(2-block, resonant, `e'=28`) — per-offset is **not** `≤1`; it is bounded by an absolute constant
`≈1.6`. The `≤1` was imprecise shorthand for "each resonance-sum term `≤1`." The **robust** bound is
the *total* `|S| ≤ 0.61·R`, `R=Σ_{e'≥2}\min(1,1/(e'\|w/e'\|))` (max ratio `0.535` over the sweep) —
it uses inter-offset cancellation (`|Σ S_{e'}| < Σ|S_{e'}|`). Either way `|S|=O(k)` and the row closes.

## Files
`04-computation/lrc14_per_offset_decomp_klein_S277.py`, `lrc14_per_offset_C_klein_S277.py` (+ outs).

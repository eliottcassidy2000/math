# The per-offset bound is not ≤1 (it reaches ~1.6) — but the reduction to endpoint exponential sums is exact, and the remainder is one k-dimensional Weyl estimate

*klein-2026-07-13-S277. Owner: prove the coupled per-offset ≤1 bound (the last rigor step from S276).
Two outcomes: the target as literally stated is **false** (per-offset reaches ~1.6, not ≤1), and the
genuine object — the *total* `|S| ≤ 0.61·R` — reduces **exactly and rigorously** to endpoint
exponential sums, leaving one standard-but-uncarried Weyl estimate. Progress plus a correction.*

---

## The correction: `per-offset ≤1` is false

`S = w·Error_cover = Σ_{e'} S_{e'}`, `S_{e'} = Σ_{R-endpoints owned by e'} ε_p G_{s(p)}(wp)`. Directly
decomposing (30-case sweep, prime grid), the per-offset contributions reach

> `max_{e'} |S_{e'}| = 1.588` (2-block cluster, resonant `w=lcm`, offset `e'=28`),

so `|S_{e'}| ≤ 1` is simply **wrong**. My S276 "per-offset ≤1" conflated `|S_{e'}|` with the
resonance-sum term `min(1, 1/(e'‖w/e'‖)) ≤ 1` — a different quantity. The per-offset contribution is
bounded by an absolute constant `≈1.6`, and it does **not** track its own `R`-term (e.g. `S_{30}=1.20`
while that offset's `R`-term is `0.14`). The clean, robust statement is the **total**:

> `|S| ≤ 0.61·R`, `R = Σ_{e'≥2} min(1, 1/(e'‖w/e'‖))` — max ratio `0.535` over the sweep,

which holds because the `S_{e'}` carry signs and partially cancel (`|Σ S_{e'}| < Σ|S_{e'}|`). Either
form gives `|S| = O(k)` and closes the row; only the per-offset "≤1" phrasing was false.

## What *is* rigorous: the exact Fourier reduction (THM-727)

The whole problem collapses cleanly. Fourier-expanding `1_{R_s}` and `g_s(wx)`, orthogonality gives —
with the factor `w` **cancelling** —

$$S = \sum_s\sum_{\ell\ne0}\frac{-1}{2\pi i\ell}\,U_s(\ell w)\,\hat g_s(\ell),\qquad
U_s(N)=\sum_{p\ \text{endpt of }R_s}\varepsilon_p e(-Np),$$
$$|S|\le\frac1{2\pi^2}\sum_s\sum_{\ell\ne0}\frac{|U_s(\ell w)|\,|\sin(\pi\ell/7)|}{\ell^2}.$$

This is an **exact identity + triangle inequality** — no heuristics. It replaces the whole two-scale
constant with the endpoint exponential sums `U_s(\ell w)`. The trivial `|U_s|\le 2\rho_s` recovers
THM-700's `O(Σe')`; the entire game is beating it.

**The uncoupled diagonal is rigorous and clean.** Grouping `U_s` by owning offset, and using that for
*integer* `N=\ell w` the full-period geometric sum is exactly `Σ_{j<e'}e(-Nj/e')=e'[e'\mid N]`, the
uncoupled contribution is `7e'\,[e'\mid\ell w]` — supported only on `\ell\in q\mathbb Z`
(`q=e'/\gcd(w,e')`) — and

$$\sum_\ell\frac{|U_s^{e'}(\ell w)||\sin(\pi\ell/7)|}{\ell^2}=O\!\big(\gcd(w,e')^2/e'\big),$$

`O(1/e')` at clean `w`, `O(e')` at resonance. This is the Denjoy–Koksma per-offset law, now on exact
footing.

## The one remaining step: a k-dimensional Weyl sum

The genuine `U_s^{e'}` is not the full geometric sum: only the `j` that yield *actual* `R_s`-endpoints
count, weighted by `ε_j=±1`, and both are decided by the **other** offsets — each an arithmetic
progression `e''/e'·j` in `j`. So `ε_j = χ(αj)` for a fixed-complexity indicator `χ` on the torus
(`k−1` frequencies), and

$$U_s^{e'}(N)=\sum_\sigma e(-N\sigma/7e')\sum_{j}\chi(\alpha j)\,e(-Nj/e')$$

— a `k`-dimensional **Weyl sum**. For bounded `k` (the row needs `k≤8`), Erdős–Turán plus `χ`'s
Fourier decay give `Σ_j χ(αj)e(-Nj/e') = O_k(\min(e',1/\|N/e'\|))`, which yields the coupled per-offset
bound and hence the full `O(k)`. This is a **standard multi-dimensional Weyl / discrepancy estimate**;
carrying out the constant (and its `k`-dependence) is the uncompleted rigor.

## Why the reduction, not the per-offset, is the right target

The session's lesson: the per-offset decomposition is a *lossy* device — it breaks R-intervals across
owners, so `S_{e'}` neither is `≤1` nor tracks `R`. The clean object is the **endpoint set as a whole**
and its exponential sum `U_s(\ell w)`; the cancellation lives *between* offsets (the signs), which the
per-offset split discards. The Fourier reduction keeps it. So the proof should target `Σ_s Σ_\ell
|U_s(\ell w)|/\ell^2 = O(k)` directly — a one-dimensional-per-`\ell` sum of `k`-dim Weyl sums — rather
than the per-offset `≤1`.

## Honest scope

- **RIGOROUS:** the exact reduction `S=Σ_sΣ_\ell(-1/2\pi i\ell)U_s(\ell w)\hat g_s(\ell)`; the uncoupled
  diagonal `O(\gcd^2/e')`.
- **CORRECTED:** per-offset is `≤1.6` (an absolute constant, not `≤1`); the robust bound is total
  `|S| ≤ 0.61·R`.
- **OPEN (one clean step):** the coupled endpoint Weyl bound `Σ_j χ(αj)e(-Nj/e')=O_k(\min(e',1/\|N/e'\|))`
  — a standard `k`-dim Erdős–Turán estimate. That is the last rigor step for the entire density-row
  tail, now isolated to a single named inequality.

This connects to the covering side, which S260 just reduced to a *mollified discrepancy* of the
coprime core (Beurling–Selberg) — both LRC(14) routes now bottom out on a discrepancy/Weyl estimate.

*Files: `04-computation/lrc14_per_offset_decomp_klein_S277.py`, `lrc14_per_offset_C_klein_S277.py`
(+ outs). THM-727, HYP-6380. Corrects [[the-van-der-corput-Ok-bound-is-confirmed-err-times-w-at-most-c-times-resonance-sum-klein-S276]]. Consumes THM-725/700.*

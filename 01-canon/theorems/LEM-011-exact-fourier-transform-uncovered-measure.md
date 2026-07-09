---
id: LEM-011
title: The EXACT Fourier transform of the uncovered-measure function 𝒲 — 𝒲̂(n) = (−1)^r (6/7)^{k−1−r} ∏_{n_i≠0} b₀(n_i) · (1[σ=0] − c(σ)) in closed form, where b₀(m)=(e(m/7)−1)/(2πim) is the 1/7-arc coefficient, c(σ)=(1−e(−σ/7))/(2πiσ), r=#nonzeros, σ=Σn_i; this is the shared "𝒲̂-decay constant" (mac-mini THM-664 handoff) made EXACT and a-priori — it pins 𝒲̂(0)=(6/7)^k, satisfies Parseval Σ|𝒲̂|²=E[W²], has geometric per-coordinate decay (7/6)/π=0.371, and expresses BOTH the density-floor decorrelation E[W]−(6/7)^k=Σ_{n·e=0}𝒲̂(n) and the THM-664 grid residual E_grid[W]−(6/7)^k=Σ_{Vmax|n·e}𝒲̂(n) as explicit convergent sums
status: PROVED (elementary Fourier computation) + VERIFIED (FFT match to grid-discretization ~1e-4 at k=3,4; Parseval exact; both resonance sums converge to the directly-computed E[W] and E_grid[W]). Replaces the "numerically-certified 𝒲̂ / mixed-variation" placeholder of opus-S157 and THM-664 with an exact closed form. Does NOT by itself give a uniform |Σ_resonance 𝒲̂| < (6/7)^k over all clusters (the signed resonance sum still depends on the cluster's arithmetic); it reduces that to an explicit low-height computation (geometric tail). The LRC(14) large-spread EXISTENCE does not need it (LEM-010 is elementary) — this is the unification + abundance
source: klein-2026-07-08-S194 (mac-mini-S59 THM-664 handoff (a): the shared 𝒲̂-decay)
depends_on:
  - THM-657   # W = uncovered measure = Σ(g_i − 1/7)_+
  - THM-664   # E_grid[W] = (6/7)^k + Σ_{Vmax|n·e}𝒲̂(n); this makes 𝒲̂ exact
related:
  - LEM-009   # density-floor tail: the n·e=0 decorrelation sum, same 𝒲̂
  - LEM-010   # the elementary (Dirichlet) large-spread closure (why 𝒲̂-decay is optional for existence)
  - THM-518   # Weyl decorrelation
external: Fourier analysis on the torus; Parseval.
---

# LEM-011 — The exact Fourier transform of the uncovered-measure function

## Statement

Let `𝒲 : T^{k−1} → [0, 6/7]` be the **uncovered-measure function**
`𝒲(φ₁,…,φ_{k−1}) = ` (Lebesgue measure of the uncovered set of the phase configuration
`{0, φ₁, …, φ_{k−1}}`), i.e. `𝒲 = Σ_gaps (gap − 1/7)_+ = W` with the first phase pinned at `0`
(THM-657). For `n ∈ ℤ^{k−1}`, write `r = #{i : n_i ≠ 0}` and `σ = Σ_i n_i`. Then

> **`𝒲̂(n) = (−1)^r (6/7)^{k−1−r} \big[∏_{i:\,n_i≠0} b₀(n_i)\big]\,\big(𝟙[σ=0] − c(σ)\big)`**

where `b₀(m) = (e(m/7)−1)/(2πim)` is the Fourier coefficient of the `1/7`-arc
(`b₀(0)=1/7`, `|b₀(m)| = |sin(πm/7)|/(π|m|)`), and `c(σ) = ∫₀^{1/7} e(−σt)\,dt = (1−e(−σ/7))/(2πiσ)`
(`c(0)=1/7`, `|c(σ)| = |sin(πσ/7)|/(π|σ|)`), `e(x):=e^{2πix}`.

## Proof

A point `t ∈ [0,1)` is **uncovered** iff the arc `(t−1/7, t)` contains no phase, so
`𝒲(φ) = ∫₀¹ ∏_{i=0}^{k−1} \big(1 − 𝟙_{(t−1/7,t)}(φ_i)\big)\,dt` (with `φ₀ ≡ 0`). Expand the product over
subsets `S ⊆ {0,…,k−1}` and Fourier-transform in `φ₁,…,φ_{k−1}` (`φ₀=0` fixed). For `i ≥ 1`:
`∫_{T} 𝟙_{(t−1/7,t)}(φ_i)\,e(−n_iφ_i)\,dφ_i = b_t(n_i) = e(−n_i t)\,b₀(n_i)` while an absent factor
contributes `𝟙[n_i=0]`; for `i=0`, `𝟙_{(t−1/7,t)}(0) = 𝟙[t∈(0,1/7)]`. Only `S ⊇ supp(n)` survive.
Summing the zero-coordinate subset choices gives `∏_{i:n_i=0,i≥1}(1 − 1/7) = (6/7)^{k−1−r}`
(the `(−1)^{|S|}` and the `b_t(0)=1/7` combine to `(1−1/7)` per zero coord), the `∏_{i∈supp(n)} b_t(n_i)
= e(−σt)∏b₀(n_i)` factors out the common phase, and the two `i=0` cases give the `t`-integral
`∫₀¹ e(−σt)(𝟙[t∈(0,1/7)]^{ε₀})dt = 𝟙[σ=0] − c(σ)` (the `ε₀=0` term yields `𝟙[σ=0]`, the `ε₀=1` term
`−c(σ)`). Collecting the signs gives `(−1)^r`. ∎

**Sanity / verification.** `n=0`: `(6/7)^{k−1}(1 − c(0)) = (6/7)^{k−1}(6/7) = (6/7)^k` — mac-mini's exact
main term (THM-664). FFT of `𝒲` on a grid matches the closed form to grid-discretization (`≤1.3·10⁻⁴`
at `k=3,4`), and Parseval `Σ_n|𝒲̂(n)|² = ∫_{T^{k−1}}𝒲² = E[W²]` holds exactly (`0.400058` at k=3,
`0.296727` at k=4). Files: `lrc14_What_exact_klein_S194.{py,out}`.

## The decay (a-priori)

Immediately from the closed form,

> **`|𝒲̂(n)| ≤ (6/7)^{k−1−r} \big[∏_{i:n_i≠0} \tfrac{1}{π|n_i|}\big]·\min\!\big(\tfrac{6}{7},\,\tfrac{1}{π|σ|}\big)`**, and `𝒲̂(n)=0` whenever `7 ∣ n_i` for some `i`, or `7 ∣ σ` with `σ≠0`.

Each additional nonzero coordinate multiplies the bound by `≤ (7/6)/π = 0.3714 < 1` (the `(6/7)^{−1}`
from the exponent times `1/(π|n_i|) ≤ 1/π`), so the transform is **geometrically concentrated on
low-weight, low-height `n`**. The `σ`-factor adds `1/(π|σ|)` decay in the frequency sum `σ=Σn_i`
except on the balanced set `σ=0`.

## The two resonance sums (the shared object, made explicit)

Let `E={e_0=0<…<e_{k−1}}` be the cluster co-offsets, `n·e := Σ_i n_i e_i`. Averaging
`𝒲(frac(e_1y),…)` gives (THM-664 identity):

- **Density-floor decorrelation (continuous first moment):**
  `E[W] − (6/7)^k = Σ_{n≠0,\ n·e=0} 𝒲̂(n)`. This is the `n·e=0` (exact-relation) resonance sum — the
  LEM-009/opus-S157 tail object.
- **THM-664 grid residual (finite-Vmax):**
  `E_grid[W] − (6/7)^k = Σ_{n≠0,\ Vmax\mid n·e} 𝒲̂(n)` (superset: `n·e=0` **plus** the pure grid
  resonances `n·e = mVmax, m≠0`).

Both are now **explicit convergent sums** of the closed-form `𝒲̂(n)`, dominated by low-height terms
(geometric tail `0.371` per coordinate). Verified: the truncated sums (`‖n‖≤22`) converge to the
directly-computed `E[W]−(6/7)^k` and `E_grid[W]−(6/7)^k` on `k=3,4` clusters (agreement `≤10⁻³`).
Files: `lrc14_resonance_from_What_klein_S194.{py,out}`.

## What this closes, and what it does not (honest)

**Closes:** the "`𝒲̂` is numerically-certified, not a-priori" caveat in opus-S157 / THM-664 / LEM-009.
`𝒲̂(n)` is now an **exact elementary closed form**; the density-floor decorrelation constant `D3_∞`
and THM-664's residual `R` are explicit convergent sums with a proven geometric per-coordinate decay,
computable to any accuracy from finitely many low-height terms. This is the shared `𝒲̂`-decay the
fleet flagged (mac-mini THM-664 handoff (a), opus-S157).

**Does not close (and is not needed for LRC(14)):** a *single uniform* bound
`|Σ_{Vmax\mid n·e}𝒲̂(n)| < (6/7)^k` over **all** `(E,Vmax)` — the signed resonance sum still depends
on the cluster's additive relations (which small `n` satisfy `n·e=0` or `Vmax∣n·e`). But the exact
formula reduces even this to a **finite** low-height check per relation-type (the tail beyond height
`H` is `≤ Σ_{r,‖n‖>H} (0.371)^r/… ` explicitly bounded). And the large-spread **existence** — the only
thing LRC(14) needs — is already unconditional via **LEM-010** (elementary Dirichlet), so this lemma
is the *unification + abundance* (`#good ≈ (6/7)^k·Vmax`), not a proof-critical gap.

## Consequence for the covering case

The finite-`Vmax` glue's large-spread half is, by THM-664 + LEM-011, *exactly* the density-floor
decorrelation read on the `Vmax`-grid — one `𝒲̂`-resonance sum, now in closed form. Combined with
LEM-010's elementary closure and the density floor (THM-663), the covering case rests only on the
bounded finite check `{Vmax ≤ 3^{12}, spread ≥ 6Vmax/7}` and Lean.

## Independent cross-validation (mac-mini-2026-07-08-S59)

Derived independently the same day by a different route — direct swap-of-integrals giving
`𝒲̂(n) = (6/7)^z ∏_{n_i≠0} b(n_i)·Q(N)` with `b(m)=(1−e(m/7))/(2πim)`, `Q(N)=𝟙[N=0]·(6/7) +
(e(−N/7)−1)/(2πiN)` — identical to the statement above (`b = −b₀` sign convention; `Q(N)=𝟙[σ=0]−c(σ)`).
VERIFIED against **direct `T²`/`T³` Fourier integrals** (a check independent of klein's FFT+Parseval):
`|formula − direct| ≤ 2e-5` at `k=3,4`, `𝒲̂(0)=(6/7)^k` exact, and the `7`-commensurate zeros
`𝒲̂(7,0)=𝒲̂(1,6)=0` exact. Two independent derivations + two independent numerical checks agree.
File: `04-computation/lrc14_What_closedform_macmini_S59.{py,out}`. (This file superseded a duplicate
`LEM-011-uncovered-measure-fourier-closed-form.md` — same result, ceded to this one.)

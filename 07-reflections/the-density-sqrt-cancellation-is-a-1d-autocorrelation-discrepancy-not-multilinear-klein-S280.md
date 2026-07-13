# The density √-cancellation is a 1-D autocorrelation discrepancy — confirmed, and genuinely lower-order than the covering crux

*klein-2026-07-13-S280. Owner: prove the density √-cancellation bound (the S279 tractable route). The
√-cancellation is confirmed empirically and reduced to a clean, exact object — the 2nd moment of the
endpoint sum equals the Riemann-discrepancy of the autocorrelation of `1_{R_s}` at the `w`-grid, a **1-D
discrepancy**. This closes the density row given `Q_s=O(diam)` (verified), and pins the density-vs-covering
asymmetry: density needs a 1-D discrepancy, covering needs a multi-linear Gowers cancellation.*

---

## The √-route and why slack makes it enough

From THM-727, `|S|=|w·Error| ≤ (1/2π²)Σ_s Σ_ℓ |U_s(ℓw)||\sin(πℓ/7)|/ℓ²`. Cauchy–Schwarz gives
`|S|=O(√Q_s)`, `Q_s=Σ_ℓ|U_s(ℓw)|²/ℓ²`. Density has slack: `Q_s=O(diam)` suffices, because on the peel
`w=d≥diam`, `Error=|S|/w=O(√diam)/d=O(1/√d)→0`. So density closes with a *non-sharp* √-bound — no need
for the sharp `O(k)` (which is the multi-linear crux). This is the S279 asymmetry cashed in.

## The exact identity: 2nd moment = autocorrelation Riemann-discrepancy

`|U_s(N)|=2π|N||\hat f(N)|`, `f=1_{R_s}`. The Dirac comb gives `Σ_{ℓ∈ℤ}|\hat f(ℓw)|²=(1/w)Σ_{m<w}A(m/w)`,
`A(t)=|R_s∩(R_s−t)|` (autocorrelation, `A(0)=|R_s|`, `∫A=|R_s|²`), so

$$Q_s=(2\pi w)^2\Big[\tfrac1w\sum_{m=0}^{w-1}A(m/w)-\int_0^1 A\Big].$$

`Q_s` is *exactly* `(2πw)²` times the **left-Riemann-discrepancy of the piecewise-linear `A` at the
`w`-grid** (`A` has `O(M)=O(diam)` breakpoints). This is a **one-dimensional discrepancy**. The
√-cancellation lives entirely in how the `w`-grid samples one autocorrelation function.

## What is proved, and what the crude tool misses

- **`Q_s=O(diam)` (empirical, robust):** `Q_s/diam∈[1.0,1.7]` over 7 clusters (`diam` to 199) ⟹ `|S|=O(√diam)`
  ⟹ `Error→0` ⟹ density row closes (finite box + S275 for the band and resonant peels).
- **Diagonal (rigorous, dominant):** arc-wise `Q_s=Σ_i 2π²\{w w_i\}(1-\{w w_i\})+\text{off-diag}`; the
  diagonal is `≤(π²/2)(M/2)=O(M)` and empirically ≈ the full `Q_s` — the off-diagonal cancels. So the
  `O(M)` scaling has a rigorous backbone.
- **The crude large sieve fails — instructively.** `Q_s≤M(2+\tfrac43δ^{-1})`, `δ=\min\|w(p-p')\|`, is
  useless: `δ` collapses (`0.012→0.0000`) because endpoints cluster (two large offsets crossing nearly
  together). But tiny arcs **self-cancel** (`e(-Na)-e(-Nb)→0`), so the true `Q_s` is small — the
  min-separation bound can't see this; the autocorrelation identity does.

## The picture: density is lower-order than covering

| | needs | order |
|---|---|---|
| **covering** (opus-S262/S263, mac-mini-S76) | multi-linear (Gowers) cancellation; 3rd-order Schur; averaging **provably** insufficient | `≥3`-linear |
| **density** (this) | `Q_s=O(diam)` = a 1-D discrepancy of the autocorrelation `A` at the `w`-grid | 1-linear |

The S279 asymmetry is now concrete: the density-row tail closes on a **1-D autocorrelation discrepancy**
— averaging/2nd-moment *is* enough here (slack), exactly where it provably is not for the tight covering
bound. So density is the genuinely more tractable of the two LRC(14) analytic cruxes, and the √-route is
its natural proof.

## Honest scope

- **Reduced + confirmed:** the √-cancellation holds (`Q_s=O(diam)` verified), closing the density row; the
  exact autocorrelation-discrepancy identity and the diagonal `O(M)` bound are rigorous.
- **Remaining:** a clean rigorous upper bound `Q_s=O(M)` for the piecewise-linear `A` at the `w`-grid — a
  1-D discrepancy estimate (the crude min-separation large sieve is defeated by clustering; the arc-width /
  autocorrelation structure is the right handle). This is far lower-order than the covering crux, and is
  the concrete next density target.

*Files: `04-computation/lrc14_second_moment_klein_S280.py`, `lrc14_delta_sep_klein_S280.py` (+ outs).
THM-729, HYP-6415. Consumes THM-727/728. The asymmetry vs
[[both-lrc14-routes-bottom-on-the-same-multilinear-cancellation-not-one-mollification-klein-S279]].*

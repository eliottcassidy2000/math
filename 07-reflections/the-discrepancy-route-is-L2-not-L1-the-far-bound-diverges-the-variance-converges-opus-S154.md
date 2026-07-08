---
source: opus-2026-07-08-S154
status: the discrepancy route (LEM-005) for Var(W) resolved STRUCTURALLY -- the L^1 (far-correction
  absolute) bound provably DIVERGES (sum_m|ahat(m)| ~ (2/pi^2)ln M = infinity), the L^2 (variance)
  bound CONVERGES (Parseval: sum|ahat|^2 = theta(1-theta)=6/49; Var=sum_nu|What(nu)|^2 -> Var_exact
  verified all families). Sharpens mac-mini LEM-007 ("extreme cancellation") to "no absolute bound
  exists"; delineates klein LEM-009's Koksma-Hlawka scope (converges for few-outlier, needs L^2 for
  general spread). Exact far Fourier formula (star) derived + verified (= mac-mini's doubly-balanced
  support, explicit integral form).
tags:
  - lrc14
  - covering-floor
  - variance
  - discrepancy
  - fourier
  - L2-not-L1
  - lem-005
---

# The discrepancy route is L^2, not L^1: the far bound diverges, the variance converges

**opus-2026-07-08-S154.** Owner: work the discrepancy route (LEM-005) for `Var(W)`. LEM-005 reduces
the k=11 spread tail to `far <= E[W]^2` (<=> `Var(W) <= near <= (2/7)E[W]` => `PZ >= bar`) and proves
it AS A REDUCTION to the phase-vector discrepancy, leaving "the explicit discrepancy rate" as the
open analytic step. This note resolves the STRUCTURE of that step: it must be done in `L^2`, and the
naive `L^1` route is not merely lossy but literally divergent.

## The exact far Fourier formula (the doubly-balanced lattice, explicit)

With `theta=1/7`, the arc-just-before indicator `phi(u)=1(u in [-theta,0])`,
`ahat(m)=phihat(m)=(e(m theta)-1)/(2 pi i m)` (`ahat(0)=theta`), the two-disjoint-arc avoider
`Psi=1-phi(.)-phi(.-t)` (`Psihat(0)=5/7`, `Psihat(m)=-ahat(m)(1+e(-mt))`), the far correlation is

> **`far = (5/7)^{k+1} + sum_{m in L, m!=0} (5/7)^{k-|S|} (-1)^{|S|} (prod_{i in S} ahat(m_i)) J(m)`**,
> `L = { m in Z^k : sum m_i = 0 AND sum m_i e_i = 0 }` (the DOUBLY-BALANCED relation lattice),
> `J(m) = int_{1/7}^{6/7} prod_{i in S}(1+e(-m_i t)) dt = sum_{T subset S} K(sum_{i in T} m_i)`.

Verified to machine precision against the exact Farey `far = E[W^2]-near` (imag parts `~1e-85`). The
`sum m_i=0` constraint comes from the `y`-average, `sum m_i e_i=0` from the `x`-average -- so the
support is EXACTLY mac-mini LEM-007's "doubly-balanced resonances" (independent, same session). This
`(star)` is the explicit integral form; support `>=3` (no support-2: `m_a e_a + m_b e_b=0`,
`m_a+m_b=0` forces `e_a=e_b`), leading term support-3 3-APs `m=(1,-2,1)` (mac-mini's correction).

## The L^1 obstruction: the absolute bound DIVERGES (rate (2/pi^2) ln M)

Taking absolute values in `(star)` with `|J(m)| <= (5/7)2^{|S|}`, `|ahat(m)| <= 1/(pi|m|)` gives
`|far-(5/7)^{k+1}| <= (5/7)^{k+1} PHI`, `PHI = sum_{m in L} prod (14/5)|ahat(m_i)|`. **PHI diverges:**
numerically `1.15 -> 2.51 -> 5.91 -> 10.70 -> 23.54 -> 32.93` as the (support, frequency) cutoff grows.
The rigorous reason is a ONE-LINE fact about the arc:

> **`sum_{m=1}^{M} |ahat(m)| = sum_{m=1}^M |sin(pi m theta)|/(pi m) ~ (2/pi^2) ln M -> +infinity`.**

(Verified: `0.59, 1.04, 1.50, 1.96, 2.42` at `M=10..1e5`, tracking `(2/pi^2)ln M`.) The `theta`-arc
indicator is bounded-variation but NOT absolutely Fourier-summable; a PRODUCT of `k` arcs inherits
this, so ANY term-by-term absolute bound on the far correction is infinite. This is the exact
mechanism of LEM-005's "the 2/7-arcs are too full, `S_1 = 2k/7 > 1`" and it SHARPENS mac-mini
LEM-007's "extreme cancellation across all orders": it is not that a large convergent sum cancels to
something small -- the underlying absolute sum is `+infinity`, so cancellation is not optional, it is
the only thing standing between us and a divergent series. **There is no `L^1` bound to be had.**

## The L^2 resolution: the variance CONVERGES (Parseval, constant 6/49)

The companion of the divergent `L^1` fact is an exact FINITE `L^2` identity:

> **`sum_{m!=0} |ahat(m)|^2 = theta(1-theta) = 6/49`** (Parseval for `1_arc`: `||1_arc-theta||_2^2`).

and the variance itself is the manifestly-convergent square-sum

> **`Var(W) = || W - E[W] ||_2^2 = sum_{nu != 0} |What(nu)|^2`,
> `What(nu) = sum_{m: sum m_i=0, sum m_i e_i = nu} (6/7)^{k-|S|} prod_{i in S}(-ahat(m_i))`.**

Truncating the `m`-cutoff, `sum_nu|What(nu)|^2 -> Var_exact` for EVERY family (ratios `1.00-1.02` at
`Mmax=4,Smax=6`), INCLUDING the compact block (`ratio 1.03`) where the `L^1`/support expansion is
non-perturbative. `Var` is dominated by `nu = ` small differences (block: `nu=+-2` is `16%` each);
the difference multiplicities are the additive energy -- this is mac-mini LEM-007 / klein THM-656
`Var ~ R2` from the `nu!=0` Fourier side. **The variance is the correct convergent home of the
discrepancy route; `far` must be reached through `Var` (via the rigorous `far <=> Var <= near`
equivalence, mac-mini LEM-007), never bounded directly.**

## Delineating klein LEM-009: why Koksma-Hlawka converges for few-outlier, needs L^2 for spread

klein LEM-009 proves `D3(block_10 u {D}) -> 0.4646 >= bar` at rate `O(1/D)` via Koksma-Hlawka, using
the SAME per-entry factor `(14/5)|ahat(n)| <= 0.891/|n| < 1` and "geometric decay in support." This
converges because block+outlier has ONE far point: every resonance that MOVES the outlier sits at
frequency `>= ~D`, and there is a single such generator, so the `L^1` sum over outlier-resonances is
effectively geometric (finite). The present divergence says precisely where this stops: for a GENERAL
spread family ALL `k` points are far, the doubly-balanced lattice `L` is dense, and the same
`0.891/|n|` per-entry bound summed over `L` is my divergent `sum|ahat|`. So:

> **klein's `L^1` Koksma-Hlawka rate is available exactly when the resonance lattice is generated by
> FEW far frequencies (block + `O(1)` outliers); the general-spread tail is genuinely `L^2` (variance
> cancellation), not `L^1`.**

This is consistent with, and explains, klein-S186's reduction of the tail to CLUSTER-monotonicity:
the tractable (Koksma-Hlawka) families are exactly cluster + few outliers; the residue is the
cluster's own internal variance, which is bounded by the cluster limit (`D3_c`) / the exhaustive.

## Ledger

- DERIVED + VERIFIED: the exact far Fourier formula `(star)` (explicit integral form of mac-mini
  LEM-007's doubly-balanced support); `sum_m|ahat(m)| ~ (2/pi^2)ln M = infinity` (L^1 obstruction);
  `sum|ahat|^2 = 6/49` and `Var = sum_nu|What(nu)|^2 -> Var_exact` (L^2 convergence, all families).
- SHARPENS: mac-mini LEM-007 "extreme cancellation" -> "the absolute sum is `+inf`, no `L^1` bound
  exists"; the discrepancy route is `L^2`.
- DELINEATES: klein LEM-009 Koksma-Hlawka converges for few-outlier (lowest resonance `~D`), diverges
  for general spread (dense `L`) -> supports klein-S186 cluster-monotonicity (tractable = cluster +
  few outliers).
- CORROBORATES (independent): mac-mini LEM-007 doubly-balanced `far <=> Var <= near` ONE wall;
  klein THM-656 / LEM-007 `Var ~ R2`; the `0.891/|n|` factor (klein LEM-009).
- Files: `lrc14_far_fourier_discrepancy_opus_S154.py` (`(star)` verify), `lrc14_far_discrepancy_
  criterion_opus_S154.py` (PHI divergence), `lrc14_discrepancy_L2_convergent_opus_S154.py` (L^1 div
  rate + L^2 convergence). Builds on LEM-005 (klein), LEM-007 (mac-mini), LEM-009 (klein), opus-S152
  (`c_j`)/S153 (relation lattice `E_2=R2`).

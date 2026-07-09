---
id: THM-670
title: The Lipschitz threshold transfer — for θ ≥ 1/7 at most 6 gaps pay, so W_θ is 6-Lipschitz in θ and E[W_{θ₂}] ≥ E[W_{θ₁}] − 6(θ₂ − θ₁); the closed 1/7-leg first moments transfer to every lifted threshold with an explicit slope
status: PROVED (four lines, below). Machine-verified: the parametric ledger (companion .out) never violates the transfer, headrooms reported per (n, r).
source: monad-explorer-2026-07-09-S5 (HYP-5747).
depends_on: []
related:
  - THM-669   # consumes E[W_{1/7+r}] floors; this transports them from θ = 1/7
  - THM-661   # the 1/7-leg moment machinery whose outputs transfer
---

# THM-670 — the Lipschitz threshold transfer

## Statement

For any gap configuration (finitely many circular gaps summing to 1) and any
`1/7 ≤ θ₁ ≤ θ₂`:

> `W_{θ₂}(x) ≥ W_{θ₁}(x) − 6·(θ₂ − θ₁)`   pointwise, hence
> **`E[W_{θ₂}] ≥ E[W_{θ₁}] − 6·(θ₂ − θ₁)`.**

## Proof

`W_θ(x) = Σ_i (g_i(x) − θ)₊` is convex and nonincreasing in `θ`, with
`∂W_θ/∂θ = −#{i : g_i > θ}` at every non-kink `θ`. If `m` gaps exceed `θ ≥ 1/7`
then `m·(1/7) < Σ g_i ≤ 1`, so `m ≤ 6`. Integrating from `θ₁` to `θ₂`:
`W_{θ₁} − W_{θ₂} ≤ 6(θ₂ − θ₁)`. Average over `x`. ∎

## Remarks

1. **Use.** THM-669's availability floors need `E[W_{1/7+r}]`; any proved 1/7-leg
   first-moment floor `E[W_{1/7}] ≥ c_n` becomes `E[W_{1/7+r}] ≥ c_n − 6r` for free.
   The parametric ledger (companion) shows the transfer is conservative by roughly a
   factor 2–3 in slope on compact minimizers (the true count of paying gaps at the
   minimum is 2–3, not 6) — the direct parametric scan is sharper where it matters.
2. **Scope.** The same argument gives slope `⌈1/θ₁⌉ − 1` for any `θ₁ > 0`; at
   `θ₁ = 1/7` this is 6. No transfer is claimed for `μ_θ` (a level-set measure can
   drop abruptly); the transfer is a first-moment statement, which is exactly what
   THM-669 consumes.

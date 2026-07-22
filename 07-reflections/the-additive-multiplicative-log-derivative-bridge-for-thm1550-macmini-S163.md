---
source: mac-mini-2026-07-22-S163
status: reformulation — THM-1550 obstacle (iii) (the additive<->multiplicative bridge, the last DvdK gap) recast as a LOG-DERIVATIVE identity, avoiding the formal log; hands death-star a simpler formalization target
tags:
  - GMC2
  - DvdK
  - THM-1550
  - additive-multiplicative
  - log-derivative
  - lonely-runner
---

# The additive↔multiplicative bridge for THM-1550, without the log

Owner: *work to close DvdK or bypass it, think additive vs multiplicative.* The whole GMC(2)
formalization is now sorry-free and reduced (via THM-2067: Galois wrapper + Vieta + irreducibility,
all landed) to a single analytic input — **THM-1550**, whose deep core is death-star's **obstacle
(iii)**: the additive series `Σ D_m t^m` versus the multiplicative small-root product `Π(t)`.
death-star's route bridges them with "a formal CT_u of `log` of the Wiener–Hopf factorization."
This note replaces the log with a **log-derivative**, which is cleaner and more Lean-shaped.

## Setup

`Φ(x) = x^M − t·R(x)`, `R` of degree `d = M+N`, `R(0) = r₀ ≠ 0`, `M ≥ 1`. `Π(t)` = product of the
`M` small roots (valuation `> 0`); `D_m = [x^{Mm}] R^m` (so `D₀ = 1`, and `D_m = CT(f^m)` for
`f = R·x^{−M}`). These `D_m` are the DvdK constant terms; the DvdK hypothesis being contradicted
is `D_m = 0 ∀m ≥ 1`.

## The identity

In `ℂ[[t]]`:
```
    t · Π'(t) / Π(t)  =  F(t) := Σ_{m≥0} D_m t^m  =  [x⁰]  x^M / (x^M − t R)      (root-free).
```
Verified to ~1e-11 across `M,N ∈ {1,2,3}` (script + .out). **Immediately**:
`D_m = 0 ∀m ≥ 1  ⟹  t Π'/Π = 1  ⟹  Π = c·t` — which is exactly obstacle (iii) (`Π` rational of
`t`-valuation 1, hence Galois-fixed, feeding boxeph's orbit-product/valuation contradiction).

## Why it is simpler: it factors into two elementary pieces

- **(a) Per-root — pure calculus.** Each small root `u_j` satisfies `u_j^M = t·R(u_j)`. Differentiate
  in `t`: `M u_j^{M−1} u_j' = R(u_j) + t R'(u_j) u_j'`, solve, and substitute `u_j^M = t R(u_j)`:
  ```
        t · u_j'/u_j  =  R(u_j) / S(u_j),        S := M·R − x·R'.
  ```
  (One line; no analysis. Verified.)
- **(b) Sum — a symmetric function / residue.** `Σ_{small} R(u_j)/S(u_j) = F(t)` (verified to ~1e-16).
  Equivalently, using `S(ξ) = ξ Φ'(ξ)/t` at a root: `F(t) = t·Σ_{small} Res_{u_j}[R/(x Φ)]`, with the
  **root-free anchor** `Res₀[R/(x Φ)] = R(0)/Φ(0) = −1/t`.

Then `t Π'/Π = Σ_j t u_j'/u_j = Σ_j R(u_j)/S(u_j) = F(t)`. So obstacle (iii) needs only:
*differentiate the defining relation* (a) + *a symmetric-function sum* (b) — **no formal log, no
Wiener–Hopf factorization over a valued field.**

## Scope (honest)

This simplifies obstacle (iii) only. **Obstacle (ii)** — constructing the `M` small roots `u_j`
(Puiseux in `t^{1/M}`) / the Weierstrass factor — remains the substantial piece (Mathlib lacks a
Henselian *factorization* theorem; death-star's monic M-th-root Hensel + fixed-point route stands).
The `F(t) = Σ D_m t^m = [x⁰] x^M/(x^M−tR)` half is fully root-free and formal (a `PowerSeries`
identity, no Puiseux) — a clean reusable lemma; the per-root formula (a) needs the roots but is
trivial once they exist. Does **not** by itself close THM-1550 / DvdK.

## The additive↔multiplicative moral

`log Π = log(c·t) + Σ_{m≥1} D_m t^m/m`: the **multiplicative** root product `Π` exponentiates the
**additive** DvdK series `Σ D_m t^m/m`; the log-derivative turns the product into the sum
`Σ u_j'/u_j` and the per-root relation `u^M = tR(u)` turns *that* into the root-free `F(t)`. The
DvdK vanishing `D_m = 0` is precisely "the multiplicative packet `Π` collapses to the monomial
`c·t`," which is the valuation-1 rational element the Galois argument needs.

→ hands death-star (THM-1550 owner) a simpler obstacle-(iii) target; builds on death-star-S106/S111
(Wiener–Hopf/log, Hensel), boxeph-S235/S236 (THM-2067 wrapper, orbit-product, done), mac-mini-S162
(Φ irreducibility). Script: `04-computation/dvdk_thm1550_logderivative_bridge_macmini_S163.py` (+out).

---
id: THM-1795
title: "Four-cell conjecture for generic GMC(2) recurrence order"
status: >
  CONJECTURAL / FINITE-EXACT ON FOUR GENERIC CELLS. Computations observe
  recurrence orders M+N for d=0 at spans 2,3,4 and order 6 at
  (M,N,d)=(1,1,1), matching (M+N)(2d+1). The general s-degree claim,
  exact creative-telescoping order, and minimality are unproved. More
  importantly, recurrence order is not automatically a uniform detection
  depth: leading coefficients may vanish at early indices or on parameter
  loci. Therefore this file supplies neither an exact depth formula nor an
  explicit Gröbner-test size. THM-1790 proves the actual unbounded
  degree-cap lower bound; THM-2022 proves NC2 by a different mechanism.
source: kind-pasteur-2026-07-20-S128c130 (owner: work the GMC(2))
depends_on:
  - THM-1770    # klein: detection depth grows with radial degree, no uniform bound
  - THM-1710    # TNC detection depth = span (the d=0 base)
related: [THM-1740, THM-1670, THM-1645, THM-1790, THM-2022, MISTAKE-247]
script: 04-computation/gmc2_depth_vs_degree_kps_S128c130.py, gmc2_depth_formula_kps_S128c130.py (+ .out)
---

# THM-1795 — a four-cell recurrence-order conjecture

> **CORRECTION (2026-07-24; MISTAKE-247).** The displayed formula is a
> conjecture about generic minimal recurrence order, verified on four cells.
> The original file incorrectly promoted it to an exact detection-depth
> theorem. A P-recurrence yields a zero-propagation bound only after its
> leading coefficients are proved nonzero throughout the required initial
> range and uniformly on the coefficient stratum. Those hypotheses are not
> established here. All theorem-like consequences below are replaced by the
> conditional statements in this repaired version.

The candidate formula is

> **`D_rec^gen(M,N,d) ?= (M+N)(2d+1)`.**

## What the numbers mean

Let `P ∈ ℂ[Z, Z̄]` have charge span `[−M,N]` and radial degree at most
`d`. The experiments measured generic recurrence orders for specialized
moment sequences. They did not compute a uniform symbolic recurrence on the
whole coefficient stratum and did not prove that recurrence order equals
detection depth.

| `(M,N)` | `d` | `D` | `(M+N)(2d+1)` | status |
|---|---|---|---|---|
| (1,1) | 0 | 2 | 2 | ✓ verified |
| (1,2) | 0 | 3 | 3 | ✓ verified |
| (1,3) | 0 | 4 | 4 | ✓ verified |
| (1,1) | 1 | 6 | 6 | ✓ verified |
| (1,1) | 2 | 10 | 10 | derivation; compute-limited |
| (1,2) | 1 | 9 | 9 | derivation; compute-limited |

The `d = 0` row is `order = span`, matching THM-1710 (TNC depth) and klein's "depth = span at
`d = 0`". The jump `2 → 6` from `d = 0` to `d = 1` at `(1,1)` is the `(2d+1)` factor.

## The proposed derivation and its missing gates

`E[P^m] = L_s(CT_u[Λ_s^m])` — the Laplace transform in the radius `s = |Z|²` of the toral
diagonal (klein THM-1770). The heuristic requires:

1. the toral diagonal `CT_u[Λ_s^m]` to have order `M+N` in `m`, with recurrence
   coefficients that are polynomials in the symbol coefficients `g_c(s)`.
2. those coefficients to have `s`-degree at most `2d`. Each `g_c` has `s`-degree `d`; the recurrence
   coefficients are quadratic in the `g_c` (the shape is visible at `(1,1)`, where the trailing
   coefficient is `disc(R) = g₁² − 4g₀g₂`, degree 2 in the `g`'s). This does
   not prove the same degree law at arbitrary span.

Integrating `s` out against `e^{−s}` (the Gaussian radial average) is a **creative-telescoping**
of one variable. Even if it produces an order bound `ρ(1+σ)`, exact equality
requires a minimality argument. Converting the resulting recurrence into a
uniform moment cutoff additionally requires its leading coefficient to avoid
zeros at all initial indices and on all parameter loci under consideration.
Neither gate is proved. Thus the calculation

> `ρ(1+σ) = (M+N)(2d+1)`

is a candidate scale, not a theorem.

## Why it matters

- THM-1740's bounded-stratum finite-elimination statement remains independent
  of this conjecture, but this file does not give its test size.
- THM-1790 now gives the rigorous scale information: at fixed charge span two
  and radial degree cap `d`, detection depth is at least `2d+2`. Thus a
  degree-uniform cutoff is impossible. The candidate formula, if repaired,
  would instead propose an upper scale of the same linear order.

## Named next

- **Break the compute wall** to verify `(1,1,2) = 10`, `(1,2,1) = 9`, `(2,2,1) = 12` directly:
  the moment engine needs the toral+Laplace factorization (compute `CT_u[Λ_s^m]` as an
  `s`-polynomial, then apply `Σ_k k!·[s^k]`) rather than full `P^m` dict-convolution, which
  reaches larger `m` per unit cost.
- **Prove `σ = 2d`** (the coefficient `s`-degree of the toral recurrence) in general — the one
  ingredient of the derivation argued only from the `(1,1)` shape. It should follow from the
  `s`-formula (THM-1690) applied with `R`-coefficients of `s`-degree `d`.
- **Prove a uniform nonsingular recurrence.** Track its leading coefficient
  in both `m` and the polynomial coefficients, including exceptional loci.
- **Separate recurrence order from detection depth.** Establish a
  zero-propagation lemma with the exact initial interval before using an order
  as a Gröbner cutoff.

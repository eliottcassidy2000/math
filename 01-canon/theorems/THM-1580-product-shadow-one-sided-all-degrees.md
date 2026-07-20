---
id: THM-1580
title: "THE MOMENTS SEE ONLY THE PRODUCT: one-sided PROVED AT EVERY DEGREE (charge span {-1,0,+1}, B=0), upgrading THM-1570's degree<=1 elimination to a structural theorem. HYP-8390 RESOLVED: the cross-term matrix is M_{ij} = (i+j+1)!, the HANKEL MOMENT MATRIX OF THE GAMMA(2) MEASURE s e^{-s} ds, hence POSITIVE DEFINITE at every size (2M = [[2,4],[4,12]] at N=2, matching THM-1570's observed cross term exactly) — but nonsingularity ALONE is NOT enough, since A^T M C = 0 does not force A=0 or C=0 (take M=I, A=(1,0), C=(0,1)), so the second moment can never finish. THE STRUCTURAL REDUCTION does finish it: a product of m factors has total charge 0 iff it uses p copies of zA, p of zbar C and r = m-2p of B, giving E[P^m] = sum_{2p+r=m} m!/(p!p!r!) E[h^p B^r] with h := s*A*C (verified symbolically through m=6). SO THE NULLCONE CONDITION SEES A AND C ONLY THROUGH THEIR PRODUCT — it is invariant under A -> lambda A, C -> lambda^{-1} C. With B = 0 the sum collapses to E[P^{2p}] = C(2p,p) E[h^p], so the nullcone forces E[h^p] = 0 for all p, whence h = 0 by THM-1570 sB (the Laplace layer), whence s*A*C = 0 and — C[s] being a DOMAIN — A = 0 or C = 0, AT EVERY DEGREE. This is the same algebraic move as the Vieta e/pi argument, and it realises the repo's multiplication/addition duality: h is the MULTIPLICATIVE shadow, B the ADDITIVE one, and they DECOUPLE at B = 0"
status: PROVED (the B=0 case at all degrees; HYP-8390 resolved) + VERIFIED (structural reduction symbolic through m=6). The B != 0 case remains.
author: opus-2026-07-20-S414
depends_on: [THM-1570 (the Laplace layer — the engine here), THM-1540, THM-1535]
resolves: HYP-8390
---

# THM-1580 — The moments see only the product

## 1. HYP-8390 resolved: the cross-term matrix

For charge span `{−1,0,+1}`, `P = z·A(s) + B(s) + \bar z·C(s)` with `s = z\bar z`. The cross
term of the second moment is

```
2·E[ (zA)(\bar z C) ] = 2·E[ s·A(s)·C(s) ] = 2·AᵀMC ,     M_{ij} = E[s^{i+j+1}] = (i+j+1)!
```

`M` is the **Hankel moment matrix of the Gamma(2) measure `s·e^{−s}ds`**, since
`∫₀^∞ s^k·s e^{−s}ds = (k+1)!`. A Hankel matrix of a positive measure with infinite support
is **positive definite**, hence nonsingular, at every size (min eigenvalues checked to
`N = 7`). At `N = 2`, `2M = [[2,4],[4,12]]` — **exactly** THM-1570's observed cross term
`2a₀c₀+4a₀c₁+4a₁c₀+12a₁c₁`.

> **But nonsingularity alone is not enough.** `AᵀMC = 0` with `M` nonsingular does *not*
> force `A = 0` or `C = 0` — take `M = I`, `A = (1,0)`, `C = (0,1)`. **The second moment can
> never finish this proof.** HYP-8390 answered, and answered *negatively as a strategy*.

## 2. The structural reduction — and what the moments can see

A product of `m` factors of `P` has total charge `0` **iff** it uses `p` copies of `zA`, `p`
copies of `\bar zC`, and `r = m−2p` copies of `B`. That product is
`(zA)^p(\bar zC)^p B^r = s^pA^pC^pB^r = h^p B^r` with **`h := s·A·C`**. Hence

```
E[P^m]  =  Σ_{2p+r=m}  m!/(p!·p!·r!) · E[ h^p B^r ]
```

*(verified symbolically through `m = 6` at degree 2 — identity holds exactly).*

> **The nullcone condition sees `A` and `C` only through their PRODUCT.** It is invariant
> under `A ↦ λA`, `C ↦ λ⁻¹C` for any `λ ∈ ℂ*`, which leaves `h` fixed. The individual
> factors are invisible to every moment.

## 3. The all-degree theorem (`B = 0`)

With `B = 0` only `r = 0` survives:

```
E[P^{2p}] = (2p)!/(p!p!) · E[h^p] ,        E[P^{odd}] = 0.
```

So the nullcone forces `E[h^p] = 0` for **every** `p ≥ 1`. By **THM-1570 §B** (the Laplace
layer: the only `g ∈ ℂ[s]` with `∫₀^∞ g^p e^{−s}ds = 0` for all `p` is `g = 0`), we get

```
h = 0  ⟹  s·A(s)·C(s) = 0  ⟹  A·C = 0  in ℂ[s]  ⟹  A = 0  or  C = 0.
```

the last step because **`ℂ[s]` is an integral domain**. ∎

> **One-sided, at EVERY degree** — replacing THM-1570 §A's degree-≤1 Gröbner result, which
> stalled at 9 unknowns. The elimination was solving the wrong problem; the grading solves
> it outright.

## 4. Multiplication/addition duality — the connection is structural, not decorative

The repo already carries this duality: *"these ARE the multiplicative and additive poles of
the sum-product structure"*, and oracle-S555o's *"addition-shadow and multiplication-shadow
COINCIDE on the rationals; the FINE pinch is where addition outruns divisibility"*.

Here it appears cleanly:

| | object | role |
|---|---|---|
| **multiplicative shadow** | `h = s·A·C` | the product of the two charged parts |
| **additive shadow** | `B` | the charge-0 (uncharged) part |

§2 says the moments are a function of `(h, B)` alone, and §3 says **the two shadows decouple
at `B = 0`**: with no additive part, the multiplicative shadow carries everything, and its
vanishing kills a factor.

**The Vieta rhyme is exact.** The classical argument — if `e+π` and `eπ` were both algebraic,
then `e, π` are roots of `x² − (e+π)x + eπ` over `ℚ̄`, hence algebraic, contradicting
Lindemann — controls the **sum and the product** and derives a contradiction. Here we control
**only the product** `AC`, and the integral-domain step `AC = 0 ⟹ A = 0 or C = 0` is the same
algebraic move: *the elementary symmetric data determines the factors up to the domain
property.* §3's proof is, in that sense, a one-sided Vieta argument.

**Schanuel: no connection found, and I looked.** Schanuel's conjecture concerns the
transcendence degree of `ℚ(x₁,…,x_n, e^{x₁},…,e^{x_n})`; nothing in the moment/charge setup
produces exponentials of independent transcendentals, and the Vieta step above is elementary
field theory, not Schanuel-strength. Recording the negative rather than manufacturing a link
(HYP-8230 rule).

## 5. What remains: `B ≠ 0`

The generating function of §2 sums to

```
E[e^{tP}] = E[ e^{tB(s)} · I₀( 2t√(h(s)) ) ]  =  1        (I₀ = modified Bessel, Σ x^p/(p!)²)
```

As `t → ∞`, `I₀(x) ∼ e^x/√(2πx)`, so the integrand behaves like `e^{t(B + 2√h)}`. If either
`B` or `h` is non-constant then `|B + 2√h| → ∞` as `s → ∞` while its argument stabilises to
that of the leading coefficient — **the same phase-stabilisation mechanism as THM-1570 §B** —
so the integral cannot stay bounded, let alone equal `1`. That forces `B` and `h` constant,
and then expanding `e^{tb₀}I₀(2t√h₀) = 1` in `t` gives `b₀ = 0` and `h₀ = 0`.

**This is a proof sketch, not a proof:** the uniformity of the saddle estimate for complex
`B, h` is exactly the step THM-1570 §B established in the single-function case and which
needs redoing for the `e^{tB}I₀(2t√h)` combination. Filed as HYP-8405.

## 6. Ledger

| locus | status |
|---|---|
| charge-definite, any `n` | PROVED (THM-1535 §1) |
| sign-coherent over `ℝ`/`ℂ` | PROVED (THM-1535 §2, THM-1540 §3, THM-1570 §B) |
| span `{−1,0,+1}`, deg ≤ 1 | PROVED by elimination (THM-1570 §A) |
| **span `{−1,0,+1}`, `B = 0`, ALL degrees** | **PROVED** (§3) |
| span `{−1,0,+1}`, `B ≠ 0` | sketch only (§5) — HYP-8405 |
| wider charge spans | open |

## Verification

`04-computation/cross_matrix_and_all_degrees_opus_S414.py` — `M`'s positive-definiteness to
`N = 7`; the structural reduction verified symbolically through `m = 6`; the `B = 0` collapse
checked term by term. Output in `05-knowledge/results/`.

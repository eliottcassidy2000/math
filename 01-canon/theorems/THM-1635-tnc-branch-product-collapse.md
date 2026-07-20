---
id: THM-1635
title: "THE GENERATING FUNCTION COLLAPSES TO A BRANCH PRODUCT: G(t) = -t (log Pi)' EXACTLY, where Pi(t) = product of the N small branches of u^N = t R(u); so TNC <=> Pi(t) = c t, rederiving klein's criterion with Pi identified explicitly. (1) THE COLLAPSE, exact, no asymptotics: G(t) = sum_m CT(Lambda^m) t^m = sum over small branches of R(u_i)/S(u_i), and along a branch t = u^N/R(u) gives t'(u) = -u^{N-1}S(u)/R^2, so R(u_i)/S(u_i) = -t d(log u_i)/dt and G = -t d/dt log Pi. (2) THE VIETA CONSTANT: for M >= 1 the product of ALL M+N roots of u^N - tR(u) is (-1)^{M+N} r_0/r_{M+N}, CONSTANT in t (verified -2, -1 on the M=2 cases; NON-constant exactly when M=0, the extreme-weight case THM-1530(B) handles separately). So Pi(t)*Pi_large(t) = const, and TNC <=> Pi = c t <=> t*Pi_large = const. (3) THE SINGULARITY CONFIRMATION: G is algebraic with radius 1/rho set by the dominant saddle (THM-1615), and on every collision case the dominant saddles are NONDEGENERATE (g'' = 2.31, 2.83, 2.67 != 0) even though their VALUES collide, so G is genuinely singular at 1/rho and CT != 0 -- the THM-1625 collision cases all satisfy TNC robustly. (4) THE RESIDUAL, now EXACT-ALGEBRAIC not asymptotic: is there a NON-monomial R with R(0) != 0 whose small-branch product Pi(t) is exactly linear in t? Equivalently t*Pi_large(t) constant. No asymptotic tower, no second-order Vandermonde needed -- a single algebraic condition on the branch structure of u^N = tR(u)"
status: PROVED (1),(2); VERIFIED (3); (4) is the exact residual, replacing THM-1625's asymptotic framing
author: opus-2026-07-20-S418
depends_on: [THM-1625 (the saddle split), THM-1615 (genuine saddle), THM-1530 (klein: extreme-weight + Pi(t)=ct), THM-1550 (klein: TNC <=> Pi=ct)]
---

# THM-1635 — The generating function collapses to a branch product

## 0. Why this supersedes the asymptotic route

THM-1625 attacked the TNC residual through the `m`-asymptotics of `CT(Λ^m)`
("second-order Vandermonde"). Attempting that numerically **failed** — the dominant saddles
are complex-conjugate, so `CT` oscillates and no clean power `m^p` is visible. That failure
is a signal: the right object is not the asymptotic series but the **exact generating
function**, which is algebraic and collapses to a single branch product. This note replaces
the asymptotic residual with an exact algebraic one.

## 1. The collapse (PROVED, exact)

`G(t) := Σ_{m≥0} CT(Λ^m) t^m = Σ_m [u^{Nm}]R^m t^m`. Summing the geometric series inside the
constant-term contour and taking residues at the `N` **small branches** `u_i(t) → 0` of
`u^N − tR(u) = 0`:

```
G(t)  =  − Σ_i  R(u_i) / S(u_i) ,      S(u) = uR′(u) − N R(u).
```

Along a branch, `t = u_i^N/R(u_i)`, so `t′(u) = −u^{N−1}S(u)/R(u)²`, whence
`R(u_i)/S(u_i) = −t/(u_i t′(u_i)) = −t·d(log u_i)/dt`. Summing,

```
G(t)  =  − t · d/dt  log Π(t) ,      Π(t) := ∏_{i=1}^{N} u_i(t)   (product of small branches).
```

**No asymptotics.** `Π` is an algebraic function of `t` with `Π(0) = 0`.

## 2. The Vieta constant, and the criterion (PROVED)

For `M ≥ 1` the polynomial `u^N − tR(u)` has degree `M+N` in `u` with leading coefficient
`−t·r_{M+N}` and constant term `−t·r_0`, so the product of **all** `M+N` roots is

```
∏_{all} = (−1)^{M+N} · (−t r_0)/(−t r_{M+N}) = (−1)^{M+N} r_0/r_{M+N}  —  CONSTANT in t.
```

*(Verified: `−2` for `R=u⁴−2u²−2`, `−1` for `R=1+u+u²+u³`. And for `M = 0` — the
extreme-weight case — the product is **not** constant, e.g. `−t/(t−1)`; that case is proved
separately by THM-1530(B).)*

Splitting the roots into `N` small and `M` large, `Π(t)·Π_large(t) = const`. Since
`G(0) = 1` and `G = −t(log Π)′`,

> **TNC  ⟺  G ≡ 1  ⟺  Π(t) = c·t  ⟺  t·Π_large(t) = const.**

This **is** klein's `Π(t) = ct` (THM-1550), now with `Π` pinned down as the **product of the
small branches** and dually as `const/Π_large`.

## 3. Singularity confirmation: the collision cases satisfy TNC (VERIFIED)

`G` is algebraic and analytic at `0`; its radius of convergence is `1/ρ`, `ρ = max|w_j|` the
dominant saddle value (THM-1615 guarantees a genuine saddle, so `ρ > 0` and `G` is
non-constant unless every dominant contribution cancels to all orders). The question is
whether the dominant singularity is genuine.

**On every THM-1625 collision case the dominant saddles are NONDEGENERATE:**

| `R` (N=2) | dominant `|g″(u_j)|` | degenerate? |
|---|---|---|
| `u⁴−2u²−2` | 2.309 | no |
| `u⁴−2` | 2.828 | no |
| `u⁴+u²−2` | 2.667 | no |

So even where the saddle **values** collide (THM-1625), the saddles themselves are simple.
A nondegenerate saddle contributes a genuine `√(t_j−t)` (or higher Puiseux) singularity that
cannot be removed unless its whole germ vanishes — impossible for `g″ ≠ 0`. Hence `G` is
singular at `1/ρ`, `CT ≠ 0`, and **all collision cases in the sweep satisfy TNC**. The
ratio test confirms it: `CT(m+1)/CT(m) → ρ` numerically.

> **Correction to the natural guess.** Collision of saddle *values* is **not** degeneracy of
> the saddle. THM-1625 conflated the two as "the residual"; §3 separates them, and only true
> degeneracy (`g″ = 0`) could threaten TNC.

## 4. The residual, now exact-algebraic

> **HYP-8465′ (sharpened).** Is there a **non-monomial** `R` with `R(0) ≠ 0` whose small-branch
> product `Π(t)` is **exactly linear** in `t` (equivalently `t·Π_large(t)` constant)?

This is a single algebraic condition on the branch structure of `u^N = tR(u)` — **no
asymptotic tower, no second-order Vandermonde, no case split on symmetric vs asymmetric
collisions.** Every prior framing (prefactor cancellation, degenerate saddles, symmetry
descent) is subsumed: they are all ways `Π` could or could not be linear.

**Evidence it cannot:** `Π = c·t` forces `t·Π_large` constant, but `Π_large ~ (t r_{M+N})^{−1}`
with `M` large branches each `~ (t r)^{−1/M}` — their product is `c/t` only to leading order,
and matching *all* orders is what a non-monomial `R` cannot do in any tested case (`CT`
nonzero for `1+u+u²`, `u⁴−2u²−2`, etc.). A proof that the higher-order terms never conspire
is the last step.

## 5. Status of the TNC (updated)

| case | status |
|---|---|
| `M = 0` (extreme weight `−1`) | **PROVED** (THM-1530 B) |
| `min(M,N)=1`, `(2,2)`, `(2,3)` | **PROVED** (THM-1595) |
| `(2,4)`, `(3,3)` | elimination (THM-1595) |
| dominant values distinct | **PROVED** (THM-1625 §1) |
| value-collision, nondegenerate saddle | **holds** (§3) |
| general `M,N` | ⟺ `Π(t)` linear for no non-monomial `R` (§4) — exact residual |

## 6. Next

1. **Prove `Π(t)` linear ⟹ `R` monomial.** `Π` and `Π_large` are the two factors of
   `u^N − tR` under the small/large Newton-polygon split; `t·Π_large = const` is a statement
   about the `M` large branches, governed by the **top** part of `R` (degree `M+N` down to
   `N`). This is a genuinely smaller object than `R` — the reduction the ladder lacked.
2. **The Newton-polygon factorisation of `Π_large`.** The large branches are the Puiseux
   solutions of `u^M ≈ 1/(t r_{M+N})·(1 + …)`; their product is an `M`-fold resultant. Its
   `t`-dependence beyond `1/t` is the exact obstruction.

## Verification

`04-computation/tnc_branch_product_opus_S418.py` (the collapse identity; the Vieta constant,
`−2`/`−1` for `M=2`, non-constant for `M=0`; the `Π = ct` test on non-monomial `R`),
`tnc_exact_gf_singularity_opus_S418.py` (nondegeneracy of dominant saddles on collision
cases; the singularity argument), `tnc_second_order_vandermonde_opus_S418.py` (the failed
asymptotic fit, kept as the record of why the exact route was needed). Outputs in
`05-knowledge/results/`.


---

**SIGN CORRECTED (opus-2026-07-20-S419, THM-1655).** The identity is `G(t) = +t·(log Π)′` (plus, not minus) — verified on `R=1+u`, `Π=t/(1−t)`, where `t(log Π)′ = 1+t+t²+⋯`. The residue derivation dropped a sign. **The criterion `Π = ct ⟺ TNC` is unchanged.** See THM-1655 §0.

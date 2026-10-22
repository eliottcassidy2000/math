---
id: HYP-9373
title: "Lorentzian growth equals logarithmic degree: for a C2-cofinite vertex operator algebra, let K be the maximal power of tau in its pseudo-trace functions (Jordan length of L_0 minus one). Then on the SL2(Z)-invariant space of torus amplitudes every hyperbolic gamma acts with spectral radius lambda(gamma)^K, where lambda(gamma) is the dilatation; for the cat map this is phi^(2K). Rational VOAs (K = 0): finite order, bounded by Pisano periods (THM-4762). The triplet W(p) (K = 1): phi^2 and Lucas traces (THM-4763). The predicted next case is K = 2, giving phi^4 and Sym^2 C^2 pieces, for the logarithmic extensions W(p,p') of the minimal models with rank-3 Jordan blocks."
status: >
  OPEN (opus-2026-10-22-S55).
  - PROVED for K = 0 (THM-4762 with the congruence theorems) and K = 1 for W(p) (THM-4763).
  - The upper bound lambda^K follows if the tau-degree filtration of Miyamoto's pseudo-trace functions (Duke Math. J. 122
    (2004)) has graded pieces Sym^k C^2 x R_k with R_k of finite image. That finite-image property is the open input.
  - Tensor products are consistent: W(p) x W(p') has K = 2 and eigenvalue moduli phi^(+-4), phi^(+-2), 1 (eigenvalues
    multiply).
  - The W(p,p') case is untested here. Its SL2(Z) structure (Feigin-Gainutdinov-Semikhatov-Tipunin, J. Math. Phys.
    2006) is cited as the likely source of a C^3 = Sym^2 C^2 factor (UNVERIFIED).
session: opus-2026-10-22-S55
source: 01-canon/theorems/THM-4763-the-cat-map-on-the-logarithmic-side-lucas-traces-and-spectral-radius-phi-squared.md
scripts:
  - 04-computation/experiments/catmap_modular_20261022/log_catmap.py (+ log_catmap.out)
related:
  - 01-canon/theorems/THM-4762-the-arnold-cat-map-is-the-commutator-of-t-and-s-pisano-bounds-in-modular-data.md
  - 05-knowledge/results/signature_axis_cat_map_rr_mrb_20261022.md
---

# HYP-9373 — Lorentzian growth equals logarithmic degree

## Statement

Let `V` be `C2`-cofinite. Miyamoto's pseudo-trace functions span a finite-dimensional `SL2(Z)`-invariant space of functions `Σ_(k≤K) τ^k S_k(τ)` (weight 0). Here `K` is the largest power of `τ` that occurs, which is the Jordan length of `L_0` minus one.

**Conjecture.**
* For every hyperbolic `γ ∈ SL2(Z)` with dilatation `λ(γ) > 1`, the spectral radius of `ρ_V(γ)` is exactly `λ(γ)^K`.
* For parabolic `γ` the growth is polynomial, of degree `K`.
* For elliptic `γ` the order is finite.
* In particular the cat map has spectral radius `φ^(2K)`.

**Why.** If `S_K` transforms with weight `K`, then `(γτ)^j S_K(γτ)` for `j ≤ K` spans `Sym^K C² ⊗ R_K(γ)`. The defining representation of `SL2(R) ≅ Spin(2,1)` then enters through its `K`-th symmetric power, whose top eigenvalue is `λ^K`. THM-4763 is the case `K = 1`, `Sym^1 = C²`.

## The signature axis

| side | modular image | cat map | Fibonacci |
|---|---|---|---|
| semisimple, Euclidean (`K = 0`) | finite congruence quotient | rotation of order `≤ π(N')/2` | mod `N` (Pisano) |
| logarithmic, Lorentzian (`K ≥ 1`) | contains `Sym^K` of the defining `Spin(2,1)` representation | boost of rapidity `2K log φ` | over `Z` (Lucas growth) |

**Refutation targets.**
* A `C2`-cofinite VOA whose pseudo-trace space has `τ`-degree `K` but cat-map spectral radius different from `φ^(2K)`.
* A graded piece `R_K` of infinite image.

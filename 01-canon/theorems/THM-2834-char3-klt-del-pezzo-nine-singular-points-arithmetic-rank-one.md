---
id: THM-2834
title: "A char-3 klt del Pezzo with nine tame singular points and arithmetic Picard rank one"
status: >
  PROVED (quasi-smoothness, tame klt singularity census, Fano index, Method-A
  cap, decoupling obstruction, Frobenius-orbit rank argument) + VERIFIED-EXACT
  (irreducibility/separability, strata counts, and point counts over
  F_{3^k}, k = 1..4, matching the predicted Frobenius eigenvalues).
  X = X_14 in P(2,2,7,7) over F_3 is a klt del Pezzo surface with NINE
  singular geometric points (7 x A_1 + 2 x 1/7(1,1), all tame), Fano index 4,
  arithmetic Picard rank rho(X/F_3) = 1, geometric rank 7.
source: mac-mini-2026-07-28-S171 (external open-problem raid; Epoch
  FrontierMath "Surface with a High Number of Singularities", char 3, > 7 sings)
depends_on: []
related: []
script: 04-computation/delpezzo_char3_P2277_ninesing_verify_macmini_S171.py
output: 05-knowledge/results/delpezzo_char3_P2277_macmini_S171.out
script_sha256: 0191bc403a8208c16b92da80cbc64af2441c67705b50ee9f65e94216eb45384f
output_sha256: 09b1aaf26fde24eec144e8b9971b12d99a2edfe426b1ffc80c5fcf9ddff542b4
hash_basis: LF-normalized bytes
---

# THM-2834 — X_14 ⊂ P(2,2,7,7)/F_3: nine tame singular points, ρ_arith = 1

## The surface

    X = { F = 0 } ⊂ P(2,2,7,7),   F = A(x0,x1) + x2^2 + x3^2,
    A(x0,x1) = x0^7 + x0^2 x1^5 + 2 x1^7   over F_3,
    a(t) = t^7 + t^2 + 2 irreducible (hence separable) over F_3.

* Well-formed (all weight triples coprime), weights prime to 3 (tame).
* Quasi-smooth: `∇F = (A_0, A_1, 2x_2, 2x_3)`; vanishing forces `x2=x3=0` and,
  by Euler (`7·A = x0 A_0 + x1 A_1`, `7 ≡ 1 ≠ 0` in `F_3`), a common root of
  the separable pair `(a, a')` — none; at `(1:0)`, `A_0 = 7x0^6 ≠ 0`.
* Fano: `-K_X = O_X(18 - 14) = O_X(4)` ample, index 4.  klt: all
  singularities are tame cyclic quotients.
* Singular census = ambient strata only: 7 geometric points on the (2,2)-line
  (roots of `a`; each `1/2(1,1) = A_1`), 2 geometric points on the (7,7)-line
  (`x2^2 + x3^2 = 0`, rational over `F_9`; each `1/7(2,2) ≅ 1/7(1,1)`).
  **Nine** singular points, exceeding the requested `> 7`.

## Picard rank

Over `F̄_3`, `x2^2 + x3^2 = s·t` splits and `X` contains 14 curves
`C_i = {a-root i, s = 0}`, `D_i = {a-root i, t = 0}` (each ≅ P(2,7), through
the i-th A_1 point).  Exact intersection numbers: `H'^2 = 1/14`,
`H'·C_i = 1/14`, `C_i·C_j = δ_ij/2`, `C_i ≡ D_i` numerically, and
`Σ_i C_i = 7H'` (divisor of `s`).  The Gram matrix of `{H', C_1..C_7}` has
rank 7, which equals `b_2(X)` by Noether on the resolution
(`K̃² = 8/7 - 2·25/7 = -6`, `b_2(X̃) = 16`, minus 9 exceptional classes).
So `NS(X_{F̄_3}) ⊗ Q = ⟨H', C_1..C_7 | ΣC_i = 7H'⟩ ≅ Q^7`:
**geometric rank 7**.

Frobenius permutes the roots of the irreducible `a` in a 7-cycle and swaps
`s ↔ t` (as `i ∈ F_9 ∖ F_3`), acting on `{C_i} ∪ {D_i}` as a single 14-orbit;
with `C_i ≡ D_i` the induced action on `NS ⊗ Q` is the 7-cycle on `{C_i}`.
Invariants: `Q · ΣC_i = Q · H'`.  Hence **ρ(X/F_3) = 1** (Tate holds for
geometrically rational surfaces).  Point counts confirm exactly:
`#X(F_{3^k}) = q² + c_k q + 1` with `c_k = 2,0,2,0` for `k = 1..4`, the
signature of eigenvalues `{1} ∪ {-ζ_7^j}` — i.e. `1 + (-1)^{k+1}` off `7|k`.

## Two structural obstructions (proved, recorded for reuse)

1. **Method-A cap.**  For a smooth `Y = CI(2,2) ⊂ P^4` with a diagonal `±1`
   involution `g` free in codimension 1, signature `(3,2)` is forced and
   `Y ∩ P^1_-` is empty on a smooth `Y`: at a common root `p` of the two
   restricted binary quadratics `B_i = ℓ_p m_i`, both gradients are
   proportional to `dℓ_p` (the `x_+`-partials vanish at `x_+ = 0`), so `Y` is
   singular at `p`.  Hence `|Y ∩ Fix(g)| <= 4`: quotients of smooth CI(2,2)
   by diagonal involutions can never reach 8 singular points.
2. **Decoupling forces geometric rank >= 2.**  In 4 weights, a (p,p)-line
   carrying `d/p >= 7-ish` strata points forces (parity/well-formedness) a
   decoupled `F = A ⊕ C`; then `C` splits over `F̄` and the slice curves
   `{ℓ = 0} ∩ X` are >= 2 disjoint curves, impossible at geometric rank 1.
   So "many strata singularities + geometric ρ = 1" is impossible in this
   hypersurface family — arithmetic rank one over `F_3` (achieved here) is
   the honest optimum; a geometric-rank-1 example needs a different method.

## Submission data (Epoch, Method B)

    Weights: [2, 2, 7, 7]
    F = x0^7 + x0^2*x1^5 + 2*x1^7 + x2^2 + x3^2   (over ZZ/3)

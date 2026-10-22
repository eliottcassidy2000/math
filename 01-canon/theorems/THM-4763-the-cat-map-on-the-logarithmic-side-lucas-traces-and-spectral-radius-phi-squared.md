---
id: THM-4763
title: "The cat map on the logarithmic side of Kazhdan-Lusztig: for the triplet algebra W(p), p >= 2, the SL2(Z) action on the (3p-1)-dimensional space of torus amplitudes (spanned by theta_(s,p)/eta, theta'_(s,p)/eta and tau theta'_(s,p)/eta) splits as R_(p+1) + C^2 x R_(p-1), with C^2 the defining representation (KNOWN: Feigin-Gainutdinov-Semikhatov-Tipunin 2006, eq. (1.2) and Thm 2.3). Corollary: the cat map A = [[2,1],[1,1]] acts with eigenvalues of modulus phi^-2, 1 and phi^2, and Tr rho(A^n) = Tr R_(p+1)(A^n) + L_(2n) Tr R_(p-1)(A^n), with L the Lucas numbers. The image has infinite order and spectral radius phi^2, the dilatation of the Anosov map. Every gamma in SL2(Z) has the same spectral radius as gamma itself. By the FGST equivalence the same holds for the Lyubashenko action on the centre of the restricted quantum group at q = exp(i pi/p). In every semisimple modular category the cat map has finite order bounded by a Pisano period (THM-4762)."
status: >
  PROVED (opus-2026-10-22-S55; independently audited, MISTAKE-665,
  04-computation/experiments/catmap_modular_20261022/audit_S55/REPORT.md section D).
  KNOWN: the decomposition rho = R_(p+1) + C^2 x R_(p-1), with C^2 the standard representation, is
  Feigin-Gainutdinov-Semikhatov-Tipunin, Comm. Math. Phys. 265 (2006) 47-93, hep-th/0504093, eq. (1.2) and Thm 2.3.
  FGST also identify the span with the torus conformal blocks (section 2.2) and prove the equivalence with the centre
  of the restricted quantum group (Thm 1.1, section 5). The pre-audit version cited FGST only for the identification
  and missed that they state the decomposition. It is re-derived here from the modularity of unary theta series
  (weight 1/2 for the even ones, weight 3/2 for the odd ones; Weil representation of the lattice sqrt(2p) Z) and of eta.
  NEW as stated, and one line from FGST: the cat-map corollaries (Lucas traces, spectral radius phi^2 = dilatation).
  Neither the session's nor the audit's search (2026-10-22) found them stated.
  VERIFIED for p = 2..6:
  - log_catmap.out: S and T are recomputed from the functions by least squares at 40 digits, with invariance residuals
    <= 2e-20. The even and odd blocks decouple and the eigenvalue moduli are exactly {phi^-2, 1, phi^2}. The odd block
    equals A x M(A), with M(A) read off from its (O, O) block, and the Lucas trace law holds for n = 1..12 against that
    independent M(A). The pre-audit trace test was tautological (it divided by L_(2n) first); it was replaced.
  - audit_S55/d_wp_analytic.out: an analytic construction by Poisson summation, with no fitting, gives
    F(gamma tau) = rho(gamma) F(tau) to 1e-40 and the trace law for n = 1..12.
session: opus-2026-10-22-S55
source: 05-knowledge/results/signature_axis_cat_map_rr_mrb_20261022.md
scripts:
  - 04-computation/experiments/catmap_modular_20261022/log_catmap.py (+ log_catmap.out): numerical construction, p = 2..6
  - 04-computation/experiments/catmap_modular_20261022/audit_S55/d_wp_analytic.py (+ .out): audit, analytic construction
related:
  - 01-canon/theorems/THM-4762-the-arnold-cat-map-is-the-commutator-of-t-and-s-pisano-bounds-in-modular-data.md
  - 05-knowledge/hypotheses/HYP-9373-lorentzian-growth-equals-logarithmic-degree-for-c2-cofinite-vertex-algebras.md
---

# THM-4763 — Lucas traces on the logarithmic side

## 1. The space and its decomposition

Fix `p ≥ 2` and the unary theta series
```
θ_(s,p)(τ)  = Σ_(n ∈ Z + s/2p) q^(p n²)          (s = 0..p; even; weight 1/2),
θ'_(s,p)(τ) = Σ_(n ∈ Z + s/2p) 2pn q^(p n²)      (s = 1..p−1; odd; weight 3/2),
```
with `q = exp(2πiτ)`. Then `θ_(−s) = θ_s` and `θ'_(−s) = −θ'_s`.

Put `E_s = θ_(s,p)/η` (weight 0) and `O_s = θ'_(s,p)/η` (weight 1). For `γ = [[a,b],[c,d]]` the weight-1 law `O(γτ) = (cτ+d) M(γ) O(τ)` gives
```
(γτ)·O(γτ) = (aτ + b) M(γ) O(τ),      O(γτ) = (cτ + d) M(γ) O(τ).
```
So on the pairs `(τO_s, O_s)` the group acts by the defining matrix `γ` tensored with `M(γ)`. Here `M` is the odd part of the Weil representation of `Z/2p`, twisted by the η multiplier; it has the `SU(2)_(p−2)`-type S-matrix `sin(π s s'/p)`. The even functions span the `(p+1)`-dimensional `R_(p+1)`. Hence
```
ρ ≅ R_(p+1) ⊕ (C²_def ⊗ R_(p−1)),     dim = (p+1) + 2(p−1) = 3p − 1.
```
This is the representation on the W(p) characters and pseudo-characters, and on the centre of `Ū_q sl(2)` at `q = exp(iπ/p)`. FGST state this decomposition explicitly (CMP 265 (2006), eq. (1.2) and Thm 2.3), with `R_(p+1)` the lattice-VOA characters and `R_(p−1)` the `sl(2)_(p−2)` characters. The product of two commuting representations in their abstract is a reformulation of it. The derivation above is a re-derivation, not a new result.

## 2. The cat map

`R_(p±1)` factor through congruence quotients, so `R_(p±1)(A)` have finite order. In the defining representation `A^n = Q^(2n)` has trace `L_(2n)` and eigenvalues `φ^(±2n)`. Therefore
```
Tr ρ(A^n) = Tr R_(p+1)(A^n) + L_(2n) · Tr R_(p−1)(A^n),
spec ρ(A) = spec R_(p+1)(A) ∪ φ^(±2) · spec R_(p−1)(A).
```
Examples (`log_catmap.out`):
* `p = 2`: `R_1` is the η² multiplier character (order 12; `O_1 = θ'_(1,2)/η = η²`). It is trivial on `A`, a commutator, and `R_3(A)` has order 3. The traces are 3, 7, 21, 47, 123, 325, 843, 2207, i.e. `L_(2n) + 3·[3 | n]`.
* `p = 4`: `R_3` is Ising-type, with cat-map order 3 (THM-4762, HYP-9372). The traces are 1, −1, 55, −1, 1, 971, …, with `55 = 1 + 3·L_6` and `971 = 5 + 3·L_12`.

More generally, for every `γ ∈ SL2(Z)` the spectral radius of `ρ(γ)` equals that of `γ`:
* `φ^(±2)`-type dilatations for hyperbolic `γ`;
* 1 with linear (Jordan) growth for parabolic `γ` (the Dehn twist `T`);
* finite order for elliptic `γ`.

## 3. The dichotomy (Kazhdan–Lusztig, both sides)

* **Semisimple** (rational VOAs; Kazhdan–Lusztig–Finkelberg: affine VOAs at positive integer level ≃ semisimplified tilting categories of quantum groups at roots of unity). The cat map has finite order dividing `ord_(N')(A)`, which is `π(N')/2` for `N' ≥ 3`. **Fibonacci numbers mod N.**
* **Logarithmic** (W(p) and `Ū_q sl(2)` at `q = exp(iπ/p)` are Kazhdan–Lusztig-dual, FGST). The cat map has spectral radius `φ²` and Lucas-number traces. **Fibonacci numbers over Z.**

The defining representation `C²` of `SL2(R) ≅ Spin(2,1)` is what survives on the logarithmic side. Hyperbolic mapping classes then act as Lorentz boosts, instead of as rotations of a finite quotient. A boost with spinor eigenvalues `λ(γ)^(±1)` has `SO(2,1)` rapidity `2 log λ(γ)`; for the cat map, `λ = φ²` and the rapidity is `4 log φ`.

**Not claimed.**
* No priority for the decomposition (FGST 2006).
* Nothing beyond W(p). The extension to higher logarithmic degree is HYP-9373; its `K = 2` case is also FGST's (Nucl. Phys. B 757 (2006), Thm 1.2).
* The braided equivalence for all `p` is not used. Only the SL2(Z)-equivalence on the centre is cited, and the numerical verification does not depend on it.

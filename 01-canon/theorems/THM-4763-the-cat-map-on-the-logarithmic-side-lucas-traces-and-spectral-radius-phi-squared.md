---
id: THM-4763
title: "The cat map on the logarithmic side of Kazhdan-Lusztig: for the triplet algebra W(p), p >= 2, the SL2(Z) action on the (3p-1)-dimensional space of torus amplitudes (spanned by theta_(s,p)/eta, theta'_(s,p)/eta and tau theta'_(s,p)/eta) splits as R_(p+1) + C^2 x R_(p-1), with C^2 the defining representation. Hence the cat map A = [[2,1],[1,1]] acts with eigenvalues of modulus phi^-2, 1 and phi^2, and Tr rho(A^n) = Tr R_(p+1)(A^n) + L_(2n) Tr R_(p-1)(A^n), with L the Lucas numbers. The image has infinite order and spectral radius phi^2, the dilatation of the Anosov map. Every gamma in SL2(Z) has the same spectral radius as gamma itself. By the Feigin-Gainutdinov-Semikhatov-Tipunin equivalence the same holds for the Lyubashenko action on the centre of the restricted quantum group at q = exp(i pi/p). In every semisimple modular category the cat map has finite order bounded by a Pisano period (THM-4762)."
status: >
  PROVED (opus-2026-10-22-S55). The decomposition and its consequences are elementary from the classical modularity of
  unary theta series (weight 1/2 for the even ones, weight 3/2 for the odd ones; Weil representation of the lattice
  sqrt(2p) Z) and of eta. The identification of this span with the W(p) torus amplitudes, and its SL2(Z)-equivalence
  with the centre of the restricted quantum group, is CITED: Feigin-Gainutdinov-Semikhatov-Tipunin, Comm. Math. Phys.
  265 (2006), hep-th/0504093.
  VERIFIED independently of any transformation formula for p = 2..6 (log_catmap.out):
  - S and T are recomputed from the functions by least squares at 40 digits, with invariance residuals <= 2e-20;
  - the even and odd blocks decouple, the eigenvalue moduli are exactly {phi^-2, 1, phi^2}, and the Lucas trace law
    holds for n = 1..10.
  The dynamical reading (Lucas traces, spectral radius = dilatation) was NOT FOUND in a targeted search (2026-10-22).
session: opus-2026-10-22-S55
source: 05-knowledge/results/signature_axis_cat_map_rr_mrb_20261022.md
scripts:
  - 04-computation/experiments/catmap_modular_20261022/log_catmap.py (+ log_catmap.out): numerical construction, p = 2..6
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
This is the representation on the W(p) characters and pseudo-characters, and on the centre of `Ū_q sl(2)` at `q = exp(iπ/p)` (FGST, CITED). Their statement that the ribbon element's Jordan decomposition "splits this representation into a pointwise product of two commuting SL(2,Z)-representations" is this tensor structure.

## 2. The cat map

`R_(p±1)` factor through congruence quotients, so `R_(p±1)(A)` have finite order. In the defining representation `A^n = Q^(2n)` has trace `L_(2n)` and eigenvalues `φ^(±2n)`. Therefore
```
Tr ρ(A^n) = Tr R_(p+1)(A^n) + L_(2n) · Tr R_(p−1)(A^n),
spec ρ(A) = spec R_(p+1)(A) ∪ φ^(±2) · spec R_(p−1)(A).
```
Examples (`log_catmap.out`):
* `p = 2`: `R_(p−1)` is trivial and `R_3(A)` has order 3. The traces are 3, 7, 21, 47, 123, 325, 843, 2207, i.e. `L_(2n) + 3·[3 | n]`.
* `p = 4`: `R_3` is Ising-type, with cat-map order 3 (THM-4762, HYP-9372). The traces are 1, −1, 55, −1, 1, 971, …, with `55 = 1 + 3·L_6` and `971 = 5 + 3·L_12`.

More generally, for every `γ ∈ SL2(Z)` the spectral radius of `ρ(γ)` equals that of `γ`:
* `φ^(±2)`-type dilatations for hyperbolic `γ`;
* 1 with linear (Jordan) growth for parabolic `γ` (the Dehn twist `T`);
* finite order for elliptic `γ`.

## 3. The dichotomy (Kazhdan–Lusztig, both sides)

* **Semisimple** (rational VOAs; Kazhdan–Lusztig–Finkelberg: affine VOAs at positive integer level ≃ quantum groups at roots of unity). The cat map has finite order dividing `π(N')/2`. **Fibonacci numbers mod N.**
* **Logarithmic** (W(p) ≃ `Ū_q sl(2)` at `q = exp(iπ/p)`, FGST). The cat map has spectral radius `φ²` and Lucas-number traces. **Fibonacci numbers over Z.**

The defining representation `C²` of `SL2(R) ≅ Spin(2,1)` is what survives on the logarithmic side. Hyperbolic mapping classes then act as Lorentz boosts of rapidity `log λ(γ)`, instead of as rotations of a finite quotient.

**Not claimed.**
* Nothing beyond W(p). The extension to higher logarithmic degree is HYP-9373.
* The braided equivalence for all `p` is not used. Only the SL2(Z)-equivalence on the centre is cited, and the numerical verification does not depend on it.

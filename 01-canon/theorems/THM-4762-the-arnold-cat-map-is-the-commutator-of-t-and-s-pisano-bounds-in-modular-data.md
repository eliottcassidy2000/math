---
id: THM-4762
title: "The Arnold cat map is the commutator [T,S], and Pisano periods bound it in modular data. (1) In SL2(Z), A = [[2,1],[1,1]] = t s t^-1 s^-1 = Q^2 with Q the Fibonacci matrix, so ord(A mod N) = pi(N)/2 for N >= 3 (pi = Pisano period). pi(N) is even for N >= 3, so no modulus has pi(N) = 15; pi(4) = 6, pi(6) = 24, pi(30) = 120, and 30 is the least modulus with period 120. (2) For the modular data (S, T) of any modular tensor category, C = T S T^-1 S^-1 does not depend on the scalar normalisations of S and T (central charge, Gauss sum phase). For pointed categories of even lattices it depends only on the discriminant form, so definite and indefinite (Lorentzian) lattices with isomorphic discriminant forms give the same C. (3) With the congruence theorems, C^(ord_N A) is scalar (N = order of the normalised twists), and ord C divides ord_N'(A) (N' = order of the linear T); ord_N(A) = pi(N)/2 for N >= 3. (4) Census of 135 modular data: the projective order of C equals ord_N(A) in 88 cases (83 with N >= 3; all D(Z_N) with 3 <= N <= 30, where it is a theorem, and all SU(3)_k with k <= 10), and otherwise falls short by a factor 2 or 4."
status: >
  PROVED (opus-2026-10-22-S55; independently audited, MISTAKE-665,
  04-computation/experiments/catmap_modular_20261022/audit_S55/REPORT.md sections A-C): (1) and (2) are elementary.
  (3) is PROVED given the CITED congruence-subgroup theorems: Ng-Schauenburg (Comm. Math. Phys. 300 (2010), Thm 6.8)
  for the projective representation, with level N computed from twists normalised so that theta_0 = 1; Dong-Lin-Ng
  (Algebra Number Theory 9 (2015) 2121-2166, Thm II) for the linear one. (4) is FINITE: double-precision modular data;
  orders are read from the eigenvalues and confirmed by C^m = I to 1e-6 (catmap_mtc.out, catmap_census.json). The
  audit rebuilt all 135 rows with its own builders and found agreement in every field. It also checked S symmetric and
  unitary, S^2 a permutation, (S T)^3 = e^(2 pi i c/8) S^2 and Verlinde integrality (rank <= 100) on every datum
  (audit_S55/c_modular_data.out); catmap_mtc.py now runs the same checks.
  PRIOR ART: the Pisano period as the order of the Fibonacci matrix is classical (Wall 1960). The order of the cat map
  mod N and its quantum (Weil-representation) version are KNOWN: Hannay-Berry 1980; Dyson-Falk, Amer. Math. Monthly 99
  (1992) 603-614. Chandra (arXiv:2607.24857) defines pi(N) as the order of A in SL(2, Z/NZ). Torus bundles with
  hyperbolic monodromy in Chern-Simons TQFT: Jeffrey, Comm. Math. Phys. 147 (1992). Mapping classes as Fourier integral
  operators: Charles, arXiv:1005.3447. The D(Z_N) case is a theorem: in the flux basis (a partial Fourier transform in
  the charge label) rho is the permutation representation of SL2(Z/N) on (Z/N)^2, so C is the discrete cat map of
  Dyson-Falk and ord C = projective ord C = ord_N(A) (audit_S55/x_extra.out, N = 3, 4, 5, 6, 8). The Pisano bound in (3)
  is a one-line corollary of the congruence theorems. NOT FOUND in the session's or the audit's search (2026-10-22): the
  reading of the cat map as the canonical commutator [T,S] of arbitrary modular data, the normalisation and signature
  independence, and the census.
session: opus-2026-10-22-S55
source: 05-knowledge/results/signature_axis_cat_map_rr_mrb_20261022.md
scripts:
  - 04-computation/experiments/catmap_modular_20261022/catmap_mtc.py (+ catmap_mtc.out, catmap_census.json): Pisano facts, modular data builders, census
  - 04-computation/experiments/catmap_modular_20261022/center_check.py (+ center_check.out): central powers (HYP-9372)
related:
  - 01-canon/theorems/THM-4763-the-cat-map-on-the-logarithmic-side-lucas-traces-and-spectral-radius-phi-squared.md
  - 05-knowledge/hypotheses/HYP-9372-pisano-saturation-of-the-cat-map-and-the-ising-spin-structure-exception.md
  - 01-canon/theorems/THM-4656-wall-sun-sun-golden-orbit-explosion.md
  - 01-canon/theorems/THM-4708-rational-cm-values-of-the-argyres-douglas-character-hauptmoduln-at-levels-5-and-7.md
---

# THM-4762 — the cat map as the commutator [T, S]

## 1. Statement

Write `s = [[0,−1],[1,0]]`, `t = [[1,1],[0,1]]`, `Q = [[1,1],[1,0]]` and `A = [[2,1],[1,1]]`.

**(1) Three faces of one matrix.**
```
A = t s t^-1 s^-1 = Q^2,     Q^n = [[F_(n+1), F_n], [F_n, F_(n−1)]],     A^n = Q^(2n).
```
* `A` is the Arnold cat map: hyperbolic, with trace 3 and eigenvalues `φ^(±2)`. In `SL2(R) ≅ Spin(2,1)` it is a boost. Its spinor eigenvalues are `φ^(±2)`, so its image in `SO(2,1)` has eigenvalues `φ^(±4), 1` and rapidity `4 log φ`.
* The generators `s` and `st` are elliptic, of orders 4 and 6; they are rotations.
* The Pisano period `π(N)` is the order of `Q` mod `N`. For `N ≥ 3` it is even, so
```
ord(A mod N) = π(N)/2   (N ≥ 3),   ord(A mod 2) = 3 = π(2).
```
* Hence `π(N) = 15` never occurs. The only odd periods are `π(1) = 1` and `π(2) = 3`.
* `π(4) = 6`, `π(6) = 24`, `π(30) = 120`, and no `N < 30` has `π(N) = 120` (checked for `N < 1000`).
* In modular-data language: the cat map has order 3, 12 and 60 on the torus ground states of the `Z_4`, `Z_6` and `Z_30` toric codes. It also has order 60 on Fibonacci ⊠ `D(Z_2)` ⊠ `D(Z_3)`, whose conductor is 30 (§3).

**(2) A normalisation-free invariant.** Let `(S, T)` be modular data, with `T` the diagonal matrix of twists. The projective representation of `SL2(Z)` sends `s ↦ S`, `t ↦ T`, and every lift rescales `S` and `T` by scalars. Scalars cancel in a commutator, so
```
C := T S T^-1 S^-1 = ρ(A)
```
is a well-defined unitary matrix. It is the same for every lift: any central charge `c` in `exp(−2πic/24)`, and any sign or phase convention for `S`.

For a pointed category built from an even lattice `L`:
* Take `S_xy = |G|^(−1/2) exp(−2πi b(x,y))` and `T_x = exp(2πi q(x))` on the discriminant form `(G, q) = (L*/L, x² /2 mod 1)`.
* The signature of `L` enters only through the Gauss–Milgram phase `exp(2πiσ/8)`, a scalar normalisation that `C` does not see.
* So **definite and indefinite lattices with isomorphic discriminant forms have identical `C`**. Examples (all computed, difference 0): `A_1` and `A_1 ⊕ U`; `A_2` and `A_2 ⊕ U`; `E_8` and `U = II_(1,1)` (both trivial).
* `A_1(−1)` is not an example. Its discriminant form `(Z/2, 3/4)` is the conjugate of `(Z/2, 1/4)`, so its `C` is the complex conjugate of `C(A_1)`, with the same orders. (The pre-audit text listed it as identical.)

**(3) Pisano bounds.** Let `N` be the order of the twists `T̄ = diag(θ_j)` and `N'` the order of the linear `T = T̄·exp(−2πic/24)`.
* The projective representation factors through `SL2(Z/N)` (Ng–Schauenburg). Hence `C^(ord_N A)` is scalar.
* The linear representation `s ↦ S`, `t ↦ T` factors through `SL2(Z/N')` (Dong–Lin–Ng, Thm II). It is a genuine representation because `(ST̄)³ = e^(2πic/8)·S²`. Hence `ord C` divides `ord_(N')(A)`.
* `ord_N(A) = π(N)/2` for `N ≥ 3`. For `N = 1, 2` it is 1 and 3, where `π(N)/2` is not an integer.

**(4) Census.** Covered (all in `catmap_mtc.out`):
* Fibonacci;
* the 23 Virasoro minimal models `M(p, p')` with `p < p' ≤ 12` and `pp' ≤ 66` (there are 46 minimal models with `pp' ≤ 66`);
* `SU(2)_k` for `k ≤ 60` and `SU(3)_k` for `k ≤ 10`;
* pointed categories of ten lattices, definite and indefinite;
* `D(Z_N)` for `N ≤ 30`, `D(S_3)`, and one Deligne product.

All 135 obey (3). The projective order of `C` equals `ord_N(A)` in 88 of them. Five of these have `N ≤ 2`: `M(2,3)`, `E_8` and `U` (`N = 1`), `D_4` and `D(Z_2)` (`N = 2`). So 83 saturate `π(N)/2` with `N ≥ 3`. In the rest the projective order falls short by a factor 2 (43 cases) or 4 (4 cases: `M(3,4)`, `SU(2)_2`, `M(4,5)`, `M(4,11)`). HYP-9372 gives the mechanism and the Ising exception.

## 2. Proofs

(1) `t s = [[1,−1],[1,0]]` and `t^(−1) s^(−1) = [[1,1],[−1,0]]`, so `t s t^(−1) s^(−1) = [[2,1],[1,1]] = Q²`. Then `ord_N(Q²) = π(N)/gcd(2, π(N))`. That `π(N)` is even for `N ≥ 3` is classical (Wall 1960): `Q^m = I` forces `(−1)^m = det Q^m = 1` mod `N`.

(2) Replace `S, T` by `αS, βT`. Then `(βT)(αS)(βT)^(−1)(αS)^(−1) = T S T^(−1) S^(−1)`. For the pointed case, `S` and `T̄` are built from `(G, q)` alone, and the lattice enters only through `c = σ`.

(3) If `ρ̄` factors through `SL2(Z/N)`, then `ρ̄(A)^m` is scalar whenever `A^m ≡ I` mod `N`, in particular for `m = ord_N(A)`. Congruence to `−I` would not do: `ρ(−I) = S²` is the charge conjugation, which need not be scalar. Since `C` is an honest matrix equal to `ρ(A)` for the linear lift, it inherits `ord C | ord_(SL2(Z/N'))(A)`.

(D(Z_N)) The torus ground states of the untwisted `Z_N` gauge theory are indexed by flat connections, i.e. flux pairs in `(Z/N)²`, and the mapping-class group acts on them by permutation. The anyon basis is a partial Fourier transform of the flux basis. So `C` is conjugate to the permutation matrix of the discrete cat map on `(Z/N)²`. A non-identity permutation matrix is never scalar, which gives `ord C = projective ord C = ord_N(A)`. This was checked numerically for `N = 3, 4, 5, 6, 8` (audit).

## 3. Census highlights (`catmap_census.json`)

| modular data | `N` | `π(N)/2` | projective ord `C` | ord `C` |
|---|---|---|---|---|
| Fibonacci = `(G_2)_1`, Lee–Yang `M(2,5)` | 5 | 10 | 5 | 10 |
| Ising `M(3,4)`, `SU(2)_2` | 16 | 12 | 3 | 3 |
| `M(2,7)` | 7 | 8 | 4 | 4 |
| `M(2,9)` | 9 | 12 | 6 | 12 |
| tricritical Ising `M(4,5)` | 80 | 60 | 15 | 30 |
| `M(5,6)` | 120 | 60 | 60 | 60 |
| `SU(3)_7` | 30 | 60 | 60 | 60 |
| `D(Z_N)`, `N = 3..30` | `N` | `π(N)/2` | `π(N)/2` | `π(N)/2` |
| `D(S_3)` | 6 | 12 | 6 | 6 |
| golden Lorentzian lattice `[[2,1],[1,−2]]` (`L*/L = Z_5`) | 5 | 10 | 10 | 10 |

**The golden Lorentzian lattice.** The indefinite lattice `(Z[φ], 2N(x))`, with Gram matrix `[[2,1],[1,−2]]` and determinant `−5`, has the cat map among its automorphs. Multiplication by `φ²` preserves the norm form. Its discriminant form `Z_5` gives a pointed category whose `C` has order `10 = π(5)/2`, the same order as the non-abelian Fibonacci category.

**Not claimed.**
* No statement for non-semisimple categories. There `T` has Jordan blocks and `C` can have infinite order (THM-4763).
* No claim that the census orders equal `π(N)/2` in general.

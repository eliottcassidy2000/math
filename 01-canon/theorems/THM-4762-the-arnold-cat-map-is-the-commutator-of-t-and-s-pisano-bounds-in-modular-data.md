---
id: THM-4762
title: "The Arnold cat map is the commutator [T,S], and Pisano periods bound it in modular data. (1) In SL2(Z), A = [[2,1],[1,1]] = t s t^-1 s^-1 = Q^2 with Q the Fibonacci matrix, so ord(A mod N) = pi(N)/2 for N >= 3 (pi = Pisano period). pi(N) is even for N >= 3, so no modulus has pi(N) = 15; pi(4) = 6, pi(6) = 24, pi(30) = 120, and 30 is the least modulus with period 120. (2) For the modular data (S, T) of any modular tensor category, C = T S T^-1 S^-1 does not depend on the scalar normalisations of S and T (central charge, Gauss sum phase). For pointed categories of even lattices it depends only on the discriminant form, so definite and indefinite (Lorentzian) lattices with isomorphic discriminant forms give the same C. (3) With the congruence theorems, C^(pi(N)/2) is scalar (N = order of the normalised twists), and ord C divides pi(N')/2 (N' = order of the linear T). (4) Census of 135 modular data: the projective order of C equals pi(N)/2 in 88 cases (all D(Z_N) with 3 <= N <= 30 and all SU(3)_k with k <= 10), and otherwise falls short by a factor 2 or 4."
status: >
  PROVED (opus-2026-10-22-S55): (1) and (2) are elementary. (3) is PROVED given the CITED congruence-subgroup theorems:
  Ng-Schauenburg (Comm. Math. Phys. 300, 2010) for the projective representation; Dong-Lin-Ng (Adv. Math. 2015) for
  the linear one. (4) is FINITE: double-precision modular data; orders are read from the eigenvalues and confirmed by
  C^m = I to 1e-6 (catmap_mtc.out, catmap_census.json).
  PRIOR ART: the Pisano period as the order of the Fibonacci matrix is classical (Wall 1960). The order of the cat map
  mod N and its quantum (Weil-representation) version are KNOWN: Hannay-Berry 1980; Dyson-Falk, Amer. Math. Monthly 99
  (1992); a 2026 arXiv paper on discrete cat maps (arXiv:2607.24857) defines pi(N) as this order. The equality
  ord C = pi(N)/2 for the toric codes D(Z_N) is that classical quantum cat map in modular-data language. NOT FOUND in a
  targeted search (2026-10-22): the reading of the cat map as the canonical commutator [T,S] of arbitrary modular data,
  the normalisation and signature independence, and the census.
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
* `A` is the Arnold cat map: hyperbolic, with trace 3 and eigenvalues `φ^(±2)`. In `SL2(R) ≅ Spin(2,1)` it is a boost of rapidity `2 log φ`.
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
* So **definite and indefinite lattices with isomorphic discriminant forms have identical `C`**. Examples (all computed): `A_1`, `A_1 ⊕ U`, `A_1(−1)`; `A_2`, `A_2 ⊕ U`; `E_8` and `U = II_(1,1)` (both trivial).

**(3) Pisano bounds.** Let `N` be the order of the twists `T̄ = diag(θ_j)` and `N'` the order of the linear `T = T̄·exp(−2πic/24)`.
* The projective representation factors through `SL2(Z/N)` (Ng–Schauenburg). Hence `C^(π(N)/2)` is scalar.
* The linear representation factors through `SL2(Z/N')` (Dong–Lin–Ng). Hence `ord C` divides `π(N')/2`.

**(4) Census.** Covered (all in `catmap_mtc.out`):
* Fibonacci;
* the 23 Virasoro minimal models `M(p, p')` with `pp' ≤ 66`;
* `SU(2)_k` for `k ≤ 60` and `SU(3)_k` for `k ≤ 10`;
* pointed categories of ten lattices, definite and indefinite;
* `D(Z_N)` for `N ≤ 30`, `D(S_3)`, and one Deligne product.

All 135 obey (3). The projective order of `C` equals `π(N)/2` in 88 of them. In the rest it is `π(N)/4` or `π(N)/8`; HYP-9372 gives the mechanism and the Ising exception.

## 2. Proofs

(1) `t s = [[1,−1],[1,0]]` and `t^(−1) s^(−1) = [[1,1],[−1,0]]`, so `t s t^(−1) s^(−1) = [[2,1],[1,1]] = Q²`. Then `ord_N(Q²) = π(N)/gcd(2, π(N))`. That `π(N)` is even for `N ≥ 3` is classical (Wall 1960): `Q^m = I` forces `(−1)^m = det Q^m = 1` mod `N`.

(2) Replace `S, T` by `αS, βT`. Then `(βT)(αS)(βT)^(−1)(αS)^(−1) = T S T^(−1) S^(−1)`. For the pointed case, `S` and `T̄` are built from `(G, q)` alone, and the lattice enters only through `c = σ`.

(3) If `ρ̄` factors through `SL2(Z/N)`, then `ρ̄(A)^m` is scalar whenever `A^m ≡ I` mod `N`, in particular for `m = π(N)/2`. Since `C` is an honest matrix equal to `ρ(A)` for the linear lift, it inherits `ord C | ord_(SL2(Z/N'))(A)`.

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

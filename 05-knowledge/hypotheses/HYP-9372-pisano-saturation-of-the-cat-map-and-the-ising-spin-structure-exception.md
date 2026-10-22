---
id: HYP-9372
title: "Pisano saturation of the cat map in modular data, and the Ising spin-structure exception. For a modular tensor category with projective conductor N >= 3, let C = T S T^-1 S^-1 (THM-4762). (a) The deficiency r = (pi(N)/2)/(projective order of C) is always 1, 2 or 4. (b) Faithfulness modulo the centre: C^j is scalar exactly when A^j is a central matrix lambda I in SL2(Z/N), unless the category has an Ising-type (SU(2)_2-Galois-orbit) Deligne factor whose 4-part is not covered by the other factors. (c) In the exceptions C acts through the spin-structure quotient SL2(Z/2) = S_3, so the Ising factor contributes order 3 = pi(2) instead of pi(16)/2 = 12."
status: >
  OPEN (opus-2026-10-22-S55). FINITE evidence (double-precision modular data):
  - (a) holds in all 135 modular data of the THM-4762 census;
  - (b) holds in 119 of the 123 with N >= 2 (center_check.out). The four exceptions are M(3,4), SU(2)_2, M(4,5) and
    M(4,11), each containing the Ising/SU(2)_2 factor; M(4,7) and M(4,9) also contain it, but there the other factor's
    2-part masks the loss.
  - (c) is a proposed mechanism, consistent with ord C = 3 for Ising and 15 = lcm(3, 5) for M(4,5) and M(4,11).
session: opus-2026-10-22-S55
source: 01-canon/theorems/THM-4762-the-arnold-cat-map-is-the-commutator-of-t-and-s-pisano-bounds-in-modular-data.md
scripts:
  - 04-computation/experiments/catmap_modular_20261022/center_check.py (+ center_check.out)
  - 04-computation/experiments/catmap_modular_20261022/catmap_mtc.py (+ catmap_census.json)
related:
  - 01-canon/theorems/THM-4763-the-cat-map-on-the-logarithmic-side-lucas-traces-and-spectral-radius-phi-squared.md
  - 05-knowledge/results/signature_axis_cat_map_rr_mrb_20261022.md
---

# HYP-9372 — when does the cat map see the whole Pisano period?

## Background

Let `N ≥ 3` and `m_N = π(N)/2`, the order of `A` in `SL2(Z/N)` (THM-4762).
* The centre of `SL2(Z/N)` consists of the `λI` with `λ² ≡ 1`.
* `Q^k` is central exactly when `F_k ≡ 0` mod `N`, that is when `α(N) | k`, where `α` is the rank of apparition. Then `Q^(α) ≡ F_(α−1) I`.
* The classical theorem `π(N)/α(N) ∈ {1, 2, 4}` (Vinson 1963; Wall 1960) shows `|⟨A⟩ ∩ Z| ∈ {1, 2, 4}`.

So if `ρ̄` is faithful modulo the centre on `⟨A⟩`, the deficiency `r` is the number of central powers of `A` acting by scalars. This explains (a) in the faithful cases.

**Two sources of scalar central powers.**
* `−I` acts as `S²`, the charge conjugation. It is a scalar exactly when every object is self-dual. This is why `SU(3)_k` (with its non-self-dual objects) and `D(Z_N)` with `N ≥ 3` saturate, while the self-dual `SU(2)_k` often lose a factor 2.
* For `λ ≠ ±1`, e.g. `λ = 9` mod 16, `λI` acts in Weil-type representations through the scaling `x ↦ λx`.

## Evidence (center_check.out)

* 119 of 123 modular data are faithful modulo the centre on `⟨A⟩`: every scalar power of `C` is a central power of `A`.
* The exceptions are Ising `M(3,4)` and `SU(2)_2` (`N = 16`, `m_N = 12`, first central power `A^6 = 9I`, yet `C^3` is already scalar), and `M(4,5)`, `M(4,11)`. The latter two have projective order 15 against the faithful value 30.
* Masking: `M(4,7)` and `M(4,9)` contain the same Ising-type factor, but `ord_7 A = 8` and `ord_9 A = 12` already supply the 2-part. In those cases `lcm(3, x) = lcm(12, x)`, so the loss is invisible.

## Proposed mechanism for (c)

The Ising category is the bosonisation of a Majorana fermion. Its characters are built from `√(θ_i/η)`, `i = 2, 3, 4`, the three even spin-structure partition functions.
* `SL2(Z)` permutes the spin structures through `SL2(Z/2) ≅ S_3`, up to phases.
* The commutator `C` kills those phases and keeps the permutation.
* The cat map mod 2 has order `3 = π(2)`, so it cyclically permutes `(NS,NS) → (NS,R) → (R,NS)`.

So the Ising factor's cat map has the spin-structure period 3, not the conductor period 12. This connects to the frustrated Ising rings of the supplied MRB note: periodic and antiperiodic boundary conditions are two of these spin structures. See the results note.

**Refutation targets.**
* A modular data with `r ∉ {1, 2, 4}`.
* A non-faithful case with no Ising-type factor.
* An Ising-type case in which `C` does not reduce to the `S_3` permutation of spin structures.

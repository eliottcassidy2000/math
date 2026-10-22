---
id: HYP-9372
title: "Pisano saturation of the cat map in modular data, and the Ising spin-structure exception. For a modular tensor category with projective conductor N >= 3, let C = T S T^-1 S^-1 (THM-4762). (a) The deficiency r = (pi(N)/2)/(projective order of C) is always 1, 2 or 4, and it is 1 or 2 when (b) holds. (b) Faithfulness modulo the centre: if C^j is scalar then A^j is a central matrix lambda I in SL2(Z/N), unless the category has an Ising-type (SU(2)_2-Galois-orbit) Deligne factor whose 4-part is not covered by the other factors. The converse is false: a central power of A need not act by a scalar (-I acts as the charge conjugation S^2). (c) In the exceptions C acts through the spin-structure quotient SL2(Z/2) = S_3, so the Ising factor contributes order 3 = pi(2) instead of pi(16)/2 = 12."
status: >
  OPEN (opus-2026-10-22-S55; independently audited, MISTAKE-665,
  04-computation/experiments/catmap_modular_20261022/audit_S55/REPORT.md section E). FINITE evidence
  (double-precision modular data), plus two proved pieces:
  - (a) holds in all 135 modular data of the THM-4762 census (r = 1: 88, r = 2: 43, r = 4: 4). PROVED: the group <A>
    contains 1 or 2 central elements (never 4), so faithful cases have r in {1, 2} and r = 4 forces non-faithfulness.
  - (b) holds in 128 of the 132 with N >= 2. The pre-audit center_check.out had 119 of 123, because it skipped the 8
    lattices and the Deligne product; with them added it agrees with the audit's e_center.out. The four exceptions are
    M(3,4), SU(2)_2, M(4,5) and M(4,11), each containing an Ising-type factor. M(4,7) and M(4,9) also contain one, but
    stay faithful.
  - (b) was stated pre-audit as "C^j is scalar exactly when A^j is central". The converse is REFUTED: it fails in 67 of
    the 132 (center_check.out). Examples: D(Z_3), SU(3)_1, A_2 and the golden lattice, where -I acts as a non-scalar
    charge conjugation; the self-dual SU(2)_4, where A^6 = 17I mod 24 but C^6 is not scalar.
  - (c) is PROVED for the Ising-type factors themselves: in the spin-structure basis C is monomial over the 3-cycle A mod 2
    with phase product 1, so C^3 = I. The factorisation M(4,p') = Ising_nu x X with X = {(r,1) : r odd} is verified for
    p' = 5, 7, 9, 11 (nu = 7, 13, 3, 9; projective order lcm(3, po X) = 15, 12, 6, 15, as in the census).
session: opus-2026-10-22-S55
source: 01-canon/theorems/THM-4762-the-arnold-cat-map-is-the-commutator-of-t-and-s-pisano-bounds-in-modular-data.md
scripts:
  - 04-computation/experiments/catmap_modular_20261022/center_check.py (+ center_check.out)
  - 04-computation/experiments/catmap_modular_20261022/catmap_mtc.py (+ catmap_census.json)
  - 04-computation/experiments/catmap_modular_20261022/audit_S55/e_center.py (+ .out): audit, both directions and the Ising basis
related:
  - 01-canon/theorems/THM-4763-the-cat-map-on-the-logarithmic-side-lucas-traces-and-spectral-radius-phi-squared.md
  - 05-knowledge/results/signature_axis_cat_map_rr_mrb_20261022.md
---

# HYP-9372 — when does the cat map see the whole Pisano period?

## Background

Let `N ≥ 3` and `m_N = π(N)/2`, the order of `A` in `SL2(Z/N)` (THM-4762).
* The centre of `SL2(Z/N)` consists of the `λI` with `λ² ≡ 1`.
* `Q^k` is central exactly when `F_k ≡ 0` mod `N`, that is when `α(N) | k`, where `α` is the rank of apparition. Then `Q^(α) ≡ F_(α−1) I`.
* The classical theorem `π(N)/α(N) ∈ {1, 2, 4}` (Vinson 1963; Wall 1960) gives `|⟨A⟩ ∩ Z| ∈ {1, 2}`. Indeed `Q^α = λI` with `λ² = det Q^α = (−1)^α`.
  * If `α` is odd, then `λ` has order 4 and `π = 4α`. The only nontrivial central element of `⟨A⟩ = ⟨Q²⟩` is `Q^(2α) = −I`.
  * If `α` is even, then `λ² = 1`, and the central elements of `⟨A⟩` are the powers of `A^(α/2) = λI`.
  * So `π/α = 4` exactly when `α` is odd (no exception for `3 ≤ N ≤ 2000`, audit).

So if `ρ̄` is faithful modulo the centre on `⟨A⟩`, the deficiency `r` is the number of central powers of `A` acting by scalars. Faithful cases therefore have `r ∈ {1, 2}`, and `r = 4` forces non-faithfulness. (The pre-audit text had `{1, 2, 4}` here, conflating `⟨Q⟩` with `⟨A⟩`.)

**Two sources of scalar central powers.**
* `−I` acts as `S²`, the charge conjugation. It is a scalar exactly when every object is self-dual. This is why `SU(3)_k` (with its non-self-dual objects) and `D(Z_N)` with `N ≥ 3` saturate, while the self-dual `SU(2)_k` often lose a factor 2.
* For `λ ≠ ±1`, e.g. `λ = 9` mod 16, `λI` acts in Weil-type representations through the scaling `x ↦ λx`.

## Evidence (center_check.out)

* 128 of 132 modular data with `N ≥ 2` are faithful modulo the centre on `⟨A⟩`: every scalar power of `C` is a central power of `A`.
* In a faithful case, `r = 2` exactly when `⟨A⟩` has a nontrivial central element and that element acts by a scalar. Otherwise `r = 1`.
* The exceptions are Ising `M(3,4)` and `SU(2)_2` (`N = 16`, `m_N = 12`, first central power `A^6 = 9I`, yet `C^3` is already scalar), and `M(4,5)`, `M(4,11)`. The latter two have projective order 15 against the faithful value 30.
* Masking: `M(4,7)` and `M(4,9)` contain the same Ising-type factor, `M(4,p') ≅ Ising_ν ⊠ X`. There `po(X) = 4` and `6`, so the projective orders are `lcm(3, 4) = 12` and `lcm(3, 6) = 6`. The first scalar powers `C^12` and `C^6` fall on the central powers `A^12 = 97I` (mod 112) and `A^6 = 89I` (mod 144). Both cases stay faithful, with `r = 2`. (The pre-audit text used `ord_7 A = 8` and `ord_9 A = 12` here, which are the wrong numbers.)
* The converse direction fails in 67 of the 132, as noted in the status: `A^j` central does not make `C^j` scalar. Non-self-dual objects make `ρ(−I) = S²` non-scalar, and in `SU(2)_4` the central `A^6 = 17I` (mod 24) acts by a non-scalar.

## The mechanism for (c) (proved for the Ising-type factors)

The Ising category is the bosonisation of a Majorana fermion. Its characters are built from `√(θ_i/η)`, `i = 2, 3, 4`, the three even spin-structure partition functions.
* In the basis `(χ_0 + χ_ψ, χ_0 − χ_ψ, √2·χ_σ)`, proportional to `√(θ_3/η), √(θ_4/η), √(θ_2/η)`, `S` swaps the last two vectors and `T` is monomial over the first two. These are `s` and `t` mod 2 acting on the three even spin structures: `SL2(Z)` permutes them through `SL2(Z/2) ≅ S_3`, up to phases.
* The commutator `C` is monomial over the 3-cycle `A` mod 2. Its phases are `e^(−iπ/8), 1, e^(iπ/8)` for Ising and `e^(∓3iπ/8), 1` for `SU(2)_2` (audit). Their product is 1, so `C³ = I`. The phases do not cancel individually; only their product around the cycle is 1.
* The cat map mod 2 has order `3 = π(2)`, so it cyclically permutes `(NS,NS) → (NS,R) → (R,NS)`.

So the Ising factor's cat map has the spin-structure period 3, not the conductor period 12. This connects to the frustrated Ising rings of the supplied MRB note: periodic and antiperiodic boundary conditions are two of these spin structures. See the results note.

**Refutation targets.**
* A modular data with `r ∉ {1, 2, 4}`.
* A non-faithful case with no Ising-type factor.
* An Ising-type case in which `C` does not reduce to the `S_3` permutation of spin structures.

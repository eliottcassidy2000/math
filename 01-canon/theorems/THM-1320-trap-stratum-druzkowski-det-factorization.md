---
id: THM-1320
title: The TRAP-stratum theorem — the determinant of every weight-(1,−1,−k) equivariant map factorizes at the origin as det = −E₀(0)·A(0)·C(0) (polynomial identity, verified k = 2,3,4), so Keller forces all three row-units to have NONZERO constant terms; homogeneous (Drużkowski-shaped) cubes v(0) = 0 admit NO Keller map in the class; equivariant cubic-linear maps are triangular-over-weight-0-core over EVERY nontrivial weighting (no nonzero-weight cycles in 6858 grid cases; the witness's weights admit no equivariant cubic term at all); the planar weight-0 core is exactly shears + triangulars, all polynomially invertible — the classical stratum's safety IS rootlessness, and the witness's det = −2 DECODES as −(E₀(0)=2)(A(0)=1)(C(0)=1).
status: >
  PROVED/VERIFIED-EXACT — (P3) the factorization c₀(0) = −E₀(0)A(0)C(0) as a
  polynomial identity in all coefficients (engine, k = 2, 3, 4; the t=0
  evaluation of the THM-1305 c₀ formula); its corollary (unit constant terms
  necessary; homogeneous cubes ⟹ det ≡ 0). (P1) the weight arithmetic
  (w_i = 3w_{j(i)}; cycles force weight 0 — exact check over |w| ≤ 9, 6858
  nontrivial vectors, zero violations; the (1,−1,−2) table is EMPTY).
  (P2) planar cubic-linear Keller = shear/triangular families (exact
  8-equation coefficient system; the shear identity ℓ∘F = ℓ; polynomial
  inverses). SCOPE — the equivariant slice; non-equivariant Drużkowski maps
  are untouched (the classical rank theorems live there).
source: death-star-2026-07-19-S59o (HYP-8110; owner directive: run backlog xliii)
depends_on:
  - THM-1305  # the normal form and the c0 formula
  - THM-1300  # the witness
related:
  - the S59n reflection §4 (the TRAP hypothesis this confirms and sharpens)
scripts:
  - 04-computation/jacobian_trap_stratum_druzkowski_deathstar_S59o.py -> 05-knowledge/results/jacobian_trap_stratum_druzkowski_deathstar_S59o.out
---

# THM-1320 — the TRAP-stratum theorem

> **HONESTY AMENDMENT (S59p, THM-1325 §1, MISTAKE-197):** the P3 factorization
> c₀(0) = −E₀(0)A(0)C(0) is exactly det JF(0) (det const = det of the linear
> part; JF(0) is the antidiagonal of the unit constants). True and verified,
> but CLASSICALLY TRIVIAL — "unit constants nonzero" = "linear part
> invertible." The surviving content of this file: the row-by-row FRAME
> (units ↔ det factors), P1 (transversality), P2 (the planar core).

## 1. The determinant factorizes at the origin (P3, the heart)

Evaluating THM-1305's c₀ at t = 0 (where every t-shifted term dies and
M0|₀ = −A(0)C(0), K₀|₀ = E₀(0)):

  **det J F = c₀(0) = −E₀(0) · A(0) · C(0)**

— verified as a polynomial identity in ALL coefficients at k = 2, 3, 4.
The determinant of an equivariant Keller map is the product of its three
row-units evaluated at the origin. Consequences:

- **Keller ⟹ A(0) ≠ 0, C(0) ≠ 0, E₀(0) ≠ 0.** The constant terms — the
  "unit-ness" of the cascade — are NECESSARY. Since c₂ = 0 forces A = v^{k+1}
  (THM-1305), this reads: **v(0) ≠ 0**, i.e. v must be an AFFINE (unit) form,
  never a homogeneous one.
- **The Drużkowski shape is banned:** cubic-linear deviations use cubes of
  linear forms WITHOUT constant terms; in-class that means v(0) = 0, hence
  det ≡ 0 — not merely "no counterexample" but **no Keller map at all**.
- **The witness's −2 decodes:** −(E₀(0)=2)(A(0)=1)(C(0)=1) = −2. With this,
  EVERY constant in the owner's expression is derived (completing THM-1305
  §3's table, whose −2 row was the one unexplained datum). The W2 conjecture
  ("|det| = orbit degree k") becomes: Keller-consistency forces
  |E₀(0)A(0)C(0)| = k — true at the only existing rung (k = 2, E₀(0) = 2).

## 2. Equivariant cubic-linear = triangular over the weight-0 core (P1)

H_i = ℓ_i³ homogeneous under a weighting needs w_i = 3·w(ℓ_i); along any
dependency cycle this iterates to w = 3^m·w, forcing weight 0. Exact check
over all nontrivial |w| ≤ 9 (6858 vectors): zero nonzero-weight cycles —
every equivariant Drużkowski map is triangular over its torus-fixed core.
In the witness's weights (1,−1,−2) the relation w_i = 3w_j has NO solutions
at all: **the witness's symmetry class and the Drużkowski stratum intersect
only in the identity**. The two strategies of the century — reduce to cubes
vs. exploit symmetry — are transverse.

## 3. The planar core is shears (P2)

The weight-0 core in dim 3 has dimension ≤ 2, and 2-variable cubic-linear
Keller maps solve exactly (8-equation coefficient system, verified): either
triangular or the shear family ℓ₂ = μℓ₁, a = −μ³b — for which ℓ∘F = ℓ is an
exact identity, so F is a shear along ℓ-level-sets with polynomial inverse
(y − ℓ³, z − μ³ℓ³). No collision can form: fibers of shears are singletons.

## 4. The verdict on the TRAP hypothesis

CONFIRMED, in a form stronger than proposed. The hypothesis was "classical
partial results = strata where the trap polynomial Φ is rootless." The
theorem: on the classical (homogeneous-cube) stratum the trap's raw material
— units with nonzero constant terms — is structurally prohibited, and the
prohibition is visible at the coarsest possible invariant (det ≡ 0 ⟹ not
even Keller). The century was safe not because its maps dodged the root but
because its normal form could not build the unit cascade that makes Φ's
root unavoidable. **The entire distance between the classical stratum and
the counterexample is the constant term of one polynomial** — v = t versus
v = 1 + t. The affine "+1" (the repo's oldest friend: the observer, the
+1 of Rédei = LRC + 1, the 1 in u = 1 + xy = the formal-group denominator)
is exactly what the classical reduction discards and exactly where the
counterexample lives.

## 5. Scope, honestly

All statements are about the EQUIVARIANT slice. Non-equivariant Drużkowski
maps (where the classical rank/nilpotency theorems operate) are not covered:
there the Φ-calculus has no invariant plane to live on. The bridge statement
worth one future session: stable equivalence adds variables — track what the
witness's unit cascade becomes under the Drużkowski REDUCTION (dimension
goes up, the torus extends trivially) and whether the constant terms are
forced into the new variables' homogeneous parts (they must go somewhere;
"where does the +1 hide in cubic-linear form?" is well-posed and computable).

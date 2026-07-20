---
id: THM-1325
title: Where the +1 goes under Drużkowski reduction — (D1) the unit constants ARE the Yagzhev identity part (JF(0) = antidiag(A(0), C(0), E₀(0)); G = L⁻¹F is X + H with det 1 and the collision transported; THM-1320's P3 factorization = det JF(0), true but classically trivial — honesty amendment); (D2) TORUS-BLINDNESS: in any equivariant cubic-linear model every collision sits in the torus-fixed core (triangular back-substitution injective), and in ambient dim ≤ 3 equivariant cubic-linear Keller maps are ALL injective; (D3) THE SURGERY WALL: no single-new-variable row substitution stabilizes the witness — every pattern is obstructed by a rank-drop locus (∇F̃₁ ≡ 0 on {y = 0, w = −1}: the +1's OWN zero-hyperplane), det = ⟨∇P, cofactor-field⟩ forced ≡ 0; stabilization is composition-shaped, and the +1 is absorbed into the TAME automorphism data.
status: >
  PROVED/VERIFIED-EXACT — D1 (JF(0) computed; G = L⁻¹F: G(0) = 0, JG(0) = I,
  det ≡ 1, triple collision transported to image (0,0,−1/4)); D2 (the 5-line
  back-substitution argument + P1/P2 of THM-1320; spot-checked 4000 points on a
  triangular-over-shear-core dim-4 sample, zero violations); D3 (ten
  substitution patterns × the full deg-≤5 P-space, 70 unknowns: NO Keller
  completion, det-const functional ≡ 0 on every solution space; the
  obstruction point (5, 0, 7, −1) verified: ∇F̃₁ = 0 and all four cofactors
  vanish exactly). INTERPRETIVE (labeled) — the conservation-of-+1 reading.
source: death-star-2026-07-19-S59p (HYP-8120; owner directive: backlog xlvi)
depends_on:
  - THM-1300, THM-1305, THM-1320
scripts:
  - 04-computation/jacobian_plus_one_tracking_deathstar_S59p.py -> 05-knowledge/results/jacobian_plus_one_tracking_deathstar_S59p.out
  - 04-computation/jacobian_partial_substitution_deathstar_S59p.py -> 05-knowledge/results/jacobian_partial_substitution_deathstar_S59p.out
---

# THM-1325 — the +1 under reduction: it is the X; the reduction hides the torus; surgery hits a wall

## 1. D1 — the +1 does not hide: it is the Yagzhev X

JF(0) = antidiag(A(0), C(0), E₀(0)) = [[0,0,1],[0,1,0],[2,0,0]] — the six-function
form's unit constants ARE the witness's linear part L. Normalizing,
G := L⁻¹∘F satisfies (verified): G(0) = 0, JG(0) = I, **det JG ≡ 1**, and the
triple collision transports to G with image (0, 0, −1/4). So in Yagzhev normal
form G = X + H, the +1's sit in plain sight as the identity X, and every
stabilization that preserves the X + H shape carries them there forever.

**Honesty amendment to THM-1320 (P3).** The factorization
c₀(0) = −E₀(0)A(0)C(0) is exactly det JF(0) — for ANY Keller map det ≡ det of
the linear part, and "unit constants nonzero" is just "linear part invertible."
The identity is true and engine-verified, but classically trivial; THM-1320's
real content is the FRAME (the identification of det's factors with the unit
cascade's constants, row by row) and P1/P2 (transversality; the planar core).
Logged as MISTAKE-197.

## 2. D2 — what the reduction actually hides: the torus

In any equivariant cubic-linear model (any dimension), the dependency order of
THM-1320 P1 makes the nonzero-weight layer triangular: if two points collide,
back-substitution up the triangle (F_i = x_i + ℓ_i³ with ℓ_i depending only on
later variables and the core) forces their triangular coordinates equal once
their core coordinates are — so **every collision lies in the torus-fixed
core**. In ambient dim ≤ 3 with a nontrivial torus the core has dim ≤ 2, where
cubic-linear Keller maps are shears/triangulars (THM-1320 P2) — injective. So:

- dim ≤ 3: equivariant cubic-linear Keller maps are ALL injective;
- any dim: the torus direction is collision-free — **a cubic-linear model can
  carry the witness's collision only where its torus sees nothing.**

The witness's own presentation has the collision ON a torus orbit (λ = ±1 of
the doubling); a Drużkowski presentation must bury it in the fixed locus. The
two normal forms are informationally complementary: cubes-form shows the
algebra (X + ℓ³) and blinds the symmetry; the equivariant form shows the
geometry (fixed + doubled orbit) and cannot reach cubes (P1: transversality).

## 3. D3 — the surgery wall (a genuine discovery of this session)

Attempting the obvious stabilization — one new variable w, substitute
u → 1 + w in various occurrences, solve the 4th row P for det ≡ const — fails
for EVERY pattern (ten patterns × the full 70-dimensional deg-≤5 P-space):
the det-const functional vanishes identically on each solution space.

The mechanism, exact: det J F̃ = ⟨∇P, c⟩ where c is the (divergence-free)
cofactor field of rows 1–3 — the derivative of P along the fiber direction of
(F̃₁, F̃₂, F̃₃). Keller needs this ≡ const ≠ 0, impossible if c has a zero. And
it does: substituting the affine form 1 + w manufactures one — on
{y = 0, w = −1}, ∇F̃₁ vanishes IDENTICALLY (verified at (5,0,7,−1): gradient
zero, all four cofactors zero). **The +1's own zero-hyperplane w = −1 is the
saboteur**: in the original map the unit u = 1 + xy never annihilates a row
(Keller certifies rank everywhere), but once its argument is a bare new
variable, the unit's vanishing locus becomes a coordinate hyperplane that
meets everything and kills the cofactor field. The unit-ness that BUILDS the
trap (THM-1320) also FORBIDS casually linearizing it.

Consequence: stabilization of the witness is composition-shaped —
F̃ = φ∘(F ⊕ id)∘ψ with φ, ψ tame automorphisms (padding first, surgery never)
— and under composition the unit constants are absorbed into the linear parts
of φ, ψ: **in cubic-linear normal form the +1 lives in the TAME conjugating
data, while H carries the wild residue.**

## 4. The conservation-of-+1 law (interpretive, labeled)

Across every presentation move computed in S59n–S59p, the +1 relocates but
never vanishes: unit constant in the cube (equivariant form) ↔ the Yagzhev
identity X (normalized form) ↔ affine part of ℓ / base-point constant
(attempted linearizations) ↔ tame automorphism data (composition-shaped
stabilization). The repo's oldest symbol — the observer's +1, the 1 of
u = 1 + xy = the formal-group denominator (HYP-1992) — is, in this corner of
affine geometry, a conserved quantity of normal-form changes. What the
century's reduction gained in algebraic simplicity (cubes of linear forms) it
paid for by shattering the torus and exiling the +1 into bookkeeping it never
looked at.

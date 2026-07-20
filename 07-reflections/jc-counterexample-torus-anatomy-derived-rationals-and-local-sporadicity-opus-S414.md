# The JC counterexample's anatomy: torus normal form, derived rationals, quotient contraction, and local sporadicity

**Instance:** opus-2026-07-19-S414 (owner: deepest possible understanding of JC
counterexamples — sporadic or families; trace every rational through our threads; be
wild but grounded). **HYP-8070.** Scripts + frozen outs:
`jacobian_torus_structure_and_family_opus_S414.py` (structure; long symbolic run),
`jacobian_family_tangent_space_opus_S414b.py` (deformation theory — the decisive one).
**Priority note:** mac-mini-S130 pushed a collision-family closed form (xy = −3/2,
identity in s) minutes before my checkpoint — their priority on the λ-family/derived-
rationals slice; everything here was derived independently and the torus normal form,
quotient contraction, and moduli computation below are complementary to their headline.

## 1. The torus anatomy (the deepest structural layer found)

F is equivariant for the FULL one-parameter torus (x,y,z) ↦ (λx, λ⁻¹y, λ⁻²z), with
target weights (−2,−1,1). The weight-0 invariants are **s = xy, w = x²z**, and F has the
normal form

> **F = ( z·P + y²·Q , y·R + xz·S , x·T )**, P,Q,R,S,T ∈ ℂ[s,w]:
> P = (1+s)³, Q = (1+s)(4+3s), R = 1+12s+9s², S = 3(1+s)², T = 2 − 3s − w.

Consequences, all verified:
- **The published points are the λ = 1 slice of a collision continuum:** for EVERY λ ≠ 0,
  F(λ, −3/(2λ), 13/(2λ²)) = F(0, 0, −1/(4λ²)) = (−1/(4λ²), 0, 0). The S413 "second
  fiber" mirror pair (±i/2, ±3i, −26) is the same orbit at λ = ±i/2. One exotic
  ℂ*-orbit double-covers the target axis (λ ↦ λ⁻² is 2:1 — the det's "2" is the torus
  weight), the rigid fixed line supplies the third preimage: the whole counterexample is
  ONE orbit + one line.
- **The quotient contraction:** on invariants (target invariants YZ, Z²X), F descends to
  q(s,w) = (T·(sR+wS), T²·(wP+s²Q)) — and q CONTRACTS the entire curve {T = 0} to the
  origin (both components carry T-factors). det J(q) is non-constant: the plane shadow
  is a contraction, not a plane Keller map — JC₂ is untouched by this route.
- **The rationals are DERIVED, not data:** the exotic orbit = {T=0} ∩ {sR+wS=0};
  substituting w = 2−3s collapses sR+wS to the LINEAR equation **4s + 6 = 0**, so
  s = −3/2 and w = 2 − 3s = **13/2**. The −1/4 is nothing but the λ=1 slice coordinate
  (the collision runs over the whole axis); det −2 is the weight; "3" is the cube/weight
  structure; "4" is Q's constant term, whose only job is to make 4s+6 have a root.
  The owner's number-archaeology closes: within the example every rational is
  structural. Across the repo, the echoes (13 ↔ our 13 speeds; 1/4 ↔ the mod-4 lore;
  3/2, 26 = 2·13) grade as CHARM, not structure — recorded so nobody mines them as
  more.

## 2. Sporadic or family? Locally sporadic, modulo equivalence (evidence-grade)

Deformation theory in the equivariant coefficient space (15 coefficients, det ≡ const):
**16 constraint equations, rank 10 at F ⟹ tangent space of dimension 5.** Of the five
kernel directions: **dir0 integrates exactly — and is the trivial target rescaling**
(it scales (P,Q) jointly; det moves linearly −2 − 2t/3; the "family" is F composed with
a diagonal linear automorphism). **The four nontrivial directions are all OBSTRUCTED at
second order.** Together with the trivial reparametrization group (source/target
weighted rescalings) accounting for the remaining tangent dimensions, the computation
says: **within this weight type and degree box, F is locally rigid modulo equivalence —
a SPORADIC object, whose apparent families are reparametrization orbits** (the
λ-continuum of collisions, the rescaling line). Honest scope: local, ansatz-bounded,
second-order; a full rigidity proof needs the complete trivial-group quotient and
higher-order obstruction theory — the named follow-up. The complementary route to
genuinely NEW examples is therefore not deformation but NEW WEIGHT TYPES: the
construction template is now explicit — pick weights, write the invariant-module normal
form, demand det ≡ const, make T a nonunit, and intersect {T=0} with {sR+wS=0} — a
MACHINE whose other instances (other tori, other invariant rings) are the real hunt.

## 3. What our mathematics contributed (the two-way flow)

The mirror-parity lemma (S413) generalizes to the full torus: fibers over the axis are
{fixed-line point} ⊔ {orbit pairs} — the 1+2k law is the ℤ/2 ⊂ ℂ* shadow. The
contraction-of-{T=0} is a blowdown seen through invariant theory — the repo's
observer-lens habit (mark the origin, read the local view) is exactly what surfaced
s = xy, w = x²z. And the duty-calculus instinct transfers: T's constant term 2 pays for
the collision the way K pays the ghost duty — "the map affords one contracted curve,
priced at det = −2." (Rhyme-grade, so labeled.)

## 4. Handoffs

(a) The weight-type hunt: run the machine at other weight systems (the template is in
§2) — each success is a genuinely new counterexample, each failure maps the rigidity.
(b) Complete the sporadicity proof (trivial-group quotient + obstructions). (c)
mac-mini: reconcile/merge with your S130 closed form — priority yours on the λ-family;
the normal form and moduli here compose with it. (d) The long symbolic transcript
(part-A script) lands asynchronously; its W1–W4 sections re-verify what is hand-proved
and ansatz-confirmed above.

Cross-links: HYP-8065/S413 (verification + parity) · mac-mini-S130 (priority, collision
closed form) · THM-378 (the ×2-parity functor rhyme) · S410 rung-theory duty calculus
(the pricing rhyme) · scripts + frozen outs (S414, S414b).

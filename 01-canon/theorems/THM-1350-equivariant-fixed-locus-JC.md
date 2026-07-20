---
id: THM-1350
title: A REDUCED JACOBIAN CONJECTURE THAT DOES HOLD — THE EQUIVARIANT FIXED-LOCUS JC — AND IT FORCES THE COLLISION MULTIPLICITY TO BE ODD, HENCE ≥ 3 — for a σ/τ-equivariant F with det JF ∈ ℂ*, F maps Fix(σ) → Fix(τ) and the restriction is again a constant-Jacobian polynomial map in dimension dim Fix(σ); when dim Fix(σ) ≤ 1 this is JC₁, which is TRUE, so **F|_Fix is injective unconditionally** and every counterexample's non-injectivity is supported entirely on FREE orbits. On the owner's dim-3 counterexample this is exact: σ = (−x,−y,z), τ = (a,−b,−c), F(σp) − τ(Fp) ≡ 0, and **F(0,0,z) = (z,0,0)** — the restriction to Fix(σ) is the IDENTITY. Consequently every fiber over a τ-fixed target contains EXACTLY ONE σ-fixed preimage, so |fiber| ≡ 1 (mod 2): an equivariant counterexample collision can never be a DOUBLE, and the minimum is a TRIPLE — which is precisely what the owner's example realises (1 fixed + 1 free orbit)
status: VERIFIED-EXACT on the counterexample (σ-equivariance symbolic, F|_Fix = identity, det JF = −2, the triple collision in exact rationals, the two non-fixed preimages a single σ-orbit). The reduced JC is unconditional for dim Fix(σ) ≤ 1 (JC₁) and conditional on JC₂ for dim Fix(σ) = 2. CREDIT: the σ-equivariance and the "Rédei-shaped 1+2 fiber" were already identified by kind-pasteur-S128c97; what is added here is the FORCING — F|_Fix bijective ⟹ odd fiber ⟹ collision ≥ 3 — and the reading of it as a surviving reduced JC
source: opus-2026-07-20-S399 (owner: creatively produce a reduced Jacobian Conjecture that does hold)
depends_on: [THM-1300 (the counterexample; JC false for n ≥ 3), kind-pasteur-S128c97 / HYP-8070 (σ-equivariance and the 1+2 fiber), klein PROBLEM-LEDGER-S332 ("Aut-group conjectures SURVIVE — mark the live boundary")]
scripts: 04-computation/jc_fixed_locus_opus_S399.py -> 05-knowledge/results/jc_fixed_locus_opus_S399.out
---

# THM-1350 — what still holds after JC falls

> **PROOF GAP REPAIRED (opus-S400).** As filed, this asserted that F|Fix has constant nonzero Jacobian without justifying it -- det JF constant bounds only the PRODUCT of the two sigma-blocks (the +1 eigenspace tangent to Fix, and the -1 eigenspace normal to it). THE MISSING STEP: along Fix, det JF = det(tangential) * det(normal), both POLYNOMIALS in the Fix coordinates; their product is a nonzero constant, and in a polynomial ring (a domain) a product equal to a nonzero constant forces EACH factor to be a nonzero constant, since degrees add. Hence the tangential block is Keller and the claim stands. Verified on the counterexample: tangential det = 1, normal det = -2, product = -2 = det JF. Prompted by kind-pasteur-S128c105's self-refutation of a DIFFERENT reduced JC (the weight-twist version, killed because sigma.F is untwisted, Keller and non-injective) -- that refutation does NOT apply here, since this claims only injectivity of the restriction, not of F.

JC is false for n ≥ 3 (THM-1300). The useful question is then which reduced
form survives. Here is one that does, and it constrains the counterexamples.

## The reduced conjecture (a theorem for dim Fix ≤ 1)

> **Equivariant fixed-locus JC.** Let σ, τ be involutions of source and target
> and F a σ/τ-equivariant polynomial map with det JF ∈ ℂ*. Then
> F(Fix σ) ⊆ Fix τ, and F|_{Fix σ} is a polynomial map with constant nonzero
> Jacobian in dimension dim Fix(σ). Hence **if dim Fix(σ) ≤ 1, F|_{Fix σ} is
> injective** — unconditionally, since JC₁ is trivial. (For dim Fix(σ) = 2 it
> is exactly JC₂, still open.)

So all non-injectivity of an equivariant counterexample lives on **free
orbits**; none of it is visible on the fixed locus.

## Verified exactly on the owner's counterexample

With F₁ = (1+xy)³z + y²(1+xy)(4+3xy), F₂ = y + 3x(1+xy)²z + 3xy²(4+3xy),
F₃ = 2x − 3x²y − x³z, σ(x,y,z) = (−x,−y,z), τ(a,b,c) = (a,−b,−c):

| check | result |
|---|---|
| F(σp) − τ(F p) | **≡ 0** (symbolic) |
| Fix(σ) = {x=y=0} → Fix(τ) = {b=c=0} | **yes** |
| **F(0,0,z)** | **= (z, 0, 0)** — the IDENTITY on Fix(σ) |
| det JF | −2 |
| d/dz of F\|_Fix | 1 |
| triple collision | exact, over (−1/4,0,0) |
| the two non-fixed preimages | a single σ-orbit |

## The forcing: collisions are odd, hence ≥ 3

Over a τ-fixed target, σ acts on the fiber, so

> |fiber| ≡ #Fix(σ|fiber) (mod 2).

Since F|_{Fix σ} is **bijective** here (it is the identity), each such fiber
has **exactly one** σ-fixed preimage. Therefore

> **|fiber| is ODD.**

An equivariant counterexample can therefore never exhibit a *double*
collision — the minimum is a **triple**, realised as 1 fixed point + 1 free
orbit. That is exactly the shape of the owner's example, and it explains why
it had to be a triple rather than the simpler double one might have searched
for first.

## Credit and scope

The σ-equivariance and the "Rédei-shaped 1+2 fiber structure" were already
recorded by **kind-pasteur-S128c97**; this file does not claim them. What is
added is (i) the reading of the fixed-locus restriction as a **surviving
reduced JC**, unconditional at dim Fix ≤ 1, and (ii) the **forcing argument**
— bijectivity on the fixed locus ⟹ odd fiber ⟹ collision multiplicity ≥ 3.

This also sits on the boundary klein's PROBLEM-LEDGER-S332 flagged
("Aut-group conjectures survive"): the surviving statement is precisely the
one about a map that *is* bijective (on the fixed locus), while the failure
is confined to the endomorphism behaviour off it.

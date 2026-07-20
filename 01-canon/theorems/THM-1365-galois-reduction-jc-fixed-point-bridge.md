# THM-1365: The Galois Reduction — a Jacobian Conjecture that holds (polynomial-deck case), and the fixed-point bridge

**Status:** VERIFIED (modulo three classical citations, named below)
**Author:** boxeph-2026-07-20-S150 (HYP-8210)
**Companion:** `07-reflections/the-galois-reduction-and-deck-poverty-boxeph-S150.md`

## Statement

Call a Keller map F: C^n → C^n (det JF ∈ C*) **polynomial-Galois** if its deck
transformations are polynomial automorphisms of C^n acting transitively on
generic fibers (the extension C(x)/F*C(u) is Galois and its deck group G is
realized inside Aut(C^n)).

**(A) Cartan freeness (all n).** For any Keller F, the polynomial deck group G
acts FREELY on C^n.
*Proof.* Let g ∈ G fix p. Differentiate F∘g = F at p: dF_p·dg_p = dF_p, and
dF_p is invertible (Keller), so dg_p = id. A finite-order automorphism with
identity differential at a fixed point is the identity in a neighborhood
(Cartan's uniqueness lemma: average h = (1/|g|)Σ_k (dg_p)^{-k}∘g^k conjugates
g to its linear part = id), hence g = id by Zariski density. ∎

**(B) Dimension 2: polynomial-Galois JC₂ HOLDS.** Every polynomial-Galois
Keller map F: C² → C² is an automorphism.
*Proof.* By (A), G acts freely. Every finite subgroup of Aut(C²) is
linearizable [CITATION 1: Kambayashi 1979, via van der Kulk amalgam + Serre's
tree fixed-point theorem], so a nontrivial G would be conjugate to a linear
group, which fixes the origin — contradicting freeness. Hence G = 1; Galois
transitivity forces generic degree d = |G| = 1, so F is birational. A
birational Keller map is an automorphism [CITATION 2: Keller 1939 — the
birational case of JC is classical; alternatively étale + ZMT open immersion +
Ax surjectivity]. ∎

**(C) Dimension n: the fixed-point bridge.** Define FP(n): "every nontrivial
finite subgroup of Aut(C^n) has a fixed point." Then FP(n) ⟹ polynomial-
Galois JC_n (same proof; (A) needs no hypothesis). FP(2) is true [CITATION 1].
FP(3) is known in substantial cases (linearizable classes) and open in
general; FP(n) for large n is related to the linearization problem [CITATION
3: Kraft's linearization problem survey]. Polynomial-Galois JC_n is thus NOT
weaker than a fixed-point theorem — it IS one.

**(D) The evasion certificate for the kernel.** The external counterexample K
has monodromy S₃ on d = 3 sheets (in-repo Chebotarev census, S140). Its deck
group is N_{S₃}(S₂)/S₂ = 1 (S₂ is self-normalizing): K is maximally
non-Galois. The counterexample evades (B)/(C) with the smallest possible deck
group — pure group theory, no computation needed.

## The Deck-Poverty Conjecture (new, named)

**Every Keller counterexample (any n) has trivial polynomial deck group.**
Equivalently, its monodromy point-stabilizer is self-normalizing. This is
implied by the full-symmetric-monodromy conjecture (S144: monodromy = S_d
always; then N_{S_d}(S_{d−1}) = S_{d−1}). It sits with the Ghost Theorem
(S145: no regular-monodromy all-ghost counterexamples) in one structure law:
**counterexamples are monodromy-rich and deck-poor.** The three known
salvage theorems (equivariant JC₂ = THM-1345, Ghost, and this Galois
reduction) each cut from a different side: symmetry, topology, Galois theory.

## Degree-2 remark (route to unconditional JC₂(d=2))

A degree-2 Keller map is automatically Galois with deck involution
σ = t∘F − id (t = fiber trace), regular exactly off F^{-1}(A(F)); escaping
sheets force genuine poles of t along A(F), so the polynomial-deck hypothesis
fails precisely through the asymptotic set. σ is fixed-point-free on the étale
locus (deck fixed points = ramification = ∅). Completion blueprint: run σ
through the Bayle–Beauville classification of birational involutions of P²
(linear / de Jonquières / Geiser / Bertini — all carry fixed curves); locate
the forced fixed locus against the boundary divisor + Euler-ledger profile
(the only surviving d=2 configuration is one asymptotic component with
(k,χ) = (1,1)). This is filed as the S150 handoff toward "JC₂ holds in degree
2" unconditionally, complementing the equivariant and smooth-caustic kills.

## Dependencies and honesty

- CITATION 1 (Kambayashi linearization / Serre tree argument): classical,
  stated without in-repo proof.
- CITATION 2 (birational Keller ⟹ automorphism): classical (Keller's paper).
- CITATION 3 (Kraft linearization survey): context for FP(n), no claim used.
- (A) and (D) are complete in-repo arguments. (B) is complete modulo 1+2.
- No computation was needed; the S140 census supplies (D)'s S₃ input.

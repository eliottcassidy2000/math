# A genuine 57-to-Heisenberg map, and the E7 representation obstruction

Session35, 2026-10-23. **PROVED**, conditional only on the explicitly cited standard E7 restriction datum; finite certificates are separate. These are elementary constructions/applications, not a claim of new finite-subgroup classification.

## 1. Recovering the missing action

The shared integer57 alone did not connect the Singer projective plane, Blum arithmetic and E7's Heisenberg radical. A more faithful source is the **57-point cyclic action**, not its cardinality.

Earlier `01-canon/theorems/THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate.md` shows that retaining an actual Singer/Heisenberg action can expose a lost fixed-point constraint. That theorem concerns finite characteristic-p affine actions, not the characteristic-zero Lie algebra here. Its method transfers, not its conclusion. The current session's Singer difference-set note supplies a C57 action on PG(2,7), with Frobenius multiplication by7.

Let W be the augmentation subspace of Q^57, the vectors whose coordinate sum is zero. Let P be the cyclic permutation matrix. Then dim W=56. The alternating matrix P-P^(-1) has kernel exactly the constant vectors: on a complex Fourier character its eigenvalue is zeta^j-zeta^(-j), zero only when 2j=0 mod57. Since57 is odd, only j=0 vanishes. It therefore gives a nondegenerate alternating form on W.

Define H=W⊕Qz with [u,v]=omega(u,v)z. This is a57-dimensional Heisenberg Lie algebra, with the distinguished C57 action acting symplectically on W and fixing z. It is not a group homomorphism from C57 into a characteristic-zero unipotent Heisenberg group; it is an action by Lie algebra automorphisms. The abstract Heisenberg Lie algebra is isomorphic to E7½'s radical. Whether the distinguished finite action extends to E7 is a further question.

## 2. Retaining Frobenius as well

The simple form P-P^(-1) is not automatically preserved by the point permutation i->7i. Repair it by taking

    Omega = sum_(a in {1,7,49}) (P^a-P^(-a)).

Multiplication by7 cyclically permutes these three exponents modulo57, so Omega is invariant under both translations and Frobenius. It is alternating. The exact polynomial certificate

    gcd(sum_(a in {1,7,49})(x^a-x^(57-a)), x^57-1)=x-1

proves that its kernel remains just the constants. A rational Bezout certificate is printed by the companion script, so this is not an inferred numerical rank. Thus the whole semidirect product C57⋊C3, where C3 acts by exponent7, acts symplectically on W and by automorphisms of H. The order171 action preserves the specified Singer line set as well.

This is a genuinely operation-preserving chain:

    Singer-labelled projective points -> augmentation module W -> invariant alternating form -> Heisenberg bracket.

It requires the cyclic labelling and the chosen Frobenius-stable exponent orbit. It forgets the actual incidence relation when passing just to W and Omega; any claim to reconstruct the plane needs the difference-set sidecar.

## 3. The action cannot extend to the E7 minuscule representation

**Theorem.** This C57 module W is not the restriction of the56-dimensional minuscule representation along any homomorphism C57->E7(C), where E7 is the simply connected group acting on that representation. Consequently the displayed C57⋊C3 action cannot so extend either.

Proof. In the regular permutation representation of C57, the element P^19 has order3 and19 cycles. Removing the constant line gives fixed dimension18 and trace−1 on W.

Use the standard maximal-rank subgroup SL8/{±I} in E7, under which the minuscule representation restricts to exterior^2(C^8)⊕exterior^2(C^8)*. Every order3 element of the complex group is semisimple and is conjugate into a maximal torus of this subgroup. Its lift to SL8 can be chosen of order3: any lift with cube−I becomes cubeI after multiplication by−I.

Let a,b,c be multiplicities of1,omega,omega² in the8-dimensional lift. Then

    a+b+c=8,  b+2c=0 mod3.

The fixed dimension on exterior^2 plus its dual is

    d0=2[binomial(a,2)+bc].

Enumerating the finitely many triples gives precisely

    d0 in {2,14,20,32,56}.

The corresponding traces (3d0−56)/2 are −25,−7,2,20,56. None has fixed dimension18 or trace−1. Contradiction.

The finite enumeration is fully exact and is also the familiar order3 trace restriction in the finite-subgroup literature. The proof supplied here uses the standard SL8 branching, not an unaudited table lookup.

This obstruction is stronger and more informative than observing that3·19 is a Blum integer. It identifies exactly what fails: the representation character of the preserved order3 subgroup. There can still be other C57 subgroups of E7 with different56-dimensional character; no blanket prohibition of elements of order57 is asserted. No obstruction is claimed to forgetting the C57 action and identifying the two abstract Heisenberg algebras.

## 4. Scope and tests

The construction and E7 obstruction above are asserted in characteristic zero. In characteristics3 and19, permutation splitting and semisimplicity can fail; no positive-characteristic E7 extension claim is made. Neither a symplectic action nor a shared dimension constructs the E7-invariant quartic tensor. The order3 character already prevents that particular extension.

`check_singer_heisenberg.py` prints the Bezout identity, all allowed SL8 multiplicity triples, the fixed-space obstruction, and verifies the explicit Singer difference set D={8,15,19,38,48,50,51,56} and Frobenius stability. Normal and optimized outputs agree. `singer_heisenberg.json` is the exact receipt.

Literature boundary: the maximal-rank SL8 restriction is a standard structural input. Source pin and independent audit accompany the final publication; no inference from an inaccessible PDF abstract is needed for the finite calculation. The initially located Griess--Ryba university PDF link redirected unsuccessfully; that failure is not represented as a full-text reading.

Independent source recovery: the coauthor-uploaded full text of Griess--Ryba, *Classification of finite quasisimple groups which embed in exceptional algebraic groups*, printed page39, states the four nonidentity order3 traces. The audit inspected this passage and independently derived the SL8 calculation:
https://www.researchgate.net/publication/243106241_Classification_of_finite_quasisimple_groups_which_embed_in_exceptional_algebraic_groups
See the separate Singer addendum audit for source retrieval and exact matrix checks.

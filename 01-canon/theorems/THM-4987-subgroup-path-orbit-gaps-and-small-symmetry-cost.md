---
id: THM-4987
title: "Every automorphism subgroup retains the 7/21 path-orbit gaps; small symmetry has a sharp threefold cost"
status: >
  PROVED, independently audited in Session232. FINITE-EXACT controls have
  their stated universes; inherited all-order theorem certificates remain
  explicit dependencies, not newly rerun claims. Chan-Pak completeness OPEN.
session: Session232-automorphisms-20261109
source: 05-knowledge/results/session232_automorphisms_20261109/subgroups/PROOF.md
audit: 05-knowledge/results/session232_automorphisms_20261109/subgroup_transfer_audit/AUDIT.md
---

# THM-4987

[Full proof](../../05-knowledge/results/session232_automorphisms_20261109/subgroups/PROOF.md) · [Independent audit](../../05-knowledge/results/session232_automorphisms_20261109/subgroup_transfer_audit/AUDIT.md)

Let T be a finite nonempty tournament, H(T) its directed Hamilton-path count,
G=Aut(T), and A<=G. Put q_A=H(T)/|A|. This integer counts A-orbits of ordered
spanning Hamilton paths. The following are PROVED, relying on the explicitly
cited inherited THM-4812 finite-certificate theorem, not a fresh replay of it.

1. q_A is never 7 or 21, for every A, not just G or the trivial group.
2. If q_A<23, primes dividing |A| belong to {3,5,7}. Set
   b=3^(v_5(|A|)+2*v_7(|A|)). Then q_A=b or q_A>=3b.
3. Equality q_A=b holds precisely when A=G and T is an ordinal sum of
   K1,C3,R5,P7, with v_5(|A|) R5 blocks and v_7(|A|) P7 blocks. Here
   R5=Cay(Z/5,{1,2}), P7=Cay(Z/7,{1,2,4}); the scope q_A<23 remains.
4. Therefore 7 divides |A| and q_A<23 forces q_A=9 and exactly one P7 among K1/C3
   blocks; 25 divides |A| and q_A<23 forces q_A=9 and exactly two R5 blocks;
   5 divides |A| and q_A<9 forces q_A=3 and exactly one R5.

## Proof mechanism

The full-G vertex orbits O_i induce vertex-transitive tournaments T_i. Put
P=product Aut(T_i), d=[P:G]. Faithful restriction gives G<=P, and positivity
of the odd-cycle collection formula gives

    q_A=[G:A]q_G >= [G:A] d product_i q(T_i).

THM-4812(b),(c), including its later audit repairs, classify the possible
vertex-transitive pieces when q_A<23 as K1,C3,R5,P7 with quotient counts
1,1,3,9 and group orders1,3,5,21. If r,p count R5/P7 pieces, then
q_A>=[G:A]d*3^(r+2p) and r>=v_5(|A|), p>=v_7(|A|).
This proves the prime support and price bound q_A>=b.

If q_A<3b, the odd indices [G:A],d must both be1 and both valuation bounds
must be equalities. Thus A=G=P. Independent transitive actions on the pieces
make them modules. If their quotient tournament contained a directed triangle,
its three whole-block Hamilton-word orders give three disjoint families,
contributing at least three times the product of those block counts. OCF
positivity with the remaining blocks would give q_A>=3b, a contradiction.
The quotient is transitive, so the pieces are the ordered SCCs. Their path
counts and groups multiply, giving q_A=b and the stated classification.
Conversely those ordinal sums attain equality.

For the global gaps, q_G divides q_A. The inherited exclusion q_G!=7,21
leaves, if q_A were7 or21, only q_G in{1,3} and a subgroup index divisible7.
But the price bound with A=G makes 7 divides |G| imply q_G>=9 in this low range.
Contradiction. This is an index obstruction, not a claim that allowed odd
integers are closed under arbitrary multiplication without producing7/21.

## Boundaries and verification

Automorphism orders7 and21 themselves occur: the carousel Q7 has
(H,|G|,q_G)=(175,7,25), and P7 has (189,21,9). The Q7 hostile prevents extending
the seven-symmetry classification to unrestricted quotient counts. The
strict threshold23 comes from the inherited vertex-transitive lemma; no
classification at or beyond23 is asserted.

The source packet's optional module-interaction lemma gives normalized costs
27 for a triangle through an R5 module and105 through P7, under the explicit
full-product orbit-module hypothesis. These are not bounds for arbitrary
seven-symmetric tournaments. They are not needed above.

The finite controls enumerate all532 classes through order7 and all640
literal subgroup/class pairs. Separate audit code checks all1099 labelled
tournaments through order5, all1301 subgroup pairs there, and named larger
sharp controls. These checks do not replace THM-4812's inherited all-order
certificates. Its unaudited part(g) is not used. See the full proofs and audits.
Chan–Pak completeness remains OPEN; only omission and symmetry strata are proved.

Inherited routes: THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md;
THM-002-ocf.md; LEM-003-aut-acts-freely-on-ham-paths.md, all in this directory.

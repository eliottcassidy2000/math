# All four odd eighths from affine-cube symmetry

**Session 115. PROVED scoped classifications; finite exact controls.** No literature-priority claim. This packet extends the weighted group carrier of THM-4726, not the enumeration or realization of arbitrary tournament switching classes.

Write rho(G) for the proportion of odd-order elements of a finite group G, and E(G)=1−rho(G). The same observable produces all four odd eighths in one completely classified family. They are not being identified with LRC times or Collatz densities.

## 1. A uniform formula in every binary dimension

Let V=F_2^d and let H be any finite odd-order subgroup of GL(V). For G=V semidirect H,

    rho(G) = (1/|H|) sum_(h in H) 2^(−dim Fix_V(h)),
    E(G) = 1−(1/|H|) sum_(h in H) 2^(−dim Fix_V(h)).       (1)

**Proof.** Fix h of odd order l. The affine map (v,h) raised to l is translation by

    N_h v=(1+h+...+h^(l−1))v.

It has odd order exactly when N_h v=0. In characteristic two the norm N_h is the projection onto Fix(h): on the fixed space it is multiplication by the odd integer l, hence identity, and its image is fixed by h. Therefore exactly 2^(d−dim Fix(h)) translations v work. Average over h. No independence approximation or uniform resampling after quotienting is used. QED.

A useful arbitrary-dimensional subfamily is a fixed-point-free odd linear complement H: if every h!=1 fixes no nonzero vector, then

    E(V semidirect H)=(1−2^(−d))/|H|.                     (2)

For any divisor h of 2^d−1, the order-h multiplicative subgroup of F_(2^d)^* gives such an action. At h=2^d−1 the even-order mass is exactly 2^(−d); at h=1 it is 1−2^(−d). These are endpoint formulas for a specified action family, not all groups of the same order.

## 2. Complete dimension-three classification

Up to conjugacy, an odd-order subgroup H of GL_3(2) has exactly one of four types:

| H | order | numbers of elements by fixed-space dimension | E(F_2^3 semidirect H) |
|---|---:|---|---:|
| trivial |1| dimension3:1 |7/8|
| C_3 |3| dimension3:1, dimension1:2 |5/8|
| C_7 |7| dimension3:1, dimension0:6 |1/8|
| C_7 semidirect C_3 |21| dimension3:1, dimension0:6, dimension1:14 |3/8|

This classifies every odd linear complement in dimension three, not just four chosen examples.

**Proof of completeness and conjugacy.** The order of GL_3(2) is (8−1)(8−2)(8−4)=168. Thus |H| is one of1,3,7,21. An order-three matrix has rational canonical form with one trivial one-dimensional block and one irreducible quadratic block, so all C3 subgroups are conjugate and their nonidentity elements fix exactly one dimension. Order-seven subgroups are Sylow7 subgroups and hence conjugate; every nonidentity element acts irreducibly on V and fixes zero dimensions.

If |H|=21, its Sylow7 subgroup is normal, since its number divides3 and is1 modulo7. Hence H lies in the normalizer of a Singer subgroup C7. Identify V with F8. A normalizing linear map induces an automorphism of the multiplication algebra F8; its field-automorphism component has order dividing3, and its remaining component is multiplication by a nonzero field element. Thus the normalizer is F8^* semidirect Gal(F8/F2), of order21. H equals this normalizer and is unique up to conjugacy. Its six nonidentity C7 elements have order7, and its fourteen remaining nonidentity elements have order3 and one fixed dimension. Formula(1) gives the displayed values. QED.

The direct finite subgroup census has45 odd subgroups: one trivial,28 C3, eight C7 and eight order21 subgroups. This is a check of the proof, not its substitute.

### Complete decoration spectra and a centralizer hostile

Every nonidentity2-element of V semidirect H is a nonzero translation, because its image in the odd quotient H must be trivial. Its conjugacy class is its H-orbit, and its centralizer is V semidirect H_v. Thus the inherited weighted identity has the following exact decoration spectra:

| Complement order | nonzero-vector orbit sizes | centralizer odd-order weights |
|---:|---|---|
|1|1,1,1,1,1,1,1|seven copies of1/8|
|3|1,3,3|3/8,1/8,1/8|
|7|7|1/8|
|21|7|3/8|

For H=C3 the fixed nonzero vector is central, so its centralizer is the entire group, whose odd proportion is3/8 by(1). The other stabilizers are trivial. For H=C7 or order21, transitivity on nonzero vectors gives stabilizers1 orC3 respectively. These facts prove the table.

This gives an explicit hostile to replacing rho(C_G(x)) by the reciprocal Sylow2 order: in V semidirect C3, the central translation has centralizer with Sylow2 order8 but weight3/8, not1/8. The normal2-complement hypothesis cannot be dropped. The same total7/8 can also conceal different decoration spectra: V has seven weights1/8, whereas D8 has weights1/8,1/4,1/4,1/4.

### Precisely what the quotient retains

The source object is an affine transformation (v,h). The map to h forgets v; the preserved counting predicate is **odd order**, after retaining the norm-kernel size 2^(d−dim Fix(h)). The needed sidecar is the fixed-space dimension, not only the order of H. Formula(1) sums over actual h with their multiplicities. The four eighths arise from different fixed-space profiles in the same eight-point translation geometry.

## 3. A second complete classification: every Sylow-D8 group

Here D8 means the dihedral group of order8. Let G be ANY finite group with Sylow2 subgroup D8. Let c be the number of G-conjugacy classes of involutions. Then

    c in {1,2,3},                  E(G)=(2c+1)/8.        (3)

Thus the complete even-order spectrum for this Sylow family is {3/8,5/8,7/8}. The examples GL3(2), S4 and D8 realize c=1,2,3 respectively. A7 also has c=1. Odd direct factors do not change the result.

The precise decoration spectrum is

- one involution class of weight1/8;
- c−1 further involution classes, each of weight1/4;
- one order-four class, of weight1/4.

A weight means rho(C_G(x)), not the fraction of G occupied by the conjugacy class of x. In particular1/8 is a decoration weight in this theorem, whereas it is a total group even-order mass in the affine C7 example above.

### The only imported group theorem

We use Burnside's normal-complement criterion: if an abelian Sylow p subgroup Q is centralized by its normalizer, the group has a normal p-complement. Its hypotheses, including the normalizer condition, are always checked below. This is the usual transfer theorem, not the solvability theorem often called by the same name. A primary university seminar statement is available from [Bilkent](https://math.bilkent.edu.tr/assets/data/Old_seminar_data/seminars/2017Fall/08112016.htm); a proof is in Kurzweil–Stellmacher, *The Theory of Finite Groups*, Theorem7.2.1, printed page169, [university-hosted text](https://theory.uchicago.edu/~sethi/Teaching/P342-W2017/Kurzweil-Stellmacher_Theory%20of%20finite%20groups.pdf).

For clarity, the standard short transfer argument is as follows. If Q is abelian, two of its elements conjugate in the whole group are conjugate under N(Q): the two Sylow subgroups containing the conjugated element lie in its centralizer and can be aligned there. If N(Q) centralizes Q, this fusion is trivial. The transfer homomorphism to Q, evaluated on x in Q by the coset-cycle formula, is then x^[G:Q]. Since the exponent is coprime to p, transfer maps Q isomorphically onto itself; its kernel is the required normal p-complement.

### Centralizer weights, with the missing hypothesis repaired

Fix P=D8 and let z be its central involution.

**Central type.** Put J=C_G(z). Its Sylow2 subgroup is P and z is central in J. The quotient J/<z> has Sylow V4. Its three nontrivial elements have preimages in P that are one C4 and two V4 groups. Every normalizer automorphism preserves the unique C4 preimage, so its action on V4 fixes a nontrivial element and has order dividing2. But the induced automorphism group has odd order: the abelian Sylow V4 lies in its centralizer, so the normalizer-centralizer index is odd. Therefore that action is trivial. Burnside applies, and the quotient has odd-order proportion1/4.

Each odd-order element of J/<z> has exactly one odd-order lift to J: of its two lifts g,gz, exactly one has odd order. Even-order quotient elements have no odd-order lifts. Consequently rho(J)=1/8. This proof computes the weight directly; it does not merely assume J is2-nilpotent.

**Other involution types.** Let x be an involution not G-conjugate to z, moved into P. Its centralizer contains C_P(x)=V4. Its Sylow2 order cannot be8, because x would then be central in a Sylow D8 and hence conjugate to z. Thus its Sylow subgroup Q is V4. The normalizer action on Q fixes x, so again it has order dividing2 and also odd order, hence is trivial. Burnside gives a normal odd complement, and rho(C_G(x))=1/4.

**Order four.** Every order-four element is conjugate into P, and the two such elements of P are conjugate there. Its centralizer has cyclic Sylow C4 (it cannot contain a Sylow D8 with that element central). The normalizer action on C4 has order dividing2 and also odd order; hence it is trivial. Burnside gives rho(C_G(x))=1/4.

There are three involution classes in P, so G can have only1,2 or3, with exactly one containing z. Every2-element is an involution or has order4. The inherited THM-4726 identity

    E(G)=sum_(nontrivial 2-element classes [x]) rho(C_G(x))

now proves(3), including its exact decoration spectrum. QED.

## 4. Switching-class consequence and scope

THM-4726 already proves that each tournament switching class C contributes E(Aut C) to its level branch. Therefore, **whenever** Aut C has Sylow D8, its contribution is exactly(2c+1)/8, with c the number of involution classes in Aut C. This turns the old finite order-eight examples into a uniform conditional classification over every such automorphism group.

The affine-cube theorem classifies group masses for V semidirect H, but it does not prove those groups occur as full automorphism groups of tournament switching classes. Conversely, knowing only a group order or only an unweighted list of decorations does not determine the mass. The old D8 hostile remains essential: conjugacy-class weights are1/4,1/4,1/4,1/8, not four equal1/8 weights.

Neither theorem supplies the still-open denominator-clearing integral packets for the whole switching branch. Fractions from distinct mathematical carriers may only be compared after identifying their observables and weights. There is no asserted map from these group masses to Collatz content or LRC maximizing times.

## 5. Reproduction

Run `python check_eighth_groups.py` and `python -O check_eighth_groups.py`. The checker uses explicit exception gates and direct permutations, not a group-theory package. It exhausts all odd subgroups of GL3(2), forms every affine group, counts actual element orders, and compares the fixed-space formula. Separate controls evaluate centralizer odd-order proportions for D8, S4, GL3(2), A6, A7 and odd direct products. `results.json` records the actual finite census. These tests support the all-family proofs and do not replace the cited transfer theorem.

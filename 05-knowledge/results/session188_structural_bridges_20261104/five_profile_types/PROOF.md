# Five ordered-cover residue types and a complete cyclic coarse-response certificate

## 1. Canonical reduction for every nonempty child

Write H(A)=q_1(A) for a tournament's Hamilton-path count and q_k(A) for its ordered spanning covers by k nonempty directed paths. H(A) is odd. Let z(A) be the number of directed Hamilton cycles up to rotation when |A| is odd, with z=0 for order one.

The exact modulo-eight profile formulas, proved in `../mod8_transfer/PROOF.md`, imply that q(A)/H(A), coordinatewise modulo eight, is one of the following five types. Division means multiplication by the inverse of an odd residue in Z/8Z.

| Type | Exact condition | Canonical actual tournament | Normalized (q1,q2,q3) |
|---|---|---|---|
| S | order one | TT1 | (1,0,0) |
| P | order two | TT2 | (1,2,0) |
| E | even order at least four | TT4 | (1,6,4) |
| O | odd order at least three and z even | TT3 | (1,6,6) |
| C | odd order at least three and z odd | directed C3 | (1,2,2) |

Every q_k for k>=4 is zero modulo eight. All five types are realized, and their displayed normalized profiles are distinct.

Here is the complete derivation of the boundary distinctions. For N>=3, 2^N vanishes modulo eight. The inherited formulas are

    q2 = 2H+2^N-4-4z (mod8),
    q3 = 6S(N,3)                         if N is even,
    q3 = 6S(N,3)+2H-2-4z                if N is odd.

For N>=3, the integer formula 6S(N,3)=3^N-3·2^N+3 gives S(N,3)=1 modulo4 when N is odd and S(N,3)=2 modulo4 when N is even. Indeed the numerator is respectively6 or12 modulo24. Consequently even N>=4 has normalized q2=6,q3=4. Odd N>=3 has q2=q3=2H+4-4z modulo8, and normalizing gives6-4z for both. Orders1 and2 are checked directly. No value of H modulo8 is discarded: it is the odd scalar that is factored out.

If R_i is the actual canonical tournament for A_i, define the odd unit

    lambda_i = H(A_i) / H(R_i) (mod8).

Then q_k(A_i)=lambda_i q_k(R_i) modulo8 for every k. In particular the C type uses an actual C3 with H=3; the normalized sequence(1,2,2) is not being represented as a literal tournament profile.

## 2. Complete block-uniform response transfer

Let Q be any quotient tournament and let its children A_i be arbitrary nonempty named tournaments. Give each of two actual ports x->y uniform incidences to each complete child, encoded by quotient masks(a,c). In each maximal-child-run decomposition of a Hamilton word, child i contributes exactly one factor q_(k_i)(A_i). The run skeleton, port order and empty prefix/suffix cases depend only on Q and(a,c). Therefore the exact response is multilinear in the complete child profiles, one profile factor for every child. Reducing modulo8 gives

    F_(Q[A_i])(a,c)
      = (product_i lambda_i) F_(Q[R_i])(a,c) (mod8).       (1)

On the right the same quotient bits are expanded over the complete canonical child sets. The original physical names are not relabelled or identified; the smaller model is an arithmetic device for selecting the block masks, which are then restored to the actual name sets.

The scalar product in(1) is an odd unit, independent of the query and the quotient orientation. Thus the entire coarse contrast bank between Q[A_i] and Q^op[A_i] is zero modulo8 exactly when the canonical bank is zero. The same statement holds when target child interiors are changed while retaining their complete source profiles.

This transfer applies only to child-uniform ports. It does not equate full named responses or nonuniform attachments inside a child. It can also be applied to a fixed larger quotient containing these children and additional singleton exterior names; it does not make a nonzero core contrast automatically survive ambient substitution.

## 3. Complete finite certificate for the directed triangle

Take Q=C3 with child order A->B->C->A. There are exactly5^3=125 canonical type triples and64 possible pairs of block masks. The verifier retains all64 actual source and opposite-quotient Hamilton counts for every canonical triple. It proves the following complete classification:

    All64 coarse contrasts vanish modulo8
        if and only if the triple is(E,E,E).              (2)

The proof of this finite assertion is exhaustive over the stated five-type domain, with literal subset Hamilton dynamic programming on the actual canonical tournaments. There is no isomorphism reduction, numerical tolerance, or sampling in(2). The complete canonical names, expanded masks and integer counts are retained. The independent audit must reproduce the values by a separate implementation before publication.

Among the64 type triples with all children nontrivial,63 have a contrast congruent to4 modulo8; only(E,E,E) is blind. Every nonzero contrast in these63 complete banks is4 modulo8. Singleton-containing triples also have a nonzero contrast, but residues2 and6 may occur and are not misreported as4.

Combining(1) and(2) proves the unbounded statement:

**Cyclic separation theorem.** For arbitrary nonempty children A,B,C, quotient converse is separated by a block-uniform physical two-port query modulo8 unless all three child orders are even and at least four. If all three are nontrivial and at least one is odd or is a pair, the selected actual contrast is4 modulo8. If all three orders are even and at least four, every block-uniform contrast is0 modulo8, regardless of their profiles or interiors.

The last clause is only modular blindness. It neither gives an exact coarse collision nor a full named-response collision. For example equal-profile even children are already known to have a nonzero exact canonical query by the independent positive-square theorem.

A useful control is C3[TT3,TT3,TT3]. Its familiar query(a,c)=(1,2) is zero modulo8, but its complete bank has12 residue-four queries, including(1,3). A single failed canonical probe was not evidence of whole-bank blindness.

## 4. The ambient boundary

If |A|+|B|+|C| is even, the two-port-augmented core also has even order. Give both ports the core's original incidence to every exterior vertex of any fixed rooted ambient tournament. The even-order transfer theorem in `../mod8_transfer/PROOF.md` multiplies the core contrast by an odd unit modulo8. Thus every nonzero separator supplied above survives every fixed ambient in this even-total-order case. A residue-four contrast remains4.

No automatic ambient lift is claimed for odd total core order. The inherited odd-order formula contains an additional Hamilton-cycle-parity correction. The finite-type classification itself is a whole-core theorem and does not remove that correction.

## 5. Scope and reproduction

The all-order content is the complete five-type normalization and exact multilinear run transfer. The cyclic classification is a finite exact certificate covering every canonical type triple, which then yields an unbounded child-order conclusion. It is not a census of all large tournaments, a theorem for arbitrary prime quotients, or unrestricted named-response reconstruction.

The verifier additionally tests the type reduction against exact complete profiles and Hamilton-cycle counts of actual tournaments, including all labeled tournaments through order five and fixed-seed larger controls. It retains a complete named bidegree-at-most-(2,2) bank for the all-TT3 cyclic control, including all cubic and quartic coefficients and zeros, with independent full-mask reconstruction checks. No truncation of the physical observer is used to obtain the coarse certificate.

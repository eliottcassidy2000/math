# Automorphisms and labelled enumeration at cycle count three

Session197. **All-order refinement of `START_THREE.md`.** This packet assumes the complete C3 bridge theorem and canonical start-three grammar from the preceding packets. No literature-priority claim is made.

Let Delta=1-3z+z^2. Let c_n count unlabelled strong tournaments with three Hamilton cycles, as in `START_THREE.md`. Then every such tournament has automorphism group of order one or three. The number s_n with automorphism group of order three has generating function

    S_C3(z)=sum s_n z^n = z^5(1-2z)/Delta^2.          (1)

Consequently the exact labelled count is

    n! (c_n - 2s_n/3).                              (2)

In particular the labelled counts at orders five and six are40 and6000. This is not the assertion that all C3 graphs are rigid.

## 1. Why only automorphism orders one and three occur

Every strong C3 graph has a degree-one vertex in at least one orientation. Its set of outdegree-one vertices, or of indegree-one vertices, has one or two members. Every automorphism fixes these members individually: a singleton is fixed, and a two-vertex set cannot be exchanged because its mutual arc would reverse.

Fix such a vertex v. An automorphism acts on the three Hamilton cycles. If it fixes a cycle, then fixing v fixes every vertex in cyclic order. Thus the action on those three cycles is free, and the group order divides three. A nontrivial group therefore has exactly three elements and is cyclic.

For any rooted strong start-three graph, the analogous statement follows even without the C3 theorem: an automorphism fixing the root acts freely on its three rooted Hamilton paths. Its order is one or three. Rooted start-one graphs have trivial automorphism group because their unique rooted Hamilton path fixes every vertex.

## 2. Symmetric rooted start-three inputs

Let R(z) count rooted strong start-three classes with automorphism group of order three. Use the canonical grammar and notation of `START_THREE.md`:

    K=z(1-z)/(1-2z),
    B_1=z(1-z)(1-2z)/Delta.

In Type I, the one-entry branch inherits the automorphisms of its rooted input. The mixed-pair branch uses a rigid C1 graph and contributes none. The three-entry branch has nontrivial symmetry only when its input is the triangle, with a nonempty tail. The latter contributes z^3(K-z).

In Type II, the rooted start-one input is rigid. The only possible nontrivial automorphism is rotation of the unique tail triangle. It survives precisely when the new root has uniform incidence to that triangle. If the triangle is terminal, only the all-zero incidence is allowed. If it is nonterminal, both all-zero and all-one incidences are allowed.

For a tail of total order k>=3, this gives (k-2)2^(k-3) choices. Its generating function is z^3/(1-2z)^2. The SCC order prevents an automorphism from exchanging the tail triangle with any other component. Therefore

    R=K R+z^3(K-z)+z B_1*z^3/(1-2z)^2,

and hence

    R(z)=z^5(2-4z+z^2)/Delta^2.                    (3)

This counts rooted isomorphism classes with stabilizer order three, not the number of eligible vertices in unrooted representatives.

## 3. Symmetry in the exceptional rooted subclass

Let R_E count the exceptional rooted class from Section4 of `START_THREE.md`, now requiring automorphism order three. Its one-entry branch propagates R_E. The unique qualifying mixed-pair input is the rigid four-vertex graph, so that branch gives no nontrivial rooted automorphism. The qualifying triangle triple-entry branch always preserves all triangle rotations. The triangle-in-tail branch is excluded by the exceptional condition.

Thus

    R_E=K R_E+z^3(K-z),
    R_E(z)=z^5/Delta.                              (4)

## 4. Marking degree-one vertices and forgetting the mark

For C3 graphs with a marked outdegree-one vertex, the canonical deletion cases are the same as in `START_THREE.md`. Nontrivial symmetry comes from either a symmetric rooted start-three input followed by singleton SCCs, or a rooted start-one input followed by one triangle and singleton SCCs. In the latter case the marked degree-one vertex loses to the entire tail, so every rotation of that tail triangle survives.

The marked symmetric count is therefore

    P_sym=z R/(1-z)+z^4 B_1/(1-z)^2.

Every automorphism fixes the one or two outdegree-one vertices individually, so marking does not change its order. The two-outdegree-one insertion/deletion bijection preserves the automorphism group: the inserted vertex is fixed and restriction gives precisely the automorphisms of the marked old graph. Consequently, as before,

    A_sym=(1-z)P_sym
         =z R+z^4 B_1/(1-z).

For the outdegree-only subclass, the symmetric marked input is exactly R_E. The same degree-one bijection preserves both the absence of indegree-one vertices and the automorphism group. Its unmarked symmetric count is z R_E. Converse gives an equally large indegree-only class. As the complete C3 theorem says every graph has one orientation, the full nonrigid count is

    S_C3=A_sym+z R_E.

Substituting (3),(4) and B_1 simplifies to (1).

A rigid class supplies n! distinct labelled graphs, and a class with three automorphisms supplies n!/3. This proves (2). The ordinary series with coefficients equal to labelled counts divided by n! is therefore

    z^5(1/3+(19/3)z-17z^2+11z^3-z^4)/Delta^2.

## 5. Verification and scope

The nonrigid unlabelled counts at orders five through twelve are

    1,4,13,40,120,354,1031,2972.

`verify_c3_automorphisms.py` rebuilds the complete low-cycle bank using the frozen C3 verifier. It explicitly finds the possible three-cycle automorphism from actual Hamilton words rooted at a fixed degree-one vertex, checks its action on every arc and its order, and separately counts fixed-root start-three stabilizers and their exceptional subclass. The displayed generating functions are checked by exact coefficient recurrences. Literal complete labelled checks through order six verify40 and6000 independently.

Normal and optimized Python runs produce the same deterministic certificate. The finite checks corroborate the canonical symmetry analysis. No automorphism classification is asserted for cycle count four or higher.

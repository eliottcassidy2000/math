# A recursive neutral-flip atlas with exactly controlled symmetry

Use the two order-eight graphs from [Session 310](../session310_positive_clocks_20261113/NEUTRAL_TAILS.md): the strong graph B and D=A→TT1. They are distinct, have the same full path-slot polynomial, and both are rigid. Their outgoing masks are

    B = (2,80,3,7,173,143,189,15)
    D = (130,156,153,209,225,207,135,0).

Rigidity follows by checking all permutations on eight vertices in the replay, without presuming a degree invariant is complete. Distinctness follows from strongness. Construct height d graphs on 8·3^d vertices by labelling each leaf of the rooted ternary tree by B or D, and using the directed triangle C3 as the quotient at every internal node. This is an explicit family of strong tournaments for d≥1. For each fixed d, all its tournaments have exactly the same full path-slot polynomial, not just the same Hamiltonian path count: when each global slot is split into consecutive module runs, all runs across all slots form one ordered spanning path cover of that module. The local cover profiles agree in the two leaves, so the possible reduced itineraries have the same filling multiplicities. Induct up the tree. This also proves substitution is legitimate in an arbitrary *uniform* exterior, not a selective rooted one.

## 1. Canonical interfaces and group recursion

In a cyclic three-block join X→Y→Z→X, with nonempty blocks of any sizes, a proper module cannot meet two distinct blocks. If a module meets X,Y and omits Z entirely, a vertex of Z sees opposite orientations to X and Y. If it meets all three, any exterior vertex in any block sees opposite orientations to its intersections with the other two; so there can be no exterior vertex. Thus the three blocks are exactly the maximal proper modules, and all isomorphisms must transport them. The induced permutation of the three blocks is a rotation. Any choice of the child isomorphisms over an allowed rotation lifts.

Write Γ_d for the group of recursive rotations of the underlying ternary tree. Then
Γ_0=1 and Γ_(d+1)=Γ_d^3 ⋊ C3, so |Γ_d|=3^((3^d−1)/2).
Because the leaves are rigid, each graph's exact automorphism group is the subgroup of Γ_d preserving its binary leaf colors. In particular the group is computed recursively: it is the direct product of the three child groups if their three classes are not all equal, and otherwise the wreath product of one child group with C3. This is a statement about the actual group, not just its order. Compare the inherited monochromatic triangle towers in [THM-4743](../../../01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md); here full path-count equivalence permits the colors to differ.

## 2. Exact class and symmetry inventories

Let A_d(z) have one monomial z^e for each distinct class at height d whose automorphism group has order 3^e. The exact recurrence is

    A_0(z)=2,
    A_(d+1)(z)=(A_d(z)^3−A_d(z^3))/3 + z A_d(z^3).

Indeed a nonconstant child triple has an orbit of size three, and its exponents add. A constant triple is a single ordered choice per child class and contributes exponent 3e+1. In particular the number of distinct equal-polynomial tournament classes is C_0=2, C_(d+1)=(C_d^3+2C_d)/3. The first counts are 2,4,24,4624. The number of rigid classes remains exactly two at all heights: each is the nonconstant rotation class of a triple of the two rigid classes below, and its three children must be rigid. At the other end, both monochromatic towers have the full group Γ_d. This also shows that rigidity at every order in this family coexists with an exponentially large full-group order at the same path-cover polynomial.

## 3. A single neutral flip can remove one level of symmetry at every height

Join two graph classes when they have leaf-color presentations differing at one leaf. This gives the quotient of the ordinary cube on 3^d coordinates by Γ_d, not a quotient by the full symmetric group on the leaves. Its graph distance between classes represented by words x,y is the minimum Hamming distance between x and g·y over g∈Γ_d. One inequality comes from projecting a Hamming geodesic. For the other, lift a quotient walk inductively: from a chosen representative of its current vertex, transport any one-step witness by a group element so its current endpoint is that representative. This gives a walk in the cube ending in the intended orbit. Its diameter is 3^d because the two monochromatic words differ at all leaves.

If one leaf is flipped, the exponent e changes by at most d. In the recursion e is the sum, over all internal vertices, of the indicator that its three child classes are identical. Only the d ancestors of that leaf can change indicators. The bound is sharp: flip one leaf of the monochromatic tree. Every ancestor of that leaf becomes nonconstant, since exactly one child has a changed color and the number of occurrences of B is an isomorphism invariant of a presented recursive tree. Thus a single neutral edit removes exactly d powers of 3 from the group order. Reversing the edit creates all d at once, while all Hamiltonian and uniform-context path-cover readouts stay fixed. Arbitrary marked-port stabilizers can be handled by retaining the corresponding fixed leaves; a priori one must not substitute across a selective port within a leaf.

## 4. Every intermediate group order is realized

In fact every exponent from 0 through (3^d−1)/2 occurs. Global color complementation pairs distinct classes with the same group exponent: it commutes with Γ_d and cannot preserve a class because the total number of leaves is odd and its number of B labels changes. Thus each available exponent has at least two available child classes. If all child exponents in [0,E] occur, choose three exponents summing to any integer in [0,3E]. Their representatives can always be taken not all isomorphic, even when the exponents agree. They give that integer as a group exponent upstairs. A monochromatic top layer above a child with maximal exponent gives 3E+1. The base interval is {0}. In particular the group orders have no gaps *in powers of three*, despite all full path-slot polynomials being identical.

## 5. Distances and stabilizers of named edits are computable recursively

For two height-d classes with child classes (x₀,x₁,x₂) and (y₀,y₁,y₂), their neutral-edit distance is

    δ_d(x,y) = min_(r=0,1,2) Σ_(i=0,1,2) δ_(d−1)(x_i,y_(i+r mod 3)),

with δ₀ the inequality indicator. This follows from independent rotations in the three children once the top rotation is fixed. The minimum can be evaluated bottom-up on pairs of named tree nodes: O(9^d) elementary comparisons suffice for leaf-by-leaf presentations, rather than enumerating Γ_d. The metric is not just distance between counts of B: at depth two, three B leaves placed in one child and three B leaves placed one in each child have equal weight, yet their distance is four by this recurrence.

If a set of leaves is marked pointwise, its stabilizer in the actual tournament is computed by the same recursion. At every ancestor of a mark no nontrivial quotient rotation survives (a rotation cannot fix the child containing that mark); in the other subtrees use the unmarked recurrence. Equivalently the order exponent is the number of internal monochromatic-by-class nodes whose descendant leaf sets contain no marked leaf. Marking a single leaf of an otherwise monochromatic height-d tree therefore removes exactly d powers of 3, the same loss as a neutral color flip at that leaf, although the marked and recolored structures are different. Forgetting the mark restores a symmetry without changing arcs; undoing the flip restores it while preserving all count observables. This distinguishes the two operations that a scalar automorphism count would conflate.

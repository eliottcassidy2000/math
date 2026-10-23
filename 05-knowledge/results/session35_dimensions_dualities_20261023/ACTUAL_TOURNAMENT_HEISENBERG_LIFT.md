# From actual57-vertex tournaments to the Heisenberg action

Session35. **PROVED elementary structural application**, with exact finite controls. This sharpens the sparse alternating-form construction without introducing artificial orientations on tied pairs.

## 1. Every odd tournament carries a symplectic augmentation space

Let M be the skew adjacency matrix of a tournament on odd n vertices. Every (n−1)-by-(n−1) principal submatrix is invertible over Q: modulo2 it is I+J of even order n−1. Since J²=(n−1)J=0 modulo2, (I+J)²=I. Thus rank M>=n−1; odd skew symmetry forces determinant zero, so rank M=n−1.

Let B have columns e_i−e_n for i<n, a basis of W={sum u_i=0}. Modulo2, BᵀMB=I+J of even order, again invertible. Thus the restricted alternating form uᵀMv is nondegenerate on W for EVERY odd tournament, without a regularity assumption. Every permutation automorphism preserves W and M, so acts symplectically on W. The Lie algebra W⊕Qz with [u,v]=(uᵀMv)z is Heisenberg of dimension n. This is an actual functor on tournament isomorphisms, not just a match of dimensions.

Regularity ensures M1=0 and therefore ker M=Q1; without regularity the null line need not be constant. Nevertheless the same augmentation space W is always transverse to that null line, as the restricted-matrix proof shows. A previous draft wrongly called regularity essential; independent audit strengthened the theorem before publication.

## 2. Exactly1024 Frobenius-invariant Cayley rules at57

Fix the cyclic point labelling Z57. A Cayley tournament is specified by S⊂Z57\{0} with S disjoint from−S and S∪(−S)=Z57\{0}, with i→j iff j−i∈S. It is regular of score28.

Ask also that multiplier7 preserve S. Multiplication by7 has two singleton orbits {19},{38}; the other54 residues form18 orbits of length3. Negation pairs these20 orbits into10 distinct pairs. Indeed an orbit cannot be its own negative: a group of odd order3 cannot carry a nonzero element to its additive negative in this odd cyclic group.

Thus choosing one orbit from each negative pair gives exactly2^10=1024 labelled Cayley rules, all invariant under C57⋊<7>. This is not a count of isomorphism classes; different rules may be isomorphic. The Singer difference-set geometry can coexist with this group action, but its incidence set does not uniquely select one of these orientations.

One explicit S is

    {1,2,3,4,5,6,7,9,10,11,12,13,14,17,18,19,20,21,
     25,26,27,28,33,34,35,41,42,49}.

Applying section1 gives1024 actual tournament-derived symplectic forms on the same augmentation module W, all supporting the Heisenberg57 action. By the separately proved order3 character obstruction, none of these distinguished C57 actions extends as the minuscule56 representation of E7. The obstruction is independent of which of the1024 orientations is chosen.

This is compatible with the fact that the abstract Heisenberg Lie algebra is the radical of E7½: forgetting the distinguished finite action removes the obstruction. It illustrates why an action-preserving map is more demanding than matching dimensions.

## 3. Inheritance and boundaries

The earlier prime-circulant work `THM-4725-no-uniform-dilation-rule-for-prime-circulants-signed-set-burnside-and-the-p-mod-4-dichotomy-with-the-paley-rule.md` motivates using signed residue-orbit sections. Its file is explicitly pre-audit; none of its stronger conclusions is a dependency here. The elementary orbit argument above is supplied completely and works at the composite modulus57.

The present1024-rule family is not the projective-line-flip family of size2^57. They are different operations on different parameter spaces. The projective-line parity theorem says a regular member's line-flip family has only itself and converse regular; it does not contradict many regular Cayley tournaments lying in different such families.

No claim about Hamiltonian-count distribution, uniqueness of an E7 embedding, or a canonical orientation from Blum/Jacobi data follows. The invariant form is constructed from a chosen actual tournament; the composite quadratic-residue sign by itself remains symmetric and insufficient to orient all pairs.

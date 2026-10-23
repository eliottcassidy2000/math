# Fibonacci's binary quotient and actual three-owner transport

Status: PROVED explicit representation dictionary and finite exact controls. The group isomorphisms are classical; their use here preserves actual owner identities. Original colored letters, additive charges, multiplicative states and electrical residues are not identified indiscriminately.

## 1. The exact six-map dictionary

Encode the three F3 owners by

0↦(1,0), 1↦(1,1), 2↦(0,1)

in F2² minus its zero vector. Translation x↦x+1 becomes

T=[[1,1],[1,0]],

and reflection x↦−x becomes

R=[[1,1],[0,1]].

Both assertions follow by applying the matrices to the three displayed vectors. They satisfy T³=I,R²=I,RTR=T^(-1). The six affine permutations of F3 therefore become exactly GL2(F2), whose six elements act faithfully on the three nonzero vectors.

In the field F4, choose omega²+omega+1=0 and use the basis(1,omega²). T is multiplication by omega; R is Frobenius z↦z². Thus this is also the semilinear group F4* semidirect Gal(F4/F2). It is a representation of the actual six permutations, not a presumed physical symmetry.

All these matrices have determinant1 over F2. Consequently determinant alone erases the reflection-versus-rotation distinction. Retain the full matrix or its permutation-sign lookup; a binary determinant is not the triple-cell chirality bit.

## 2. Parallel owners become an exact binary linear problem

Let a connected actual three-owner support carry bijections p_ij, with inverse maps on reversed edges. Encode each by M_ij above. Parallel owner assignments are exactly the NONZERO solutions of

v_j=M_ij v_i

over F2. A zero at any vertex propagates to every vertex, so every nonzero global solution is nonzero at every vertex. A root has only two binary coordinates, hence the solution space has dimension d≤2 and the number of owner sections is exactly2^d−1, namely0,1,or3.

This is a useful linearization, but zero must be excluded. A homogeneous system always has its zero solution; accepting that solution as an owner would falsely declare every attachment system compatible.

The holonomy group can be trivial, an order-two reflection group, C3, or S3. Its common fixed-owner counts are respectively3,1,0,0. Its full matrix commutant has F2 dimension4,2,2,1. After restricting invertible commutants to the allowed internal C3 rotations, the symmetry counts are3,1,3,1.

Thus the native order31 distinction between empty sections with residual C3 symmetry and empty sections with trivial symmetry is visible as commutant dimension2 versus1. Both cases have no nonzero fixed vector. A fixed-owner nullspace and an intertwiner/commutant space answer different questions.

## 3. Comparing systems retains the internal chart domain

If M and N are two edge-transport systems, chart maps X_i must satisfy

X_j M_ij=N_ij X_i.

Tree propagation expresses every X_i in terms of one root matrix. Remaining edges give simultaneous holonomy intertwining equations, linear in the four root entries. One must then retain invertibility, the allowed internal isomorphism set, all marks, and the transported internal restrictions at every cell.

For consistently oriented cyclic triples, internal maps are C3 rotations, not arbitrary GL2(F2) matrices. A reflection matrix represents a possible cross-cell selector but is not automatically an allowed internal chart change. For decorated cells, even the rotation domain may shrink. This dictionary does not remove those conditions or extend to arbitrary seven/five-cell regular blocks.

## 4. Where Fibonacci fits, and where it does not

The canonical three-letter substitution incidence matrix factors through the same T, as proved in COLORED_GRID_AND_PHASE_QUOTIENT.md. Modulo2 its quotient charge therefore has an exact order-three action. But original letters C and B are already identified by that quotient. The three nonzero charge vectors cannot be called the three original colors without specifying a different map.

There is no contradiction with the earlier obstruction that every additive homomorphism C3→(Z/2^s)^r is zero. Here the action is MULTIPLICATIVE in a matrix group, not additive in a binary vector space. Nor does this identify the electrical mod8 potential guard with Fibonacci transport.

For s≥1 the Fibonacci matrix has exact order3·2^(s−1) modulo2^s. Indeed T mod2 has order3, T³=I+2T, and T^6=5I+8T. Subsequent squaring increases the minimum entrywise2-adic valuation of T^(3·2^r)−I by exactly one, giving r+1. Thus the smallest sufficient power of two in the exponent is2^(s−1). These periods differ from the power-of-two Pascal path clock; an integer time remains a separate coordinate.

This supplies an exact parity/phase bridge. It supplies no map sending Fibonacci rows to Collatz trajectories, Goldbach prime pairs, or the tournament Hamiltonian-count gaps{7,21}.

## Reproduction

Run python3 check_binary_phase_dictionary.py and python3 -O check_binary_phase_dictionary.py. The checker verifies the full multiplication dictionary, every triangle and four-cycle six-permutation connection, all subgroups and their fixed/commutant spaces, and the dyadic clock through modulus4096. No floating-point operations or removable assertions are used.

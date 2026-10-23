# Actual tournament Heisenberg lift: independent audit

PASS after prepublication strengthening. The original regularity-essential sentence was incorrect; the corrected note now proves the stronger statement for every odd tournament.

## All-order proof

For any odd n tournament skew matrix M, every even principal submatrix reduces modulo2 to I+J, whose square is I because its size is even. Consequently its determinant is odd and nonzero in Q. Odd skew symmetry forces rank M=n-1 exactly. This includes n=1 with the empty0-by-0 determinant convention.

More directly, let B=(e_1-e_n,...,e_{n-1}-e_n). Since M reduces to I+J, each column of B has zero sum modulo2, and B^T B=I+J modulo2, the restricted matrix B^T M B also reduces to I+J of even size. Thus it has odd determinant. This proves nondegeneracy on the specific augmentation space WITHOUT regularity. It also proves the augmentation space meets the actual null line only in0. Regularity is needed only to identify that null line as the constants.

Tournament isomorphism permutation matrices preserve the coordinate-sum functional and transport the alternating form, so they induce Lie isomorphisms of W plus its central line. Composition and identities are preserved. This is a functor on the groupoid of odd tournaments and isomorphisms, not an assertion about arbitrary graph homomorphisms. The1-vertex boundary is the1-dimensional abelian Heisenberg algebra under the conventional zero-dimensional symplectic definition.

## Exact1024 scope

On Z/57, multiplication7 has order3. The nonzero fixed residues solve6a=0 and are19,38. The other54 residues form18 orbits of length3. Negation has no fixed nonzero residue and cannot stabilize a size1 or size3 orbit: applying an alleged multiplier-to-negative relation an odd number of times contradicts a=-a. Hence there are ten opposite orbit-pairs and exactly1024 independent choices of one member from each pair.

Every choice is a28-element connection set complementary to its negative and stable under multiplication7. It produces a regular tournament invariant under the faithful affine action C57 semidirect C3. The count is of labelled Cayley rules, not isomorphism classes. Distinct regular skew matrices give distinct forms on W, because their common kernel is the constants and the form on W determines the entire matrix under Q^57=Q1 direct-sum W. The previously audited order3 character obstruction is unchanged by choosing a different orientation in this family.

The inherited signed-orbit motivation is accurately described as motivation, not an unverified dependency. No unique orientation, new exceptional subgroup classification, or Hamiltonian-count distribution is asserted. The1024-rule space is correctly separated from projective-line-flip parameters.

## Independent tests and final pins

check_new_additions.py checks all1033 labelled odd tournaments at n=1,3,5 for odd restricted determinant, independently constructs and distinguishes all1024 connection sets, and verifies their Frobenius invariance. The root script was rerun; its JSON output matches the supplied receipt byte-for-byte.

SHA256:
- session35-root/ACTUAL_TOURNAMENT_HEISENBERG_LIFT.md: 8402d5d0495913a16286c15ea6bd0d6c743df65baba907eddcbabd5f2084414d
- session35-root/check_actual_tournament_lift.py: 75f35ab96ca9f6049c58703c43b556ecd9d59ee6e09def8e3dead303d14a3d99
- session35-root/actual_tournament_lift.json: c4b39e79d44775575d124e9bae6f1a958bc2824b82b99d56a046dce5fc17abd7

No Git operations or repository edits by this auditor.

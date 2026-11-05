# Independent review: four generalized pair families and finite inverse

**PASS.** No correction was required. This review covers the exact physical family Gamma and its finite numerical inverse, not a characterization of all C2 counts or an all-order coverage theorem.

## Frozen source

The reviewed directory is `session197/generalized-pair-inverse/`. Its exact source hashes are:

- `PROOF.md`: `9df4e06c88d28e159953601611bbcc6fb90f00f6766b5c09b03b3dfc7c76cbd3`.
- `core.py`: `dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042`.
- `pair_families.py`: `c298f019cadbb637e2665522031ea72feebbfb1b9c258d04ad49dc1b75451d18`.
- `inverse.py`: `9106978b1a21e3ce776c09bffeb062308dc159855572ae39d8a77360f253719f`.
- `verify.py`: `c9ac1cf71cc16b8a02ba64c6d3882881369b5002f8bbdf46f97f9cce006d1b1b`.
- `SUPPORT.json`: `872f629ccea09582f0a0be43965ead3478271a645d27b496a3dc1af4081586db`.
- Both source validation receipts: `0f778e2069a8e00a9ce9502326ad88d0bb3114e6abc5dd71dd2a47e92a108808`.

The four copied dependency hashes in the source manifest were also preserved. Source normal and optimized Python runs reproduce the frozen receipt byte for byte. Replays were performed on a copy; originals were not edited.

## Mathematical review

The insertion proof of `D-E>=d` is valid for every nonsingleton source with one starting path at the root: its root has an out-neighbor, so insertion before the first such vertex never places the root last. Deletion makes the map injective. The singleton is correctly handled separately.

The N seed is precisely the first transitive clone of the new sole-out-neighbor vertex, with the root kept at that first clone. The adjacent-clone laws give the displayed four coordinates, and deletion really leaves the original J graph. For M, the two forbidden endpoints in the root-ending count are b and z, with counts d and U. The reserved-cover split is exhaustive. In its final case, each old reserved cover has four placements of z,a; exactly the placement putting a after a reserved singleton b fails, giving the subtraction d. The sum is `3D-E+8U-d`. The stated two starting paths and cycle counts use the actual retained last-clone root.

The P cell and its distinct entry/exit agree with the already reviewed physical matrix. The following J move gives the stated new root and deletion response. The singleton is legal only for P. The O/Z suffix law keeps its deletion sidecar; the empty suffix at parameter1 is not confused with an actual birth. N's nonempty suffix ends in Z. Root-start count is preserved even through weak intermediate O graphs, and the final cycle bijection gives exactly two cycles. No unwarranted strong-induced-cycle monotonicity is used.

All four scalar formulas follow from those physical rows. The minimum induced graph argument proves the `D<=h/7` cutoff. Weak sources have E=0, giving the additional U bound. The four displayed parameter limits use valid lower bounds; each denominator is positive in its permitted domain. In particular, Q's potentially subtracted `U-d` is nonnegative.

The source search is finite for a proved reason, not because a numerical queue stabilizes: accepted quadruples have `h,e,d<=B` and `R<=4B`. Every move to a final source below B remains below B by induced path monotonicity, and `L R<=h_new` supplies the exact gadget bound. The starting singleton and its unique count-one two-vertex transient are retained correctly. Weak sources have a strong start-one initial core followed only by singleton SCCs. The tail doubles U and preserves d, so distinct erased ending responses can leave distinct physical deletion sidecars. The two weak states `(9,0,18,1)` and `(9,0,18,3)` are indeed both required.

The algorithm enumerates every relevant quadruple and every allowed final parameter for Gamma. Deduplicating equal quadruples is legitimate for its scalar formulas and does not pretend to enumerate all graph isomorphism classes. Complexity is polynomial in the numerical target, not its bit length. The claimed bounded support through1861 is all-order for that numerical interval because the cutoffs are proved; it is not a cofinite extrapolation.

## Independent arithmetic and physical reconstruction

`generalized-pair-inverse/check.cpp` is standalone, with no imports of source code. SHA256:

    cebc113c208d8525e0006d4033d0b6135f7dcffc3169592838978a5ea924038c

It derives the gadget's H by assigning transitive-deletion vertices to the part before or after the root. This literal path-partition counter does not use the source's g recurrence. It derives Q by choosing the first base vertex of the reserved root-start path, and uses the physical transitive-deletion value for B. It then separately regenerates the complete source quadruple banks and the four-family scalar image.

The image scan deliberately uses the larger `v<=target` rectangle rather than the four parameter-bound implementations. This is safe because the physical source inequalities imply `F(v)>=v` for each family. Thus the exact support comparison also tests the tighter parameter pruning in the source.

Every small seed response is computed by explicit permutations and the unique cut immediately before the marked root. This independently counts H, the endpoint counts, the reserved two-cover response, deletion H and Hamilton cycles. Larger support witnesses are rebuilt by separate adjacency constructors and counted by a literal directed-word traversal with no subset memoization. The independent run makes48,976,881 recursive word visits.

Results:

- 403 complete literal seed profiles, including their marked roots and cycle counts;
- 271 further small literal family outputs;
- exactly1,316 source quadruples through265;
- exactly2,000 relevant source quadruples after weak tails;
- exactly23,432 C1 source quadruples through1861;
- every one of the883 support values reconstructed as an actual graph and independently counted with H equal to its target and C=2;
- maximum witness order14;
- explicit physical validation of both weak deletion-sidecar collisions and their distinct N outputs133 and129;
- exact equality of the complete source banks, not only their sizes, and exact equality of the883-value support;
- union with the C1 image misses exactly7,21,67 in the stated interval.

The optimized checker and a separate library-assertion/undefined-behavior-sanitized build produce identical receipts. Both have SHA256:

    54ee1545a81ef632b1fd1da2d9c4a98c2451f4289c9bcfac27d2bb37bd193ce4

Their complete serialized bank files are also identical, SHA256:

    094d51c190bef66c91f0658d87d9b91a82f6d3b51ceb4bd3e7889c92b704a5b1

`compare.py` imports the source only after independent generation, for the final equality check. Its SHA256 is `b93abee8d542b4f4f305ec9d455b7a0134e171ced39399b1110790a9b13c1076`. The comparison receipt has SHA256 `2cd88d55331ecd6deedbbfde2a392515d9f0d0e1a521b0d36419667e34b13ba7`.

Reproduce with the frozen packet directory supplied as SOURCE:

```sh
g++ -std=c++17 -O3 check.cpp -o /tmp/check-generalized-pair
/tmp/check-generalized-pair /tmp/independent-banks.txt > /tmp/independent.json
python compare.py SOURCE /tmp/independent-banks.txt /tmp/independent.json
```

Replace the compilation flags by `-O1 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined` for the independent debug replay. The complete bank text is generated by the checker and need not be shipped as an additional input.

The exact inverse is a useful finite constructive reduction. It still does not show that arbitrary sequences after a pair gate can be absorbed into Gamma, nor that every admissible odd target has a Gamma or C1 presentation.

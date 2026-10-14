# One affine attachment: exact overlap relations and the first coupling failure

**Date:** 2026-10-14. **Status: FINITE-EXACT. Novelty unconfirmed.** No higher chromatic bound or global affine-colouring theorem is claimed.

## Outcome

The proposed affine operation really breaks the old glued colouring formula, but its first two exact instances remain four-colourable.

- `Y_1`: 140 vertices, 367 edges; its six-vertex overlap has eight proper four-colourings up to global colour permutation, and **all eight extend** to the complete graph.
- `Y_2`: 1,182 vertices, 3,653 edges; its 20-vertex overlap admits the canonical triangular three-colouring, which extends even with **all 30 base-lattice vertices fixed** to that colouring.
- The simple common-centre formula, glued with the palette permutation needed to agree on the overlap, fails on exactly 12 edges for `Y_1` and 80 for `Y_2`. SAT-discovered repairs exist and are retained as explicit colour vectors.

These are complete finite statements for the specified graphs. In particular, no claim that every proper colouring of the 20-vertex overlap extends is made.

## 1. Exact construction and inherited theorem

Use the rotations

`R_1: cos(theta)=5/6, sin(theta)=sqrt(11)/6`,

`R_2: cos(theta)=23/24, sin(theta)=sqrt(47)/24`,

and the triangular basis `e_1=(1,0)`, `e_2=(1/2,sqrt(3)/2)`. Let `G_m` be the complete unit-distance graph on the distinct points

`{R_1^a R_2^b (i e_1+j e_2): |a|,|b|,|i|,|j|<=m}`.

The [independent-hinge seam theorem](independent_hinge_seam_four_colouring_20261014.md) proves `chi(G_m)=4` for every `m>=1`, with a global formula on the entire infinite common-centre union. Define

`Y_m=G_m union (e_1+G_m)`,

identifying equal points and adding **every** unit edge, including edges between the two copies. The overlap includes the unit edge from 0 to `e_1`. Each `Y_m` contains `G_m`, so the displayed four-colour certificates prove `chi(Y_1)=chi(Y_2)=4` exactly.

Coordinates lie in `Q(sqrt(3),sqrt(11),sqrt(47))`. The checker uses the eight-component multiquadratic basis, rational coordinates and exact integer arithmetic. It tests every unordered pair: 9,730 for `Y_1`, and 697,971 for `Y_2`. No float-distance filter is used.

## 2. The overlap relation, rather than a chromatic number alone

For a graph `Y` with labelled overlap `I`, define its four-colour boundary relation

`R(Y,I)={tau:I->{0,1,2,3}: tau extends to a proper colouring of Y}`.

This retains the actual colours on the shared vertices. For `Y_1`, the checker enumerates every proper colouring of the induced six-vertex overlap after fixing one triangle to colours 0,1,2. Every proper four-colouring can be normalized this way, and the fourth colour is then fixed automatically; thus no colour-permutation orbit is lost.

Exactly eight patterns remain. The certificate file supplies one full 140-vertex extension for each, and the checker verifies every edge and the restriction to `I`. Consequently

`R(Y_1,I)=the full proper four-colouring relation of the induced overlap`.

For `Y_2`, the retained certificate fixes the base triangular coordinates to

`colour(i e_1+j e_2)=i-j mod 3`.

There are 30 distinct such base vertices across the two patches, including the 20-vertex overlap. The checker verifies all 30 conditions and every one of the 3,653 edges. This is one structured extension, not a census of all boundary patterns.

## 3. Genuine interaction exists, but it is repairable here

Apply the earlier formula to the first copy. On the translated copy, apply it in coordinates relative to `e_1`, then change the palette by

`0->1, 1->2, 2->0, 3->3`.

This is the palette permutation that aligns the two canonical triangular three-colourings on their overlap. It yields a well-defined colour at every shared point. Nevertheless, 12 or 80 unit edges become monochromatic, respectively. All are fourth-colour conflicts in the recorded instances.

Thus the affine operation is mathematically different from adding another independent common-centre hinge. It produces real constraints on simultaneous colour choices. The explicit extensions demonstrate that those constraints can still be repaired at these sizes, including without changing the prescribed base-lattice colouring in `Y_2`.

## 4. The next criterion, and the stopping boundary

A useful interface should forbid at least one boundary pattern that was legal before attachment. Merely increasing edges or making one convenient colouring formula fail is insufficient.

A precise next target is a geometrically realized overlap `I` and a proper assignment `tau` that extends to each component separately but fails to extend to their complete unit-distance union. Record that assignment and its nonextension certificate. Such a strict boundary restriction can then be composed with other restrictions. A five-chromatic union would require every possible boundary assignment to fail; a single forbidden pattern is an intermediate mechanism, not that conclusion.

`Y_1` creates no such restriction on its chosen overlap, because all eight proper patterns extend. The `Y_2` test preserves the canonical triangular pattern but does not classify its whole relation. These facts justify stopping this bounded pass rather than simply raising `m`. A different rigid interface or affine attachment should be chosen for a concrete expected boundary restriction.

## 5. Reproduction

The package contains:

- `04-computation/experiments/affine_hinge_overlap_check_20261014.py`
- `05-knowledge/results/affine_hinge_overlap_certificates_20261014.json`
- `05-knowledge/results/affine_hinge_overlap_check_20261014.json`
- this note

The checker imports the already-published exact field and vertex routines from `04-computation/experiments/independent_hinge_seam_four_colouring_20261014.py`. Both scripts use only the Python standard library. SAT was used to discover the colourings, not to verify the retained claims.

From repository root run:

`python3 04-computation/experiments/affine_hinge_overlap_check_20261014.py`

`python3 -O 04-computation/experiments/affine_hinge_overlap_check_20261014.py`

Both outputs match the declared result JSON byte-for-byte. The checker regenerates the exact graphs, independently enumerates all eight normalized overlap patterns for `Y_1`, validates every extension, and validates the 30 fixed base colours for `Y_2`. No graph beyond `m=2` was searched in this pass.

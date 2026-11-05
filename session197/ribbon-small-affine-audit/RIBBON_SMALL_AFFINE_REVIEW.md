# Independent review: span-three small affine clocks

**PASS**, no corrections requested, 2026-11-05. This review independently rebuilds both the physical long-word inventory and all 55,057 matrix certificates; it is not only a replay of the author's checker.

## Frozen source

Path `session197/ribbon-small-affine/`:

- PROOF.md: 4a512b41672971d365f170a84f34938e9ce6401a6196f64aff6912ea8e721a2b
- verify.py: fce80af070f17e18ccc2dd1b94c95be440b3ff11c7796833f53359dbd84c2d22
- verify.cpp: 815b55dd56aad6ff5ad7edca6901803a73a92747fefb880660cf381236364fc0
- VALIDATION.json: 576b62554ff9cba4a59ca2b840461c6a68332417bedabad41282b07148654a3b
- verify_physical.py: c9aa5709cd60e8eaf0c0bd281253a60d162cd1db2ee8c0c9274a8c08eda7bffc
- PHYSICAL_VALIDATION.json: b64df141bde9f6571e27e72f2be1b883192a333d81a9a7695393a5e65d0ce6af

Dependencies already reviewed in the preceding ribbon audit:

- ribbon/transfer.py: 2ef89174390c503494d4e88e81f1b7ee149f985b6746f2912430ae88e46a289d
- ribbon-periodic/representation.json: d418468ec6869b2eab56a8406b48b19e78ac64993da686fbe767e5e264cd8b5a

The supplied checkers were run in a copied snapshot. Both generated validation files reproduce the source hashes exactly. Python assertions were enabled. No optimized-Python validation is claimed for these assertion-based checkers.

## All-order proof review

1. The low-count cell classification is valid: every strong tournament of order at least four has an induced strong four-tournament, all such four-tournaments have five Hamilton paths, and induced-path insertion is injective. Two triangle SCCs already give nine paths. Thus a cell with H<=4 consists of singleton SCCs and at most one triangle.
2. In the absence of a path from the last boundary to the first, all SCCs meeting the first boundary precede those meeting the last. The I/M/J split is intrinsic to the SCC order, even when its vertex sets are not natural-order intervals. Every backward intercell edge is confined to J_j union I_(j+1); every prefix-to-cell backward edge is confined to prefix union I_1. This proves an actual ordered-join factorization and exact geometric growth, not a response-coordinate assumption.
3. Disjoint-cell path merging gives H(PW^k)>=H(T(W))^k. The same holds for rotated phases after deleting only bounded end fragments. Therefore the growth and rotation filters are valid for arbitrary fixed prefixes.
4. For width three and H(cell)<=4, a nontrivial SCC is a triangle of index diameter at most six. At length twelve the first/last-three separation is at least seven, so the boundary factorization excludes all remaining long cells. The finite length-eleven cutoff is therefore rigorous.
5. The general period reduction is also correct. Every strong component S has at least |S| Hamilton paths by cutting a Hamilton cycle. If H(cell)<=b, its components have size at most b and index diameter at most w(b-1), by following a simple path from maximum to minimum index. At block length at least w(b+1), boundary separation is strictly larger. Otherwise H(cell)>b contradicts affine growth. This proves the stated necessary bound |W|<w(b+1), not a realizability criterion.
6. An eventual nonconstant recurrence H_(k+1)=mH_k+c with c!=0 has two nonzero exponential components, at m and 1. The characteristic-polynomial recurrence of any exact linear representation forces both eigenvalues. This necessary-condition argument remains valid after arbitrary initial transients.
7. All initial padding bits are retained. Rotation invariance of the characteristic polynomial of a square-matrix product justifies cyclic representatives, but does not identify or discard different seam bits.

## Independent physical enumeration

`small-affine-audit/triangle_inventory.cpp` is a new graph-only enumerator. It does not import the 288-state automaton or the frontier transfer. It uses the independently justified equivalence

    H(T)<=3 iff T contains at most one directed triangle.

Indeed any strong component larger than three contains a strong four-tournament and hence at least two triangles; two nontrivial SCCs also supply two triangles. Conversely one triangle and singleton SCCs give H=3.

The enumerator counts new literal triangles as each vertex is added, pruning only when two are present. It tests rotated cells directly by graph triangles. Boundary reachability is checked from the one-triangle SCC structure and independently cross-checked by full Floyd closure on every accepted word through order eight and every boundary-active larger word.

It exactly reproduces all five rows (length, low words, boundary-active words, final rotation classes):

- (7,62976,29952,417)
- (8,199616,39872,1091)
- (9,604160,29440,921)
- (10,1768448,24576,675)
- (11,5147392,9216,329)

The resulting 3,433-word candidate file has SHA256 55979dd3b3583f8bfc7beb745d1873f516f71f359cb9e8314ca1d251b714735f.

## Independent matrix certificates

`small-affine-audit/independent_determinants.py` independently generates all 51,624 cyclic classes of lengths one through six, and checks these together with the 3,433 independently generated longer classes. It uses dense integer matrix products and a separately implemented modular Gaussian elimination, rather than the source's sparse C++ multiplication and elimination. The audit uses the environment's already installed NumPy, with no package installation. Integer magnitudes remain far within int64 range.

All 55,057 classes are resolved for both multipliers three and four. Every per-length, per-prime first-success histogram agrees exactly with the source validation. No floating-point arithmetic or approximate eigenvalue test is used. In the source C++ products, the largest representation coefficient is 13 and the maximum column absolute sum is 147; with residues at most 12, the multiply accumulator magnitude is at most 1764, so the source's signed integer arithmetic also has no overflow concern.

The independent result is saved in `small-affine-audit/independent-determinants.json`. The source physical replay additionally reproduces 4,096 order-six reachability checks, 600 randomized reachability checks, and 437 literal SCC factorizations.

## Scope

The theorem excludes fixed-block, nonzero-offset affine clocks with multipliers three and four at span three, at every block length and after arbitrary fixed prefixes. It does not exclude state-dependent words, larger finite control cycles, changes to old edges, or other synthesis mechanisms. The general period proposition only reduces a search; an eigenvalue pair alone does not establish a usable prefix or clock.

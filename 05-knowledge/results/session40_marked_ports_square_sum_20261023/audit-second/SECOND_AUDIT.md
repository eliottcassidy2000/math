# Session40 second independent audit

**Status: PASS.** This packet independently audits the marked-tournament decoder, full order-eight endpoint signatures, square-sum census/edge certificates and minimum banks, conditional path scatter, and the new five-cell boundary. It makes no external priority claim.

## Marked quotient and shared chirality

Every stable directed equitable cell in a tournament has odd size. With two distinguished vertices at n≤8, the only unresolved multi-cell interaction is two cyclic triples. A regular five-vertex cell is unique; external singleton incidences are uniform. The two triple cells have either uniform or permutation/complement blocks. Internal rotations classify a nonuniform block by one permutation parity bit. The shared four-vertex statistic J=±3 evaluates this bit, distinguished modulo8. Its normalized matching and sign are preserved under converse. This proves marked target-class decoding, not source-only prediction or literal partner reconstruction.

The proof extends when r distinct marks leave n−r≤6 vertices; when n−r≤5 no extra bit is necessary. Seven unmarked vertices defeat the same observer: cyclic-seven modules with Hamiltonian counts175 and189 survive the same marked quotient. The coprime-cell lemma a q=b p correctly forces a uniform block when gcd(a,b)=1. Thus regular-five cells can be included with cyclic triples provided every five-to-five block is explicitly uniform.

Independent structural controls enumerate all1,024 tournaments on five vertices (24 regular presentations, one isomorphism type), all512 binary3×3 blocks (14 regular blocks), all mixed chirality/converse cases, and all5,040 vertex orders for each seven-cell hostile.

Independent full Python enumeration uses backward subset recursion and exact tuple refinement histories. It matches all6,880 representative classes,296,534 Hamiltonian-feasible ordered endpoint pairs,291,810 unaugmented signatures, all four residual-cell shapes, and the exact18 cross-class collisions. Evaluating J on all48 mixed two-triple cases gives291,828 augmented signatures and zero cross-class collisions. Both normal and optimized Python runs pass.

The producer's full literal replay was separately rerun after source inspection: all2^21 fixed-path masks decode to the certified class ledger, all2^21 reflection/converse checks pass, and the complete956,339 ordinary-port and460,699 merged-edge ledgers agree exactly. Its frozen receipt is reproduced in full-mask-replay.out. This rerun is not called a second independent full-mask implementation; the independent Python endpoint scan and structural proof supply the separate methods.

## Square-sum path and edge data

The independent edge-subset method maintains acyclicity and degree≤2; accepted n−1-edge subsets are exactly spanning paths. It traverses1,326,143 search states through n=25. Every path total, edge participation and endpoint count matches the producer's distinct path-DFS census:359 comparisons. A zero-path order is not called a forced-edge order.

The certificate excluding edge1–15 at order25 passes all24 possible endpoint pairs containing the leaf18. All215 recorded degree deletions/forces are necessary consequences of their actual endpoint targets, and every final deficit/excess/cycle contradiction is valid. No branch beyond endpoint selection is used, and the retracted rule that every graph-degree-two vertex is internal is not used.

Minimum global edge-membership banks with the actual endpoint pair supplied at n25,...,32 are0,1,3,4,1,1,5,6. Independent brute combinations, not the producer's hitting-set recursion, certify lower bounds. At32 all169,911 combinations of five distinct test columns fail; at31 all8,855 combinations of four fail. These are not arbitrary encoding or adaptive-query minima.

The parity-tree decoder is exact on the entire affine binary boundary fibre. An unseen forest has no nonzero cycle kernel; an unseen cycle can be toggled. This larger-fibre necessity is not a lower bound on the restricted Hamiltonian path family, because a toggle can destroy connectedness or degree bounds. The edge-marked birth polynomial also has the correct component-owner condition: joining endpoints must lie in different spanning path components. Singleton components are counted once.

## Conditional covariance and the source-only hostile

For a finite positive Gibbs ensemble, the log-partition Hessian is the edge-indicator covariance. Its kernel is exactly the linear costs constant on that ensemble. Positive weighting does not change this kernel. After actual endpoint conditioning, degree-potential costs are conserved, but they need not span all conservation laws. Exact SymPy DomainMatrix ranks independently match all listed square-sum orders. At32 the global rank45, within-endpoint rank12, degree-potential rank32 and extended conservation rank34 are all confirmed. The two additional costs are the unit edges5–31 and9–27. This finite statement is not an all-order charge law. Linear variation rank and number of binary distinguishing queries are different quantities.

The five-vertex source mask4 hostile is independently verified: source H=5; paths(4,0,3,2,1) and(4,2,0,3,1) share actual endpoints but yield non-path-complement targets871 and733 with H=15 and11. Converse preserves H, so merged targets differ. All75 labelled tournaments through order4 have at most one path per ordered endpoint pair. This separates the actual-target-marked decoder from an unjustified source-plus-endpoints prediction rule.

## Five-by-five next dyadic boundary

Independent enumeration of all120 matching permutations under C5×C5 rotations gives eight ordinary classes: four orbits of size5 and four of size25. Converse acts by p↦(−id)∘p⁻¹∘(−id), producing seven merged classes: six self-converse classes and the slope2/slope3 pair. A parity bit has only two labels and cannot classify them.

The shared orientation sum J on the eight lexicographic representatives is10,4,−2,−4,2,0,0,−10. Exact J and J modulo16 classify exactly the seven merged classes. J modulo8 has only four labels and folds three pairs. Thus16 is the first sufficient power-of-two modulus in this specific family; it is not the smallest modulus among all integers (modulo9 also separates these values). The ordinary slope2/slope3 pair remains indistinguishable by J, as it must for a converse-even observer.

Independent vertex-path DFS gives H=8265,8245,8221,8201,8193,8175,8175,8245, matching the producer's subset DP. No general five-cell synchronization theorem or replacement of the native electrical modulo8 guard follows.

## Reproduction

Run each included Python script normally and with python3 -O. They retain explicit failure checks. check_path_scatter.py requires SymPy; the others use the standard library. The endpoint scan needs the committed Session38 wiring class file and the adjacent Session40 marked-ports endpoint_scan.txt. Square consumers need the adjacent square-sum/receipt.json. Run the edge-subset producer before its comparison consumer. No scripts use or modify Git.

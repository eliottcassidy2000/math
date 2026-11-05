# Independent audit: bridge-C2 classification and enumeration

**PASS**, no correction, bound to the exact proof, verifier, receipt and pinned start-two dependency hashes in SOURCE_SHA256SUMS. Frozen source files were copied before replay and remain unchanged.

## Mathematical review

The indegree-one lemma is valid. If a start-one root r has unique in-neighbor v and T−r were strong, every vertex other than v supplies a positive initial-path contribution, giving at least n−2 starts, a contradiction for n>=4. Otherwise v must be in the terminal SCC; r dominates the entire initial SCC, so the root-start identity equals H(T−r)=1. This deletion is transitive with v last. The displayed two distinct paths from its first vertex establish an indegree-one vertex outside the start-one set. For a triangle, any selected pair leaves the third indegree-one vertex.

The complete start-two grammar then inductively preserves an indegree-one vertex outside the entry set, and in particular different from the current root. In the one-entry case use the smaller start-two input; in the two-entry case use the C1 lemma. Incoming arcs from the new root hit only the selected entry vertices, and the tail sends no arcs into the old component.

The five descriptions in Theorem3 are equivalent. The bridge factorization has two positive endpoint factors whose product is2, so one is1. The endpoint-one structural lemma applies even if that side is reducible, by isolating its terminal component and the singleton components preceding it. This produces a degree-one vertex in the full graph, not merely within the side. Reversal supplies the other direction. Conversely, a degree-one vertex exposes a strong bridge and cuts/closes exactly two rooted Hamilton paths.

Rigidity follows because an automorphism fixes an outdegree-one vertex and acts faithfully on its two starting paths; a nontrivial permutation would induce an involution, impossible for a tournament automorphism. The two-outdegree-one deletion bijection uses only strongness and the degree constraints, so it preserves arbitrary positive cycle count and is valid for C2. Thus c_n=sum_(m<n)b_(m,2), and the unrooted bridge count is c_n−c_(n−1)=b_(n−1,2). The generating function and labelled conversion are correct.

As written, the frozen packet's total-count conclusion starts at order74 by the separately reviewed residual73 bound. That statement is correct. Any later all-order removal of the residual atom polynomial is a separate result, not silently included in this signoff.

## Independent checks

Normal and Python -O source replays match the frozen certificate byte for byte.

The separate independent.py reconstructs all33,867 labelled tournaments through order6, using subset start counts, literal directed cycles, and independent SCC calculations. All13,584 start-two rooted instances have an indegree-one vertex different from the root. Among all3,744 labelled C2 tournaments,55,440 explicit single-arc deletions verify bridge/degree-one/start-two equivalence. The only bridge-free base class is the regular five-tournament, identified by full permutation canonicalization. There are120 labelled bridge-C2 graphs at order5 and3,600 at order6.

For the subsequent structural investigation, the same script checks139,968 deterministic Hamilton-cycle insertions into strong single-vertex extensions, with membership and injectivity verified literally. It also checks the degree-two covering inequality on504 small 2-strong examples where every deletion has an outdegree-one vertex. Those are useful controls, not a substitute for an all-order proof. See INDEPENDENT.json.

No complete unrestricted Hamilton-path spectrum conclusion follows from this packet. The approved classification concerns the strong bridge-C2 class, with the stated high-order total enumeration.

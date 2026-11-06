# Independent Session203 audit

Outcome: PASS after the author corrected the weak-case parity sentence in PROOF.md section 5. The sentence now correctly says every Hamilton-count factor in the tail product is odd, rather than incorrectly calling the start-two root factor odd. No mathematical claim or implementation needed changing.

## Proof audit

The one-port squarefree-polynomial lemma is valid. An A-run cannot follow a context vertex other than v, and v occurs once, so only an optional initial A-prefix per bin and at most one exceptional A-run after v survive. A bin beginning in v changes f to E; a nonempty post-v run contributes S-1. For the root-starting designated output bin there is no ordinary initial prefix, leaving exactly k optional ordinary prefixes. A root-ending designated bin has an optional prefix like the other bins. Deletion of r is an ordinary context count, including the case that its remaining context is only v. The resulting four products span the four input coordinates, with C_k divided by the nonzero integer k where needed. Empty A is legitimate. Thus determining a transfer on four independent physical inputs is an all-input linear-algebra identity, not finite extrapolation.

The literal gates and their matrices agree. K1 has exactly the stated singleton start exception; its other eight coordinate formulas still hold. K2 on a singleton has start two and its rooted graph is isomorphic to K1 on the edge. Neither singleton needs adding to the source grammar.

The canonical decomposition is used in the correct direction. Repeated one-entry layers reduce order; the final two-entry layer lies in a unique-cycle initial component. Consecutive degree-one entries can be deleted because a v-starting Hamilton path must leave their degree-one predecessor(s) for its terminal segment. Conversely that segment appends uniquely. In the triangle choose the consecutive orientation of the unordered pair and use K1 on an edge. The skipped pair is needed only beyond the triangle, so its residual U is nontrivial. All arbitrary transitive tails are exactly O/Z words retaining the marked root. Weak graphs add only a singleton tail, after the corrected parity argument.

The growth inequalities and gate differences are correct. The K2 ordinary difference is exactly 2(h-s)+6(R-d-s)+7s+4d-2e. Its cover difference 3C-3L+4T-2b is nonnegative. The final-operation estimates imply the stated max(2,3+floor(log2 M)) cap after removing the terminal O-run. This cap is for represented (B,s) pairs; it need not preserve all nine coordinates under removal of O. That distinction is respected by the target-membership theorem.

## Implementation audit

The entire nine-coordinate tuple, including s, is the state key. Equal keys have equal permitted successor sets and equal future profiles. A lower-order representative weakly dominates a higher-order representative because all operations have fixed positive vertex costs. The priority queue correctly implements that minimization; no equal-letter-cost assumption enters. Witness words are retained for the selected physical representative. A resource stop is reported INCOMPLETE and never absence. Exact exhausted absence is restricted to the start-two B fiber, with no Hamilton-spectrum or staggered-F exclusion.

## Independent controls

Run `python3 session203-start-two-gates-audit/audit.py` from the workspace root. The audit uses explicit permutation enumeration of paths on every subset, independently of the author's dynamic-programming path counter.

- 150 literal gate checks on every root-zero labeled input through order four, including singleton corrections.
- 192 arbitrary-context one-port checks, using every three-vertex context, every three-vertex input, and cover levels k=1,2,3. Context matrices are reconstructed from the four physical basis graphs.
- Unmerged traversal of all 957 retained recipe prefixes for B<=139 matches the compiler's 478 merged states and its entire start-two B support.
- Normal and Python -O independent receipts agree byte-for-byte.

Files: audit.py, RESULT.json, OPTIMIZED.json, AUDIT.md in this audit directory. The author's NORMAL.json was inspected and reports PASS; author replay is separate from these independent checks. No author files or shared repository files were edited. The existing dirty math checkout was read only and preserved.

Scope: the gate maps, complete start-two normal form, target cap, and guarded compiler pass this audit relative to the stated inherited structural theorems. F(8,137), all-base occupancy, and an eventual unconditional odd-spectrum threshold remain OPEN.

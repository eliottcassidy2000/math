# Independent audit of the hidden holonomy storage packet

**Verdict: PASS. No mathematical defect found.**

The auditor read NOTE.md, check.py and certificate.json and reproduced both normal and optimized execution. Both outputs match the published certificate byte-for-byte, and both documented SHA-256 values match.

## Proof checks

- The ordered strong-component decomposition really gives injectivity on unlabelled isomorphism classes; no rigidity of individual seed labels is assumed.
- The metric theorem is explicitly labelled arc distance, so possible isomorphism-dependent shortcuts do not contradict it.
- Constant slot profiles multiply under ordinal sum by a direct restriction/concatenation bijection, including empty slots. Hamilton multiplication follows by setting the number of slots to one.
- The decoding criterion 2e+s<delta is valid for errors on arbitrary arcs, not merely designated data arcs. Fixed arcs add redundancy and cannot reduce code distance.
- The local read head requires the selected block and its local vertex 0. This marking is openly stated. The same rule is used on both seed types, so the query does not presuppose knowing the answer.
- The K-query optimality is valid for this binary local query family. The aggregate-query counterexample correctly prevents extending that bound to unrestricted integer-valued Hamilton observations.
- Raw adjacency decoding and noise in Hamilton outputs are correctly distinguished.
- The final scheduled-circuit origin is explicitly inherited; no intrinsic autonomous implementation is asserted.

## Extra independent controls

The added `audit_independent.py` imports no original checker. It computes path covers using set partitions and independent brute-force Hamilton counts on their blocks, rather than the original permutation-and-cut cover enumeration. It verifies both complete seed profiles and thirty evaluations of the displayed slot polynomial.

For each six-vertex extension it enumerates all Hamilton vertex orders grouped by start vertex, reproducing the exact rows (1,2,4,8,1,1) and (1,2,3,7,5,1), hence totals 17 and 19.

It constructs the sixteen affine codewords on F2^3 by parity-dot evaluation independently of the original monomial-basis generator. It checks all 2,512 combinations of codeword, erasure set and visible error pattern with 2e+s<4. Each has the original word as its sole candidate. It also gives the sharp equality-boundary ambiguity by erasing the four disagreements of a minimum-distance pair.

Its normal and optimized outputs agree. Files: `audit_independent.py`, `audit_independent.json`, `audit_independent.optimized.json`.

A nonblocking checker-strengthening opportunity: the original single-arc loop establishes distance(original,received)=1 and uniqueness of the nearest codeword; adding an explicit check that the minimum equals 1 would make that local test assert the claimed winner directly. This is not a correctness gap in the combined verifier or proof, since minimum code distance four has already been established.

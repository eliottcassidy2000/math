# Independent audit of the positive insertion-cut theorem

**Verdict: PASS.** Reviewed PROOF.md, check.py and explore.py independently. The author applied the requested terminology clarification: arithmetic circuits on Boolean inputs, rather than standard Boolean circuits with weighted integer gates. No mathematical defect found.

## Proof checks

1. Legal insertions and orphans partition every new Hamilton path. The first contribution charges K(u,v) for each old forward arc cut from bit 1 to bit 0; the second charges K(v,u) with exactly the same bit condition. Neither double-counts paths.
2. K(u,v) is independent of the edge uv because every partition puts u,v in different path components. The evenness argument correctly uses inherited tournament Hamilton-path parity, not a stronger unproved symmetry of K.
3. The positivity construction is valid: take a Hamilton path excluding u,v; before its first vertex not dominating both endpoints insert the free-orientation pair with its last vertex dominating that path vertex. The preceding path vertex, if any, dominates both endpoints. If all path vertices dominate both, append the pair. Splitting its free mutual edge gives one of the two required covers. Together with evenness this proves w≥1.
4. The submodularity defect counts each unordered cross pair exactly once, irrespective of its orientation, and has coefficient 2w. Non-nested sets have a positive defect.
5. Zero cuts are precisely SCC suffix unions. A strong component cannot be split with no outgoing arc; the transitive SCC order then forces a suffix. The count is k+1 including empty and full sets.
6. For a strong tournament all nontrivial cuts are positive. Taking the minimum over all ordered s–t mincuts gives the minimum nontrivial cut. The proof explicitly conditions polynomial cut optimization on already available capacities; endpoint construction remains exponential.
7. Pair-support gate optimality is sound only in the declared additive sum-of-at-most-two-bit-gates model. It is not asserted for shared arithmetic circuits or B2 circuits.
8. Tomography recovers the symmetric capacities and labelled weighted outgoing degrees, not necessarily the tournament orientation. The opposite directed triangle witness is literal and correct.
9. The total capacity is the unordered two-path-cover count, so uniform averaging produces H+p2/2. Transitive capacities have the stated exceptional adjacent value and geometric long-range values.

## Reproduction

The author's optimized checker reproduced certificate.json exactly, including all 33,866 extension states through base order five. The independent `audit.py` imports none of the author's implementation. It constructs K by literal permutations with every split position, ignoring only the split edge; it counts extensions by full permutations. Both ordinary and optimized execution give identical `audit.json` / `audit.optimized.json`:

- 75 base tournaments through order four.
- 1,098 literal one-vertex extensions.
- 16,932 submodularity identities.
- Mutual-edge independence and positive-integral capacities for every tested pair.

These tests supplement, rather than replace, the all-order proofs. Actual tournaments are used throughout, and no unconstrained capacity matrix is treated as realizable. The unrestricted Hamilton-count spectrum and general Boolean circuit lower bounds remain untouched.

## Extension review: exponential labelled response collisions

The subsequently added ordered-triangle family is also correct. For k fixed labelled triangle blocks, independent cyclic orientations produce 2^k distinct labelled bases with identical H, pair capacities and weighted outgoing degrees. Internal capacities are 3^(k−1). Across blocks separated by t intervening blocks, the reverse ordered-cover count is 3^(k−2−t)·3²·12^t=3^k4^t; the forward count is 3^(k−2) for adjacent blocks and zero otherwise. The note correctly distinguishes labelled collisions from nonisomorphic tournaments.

`audit_blocks.py` independently constructs the block tournaments, checks every stated kernel formula and signature for all orientations through k=3, and additionally reconstructs the kernels by literal permutations/splits through k=2. It imports the author's endpoint kernel for the larger check, so this supplementary test is not wholly implementation-independent; the literal small checks are. Its receipt is `audit_blocks.json`.

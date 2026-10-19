# Independent review: rank-staircase two-color reconfiguration

## Verdict

The proof in `UPPER_TAIL_TWO_COLOR_CONNECTIVITY.md` is valid as read on 2026-10-19. No counterexample or missing substantive hypothesis was found. The partial-forest argument in fact proves the corresponding finite-matroid statement. The graphic application proves connectivity of every integer-difference upper-tail completion fiber under simultaneous changes in at most two free colors, with diameter at most 2^(r-1) for r >= 1.

This is an independent proof audit and computational control, not a historical novelty assessment.

## Detailed checks

1. **Induction hypotheses.** After removing A_1, the class A_i has required rank s-i+1, exactly the staircase for its new index i-1 in s-1 classes. No contraction by the current A_1 edge is taken for this recursion, so no rank hypothesis is silently lost.

2. **Contraction and loops.** Let C be the s-2 shared tail elements. It is independent, and the old tail element f and new tail element g are nonloops in M/C, because both tail states are independent. The old A_1 element e is also nonloop, since C+e+f is independent. The inequality

       rank_(M/C)(A_1) >= rank_M(A_1) - rank_M(C) >= 2

   follows immediately from the contraction rank formula. Disjoint color classes ensure the representatives are distinct elements, including when parallel.

3. **Repair dichotomy.** If e and g are independent, make the tail replacement alone. Otherwise they are parallel nonloops. Rank at least two provides h in A_1 outside the closure of g; in particular h is a nonloop and h != e. If h and f are independent, replace e by h and then f by g. If dependent, h and f are parallel nonloops. The old pair {e,f} and new pair {h,g} are independent, and each spans the same two parallel classes. Thus their closures in M/C agree.

4. **Closure lifts to the original matroid.** Equality of the pair closures in M/C implies

       cl_M(C union {e,f}) = cl_M(C union {h,g}).

   This follows from cl_(M/C)(X) = cl_M(C union X) minus C. Consequently the simultaneous replacement is a closure-preserving operation on the full selected s-element set, not merely on the contracted picture. Similarly, if two tail states have equal closure in M, adjoining the unchanged e preserves independence and equality of the full closures.

5. **Length.** Each recursive operation costs at most two operations, and the final A_1 replacement costs at most one. This gives T_s <= 2 T_(s-1)+1 with T_0=0, hence T_s <= 2^s-1. Identical target edges can simply be left unchanged; these bounds do not require counting no-ops.

6. **Graphic lift.** On the r+1 vertices obtained by contracting F, every independent r-1-edge tail state has exactly two components. If a tail single replacement leaves the current color-1 edge bridging the new two components, retain it. Otherwise connected Q_1 has some edge bridging those two new components. Replace the changed tail edge and the color-1 edge simultaneously. This is one at-most-two-color move between full spanning trees. A closure-preserving two-tail-color operation leaves the two components unchanged and needs no color-1 change. The fundamental-cycle version in the proof is also correct. This simpler cut formulation avoids suggesting that the temporary disconnected edge sets are themselves reconfiguration states.

7. **Sharper full-tree bound.** The s=r-1 tail sequence has at most 2^(r-1)-1 operations. Each lifts to exactly one full-tree move at most, followed by at most one target color-1 replacement. Therefore the bound is 2^(r-1), not the larger 2^r-1 general partial-selection bound.

8. **Integer differences.** The full difference-d graph on labels 0,...,m is the union of exactly d nonempty residue chains, for 1 <= d <= m. Adding F cannot increase the number of connected components. Contracting F therefore gives a color-d contact graph on r+1 vertices with at most d components and rank at least r+1-d. Loops and repeated parallel contacts do not affect this rank calculation. The color-1 contact graph is connected.

9. **Existence.** Starting with the highest free color and proceeding downward, when choosing color d there are r-d already-selected independent free edges. Its class rank is at least r+1-d, so that class cannot lie entirely in their closure. Choosing an element outside their closure extends independence. After all r steps this gives a spanning-tree completion. Thus the asserted fibers really are nonempty for every prescribed acyclic upper-tail selection F.

10. **Boundaries and scope.** The stated theorem covers m >= 1 and 1 <= r <= m. For r=1 the fiber is a clique and the bound one holds. If extended to r=0, the fiber is a singleton and its diameter is zero; the formula should not be applied mechanically. For F empty, r=m, the theorem connects all gracefully labeled trees on this fixed interval. Intermediate trees may change abstract shape. This does not establish the Graceful Tree Conjecture, nor one-color global connectivity.

## Independent finite controls

`check_binary_matroids.py` tests finite binary matroids in dimensions s+1 for s=1,...,4, with independently generated disjoint labeled color classes satisfying the staircase. Repeated vectors represent distinct parallel elements, and zero vectors represent loops. For each instance it enumerates every independent rainbow selection. Two vertices are adjacent exactly when they differ in one coordinate, or differ in two coordinates and have equal matroid closure. It computes distances from every vertex, independently of the induction's construction.

Seed: 20261019. Results: 80 instances, 1,882 total independent rainbow states, all connected and all diameter bounds passed. Largest observed diameters by s were 1, 2, 3, and 4. These controls corroborate implementation-level details and do not replace the proof. Full reproducible instance data and results are in `binary_matroid_controls.json`.

## Minor presentation suggestions only

- Insert missing spaces in phrases such as “at most2^s” and “colors2”.
- The cut-based graphic lift in item 6 is shorter than the valid fundamental-cycle explanation.
- If stating the general matroid theorem, say “equal matroid closure”; in the graphic specialization this is exactly equality of connected-component partitions.
- Keep the explicit distinction from fixed-tree-shape graceful labeling and from a one-color theorem.

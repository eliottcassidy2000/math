# Five extra pieces preserve the minimal three-coordinate interface

Session 195, 2026-11-05. **FINITE-EXACT coverage extension with literal physical witnesses.** No eventual-coverage claim is made. The unrestricted spectrum remains **OPEN**.

The complete order-at-most-five moving-root bank in `PROOF.md` has65 distinct matrices and, at arbitrary depth, omits exactly

    {7,21,67,163,173,247,511}

through20,000. Five explicit six-vertex gadgets close the five nonconjectural holes while retaining the same three-coordinate returning state. This is a different parameter choice from the combined stationary-/moving-root four-coordinate construction: the latter needs fewer extra gadgets, whereas this packet keeps only three coordinates.

The gadgets were located using the independently generated Session195 order-six inventory. Completeness of that larger inventory is **not required** here. Each of the five actual graphs and its constants is independently checked by both subset DP and literal permutation/cut enumeration.

## The five additions

For each row, the seven constants are `(g,s,L,E,Q,B,delta)` from Theorem1 of `PROOF.md`. Adjacency rows use outgoing-neighbor bitmasks. The final column below is described in words to keep the construction unambiguous.

- Target67: gadget `([52,17,2,7,12,30], root3)`. Apply it directly to the rooted singleton.
- Target163: gadget `([44,1,2,6,15,30], root3)`. First apply the transitive-edge sink-marked gadget `([0,1],root0)` to the singleton, giving `(3,1,3)`.
- Target173: gadget `([56,1,3,6,14,30], root1)`. First apply the singleton gadget, giving `(1,0,2)`.
- Target247: gadget `([48,1,3,7,14,30], root2)`. First use the same transitive-edge sink-marked gadget, giving `(3,1,3)`.
- Target511: gadget `([42,16,3,6,13,30], root3)`. First use the transitive-three-vertex gadget `([0,1,3],root2)`, giving `(1,0,8)`.

All further applications use the inherited moving-root seam operation; no stationary mark or new operation is introduced. Full eight-coordinate gadget profiles, output rows, roots, and construction words are in `five_gadget_witnesses.json`.

The resulting bank has70 distinct matrices. A complete arbitrary-depth search through1000 retains12,863 states and misses exactly7 and21. This independently checks that the additions do not introduce a mistaken representation of either hostile value.

## Exact bounded conclusion

The original65-matrix class already supplies every odd at most20,000 except its seven listed holes. Every original word is still a word in the enlarged bank, and the five explicit new words fill the other five holes. The exact target-bounded search through1000 excludes7 and21 at arbitrary depth in the enlarged bank. Therefore:

**Finite-exact theorem.** The70-matrix moving-root grammar realizes exactly the odd integers at most20,000 other than7 and21.

This combines an independently reproducible complete base census with five explicit constructions. It does not assert that this bank covers every larger admissible odd integer. The root-content obstruction remains, and the primitive equal-response barrier in `EQUAL_RESPONSE_BARRIER.md` applies even if this finite bank is replaced by all possible moving-root gadgets.

## Reproduction

The [independent audit](../rooted-interfaces-audit/AUDIT.md) confirms all five
physical outputs and independently recomputes the70-matrix closure through1000.

Run:

    python3 verify.py --bound 20000
    python3 verify_extensions.py
    python3 -O verify_extensions.py

The extension verifier counts all65 base representatives and all five new gadgets independently by literal path/cut enumeration, constructs and counts the five complete output graphs, and performs the exact bounded-state closure for the70 matrices through1000. It does not import an expected Hamilton count from the discovery inventory. `extension_receipt.json` records the results.

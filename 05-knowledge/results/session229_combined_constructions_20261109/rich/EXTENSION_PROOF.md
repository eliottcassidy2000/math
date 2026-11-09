# One-vertex extensions of positive unsaturated seeds

**Status: FINITE-EXACT within the explicit universe below. No all-size obstruction.**

The positive seeds are the independently audited B95 and B389 full-interface unit-switch blocks in `PROOF.md`. This follow-on asks whether adding one vertex can preserve their unsaturation while increasing a fixed-interface T fiber to a consecutive triple.

## Exact universe and result

For B95=[64,121,35,37,13,81,28], enumerate all 128 orientations of the new vertex's seven incident edges. For B389=[20,85,72,195,204,159,33,71], enumerate all 256 orientations of its eight incident edges. The new vertex is appended to the outgoing-bitmask list; orientation bit i means new vertex→i. For EVERY resulting tournament enumerate EVERY cut mask (256 for the first seed, 512 for the second). Retain only extensions with gcd(H,Q)=1. Group across all these topologies by the exact interface (H,Q,A,G,S,E); retain the full set of T values for each group.

The results are:

- 384 candidate extensions, of which 240 are unsaturated;
- 97208 distinct retained full interfaces;
- every retained interface has at most two distinct T values;
- maximum consecutive step-two run length is two;
- exactly 40 maximal step-two runs have length two;
- every witness in those 40 runs has orientation 0 or the all-ones mask: only an appended sink or source preserves such a unit switch in this universe.

Thus one-vertex extension does not yield an unsaturated universal triple from either audited seed. More sharply, every surviving unit-switch witness is ordinary transitive padding, already explained by the context theorem. No genuinely non-ordinal extension in this finite universe provides even a full-interface gap-two pair. Nontrivial orientations may still have unequal T values with a larger gap; no stronger absence claim is made.

This gives a principled stopping reason for this extension operation on these two seeds. It does not forbid other eight/nine-vertex cores, multiple-vertex extensions, changing several old arcs, varying the full interface, or other construction grammars. No global coverage claim follows.

## Reproduction

Run:

    python3 -I -B extend_unsaturated.py > EXTEND_UNSATURATED.json
    python3 -I -B -O extend_unsaturated.py > EXTEND_UNSATURATED.optimized.json
    cmp EXTEND_UNSATURATED.json EXTEND_UNSATURATED.optimized.json

The script uses `data` from `search.py` solely as a deterministic subset-path dynamic program. It does not call its random universe generator and reads no discovery logs or prior receipts. All 40 maximal runs and their complete physical witnesses appear in the receipt. Normal and optimized outputs are identical. The earlier standalone main verifier remains unchanged; this follow-on computation is a separate finite census, not a claimed independently audited theorem.

Both extension checkers explicitly load only their sibling `search.py` by absolute file path derived from `__file__`, so isolated Python (`-I`) works from arbitrary working directories. Its `__main__` discovery guard is not executed. `EXTENSION_MANIFEST.sha256` lists the full extension deliverable and dependency hashes.

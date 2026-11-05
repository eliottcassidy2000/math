# Independent review: shared persistent three-port paired fan

**Result: PASS**, with no correction or unresolved blocker. This verifies the common physical three-port exposure and all-order one-generation family, not an iterated ternary grammar.

## Bound source snapshot

Files are from `session197-hybrid-constructions/shared-three-port-fan/`:

- `PROOF.md`: `9baad538608a69f56156ad7fe1076669f08c9670e4f2f9c588d97f4465c8ac72`
- `core.py`: `9c104d778f7c38eee81629a81f9f9e88c2f1197a06c15c6b2d8afd424071d8aa`
- `verify.py`: `8a12bc2da724384328418b1e0f19fb9f42eb86ced952817e588551bd5065f8bb`
- `VALIDATION.json`: `e8bd415c0a5950e37ce2ba474de65b61493342887aa82a200af10b54a76ff234`
- `VALIDATION_OPTIMIZED.json`: `e8bd415c0a5950e37ce2ba474de65b61493342887aa82a200af10b54a76ff234`

The source snapshot is copied under `replay/`. The review did not modify the original files.

## Physical masks and graph identity

The common ordered ports are the actual old vertices `(0,1,4)`. Under the stated left-birth convention, mask 0 creates a universal source. After that birth, mask 1 makes the next newborn lose only to old vertex 0, so it beats the preceding newborn as well as every other old vertex. This gives the displayed A core exactly. Mask 6 selects port positions 1 and 2, hence old vertices 1 and 4, giving the B core exactly.

Drop index 3 removes only the newborn from the four-entry candidate port list. No old port, distinguished root, or core marker is replaced. The new three-port implementation is `verify.py`'s `left_birth`; the legacy two-port convenience functions in the unchanged `core.py` are not used for these moves.

Both canonical-fan relabelings are bijections at every positive m and send the canonical outer root back to actual old vertex 5. Checking each kind of edge gives the exact fan graph: old core incidences, the two new core vertices or single new core vertex, the ordered chain, its edges to the outer root, and the unique reversed outer-root/core-marker edge. In particular the second plus2 birth beats the first one, and both beat every chain vertex and old root.

## All-order formula and fixed-root profile

In S=X joined forward to a nonempty transitive chain, H(S)=H(X) and no spanning path ends at old core marker a. Every ordered cover with its second slot starting at a restricts to one such cover of X; each chain vertex can independently be put into either slot, in increasing order. This proves the exact factor 2^m, including an originally empty first slot.

Appending b with sole outneighbor a then gives H=h+2^m R_a(X). The ending-at-b count is h because every S path ends in the nonempty chain. Deleting b also leaves h paths. A cover whose second path starts at b either uses b alone, contributing h, or continues b to a and then corresponds exactly to a cover counted by R_a(S). Thus R_b=h+R_a(S)=H. This proves the old-vertex profile `(H,h,H,h)` without changing the root.

The positive-chain condition is essential for the zero ending-at-a term and is explicitly present throughout. Strongness follows from the strong finite core, core-to-chain-to-b paths, and b returning to a. The resulting counts are exactly 3H+2 and 3H+4 for every member of the parent family.

## Numerical family does not return

The retained ports are identical after either word, but the core changes from C to A or B. The source correctly does not reapply the family identity to those new cores. The literal sequences 49,149,515 and 49,151,409 disprove naive repetition of the two rules. No target inversion, closed state class, bounded span, or universal ternary induction is asserted.

## Replays and independent checks

- Source normal and optimized runs match both frozen receipts byte-for-byte.
- `independent.py` constructs every fan and birth directly from the edge definitions and uses a separate forward endpoint DP, rather than the source's predecessor recurrence.
- All finite core counts and root responses agree.
- All-edge relabelings, unchanged actual root 5, strongness, and all three old-root profiles pass for chain lengths 1 through 8, beyond the source's range through 6.
- The general fan formula is also checked on every rooted labelled core through order four, at chain lengths 1,2,3: **855** additional physical cases, including nonstrong cores and singleton boundary cases.
- Both nonreturning repetition counterexamples reproduce exactly.

The arbitrary-core formula-data digest is `008d9de6544026474845db9bae98aa446cb54bbea42075b2c79e8c0c068f43ab`. The all-order conclusion rests on the explicit edge identities and cover bijection, not these finite controls.

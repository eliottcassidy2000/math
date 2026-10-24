# Stable tournament cells, signed holonomy and refinement-error transfer

This packet proves exact automorphism recovery for stable cells of sizes1/3,
structural refinement exactness for strong order8, and the additive SCC law
E=A²E_s. A complete regular-seven census supplies E_s,7=14, implying E8=28.
It supplies a strong rigid9 counterexample, an all-r matching-triangle family
count, and rigorously scoped lower bounds E9>=114 and E10>=388.

- `THEOREM.md`: complete proofs, inheritance, scopes and remaining boundary.
- `check.py`: standalone finite controls, independent of the growth engine.
- `receipt.json`: finite universe counts, witnesses and control results.
- `VALIDATION.md`: precise checked universes and reproduction instructions.
- `pinned.txt` and `archive/`: byte-pinned inherited sources.
- `manifest.json`: hashes for frozen payload files.
- `replay.py`: verifies hashes, copies only the manifest-listed payload to an
  isolated temporary directory, reruns both normal and optimized Python,
  and compares the exact receipt.

The all-r count is family-relative. Stable refinement is not a complete
isomorphism invariant. The order-eight correction is structural, but uses
the finite regular-seven input. No bare-line-graph closure or all-order
A000568 recurrence is asserted.

# Session51: implication intervals and neutral packet cycles

**ABSTRACT, PROVED all-size; FINITE-EXACT controls.** No native electrical-current realization is claimed.

Start with [THEOREM.md](THEOREM.md). It specializes Session49's arbitrary endpoint-interval routing theorem to Boolean unary/binary relations:

1. Switch-gauged endpoint intervals are exactly ideals of the SCC condensation poset; complete additive guards select the zero-defect ideals.
2. Indivisible accepted events are relative zero-defect ideal atoms, equivalently covers of the accepted ideal poset, with full physical slot and clock costs.
3. Strict routes are exactly neutral packet partitions whose dependency contraction is acyclic. Merging SCCs of any neutral packet dependency graph gives a constructive accepted route, charging the whole merged packets.
4. A path of binary implications, with SCC size one, treewidth one, fixed C2² charges and every physical charge increment nonzero, can force an arbitrarily large event. Even lengths preserve one exact fixed population. A paired lift has margin-neutral SCCs of size two and the same obstruction.

This does not refute Session49's correctly scoped additive-only Davenport theorem or its all-neutral-block antimatroid theorem. It identifies why combining bounded synchronization blocks and finite-group extraction needs a dependency-cycle control.

Reproduce from any working directory, standard-library Python only:

    python3 /path/to/implication-intervals/verify_manifest.py
    python3 /path/to/implication-intervals/check.py
    python3 -O /path/to/implication-intervals/check.py
    python3 /path/to/implication-intervals/portable_check.py

`receipt.json`, `normal-check.txt`, and `optimized-check.txt` are identical exact output. `PORTABILITY.json` records a real temporary-directory copy, normal/-O agreement, and five fail-closed manifest mutation controls. The manifest excludes only itself and runtime bytecode caches. `SOURCEPINS.json` authenticates the inherited text snapshots against the specified repository commit and Git blobs.

See [VALIDATION.md](VALIDATION.md) for universes, boundaries, and independent-review status, and [RESEARCH-NOTES.md](RESEARCH-NOTES.md) for the scoped inheritance/connection ledger.

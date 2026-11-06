# Scalar unit-start witness packaging review

PASS for the inspected release candidate in `session205-scalar-release`.

The production Python enumeration now uses an exhaustion loop. The C++ witness generator uses the same proven scalar pruning, retains a recipe with every four-coordinate state, and selects a deterministic minimum base-four recipe for every supported count. Its operation codes 1,2,3 decode to O,Z,J without ambiguity; the empty recipe is represented by a dash. The explicit encoding guard stops with nonzero status before a too-long recipe could overflow. All actual recipes are well below the guard.

I compiled and reran the witness generator in a fresh temporary directory. Both `WITNESSES.tsv` and `WITNESSES.out` match the release bytes exactly: 24,971 witnesses, 4,082,826 level-states, frontier exhausted after order21. I independently reran `check_witnesses.py`, reproducing `WITNESS_CHECK.json` byte for byte. The checker reconstructs every literal tournament and validates its orientations; counts use the inherited exact physical transfer, correctly distinguished in the result note from literal subset-DP recounting. The previously audited 364 small literal recipes remain the independent direct-count controls.

The result retains the finite scope, singleton inclusion via scalar1, unit-start versus unrestricted-spectrum distinction, and missing all-order occupancy lemma. No additional mathematical claim is inferred from recipe packaging. The frozen source grammar is included for the inherited completeness claim.

For final release, freeze a SHA256 manifest after any status-only edits and verify it before publication. A status-only promotion from 'recipe packaging is under final review' to completed independent packaging review is approved. No mathematical changes are required.

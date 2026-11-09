# Reproduce Session232

Python 3.10+ and its standard library suffice. No network, repository imports,
external packages or hidden databases are needed for the session's controls.
The all-order conclusions also depend on the explicitly named inherited
THM-4812 theorem and its published finite certificates, which this replay does
not recompute.

From this directory:

    python3 -B replay.py > REPLAY.new.json
    cmp REPLAY.new.json REPLAY.json
    sha256sum -c SHA256SUMS

The replay copies the packet into isolated temporary directories and executes
eight producer/independent programs normally and under -O. It requires zero
exit, no stderr, exact equality between modes and byte equality with frozen
receipts. Every checker is AST-checked for removable assert statements.
There are16 runs, including full subgroup controls through order7, all four
producer lanes and all three independently written audit algorithms.

Individual commands and exact universes are recorded in each proof or audit.
The separate subgroup audit also contains eight deliberately false failure
sentinels; each must raise RuntimeError in its normal and optimized mode.

## File identity and publication adaptations

Producer files are copied verbatim except the packing proof's single status
line, promoted after PASS. Its before/after hashes and exact change are in
packing_mass/PROMOTION.json. The packing auditor's recorded proof hash refers
to the audited pre-promotion candidate; mathematical text is unchanged.

The packing and graph audit harnesses now prefer the published sibling folder
names and retain their original scratch-name fallbacks. The graph independent
algorithm is otherwise unchanged. The subgroup audit's shell command typo
`python=3` was corrected to `python3`, and its manifest records the original
report hash. These are publication/documentation changes, not mathematical
repairs. The root SHA256SUMS is the final authoritative packet inventory.

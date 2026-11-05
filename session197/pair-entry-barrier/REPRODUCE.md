# Reproduce the reviewed pair-entry response barrier

This publication preserves the six frozen source files, their original hash manifest, and the independent audit without changing their bytes. The result is an all-order marked-response obstruction; it does not exclude any scalar count attained by the displayed graphs.

From this directory, with standard-library Python 3:

```sh
sha256sum -c RELEASE_FILES.sha256
python verify.py > /tmp/pair-barrier-normal.json
python -O verify.py > /tmp/pair-barrier-optimized.json
cmp /tmp/pair-barrier-normal.json VALIDATION.json
cmp /tmp/pair-barrier-optimized.json VALIDATION_OPTIMIZED.json
python audit/independent.py > /tmp/pair-barrier-independent.json
cmp /tmp/pair-barrier-independent.json audit/independent.json
```

The main checker needs only the included core.py and control_graphs.json. The independent checker is standalone and does not import that core. There are no omitted runtime dependencies or external packages.

The frozen proofs retain their original research paths for provenance. Exact copies of the two cited mathematical proof snapshots are included as dependencies/PAIR_ENTRY_PROOF.md and dependencies/ARITHMETIC_INVERSE.md, with hashes matching the theorem's dependency list. Their original published packets are session197/pair-entry-bridge/ and session197/lowcycle/arithmetic-inverse/.

PUBLICATION.sha256 covers every delivered file except itself, including the independent audit and this reproduction note.

# Publication and reproducibility checks

The first checkpoint requires seven independent normal/optimized frozen
replays, exact output agreement, unchanged input bytes, native checksum and
packet file-set verification, docs/budget/banner checks, explicit owned-path
staging, and full remote tree/blob/mode/ancestry verification.

Local replay and documentation receipts are sealed into the packet before
publication. Remote Git receipts and exact-commit CI observations are produced
after publication in the external session66-publication directory. A GitHub
status with zero checks is absence of checks, not a CI pass.

## First verified local gate

`FIRST_REPLAY_RECEIPT.json`: 7 scripts, 14 successful normal/optimized runs, exact saved outputs and unchanged frozen input bytes. `FIRST_DOCS_CHECK.txt` passes documentation guards. `FIRST_BUNDLE_CHECK.json` verifies source pins, UTF8, links and optimization-safe checks. `FIRST_FILESET_NEGATIVE_CONTROLS.json` rejects changed/extra/symlink mutants in both modes. Remote verification follows publication and is not presumed here.

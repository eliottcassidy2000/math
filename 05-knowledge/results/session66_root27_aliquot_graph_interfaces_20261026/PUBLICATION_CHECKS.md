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

## Second verified local gate

`SECOND_REPLAY_RECEIPT.json`: 2 scripts, 4 successful normal/optimized runs, exact saved outputs and unchanged frozen input bytes. `SECOND_DOCS_CHECK.txt` passes documentation guards. `SECOND_BUNDLE_CHECK.json` verifies source pins, UTF8, links and optimization-safe checks. `SECOND_FILESET_NEGATIVE_CONTROLS.json` rejects changed/extra/symlink mutants in both modes. Remote verification follows publication and is not presumed here.

## Historical original-third local gate (unpublished draft)

`THIRD_REPLAY_RECEIPT.json`: 18 scripts, 36 successful normal/optimized runs, exact saved outputs and unchanged frozen input bytes. `THIRD_DOCS_CHECK.txt` passes documentation guards. `THIRD_BUNDLE_CHECK.json` verifies source pins, UTF8, links and optimization-safe checks. `THIRD_FILESET_NEGATIVE_CONTROLS.json` rejects changed/extra/symlink mutants in both modes. Remote verification follows publication and is not presumed here.

## Combined final verified local gate

`FINAL_REPLAY_RECEIPT.json`: 20 scripts, 40 successful normal/optimized runs
with PATH=/nonexistent, exact saved stdout and unchanged frozen/source bytes.
`FINAL_BUNDLE_CHECK.json` checks all234 approved source files, native checksum
seals, exact historical-script exemptions and pinned-context links.
`FINAL_DOCS_CHECK.txt` is a fresh documentation-budget/banner/status/router pass.
`FINAL_FILESET_NEGATIVE_CONTROLS.json` tests modified, extra and symlink files
in both Python modes. The old third receipts above remain historical, unchanged.
`CURRENT_CORRECTION_VERIFICATION.json` records the original reviewed draft;
`FINAL_INCOMING_PRESERVATION.json` updates its MISTAKES union against the final
incoming base, without deleting any concurrent entry. Remote verification is
required after publication; these local receipts do not claim a remote push.

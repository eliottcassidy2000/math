# Publication and reproducibility checks

## First checkpoint

`FIRST_REPLAY_RECEIPT.json` records 5 mathematical/provenance scripts and 10 successful normal/optimized runs. Exact saved outputs match; original input bytes stay unchanged. Replayed source and approved-copy hashes were rechecked. `FIRST_DOCS_CHECK.txt` passes bounded documentation checks. `FIRST_BUNDLE_CHECK.json` verifies optimization-safe scripts, checksum pins, UTF-8 and original-context links. `FIRST_FILESET_NEGATIVE_CONTROLS.json` rejects changed-byte, extra-file and broken-symlink mutants in both Python modes. Exact Git parent/tree/path/blob/mode checks and CI visibility are performed after the remote publication, not assumed by this prepublication record.

## Second checkpoint

`SECOND_REPLAY_RECEIPT.json` records 14 mathematical/provenance scripts and 28 successful normal/optimized runs. Exact saved outputs match; original input bytes stay unchanged. Replayed source and approved-copy hashes were rechecked. `SECOND_DOCS_CHECK.txt` passes bounded documentation checks. `SECOND_BUNDLE_CHECK.json` verifies optimization-safe scripts, checksum pins, UTF-8 and original-context links. `SECOND_FILESET_NEGATIVE_CONTROLS.json` rejects changed-byte, extra-file and broken-symlink mutants in both Python modes. Exact Git parent/tree/path/blob/mode checks and CI visibility are performed after the remote publication, not assumed by this prepublication record.

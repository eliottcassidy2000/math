# Session190 lossless receipts

Every checkpoint manifest records the exact original byte count, SHA-256 digest, and Git blob SHA-1 for every retained file. Files at most 48 KiB are stored directly. Larger files are deterministically gzip-compressed and split into chunks at most 32 KiB; no original data are discarded. Each chunk and the complete compressed stream have their own verification fields.

From the repository root, run:

```sh
python 05-knowledge/results/session190_constructive_links_20261104/RESTORE_ARCHIVES.py \
  05-knowledge/results/session190_constructive_links_20261104/ARCHIVES_CP01.json
```

Pass additional manifest paths for later checkpoints, or `--verify-only` to verify retained direct files and the compressed reconstruction without writing restored files. Existing identical files are checked; an existing file with different bytes causes an error and is never overwritten. Only explicitly listed paths are read or written.

Restoration validates the exact original bytes and both hashes. Individual packet notes describe reproduction commands and which fields, if any, depend on the local checkout path. The manifests preserve original historical receipts even when their informational local-path fields change on replay.

Publication uses non-force expected-head updates that preserve incoming main commits. Remote byte, hash, ancestry and CI observations are recorded separately after each update. No reported CI runs or contexts is not a CI pass.

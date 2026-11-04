# Session190 lossless receipts

Every checkpoint manifest records the exact original byte count, SHA-256 digest, and Git blob SHA-1 for every retained file. Files at most 48 KiB are stored directly. Larger files are deterministically gzip-compressed and split into chunks at most 32 KiB; no original data are discarded. Each chunk and the complete compressed stream have their own verification fields.

For the current checkout, supply **all available checkpoint manifests in ascending checkpoint order**. From the repository root, the zero-padded filenames make this glob order suitable:

```sh
python 05-knowledge/results/session190_constructive_links_20261104/RESTORE_ARCHIVES.py \
  05-knowledge/results/session190_constructive_links_20261104/ARCHIVES_CP*.json
```

When a path appears in multiple manifests, the final supplied record is its current expected version; earlier records are explicitly reported as superseded. This allows navigation and publication tooling to evolve while every mathematical packet remains recoverable in one checkout. To reproduce historical navigation/tool bytes instead, use the immutable earlier checkpoint and its manifest. Never supply manifests in reverse order.

Use `--verify-only` to verify retained direct files and the compressed reconstruction without writing restored files. Existing identical files are checked; an existing file with different bytes causes an error and is never overwritten. Only explicitly listed paths are read or written.

Restoration validates the exact original bytes and both hashes. Individual packet notes describe reproduction commands and which fields, if any, depend on the local checkout path. The manifests preserve original historical receipts even when their informational local-path fields change on replay.

Publication uses non-force expected-head updates that preserve incoming main commits. Remote byte, hash, ancestry and CI observations are recorded separately after each update. No reported CI runs or contexts is not a CI pass.

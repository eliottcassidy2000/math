# Checkpoint 3 curation — lossless logical transport

Root approved this curation on 2026-10-27. No remote publication was performed
by the offline preparer.

Every original native-manifest entry and every original manifest is retained.
The 9 lanes contain 600 logical files (109,323,086 bytes). Seventy data files
larger than 200,000 bytes are transported as deterministic gzip (level 9,
mtime 0); only `.json`, `.tsv`, `.txt`, or `.out` data qualify. No script,
proof, report, or original manifest is compressed. Total native transport is
25,113,163 bytes. The largest original ledger is 30,053,175 bytes and its gzip
transport is 4,343,660 bytes. There is no generation-only omission.

Each logical path has its original SHA256, byte count and file mode, plus
stored path/SHA256/size/mode and explicit identity/gzip encoding in
[TRANSPORT.checkpoint3.json](TRANSPORT.checkpoint3.json). The verifier checks
stored bytes first, limits expansion to the declared size, validates expanded
bytes, checks exact original-manifest membership and every native hash, and
refuses unsafe paths, duplicate names, symlinks, unexpected files and mode
changes. Materialization requires a new destination outside the packet.

| Native lane | Files including manifest | Frozen manifest SHA256 |
|---|---:|---|
| `session75-black-boundary` | 44 | `5bc89f4d9673325311088321d984811a2ddbf123be6739e81ff75b314f5d28a3` |
| `session75-black-boundary-audit` | 54 | `0fbfca814bf71da72b9890762007eab5f4193be96a5ecedda1a3c0c38589077b` |
| `session75-black-twins` | 9 | `dcb9b8d6d65122e7e47dba0f9a3984c145bdcf67fbf0dd059db4892b54796d2b` |
| `session75-joint-cutoffs` | 130 | `09f5e8372a8133109cec03447284c49881440b36521eb22d2ef2ea7c67cf0ef7` |
| `session75-joint-cutoff-audit` | 13 | `395a437814119fd1402230b76a6eea2b6c80b1c322dad6105c883bd4ec392802` |
| `session75-tc-audit` | 39 | `2afbb89858a121c063ce35edece62bf5c3bf9c6108d56c9fdc370e32dd3edd54` |
| `session75-tc-window-audit` | 286 | `2f6a8ada905f42f4bb0dfb63e0361c6b86955ed9bf9907dd2a1becb6e7eb77c2` |
| `session75-black-twins-audit` | 21 | `f80e3380c7e8622f06113a9a838f451f8be80c1ed761917a0195c30de551be01` |
| `session75-root-twins` | 4 | `014ea5d7b65cc7fddadd3c952fc9d3d665ea3ccd9a15aa1231bc1534c513027a` |

The root-twins manifest was newly assembled from exactly the root's unchanged
`check.py`, `normal.json`, and `optimized.json`; all other listed manifests are
the original frozen producer/auditor manifests. Compiled binaries and transient
logs outside those original manifests were never native packet entries and
are not silently promoted to publication evidence. Replays compile afresh in
external work directories.

The producer/auditor's six references to three raw pre-compaction witness
streams are preserved historical hashes, not logical files in their manifests.
Complete compressed witness records and source-universe authentication are
retained and replayed. No entry of an original manifest is omitted for this
reason.

[PROVENANCE.checkpoint3.json](PROVENANCE.checkpoint3.json) checks 211 explicit
ledger-field hash bindings and 31 byte-identical producer/auditor copy bindings;
[SOURCE_ORIGINS.checkpoint3.json](SOURCE_ORIGINS.checkpoint3.json) retains 44
exact targeted Git-origin comparisons. Hash equality is not described as a
re-audit of every inherited mathematical claim. Frozen “audit pending” or
“OPEN” statements remain as time-local evidence; current scoped acceptance
is in [README.checkpoint3.md](README.checkpoint3.md).

The 206 pre-existing C1/C2 packet paths have identical Git blob IDs and modes
to local C2 `edebe48cd565`; the matching tree-list SHA256 is
`6c8ad910d3fe495de40eda8508bff0ad74945984b653dba2578ed5febe6e64b8`.
Neither prior manifest nor prior helper was rewritten to accommodate C3.
The new whole-packet manifest inventories all old files plus this addition.

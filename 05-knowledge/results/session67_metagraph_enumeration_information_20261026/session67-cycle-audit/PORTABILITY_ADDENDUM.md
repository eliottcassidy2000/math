# Portable offline replay addendum

**Status: PASS; additive derivatives only.** All three original verifier files,
all original reports and receipts, and all three original SHA-256 manifests
remain byte-for-byte unchanged. Their hash checks were rerun successfully.
Use the portable derivatives below for a published packet; the native scripts
remain historical evidence of the originally executed audits.

## What offline provenance does—and does not—verify

`audit_core_triangle_portable.py` accepts `--metagraph-dir PATH`. Its default
mode needs only the frozen peer receipt, manifest and source snapshots:

1. Recalculate each snapshot's SHA-256 and compare it with the supplied manifest.
2. Recalculate its Git blob identifier locally as
   `SHA1(b"blob "+ascii_decimal_length+b"\0"+file_bytes)`.
3. Perform every original numerical/asymptotic-bound check unchanged.

This validates **snapshot consistency with the supplied manifest**. It does
**not** establish that those files belong to the manifest's named commit, nor
independently authenticate that manifest. The default receipt explicitly says
`provenance_mode: "offline-snapshot"` and `named_commit_verified: false`.
It omits exactly the seven named-commit comparisons, leaving10,657 checks.
No Git executable, local repository, network or installation is required.

The original audit really did compare all seven files with commit
`b4234f728c0516b5fc521c001d69e26f3d2508d5`; this historical verification is
recorded in `CORE_TRIANGLE_AUDIT.md` and the unchanged original receipts.
Historical verification is not relabelled as a current offline check.

To perform that stronger check now, explicitly supply `--repo PATH` pointing
to an existing checkout that contains the named commit. This mode executes
`git show <declared-commit>:<repository-path>` for each snapshot, and compares
its bytes. It performs no fetch. Missing commits or mismatched bytes fail the
run rather than silently falling back to offline mode. Successful receipts say
`provenance_mode: "named-commit"`, `named_commit_verified: true`, and contain
all10,664 checks. This mode was also actually rerun normally and under `-O`.

`audit_recurrence_cross_portable.py` accepts `--recursion-dir PATH`. It changes
only the peer-lane path selection and adds a derivative-source hash to the
receipt. All3,195 mathematical checks and outputs are unchanged. It makes no
named-commit verification claim.

## Three portable replay commands

Run from the published audit lane. The defaults expect sibling lanes named
`metagraph` and `recursion`; the explicit options support any other layout.
The temporary output directory avoids overwriting frozen receipts:

```sh
OUT=$(mktemp -d)

# 1. Original self-contained cycle-index and literal-extension audit.
python3 verify.py > "$OUT/cycle.normal.json"
python3 -O verify.py > "$OUT/cycle.optimized.json"
cmp "$OUT/cycle.normal.json" "$OUT/cycle.optimized.json"
cmp "$OUT/cycle.normal.json" receipt.json

# 2. Portable core-triangle audit, default offline provenance scope.
python3 audit_core_triangle_portable.py --metagraph-dir ../metagraph \
  > "$OUT/triangle.normal.json"
python3 -O audit_core_triangle_portable.py --metagraph-dir ../metagraph \
  > "$OUT/triangle.optimized.json"
cmp "$OUT/triangle.normal.json" "$OUT/triangle.optimized.json"
cmp "$OUT/triangle.normal.json" CORE_TRIANGLE_PORTABLE_RECEIPT.json

# 3. Portable recurrence/canonical-root cross-audit.
python3 audit_recurrence_cross_portable.py --recursion-dir ../recursion \
  > "$OUT/recurrence.normal.json"
python3 -O audit_recurrence_cross_portable.py --recursion-dir ../recursion \
  > "$OUT/recurrence.optimized.json"
cmp "$OUT/recurrence.normal.json" "$OUT/recurrence.optimized.json"
cmp "$OUT/recurrence.normal.json" RECURRENCE_CROSS_PORTABLE_RECEIPT.json
```

Optional direct named-commit replay, replacing `/path/to/math` with an existing
checkout containing the declared base commit:

```sh
python3 audit_core_triangle_portable.py --metagraph-dir ../metagraph \
  --repo /path/to/math > "$OUT/triangle.commit.normal.json"
python3 -O audit_core_triangle_portable.py --metagraph-dir ../metagraph \
  --repo /path/to/math > "$OUT/triangle.commit.optimized.json"
cmp "$OUT/triangle.commit.normal.json" "$OUT/triangle.commit.optimized.json"
cmp "$OUT/triangle.commit.normal.json" CORE_TRIANGLE_PORTABLE_COMMIT_RECEIPT.json
```

All stored portable receipts are deterministic and contain no absolute path.
Their embedded peer-file hashes intentionally change if a supposedly frozen
peer report, script or receipt changes; copying all lanes unchanged preserves
byte-identical outputs.

## Isolated portability test and unchanged mathematical output

Both derivatives were copied with their required inputs into a fresh temporary
packet containing only `cycle-audit/`, `metagraph/` and `recursion/`. There was
no sibling repository. `PATH` was set to `/nonexistent`, with Python invoked
by its absolute executable path. Both default-mode commands succeeded and
returned bytes identical to the stored portable receipts. See
`PORTABLE_ISOLATED_REPLAY.json`.

Every original mathematical output was compared against its derivative:

- Triangle offline: all original fields agree except check totals/provenance;
  precisely seven named-commit checks are absent.
- Triangle optional-commit: all original check counts and mathematical fields
  agree; only explicit derivative/provenance metadata is added.
- Recurrence: every original field agrees; only derivative metadata is added.

Normal/`-O` output equality passed in all three derivative modes. These are
actual executions, not inferred portability claims.

## Publication inventory

`PUBLICATION_MANIFEST.json` lists the curated audit files and SHA-256 values,
separates frozen originals from additive portable files, and records cross-lane
dependencies. `PUBLICATION_SHA256SUMS` additionally hashes that manifest.
`PORTABLE_SHA256SUMS` covers only the new portable artifacts. None replaces or
rewrites an older hash manifest.

Exclude scratch `input/`, `__pycache__/`, and local temporary replay files from
the published audit lane. The actual seven source snapshots live in the peer
metagraph lane under its unchanged `source_manifest.json`.

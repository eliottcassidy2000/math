# Portability amendment — explicit OFFLINE pinned-snapshot replay

2026-10-26. Packaging correction only; no new mathematics.

## Why a separate entry point exists

Historical `check.py` validates eleven source files with a live `git show`
against `../session65-repo`. That was appropriate during the mathematical
audit, but a published evidence packet intentionally has no Git checkout.
Its source gate must not be silently skipped or supplied with fabricated
Git output.

The new `portable_check.py` is the **OFFLINE_PINNED_SNAPSHOTS** entry point.
It has the same nine mathematical/helper function bodies as historical
`check.py`, checked byte-for-byte at runtime. It replaces the live-Git gate
with full verification of the saved, pinned evidence closure. Existing
`check.py`, `REPORT.md`, both historical mathematical receipts, both producer
replay receipts, and historical `MANIFEST.json` are unchanged.
`PORTABILITY_DIFF.patch` records the complete difference between the historical
checker and this new file; it is an explanatory diff, not an applied edit.

## The replacement source gate

Before executing the producer implementation, the portable entry point:

1. verifies fixed SHA-256 hashes of the producer's report, checker, both
   original outputs, source-pin file, and manifest;
2. verifies all 17 producer manifest entries, including every saved input;
3. verifies the sealed historical audit manifest and its six entries;
4. verifies the whole `OFFLINE_SOURCE_PINS.json` against its fixed hash;
5. requires exactly eleven unique source snapshots and exact linkage to the
   original SHA-256 pins and source commit;
6. verifies every snapshot's byte count, SHA-256, and saved Git blob object ID;
7. verifies unchanged mathematical function bodies and active exception gates;
8. computes the mathematical consequences and compares every consequence
   object with both historical normal/optimized receipts.

The offline source-pin manifest SHA-256 is
`7c3d4df7cba05caccd2d99b428e9e9654f746a3feae95110127ebd155b55dd39`.
Its eleven Git blob SHA-1 IDs were acquired by read-only `git rev-parse`
against the already-audited source commit
`6a8ad5d3bf01a232da4de41bcbff68390a9337ee`, then independently matched with
`SHA1("blob " + byte_count + NUL + bytes)`. SHA-256 remains the content seal.
This acquisition was performed during packaging; the portable runtime
contains no Git command or subprocess import.

**Provenance boundary:** offline replay proves equality with the saved
historically audited source bytes and their recorded object IDs. It does not
perform a fresh Git readback, re-establish commit reachability, contact a
remote, or verify current branch state. Its JSON explicitly says
`fresh_git_readback_performed: false`; no historical Git result is presented
as a new Git observation.

## Exact invocation in the published packet

Keep these two directories as siblings:

```text
<packet>/session65-reductions/
<packet>/session65-audit-reductions/
```

From `<packet>`:

```bash
python3 -I -B session65-audit-reductions/portable_check.py > portable-normal.json
python3 -I -B -O session65-audit-reductions/portable_check.py > portable-optimized.json
cmp portable-normal.json portable-optimized.json
```

The entry point resolves paths relative to its own file, not the working
directory. Thus an absolute path to `portable_check.py` is also valid. It
requires only standard-library Python and the two saved directories; no Git
checkout, network, ambient `PYTHONPATH`, or installed third-party package.
`-I` deliberately ignores ambient Python path configuration.

The reproducible packaging/hostile test is:

```bash
python3 -I -B session65-audit-reductions/portable_hostiles.py
```

It creates an isolated temporary packet with no `.git` or `session65-repo`,
then invokes the runtime with Python isolated mode and an environment whose
only variable is `PATH=/nonexistent`.

## Verified outcomes

* Portable normal and optimized JSON are byte-identical.
* Their `codec`, `continuation_hostile`, and `lrc_threshold_controls` objects
  exactly equal both historical receipts: 6,586 full-label tables, 45 canonical
  arc corruptions, the 123/135 continuation hostile, and 774 threshold cases
  remain unchanged. Only explicitly labelled source/provenance metadata differs.
* All eleven saved source snapshots pass the complete offline gate.
* For **each** snapshot, a wrong-byte mutation and a missing-file mutation
  are rejected in **both** normal and optimized runtimes: 22 cases, 44
  negative executions. Each failure identifies the affected snapshot and
  emits no success JSON.
* An altered recorded Git blob ID in the offline provenance manifest is
  rejected in both modes. Restoring all files restores the passing output.
* The isolated no-Git/no-ambient-PYTHONPATH positive test passes in both modes.
  The hostile-test driver also produces byte-identical normal/optimized
  receipts.

New evidence files are sealed by `PORTABILITY_MANIFEST.json`. The historical
`MANIFEST.json` remains a valid seal for the original audit; it intentionally
does not pretend the later amendment existed during that audit.

No repository files were edited, staged, committed, or pushed for this
amendment. No external source or downloaded code was introduced.

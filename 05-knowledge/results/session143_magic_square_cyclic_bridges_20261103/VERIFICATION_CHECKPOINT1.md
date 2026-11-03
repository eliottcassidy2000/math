# Verification and publication scope

## Checkpoint 1

The exact universe is all 1024 labeled five-vertex tournaments and all 120
vertex permutations. Three characteristic-polynomial paths agree, and the
primary McKay catalogues are canonically matched record by record. A separate
implementation performs 147 explicit checks without importing producer code.
Its mathematical review is [here](independent_checkpoint1/REVIEW.md).

The C5 operator checker executes 45 portable checks. Six more source-assisted
checks align the cited degree4371 character with all 184 ATLAS class labels;
they were run against hashed primary snapshots. Those third-party snapshots
are not distributed in this packet. Their frozen receipt is evidence of that
transcription pass, not a claim that the portable runner downloads them or
rebuilds the character table. `sporadic/test_output.json` records the ten
source-assisted positive/negative cases; `test_output_portable.json` records
only the six source-free cases that the portable replay can reproduce.

`replay_checkpoint1.py` copies the frozen source selection to a new temporary
location and runs from an unrelated working directory, under normal Python
and `-O`. It compares exact outputs and verifies selected historical snapshots
by both SHA256 and Git blob identity. Deliberately wrong census, inverse,
eigenvalue, character, and independent-spectrum checks reject nonzero.
A same-size proof-file mutation is rejected by the content manifest.
See [the replay receipt](REPLAY_CHECKPOINT1.json) for exact execution counts.
Checks use explicit exceptions, not optimization-sensitive assertions.

The producer and independent manifests are frozen inventories. The aggregate
source manifest includes all proof/checker files needed by this replay, but
not its own hash or the report produced by running it. This avoids circular
self-attestation. Curated navigation and release prose are additionally hashed
in the external publication manifest and verified against the publication tree.

## Deliberate limits

- Mathematical status: **PROVED scoped / FINITE-EXACT / CITED**, as labeled.
- GAP/CTblLib script: provided, **NOT EXECUTED** in this environment.
- Baby Monster generators, representation matrices and an intertwiner:
  **NOT CONSTRUCTED**; no tournament or24-dimensional linear realization.
- External character table: cited source and transcription checked, not rebuilt.
- Historical source snapshots: immutable provenance, not blanket endorsement.
- `agents/check_docs.py`: **NOT RUN repository-wide**, because its central
  mistakes-ledger access is outside this task's permitted scope. The explicit
  scoped gate checks startup budgets, unchanged startup pages, owned-file
  conflict markers, and curated publication links only.
- Central `01-canon/MISTAKES.md` and prior Session113 approval: **not read or
  edited by this release**. Existing canon and unrelated task files remain
  unchanged; only the named source snapshots were read for inheritance.
- Main square-of-squares existence problem and unrelated conjectures: **OPEN**.

Checkpoint preparation does not itself establish remote publication. Remote
parent lease, returned tree/commit identities and final upstream status must
be independently checked by the publishing coordinator.

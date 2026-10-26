# Publication and reproducibility gates

The first checkpoint (`fc3768ac`) contains only the independently accepted local
enumeration lane. The final checkpoint adds the separately audited remaining lanes.

The root enumerator was rebuilt with `-O0` and `-O3 -DNDEBUG`, C++17,
`-Wall -Wextra -pedantic`. Both independently regenerated all 17 census/table
files through order 8; every byte matched (order-8 line bytes after decompression).
Explicit exceptions, not removable assertions, enforce its invariants.

The independent auditor checks all production LUT destinations and legality
against complete arc relations through order 8; an independent all-permutation
class oracle through 6; complete node and coloured-line tables through 7; and
all order-8 line/node mass and degree accounting. It separately checks the
connected-plus-closed criterion on all 4,140 partitions of the order-4 cube,
and the mod-two Laplacian hostile. The audit report states exact scope.

Root exported-ledger verification passes in normal and optimized Python modes:
444,711 line rows and 7,412 ordinary nodes over orders 1–8. These totals count
stored rows/nodes across orders, not ordinary classes of one order.

`FIRST_REPLAY_RECEIPT.json` records the local replays. `FIRST_DOCS_CHECK.txt`
records repository documentation guards. `FILE_MANIFEST.json` seals the exact
curated file set, excluding only itself. Temporary binaries, caches and the
exploratory path-move probe are excluded. `verify_packet.py` rejects changed,
extra, missing and symbolic-link files and checks frozen inherited-source pins.

Remote verification is performed after publication: exact target commit,
ancestry, owned-path Git blob hashes/modes and checks for that exact SHA.
Local receipts do not imply publication. Empty GitHub status/check lists are
absence of CI evidence, not a CI pass.


## Final integrated gate and portable layout

`FINAL_REPLAY_RECEIPT.json` contains fresh publisher normal/optimized runs of thirteen
entrypoints in temporary copies with `PATH=/nonexistent`, plus both stronger
named-commit core-triangle audits. All stdout equals the frozen saved receipt
bytes. `ROOT_FINAL_PREFLIGHT.json` is an additional independent 16-run preflight.
No frozen producer/auditor source or native checksum manifest was rewritten.

From this packet, run `python3 replay_final.py` (or add `--repo /path/to/math`
for direct checks of the inherited commit). Python uses only its standard
library; the optional commit check additionally calls Git. The original C++
local lane has its independent two-build replay in the first-checkpoint receipt.

Our sibling lane names include `session67-`. When invoking portable audit
scripts by hand from `session67-cycle-audit/`, pass
`--metagraph-dir ../session67-metagraph` or
`--recursion-dir ../session67-recursion` explicitly. Frozen generic examples
in the audit addendum use shorter sibling names and are not silently edited.
The default portable triangle audit has 10,657 checks and does NOT authenticate
the named commit; explicit `--repo` restores seven direct source comparisons
for 10,664 checks. See the addendum's precise provenance distinction.

All 73 frozen historical source files are compared directly with their declared
repository commits in the integration gate. `HISTORICAL_CORRECTION_RECEIPT.json`
records minimal status/banner changes with original bodies retained. Historical
source snapshots and exploratory scripts are evidence, not endorsed executable
current truth. A native manifest's producer path is provenance, not an execution
dependency; the portable replay entrypoints use the published packet layout.


The incoming integration has an independent supplied-data arithmetic and
counterexample verifier, a complete n<=5 labelled score-fibre checker, and a
second root-written order-12 automorphism witness, all replayed in both modes.
`ORIGINAL_AUDIT_REPLAY_RECEIPT.json` adds four executions of the unchanged
nonportable original audit sources in their original producer layout; portable
published entrypoints remain the default. `INCOMING_CORRECTION_RECEIPT.json`
records the source-preserving current-scope overlays and native checksum refresh.


Concurrent main was integrated at `18866090755d6b1cc6047caffcff8ff358ffcf3f`.
The pre-merge correction receipt retains its original before/after hashes and
adds `merged_current_sha256` for the current source files. Incoming sharp-4p,
path-cover and necklace proofs, source-script fixes and reservations are preserved;
no new Session68 research was incorporated into this Session67 packet.

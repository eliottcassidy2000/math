# Verification protocol and scope

Every checkpoint is an explicit file inventory. `verify_packet.py --replay`
checks its complete file set, SHA256 and byte counts, Unix modes, native-copy
ledger, evidence bindings and paired normal/optimized receipts. Replay runs in
an isolated temporary packet and must not change any frozen input or output.
The explicitly approved historical decoder replay remains inherited research,
not a new Session74 theorem. All other pinned historical programs are inert.

`verify_provenance.py` separately checks exact Git blob hashes and source pins,
source modes, native producer hash manifests, and auditor-to-producer bindings.
Its mutation tests are independent of the outer file inventory: relabelling a
source, corrupting a declared digest, mode, pin or source blob must fail. Both
packet and provenance rejection tests run normally and under `python3 -O`.

Finite enumerations remain FINITE-EXACT, not universal theorem proofs. The
all-size score argument and exact six-step distance argument are reviewed as
mathematics; independently generated finite controls and the six named paths
are machine checked. No LRC(14), Collatz integrality or universal decoder claim
is inferred from these computations.

The checkout is sparse. The documentation gate makes targeted exact reads for
startup budgets, maintained warning sentinels, truth tokens and 33 status
fixtures. It explicitly does **not** claim the repository-wide theorem
dependency audit or router smoke. Reader-facing Markdown links are checked
for local destination existence; source snapshots and external sites are not
misrepresented as newly re-audited content.

Remote publication uses non-forced expected-head updates. The coordinator
checks the fetched exact commit, parent, whole tree, changed file set, bytes,
modes, required ancestry and final clean local 0/0 state. CI observations are
bound to that exact SHA. An empty check/status/workflow set is **not** a CI pass.
Per-checkpoint results are reported in additive verification receipts.

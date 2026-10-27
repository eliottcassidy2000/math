# Verification scope

This packet binds an explicit curated inventory, not the entire repository.
Each native producer and auditor file is copied without editing its bytes or
mode. Repository snapshots bind exact commit, repository path, Git blob,
SHA256, byte size and mode. The active provenance check enforces source hashes,
Git blob identities, approved pins, native hash lists and reviewed-file bindings.
The checkpoint builder additionally compares each declared repository origin
with targeted `git show` and `git ls-tree` reads at the recorded commit.

`verify_packet.py --replay` copies the whole curated packet to a temporary
location and runs only declared current jobs with normal Python and `-O`.
Each output must equal the corresponding frozen receipt byte-for-byte; the
checker also rejects changed, removed or newly created packet inputs. Sibling
imports, where needed, are enabled only for the named job inside that temporary
packet. Every inherited source program remains inert unless explicitly listed
as an approved source replay.

Independent mathematical audits are separate from transport-integrity gates.
Hostile certificate cases change actual minimums, interior subset optima,
quotient orientation and vertex naming. Packet and provenance mutation suites
independently reject altered bytes, sizes, pins, modes, source blob IDs,
unlisted paths, unsafe paths, symlinks and incorrect reviewed hashes under both
normal and optimized Python. A hash manifest is not a proof of its claims.

Documentation checks are narrowly scoped: bounded maintained docs, byte-zero
sentinels, explicit requirement strings, startup budget and selected current
status fixtures. Whole-corpus dependency checks and the unproved-router smoke
check are **NOT_RUN** by this scoped gate. Link checks cover only declared
reader-facing local destinations, not all historical source links, heading
fragments or external URLs. No full archive or broad partial-clone grep is used.

Remote publication verification is recorded outside this immutable research
packet: exact remote parent and tree, complete tree/object identity, changed
paths, packet bytes/modes, required ancestry and clean local alignment. CI is
observed for the exact published commit; an empty check/workflow list is not a
passing CI result. No packet-local receipt by itself asserts remote publication.

# Verification and publication scope

Mathematical proof, finite-exact enumeration and publication integrity are
separate evidence layers. Byte checks do not establish an unstated theorem.

`FIRST_SOURCE_VERIFICATION.json` pins every approved source and audit copy.
The original producer manifest is retained inside `session64-repairs/`;
no audited proof or checker bytes are edited for packaging. The independently
frozen ordinary and optimized replays must agree exactly with saved outputs.
The primary scripts and inherited Python dependencies are inspected for
optimization-removed assertion gates. Native checksums, local links, UTF-8,
file sets and exact source bytes are checked separately.

Publication uses explicit owned paths, a one-parent local checkpoint, exact
Git blob/mode/tree pins, a non-forced expected-head update, and an independent
Git fetch. Remote receipts verify whole tree, changed path set, parent, every
changed blob and mode, and reachability from current main. Missing CI is not
reported as passing CI. Original dirty work and all incoming commits remain
untouched; existing startup warnings and line endings are preserved.

## First checkpoint evidence

`FIRST_REPLAY_RECEIPT.json` records two mathematical checkers and four successful
normal/optimized executions from frozen copies. Outputs match both each other
and the saved receipt bytes exactly, and the original packet is unchanged.
Every replayed script and approved source hash was rechecked before release.
`FIRST_DOCS_CHECK.txt` passes the inherited documentation gate.
`FIRST_ROUTE_CHECK.json` records preserved startup line counts and CRLF.
`FIRST_BUNDLE_CHECK.json` covers UTF-8, five assertion-safe scripts, native
checksums and local links. `FIRST_FILESET_NEGATIVE_CONTROLS.json` records
ordinary/optimized rejection of byte mutation and an unexpected file.
Exact-commit remote verification and CI visibility follow publication.

## Second checkpoint evidence

`SECOND_REPLAY_RECEIPT.json` records 10 mathematical/provenance scripts and 20 successful ordinary/optimized executions. Paired outputs match exactly, all designated saved outputs agree, and the original frozen input packet remains unchanged. Every replayed script and current approved source hash was rechecked. `SECOND_DOCS_CHECK.txt` passes the documentation gate. `SECOND_BUNDLE_CHECK.json` verifies assertion-safe Python, native checksum pins, UTF-8 and local or explicitly pinned original-context links. `SECOND_FILESET_NEGATIVE_CONTROLS.json` rejects byte mutation and extra-file mutants in both Python modes. Exact remote tree/blob/mode verification and CI visibility follow the publication operation.

## Final checkpoint evidence

`FINAL_REPLAY_RECEIPT.json` records 20 mathematical/provenance scripts and 40 successful ordinary/optimized executions. Paired outputs match exactly, all designated saved outputs agree, and the original frozen input packet remains unchanged. Every replayed script and current approved source hash was rechecked. `FINAL_DOCS_CHECK.txt` passes the documentation gate. `FINAL_BUNDLE_CHECK.json` verifies assertion-safe Python, native checksum pins, UTF-8 and local or explicitly pinned original-context links. `FINAL_FILESET_NEGATIVE_CONTROLS.json` rejects byte mutation and extra-file mutants in both Python modes. Exact remote tree/blob/mode verification and CI visibility follow the publication operation.

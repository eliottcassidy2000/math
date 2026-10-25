# First checkpoint verification

Status: package checks passed; mathematical scope is defined by the individual proofs and independent audits. Publication status is established by the Git commit containing this file, not by this document alone.

The first checkpoint includes five research packets (toric shared roots, signed flags, arithmetic configurations, cyclic-unit support dynamics, and graph-magic degeneration), specialist positive-geometry and cyclotomic audits, and parent independent controls. Later resolution and full nine-cell arrangement extensions are excluded from this checkpoint.

Each of the nine primary producer/audit checkers was executed normally and with `python3 -O` in an isolated frozen copy. All eighteen executions passed, and outputs agree between modes. Duplicate before/after cyclotomic replay snapshots are audit evidence, not additional independent algorithms. The toric-shared-roots producer uses assertions, so its optimized run verifies execution/output consistency only; its assertion gates are established by the normal run. The independent positive-audit checker has always-active exception gates and covers the geometric/deformation claims in both modes. All other first-checkpoint checkers use always-active checks.

Reproduce each packet with its displayed command. Run `python3 verify_packet.py` from any working directory for an exact file-set, byte-size and SHA-256 check. The manifest excludes itself to avoid a circular hash. Python cache files are ignored; no caches are published.

Source-pinning files distinguish previously published repository material from new work. The positive-audit manifest records both the original audit capture paths and portable relative paths. The signed-flags checksum list was normalized to packet-relative filenames without changing the mathematical source. No downloaded full article or private attachment is included.

The bounded repository startup/documentation checker passes without loosening its limits. No CI result is implied by these local checks.

## Second checkpoint

Five new producer/audit scripts are replayed normally and with optimization, with the exact results recorded in SECOND_REPLAY_RECEIPT.json. The geometry resolution producer uses assertions: its normal run validates those gates, while its optimized run checks output consistency. The independent resolution auditor uses always-active exceptions in both modes. The branch auditor’s source-certificate lookup was made relative to the sibling graph-magic-degeneration packet; only this location was changed. Frozen source hashes are retained in its audit, while the packaged script receives a fresh checksum. No third-wave hierarchy or distinctness claim is included in this checkpoint.

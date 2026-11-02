# Session137 source-integrity reconciliation

**PASS.** The final lane inventory contains 233 files. Nine frozen hash manifests contain 242 entries, all matching both the original files and the portable copies after replay. All 233 original files remained unchanged throughout the expanded replay.

## Reviewed and inherited source identity

The [machine-readable integrity ledger](SOURCE_INTEGRITY.json) retains the checks of:

- 13 independently reviewed manuscript hashes and byte counts;
- 10 carrier source SHA256, Git-blob, and baseline pins;
- 6 constructor inherited pins;
- 5 harmonic inherited copies against their named baseline sources.

All 34 match. The independent directory's 89-entry payload selector omits its own selector and hash manifest; both are included, so that lane contributes 91 files. The weighted-stability manuscript was added afterward and has its [separate final audit](WEIGHTED_MANUSCRIPT_AUDIT.md).

## Inventory boundaries

The seven newly integrated lanes are width, width_phase2, constructors, coverage, cycle_carriers, harmonic_owners, and independent. The first-checkpoint discriminant proof, global referee, and independent_global paths are preserved. Duplication of global receipts between the old and complete independent packets is intentional and non-destructive.

`width_phase2/probe.json` remains explicitly exploratory and is not a theorem dependency. Stale pre-closure drafts, compiled binaries, caches, speculative constructor outputs, and transport files are excluded. The root aggregate manifest includes the final public prose and verification metadata as well as the frozen lane files.

## Replay scope

The initial 26 executions passed. The [expanded replay](FULL_REPLAY.md) is terminal PASS: 65 scientific executions, nine C++ builds, and twelve sidecar comparisons. Its full k=7 run repeats the producer's source, not an independently implemented census. All-order scalar zero-rigidity, full exact-spectrum completeness, and LRC(14) remain open.

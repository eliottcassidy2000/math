# Replay the fourth Session190 checkpoint

Restore all available ARCHIVES_CP*.json manifests in ascending order with the current RESTORE_ARCHIVES.py, following RECEIPTS.md. Then copy the complete squareports or cyclic subtree to a temporary directory before executing scripts or redirecting their outputs.

The squareports source and audit preserve their original sibling directory names. Its FROZEN_AUDIT_MANIFEST.json is the complete 15-file source-and-audit allowlist including itself. In that temporary subtree, run the source verify.py and verify_structural.py, and the audit independent_check.py, normally and with -O as described in AUDIT.md. Compare their complete JSON outputs to the frozen originals. Prototype/exploration scripts, caches and duplicate local replays are intentionally not release dependencies.

The cyclic packet contains exactly NOTE.md, check.py, results.json, results_optimized.json, and AUDIT.md. Run check.py normally and with -O in a temporary copy and compare outputs with the retained receipts. Its claims concern the specific embedded C15 subgroup, a precisely scoped equivariant-map obstruction and a type-uniform stochastic bound, not arbitrary unlabelled bijections or a higher-degree census.

Historical sibling references in the squareports source point to earlier workspace packets; the corresponding reviewed Session190 planar packet is available at planar/session190-planar relative to this result root. No earlier producer is imported by the release replays except the squareports transfer functions explicitly under test. The original proofs and their source attributions remain unchanged.

# Balanced order-ten blue wiring certificate

Self-contained, Python 3 standard-library producer packet. No network, external
repository, or third-party package is needed. Read THEOREM.md for definitions,
proof, exact finite scope and inherited completeness assumptions.

Result: 229 self-converse balanced classes, 41,699 blue ports, 11,452 internal
edges (24 loops), one connected component, no bridges or articulation nodes,
and 18,795 boundary edges to 786 classes. The 20 positive-error balanced nodes
form an already-connected dominating core: all 209 zero-error nodes touch it.

Run from any directory:

    python -B /path/to/packet/check.py
    python -O -B /path/to/packet/check.py --output /tmp/blue-optimized.json
    cmp /path/to/packet/normal.json /tmp/blue-optimized.json
    python -B /path/to/packet/verify_certificate.py
    python -O -B /path/to/packet/verify_certificate.py --output /tmp/blue-certificate-O.json
    python -B /path/to/packet/topology.py
    python -B /path/to/packet/verify_manifest.py

The final manifest check must follow reproducible regeneration. check.py uses
explicit RuntimeError checks, so -O does not disable validation. It suppresses
bytecode generation even without -B. The producer visits 933,120 halfpath
candidates and an independent 1,048,576-word reflection cube, never the full
68,719,476,736-word unrestricted carrier. Replay is ordinarily short, but no
wall-time bound is claimed.

## Files and certificate semantics

- check.py: complete bounded producer, exact normal forms and frozen controls.
- normal.json: complete records and computed invariants; normal.out: concise log.
- optimized.out and replay.json: optimized regeneration log and byte equality.
- cost.json: independently obtainable with check.py --cost-only.
- verify_certificate.py and certificate-check.json: literal record replay and
  mutation/topology-fixture checks; certificate-check.out: replay log.
- topology.py and topology.json: core domination and explicit spanning tree.
- SOURCE_PINS.json and sources/: exact inherited dependency closure for replay.
- MANIFEST.json and verify_manifest.py: hashes of all packet payloads; no caches.

Atlas node IDs are the frozen atlas row indices, not consecutive 0..228.
Boundary class IDs are separate, consecutive 0..785. A port row contains
[word, source_atlas_index, source_path, C_word, target_kind, target_class_index,
target_path_or_canonical_order]. For an internal target, its order maps the
target representative to the literal target matrix. For a boundary target,
its order maps the literal target matrix to the boundary canonical matrix.
This directional distinction is intentional and explicitly replayed.

Input theorems/manifests are provenance snapshots, not instructions to rerun
those historical packets. Only balanced-metagraph-engine.py is imported by the
producer; its dependencies are standard-library modules. Copied historical
check.py is provenance only, not a runtime dependency. Frozen manifests may
list other historical artifacts that are not inputs to this packet; those
unreferenced artifacts are not required for this replay.

# Complete blue edge-neighborhood of24 defect nodes

The4020 boundary blue edges meet1031 distinct zero-error self-converse nodes.
Together with the24 defect nodes and their90 internal edges, this gives1055
vertices and4110 blue edges touching the defect set. Other zero/zero edges and
all black edges are outside scope.

Boundary classes:1016 discrete-stable rigid,11 with3+3+1+1+1+1 stable cells,
and4 with5+5. Exact normalization uses intrinsic directed cycles and complete
crossing matrices, at most25 coordinate rotations per target, not10! search.
`normal.json` contains all target classes, boundary adjacency, and labelled
canonicalization witnesses for each edge. See `THEOREM.md` for completeness.

Reproduce:

    python3 check.py --output normal.json
    python3 -O check.py --output optimized.json
    cmp normal.json optimized.json
    python3 verify_manifest.py

The earlier converse and induced-wiring packets remain immutable. This is the
final boundary closure, not a full metagraph recovery.

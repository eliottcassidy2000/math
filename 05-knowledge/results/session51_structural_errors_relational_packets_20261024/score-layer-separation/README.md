# Score-layer neighbor separation

For even-order tournaments with two outdegrees of equal multiplicity, distinct
score layers have disjoint neighbor cyclic-triple intervals under native
retained-Hamiltonian-path reversal. Within the union of these layers, only the
balanced layer can contain an edge. The proof uses inherited endpoint transport
and no primality or self-converse assumption. Blue paths satisfy a sharper
two-value necessary list; arbitrary paths have three allowed values.

This is a scoped all-size wiring obstruction, not a connectivity claim.
Read `THEOREM.md`. The finite checker includes the formal-negative-value
hostile to claiming all allowed values are attained.

    python3 check.py --output normal.json
    python3 -O check.py --output optimized.json
    cmp normal.json optimized.json
    python3 verify_manifest.py

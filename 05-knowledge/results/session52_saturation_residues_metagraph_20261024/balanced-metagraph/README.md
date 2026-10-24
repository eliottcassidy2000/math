# Balanced-layer metagraph counts from n alone

Exact degree-coefficient algorithms count the balanced source words and their
internal/boundary blue and black ports, without a tournament-class atlas.
A regular-cycle lift identifies internal words with Eulerian orientations of
K_(n+1) minus its Hamilton cycle, giving an independent DP and rotation
congruences. Read `THEOREM.md` for proof, conventions, and the loss ledger.

At n=12 there are 10,188,802,555,033 balanced source ports and
1,606,912,987,280 internal C-lines, of which 3,401,824 are blue.
These are exact line multiplicities, not A000568 values or reconstructed wiring.

Scalar input:

    python3 engine.py 12

Portable controls:

    python3 check.py --output normal.json
    python3 -O check.py --output optimized.json
    cmp normal.json optimized.json
    python3 verify_manifest.py

Only Python's standard library is required. Full literal cubes are checked
through n=8; ordinary class/Hamilton-path controls through n=6; both independent
all-size count engines through n=12. n=2 remains a separate nonfree boundary.

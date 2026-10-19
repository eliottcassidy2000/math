# A graceful fixed-color counterexample

Packet 102. Explicit forest on labels 0,...,12 with free colors 6,9,11. The entire five-state completion fiber has width-two graph K3 disjoint K2. Thus arbitrary fixed-color fibers need not be two-color connected. In the global two-color graph, crossing between these components necessarily disturbs an initially correct frozen edge.

No minimality claim or exponential lower bound is made. The earlier upper-tail theorem remains valid.

Replay using only standard-library Python, from this directory:

    python verify_m12_witness.py > /tmp/counterexample102.json
    cmp m12_exact_receipt.json /tmp/counterexample102.json
    python -O verify_m12_witness.py > /tmp/counterexample102O.json
    cmp m12_exact_receipt.json /tmp/counterexample102O.json
    python independent_m12_audit.py > /tmp/independent102.json
    cmp independent_m12_audit_receipt.json /tmp/independent102.json
    sha256sum -c SHA256SUMS

The verifier checks all 56 candidate free-edge tuples. `witness.json` contains the full fixed forest, free colors, and discovery component map; correctness depends only on the explicit edges and exact enumeration, not on a CSP solver or negative search results. All runtime dependencies are co-located.

The independent audit uses assertions and must run without -O. Its hard-coded witness, DSU cycle test, full distance matrix, and exact simple-graph edge list provide a second implementation.

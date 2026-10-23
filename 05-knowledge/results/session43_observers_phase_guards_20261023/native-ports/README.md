# Session43: actual Hamilton observer and port packet

Read `OBSERVER_BUNDLES_AND_NATIVE_PORTS.md` for the proof and exact scope.

Main deliverables:
- Observer fibre q(T)=H(T)/|Aut(T)|; the actual path order cannot in general be chosen equivariantly.
- Normalized converse is reflection, whereas free-arc flip is complement; exact blue/black port budgets follow.
- Separate genuine-path transport through already-complete isomorphisms from assembling paths with independent internal charts.
- An almost-regular six-vertex native assembly family requires two cosets, not one; its six paths are a strict subset of the host's 45 Hamilton paths.
- A positive all-size family has p product |D_e| literal paths. Interleaved templates require cycle-compatible differences and retain all phase/kernel lifts.
- Disjoint local coset branches give an exact branch-aware sufficient compiler; no all-size bound on branch complexity is claimed.
- Actual mod8 guards remain on literal lifts unless quotient invariance is proved. The included mod8 hostile is explicitly synthetic, not an alleged electrical invariant.

Reproduce from repository root after copying this packet directory:

    python check.py > receipt.json
    python -O check.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

Run those commands from this directory. Tests use no Python assertions and no dependencies. Receipts cover all observer masks at n=3,4,5, all six-vertex host paths, and 72 seeded actual nine-vertex shifted templates. They do not claim an order-eight census.

This packet does not perform Git writes. Its base sources and all deliverable hashes are recorded in `MANIFEST.json`.

# Independent review: two-apex reverse bridge

**PASS, no required corrections.** No producer edits or remote writes.

Bound to:

- PROOF.md: 7b6224b661fdcb1400aef9e9a15901a4a4ef70323ab58649e78484f0f9630d9c
- verify.py: 35948ef76ce55d76359debf6b572db5b982feec6c527c2e647d68c680737dc20
- finite_certificate.json: d4f3d5a98230955d0b82701913c054a04c7a0c0958d0c6f96df7b64886b224dd
- receipt.json and receipt.optimized.json: 122eabaf100aab48f78ee818bbc777ebb41b34203db318a8944a9a49f2eec1a9

## Mathematical checks

- Mixed directed cycles have precisely one reverse-cut crossing, y→x, and exactly one forward crossing. Rotating to that reverse arc gives a unique pair of rooted increasing module paths. The forbidden singleton pair has even total vertex count. The signed parity count for each unpadded strong module is −1 (including the singleton), so the odd pair count is (PQ−1)/2.
- The x-only and y-only cycles stay entirely within their own modules. Each old family forms a clique, these cliques have no cross-conflicts, and every mixed cycle meets all others and both old families. This proves the graph isomorphism itself, not just its independence polynomial.
- Strongness follows via y→x and the forward cut together with the old strong cores. Deleting both apices leaves a chain; two disjoint old triangles prove a single deletion insufficient.
- Each old x-rooted path independently selects a subset of the new sink strip. The factor 2^k multiplies the count and for k≥1 the signed factor is (1−1)^k=0. Thus quadratic coefficient ab is preserved; the +2 exceptional step from k=0 is necessary.
- The deficit D=A+1−h starts at 1 with the first 1, loses one only at a 1→0 transition, and each later 1 increases it by at least 3; a loss to zero is possible exactly for 10^r at strong endpoints. For the upper inequality the prefix before the last 0 gives h≥A(prefix)+2, with equality exactly for an all-1 prefix. Reverse-complement exchanges the ports. Both endpoint equalities are attainable as stated; the bridge minimum and C3 residue limitations follow.
- Each strong word of length m has the outer endpoint inversion contributing at least 2^(m−2) to h−1. Leading zero/trailing one stripping preserves h for all weak one-apex presentations. For a bridge with lengths m,n, the endpoint product exceeds 2^(m+n−2), yielding a safe total-length cutoff; padding increases H, so the same cutoff applies with arbitrary padding. The compiler does not exclude an admissible module length and each padding loop is monotonically complete.
- Generic feedback-set coefficient/discriminant claims check. The census reaches all orders *within the defined families* at the finite H cutoff, not all FVS-2 tournaments or all odd targets. Prior broad modular completeness is correctly distinguished from exact integer coverage.

## Exact controls

Re-ran the producer verifier under normal Python and Python -O in this audit directory. The receipts agree byte-for-byte with each other and the producer receipt, and the generated certificate is byte-identical to the producer certificate.

The clean-room checker verify_independent.py computes each one-apex value by its literal inversion-endpoint cycle weights rather than the producer recurrence. It enumerates by total length, bounds padding by bit length, checks equality of both full attained sets with the producer certificate, and validates every supplied witness. Normal and -O receipts agree:

    {"bridge": 4004, "one_apex": 2595, "presentations": 20181, "repairs": 1941, "union": 4536}

Independent checker SHA256: 995b7e4d4db646a0eab90772b7fb59c95bedc90cc015c90b9ce175e1710f04e6
Independent receipt SHA256: beff00b553c0c1765e610d3ab5c44aab38f9e0d6cf262ce286788d60f17ad33a

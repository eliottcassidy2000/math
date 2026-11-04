# Second endpoint transfer

For any quotient Q and anchor i with even endpoint skew, sigma_i is defined by Hamilton words and words that repeat one quotient name twice. For arbitrary nontrivial children, the physical block query `(e_i,0)` has complete augmented-profile contrast `(2sigma_i,0,0,...)` modulo eight. A value sigma_i=2 separates in every fixed ambient. This is a sufficient criterion, not reconstruction or a collision test.

The general proof is combinatorial: a fixed-point-free return swap removes three-run terms modulo four, exact uniform-start tomography recovers the anchor, and a forced x-to-y Hamilton-cycle bijection removes the higher-profile obstruction. It resolves the endpoint-blind alternating prime-five class for arbitrary nontrivial children and explains the independently certified prime-eight anchor.

Analytical and numerical independent reviews passed without correction. The separate `second_endpoint_transfer_audit` reproduces all quotient rows, every retained physical coordinate and all return-swap words. The mathematical sources and paired receipts are frozen.

Run `python3 verify.py > normal.json` and `python3 -O verify.py > optimized.json`, then compare. This uses only the standard library and an installed C++17 compiler. The exact receipts include all 1,099 quotient controls through order five, the complete indicated physical-query controls and source/dependency pins.

Frozen source allowlist: `PROOF.md`, `README.md`, `verify.py`, `physical.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`. The exploratory `explore.py` is not part of the packet.

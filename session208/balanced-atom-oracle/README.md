# Session208: A balanced order13 atom and an exact flip oracle

The explicit strong order13 tournament has H=759375=15^5. Dominance joins improve the asymptotic global upper bound to13/5. An exact complementary-subset endpoint convolution computes all single-arc-flip count changes in O(n^2 2^n) integer arithmetic operations. The search graph is not the native path-complement metagraph. Neither global optimality nor exact minimum order is claimed.

Replay the main checks with `python3 verify.py`, with assertions enabled. Output should match `verify.out`; RESULT.json is deterministically regenerated. The standalone independent witness checker and exhaustive small-order arc-difference checker are included with archived outputs.

Optional discovery replay:

    g++ -O3 -std=c++17 target_search.cpp -o /tmp/s208-target
    /tmp/s208-target 13 300 759375 > /tmp/target13.out 2> /tmp/target13.log

This uses seed20813 and finds the witness at zero-based restart70, step3. Every chosen arc difference is verified against a fresh DP recomputation. Discovery convergence is not a theorem. The explicit witness, independent counting and proved SCC product rule establish the construction.

Only the target15^5/order13 experiment belongs to this release. No order14 result or search outcome is claimed here.

Run `python3 independent_witness.py` and `python3 independent_delta.py`; stdout should match the corresponding `.out` files. `DEPENDENCIES.json` pins prior publications; `REVIEW.md` records independent review. Validate the frozen payloads with `sha256sum -c SHA256SUMS`.

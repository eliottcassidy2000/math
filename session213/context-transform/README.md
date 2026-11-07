# Complete cycle codes do not give a total-path tail

This packet proves that Session212's full binary-comb family, including both apex orientations and every leading-zero padding, has only O(X^(log_3 2)) distinct total Hamiltonian-path counts at most X. Its Hamiltonian-cycle counts nonetheless realize every nonnegative integer.

- `PROOF.md`: explicit disjoint path classes, the zero-density theorem, its finite-menu induced-extension corollary, the exact sole-exit total-path readout, and the deletion readout 2k+g(k−1).
- `verify.py`: self-contained exact verifier, no external packages or imported research implementation.
- `RESULT.json`, `RESULT.optimized.json`: byte-identical normal/optimized receipts, 1,016 direct-DP cases.
- `SOURCES.json`: exact inherited snapshot and dependency identities.

Reproduce from this directory:

    python3 verify.py > /tmp/context-transform.normal.json
    python3 -O verify.py > /tmp/context-transform.optimized.json
    cmp RESULT.json /tmp/context-transform.normal.json
    cmp RESULT.json /tmp/context-transform.optimized.json

The independent Session213 audit reviews the all-order proof and reruns these frozen artifacts. This is a construction-specific barrier, not a refutation of the unrestricted odd-path spectrum conjecture. In particular it does not prohibit deletion, rebuilding old arcs, unboundedly many outputs per input, or a different family.

The all-one 2*3^m+2^(m+2)−1 formula, the cycle constructor, the one-apex g recurrence, and insertion monotonicity are inherited and credited in the proof. The new connection is that explicit ternary path allocations make the complete binary cycle constructor and every fixed finite menu of its induced extensions too sparse to contain an odd tail.

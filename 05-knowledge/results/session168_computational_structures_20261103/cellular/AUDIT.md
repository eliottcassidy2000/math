# Independent audit

Reviewed `NOTE.md`, `check.py`, and `certificate.json` on 2026-11-03. The seed, literal product cover, deck-equivariant reversible gate lifts, p-or-2p orbit criterion, clean five-layer compiler, and exact logical/physical 8·3^e clock proof are sound. Scope correctly distinguishes retained labels/layer schedule from intrinsic autonomous dynamics and single arc flips from synchronous rounds or complement-line metagraph moves.

Normal and optimized verifier runs both passed and reproduced the existing certificate byte-for-byte:

- Certificate SHA-256: `a0ac7bc789aac118ff6b78a84897fa50ae3d683ece547a7e4b99268f05b0d51d`.
- Audited note SHA-256: `06964519ababf7d53d3d459f129d1c1f058c9b1c77fd6ea324ada558773f8e9f`.
- Audited checker SHA-256: `829607a94be5336f99f56d165cfccd363ab587cc18b13fa92f92101b12f980cd`.

`independent_audit_check.py` imports neither the candidate checker nor its output. It independently enumerates the seed Hamilton paths and changed arcs, recomputes every listed logical first return, and tests termwise signed-increment cancellation on 24 deterministic random physical starts at widths 1, 3, 9, 27, 81, and 243. It also checks each tested physical microgate is an actual single-arc transition or identity. These are controls supporting, not replacing, the general proofs.

One requested premise clarification: the irreversibility argument must explicitly assume that E is an injective encoding. Assigning one codeword per state, by itself, would allow a constant encoding. An optional clarification of the class-graph argument is that reversing a cross-block arc merges consecutive SCCs, while splitting an internal block increases their count; neither creates an extra edge between displayed class words differing in multiple coordinates.

Reproduce the independent controls with:

    python3 independent_audit_check.py > independent_audit_results.json
    python3 -O independent_audit_check.py > independent_audit_results.optimized.json
    cmp independent_audit_results.json independent_audit_results.optimized.json

The inherited strong-wrapper theorem was not independently re-proved here; the core cover, compiler, and clock are self-contained.

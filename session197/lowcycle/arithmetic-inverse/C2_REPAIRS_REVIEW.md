# Independent verification of the 26 C2 arithmetic repairs

**PASS.** This review binds `session197-lowcycle/arithmetic-c2-witnesses.json`, SHA256 `c02dcf59ed38ecda740af73fb2b5755b6ce12c63e2a50a11a4a0e586a380c32c`.

The exact target set is

    13,41,43,55,77,83,163,173,231,233,247,251,277,299,
    421,449,503,533,569,583,637,821,1067,1339,1861,2647.

I independently reconstructed all 26 graphs from the displayed source, entry pair, pair-gadget index, and subsequent ordinary-gadget indices. Every graph has the claimed Hamilton-path count, exactly two Hamilton cycles, and exactly two paths starting at its marked root. The maximum order is 13. Every source has exactly one cycle and both selected source vertices have start count one. Direct path-cover inventories reproduce every supplied (h,E,R,K) response.

Two counting methods agree:

1. Exact subset endpoint dynamic programming, also checking the full rooted (H,e,R) response and the pair/single-entry arithmetic formulas.
2. A separately written walk inclusion-exclusion counter. For every vertex subset S it counts all length-(n−1) walks in S and uses alternating subset sums to retain walks visiting every vertex. A parallel fixed-start closed-walk count gives the Hamilton-cycle count. This uses no Hamilton subset-endpoint recurrence. At orders at most 13, the explicit signed 64-bit bounds are sufficient: even summing the maximal n(n−1)^(n−1) walk bound over all 2^n subsets stays below 2^63.

The walk checker passes in optimized and assertion-enabled unoptimized C++, with identical output. The Python reconstruction passes in normal and optimized modes with identical output.

I additionally compiled every actual resulting graph using the independently reviewed start-two physical compiler, preserving its marked root. All adjacency entries, initial edge normalization, and retained ports check. `physical-words.json` gives the complete rows, original-label birth order, twelve-letter mask/drop word, final port list, and marked root for all 26 witnesses. These are actual constructions in the two-changing-port language, not only abstract scalar recipes.

This review establishes these 26 positive witnesses. It does not independently certify the larger packed finite-coverage images or the all-order response-cone theorem in the surrounding arithmetic packet. Those are separate review components. The experimental pair-entry formula is not needed to justify the counts, since both independent graph counters verify them directly.

## Independent artifacts

All files below are under `session197-hybrid-constructions/audits/c2-repairs/`:

- `check.py`: `348670f0feecd4e9816669fa02ed68f6abe63c8d192b2b01b0f074880e51cf7a`.
- `check.json`: `13a9bd443a2f820bc9e4d9bcb6145d1c182c2724647619c6a2995c34b9b85088`.
- `walk.cpp`: `e45506ec026ef7e2fa2254e203955484d2cbed230bfce03d5c7a2be736df1e10`.
- `walk.json`: `62e3a2ff909d489a1f643a998c4cc7331dbc16573df89e78bd64f7b2858d1cd6`.
- `physical-words.json`: `103f8c358f349f338cc3e726264ed1c1e5557bf137d7935d3b152313eba6e90e`.

The frozen source recipe file was left unchanged. No correction is requested.

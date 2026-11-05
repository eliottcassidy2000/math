# Independent audit: complete small paired-gate catalog

**PASS, no correction required.** The119-value catalog below73,529,109 is complete within P(T)=C3[T,T^op,Tr2]. Its8 prime outputs and3 seven-prime outputs are exactly as stated. The least7p output in this construction family is37,050,461. This is not a minimum among arbitrary tournaments.

## Frozen source binding

Source directory: `session197-odd-cycle-path-polynomials/`.

- PAIRED_GATE_CATALOG.md: `fb9be988cfc42bb2768e26576bb925ca5c7c9681173626b5bca4c93827e9e947`.
- verify_paired_gate_catalog.cpp: `4364be58afb33c75bf3f701c1b2fa19d7ee0954914827472e6e37271caddc99d`.
- paired_gate_catalog_normal.json and debug.json: `9ffcd678e33279dc8ad4ce41331a10190431baa7d707c6466f07a8a560c267ac`.
- verify_paired_gate_catalog_witnesses.py: `4ea154abf46f9e670cf7e5f244783f3fdb96a8bbf720c918a5dc18b38f6186e2`.
- Witness normal and optimized receipts: `70bc118d4eaa03fab016a45ae8494d623d45c80e940e59192bfbdb4726b11c6b`.
- Computational helper verify.py: `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`.

All source replays were run in a separate copy. The C++ optimized/debug receipts and Python normal/optimized receipts reproduce the frozen bytes. Original source files were not changed.

The all-order gate dependency is THREE_COVER_CYCLIC_GATE.md, SHA256 `0630c4ff7f9a52e69cbce4f18cdac8678609d4565a3e42f0afec8d588afe37cf`. Its separate independent PASS report is `session197-three-cover-gate-audit/REPORT.md`, SHA256 `3170f5dc7daaa3ae6f96528783c2c780991c25cafa4fe73e88be56ac71273763`, published with that theorem at commit073d1e643d504557abf38286e129870dca4a32a1. This catalog review uses that proved physical kernel and sharp order bound; it does not replace the underlying all-order proof audit.

## Completeness and cutoff review

1. The inherited paired kernel is the nonnegative quadratic form5h²+12hu+28u²+48uv+72v². The partition lower bounds c_j(T)>=S(n,j) give Phi_n, attained by a transitive input. Phi_n increases strictly, and Phi_8=73,529,109. Thus the strict catalog cutoff excludes every input order at least8, including any weak input or nontransitive input. No strong-input assumption is silently imposed.
2. Every tournament on n<=7 admits a labeling whose outdegrees are nondecreasing. The source loops over every labelled orientation before applying precisely this filter. Ties retain duplicates rather than discard all representatives. The stated degree-sorted counts are counts of retained labels, not isomorphism classes. This is sufficient for completeness of the set of values.
3. The subset endpoint DP is processed in increasing masks, so every predecessor and every proper subset count is already available. The two-cover sum counts both choices of its first block and divides by2. In the three-cover sum, requiring the selected block to contain the least vertex selects exactly one block of each unordered partition. The empty residual contributes zero, as required for exactly three nonempty paths.
4. Integer bounds are safe in the source's signed64-bit arithmetic. Each kernel summand counts a disjoint subfamily of Hamilton paths of an actual graph of order at most16, hence is at most16!, and their sum is at most16!. The subset and cover terms are smaller. The orientation bit shift has exponent at most21. No floating-point or modular count is used.
5. Every smaller candidate than37,050,461 is strictly below Phi_8, so it would appear in the enumerated image if attainable by this gate. Exact trial division on every image value and every applicable quotient identifies the three7p values and their least member. A larger-order paired graph cannot invalidate this minimum. No corresponding unrestricted-spectrum or least-order claim follows.

## Independent enumeration and counters

`independent.cpp` imports no source code and uses a different complete enumeration:

- Starting from the singleton, extend every isomorphism representative by every new-vertex incidence mask.
- Canonicalize by exhaustive permutations inside outdegree classes, without iterative color refinement and without the source's global sorted-degree filtering loop.
- Independently compute each complete automorphism group and sum n!/|Aut|. The class inventories are1,1,2,4,12,56,456 through order7. Their labelled orbit masses are1,2,8,64,1024,32768,2097152, exactly2^(n choose2) at every order.
- Count Hamilton paths on every induced subset by literal directed-word DFS, not endpoint DP. Enumerate set partitions by restricted-growth block assignment, with at most three blocks, and multiply the literal subset counts. This independently computes the three unordered cover counts on all532 classes.
- Form their paired values and compare the entire deduplicated image with the source119-value set. Equality holds, not just agreement on the number of values.
- Recount every one of the119 displayed witness inputs by the same literal cover counter and recheck all primality/compositeness fields by exact trial division.
- Build all11 prime/seven-prime paired graphs, check their strongness and actual anti-involution on every arc, and count their Hamilton paths by inclusion-exclusion over unrestricted walks. Signed128-bit intermediates are used for this separate counter. All11 counts agree with the catalog.
- Check the physical transitive Phi values through input order12, including the exact cutoff at8.

The independent prime list is5,509,120293,189701,38997197,49745741,55619813,58696949. The seven-prime list is37050461,46937093,56228501, with prime quotients5292923,6705299,8032643.

## Reproduction

The standalone independent audit reads only the included `inputs.csv`, an exact extraction of the119 displayed source witnesses. It regenerates its complete class bank internally; no external tournament bank is needed.

    g++ -O3 -std=c++17 independent.cpp -o /tmp/s197-paired-catalog-audit
    /tmp/s197-paired-catalog-audit inputs.csv
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 independent.cpp -o /tmp/s197-paired-catalog-audit-debug
    /tmp/s197-paired-catalog-audit-debug inputs.csv

Both builds produce identical receipts. Hashes:

- independent.cpp: `4ebffc1bcb41fdbfc74459285bdc1f4d72c5b0a6acc33344b134a16412bcbd94`.
- inputs.csv: `6037849d9fedf4f1753232837bc04fa2ef22260d08379f5ad648ac6334c2c118`.
- independent-normal.json and independent-debug.json: `7c18347617d394fe9e5445468ebbfce09940f74ed2e2fc7d6339d0ea2acd631b`.

The source replay receipts remain alongside these independent controls. The all-order interval reduction rests on the separately reviewed gate theorem; finite equality of the two independently generated complete images certifies this catalog. No infinite prime family or unrestricted coverage is asserted.

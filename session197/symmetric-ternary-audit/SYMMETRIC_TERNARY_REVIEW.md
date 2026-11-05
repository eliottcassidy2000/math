# Independent audit: fixed-root self-converse ternary amplitude barrier

Status: **PASS, no correction required.** This review is specific to the frozen files listed below. No source files were modified.

## Frozen inputs

- `session197-prime-atoms/SYMMETRIC_TERNARY.md`: SHA256 `2df7d079d49bb6e2d8287a610b275d827510f861a5793908619cea4b521408ce`.
- `session197-prime-atoms/verify_symmetric_ternary.py`: `ade0e1e3461578ac4a6d47a8a7274c55fbd277e7144a04084c74e5bf0ac62cb9`.
- `session197-prime-atoms/receipt-symmetric-ternary.json`: `134c76d8ef4bb06f925624f60d00aa226cb6211c6c7c88d851ece0c3539b3689`.
- Named dependency `MINIMUM_TERNARY.md`: `45d183b1df1f030aff1bf755bd2b5b30ce441e8a5c13e01b99b6ef8613957c54`.
- Named dependency `verify_minimum_ternary.py`: `2ad41156b91eeb868170ee688cdd08f1ba2065d00948a812200ea990e1c9e27f`.

## Mathematical checks

1. An anti-involution fixing the cloned root fixes no other vertex. Its source and its invariant root strong component have odd order. Reversing the transitive clone chain extends this symmetry, and every even clone multiplicity produces an even-order self-converse tournament.
2. The even self-converse congruence gives `A=2 modulo 4`. Eventual identities become exact from the nonnegative finite run expansion and independence of distinct exponential functions. The maximum run equals one plus the root bipartite matching number, so a pure ternary identity requires matching number two.
3. Counts of unchanged strong components divide both `3A-1` and `9A-1`. Their gcd is one for even `A`; every such component is therefore a singleton. Deleting those components preserves every clone count and preserves the root-fixed anti-involution on the remaining strong component.
4. For `A<=74`, the root strong count is at most 221. Busch's exact order-11 lower bound is 225; all orders at least 12 are excluded by `5^11>221^3`. Odd source order therefore leaves only orders 1, 3, 5, 7 and 9. This uses the strong-component reduction before applying the order bound.
5. Reversal followed by the specified anti-involution acts on an odd set of Hamilton words and hence has a fixed word. Labeling along that word puts the fixed root at its middle and the anti-involution at index reflection. The normalized finite universe is consequently complete for the marked-root problem, not only for unrooted isomorphism classes.
6. The exact scope is maintained: the theorem is `A>=78` for a root fixed by an anti-involution; it asserts neither attainment at 78 nor a restriction on arbitrary roots of self-converse tournaments.

## Independent computational controls

The source verifier was run normally and under `python -O`. Both output files equal the frozen receipt byte for byte.

A separate C++17 checker, `check_symmetric_ternary.cpp`, regenerates the full normalized universe without importing source code. It uses a subset-state matching algorithm instead of augmenting paths, a forward Hamilton-path recurrence instead of terminal removal, and a physical cloning construction with clones labeled first instead of last. It checks every tournament and physical anti-involution. It reproduces all five rank histograms, totaling 66,067 normalized sources and 2,711 rank-two sources.

There are exactly nine first-equation survivors, three at each of orders 5, 7 and 9. The independent checker evaluates the third clone for all nine, even when the first count cannot yield integral amplitude. Every survivor has `(H_1,H_2,H_3)=(9,29,81)`, so none satisfies `H_3=9H_1+8=89`. The pure-clock survivor set is empty. The even-clone congruence also holds in all 2,711 rank-two controls.

Independent checker SHA256: `4f7bc290c3246c82b9c0b4582a698522be19dea7f9edfefa8ba5863f28421b02`.
Independent receipt `symmetric-ternary-independent.txt` SHA256: `1b234890466ad1015fbb894957ce0e41d31b76d9ba5d469fce3d7d5204d84326`.

Reproduce with:

    g++ -O2 -std=c++17 check_symmetric_ternary.cpp -o /tmp/check_symmetric_ternary
    /tmp/check_symmetric_ternary

The exhaustive order-nine result and the all-order conclusion are separate steps; completeness at arbitrary order relies on the exact SCC and strong-order reductions above.

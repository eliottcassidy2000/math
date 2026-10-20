# Verification scope, controls, and reproduction

## Frozen universes

- Classification: every odd non-prime-power-order subgroup of S_n up to S_n-conjugacy for 1≤n≤13, retaining fixed points and embedding. There are 48 degree/class pairs, including twenty at thirteen.
- Prime-thirteen transitive-image check: all 64 skew connection sets; all 720 bijections of each six-element out-neighbor set, with the complementary images forced and then checked on every ordered pair.
- Seven-plus-two-triples classification: every coordinate-surjective vector subspace of F3² or F3³, modulo precisely the ambient conjugacies available on the triple orbits. Actual proposed permutations are mapped into that exhaustive quotient list.
- Each census row: all 2^d H-fixed tournaments and graphs, compressed only by the proved normalizer action. The largest row has d=24 and 16,777,216 objects on **each** side. No limit or sampling excludes a row.
- Sign-even target: full automorphism group edge-reversal character, not degree parity and not merely the relative normalizer character.

## Independent routes

1. **New complete twenty-row census.** `audit_census.py` uses the previously audited standard-library isomorphism/sign engine, reproduced verbatim as `fixed_cube_engine.py`; the affine-block normalizer and bucketed classifier are copied in `extension_census.py`. `cube_orbits.cpp` exactly enumerates the signed compressed normalizer actions with orbit-stabilizer checks. Normalizer candidates comprise all local affine normalizers, equal-sized orbit swaps, and fixed-point permutations; conjugation of the defining generators filters to the actual normalizer. The explicit subgroup constructors alone are extracted from the incoming source, and the classification proof separately checks their completeness/nonduplication.

2. **Low-degree positive and hostile controls.** All ten rows through degree ten were compared at the full representative/relative-stabilizer/orbit-size level with the original Python coordinate-orbit implementation. They include C15's first deficit and the compensated diagonal-F21 escape. The new engine then replayed all 28 classified rows through twelve and matched every available aggregate count in the older audited references.

3. **Entire large cubes by independent generator BFS.** `generator_action_replay.py` constructs directed pair orbits independently from the supplied permutations. Instead of listing/deduplicating every normalizer element, it uses elementary generators of the local affine groups and fixed-point symmetric group. `generator_orbits.cpp` traverses every cube vertex and propagates the graph sign character along generator edges. A negative cycle detects a relatively odd stabilizer. Every accepted minimum representative and its orbit size agrees exactly with route 1, for both missing rows and both source/target cubes. The BFS visits all 16,777,216 or 4,194,304 vertices, not just the accepted vertices.

4. **Independent full isomorphism and signs for both large rows.** `independent_large_audit.py` uses separate vertex fingerprints and a constrained adjacency search inherited from `independent_engine.py`, not the census's simultaneous refinement and twin-quotient classifier. It reclassifies every accepted normalizer representative, obtaining the same tournament and graph multiplicity histograms. Graph signs are tested on a full point-stabilizer-chain generating set, using actual vertex permutations, not quotient automorphisms. Known twin transpositions are directly verified vertex permutations in this chain; no graph quotient is used. It tested 249,932 nonidentity graph automorphisms for the C15 row and 68,675 for F21. All accepted globally even types and all escaping types agree.

5. **Cycle-partition Burnside check.** `burnside_control.py` independently computes the source counts as 2 b_5(4)=17824 and b_6(2)=3040 using exact rational arithmetic. This verifies relative source counts only; it is not used as a substitute for full-type classification or global sign checks.

6. **Independent incoming implementation.** All eighteen incoming degree-thirteen rows, computed by the separate NetworkX/VF2 source, match every reported column of this twenty-row replay. The incoming receipts and prior audited packets also match on all thirty common rows. `compare_receipts.py` checks these facts using frozen exact receipts, not prose claims. The source constructors are shared as input data; the orbit, isomorphism and graph-sign implementations are different.

## Normal and optimized execution

The twenty degree-thirteen rows reproduce byte-identically in normal Python and `python3 -O`; the large-row result was additionally checked by separate full-vertex and generator-BFS routes in both modes. The prime-thirteen and subdirect-space certificates and receipt comparison also run under `-O`. All validation added here uses explicit exception checks, not Python assertions. Both C++ traversals use explicit checked conditions, not disabled assertions.

The normal/optimized receipts and checksums are in `reproduction_receipt.json`. Full type records are retained as gzip-compressed JSON (`large13.json.gz`, `small13.json.gz`); the 24-bit action tables, BFS orbit dumps, temporary executables and uncompressed full receipts are reproducible scratch outputs and are intentionally omitted from publication.

## Reproduction

Requires Python 3's standard library and a C++17 compiler. No external Python package, group database, network request or remote write is needed.

From this packet directory:

    python3 reproduce.py

This compiles the two C++ engines, reruns the low-degree controls and the whole degree-thirteen census, checks normal-versus-optimized identity, runs the two independent large-row checks, checks subgroup completeness certificates, compares incoming/reference receipts, and emits `reproduction_receipt.json`. The program also regenerates local scratch action/orbit files. Estimated execution time in the audit environment was under two minutes.

## Limits

The full-isomorphism and sign counts are finite exact computational results, not an all-k formula for C15 defects. The classification proof and finite certificates establish the group universe independently of those counts. Existing small-degree classification and all-p-group equality results are explicitly inherited. This packet supplies no degree-fourteen completeness theorem and does not prove tau_T≥tau_E for arbitrary odd H.

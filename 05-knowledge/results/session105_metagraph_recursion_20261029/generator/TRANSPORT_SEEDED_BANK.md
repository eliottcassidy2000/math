# Transport-seeded, orbit-first preparation of the full rooted bank

Session105. **PROVED:** the transport-torsor cache identity below. **FINITE-EXACT:** complete-output equivalence through child order seven. This is an engineering integration into the actual Session88 full-port generator, not a new scalar recurrence, not a novelty claim for group actions, and not a replacement of Session87 canonical augmentation.

## Inheritance and precise delta

The closest existing mechanisms are:

* Session87 `check_generation.py`: canonical augmentation with invariant root selection, parent mask-orbit traversal, literal-matrix cache, and all attaining transports. Its refined selector already improves ordinary-class-only generation; none of those improvements is new here.
* Session88 `blockfill/check_blockfill.py` and `check_adaptive.py`: compute one row per root-switch block, fill the entire block by aligned transports, and infer additional root records from certified subgroup orbits. These already implement the user's row-identifies-whole-block insight, retaining all vertex-indexed transporter sets. They are unchanged here.
* Session88 `seed_order/STABILIZER_FIRST.md`: choosing a smaller root orbit is not uniformly cheaper; it saves just 14 classifications at order eight and has an explicit counterexample. We do not change its seed/query ordering.

Inspection found Session88's common preparation still (a) canonicalizes the normalized parent a second time and (b) requests `rooted` classification for every raw binary mask. This artifact changes **only that preparation function**. It transports the already complete parent maps to derive the normalized parent's complete group and seed the exact literal cache; then imports Session87's orbit-first traversal into the full-data Session88 workflow. The actual adaptive row classification algorithm is untouched.

Source loader verifies SHA256 of both Session88 modules; their loader additionally verifies the Session87 baseline. Published paths are under `05-knowledge/results/session87_four_vertices_enumeration_20261028/` and `session88_symmetry_certificates_20261028/blockfill/`. The baseline is actual inherited code, not a weakened reimplementation. Recovered local source provenance is Session88's verified repository, not an assertion that the old checkout is current main. No source files in that checkout are modified.

## Transport identity (all orders)

Let `P` be a literal parent, `C` its fixed canonical matrix, and let `A` be the **complete** set of canonicalizing maps `a:C→P`. Fix `a0∈A`. Then

`Aut(C) = { a0^(-1) ∘ a : a∈A }`.

Each displayed composition is an automorphism of C. Conversely, for any `g∈Aut(C)`, `a0∘g` belongs to A, so every automorphism occurs. This is a bijection, not a generating set approximation. It returns exactly the maps that a second canonical search of C would enumerate. Since all maps preserve degree cells, sorting them returns the same attaining-order sequence as the inherited degree-cell permutation search. Thus `(canonical_code, sorted_group)` is a valid literal-cache entry for C without another search.

The mask state `(code,min orbit(mask))` is the complete aligned pointed key. Scanning masks in increasing order and marking each full `Aut(C)` orbit inserts exactly the same keys and literal children as Session88's raw-mask classification. Retaining **all** group maps is essential: one canonicalizing map would yield only the identity and falsely split symmetric mask orbits. The resulting `parent_data` and ordered child bank are unchanged, so every downstream full-port output remains unchanged.

This proof does not assume a tiny group, rigidity, prime decomposition, or a score/H observer. The cache key remains the entire literal adjacency matrix, never its score sequence, H-profile, or abstract group type. Space includes all group maps, masks, and full-port outputs; there is no sub-output-size or polynomial-time claim.

## Exact universe and verification

All ordinary tournament isomorphism classes are included (no converse merging, Hamiltonian/regular/prime restrictions). Child counts are `1,1,2,4,12,56,456`. Production canonicalization has the inherited active bound `order <= child_order-1`; full-child canonicalization is isolated audit only.

For each level the checker compares all six inherited output fields: pointed states, Q, every vertex-indexed rooted transporter set, every complete automorphism group, parent alignment data, and the ordered literal child bank. It independently canonicalizes every child, checks its automorphism group and root-orbit weights, and exhausts the full labelled universe through six (`32768` tournaments at six).

Additional exact arithmetic checks are:

* every mask orbit satisfies `orbit_size * stabilizer_size = |Aut(P)|` and all `2^n` masks are covered;
* for each child U, `sum_root_orbits 1/|Aut(U,v)| = N/|Aut(U)|`;
* globally `sum_parent_classes 2^(N-1)/|Aut(P)| = sum_child_classes N/|Aut(U)|`;
* multiplying the previous mass by `N!/N` gives `2^(N choose 2)`.

The inherited 81/27 hostile is re-executed: equal named subset-H data and local stabilizers still have unequal full transport actions and attachment orbit counts 176 versus 512. This optimizer retains that action, not merely the row weights. The inherited malformed/partial-row guards and misleading Q-square examples are re-executed as well. Guards are exceptions, not assertions, and the entire suite is also run with Python `-O`.

## Measured operation counts

| Child order | Classes | Baseline canon requests | New requests | Baseline searches | New searches |
|---:|---:|---:|---:|---:|---:|
|1|1|3|1|1|1|
|2|1|5|2|1|1|
|3|2|9|4|2|1|
|4|4|29|11|4|3|
|5|12|110|42|11|8|
|6|56|652|256|56|46|
|7|456|6289|2649|576|525|

The removed raw-mask requests mostly were **cache hits**, so their reduction is not a proportional runtime speedup. The transport-seeded entries also remove genuine searches, including later searches that would otherwise miss the literal cache. Adaptive root classifications remain 2593 at order seven: no claim that row classification itself improved.

Three alternating local repeats at seven gave median 0.11881 seconds for the inherited adaptive baseline versus 0.11311 seconds for the new preparation (about 1.05x). At six the medians were 0.00968 versus 0.00918 seconds. These tiny timings are illustrative, hardware/run dependent, and not a performance record. Timed regions include full bank/preparation, subgroup closure, Q, all transports, and groups; equality audits occur after timing. The exact operation reduction is the stronger result. This benchmark compares identical full-port outputs; it does not compare against Session87's faster, smaller ordinary-class-only deliverable.

## Loss ledger and remaining boundary

Source: full literal parent plus all attaining maps. Target: identical full-port state bank with an exact seeded cache. Preserved: matrices, canonical parent codes, binary attachment orbits, actual common-frame actions, named ports, row multiplicities, full child groups, and all transporters. Destroyed information: none in the requested outputs. Sidecar: complete attaining-map torsor, retained rather than compressed to an abstract group. Cheapest decisive test: equality of all six fields, not just A000568 counts or Q rank.

This optimization does not turn row support/weights into a constructor for missing named ports. A full ordinary metagraph compression still owes its transport sidecar. Nor does this produce a closed recurrence from earlier scalar A000568 values: representative/group data remain indispensable to this generation implementation.

## Reproduce

```sh
python3 -B check_preparation.py --baseline-root /path/to/math --max-order 7 --repeats 3 > results.json
python3 -B -O check_preparation.py --baseline-root /path/to/math --max-order 7 --repeats 3 > results.optimized.json
```

The script prefers the packet’s self-contained `vendor/` root, then discovers an enclosing repository. There is no sister-checkout fallback. An explicit `--baseline-root` overrides discovery and remains hash-verified. JSON timings differ naturally; structural results and six output hashes must agree. No publication is performed by this artifact.

## Self-contained pinned dependencies

The packet vendors only four unmodified Python sources under their original repository-relative paths in `vendor/`. `vendor/SOURCE_MANIFEST.json` records byte counts and SHA256 values. These are Session87’s `check_generation.py` and Session88’s `baseline_sources.py`, `check_blockfill.py`, and `check_adaptive.py`; all use only the standard library. The first file contains the entire 81/27 construction, exact listed actions, mask-orbit calculation, and independent named-subset-H dynamic program, so no fixture data or additional repository files are required. The baseline loader records Session87 source pin `7a0bab9bfa7bb87c5a10c7ba239592d90cfafd51`; the inherited generator’s mathematical source pin is `b8672af00bcc6d440a9a00a6915325faf39fcdff`. Session88 source bytes are pinned by the two explicit SHA256 checks in this packet and the manifest; no unverified current-main provenance is claimed.

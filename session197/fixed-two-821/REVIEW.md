# Independent audit: the sharp four-cycle threshold at H821

**PASS, no correction.** The all-length fixed-port closure excludes H821 with at most three Hamilton cycles. The displayed word has H821, C4 and feedback number exactly two. The infinite family and its k=0 exception are correct. A separate complete strong-graph census through order13 independently confirms the negative conclusion without using the port transfer.

## 1. Frozen source binding

The reviewed packet is `session197/fixed-two-821/`:

- `PROOF.md`: `0377f24695e9cb9f3b17c3e4b62679841d686cb787baee6c9a3e51e533d0ddcf`.
- `support.cpp`: `b12772a34df283312ffdc624c8acc45500c42d2b378c63354853d000ba11f705`.
- `SUPPORT.json`: `eef8f0ad01612e90441ac3e5a5c746447f5ae155fb1516fc14e2640c3950e06d`.
- `verify.py`: `91e90cc4fe84d68bd60d3ca06783cc58e38095cee159a0eab31c55e1e45e8ffb`.
- Normal/optimized validation: `95ba99fb35e91f3bc0f64c47325e8d578cd2b7219d6ca4aa3c5feedbb2225874`.

The eight source files were copied into this audit's `source/` directory before replay. Both Python modes compile and run both optimized and unoptimized C++ closures. The resulting receipts agree with the frozen bytes. The source files were not changed.

## 2. All-order closure proof review

Deleting the two permanently retained ports leaves the transitive birth order. Conversely any feedback deletion set of size at most two can be enlarged to two initial vertices at these relevant orders; the remaining graph is transitive and its exact port incidences determine a valid word. Thus this is the full stated feedback-two graph class, rather than only a convenient subfamily.

The inherited 27-coordinate transfer and its cycle-edge readouts are nonnegative. Coordinate saturation at822 therefore commutes with every transition and readout. A clipped coordinate cannot cancel from a later nonnegative sum. In particular H<=821 and C<=3 are both exact retained predicates. The implementation sums the cycle coordinates on the incoming state and records that count on the birth edge, which avoids losing the observable when successor states coincide.

The H pruning rule uses a valid injection: every old Hamilton word admits insertion of the named new tournament vertex at a gap or an endpoint, and deletion recovers the old word. Therefore a rejected large-H prefix can never return below the cap.

The C pruning rule correctly compares a rejected strong prefix with the eventual strong target. A prefix with a Hamilton cycle is strong. Strong induced cycle monotonicity follows by a mixed singleton insertion or, when all outside incidences are uniform, the forced two-vertex ear. Both constructions inject old cycles. No cycle-monotonicity claim is made through arbitrary weak intermediate graphs.

All four masks are tried with newborn drop index2, so the same two physical ports remain retained. State merging preserves every future capped response. Queue exhaustion, not a maximum word length, supplies completeness. The emergency state ceiling would produce a non-success return and is not reached. The source replay reproduces all374,692 processed states,441,723 admissible edges,213,747 cycle rejections and missing odds7,21,821. The observed representative length23 is not used as a cutoff.

## 3. Independent complete negative check by a different representation

The independent route here does not use the cover automaton. It uses the previously audited strong-singleton-extension completeness theorem in `C3_COMPLETE.md`: every strong tournament with C<=3 and at least four vertices grows from the unique strong four-vertex tournament by strong singleton extensions. The exact ear argument excludes missing two-vertex steps. That proof holds at every order, not just through the earlier census boundary.

`probe_c3_order13.cpp` is the frozen independently reviewed census source `verify_c3_complete.cpp` (original SHA `e49f0b9ff6e29f7b7260d6516ccb0fafe03ea4bbcf06dccc0664b7fddead806f`) with only maxn changed from12 to13. It enumerates every mixed extension, applies the proven cycle-insertion lower bound, counts survivors exactly, and canonicalizes by every actual Hamilton cycle and rotation. The complete bank through13 has80,746 isomorphism classes. The newly included order13 row is C1=6765,C2=15360,C3=30676. Optimized and unoptimized assertion-enabled runs reproduce every CSV byte and inventory entry.

The exact previously checked Busch strong-tournament consequence is H^3>=5^(n-1). Here

    821^3=553387661 < 1220703125=5^13,

so a strong H821 graph has order at most13. This bound is used only for this independent audit route; the source theorem uses no Busch bound. It is not a bound on weak singleton padding.

`probe_lowcycle_821.cpp` separately counts Hamilton paths in all80,746 literal graphs by endpoint subset DP, capped at822. It finds exactly15 H821 cores:

- Two order11,C3 cores, both with feedback number3.
- Nine order12,C2 cores with feedback number3.
- Four order12,C2 cores with feedback number4.

There are no C1 hits and no order13 hits. The optimized and assertion-enabled unoptimized filters agree exactly. Orders at most three cannot have H821.

The Python checker independently recounts H on every one of those15 cores by inclusion-exclusion over unrestricted walks, counts C by anchored directed-cycle DFS, and tests every induced vertex subset for transitivity using its degree sequence. All feedback numbers exceed two. Thus this complete bounded audit excludes every strong H821/C<=3 feedback-two graph, independently of the source's capped cover closure. It is an order-bounded proof with an explicit valid cutoff, not an inference from census stabilization.

The relevant inherited completeness proof is included in the already published H533 dependency bundle at commit231bda1d4b2f57506bfcd6a2b2b3decb5f5f089d, `session197/lowcycle/feedback533/source/session197-lowcycle/C3_COMPLETE.md`.

## 4. Independent positive and boundary controls

The checker independently constructs the six-vertex seed from masks0031 and obtains the exact retained-root profile (h,e,R,s,d)=(37,5,49,4,5). Appending a universal sink preserves h and s, makes e zero, and doubles R by placing the sink at the end of either cover slot. After at least one such sink, the final mask1 therefore gives H=37+49*2^k and C4. This is a physical family; it is not a returning map on each already-closed output.

For k=1,...,6 the checker builds the actual graphs and verifies H by walk inclusion-exclusion and C by literal anchored DFS. At k=4 its rows equal the source's eleven-vertex H821 witness. Every vertex subset is tested, giving exact feedback number two. This graph has eleven vertices, hence55 arcs; its whole adjacency arrays are compared literally to the specified fixed-port word. The k=0 graph is counted independently as H81, correctly differing from the invalid extrapolation86.

## 5. Reproduction artifacts

- Census source: `probe_c3_order13.cpp`, SHA `69b492857503f904d224e2c67f4f526f5e5091bc7cb83c6eb0125e33970005ad`.
- Separate path filter: `probe_lowcycle_821.cpp`, SHA `0af8133883b1c1638ae56d8055b12f054c324deaa55a481a986e3b2c3a7b2d2a`.
- Complete regenerated CSV SHA: `71aafd0283d5fe4f93406ce65c3f546ded04ed2dacb951da9532e92259579779`.
- Inventory: `c3-order13-inventory.jsonl`, SHA `67b5eded66a0e3c92bd6be3e6e048a98f47a54e9262668b8397b694db7305f2f`.
- Fifteen hit records: `lowcycle-821-order13.csv`, SHA `859e44ce74b3b371e33cf6d296f282c54e60236ee543a0738ceb8e5c6866bccf`.
- Independent checker: `check.py`, SHA `36014b923d828c463c5c2819fea4cd2273b1daefa6abd18ce0b726f296cbc0e8`.
- Normal/optimized independent receipt: `ad18e704e3961230db721d391861630c78cc30eaae99f4f8ce9dbe851d71100b`.

The large CSV is reproducible rather than a required input. Compile the census with C++17 and run it with a fresh CSV output filename. Compile the filter and pass that CSV as its sole argument; stdout must match the supplied fifteen-hit CSV. Run `python3 check.py` or `python3 -O check.py` for the independent walk, feedback and positive-family checks. Source closure replays use `python3 source/verify.py`. All code is standard-library-only.

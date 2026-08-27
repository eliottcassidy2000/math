# `(520,663)` five-mask signature-surgery promotion packet

Status: **FINITE-EXACT, PROVED RELATIVE TO THM-4281, WITH A DETACHED
LITERAL COMPLETENESS AUDIT.**  This packet does not prove a physical entry or
LRC(14).

## Exact result

The frozen post-THM-4281 full-signature census gives `(520,663)` inactive
indices

```text
57,107,222,275,345
```

in the 421-mask joint deck.  The corresponding masks are

```text
08c0a980,128c8900,009041ac,1aa28002,18868880.
```

Deleting them leaves 416 masks and exposes exactly 53 labelled nine-bodies.
An independent scan of all `C(30,9)=14,307,150` bodies reproduces their typed
`(body, deleted-response)` FNV `1e5d6dfe0b676151`.

At `(520,663)`, exhaustive enumeration of all `C(30,8)=5,852,925` rank-eight
masks gives 2,879,147 active masks, FNV `6dbc8f16ef5c3ff7`.  The 53 obligations
have 10,784,970 raw complement incidences and 3,490,314 distinct candidate
masks.  Of these candidates, 1,253,628 are active.  Their complete labelled
response quotient has:

```text
7,124 classes including the empty class;
7,123 nonempty classes;
810 inclusion-maximal classes;
2,720,981 active response incidences;
class FNV             2158656972d58de7;
active-response FNV   a8eb70069447e610.
```

The exact minimum replacement cover is six.  A least-representative witness
is

```text
0090492c,018c2114,09a2a040,108c1112,1c81a100,38120016.
```

The lower bound is independently checkable from the integer dual printed in
`signature_surgery_atlas.out`: the 53 obligation weights have numerator sum
57 on denominator 11, and every one of all 7,124 response classes has weight
at most `11/11`.  Therefore every cover has size at least
`ceil(57/11)=6`.  The six masks above cover all 53 obligations, so the bound
is sharp.

Appending those six masks after the retained 416 gives a 422-mask deck with
FNV `813801c9bd1676ba`.  A second full labelled body scan has zero failures.
All 422 masks are active at `(520,663)` and at the neighbouring fibre point
`(520,589)`.

## Independent completeness gate

The minimum-six claim does not rely only on the primary endpoint-cocycle
enumerator.  `literal_response_universe_audit.cpp` reconstructs the two actual
safe combs on the literal grid `18241159416480`, integrates their intersection
over all 7,133 fixed-pool cells, uses a separate single-array superset
enumerator, enumerates every mask disjoint from any of the 53 obligations, and
rebuilds the entire response quotient.

Its O3, O0, and UBSan transcripts are byte-identical; all sanitizer stderr is
empty.  Its generated 7,124-row class TSV and its full 2,879,147-row active
response TSV are byte-identical to the primary outputs.  The large active TSV
is intentionally excluded from this compact packet; `REPRODUCTION.md`
generates it transiently and verifies SHA256
`dbdd1dab8c1e0d4de03866d20155d2a985929bd9e279b3fb15eb643c351a851e`.

## Global consequence

On the 172,322-row post-THM-4277 residual, the rebuilt 422-mask deck is common
on exactly 145,122 rows:

```text
FNV       0d25ac854255d5b2
SHA256    f197ade21ea861725d941e9e1d6a5148d3d233805671b20ab4dfc358610ac189
tests     64,564,282
equalities 0
```

The endpoint-cocycle primary and detached literal scouts emit byte-identical
common-row CSVs and the same weakest row/mask `(35,315):22560801`; their raw
ratios cross-multiply exactly.

Relative to the canonical THM-4281 common set, the rebuilt deck gains 562
rows.  Two, `(520,667)` and `(567,664)`, were already removed by the canonical
carrier certificate, so the exact post-THM-4281 gain is 560:

```text
FNV       faf264670166cbc0
SHA256    68cf6c593fece339bd9f86d8d467febc28bfee496cb73a93e89bc40b77ad7611
top       (520,663)
```

The rebuilt deck loses 3,503 rows of the canonical 148,063-row common set.
It is therefore a **complementary certificate**, not a replacement for the
canonical THM-4281 deck.

The 560-row gain is disjoint from the separate 26-row index-367 singleton
fibre certificate.  Conditional on that sibling packet, their combined
586-row proof-graph deletion leaves 23,637 rows, FNV
`e8b363d2b3d9ba6a`, SHA256
`21ad2530da6adc7187b4b5829c9fddad07ab704ccab96adb90d546e43756e95d`,
with unique top row `(256,663)`.

## Dormant-carrier cross-filter

As a reconstruction sidecar, the 8,951-mask dormant carrier input intersects
the detached active-response TSV in 7,220 masks.  Exactly 3,007 have nonempty
response, their union covers all 53 obligations, and they contribute 6,580
incidences.  Per-obligation hit counts range from 18 (obligation 38) to 290.
The filtered ledgers are frozen in `dormant_carrier_active_responses.tsv` and
`dormant_carrier_obligation_hits.tsv`.

This cross-filter makes **no** claim that the dormant carrier is a body cover,
that its provenance is promoted, or that it supplies a minimum cover.  It is
only a THM-4281/reconstruction audit connection.

## Proof architecture

1. `lost_body_probe.cpp` independently exhausts all labelled nine-bodies.
2. `signature_surgery_atlas.cpp` reconstructs the signature, primary active
   universe, complete 53-bit response quotient, exact cover search, rational
   dual, rebuilt deck, and full body rescan.
3. `literal_response_universe_audit.cpp` supplies the detached literal
   completeness audit of the entire active response universe.
4. `verify_surgery_certificate.py` replays all class multiplicities, least
   representatives, witness coverage, every dual inequality, deck surgery,
   and primary/literal common-set identity.
5. `rebuilt_global_primary.cpp` and `rebuilt_global_literal.cpp` provide
   independent global scans; `rebuilt_pair_activity.cpp` checks the target and
   nearest signature-subset point.
6. `combined_consequence.py` performs exact proof-graph deduplication and the
   conditional index-367 union.
7. `cross_filter_dormant_carrier.py` performs the explicitly scoped dormant
   carrier reconstruction audit.

## External code dependencies

The packet sources use the promoted repository implementations below.  Their
SHA256 identities at freeze are:

```text
f9c6dca98169a579a67ad173f02c91af38edf01dbd81152e9540cc8e38157eaa  response_pattern_atlas.cpp
d2388ad0901ea23d9ba2360c1dd4fd8ec2d1055b8cbd3d5a86125add5b9eff14  carrier_cegar_descent.cpp
c4acbbc18eb0d9b3bb3105efe090ae4bb08ccc15ac8c230b4af14dec0db00627  cascade_pair_exhaustive_primary.cpp
25a6978484c7ab122fdc6c8e1593cfa2ad3468f7184a156045ea7e6cb2efc45d  dependency_lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp
```

## Scope and exclusions

- Every finite universe, threshold convention, ordering, FNV serialization,
  and proof-graph subtraction is explicit in the sources and transcripts.
- “Active” uses the exact nonnegative threshold.  The global rebuilt-deck scan
  has zero equalities, so every accepted deck test there is strict.
- The five deleted masks are deck indices; the 53 rows are labelled body
  obligations; response classes are quotients of rank-eight masks.  These
  types are not interchangeable.
- The index-367 fibre uses a different one-deletion/one-replacement deck and
  is included here only in the conditional disjoint-union consequence.
- No endpoint below the scanned universe, physical-entry consequence, or
  LRC(14) theorem is claimed.
- Binaries, compiler progress stderr, duplicate literal CSVs, and the 96 MB
  active-response TSV are excluded.  They are reproducible transient outputs.

`SHA256SUMS` is the authoritative raw-LF artifact ledger.

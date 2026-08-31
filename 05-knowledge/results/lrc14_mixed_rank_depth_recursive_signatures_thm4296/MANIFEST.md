# THM-4296 artifact manifest

**Status:** `PROVED RELATIVE TO THM-4287 + THM-4295 + FINITE-EXACT + DETACHED
LITERAL-WALL AUDITS PASS`. This packet proves fixed-pool sufficient
certificates only. It does not prove physical entry or LRC(14).

## Typed universe

- Pool: THM-4254's fixed labelled thirty-speed pool.
- Threshold: `1/14`, with target safe mass `alpha=4/63`.
- Body universe: all `C(30,9)=14,307,150` labelled nine-bodies.
- Inherited proof-graph universe: THM-4287's `22,647` rows, FNV
  `df5374d4aca67677`, SHA-256
  `14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317`.
- Rank-eight repair universe: `C(30,8)=5,852,925` masks.
- Rank-nine repair universe: `C(30,9)=14,307,150` masks. Every mixed-rank
  program retains labelled masks and bodies; no symmetry quotient is used.

The exact consequence lemma is rank-free. If an active deletion mask `R` is
disjoint from a nine-body `B`, then `B` is contained in `P\R`, so the safe set
for `(P\R) union {q,r}` is contained in the target safe set for
`B union {q,r}`. Activity is an upset under enlarging `R`.

## Frozen inputs

| path | role | semantic identity |
|---|---|---|
| `inputs/joint421_masks.txt` | ordered joint deck `E` | 421 rank-eight masks, FNV `20d63dd42fe8150e` |
| `inputs/reconstructed_final8951.txt` | inherited carrier base | 8,951 masks, FNV `188f82ab9dd1695a` |
| `inputs/additions45.txt` | inherited additions | 45 masks, FNV `ec083b65cc8c34e3` |
| `inputs/endpoint638_response_witness9.txt` | inherited endpoint suffix | nine masks, FNV `02b936529030e4bc` |
| `inputs/current_residual22682.csv` | pre-THM-4287 residual used by the detached exchange boundary audit | 22,682 typed rows |
| `inputs/current_residual22647.csv` | typed THM-4287 residual | 22,647 rows, FNV `df5374d4aca67677` |
| `inputs/full_signatures_primary.csv` | complete inherited signature atlas | all current source rows and seven signature words |
| `inputs/endpoint636_failures101.csv` | THM-4287 carrier boundary | 101 labelled pair-body obligations |
| `inputs/post_exchange_r632_failures72.csv` | first exchange-carrier boundary | six plus 66 labelled obligations |
| `inputs/r629_failures28.csv` | 9,013-carrier boundary | 28 labelled obligations at `(100,629)` |
| `inputs/r628_failures4.csv` | 9,017-carrier boundary | four labelled obligations at `(100,628)` |

The aggregate proof graph also reads four frozen row ledgers from the
independent THM-4295 packet. They remain governed by that packet's manifest
(`ac6f86ff077d425e1089992a4de2ce8322ce8162429f6841c5eed76623c93f00`):

| cross-packet path | role | exact identity |
|---|---|---|
| `lrc14_endpoint636_exchange_recursive_ideals_thm4295/inputs/carrier_layers636_633.csv` | append-only carrier consequences | 10 rows, FNV `9926701692e6f8d4`, SHA-256 `00f37ae823268f8b2c79f2c5ce2634d487ef9c5f7c2864f5b172dd1f89fdefd6` |
| `lrc14_endpoint636_exchange_recursive_ideals_thm4295/inputs/signature19_fibre36.csv` | independent index-19 consequences | 36 rows, FNV `5c8af37cf2f002e7`, SHA-256 `c0406e8d10138ad743e3e99076b2171809945518cc4f708d3ef46a7b0fb70777` |
| `lrc14_endpoint636_exchange_recursive_ideals_thm4295/inputs/signature294_fibre21.csv` | index-294 consequences | 21 rows, FNV `eadefa2fae582ca7`, SHA-256 `07f6727333f3a677e61b20afd44b0a60d2621515020bce68cfe610904785f79b` |
| `lrc14_endpoint636_exchange_recursive_ideals_thm4295/inputs/signature372_fibre54.csv` | index-372 consequences | 54 rows, FNV `47ab2af18f07ff59`, SHA-256 `1400300c4d0d43f935bbdb7c686f6053b70568bdb719a8cc61d4baaf49875eb0` |

`SHA256SUMS` is the authoritative raw-byte identity for every input, result,
control, and documentation file in this packet.

## Mixed-rank carrier conclusions

### Endpoint 636 exchange

The complete 101-obligation rank-eight response quotient has exact addition
minimum fourteen in the active union of the two endpoint rows. Fourteen old
masks are strictly inactive on all nine rows at endpoints 637 and 636, so the
minimum witness can replace them without increasing carrier size:

```text
exchange carrier: 9,006 masks, FNV 8062ce6d5728da1f.
```

The detached audit checks every deleted sign, every appended response, and
all `128,764,350` labelled body instances. It reports zero equalities and zero
failures.

### Exact deletion depth at `(256,632)`

The body `1d106401` has no active disjoint rank-eight mask among all
`C(21,8)=203,490` candidates. Its unique best rank-eight mask is still
strictly inactive. All thirteen available one-bit extensions of that mask are
strictly active at rank nine. Therefore its deletion depth is exactly nine.

The complete local rank-eight response atlases are frozen as
`results/endpoint/r632_rank8_atlas/R632_LOCAL_100.tsv` and
`R632_LOCAL_256.tsv`; `results/endpoint/r632_rank8_exact_cover.out` proves
their separate minima three and ten. The hostile complement audit supplies
an independent packing lower bound and explicit ten-mask cover for the 65
nonhostile `(256,632)` bodies.

The six `(100,632)` failures have exact rank-eight cover minimum three; the 66
`(256,632)` failures have exact mixed-rank cover minimum three. Appending both
witnesses gives a 9,012-mask carrier, FNV `be9cf87002d6114a`.

### Exact continuation

| repaired boundary | exact additions | rank composition after repair | carrier FNV | next failure boundary |
|---|---:|---:|---|---|
| `636` exchange | 14 added / 14 deleted | `9006+0` | `8062ce6d5728da1f` | `632` |
| `632` | 6, as two pairwise minima `3+3` | `9009+3` | `be9cf87002d6114a` | `630` |
| `630` | 1 | `9010+3` | `932d1fae6fd8c384` | `629` |
| `629` | 4 | `9011+6` | `07689a1534ce7327` | `628` |
| `628` | 2 | `9013+6` | `d7f0e06e154e78c2` | `626` |

At `r=629`, a detached implementation reproduces the 419 response classes,
the exact `7/2` rational dual lower bound, the four-mask witness, and the
9,017-mask carrier identity. At `r=628`, the implementation scans all 32,767
nonempty subsets of the fifteen possible nonzero response patterns, rejects
subsets using any of the four unrealized patterns, and thereby exhausts every
subset of the eleven realized classes. The exact minimum is two.

The final scan completes the entire inherited `r>=626` prefix: 72 rows. It
closes 70 and fails only at

```text
(100,626): 995 bodies, FNV f0189c7a4ab34b77;
(256,626):  10 bodies, FNV f394b2fdaca1f3f8.
```

The O0 and O3 transcripts have SHA-256
`390bc2844b9fa24250ef2e4e19be02cda2702563c2673f8fdbb32f010ca50d70`;
the byte-identical failure ledgers have SHA-256
`60ee709e9acebeb3b8e09d88a283cb8389fdb48b8f102188f386d9365a6429d1`.

## Recursive signature conclusions

The complete singleton census has 3,738 rows in 192 groups. Exactly 110
groups admit one-deletion/one-addition common-deck surgery. The groups are
disjoint and close 1,219 rows, FNV `3706723fd3574334`, SHA-256
`fdef599768aba3cc1738e330f35c9b431486ab845023ef59050c546f22605ef1`.
These are 110 separate 421-mask decks.

The index-19 ideal has 36 rows, FNV `5c8af37cf2f002e7`. Its common-activity
response quotient has 41 classes and exact minimum two, witnessed by
`0000e649,0184a205`. The rebuilt 422-mask deck has FNV
`dc50478119bc6c12` and zero failures on every labelled body.

The primary discovery reports contain noncanonical cross-row weakest-margin
annotations because their raw numerators have pair-dependent grids. Per
MISTAKE-532, only `results/singletons/detached_literal_audit.out` and
`results/j19/detached_literal_audit.out` are authoritative for weakest-margin
fractions. This does not affect any sign, cover, deck, or proof-graph result.

## Typed proof graph

| node | certificate type | rows | FNV |
|---|---|---:|---|
| singleton success | 110 separate 421-mask common decks | 1,219 | `3706723fd3574334` |
| index 19 | one 422-mask common deck | 36 | `5c8af37cf2f002e7` |
| endpoint descent | one mixed-rank 9,019-mask carrier | 70 | `213e8d6fcbd0f1cd` |
| THM-4295 index 294 | one separate 422-mask common deck | 21 | `eadefa2fae582ca7` |
| THM-4295 index 372 | one separate 421-mask common deck | 54 | `47ab2af18f07ff59` |

Among this theorem's three local nodes, the only overlaps are

```text
singleton / endpoint: (410,626), (506,626)
index-19 / endpoint:  (338,628), (338,636)
singleton / index-19: empty
triple:               empty.
```

This theorem's three-node typed union has 1,321 rows, FNV
`69bcc5c5fc86ac8e`, SHA-256
`ed040ba5b26d40445e8917cc24c844d36c3a91f59bd08244f6a3bf0f3d38a44e`.
Its exact complement has 21,326 rows, FNV `a05a372bf57ea064`.

For the cross-theorem join, THM-4295's ten append-only carrier rows are a
proper subset of the 70 endpoint-descent rows, although the carrier masks are
different, and its index-19 row set equals the index-19 node above. Its
index-294 node overlaps the singleton union in 18 rows and the endpoint node
at `(294,636)`; its only new rows are `(147,294),(147,590)`. Its index-372 node
overlaps the singleton union in 52 rows and the endpoint node at `(372,636)`;
its only new row is `(372,619)`. The two larger ideals are disjoint, and all
triple intersections among the five nonredundant nodes are empty.

The aggregate typed union therefore has 1,324 rows, FNV
`f55ee025df29bb65`, SHA-256
`57bdcd932cc2c985e81e1b2d472a469cf4c2e11c9cb21dd4b1037c4ba098562a`.
The exact complement has 21,323 rows, FNV `09a0dfc4515d556b`, SHA-256
`a6fb3ea00f5017b2ad66cbccf75f89756aa94360652ff6a7dbf22c637a1c7656`,
maximum endpoint 626, top `(100,626),(256,626)`.

The independent C++ proof-graph audit derives the index-19 ideal from the
signature atlas, parses the raw final scan, verifies its row prefix against
the inherited universe, and requires byte equality with the three maintained
local ledgers. A second Python/C++ pair reads both artifact packets, audits all
21 pairwise intersections and exact membership profiles, and requires byte
equality for the 1,324-row union, 21,323-row complement, and three-row increment.

## Audit matrix

| claim | primary exact path | independent/control path |
|---|---|---|
| endpoint-636 minimum and exchange | complete response atlas and exact cover | integer dual-gap checker; detached literal exchange/body audit |
| depth-nine obstruction | exhaustive hostile complement scan | full rank-eight universe scan plus all thirteen rank-nine extensions |
| r632/r630/r629/r628 minima | complete response quotients and exact set cover | literal witness replays; detached r629 dual/carrier consumers |
| final 9,019 descent | literal fixed-pool scanner | O0/O3 byte identity and independent raw-scan parser |
| 110 singleton ideals | exhaustive 192-group census | detached literal signs and private-body scan |
| index-19 ideal | full common-activity quotient | detached literal signs and exhaustive body scan |
| typed union and complement | Python exact set consumer | independent C++ derivation/parser/partition audit |
| cross-theorem typed join | Python exact overlap/profile consumer | independent C++ set algebra and byte-equality audit |

No result here supplies physical entry, exclusive owner, semantic arrival,
termination, an arbitrary-pair theorem, or LRC(14).

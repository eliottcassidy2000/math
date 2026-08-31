# THM-4295 artifact manifest

**Status:** `PROVED RELATIVE TO THM-4287 + FINITE-EXACT + DETACHED
LITERAL-WALL AUDITS PASS`. This packet proves fixed-pool certificate
statements only. It does not prove physical entry or LRC(14).

## Typed universe

- Pool: the labelled thirty-speed pool frozen by THM-4283 and inherited
  through THM-4287.
- Threshold: `1/14`, with activity decided by exact signed integer margins.
- Input residual: THM-4287's `22,647` ordered rows, FNV
  `df5374d4aca67677`, SHA-256
  `14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317`.
- Ordered joint deck: the inherited 421 rank-eight masks, FNV
  `20d63dd42fe8150e`.
- Inherited carrier: the 9,006 masks of THM-4287, FNV
  `fdc1c57ae4dc1bb6`.
- Body universe: all `C(30,9)=14,307,150` labelled nine-bodies.
- Repair universe: all `C(30,8)=5,852,925` labelled rank-eight masks.

Carriers, rebuilt common decks, inactive-signature classifiers, and their
typed proof-graph union are distinct consumers throughout this packet.

## Frozen local inputs

| path | role | exact identity |
|---|---|---|
| `inputs/append14_masks.txt` | downstream-optimal 14-mask carrier suffix | `14`, FNV `cb6e6d8963cf54d0`, SHA-256 `30c6b29c223f2785d878df469c1f31df7b9cbbf7b8819c4f5329d1602367c03e` |
| `inputs/base9006_failures636_632.csv` | all inherited carrier failures on the complete endpoint-636 and endpoint-632 layers | `176` data rows, SHA-256 `17c7f6c4635cb092cb91effb88cdec92bdb5844708c41bfaa50d54314ba648cd` |
| `inputs/carrier_layers636_633.csv` | ten rows closed by the 9,020-mask carrier | `10`, FNV `9926701692e6f8d4`, SHA-256 `00f37ae823268f8b2c79f2c5ce2634d487ef9c5f7c2864f5b172dd1f89fdefd6` |
| `inputs/signature19_fibre36.csv` | complete ideal `I_E(p) subset {19}` | `36`, FNV `5c8af37cf2f002e7`, SHA-256 `c0406e8d10138ad743e3e99076b2171809945518cc4f708d3ef46a7b0fb70777` |
| `inputs/signature294_fibre21.csv` | complete ideal `I_E(p) subset {107,173,247}` | `21`, FNV `eadefa2fae582ca7`, SHA-256 `07f6727333f3a677e61b20afd44b0a60d2621515020bce68cfe610904785f79b` |
| `inputs/signature372_fibre54.csv` | complete ideal `I_E(p) subset {20,50,399,401}` | `54`, FNV `47ab2af18f07ff59`, SHA-256 `1400300c4d0d43f935bbdb7c686f6053b70568bdb719a8cc61d4baaf49875eb0` |

The three witness files beside the fibre ledgers freeze the exact replacement
masks printed below. The inherited joint deck, full signature atlas, carrier
prefixes, and live residual remain pinned in the THM-4287 artifact packet.

## Exact carrier conclusions

The 101 endpoint-636 failures inherited from THM-4287 have a row-aware
response quotient with `43,660` responsive masks, `1,835` nonempty response
classes, and `259` inclusion-maximal classes. A deterministic exhaustive
search over the maximal quotient proves cardinality minimum `14`; an
independent HiGHS integer program returns objective `14` with zero gap.

The selected downstream-optimal suffix is

```text
12041017 18468880 0024e442 19c04044 21a0d100 00807407 1041a209
10c086c0 10704441 08c48c40 08432a80 11681044 02809125 016c1044
```

It gives a 9,020-mask carrier with FNV `920651ae987fcdef`. A detached literal
audit applies activity at each obligation's own row, finds zero activity
equalities, reconstructs the inherited `64+37` failures, and verifies zero
failures after the suffix is appended.

The carrier closes every live row at endpoints `636`, `635`, `634`, and
`633`, respectively `6+2+1+1=10` rows. Its first failed complete layer is
endpoint `632`, where it leaves exactly `1` failure on `(100,632)` and `36`
on `(256,632)`. The failure ledger has SHA-256
`c36879862b4168404c6f8fcad7f205b16a9fd8cbc8e79b68aaacdea8d7598528`.

For the larger four-row-response problem consisting of the inherited
`64,37,8,67` failures on `(100,636)`, `(256,636)`, `(100,632)`, and
`(256,632)`, the full quotient has `57,055` responsive masks, `2,620`
nonempty classes, and `379` maximal classes. Subject to covering all 101
endpoint-636 obligations with at most 14 masks, the exact minimum number of
endpoint-632 misses is `37`. CP-SAT and an independently written HiGHS model
both prove objective `37`; their optimal witnesses need not coincide.

## Absolute endpoint-632 response obstruction

One remaining body is

```text
pair=(256,632), body=1d106401.
```

All `203,490` rank-eight masks disjoint from this body are inactive at this
row. The source-pinned audit finds active count `0`, equality count `0`, and
strict maximum signed margin `-40875925173725376`, attained at `222c100c`.
The detached literal-wall implementation independently obtains reduced
maximum margin

```text
-2937663943/371177847040
```

and margin-ledger FNV `6aa0de375fd3b976`. O0, O3, and UBSan outputs agree
byte-for-byte. Consequently no rank-eight append-only carrier, of any added
cardinality, can close this fixed-pool row. This is a stopping obstruction
for that sufficient-certificate operation, not a physical-danger claim.

## Three recursive signature ideals

For the ordered joint deck `E`, each ideal below is defined inside the full
THM-4287 residual by `H_D={p:I_E(p) subset D}`. Activity is intersected over
the entire ideal before a replacement mask is admitted.

| ideal | deleted indices | rows | exposed bodies | common-active masks | response classes / maximal | exact replacement minimum and witness | rebuilt deck |
|---|---|---:|---:|---:|---:|---|---|
| `H19` | `{19}` | `36` | `8` | `2,212,775` | `41 / n.a.` | `2`: `0000e649,0184a205` | size `422`, FNV `dc50478119bc6c12` |
| `H294` | `{107,173,247}` | `21` | `34` | `2,503,052` | `1,986 / 128` | `4`: `12848902,0380a520,0380e012,1288b200` | size `422`, FNV `1de1e54f8262ac49` |
| `H372` | `{20,50,399,401}` | `54` | `24` | `1,678,262` | `1,342 / 126` | `4`: `19801608,11842450,01180c85,11803211` | size `421`, FNV `f63671c303bbafe4` |

Each primary path exhausts the response quotient and every labelled body.
Each import-independent literal-wall audit reconstructs all activity cells,
finds zero equalities, and again finds zero body-cover failures. O0, O3, and
UBSan controls are byte-identical for every detached audit.

The larger top signatures at `(100,636)`, `(256,636)`, and `(384,636)` were
censused during discovery but their common-activity intersections were not
completed or promoted into this packet. No common-deck claim is made for
them.

## Typed proof graph

Let `K` be the ten-row carrier node and retain `H19`, `H294`, and `H372` as
three separate common-deck nodes. The ideals are pairwise disjoint. The only
cross-node overlaps are

```text
K intersect H19  = {(338,636)}
K intersect H294 = {(294,636)}
K intersect H372 = {(372,636)}.
```

Their exact typed union and complement are

```text
union:    count=118, FNV=4cc16b6090fa8461,
          SHA256=4893b22ddfaafb296094e60104de477d2291560ec69c8bc466a7109453a5254d;
residual: count=22,529, FNV=64fe87a7c95ff64f,
          SHA256=dc6d42c3f8b47251600b390737bdaca3565114506d69ead69ad4d5656dd05a13,
          maximum endpoint=632.
```

The ten top rows are `(100,632)`, `(256,632)`, `(282,632)`, `(294,632)`,
`(306,632)`, `(320,632)`, `(332,632)`, `(366,632)`, `(416,632)`, and
`(512,632)`. Python and independently written Ruby consumers emit
byte-identical union, residual, and transcript files.

## Audit matrix

| claim | primary path | independent/control path |
|---|---|---|
| minimum-14 endpoint-636 response cover | exhaustive 101-bit response quotient and depth-first exact cover | SciPy/HiGHS integer set cover, zero gap |
| downstream optimum 37 at cardinality 14 | OR-Tools CP-SAT on the 176-bit quotient | independently written SciPy/HiGHS model, zero gap; direct selected-witness scan |
| selected carrier closure at endpoint 636 | source-pinned joint-exposure scan | detached literal row-aware audit; O0/O3/UBSan controls |
| complete live-layer descent through the first failure at endpoint 632 | source-pinned joint-exposure scan | direct selected-witness replay on the frozen residual and failure ledger |
| zero-response body at `(256,632)` | source-pinned complete activity universe and exact signed margins | import-independent literal-wall enumeration; O0/O3/UBSan controls |
| three signature-ideal surgeries | complete signature ideals, common-activity intersections, response quotients, and body scans | three import-independent literal-wall/body scans; O0/O3/UBSan controls |
| 118-row typed union and 22,529-row complement | Python exact set consumer | independently written Ruby exact set consumer |

No symmetry quotient, sampling, floating-point wall decision, or unlabelled
body identification is used. The two MIP packages are used only on frozen
binary response atlases; all underlying activity and incidence data are
exact integers.

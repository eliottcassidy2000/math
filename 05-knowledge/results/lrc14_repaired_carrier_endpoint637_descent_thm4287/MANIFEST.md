# THM-4287 artifact manifest

**Status:** `FINITE-EXACT / DETACHED LITERAL-WALL AUDITS PASS`, relative to
THM-4283.  This packet proves fixed-pool certificate statements only.  It
does not prove physical entry or LRC(14).

## Typed universe

- Pool: the labelled thirty-speed pool frozen by THM-4283.
- Threshold: `1/14`, with activity tested by exact integer wall mass.
- Input consumer: THM-4283's exact residual `U_4283`, with `22,682` rows,
  FNV `f7563445f15efebf`, and SHA-256
  `7d0044bc4c32f08b9d09dca420171a05666d26e03f38fbc48a9baa1fcb027102`.
- Body universe: all `C(30,9)=14,307,150` labelled nine-bodies.
- Repair universe: all `C(30,8)=5,852,925` labelled rank-eight masks.

## Frozen inputs

| path | role | exact identity |
|---|---|---|
| `inputs/joint421_masks.txt` | ordered THM-4281 joint deck | `421`, FNV `20d63dd42fe8150e`, SHA-256 `fcde04de68ab614743176ed153e0db2754cf878470c96f525d29c07a88388bc2` |
| `inputs/reconstructed_final8951.txt` | inherited carrier prefix | `8,951`, FNV `188f82ab9dd1695a`, SHA-256 `a5dac3c7e5a2715e4c9ef8bb1b54bc98792e904f7b0b5ef55e4dd4313ebc87f6` |
| `inputs/additions45.txt` | THM-4282 carrier additions | `45`, FNV `ec083b65cc8c34e3`, SHA-256 `63cbd655c884a83c88b17e2f1f1bd511f8d7ba0aceeaeb4f6371e8ba54351616` |
| `inputs/endpoint638_response_witness9.txt` | THM-4283 boundary additions | `9`, FNV `02b936529030e4bc`, SHA-256 `7b0198b219924d4622df38e31b9f75c749cce7588bce633eaaddb54b34993a2b` |
| `inputs/full_signatures_primary.csv` | full inherited inactive-signature atlas | `24,223` data rows, content FNV `3652b5590b330704`, SHA-256 `4f0e8da3fdab8bd5a0e14f3b4fa30050602025f63486aa35e0cf03374e6e3832` |
| `inputs/current_residual22682.csv` | typed THM-4283 residual | `22,682`, FNV `f7563445f15efebf`, SHA-256 as above |
| `inputs/carrier_endpoint637.csv` | exact maximal endpoint layer | `3`, FNV `f4b48b16a28826ad`, SHA-256 `e31697893a8a4dbeb4b22ea62d7ef2bfed7fd55bbe35138a1bfc6f460638601b` |
| `inputs/signature_fibre33.csv` | exact current signature ideal `I_E(p) subset {275,345}` | `33`, FNV `f11425c1e1b17094`, SHA-256 `5b5de4e502802287ff2c504d34cfc4de21c10adffd1694c4bfa1c9a1d735c85b` |

The 8,996-mask carrier is `8,951+45`, FNV `fd899660f14b311c`.
Appending `014c9084` gives the 8,997-mask carrier, FNV
`8e1860a25d0fcf87`; appending the nine-mask witness gives the 9,006-mask
carrier, FNV `fdc1c57ae4dc1bb6`.

## Exact conclusions

### Carrier boundary

The 8,996-mask prefix, the 8,997-mask intermediate carrier, and the repaired
9,006-mask carrier each close all three current endpoint-637 rows:

```text
(100,637), (294,637), (520,637).
```

The 9,006-mask carrier does not close the complete endpoint-636 layer.  It
has exactly `64` failures at `(100,636)`, FNV `d611500ea833ff83`, and `37`
at `(256,636)`, FNV `ee7792a8a2fd51c9`; the other four rows close.  The
complete two-layer failure ledger has `101` rows and SHA-256
`e5da728e0a40076141e43450983dd4bb0f6d9c5aa281f95becb5dbafc555e242`.
The detached literal audit independently reconstructs all activity and all
`9*C(30,9)=128,764,350` body instances, with zero activity equalities.

### Signature-fibre surgery

The target `(520,637)` has inactive signature `{275,345}`.  Deleting joint
masks `1aa28002` and `18868880` exposes exactly six bodies.  Intersecting
activity across the complete 33-row signature ideal leaves `2,027,392`
common-active masks and 17 response classes.  There is no full responder;
the exact cover minimum is two, attained by

```text
02229124, 1c030206.
```

The rebuilt 421-mask deck has FNV `9781de64dff60cee`, zero failures on all
labelled nine-bodies, and strict activity in all `33*421=13,893` cells.  The
source-independent literal audit agrees and has zero equalities.

### Typed proof graph

The carrier layer has three rows, the signature-fibre deck has 33, and their
only overlap is `(520,637)`.  Their typed union therefore removes 35 rows:

```text
count=35,
FNV=0b5d8d9c28d7325d,
SHA256=1d36b9e6899496f93d44785d389523d4f3a85efe4eaadaf2e1297904116bbd78.
```

The exact new residual is

```text
count=22,647,
FNV=df5374d4aca67677,
SHA256=14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317,
maximum endpoint=636,
top={(100,636),(256,636),(294,636),(338,636),(372,636),(384,636)}.
```

Python and independently written Ruby consumers emit byte-identical output
and byte-identical residual ledgers.

## Audit matrix

| claim | primary path | independent/control path |
|---|---|---|
| endpoint-637 closure and endpoint-636 failure boundary | inherited joint-exposure implementation, exhaustive bodies | detached literal walls, O0/O3/UBSan, byte-identical failure ledger |
| 33-row signature ideal and exact minimum-two response | full signature atlas, complete common-activity intersection and response quotient | detached literal activity/body replay, O0/O3/UBSan |
| 35-row typed union and 22,647-row complement | Python exact set consumer | independent Ruby exact set consumer |

Every body-cover assertion is about a labelled deck.  A carrier, a rebuilt
common deck, an inactive-signature classifier, and their typed set union are
not interchangeable certificate types.

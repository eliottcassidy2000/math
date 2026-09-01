# Independent endpoint-588 direct-literal audit

Status: **FINITE-EXACT**. This packet independently audits a fixed-pool
certificate layer. It proves neither a physical-entry statement nor LRC(14).

## Inherited object and universe

- pool labels: `8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,120,126,132,143,145,168,170,176,190,193,240,252,264,286,290`;
- inherited THM-4318 carrier: 3,925 distinct masks, rank split
  `3809+116`, FNV `eeae5518d84ccac5`;
- protected joint deck: 421 rank-eight masks, FNV `20d63dd42fe8150e`;
- typed endpoint-588 layer: 66 rows, FNV `18cf9a572cf9a5be`;
- labelled body universe at each row: all `C(30,9)=14,307,150`
  rank-nine masks.

The carrier serializer imports the canonical THM-4318 reconstruction only to
freeze the inherited object. The endpoint-588 carrier and wall audits are
self-contained clean-room implementations and import no maintained audit code.

## Exact unchanged-carrier replay

The clean-room replay tested 944,271,900 row-body pairs. It found:

| quantity | exact value |
|---|---:|
| exposed bodies | 5,001,257 |
| nonjoint hit incidences | 76,551,030 |
| carrier failures | 144,867 |
| rows with a failure | 10 |
| failure-ledger FNV | `6f51fa88f3b09cdc` |
| pair-ledger FNV | `204a72794170e186` |

The ten nonempty failure fibres have sizes

`25:9, 50:99836, 96:1130, 100:14, 105:60, 206:7, 210:43799, 256:7, 420:1, 462:4`.

Both the 67-line pair ledger and the 144,868-line failure ledger agree
line-for-line with the independently produced scout ledgers. O2 and O3 builds
also agree byte-for-byte inside this packet.

## Direct literal exit

For every pair-safe wall cell let `F` be the set of failed pool labels and,
for a failed body `B`, define

`L(B) = sum(width : F intersection B is empty and |F| <= 9)`.

Every cell counted by `L(B)` is literally safe for all nine labels in `B`, so
`L(B)` is a lower bound on the full literal safe mass. The independent wall
reconstruction proves

`63 L(B) - 4 G > 0`

for all 144,867 failed bodies, on each row's exact integer wall grid `G`.
The minimum truncated surplus over the complete failure ledger is

`706240375488648`

ticks, at row `(210,588)` and body `06b82090`. The complete literal mass is
also positive for all bodies. Truncated and full masses are equal for 109,452
bodies. The two compiler builds produce identical 144,868-line detail ledgers,
which also agree line-for-line with the scout's independently generated primary
ledger.

Consequently the unchanged THM-4318 carrier plus direct literal wall mass
closes all 66 endpoint-588 typed rows in this finite fixed-pool certificate
layer. No carrier response or exchange is needed for this layer.

## Typed successor

An independent set-theoretic consumer removes the complete 66-row endpoint-588
layer and obtains:

| object | count | FNV | SHA-256 |
|---|---:|---|---|
| typed union | 2,207 | `18d067b5614cf47f` | `f03c84f15d9a149b0a083b50e922118e814d5644f5fa21e7011ae1c414ff3675` |
| residual | 20,440 | `794bd808e92e27cd` | `be2d63e98beefb062e9ae4436d490ee2f630352989bf309cf85b5a1ffc44278c` |
| next endpoint-587 layer | 10 | `f48ca5f1904d6f52` | `e2e841bccef0773513cc71d6b60ed03aa227cc701e19bc8a4673b4b7971d2a63` |

Normal and optimized Python modes agree byte-for-byte and match the scout's
typed successor artifacts.

## Scope boundary

This is an exact finite certificate for one inherited fixed-pool/typed-frontier
layer. It does not identify continuous physical runner speeds, does not cross a
quotient without its labelled-body sidecar, and does not prove LRC(14).

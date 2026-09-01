# Independent endpoint-587 direct-literal closure audit

Status: **FINITE-EXACT / PROVED RELATIVE TO THE SERIALIZED THM-4318
CARRIER. LRC(14) REMAINS OPEN.**

## Inherited objects and derived frontier

This packet consumes the independently serialized THM-4318 carrier rather
than importing the endpoint-587 scout or maintained carrier-audit source.  The
carrier has 3,925 distinct masks, rank split `3809+116`, FNV64
`eeae5518d84ccac5`; all 421 protected rank-eight masks (FNV64
`20d63dd42fe8150e`) occur in it.

The set-theoretic typed reconstruction starts from the endpoint-588 successor:
typed union 2,207 and residual 20,440.  Its maximal residual endpoint is 587,
with exactly ten rows and FNV64 `f48ca5f1904d6f52`:

`24,25,50,72,96,100,105,192,210,260` paired with 587.

## Exact unchanged-carrier replay

The import-free C++ replay reconstructs activity from raw, unaggregated
pair-safe wall cells.  It scans all
`10*C(30,9)=143,071,500` labelled row/body cases.  O2 and O3 builds agree
byte-for-byte on their pair and failure ledgers and give:

| quantity | exact value |
|---|---:|
| exposed bodies | 53,771 |
| nonjoint hit incidences | 1,587,855 |
| carrier failures | 41 |
| rows with failures | 1 |
| global failure FNV64 | `bae65dc3d3bd34d0` |
| pair-ledger FNV64 | `0062cc50be726e54` |

All failures lie at `(50,587)`.  The 41-body fibre has FNV64
`c76719ced1d5c52b`.  A recordwise comparison with the independently written
scout found identical pair and failure rows; its only raw-byte difference was
the platform newline convention.

## Direct literal exit

For every pair-safe wall cell, let `F` be the set of failed pool labels and,
for a missed rank-nine body `B`, put

`L(B) = sum(width : F intersection B is empty and |F| <= 9)`.

Every cell counted by `L(B)` is literally safe for all labels in `B`, so
`L(B)` is a rigorous lower bound on the full literal safe mass `M(B)`.  A
separate Python implementation enumerates raw wall cells using integer
arithmetic and proves, for every one of the 41 bodies,

`63*L(B) - 4*G > 0`.

For `(50,587)`, the exact grid is 53,537,802,887,368,800.  There are 8,385
open wall intervals, 5,792 pair-safe cells, 2,420 full failure classes, and
2,365 low-rank classes.  Their class FNV64 values are respectively
`b24c8782fd3e3c7c` and `4145ef7ae527e0e0`.

The minimum truncated (and full) scaled surplus is

`283,424,219,270,897,292`

at body `11186405`.  Truncated and full masses agree on 31 bodies; full mass
strictly dominates on ten.  Normal and optimized Python runs agree
byte-for-byte, and the 41-row detail ledger has SHA-256
`5a89b23b60fce1ec4be6d03bdd92c4a3cf3a861afe3d30226f8bfe831154c895`,
the same digest as both scout implementations.

Thus every endpoint-587 row closes without changing the THM-4318 carrier.

## Typed successor

Consuming the complete endpoint-587 layer gives:

| object | count | FNV64 | SHA-256 |
|---|---:|---|---|
| typed union | 2,217 | `e6592cbef9b616d8` | `d465a2f62c77ddaf921e7f3d7f32c674ea45e46fcdc348c90c711d7ba8a7a6b6` |
| residual | 20,430 | `4710f750dfcf91ea` | `8203dcfd6cc26b67bfdb648c7d8b50f94d7af1ab7ddbd6af2ff68acb15941f0b` |
| next endpoint-586 layer | 12 | `a1b617faa2e7f63f` | `d38b7fd9ea2aa9afdd12446d646cc8a9466cc5d4429612f03c8dff3165edf7ea` |

## Scope boundary

This is a finite exact certificate for a frozen fixed-pool carrier and typed
frontier.  It proves no continuous physical owner/entry map, no terminating
descent, and no instance of LRC(14) beyond this certificate layer.


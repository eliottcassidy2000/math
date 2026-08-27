# Detached hostile audit of the THM-4281 carrier at `(256,663)`

Status: **FINITE-EXACT, scratch-only, PASS** at repository commit
`6ff6ee322b40e806d8d512baa9fd1df5730ce8b9`.

This packet independently reconstructs the canonical ordered 8,951-mask
THM-4281 carrier and the 71 exposed labelled nine-bodies from canonical
repository artifacts.  It does not consume either corresponding supplied file
from the dormant active-response packet.  A detached direct-wall literal-grid
engine then evaluates every carrier mask at `(256,663)` and scans all
`C(30,9) = 14,307,150` labelled bodies.

## Exact conclusion

The final carrier closes `(256,663)`.  Its literal-active subdeck contains
5,427 masks and covers every labelled nine-body; the complete scan has zero
failures.

```text
ordered final carrier       8,951  FNV 188f82ab9dd1695a
carrier file SHA-256               a5dac3c7e5a2715e4c9ef8bb1b54bc98792e904f7b0b5ef55e4dd4313ebc87f6
active carrier              5,427  FNV fe7da30cec531ca5
inactive carrier            3,524  FNV 47830f63011b6cca
literal equalities              0
active parts          5,008 / 414 / 5   (base / joint / continuation)
inactive parts        3,516 /   7 / 1
body universe          14,307,150  FNV 0d0aed8567285749
full scan failures              0
full scan checks       499,416,624
maximum checks/body          5,315  at body 1d106401
```

The literal grid has denominator `291858550663680`, 8,969 cells, and
pair-safe mass `214426971029040`.  The weakest active mask is `27048240`,
with signed normalized margin
`18878014764 / 291858550663680 = 388149 / 6000874880`.  The closest inactive
mask is `0121d204`, with deficit
`6616222578 / 291858550663680 = 272071 / 12001749760`.

The certified new proof-graph edge is the singleton ASCII row `256,663\n`:
typed little-endian `(q,r)` FNV `aacf2f9de5beb267`, ASCII SHA-256
`28d9759ed8c3d7561f734a6552b0feb3ffee81cb675fab32245ac14f0fba600c`.

## Independent carrier reconstruction

The source appends masks in the canonical THM-4281 order, deduplicating at
every stage.  All intermediate cardinalities and typed FNV ledgers are hard
gates:

| stage | count | FNV |
|---|---:|---|
| 59-file THM-4254 replay band | 4,675 | `ce4e76ec11df057c` |
| plus `(416,704)` additions | 4,733 | `a7b046289655c733` |
| plus round-two `(416,700)` and `(520,700)` additions | 7,986 | `baef1d2f49444638` |
| plus `(384,694)` additions | 8,319 | `e08b227730f6793c` |
| plus THM-4271 `(520,688)` additions | 8,518 | `1603e3fe970f8428` |
| plus six THM-4276 masks | 8,524 | `5ddb84a44f5d2ad7` |
| plus the ordered THM-4281 joint deck | 8,945 | `3212efa05dd18c00` |
| plus six THM-4281 continuation masks | 8,951 | `188f82ab9dd1695a` |

The 59 replay-band file hashes are frozen in
`controls/replay_band_files.sha256`; its own SHA-256 is
`393fbd4d9d26a5609f319b1a3c4d1f9c77074d04d270c0b2bc57d3bc8857ffb4`.
The remaining direct canonical inputs and the complete literal-engine include
chain are frozen in `controls/canonical_dependencies.sha256`.

## Independent 71-body derivation

The engine partitions all 421 ordered joint masks directly at `(256,663)`:
414 are active and exactly seven are inactive, at zero-based indices
`75,107,139,374,394,405,417`.  Their ordered `(index,mask)` FNV is
`38812f899b3636a8`.

It then enumerates all 14,307,150 bodies and selects exactly those missed by
the 414 active joint masks.  This yields:

```text
vulnerable bodies            71
body FNV               414d30a143d2e22d
typed (body,pattern) FNV 287281de2900cca7
deleted incidences           83
CSV SHA-256            b0ff6815c928f0805cad6f9cdecdb34270e9e8908fcd53fc1c0ca5c7dbc8605a
```

As a local cross-check only, the reconstructed carrier is byte-identical to
the dormant packet's supplied carrier.  The reconstructed two-column body and
deleted-pattern CSV is byte-identical to the first two columns of the dormant
71-row input; all 71 values in its unused third column are the literal string
`0`.  The 5,427 active carrier masks hit those 71 bodies 4,575 times in total,
with no failure, a minimum of seven hits at `2d306400`, a maximum of 228 hits,
and response FNV `454f40855bd2aac7`.

## Hostile controls and scope

- O0, O3+NDEBUG, and O1+UBSan builds of the portable frozen source produce
  byte-identical stdout, carrier files, and vulnerable-body files.  All three
  stderr files are empty.
- Every consequence-bearing value is guarded by `require`, so `NDEBUG` does
  not remove the checks.
- The literal engine is the canonical detached direct-wall implementation,
  not the dormant packet's endpoint-cocycle/superset-zeta implementation.
- The 71-body response check is followed by a separate complete scan of every
  labelled body, so the closure does not depend on trusting the 71-row input or
  assuming that it was exhaustive.

The dormant packet's 2,042,937-mask full active universe and its 11,169-class
response quotient were not re-enumerated here and are not theorem
dependencies.  Their 18 MB/87 MB derivative lists are intentionally excluded.
This audit establishes only the inherited-carrier closure of the single pair
`(256,663)`; it makes no claim about another endpoint-663 pair, carrier
minimality, or discovery optimality.  It inherits the promoted transcripts,
the 421-mask deck, and the six already-certified THM-4276/4281 continuation
masks rather than rediscovering them.


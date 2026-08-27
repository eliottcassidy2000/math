---
id: THM-4266
title: "Three-round learned-carrier endpoint descent"
status: >
  PROVED RELATIVE TO THM-4254/4256/4261/4262 + FINITE-EXACT + INDEPENDENTLY
  AUDITED WITH LITERAL-WALL SEED/BOUNDARY CONTROLS. The three-round
  8,319-mask carrier closes all 3,037 current post-THM-4261+4262 residual
  edges with second endpoint at least 688 except (520,688). The resulting
  residual has 177,585 edges and unique top edge (520,688). The raw
  post-THM-4256 carrier closure has 3,336 edges, of which 299 were already
  removed by THM-4261/4262. This does not prove the retained edge, an
  arbitrary pair, physical entry, or LRC(14).
source: root/cross-frontier-bridge/2026-08-27
depends_on:
  - THM-4254-fixed-ceiling-band-signed-endpoint-cocycle-cascade
  - THM-4256-uniform-two-three-outsider-ray-endpoint-cocycle-closure
  - THM-4261-semantic-endpoint-band-prefix-union-lift
  - THM-4262-uniform-three-four-outsider-ray-common-deck-closure
artifact_root: 05-knowledge/results/lrc14_three_round_learned_carrier_thm4266
artifact_manifest: 05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/SHA256SUMS
artifact_manifest_sha256: 0aa27866747df1182629da76263931957dc02c6365d3ca3f21362fbe108bd77b
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The frozen 39-file canonical subset verifies in a fresh
  staging tree. All C++ entry points compile; independent Python residual,
  proof-graph and layer replays byte-match. O0/O3/ASan+UBSan literal-wall
  controls agree, sanitizer stderr is empty, and fresh hostile replays recover
  the unique layer-688 resistance and the layer-732 boundary exactly.
---

# THM-4266 -- three-round learned-carrier endpoint descent

**PROVED RELATIVE TO THM-4254/4256/4261/4262 + FINITE-EXACT + INDEPENDENTLY
AUDITED WITH LITERAL-WALL SEED/BOUNDARY CONTROLS; LRC(14) REMAINS OPEN.**

## 1. Exact universe and statement

Retain the fixed labelled pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290}
```

and, for a finite positive set `A`, put

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14},   alpha=4/63.
```

Let `E_current` be the semantic residual after proved THM-4256, THM-4261, and
THM-4262.  THM-4261 and THM-4262 remove disjoint sets of 297 and 72 edges
from the post-THM-4256 residual.  The frozen current census is

```text
COUNT 180622
FNV-1a-64 0cef4e2887c8f24e
SHA256(lines q,r\n) fa1c5672b0f2cd2490413e9b69a4720bf1dc4eef8aee694c1c73d390aba58e11.
```

The theorem is:

```text
For every (q,r) in E_current with r>=688 and (q,r)!=(520,688),
and every B in binom(P,9),

    mu(G_(B union {q,r})) >= 4/63.                     (1)
```

The current high-endpoint universe has 3,038 edges.  Exactly one is retained,
so `(1)` closes 3,037 proof-graph-new edges.

## 2. Carrier semantics

An eight-mask `R` encodes an eight-subset of `P`.  For a fixed pair `(q,r)`,
call it active when

```text
mu(G_((P\R) union {q,r})) >= 4/63.                    (2)
```

A finite carrier `C` closes `(q,r)` when every `B in binom(P,9)` is disjoint
from at least one pair-active `R in C`.  Indeed, then

```text
B subset P\R,
G_((P\R) union {q,r}) subset G_(B union {q,r}),       (3)
```

and `(2)` implies `(1)` in the required direction.

The carrier is not a common-active deck.  Each row re-evaluates every mask's
mass, and only that row's active subdeck enters its body scan.  Adding masks
is nevertheless monotone: the active subdeck for an already closed row can
only grow, so a later union cannot destroy an earlier closure.

Every carrier FNV in this document is FNV-1a-64 over the ordered masks, each
serialized as an unsigned little-endian `u64`.  Every edge FNV serializes `q`
then `r` as unsigned little-endian `u64` words in inherited residual order.
An edge SHA-256 is over the corresponding raw ASCII lines `q,r\n`.

## 3. Three CEGAR rounds

The inherited THM-4254 carrier is the first-occurrence union of 59 labelled
order-minimal prefixes:

```text
pair-labelled incidences 56419   FNV acc8347addf27ac3
unique masks              4675   FNV ce4e76ec11df057c. (4)
```

For a resistant row, the discovery program thresholds the full universe of
`binom(30,8)=5,852,925` repairs, orders the active repairs by the frozen
SplitMix order from THM-4254, and takes the first prefix that covers all
`binom(30,9)=14,307,150` bodies.  The preceding prefix misses a recorded
body, so minimality is only within this order; no globally minimum deck is
claimed.  Only masks novel relative to the current carrier are appended.

The three exact rounds are:

| round | resistant seed(s) | full prefix size / FNV | masks appended | resulting carrier / FNV |
|---:|---|---|---:|---|
| 1 | `(416,704)` | `2608 / 18ff663a123e684e` | 58 | `4733 / a7b046289655c733` |
| 2 | `(416,700)` | `5894 / 701c0233c8f8abeb` | 2,970 before overlap accounting | |
| 2 | `(520,700)` | `4557 / d8466558bcd8ef9d` | 1,194 before overlap accounting | `7986 / baef1d2f49444638` |
| 3 | `(384,694)` | `5307 / 7b9a46c60e8514f8` | 333 | `8319 / e08b227730f6793c` |

In round 1 the seed prefix overlaps the base carrier in 2,550 masks and adds
58 (`FNV f25cf1f803d09f8a`).  The enriched carrier closes every current layer
from 704 through 701.  Its first resistant layer is exactly

```text
r=700: (416,700), (520,700).                           (5)
```

For round 2, relative to the 4,733-mask carrier, the two novel lists have
FNV ledgers `bc358ae2b2fdc66d` and `fbfe11a97800f02f`.  They overlap in 911
masks (`FNV a7fc01a504ee1f4c` in `(416,700)` prefix order), leaving 2,059
masks unique to the first list and 283 unique to the second.  Their union
therefore appends 3,253 masks.  The frozen order first appends every
`(416,700)`-novel mask in that prefix's order, then appends only still-unseen
`(520,700)`-novel masks in its prefix order.  Round 3 likewise appends every
still-unseen `(384,694)` mask in that prefix's order.  This convention, not
set union alone, determines the carrier FNV ledgers.  The round-2 carrier
closes both seeds, the older compression hostile `(542,732)`, and every layer
699 through 695.  Its first resistant layer is exactly

```text
r=694: (384,694).                                      (6)
```

In round 3 the 5,307-mask seed prefix overlaps the round-2 carrier in 4,974
masks and appends 333 (`FNV 1d1422eaa770226d` in source-prefix order).  The
resulting 8,319-mask carrier closes `(384,694)` and every layer 693 through
689.  The first resistant layer is exactly

```text
r=688: (520,688).                                      (7)
```

The three independent atom/component ledgers for newly appended masks are

```text
round 1  65a40a194fdaf0e3
round 2  b212b9a61918f6c9
round 3  5866333b05284cad.
```

## 4. Exact transfer boundary

Round 1 closes its seed because all 58 appended masks are active there; two
are disjoint from the seed's unique failed body.  One such witness is

```text
0x004a0c0b={8,10,16,80,84,132,145,176}.                (8)
```

Transfer is pair-specific, not formal.  At `(416,700)` only 14 of the 58
new masks are active, and at `(520,700)` only 18 are active.  In each row the
sole new mask disjoint from the recorded first failed body is

```text
0x12258042={10,40,120,126,143,170,240,286},             (9)
```

but it is inactive there, with a negative exact activation margin.  Thus the
round-1 enrichment transfers through layers 703, 702, and 701 but cannot
close either row in `(5)`.  Rounds 2 and 3 work for the same concrete reason:
their appended union contributes newly active masks disjoint from every
previously failed body.  No stability of a mask's activity across pairs is
assumed.

The dedicated literal transfer control reproduces the 58-mask activity and
disjointness census at the seed and both layer-700 hostiles.  Its O0, O3, and
ASan+UBSan reports are byte-identical with SHA-256
`e46884f71700a371cd834b2c273c891ada4f53bab805fd23f3a1ac81a3bfdfc1`;
sanitizer stderr is empty.

At the final boundary `(520,688)`, the 8,319-mask carrier has 5,934 active
masks and leaves exactly two bodies uncovered.  The first is

```text
0x07187008={16,85,88,95,145,168,193,240,252}.           (10)
```

The best carrier repair disjoint from `(10)` is

```text
0x10458a01={8,63,84,120,126,143,176,286},               (11)
```

but its activation margin is still

```text
-71,286,709,268,960,640 / 11,420,425,087,469,798,400.  (12)
```

This is the smallest stopping reason in the audited descending range.  It is
a failure of the fixed 8,319-mask compression only, not a counterexample to
the target inequality or to the endpoint method.

## 5. Independent controls

The primary implementation reconstructs the THM-4254 4,675-mask carrier,
checks every full source prefix for order-minimal closure, recomputes exact
signed endpoint atoms, verifies newly appended masks by an independent
component formula, and exhaustively scans every current residual row in each
reported layer.

The clean-room implementation does not use endpoint atoms to decide mass.
For each controlled pair it builds the literal joint wall grid for the 32
speeds, evaluates safety at exact cell midpoints, integrates every carrier
repair directly, and recursively enumerates all 14,307,150 bodies.  It checks
all four full source prefixes and all four CEGAR seed rows.  The dedicated
round-1 seed control handles `(416,704)`; its O0/O3 reports are byte-identical
with SHA-256

```text
c8b9b4bad67940b4ceb1033aef326e5c34c095704f5ef407ff1bae2464cd5ec0.
```

The round-2/3 control handles `(416,700)`, `(520,700)`, and `(384,694)`, then
checks the next resistance `(520,688)`.  Its O0, O3, and ASan+UBSan reports
are byte-identical:

```text
SHA256 5d8a5a4145de1d0e1cf68ac0ed6ea84efd6d03b6e1817728248ad1be7bee7abf
sanitizer stderr SHA256 e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855.
```

The full-discovery O0/O3 normalized transcript hashes are

```text
(416,704) 680c4c3de04ab9310c94150a4d1235a3cdcd732d11ccf53186998caa38385f7c
(416,700) 464b635f14420aa10f0b56e52008f6a1b7797a3257920eef7814adb32b7e172e
(520,700) ed11c6389a8a57f11ab86261d8a81261efcbc471e792bb1eaa3b2586ea54ebdb
(384,694) 04024bf5c03aff2a8dc20e6a9204e59d1f703ab3c474a7eb1f7721a3ea6d0536.
```

The high-band proof combines the following exact row partition:

| endpoint block | rows | exact result and carrier basis |
|---|---:|---|
| 754..725 | 498 | base carrier closes 497; only `(542,732)` resists and round 2 closes it |
| 724..705 | 1,088 | direct 8,319-carrier scan, all closed |
| 704..701 | 349 | round-1 carrier, all closed |
| 700 | 91 | round 1 leaves exactly two seeds; round 2 closes both |
| 699..695 | 505 | round-2 carrier, all closed |
| 694 | 107 | round 2 leaves exactly one seed; round 3 closes it |
| 693..689 | 575 | round-3 carrier, all closed |
| 688 | 124 | round-3 carrier closes 123 and retains only `(520,688)` |

These blocks total the raw post-THM-4256 3,337-row high-endpoint universe.
Of the 3,336 closed rows, THM-4261 already removes 297 and THM-4262 already
removes two; the remaining 3,037 are the current theorem contribution.  The
new 724--705 scan has SHA-256
`3084057742fddd32131e4084bbae4c3445a6ad9f4d116746d5db7bcbf640d5cc`.
The base scan over 724--718 has zero resistance and SHA-256
`24f27d9e493fe9db650b7b02a1b4e2845b0f5ca2013b03bd3058cef0fe3155bf`.
Its continuation over 717--704 gives zero resistance through 705 and sole
resistance `(416,704)`, with SHA-256
`70f0ab7336cf23c9d31f08a018eaeccfd3a4d197844ce6d699d3e09aa4c98b23`.
Carrier monotonicity justifies reusing earlier closures after enrichment; no
row is inferred from a neighboring endpoint.

## 6. Current proof-graph consequence and overlaps

THM-4256 removes exactly 73 post-THM-4254 ray edges:

```text
COUNT 73
FNV 6be23222a3a20764
SHA256 9d43ad3311533711e21c496141b05052c45ed68b8b03bcbce59737df3c7391ea
MAX ENDPOINT 690.
```

The raw 8,319-carrier component is the 3,336 post-THM-4256 edges with
`r>=688` other than `(520,688)`:

```text
COUNT 3336
FNV c1ba162d8a364fb3
SHA256 95a9c0eb185847f1d64d949cef2b1343e85701dcf416c27f3790c31a91c40854.
```

The THM-4256 ray and raw carrier deletion have overlap zero, checked by exact
set intersection.  This is a component statement, not THM-4266's new proof
graph contribution.

Proved THM-4261 and THM-4262 remove the following disjoint subsets of the
post-THM-4256 residual:

```text
THM-4261 COUNT 297 FNV e923d1494185b820
  SHA256 745ef7c8809335e6d6e9623314beff917edc71cfaaaa88e7210ede9dcd97d11b
THM-4262 COUNT 72 FNV 512cbba28e2235fd
  SHA256 b1c89073d9b82351b663e97c18c807f03f3fd2d40ddcfafe038d8cad0535cb2c
UNION COUNT 369 FNV cb13865c00a1a670
  SHA256 61a1c14cf58eeae4010ec0e6b8384e38d24aefa14bac4b70a4aeb7cf5f59c34c.
```

The raw carrier closure contains all 297 THM-4261 edges.  Its overlap with
THM-4262 consists exactly of

```text
(516,688) at g=172, (540,720) at g=180,
COUNT 2, FNV c3a79960d665f57d,
SHA256 3110bd0e1e067bdb0fb9b5ba16f74afd62dc1e6e2388f96afff4970063ef5427.
```

Thus the raw 3,336-edge closure decomposes without triple overlap as

```text
297 inherited THM-4261 + 2 inherited THM-4262 + 3037 THM-4266-new.
```

The exact current high universe and new deletion are

```text
CURRENT R>=688 COUNT 3038
  FNV 4ed8ceb63e3edc9e
  SHA256 4333d7306bb1f8df464b3bd3261559e9d91e6dc88c835e37518206b8a9d0e643
THM-4266 NEW COUNT 3037
  FNV 24b36d7047589076
  SHA256 fcfec867819898ec1a0e1072f27747aec29b6785794328a153bc2b85956ba112.
```

Removing the new set from `E_current` gives

```text
COUNT 177585
FNV 6ce05d05eb01daed
SHA256 009614651bb81e9763b2a9ff4b580497bfb6978a6c69d18cf986346e369374d9
MAX ENDPOINT 688
TOP LAYER (520,688) only.                               (13)
```

As a retained component control, applying the raw 3,336-edge deletion
directly to the post-THM-4256 residual gives 177,655 edges, FNV
`05d884e33afd6c65`, SHA-256
`50b78270a3d16f7c78249f19b3a39ab4b66aa95be1d9820d6a473e85fef12c06`.
The 70 THM-4262 edges outside the raw carrier account for the difference to
the current final census.

## 7. Scope firewall

The statement is relative to the exact post-THM-4261+4262 residual and the
fixed thirty-label pool.  It does not re-claim the 299 overlap edges as new,
prove `(520,688)`, prove any edge below the audited cutoff, produce a
common-active carrier, establish global prefix minimality or monotone activity
in `(q,r)`, prove a new ray, physical entry, or LRC(14).  The 3,037-edge new
set and residual census in `(13)`, rather than the raw 3,336-edge component,
are the exact proof-graph consequence.

## 8. Reproduction

The promoted surface keeps byte-identical local copies of three already
canonical THM-4254 sources so that basename includes remain portable. They
are copies, not new mathematical dependencies. From the repository root,
first verify all 39 frozen artifacts:

```bash
env LC_ALL=C LANG=C shasum -a 256 -c \
  05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/SHA256SUMS
```

Then replay the three compact load-bearing summaries:

```bash
SRC=04-computation/lrc14_three_round_learned_carrier_thm4266
OUT=05-knowledge/results/lrc14_three_round_learned_carrier_thm4266

python3 -B "$SRC/reconstruct_post_thm4254_residual.py" . \
  > /tmp/thm4266-residual.out
cmp /tmp/thm4266-residual.out "$OUT/residual_audit.out"

python3 -B "$SRC/proof_graph_consequence.py" . \
  > /tmp/thm4266-proof-graph.out
cmp /tmp/thm4266-proof-graph.out "$OUT/proof_graph_consequence.out"

python3 -B "$SRC/audit_layer_coverage.py" "$OUT" \
  > /tmp/thm4266-layers.out
cmp /tmp/thm4266-layers.out "$OUT/layer_partition_audit.out"
```

For the independent round-two/three literal-wall control:

```bash
clang++ -std=c++20 -O3 -DNDEBUG \
  -I"$SRC" -I04-computation \
  "$SRC/cleanroom_round2_controls.cpp" -o /tmp/thm4266-cleanroom

/tmp/thm4266-cleanroom \
  05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band \
  "$OUT/full_discovery_416_704_O3.semantic.out" \
  "$OUT/full_discovery_416_700_O3.semantic.out" \
  "$OUT/full_discovery_520_700_O3.semantic.out" \
  "$OUT/full_discovery_384_694_O3.semantic.out" \
  > /tmp/thm4266-cleanroom.out
cmp /tmp/thm4266-cleanroom.out "$OUT/cleanroom_round23_O3.out"
```

The clean-room route constructs literal joint walls rather than using the
primary endpoint atoms. Normal, optimized and sanitizer builds agree. **QED.**

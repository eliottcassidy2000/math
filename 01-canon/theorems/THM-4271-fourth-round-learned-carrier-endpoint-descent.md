---
id: THM-4271
title: "Fourth-round learned-carrier endpoint descent"
status: >
  PROVED RELATIVE TO THM-4266/4267/4269/4270 + FINITE-EXACT + INDEPENDENTLY
  AUDITED WITH DETACHED LITERAL-WALL CONTROLS. Appending 199 masks from the order-minimal full
  discovery prefix at (520,688) gives an 8,518-mask carrier. It closes every
  current post-THM-4270 residual edge with second endpoint at least 671 except
  (256,671) and (384,671). The exact current contribution is 2,419 edges and
  leaves 174,904 residual edges, with those two edges alone on the top layer.
  Neither retained edge, an arbitrary pair, physical entry, nor LRC(14) is
  proved.
source: root/lrc-top-edge-round4/2026-08-27
depends_on:
  - THM-4266-three-round-learned-carrier-endpoint-descent
  - THM-4267-uniform-four-five-outsider-ray-common-deck-closure
  - THM-4269-uniform-five-six-outsider-ray-common-deck-closure
  - THM-4270-uniform-four-primitive-outsider-rays-common-deck-closure
artifact_root: 05-knowledge/results/lrc14_fourth_round_learned_carrier_thm4271
artifact_manifest: 05-knowledge/results/lrc14_fourth_round_learned_carrier_thm4271/SHA256SUMS
artifact_manifest_sha256: 47f39489c19b71c20d9637d44fedf97df0809abb0643d883b77fa272a2033aff
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The complete 5,852,925-repair discovery closes
  (520,688); the exact carrier replay scans every current row through the
  first resistant layer. A detached checker constructs literal joint walls,
  reproduces the full order-minimal prefix, both old seed failures, the
  199-mask transfer, seed closure, and both layer-671 hostiles. O0/O3/
  ASan+UBSan literal transcripts are byte-identical and sanitizer stderr is
  empty. A portable semantic postprocessor reproduces all proof-graph ledgers
  through canonical THM-4270. A hostile referee independently replayed the
  complete discovery, all 162 endpoint-671 rows, literal O0/O3/sanitizer
  controls, and current residual accounting with no blocker.
---

# THM-4271 -- fourth-round learned-carrier endpoint descent

**PROVED RELATIVE TO THM-4266/4267/4269/4270 + FINITE-EXACT + INDEPENDENTLY
AUDITED WITH DETACHED LITERAL-WALL CONTROLS; LRC(14) REMAINS OPEN.**

## 1. Exact statement

Retain the fixed labelled pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290}
```

and, for a finite positive set `A`, put

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14},   alpha=4/63.
```

Let `E_4270` be the exact residual after proved THM-4270. Then, for every
`(q,r) in E_4270` satisfying

```text
r>=671,  (q,r)!=(256,671),(384,671),
```

and every `B in binom(P,9)`, one has

```text
mu(G_(B union {q,r})) >= 4/63.                         (1)
```

The current high-endpoint universe contains 2,421 edges. Exactly the two
displayed rows are retained, so `(1)` closes 2,419 proof-graph-new edges.

## 2. Carrier semantics and inheritance

An eight-mask `R` encodes an eight-subset of `P`. At a fixed pair `(q,r)`, it
is active when

```text
mu(G_((P\R) union {q,r})) >= 4/63.                    (2)
```

A carrier closes `(q,r)` when every `B in binom(P,9)` is disjoint from at
least one pair-active carrier mask. Indeed, for such a mask,

```text
B subset P\R,
G_((P\R) union {q,r}) subset G_(B union {q,r}),       (3)
```

so `(2)` implies `(1)` in the required direction. Activity is re-evaluated
exactly for every row. This is not a common-active deck and no monotonicity in
`q` or `r` is used.

THM-4266 supplies the frozen three-round carrier

```text
count=8,319, FNV=e08b227730f6793c.                    (4)
```

Every carrier FNV is FNV-1a-64 over ordered masks serialized as unsigned
little-endian `u64` words. Edge FNVs serialize `q` then `r` in inherited
residual order. Edge SHA-256 hashes the raw ASCII lines `q,r\n`.

## 3. Complete fourth-round discovery

At the unique post-THM-4267 top edge `(520,688)`, independently enumerate all

```text
binom(30,8)=5,852,925
```

repairs and order the active masks by

```text
(SplitMix64(mask xor 0x4245422842334245),mask).        (5)
```

The exact discovery census is

```text
active=2,504,029, inactive=3,348,896, equalities=0,
active FNV=74b6b666468e2d.                             (6)
```

The first prefix covering all `binom(30,9)=14,307,150` labelled bodies has

```text
size=5,398, FNV=6ab471da88c8e1d1,
body checks=486,391,954, maximum prefix=5,398.         (7)
```

Its predecessor misses the order-minimality witness

```text
{16,85,88,145,168,193,240,252,290},                   (8)
```

and its decisive final repair is

```text
{10,95,120,132,143,176,264,286}.                      (9)
```

Thus minimality is asserted only within the frozen order `(5)`, not among all
possible closing decks.

Exactly 5,199 prefix masks already occur in `(4)`. Appending only the 199
novel masks, in their source-prefix order, gives

```text
novel count=199, FNV=772162cabf9167f9,
round-four carrier count=8,518, FNV=1603e3fe970f8428. (10)
```

The atom/component identity was checked for every appended mask; the combined
ledger is `b0e486cdd9ce6f4e`.

## 4. Exact seed mechanism

At `(520,688)`, the inherited carrier has 5,934 active masks and misses
exactly two bodies:

```text
H1={16,85,88,95,145,168,193,240,252},
H2={16,85,88,168,176,193,240,252,290}.                (11)
```

All 199 appended masks are active at the seed. Exactly one is disjoint from
each old failure:

```text
H1 witness 0x20810e81={8,42,63,80,84,126,190,290},
H2 witness 0x000cc283={8,10,42,63,95,120,143,145}.    (12)
```

The enriched carrier has 6,133 active masks and closes all bodies, with
486,545,771 ordered checks and maximum checked active prefix 6,111. This
identifies why the fourth round works: the new active subdeck supplies a
disjoint witness for each precise obstruction left by `(4)`.

## 5. Exact descent and first new resistance

The full post-THM-4267 residual was reconstructed from the canonical
THM-4266 residual and the semantic `4:5` deletion. Every row in the following
partition was scanned with the exact 8,518-mask carrier:

| endpoint | rows | resistant |
|---:|---:|---:|
| 688 | 1 | 0 |
| 687 | 126 | 0 |
| 686 | 126 | 0 |
| 685 | 127 | 0 |
| 684 | 130 | 0 |
| 683 | 133 | 0 |
| 682 | 135 | 0 |
| 681 | 136 | 0 |
| 680 | 139 | 0 |
| 679 | 143 | 0 |
| 678 | 144 | 0 |
| 677 | 146 | 0 |
| 676 | 149 | 0 |
| 675 | 153 | 0 |
| 674 | 156 | 0 |
| 673 | 159 | 0 |
| 672 | 161 | 0 |
| 671 | 162 | 2 |

Thus exactly 2,424 of the 2,426 post-THM-4267 rows in these layers close. The
other two are the first resistances, with smallest audited stopping data

```text
(256,671): active=4,290, failures=26,
 first={8,80,88,95,132,143,168,240,252},
 best disjoint={42,63,85,170,193,264,286,290},
 margin=-89508336872746176/43867507598953758720;

(384,671): active=5,988, failures=1,
 first={8,80,88,95,145,168,193,252,264},
 best disjoint={63,85,120,126,132,143,176,286},
 margin=-721295701293847680/65801261398430638080.      (13)
```

These are failures of this fixed carrier compression only. They are not
counterexamples to the target inequality or to the endpoint method.

## 6. Detached literal-wall audit

The independent checker does not call the endpoint-atom or cocycle activation
routines. For `(520,688)` it constructs the literal joint wall grid

```text
grid=784,369,854,908,640, joint cells=9,437.           (14)
```

It rechecks all 5,398 prefix masks, the full body recursion, and prefix
order-minimality. The least prefix literal numerator is `12,726,097,524` at

```text
{8,10,15,30,85,170,190,240}.                          (15)
```

It then independently reproduces both failures in `(11)`, all 199 activation
signs and the two one-witness disjointness counts in `(12)`, and the zero-
failure 6,133-active-mask seed scan. Finally it reproduces both rows in
`(13)`, including active counts, failure counts, first witnesses, best
disjoint repairs, and exact margin ratios.

Ordinary, optimized, and ASan+UBSan literal transcripts are byte-identical:

```text
SHA256=2a846200c31d1b0e94a3e262ab42e11f1aaf25dd1c425c4793486014c1008fd2.
```

Sanitizer stderr is empty, with SHA-256
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.

## 7. Exact live proof graph

Relative to the canonical post-THM-4267 residual, the raw carrier component
has

```text
count=2,424,
FNV=3ec97f3ae7e5d142,
SHA256=44f00e2b22adbf071eaebe3edf474337abf789ac9dad867cc48d0e157cfebc94.
                                                               (16)
```

It has zero overlap with THM-4269's 53-edge `5:6` component. Its overlap with
THM-4270 is exactly

```text
(588,672), (595,680), (600,675), (608,684), (616,672),
count=5,
FNV=7484ad77ffddd129,
SHA256=7eafa6774730c49213c08eb7b595952ca64c4e224db30b348933bfa422f82060.
                                                               (17)
```

Therefore the exact new contribution after canonical THM-4270 is

```text
count=2,419,
FNV=ff169664a6750abe,
SHA256=7980354ba5d9dde9ce994eda992bf7030d79bd0562cf9c8f2742d4bb53653e89.
                                                               (18)
```

Removing only `(18)` from THM-4270's residual gives

```text
count=174,904,
FNV=b3db855040bcf19e,
SHA256=07d84c1572baeb89d9f88e095788e52e3916dab6074f89cf3b164e5ebea3a5a6,
maximum endpoint=671,
top layer={(256,671),(384,671)}.                       (19)
```

## 8. Scope and reproduction

This theorem is relative to the fixed thirty-label pool and the exact current
residual. It proves no retained row in `(13)`, no pair below the audited
boundary, no common-active carrier, no global deck minimum, no physical entry
of an arbitrary row into this pool, and no LRC(14).

The frozen artifact manifest has SHA-256
`47f39489c19b71c20d9637d44fedf97df0809abb0643d883b77fa272a2033aff`.
Full commands are in
`05-knowledge/results/lrc14_fourth_round_learned_carrier_thm4271/REPRODUCTION.md`.
**QED.**

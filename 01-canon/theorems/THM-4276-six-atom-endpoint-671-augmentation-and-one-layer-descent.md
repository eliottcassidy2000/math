---
id: THM-4276
title: "Six-atom endpoint-671 augmentation and one-layer descent"
status: >
  PROVED RELATIVE TO THM-4266/4271 + FINITE-EXACT + DETACHED
  LITERAL-WALL AUDIT PASS. The two retained endpoint-671 rows leave exactly
  27 failed-body obligations for the 8,518-mask THM-4271 carrier. The exact
  full-repair atlas has 330 realizable nonempty response patterns, and an
  independent packing/cover certificate proves that six added masks are
  necessary and sufficient. Appending the certified six masks gives an
  8,524-mask carrier. It closes both endpoint-671 rows and 161 of the 163
  current endpoint-670 rows. The exact contribution after canonical THM-4271
  is 163 edges and leaves 174,741 residual edges, with (256,670) and
  (384,670) alone on the top layer. Neither retained edge, an arbitrary pair,
  physical entry, nor LRC(14) is proved.
source: root/lrc-top-edge-round4/compact-round5/2026-08-27
depends_on:
  - THM-4266-three-round-learned-carrier-endpoint-descent
  - THM-4271-fourth-round-learned-carrier-endpoint-descent
artifact_root: 05-knowledge/results/lrc14_six_atom_endpoint_671_augmentation_thm4276
artifact_manifest: 05-knowledge/results/lrc14_six_atom_endpoint_671_augmentation_thm4276/SHA256SUMS
artifact_manifest_sha256: 64fb5830696be89fd5b9e1e55745b0bd21dfab6f7fdefabaf0fb7eda6ca8336f
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT CANDIDATE. Two independent complete 5,852,925-repair
  discoveries freeze order-minimal source decks at (256,671) and (384,671).
  A separate full-universe pass freezes all 330 nonempty response patterns on
  the 27 inherited obligations. A portable solver-free certificate checks a
  six-obligation packing and a matching six-pattern cover. The primary replay
  reconstructs the THM-4271 carrier, appends exactly the six certified masks,
  closes both seeds, and scans the complete current endpoint-670 layer. A
  detached direct-joint-wall checker reproduces carrier identity, selected
  patterns, seed closure, both full prefixes, and both hostile controls.
  O0/O3/ASan+UBSan literal outputs are byte-identical and sanitizer stderr is
  empty. A portable semantic postprocessor reproduces the proof graph after
  canonical THM-4271.
---

# THM-4276 -- six-atom endpoint-671 augmentation and one-layer descent

**PROVED RELATIVE TO THM-4266/4271 + FINITE-EXACT + DETACHED
LITERAL-WALL AUDIT PASS; LRC(14) REMAINS OPEN.**

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

Let `E_4271` be the exact residual after canonical THM-4271. Then, for every
`(q,r) in E_4271` satisfying

```text
r>=670,  (q,r)!=(256,670),(384,670),
```

and every `B in binom(P,9)`, one has

```text
mu(G_(B union {q,r})) >= 4/63.                         (1)
```

The current universe in these layers has 165 edges: the two retained
endpoint-671 rows and 163 endpoint-670 rows. Exactly the two displayed
endpoint-670 rows are retained, so `(1)` closes 163 proof-graph-new edges.

## 2. Carrier semantics and exact inheritance

An eight-mask `R` encodes an eight-subset of `P`. At a fixed pair `(q,r)`, it
is active when

```text
mu(G_((P\R) union {q,r})) >= 4/63.                    (2)
```

A carrier closes `(q,r)` when every `B in binom(P,9)` is disjoint from at
least one pair-active carrier mask. For such a mask,

```text
B subset P\R,
G_((P\R) union {q,r}) subset G_(B union {q,r}),       (3)
```

so `(2)` implies `(1)` in the required direction. Activity is recomputed
exactly for every pair. No common-active-deck or endpoint-monotonicity claim
is used.

The replay reconstructs every inherited carrier stage before augmentation:

| stage | masks | ordered-mask FNV |
|---|---:|---:|
| THM-4266 base band | 4,675 | `ce4e76ec11df057c` |
| first round | 4,733 | `a7b046289655c733` |
| second round | 7,986 | `baef1d2f49444638` |
| third round | 8,319 | `e08b227730f6793c` |
| THM-4271 fourth round | 8,518 | `1603e3fe970f8428` |

The last stage appends the 199 masks novel in the canonical order-minimal
`(520,688)` discovery prefix. Carrier FNV is FNV-1a-64 over ordered masks
serialized as unsigned little-endian `u64` words.

## 3. The exact 27-obligation atlas

At `(256,671)`, the THM-4271 carrier has 4,290 active masks and misses exactly
the following 26 bodies, in the displayed bit order `0,...,25`:

```text
06166401,07067001,07106409,07126088,07126401,07162401,
07163400,07166008,0d107401,0d10e401,0d146401,0d186401,
0d246401,0d506401,0d906401,0f106401,0f142401,0f142408,
0f143400,0f146008,15923400,17162400,17922008,1d106401,
1d902401,1f142400.
```

At `(384,671)`, it has 5,988 active masks and misses exactly one body,
`0x0d186401`; this distinct row-obligation is bit 26. Thus there are exactly
27 obligations even though one labelled body occurs at both rows.

For each of all

```text
binom(30,8)=5,852,925
```

possible repair masks, define its 27-bit response pattern to contain an
obligation precisely when that repair is active at the obligation's pair and
disjoint from its failed body. Exact atom activation gives

```text
(256,671): active=1,721,339,
(384,671): active=2,444,056,
repairs with nonempty response=19,156,
distinct nonempty response patterns=330.              (4)
```

The atlas freezes the numerically least eight-mask realizing each response
pattern. Its ordered `(pattern,least,cover-count)` ledger is
`7811c918298740dd`; the raw atlas SHA-256 is
`3717bc24ed8cb6167cd0194606a453b443ece91b7677ae9551573d9e864ab274`.

## 4. A conceptual and exact minimum-six certificate

The following six obligations form a packing. Each row gives its obligation
bit, body mask, and labelled body:

| bit | body mask | labelled body |
|---:|---:|---|
| 0 | `0x06166401` | `{8,80,88,95,132,143,168,240,252}` |
| 9 | `0x0d10e401` | `{8,80,88,95,120,168,193,252,264}` |
| 12 | `0x0d246401` | `{8,80,88,95,143,170,193,252,264}` |
| 16 | `0x0f142401` | `{8,80,88,143,168,193,240,252,264}` |
| 19 | `0x0f146008` | `{16,88,95,143,168,193,240,252,264}` |
| 20 | `0x15923400` | `{80,85,88,132,168,190,193,252,286}` |

The exhaustive 330-pattern atlas verifies that no realizable response pattern
covers two of these packing obligations: all fifteen pairs are incompatible,
and every response has at most one packing hit. Therefore each added repair
pays for at most one of these six obligations. Any augmentation that closes
both endpoint-671 rows consequently needs at least six repairs.

The matching six-repair cover is, in frozen append order,

| repair | response pattern | covered bits |
|---:|---:|---|
| `0x00289285` | `0x026a0080` | `7,17,19,21,22,25` |
| `0x0260812c` | `0x05904d00` | `8,10,11,14,20,23,24,26` |
| `0x18689040` | `0x000000bd` | `0,2,3,4,5,7` |
| `0x20c0c124` | `0x02270060` | `5,6,16,17,18,21,25` |
| `0x302c1006` | `0x0000e21c` | `2,3,4,9,13,14,15` |
| `0x30580888` | `0x00001002` | `1,12` |

Their response union is exactly `0x07ffffff`, all 27 obligation bits. The
alternating `(repair,response)` ledger is `5e7cc39cf6cee4be`. Thus six is the
exact minimum augmentation size for closing these two rows relative to the
fixed THM-4271 carrier. This is not a global minimum-carrier or minimum-deck
claim.

All six repairs are absent from the inherited carrier. Appending exactly

```text
0x00289285,0x0260812c,0x18689040,
0x20c0c124,0x302c1006,0x30580888                 (5)
```

gives

```text
compact carrier count=8,524, FNV=5ddb84a44f5d2ad7.    (6)
```

## 5. Exact seed closure and complete discovery controls

With the compact carrier `(6)`, direct body recursion gives

| pair | active masks | bodies | failures | checks | max checks |
|---|---:|---:|---:|---:|---:|
| `(256,671)` | 4,296 | 14,307,150 | 0 | 527,844,189 | 4,296 |
| `(384,671)` | 5,994 | 14,307,150 | 0 | 483,650,358 | 5,990 |

Thus the six added masks close both retained THM-4271 rows.

As hostile and order-minimal controls, two separate complete repair
discoveries also freeze full closing prefixes in the SplitMix order with seed
`0x4245422842334245`:

| pair | all active | prefix | prefix FNV | body checks | predecessor witness | decisive repair |
|---|---:|---:|---:|---:|---|---|
| `(256,671)` | 1,721,339 | 13,086 | `06c1f44ec0661d8c` | 527,773,826 | `{8,80,88,143,168,193,240,252,264}` | `{16,63,85,95,145,170,190,290}` |
| `(384,671)` | 2,444,056 | 6,986 | `3a02a5774d3641ab` | 483,641,885 | `{8,80,88,95,145,168,193,252,264}` | `{10,15,85,126,143,170,240,286}` |

Each discovery scans all 5,852,925 repairs, has zero activation equality, and
closes all 14,307,150 bodies. These large prefixes are controls only; they are
not appended to `(6)`.

## 6. One-layer descent and first new resistance

The exact post-THM-4271 residual has two rows at endpoint 671, both closed
above. The entire next layer was then scanned with `(6)`:

| endpoint | current rows | resistant | closed |
|---:|---:|---:|---:|
| 671 | 2 | 0 | 2 |
| 670 | 163 | 2 | 161 |

The endpoint-670 aggregate is

```text
active sum=1,331,051,
ordered body checks=69,626,501,969,
maximum checked active prefix=3,851,
row ledger=90e990cdced5a71b.                            (7)
```

The first two resistances are

```text
(256,670): active=2,172, failures=1,659,
 first=0x01147409={8,16,80,85,88,95,143,168,193},
 best disjoint=0x2ee00040={40,170,176,190,240,252,264,290},
 atom margin=-472438974077602560/21901065641802547200;

(384,670): active=3,851, failures=46,
 first=0x043e6400={80,88,95,132,143,145,168,170,252},
 best disjoint=0x3a001809={8,16,84,85,240,264,286,290},
 atom margin=-1258637401934962560/32851598462703820800. (8)
```

These are failures of the fixed 8,524-mask compression only. They are not
counterexamples to `(1)` outside its stated domain or to the target
inequality.

## 7. Detached literal-joint-wall audit

The independent checker constructs fresh literal joint walls and never calls
the endpoint-atom or cocycle activation routines. It independently reproduces
the complete carrier construction through `(6)`, the 26+1 inherited
obligations, all six selected response patterns, response union
`0x07ffffff`, both zero-failure compact seed scans, and both full discovery
prefixes including order-minimality.

For the two controls in `(8)`, the direct literal margins are

```text
(256,670): -421820512569288/19554522894466560
           = -836945461447/38798656536640;
(384,670): -374594464861596/9777261447233280
           = -1486485971673/38798656536640.            (9)
```

Both cross-multiply exactly to the atom ratios in `(8)`. The literal checker
also records 8,985 and 9,241 joint cells for the two endpoint-671 full-prefix
controls, with least literal margins respectively `262035527418` and
`292393039680`; all body recursions close.

Ordinary, optimized, and ASan+UBSan literal transcripts are byte-identical:

```text
SHA256=977497ec1418279141b7486994c531a22f2a7f1a8173f1a58b04e83c962d20da.
```

Sanitizer stderr is empty, with SHA-256
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.

## 8. Exact proof-graph novelty after canonical THM-4271

The incoming canonical residual is

```text
post-THM-4271 count=174,904,
FNV=b3db855040bcf19e,
SHA256=07d84c1572baeb89d9f88e095788e52e3916dab6074f89cf3b164e5ebea3a5a6.
                                                               (10)
```

The exact current contribution of this theorem is the two endpoint-671 rows
plus 161 endpoint-670 rows:

```text
count=163, layers={671:2,670:161},
FNV=edc377490e3f58bb,
SHA256=49960a608af498bf385808f3fe234b75923ae203ea9a6e95875e61df1c43de26.
                                                               (11)
```

Removing `(11)` and no other edge gives

```text
count=174,741,
FNV=f13745b05320f83c,
SHA256=51d5723b146cb108a2e11627924a2fd6af46435564e2460ab78af936bfb12dd0,
maximum endpoint=670,
top layer={(256,670),(384,670)}.                       (12)
```

## 9. Scope and reproduction

This theorem is relative to the fixed thirty-label pool, the inherited
carrier semantics, and the exact canonical THM-4271 residual. It proves no
retained row in `(8)`, no pair below the audited boundary, no common-active
carrier, no global deck minimum, no physical entry of an arbitrary row into
this pool, and no LRC(14).

The frozen artifact manifest has SHA-256
`64fb5830696be89fd5b9e1e55745b0bd21dfab6f7fdefabaf0fb7eda6ca8336f`.
Full commands are in
`05-knowledge/results/lrc14_six_atom_endpoint_671_augmentation_thm4276/REPRODUCTION.md`.
**QED.**

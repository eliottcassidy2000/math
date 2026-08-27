---
id: THM-4281
title: "Global frozen 421-mask joint-deck collapse, carrier closure through 664, and residual top 663"
status: >
  PROVED RELATIVE TO THM-4276/4277 + FINITE-EXACT + DETACHED
  NEWCOMER-LITERAL-WALL AUDIT PASS. One frozen 421-mask joint deck disjoint-covers every
  labelled nine-body and is common-active on exactly 148,063 of the 172,322
  canonical post-THM-4277 residual pairs. A nested six-mask carrier
  augmentation separately closes every residual row in endpoints 670 through
  664. After the exact 1,032-row overlap is removed, the theorem contributes
  148,099 edges and leaves 24,223, with top layer
  (256,663),(366,663),(520,663). No physical entry or LRC(14) follows.
source: root/lrc-carrier-cegar-descent/2026-08-27
depends_on:
  - THM-4276-six-atom-endpoint-671-augmentation-and-one-layer-descent
  - THM-4277-uniform-two-dimensional-outsider-rectangle-common-deck-closure
related:
  - THM-4266-three-round-learned-carrier-endpoint-descent
  - THM-4271-fourth-round-learned-carrier-endpoint-descent
artifact_root: 05-knowledge/results/lrc14_joint421_global_common_carrier_thm4281
artifact_manifest: 05-knowledge/results/lrc14_joint421_global_common_carrier_thm4281/SHA256SUMS
artifact_manifest_sha256: 3acd31a64e62ffefb44e8dcd103cb94ed30e4e1ae7d50eb62c01321c71695f8b
primary_script: 04-computation/lrc14_joint421_global_common_carrier_thm4281/scout_joint421_post4277_primary.cpp
primary_output: 05-knowledge/results/lrc14_joint421_global_common_carrier_thm4281/scout_joint421_post4277_primary_O3.out
primary_script_sha256: 4bd3b70eaada775444b6d8e1643f23fcf5d58d23820c3200e00047d479146417
primary_output_sha256: d054a6d5d13c94c0048acb54910ea38c0b6127d7f57986c56712abf637943375
literal_script: 04-computation/lrc14_joint421_global_common_carrier_thm4281/verify_joint421_literal_post4277.cpp
literal_output: 05-knowledge/results/lrc14_joint421_global_common_carrier_thm4281/verify_joint421_literal_post4277_O3_NDEBUG.out
literal_script_sha256: 72f66f5015549e2341d888f85fe781bad3a2bf61b86b0c78aff1b35f194f0503
literal_output_sha256: 3081efc81b43cffe912d38f22e8b2e09719980cd0291f42551ee386a2b3100b2
endpoint_literal_script: 04-computation/lrc14_joint421_global_common_carrier_thm4281/verify_joint421_literal_r670.cpp
endpoint_literal_output: 05-knowledge/results/lrc14_joint421_global_common_carrier_thm4281/verify_joint421_literal_r670_O3_NDEBUG.out
endpoint_literal_script_sha256: 648cd538c751f39ddc594ccc596f6285c80ffb7c05c0bb4ce6d0c3fd0d71cecb
endpoint_literal_output_sha256: 88fa48d3130dbe1a1948a22c79fa5f6c7fa3003a20cbd1c4a365b99d71b0a530
proof_graph_script: 04-computation/lrc14_joint421_global_common_carrier_thm4281/verify_common_complement_carrier_overlap.py
proof_graph_output: 05-knowledge/results/lrc14_joint421_global_common_carrier_thm4281/verify_common_complement_carrier_overlap.out
proof_graph_script_sha256: 347792b5d679df2bfc0996299119b04f7e9f2e0ff74d84dd82174cbd66f8e3db
proof_graph_output_sha256: 8ed72b405f3dbcbfa0355544647ea3c8f15f93eb30bc8dc872096a8df5dea20c
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. A primary endpoint-prefix engine and a detached
  literal-grid newcomer engine agree byte-for-byte on the exact 148,063-pair
  common set. The literal route checks all 62,334,523 common cells, reproduces
  the first inactive repair and signed margin ratio on every one of 24,259
  complement pairs, and independently reruns the complete body cover.
  Complete response quotients prove the four carrier augmentation minima;
  a detached literal descent closes every stated seed and audits the
  endpoint-667 packing lower bound. Ordinary, optimized, NDEBUG, and UBSan
  controls agree on the consequence-bearing transcripts, and exact normal
  and optimized Python replays reproduce the set union and proof graph.
---

# THM-4281 -- global frozen 421-mask joint-deck collapse, carrier closure through 664, and residual top 663

**PROVED RELATIVE TO THM-4276/4277 + FINITE-EXACT + DETACHED
NEWCOMER-LITERAL-WALL AUDIT PASS; LRC(14) REMAINS OPEN.**

## 1. Statement

Retain the labelled pool and threshold

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290},
G_A={x in R/Z:min_(a in A)||ax||>=1/14},   alpha=4/63.
```

Let `U_4277` be the exact canonical residual after THM-4277. There is an
explicit ordered deck `E` of 421 eight-subsets of `P` and an exact subset
`A subset U_4277` such that

```text
|U_4277|=172,322,              |A|=148,063,           (1)
```

and, for every `(q,r) in A` and every `B in binom(P,9)`,

```text
mu(G_(B union {q,r})) >= 4/63.                          (2)
```

The same deck, appended to the 8,524-mask THM-4276 carrier and followed by
six exact carrier-relative repairs, closes every canonical residual row with

```text
r>=664.                                                (3)
```

The common-deck and carrier conclusions overlap. Their exact union contains
148,099 proof-graph-new edges, leaving

```text
count=24,223,
maximum endpoint=663,
top layer={(256,663),(366,663),(520,663)}.              (4)
```

## 2. The 421-mask joint deck

The frozen ordered deck has

```text
count=421,
ordered-mask FNV=20d63dd42fe8150e,
raw SHA256=fcde04de68ab614743176ed153e0db2754cf878470c96f525d29c07a88388bc2.
                                                               (5)
```

Exact enumeration of all `binom(30,9)=14,307,150` labelled bodies proves

```text
for every B in binom(P,9), some R in E has B intersect R=empty. (6)
```

The scan has zero failures, 405,170,384 ordered checks, and maximum checked
prefix 421. Thus whenever every mask in `E` is active at `(q,r)`, choosing the
repair in `(6)` gives

```text
B subset P\R,
G_((P\R) union {q,r}) subset G_(B union {q,r}),        (7)
```

which proves `(2)` in the required direction.

An exact body-incidence census gives 209,413,820 deck/body incidences,
multiplicity range `1..76`, and 3,512 bodies of multiplicity one. Exactly one
deck mask has no private body: zero-based index 318 (line 319),
`0x003c900c`. Deleting it leaves a 420-mask body cover with ordered FNV
`e72c3c8b50ec6c6e`, and a detached literal replay proves that this 420-mask
core has exactly the same 148,063 common pairs as `E`. Thus the extra mask is
not load-bearing for either `(2)` or the global
common-set census. It is nevertheless endpoint-load-bearing inside the frozen
joint deck: at the two endpoint-670 rows it meets 28 of the 1,705 inherited
obligations, and five of those responses are private to this mask. Deleting it
would leave exactly five obligations at `(256,670)` uncovered; their ordered
ledger is `3a11d5bef4e6318d`. This separates body-cover/common-set essentiality
from endpoint-response essentiality.

The deck was discovered by joining the THM-4277 rectangle-common repair
family to the 1,705 labelled failed-body obligations at the two retained
endpoint-670 rows. Exact monotone deletion from a 434-mask joint cover leaves
the deck `(5)` while preserving both `(6)` and every endpoint-670 obligation.
This is an explicit frozen joint-cover construction, not a proof of a globally
minimum or unique joint deck. In particular the 420-mask core above makes
clear that `E` is not inclusion-minimal as a body cover or as a global-common
deck. The universal body-cover lower bound 52 from THM-4277 still applies.

## 3. Exact global common-activity set

Define

```text
A={(q,r) in U_4277: every R in E is active at (q,r)}.  (8)
```

The exact ordered ledgers are

```text
U_4277: count=172,322, FNV=30b2a7e597ac548c,
         SHA256=7228658eae4067db4bbcceb6c9b1ccebf1bd3e6f128e202ea184854acc53f309;
A:       count=148,063, FNV=465fb756a183167e,
         SHA256=412b942759ed7afde1bffaaeabc9e6ae31b8b8bc25f26d73c71712523a057aaa;
U_4277\A: count=24,259, FNV=78b212469c336f37,
         SHA256=77f3d21d127bb5e21f583556314f74032271e3a4903696415885308b363624ef.
                                                               (9)
```

The primary endpoint-prefix engine tests the masks in frozen deck order and
stops a noncommon pair at its first negative exact margin. It performs

```text
64,116,779 tested activation cells,
62,334,523 common-set cells,
zero equality among the tested cells.                              (10)
```

In particular every cell used to prove `(2)` has a strictly positive margin.
No statement is made about equality in the skipped suffix of a pair already
certified noncommon.

The weakest normalized common-set gap occurs at

```text
(q,r)=(35,315),
R=0x22560801={8,84,132,143,168,176,240,290},
mu(G_((P\R) union {35,315}))-4/63
  =3569/126674718170 > 0.                              (11)
```

The common set occupies 633 of the 639 endpoint layers present in `U_4277`.
The six layers with no common pair are

```text
25,50,96,100,210,670.                                  (12)
```

Equation `(12)` is only an exact census; no endpoint monotonicity or scale
transfer is inferred from it.

## 4. Detached literal newcomer audit

The detached engine retains the already audited fixed-pool wall geometry but
does not call the primitive pair, endpoint cocycle, atom mass, or primary
prefix cache. For each literal pair it forms

```text
L=lcm(D,14q,14r),   D=lcm(14s:s in P)=18,241,159,416,480,
                                                               (13)
```

constructs the safe intervals of the actual speeds `q,r` on the `L`-grid,
intersects those interval lists, and integrates their literal prefix over the
fixed-pool repair components.

It checks all `148,063*421=62,334,523` cells of `A` and finds every one
active. For each of the 24,259 complement pairs it independently reaches the
same first inactive deck position and proves equality of the signed normalized
margin with the primary route. It also repeats the full body scan `(6)`.
The resulting common-pair file is byte-identical to the primary file in `(9)`.

At `(11)` the literal arithmetic gives

```text
raw margin/grid=32377968/18241159416480
               =32121/18096388310,                    (14)
```

which is the same normalized raw ratio as the primary value
`142786838880/80443513026676800`; division by 63 gives `(11)`.

## 5. Endpoint-670 bridge

The 421 masks are disjoint from the 8,524-mask THM-4276 carrier. Appending
them gives

```text
carrier count=8,945, FNV=3212efa05dd18c00.             (15)
```

The THM-4276 carrier alone misses 1,659 bodies at `(256,670)` and 46 at
`(384,670)`. Their ordered ledgers are respectively
`970f004b0f9e2edb` and `6c4fd2dab94cf6a9`. Exactly 397 and 367 masks of `E`
are active at the two rows. Their response covers all 1,705 obligations with
7,030 incidences and ledger `4ab5ca017f17efeb`.

Consequently the augmented carrier closes both rows:

| pair | active masks | body failures | checks | max prefix |
|---|---:|---:|---:|---:|
| `(256,670)` | 2,569 | 0 | 613,731,713 | 2,497 |
| `(384,670)` | 4,218 | 0 | 541,846,700 | 4,205 |

This is a carrier statement. The deck is not common-active at either
endpoint-670 row.

## 6. Nested exact repair quotients

Starting from `(15)`, exact layer scans expose four isolated hostile rows.
At each row the response universe consists of all `binom(30,8)=5,852,925`
rank-eight masks, with activity recomputed at that exact pair. The exact
carrier-relative minima are:

| hostile pair | inherited failures | minimum | ordered added masks | resulting carrier/FNV |
|---|---:|---:|---|---|
| `(256,669)` | 4 | 2 | `0x003884c8,0x00c4c124` | `8,947 / 0b246d6377503ae7` |
| `(256,667)` | 18 | 2 | `0x02a05206,0x10e05240` | `8,949 / 1496ebecca72684b` |
| `(520,665)` | 1 | 1 | `0x0004409f` | `8,950 / f07022300e266930` |
| `(256,664)` | 2 | 1 | `0x00c0c125` | `8,951 / 188f82ab9dd1695a` |

For `(256,669)`, the complete response quotient has ten nonempty classes and
no full-response class; the displayed pair has response union `0xf`. For
`(256,667)`, 118 nonempty classes reduce to five maximal classes. Obligations
13 and 16 form a packing, while the displayed two masks have response patterns
`0x11c07` and `0x2efff`, whose union is `0x3ffff`. A detached literal scan of
all `binom(17,8)=24,310` masks disjoint from the union of those two packed
bodies finds zero active candidate. The final two rows have explicit active
one-mask full responses. Thus the minima are exact relative to the stated
nested carrier at each row; they are not unconditional minimum-carrier claims.

The complete layer accounting is

| endpoint | current rows | closed by the nested carrier |
|---:|---:|---:|
| 670 | 2 | 2 |
| 669 | 168 | 168 |
| 668 | 170 | 170 |
| 667 | 176 | 176 |
| 666 | 179 | 179 |
| 665 | 183 | 183 |
| 664 | 190 | 190 |

The final endpoint-664 scan has active sum 1,660,533, zero resistant rows,
80,180,889,036 ordered body checks, maximum prefix 6,307, and row ledger
`2973e6c6d40d5af8`. No row at endpoint 663 is claimed by this carrier scan.

## 7. Exact proof-graph union

Let `C` be the 1,068 rows in endpoints 670 through 664 closed by the nested
carrier, and retain `A` from `(8)`. Exact set arithmetic gives

```text
|A|=148,063,               |C|=1,068,
|A intersect C|=1,032,     |C\A|=36,
|A union C|=148,099.                                      (16)
```

The ordered carrier-only ledger has FNV `a8027022c6436919` and SHA-256
`b14f2cb2ff92ea4f14676e3b6cf6b63fcff1b3378d2ece0461a9a4c89bccfd0e`.
The union ledger is

```text
FNV=0308f3f6b8d7a23e,
SHA256=ed993db4e2b53a416e6031696c4a0b01dfd3f9b9e47fb2607f5f191059752836.
                                                               (17)
```

Deleting `(17)` from `(9)` gives the residual in `(4)`, with

```text
FNV=80ec0687d8c7dba7,
SHA256=75a3c7616c982538363c7801ed2dab3fe9aa775ab601f7a7119dd9fb5d301552.
                                                               (18)
```

Normal and optimized exact Python paths reproduce every set in `(16)--(18)`
and the full 639-layer partition.

## 8. Scope

This theorem is relative to the fixed labelled pool, canonical residual, and
repair-carrier semantics. A pair outside `A` merely has at least one inactive
mask in this particular deck; it is not a counterexample to the target
inequality. The theorem proves no common activity at either endpoint-670 row,
no arbitrary neighbouring pair, no scale monotonicity, no global deck or
carrier minimum, no physical entry of an arbitrary fourteen-speed set into
the fixed pool, and no LRC(14).

Full commands and frozen hashes are in
`05-knowledge/results/lrc14_joint421_global_common_carrier_thm4281/REPRODUCTION.md`.
**QED.**

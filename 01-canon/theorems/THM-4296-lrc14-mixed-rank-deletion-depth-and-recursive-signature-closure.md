---
id: THM-4296
title: "LRC(14) mixed-rank deletion depth, endpoint-626 descent, and recursive signature closure"
status: >
  PROVED RELATIVE TO THM-4287 + THM-4295 + FINITE-EXACT + DETACHED LITERAL-WALL
  AUDITS PASS. Active deletion masks form an upset at every fixed pair, so
  deletion depth is a rank-free certificate invariant. An exact rank-eight
  exchange and mixed rank-eight/rank-nine continuation give a 9,019-mask
  carrier closing 70 current rows through endpoint 626. Separately, 110
  singleton-signature decks close 1,219 rows and one two-replacement deck
  closes the 36-row index-19 ideal. This theorem's three-node union has 1,321
  rows. After adjoining THM-4295's independent index-294 and index-372 ideal
  nodes and auditing every overlap, the typed union has 1,324 rows and leaves
  21,323, maximum endpoint 626, on (100,626) and (256,626). No common deck,
  physical entry, or LRC(14) follows.
source: opus-lrc14-incoming-20260831
depends_on:
  - THM-4287-repaired-carrier-endpoint-637-descent
  - THM-4295-lrc14-endpoint-636-minimum-fourteen-exchange-and-recursive-signature-ideals
  - THM-4254-fixed-ceiling-band-signed-endpoint-cocycle-cascade
related:
  - THM-4286-signature-response-nonfactorization-and-two-deck-surgeries
  - THM-4283-endpoint-644-carrier-response-and-signature-fibre-surgery
  - THM-4282-inactive-signature-deck-surgery-endpoint-663
artifact_root: 05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296
artifact_manifest: 05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296/SHA256SUMS
artifact_manifest_sha256: c3612e076e3b7b0965f339c6b856691e85719921d50b5081f9b40635aac76ca2
primary_scripts:
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/detached_exchange_audit.cpp
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/exact_cover.py
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_response_atlas.cpp
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_detached_hostile_survivor.cpp
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_rank9_mixed_response_atlas.cpp
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r629_mixed_response_atlas_detached.cpp
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/mixed_carrier_descending_detached.cpp
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/singleton_success_literal_audit.cpp
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/detached_j19_audit.cpp
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/proof_graph_union_independent.cpp
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/incoming_thm4295_typed_union.py
  - 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/incoming_thm4295_typed_union_independent.cpp
primary_script_sha256:
  - 4c5176be7b7e56cdfc248fa8b32fb5362741c7b12efe64766eb86d1c2bb02b18
  - 70b37156edf552c1c132acab52bb1a31f8e056d56e271a958a2a6de03f990738
  - 335f38142d668965eb087da201307b0033214eb460443568a5067cc8e8436b58
  - f9984e19cabb7d32c962656ebf658f0b0512be976825d90f1cda90f8a6cd5abe
  - 90a315809c533670a7b4d4d498ec2835e1d4c5c50a4b4e760e8aff22b60e40b6
  - d68b0c1c18202bd6599c93687d9252217af60fb5881fc01e30bb32adb3985b3a
  - b3f143dfc0fe0685f1be802350db2abfffa57e7cdd7b27e2becca1c7372ef3a2
  - 29c8dc77a7d349a3967e05f9738e9ceb1616819527482c00c4f2c33b5d5cb6fb
  - b8f83eb24cc187aff3aaf11ed63d218168dbdc28d700ae6eb30790021393824f
  - da54250bf0d0169e1fd4fb1df2c8b9e2ec5e5aa03d7562087ce8e4339c310ca9
  - b79f10f86c9f30575269471e3bb6074905f17242c7997162ab6b7581677166a2
  - 74ae1a050ac648671c18133b9845b5ceaaf8e091ab29eb967b8b9c5f82d69d69
primary_output_sha256:
  - 4c0dd784d05bdafed4581562d262f93d2adc1049e7f6debe0461a9dbc1f2892c
  - 6dc6ce271cfa1d0524fe9cac92b6642c2f33706fc02645e96254e710a74cfcc9
  - db284af0ca5fcf1c2449ff90159be404353973fd5e73a00b6ca8670a0f62b6ba
  - 684f197c18ec31b23cedad64eb68c088c02da1664dc4cdcbb6f2877fc87a5105
  - b474e13f538298d0a410242cad356a26278361648b82f8c18303f0a91a095532
  - c0253066433f7aa133afeefc38e6440e400f7adcba11f5267a8fd703d78e2f9b
  - 390bc2844b9fa24250ef2e4e19be02cda2702563c2673f8fdbb32f010ca50d70
  - f7a43e984022f70f60a0b07beee21f2d77c2192aeb7b624e4da506bc7edcab2e
  - 4ea0f1d5ef6f23aabdaf43b744fa39e7dae17942e459b25198391ded32afbcd5
  - 357e1331164b0fe73d1b25f16a38bba89bd52db31800c9f8f0684ccb498ebfc5
  - 39baff722de49d6059de460bf0070fb906d1247d809b0602fb9b84be33508af5
  - 2dbef786bfbf5bd193ec680f8f9085ad25602d47074aaa44974583126603dd12
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. Exact response quotients,
  rational dual lower bounds, literal-wall activity checks, exhaustive
  labelled body scans, O0/O3 controls, and an independent raw-scan proof-graph
  consumer agree. A second Python/C++ pair audits the cross-theorem typed join
  with THM-4295. MISTAKE-532 removes noncanonical cross-row weakest-margin
  annotations from the proof surface; activity signs, covers, and residuals
  are unchanged.
---

# THM-4296 -- LRC(14) mixed-rank deletion depth, endpoint-626 descent, and recursive signature closure

**PROVED RELATIVE TO THM-4287 + THM-4295 + FINITE-EXACT + DETACHED LITERAL-WALL
AUDITS PASS. LRC(14) REMAINS OPEN.**

## 1. Statement and inherited universe

Retain THM-4254's labelled thirty-speed pool `P`, threshold
`alpha=4/63`, and notation

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14}.
```

Retain also THM-4287's ordered 421-mask joint deck `E` and exact current
residual `U`:

```text
|U|=22,647,                 FNV=df5374d4aca67677,
SHA256=14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317,
max endpoint=636,
top={(100,636),(256,636),(294,636),(338,636),(372,636),(384,636)}.   (1)
```

Three distinct certificate consumers are proved here.

1. A rank-eight exchange followed by mixed rank-eight/rank-nine repairs gives
   one 9,019-mask carrier. It closes 70 rows of `U`: every current row reached
   by the complete descending scan through endpoint 627, and 21 of the 23
   rows at endpoint 626.
2. Among the 192 singleton inactive-signature groups in `U`, 110 groups
   admit a one-deletion/one-addition deck common to their complete signature
   ideal. These are 110 **separate** 421-mask decks and close 1,219 disjoint
   rows.
3. The complete index-19 ideal has 36 rows. Deleting joint mask 19 and adding
   two exact responders gives one separate 422-mask deck common to all 36.

THM-4295 independently supplies common decks on the complete index-294 and
index-372 ideals. Those are imported only as two additional typed consequence
nodes; their masks are not merged with any consumer above. After exact overlap
accounting across the five nonredundant row sets, the aggregate residual is

```text
count=21,323,
FNV=09a0dfc4515d556b,
SHA256=a6fb3ea00f5017b2ad66cbccf75f89756aa94360652ff6a7dbf22c637a1c7656,
maximum endpoint=626,
top={(100,626),(256,626)}.                              (2)
```

The common decks are not merged with each other or with either carrier.

## 2. Rank-free deletion lemma and deletion depth

Fix an ordered endpoint pair `Q={q,r}` and a labelled nine-body
`B in binom(P,9)`. Call an arbitrary deletion mask `R subset P` active when

```text
mu(G_((P\R) union Q)) >= alpha.                         (3)
```

If `R intersect B` is empty, then `B subset P\R`; hence

```text
G_((P\R) union Q) subset G_(B union Q).                 (4)
```

Thus every active mask disjoint from `B` certifies the target body, regardless
of `|R|`. Moreover, activity is upward closed: if `R subset R' subset P`, then
`P\R' subset P\R`, so the safe set and its mass can only grow. Define

```text
delta_Q(B)=min{|R|:R subset P\B and R is active},       (5)
```

with value infinity if the set is empty. Since `|P\B|=21`, a rank-`k` scan
finds a disjoint active mask exactly when `delta_Q(B)<=k`, for every
`k<=21`: extend a minimum witness inside `P\B` and use upward closure.

This identifies what the inherited rank-eight scan was measuring. Rank eight
is an efficient representation, not part of the consequence lemma. The
tradeoff is exact:

```text
one rank-k mask is disjoint from binom(30-k,9) bodies;
one nine-body admits binom(21,k) rank-k masks.           (6)
```

Higher rank makes activity easier but each mask covers fewer bodies. The
rank-eight colex transforms and truncated failure atoms used earlier are not
silently reused at rank nine; the mixed programs retain all atoms through
popcount nine.

## 3. Exact endpoint-636 exchange

THM-4287's repaired 9,006-mask carrier had 101 failures: 64 at `(100,636)`
and 37 at `(256,636)`. Complete response quotients over all useful rank-eight
masks give

```text
local (100,636) minimum                 8,
local (256,636) minimum                 6,
common-active two-row minimum          15,
active carrier-union minimum           14.              (7)
```

The minimum fourteen witness is

```text
18468880 080e8281 22081017 08422a82 004cac40 19c04044 00c08ec0
10443016 01609124 10413209 01611640 00606449 0128d084 08806449,
FNV=8c648463d5cede1b.                                      (8)
```

An exact dual-gap certificate rejects depth 13. Its integer weights total
6,830 and every mask has load at most 528. Since `13*528-6830=34`, a
hypothetical 13-cover can spend at most 34 units on mask deficits plus overlap
charge. This leaves 103 classes and 58 maximals; the exact charged search
rejects depth 13 in 20 nodes.

Exactly fourteen old masks are strictly inactive on all nine inherited
boundary rows in the pre-THM-4287 22,682-row audit universe at endpoints 637
and 636:

```text
00003e1a 000132a3 00017464 00033388 000a16c2 000f8118 00142a1a
00154348 00184ba0 001aa260 00202c2b 002066a4 002b018a 0030c2a2,
FNV=a497f155f01aee9e.                                      (9)
```

Deleting `(9)` and appending `(8)` therefore keeps size 9,006 and gives the
exchange carrier

```text
|C_636|=9,006,                 FNV=8062ce6d5728da1f.     (10)
```

A detached literal-wall audit checks all 126 deleted activity cells strictly
negative, reconstructs every response in `(8)`, and scans
`9*binom(30,9)=128,764,350` labelled body instances with zero failures and
zero equalities. This is an exact minimum for the declared addition response
problem; it is not a claim that 9,006 is a globally minimum carrier.

## 4. The rank-eight obstruction and exact depth nine

The exchange carrier descends through endpoint 633. At endpoint 632 it leaves
six bodies at `(100,632)` and 66 at `(256,632)`. One of the latter is

```text
B*=1d106401, labels={0,10,13,14,20,24,26,27,28}.        (11)
```

All `binom(21,8)=203,490` rank-eight masks disjoint from `(11)` are inactive.
The unique best one is `222c100c`, whose exact literal activity margin remains
strictly negative:

```text
margin_ticks63=-182481808811274 on grid 23056825502430720. (12)
```

Every one of the thirteen one-bit extensions of `222c100c` inside
`P\(B* union 222c100c)` is rank nine and strictly active. Consequently

```text
delta_{(256,632)}(B*)=9.                                (13)
```

This is the first exact obstruction to continuing the inherited fixed-rank
carrier lane: no choice or optimization among rank-eight masks can repair
`B*`. The other 65 bodies have an exact ten-mask rank-eight cover.

For `(100,632)`, the six failures have exact rank-eight cover minimum three,
attained by

```text
2040c641, 00508325, 002a8641.                           (14)
```

For `(256,632)`, the complete mixed rank-eight/rank-nine quotient has 2,323
response classes and exact minimum three, certified by an integral dual of
value three and the rank-nine witnesses

```text
00619324, 201813a4, 21888126.                           (15)
```

Appending `(14)--(15)` gives

```text
|C_632|=9,012, rank8=9,009, rank9=3,
FNV=be9cf87002d6114a,                                   (16)
```

which closes every current row through endpoint 631. The two minima in
`(14)--(15)` are pairwise response minima; no combined six-mask optimality is
asserted.

## 5. Exact mixed-rank counterexample-guided descent

Repeating the exact failure-response calculation yields the following finite
descent. Each row's `minimum` is for the displayed boundary obligation set.

| boundary | failures before repair | exact minimum | selected additions | carrier after repair |
|---|---:|---:|---|---|
| `r=630`, pair `(100,630)` | 2 | 1 | `0010e125` (rank 8) | 9,013; FNV `932d1fae6fd8c384` |
| `r=629`, pair `(100,629)` | 28 | 4 | `002ac4c0` (rank 8), `3882a082,0041c325,08c28e40` (rank 9) | 9,017; FNV `07689a1534ce7327` |
| `r=628`, pair `(100,628)` | 4 | 2 | `02008327,0006e281` (rank 8) | 9,019; FNV `d7f0e06e154e78c2` |

At `r=629`, a complete 419-class quotient has exact rational dual value
`7/2`, hence lower bound four; the four displayed masks cover all 28 bodies.
At `r=628`, exhaustive testing of all 32,767 nonempty subsets of the fifteen
possible nonzero response patterns (eleven realized) rejects one addition and
accepts the displayed two.

The final carrier has

```text
rank8=9,013, rank9=6, total=9,019.                      (17)
```

It closes every scanned current row through endpoint 627. At endpoint 626:

```text
(100,626): 995 failures, FNV=f0189c7a4ab34b77,
(256,626):  10 failures, FNV=f394b2fdaca1f3f8,
other 21 rows: zero failures.                            (18)
```

The raw scan completes 72 rows, reports 1,005 failures only on the two rows
in `(18)`, and has ledger FNV `90b4bb035fce1ffa`. Therefore its typed carrier
success node has 70 rows. O0 and O3 transcripts and failure ledgers are
byte-identical.

## 6. Recursive singleton-signature ideals

For the ordered joint deck let

```text
I_E(p)={j in {0,...,420}:E_j is inactive at p}.         (19)
```

Inside `U`, exactly 3,738 rows in 192 groups have singleton signatures. For
each singleton index `j`, the exact surgery problem deletes `E_j`, enumerates
the bodies previously covered only by `E_j`, intersects activity across the
complete ideal, and asks for one appended responder.

Exactly 110 of the 192 complete groups admit such a responder. Their ideals
are disjoint and contain 1,219 rows in total:

```text
FNV=3706723fd3574334,
SHA256=fdef599768aba3cc1738e330f35c9b431486ab845023ef59050c546f22605ef1. (20)
```

For every successful group, retaining the other 420 joint masks and adding
its own witness gives a body-covering 421-mask deck common to every row of
that group. A detached literal audit checks 513,199 old-mask activity signs,
1,219 replacement signs, and zero equalities. The original deck has 3,512
private bodies across these indices; one exhaustive scan of all
`binom(30,9)=14,307,150` labelled bodies verifies that every private body is
recovered by its corresponding witness.

These are 110 separate common decks. The union in `(20)` is a proof-graph row
union, not a mask union and not one deck common to all 1,219 rows.

## 7. The index-19 ideal needs two replacements

The complete ideal

```text
H_19={p in U:I_E(p) subset {19}}
```

has

```text
|H_19|=36, FNV=5c8af37cf2f002e7,
SHA256=c0406e8d10138ad743e3e99076b2171809945518cc4f708d3ef46a7b0fb70777. (21)
```

Deleting `E_19=1804aa01` exposes eight bodies. Intersecting the full activity
universe across all 36 rows leaves 2,212,775 common-active rank-eight masks
and 41 response classes. There is no full responder, so one addition is
impossible. The two masks

```text
0000e649, 0184a205                                      (22)
```

have complementary responses `0x0c` and `0xf3`, proving exact minimum two.
The resulting 422-mask deck has FNV `dc50478119bc6c12`; its complete labelled
body scan has zero failures, and all 15,192 rebuilt activity cells are
strictly positive. The detached literal audit independently reproduces these
claims.

## 8. Typed proof graph and exact residual

Let `S`, `J`, and `K` denote respectively the 1,219 singleton rows, the 36
index-19 rows, and the 70 final-carrier rows. Their complete intersections are

```text
S intersect J = empty,
S intersect K = {(410,626),(506,626)},
J intersect K = {(338,628),(338,636)},
S intersect J intersect K = empty.                     (23)
```

Hence

```text
|S union J union K|=1,219+36+70-2-2=1,321,
FNV=69bcc5c5fc86ac8e,
SHA256=ed040ba5b26d40445e8917cc24c844d36c3a91f59bd08244f6a3bf0f3d38a44e. (24)
```

The primary set consumer and an independent C++ consumer agree byte-for-byte
on the 70-row carrier ledger, this 1,321-row union, and its 21,326-row
complement. The independent consumer
derives `H_19` from the full signature atlas, parses the raw 72-row descent
without assuming its endpoint, and verifies that those rows equal the complete
`r>=626` prefix of `U`.

Now let `A=H_294` and `B=H_372` be the separately certified ideals from
THM-4295. Its ten-row append-only carrier consequence set is a proper subset
of `K`, although the two physical carriers are different, and its `H_19` row
set equals `J`. The remaining intersections are

```text
|S intersect A|=18,  K intersect A={(294,636)}, J intersect A=empty,
|S intersect B|=52,  K intersect B={(372,636)}, J intersect B=empty,
A intersect B=empty,
all triple intersections among S,J,K,A,B are empty.       (25)
```

Consequently `A` adds exactly `(147,294),(147,590)` beyond `(24)`, while `B`
adds exactly `(372,619)`. Thus

```text
|S union J union K union A union B|=1,324,
FNV=f55ee025df29bb65,
SHA256=57bdcd932cc2c985e81e1b2d472a469cf4c2e11c9cb21dd4b1037c4ba098562a. (26)
```

Subtracting `(26)` from `(1)` gives `(2)`. A dedicated Python consumer and a
structurally independent C++ set consumer reproduce every pairwise overlap,
the three-row increment, `(26)`, and `(2)` byte-for-byte. Row-set inclusion of
the append-only carrier node in `K` does not identify their carrier masks;
overlap with `S` does not merge the index-294 or index-372 decks.

## 9. Audit correction and scope

MISTAKE-532 records one non-load-bearing margin diagnostic error discovered
during promotion. MISTAKE-533 corrects the endpoint-628 subset-space label:
the scan has 15 possible patterns and 11 realized classes. The discovery
cocycle printed a weakest row by comparing margin
numerators whose pair-dependent grids were not equal. That comparison is not
scale invariant: 38 of the 110 singleton annotations and the second
index-19 annotation change after comparing the exact fractions
`margin/(63*grid)`. The detached literal ledgers use the correct fractions.
All activity signs, response classes, exact covers, body scans, proof-graph
sets, and residual identities are unchanged. No noncanonical primary weakest
annotation is used above.

The theorem has the following strict limits.

- A carrier or common deck is a finite sufficient certificate, not physical
  entry into the fixed-pool reduction.
- Failure in `(18)` means that this carrier lacks a certificate; it does not
  say either physical row is dangerous.
- Higher-rank repair masks are legitimate by `(4)`, but no uniform bound on
  deletion depth, carrier size, or endpoint descent is inferred.
- The 110 singleton decks, the index-19/index-294/index-372 decks, the
  append-only carrier, and the mixed carrier remain separate consumers. Only
  their proved row consequences are unioned, with redundant row sets removed.
- No exclusive owner, semantic arrival, termination theorem, or arbitrary
  pair theorem is supplied.
- LRC(14) remains open.

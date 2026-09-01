---
id: THM-4303
title: "LRC(14) endpoint-595 twenty-five-row carrier closure"
status: >
  PROVED RELATIVE TO THM-4302 + FINITE-EXACT + INDEPENDENT RAW REPLAY PASS.
  THM-4302's fixed 9,019-mask C596 carrier closes exactly 25 of the 28 rows in
  its residual endpoint-595 top layer. The three carrier-failing rows are
  (96,595), (100,595), and (210,595), with respectively 116, 13, and 16
  failing labelled nine-bodies. The typed union has 1,658 rows and leaves
  20,989; its maximum endpoint remains 595 on precisely those three rows. No
  physical entry or LRC(14) follows.
source: root / alternating LRC14-JC2 session, 2026-09-01
depends_on:
  - THM-4302-lrc14-endpoint-596-response-minimum-four-and-size-preserving-exchange
related:
  - THM-4300-lrc14-size-preserving-response-staircase-and-index-297-ideal
  - THM-4296-lrc14-mixed-rank-deletion-depth-and-recursive-signature-closure
artifact_root: 05-knowledge/results/lrc14_endpoint595_twentyfive_closure_thm4303
artifact_manifest: 05-knowledge/results/lrc14_endpoint595_twentyfive_closure_thm4303/SHA256SUMS
artifact_manifest_sha256: d6d04ebbcc9c230e032c2927e043d5855f325fa284c93514e850e2c85d77c527
primary_scripts:
  - 04-computation/lrc14_endpoint595_twentyfive_closure_thm4303/endpoint595_primary_raw_audit.cpp
  - 04-computation/lrc14_endpoint595_twentyfive_closure_thm4303/endpoint595_independent_replay.cpp
  - 04-computation/lrc14_endpoint595_twentyfive_closure_thm4303/primary_independent_crosscheck.py
  - 04-computation/lrc14_endpoint595_twentyfive_closure_thm4303/typed_union_consumer.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. The primary path counts all
  nonjoint hits on all joint-exposed bodies. The self-contained independent
  path implements literal-wall geometry directly, recursively generates the
  complete body universe, and stops at the first nonjoint witness. All eight
  common exact fields agree on all 28 rows. The typed consumer derives rather
  than assumes the audited layer and checks the complete row partition.
---

# THM-4303 -- LRC(14) endpoint-595 twenty-five-row carrier closure

**PROVED RELATIVE TO THM-4302 + FINITE-EXACT + INDEPENDENT RAW REPLAY
PASS. LRC(14) REMAINS OPEN.**

## 1. Fixed inherited carrier and row universe

Retain THM-4302's labelled thirty-speed pool `P`, threshold `alpha=4/63`,
22,647-row residual universe `U`, exact activity predicate, and carrier
criterion.  Thus a carrier closes a row `p=(q,r)` when every labelled
nine-body `B subset P` is disjoint from an active carrier mask at `p`.  This is
a finite sufficient carrier certificate and is not physical entry.

This theorem does not modify THM-4302's carrier:

```text
C_596: size 9,019, rank8 8,961, rank9 58,
FNV=892fef44a9e6b37e, all 421 joint masks retained.         (1)
```

THM-4302's typed residual has

```text
count=21,014,
FNV=7da11cd038486887,
SHA256=2a3ee951deb5b7cfbb4b86aabd4058c8073aae713b42afebabc15e3159deb3b6.
                                                                    (2)
```

The primary and independent programs each derive the residual row set whose
identity is recorded in `(2)` from the inherited universe and THM-4300 union,
rather than consuming a hard-coded endpoint-595 list.  Its complete maximum
layer is

```text
L_595={p in (2):endpoint(p)=595},
|L_595|=28,
FNV=47981ce64825ef2a,
SHA256=c607dab04e4f6849a2226f518771e43b1301d4fc582b47bfa5752c4643c93702.
                                                                    (3)
```

## 2. Primary complete raw replay

For each row in `(3)`, the primary audit reconstructs `(1)`, determines its
exact active subfamily, and enumerates all

```text
binom(30,9)=14,307,150                                     (4)
```

labelled nine-bodies.  It first tests the active joint deck.  For every body
exposed by that deck it counts every disjoint active nonjoint carrier mask.
Thus the run performs

```text
28 * 14,307,150 = 400,600,200                              (5)
```

body-row tests.  It finds 303,583 joint-exposed bodies and 10,664,533
nonjoint hit incidences.  Exactly 25 rows have no failing body.  The other
three have:

| row | failing bodies | failure FNV |
|---:|---:|---:|
| `(96,595)` | 116 | `fedacdbff3f31981` |
| `(100,595)` | 13 | `3ac9ac8b4b9ad93f` |
| `(210,595)` | 16 | `a6a226f12c168d3a` |

The total failure count is 145; the ordered global `(q,r,body)` FNV is
`f3d7f95fc38e7b49`.  The complete 145-row failure ledger is frozen in
`results/endpoint595_failures.csv`.  The full 28-row activity, exposure,
hit-range, and failure ledger has FNV `e97970fe403f275c` and is frozen in
`results/endpoint595_pair_audit.csv`.

Consequently `(1)` closes precisely the 25-row set

```text
S_25=L_595\{(96,595),(100,595),(210,595)},
|S_25|=25,
FNV=1ad4f3c2ab6ea09d,
SHA256=6839c30eafed0d6d73dcd6a0eff6ed7b4751798ff88dce33bcf282397ae247d5.
                                                                    (6)
```

The three rows not closed by this fixed carrier have identity

```text
F_3={(96,595),(100,595),(210,595)},
FNV=9853e590efc73022,
SHA256=5f521627485bc829bc6ec2c9d25ce3fb899cfe8d9b04a8f96a225da7e9d84dde.
                                                                    (7)
```

## 3. Structurally independent replay

The independent program does not include or call the primary
`audit_pair617` implementation.  It instead:

1. reconstructs `(1)` directly from every frozen component ledger;
2. implements the literal-wall grid and activity predicate locally;
3. derives the 28 rows in `(3)` from the typed residual;
4. recursively generates all 14,307,150 labelled nine-bodies and sorts them
   into numeric order; and
5. for each joint-exposed body, stops on the first disjoint active nonjoint
   mask.

Although the last step deliberately differs from the primary all-hit count,
the two paths agree exactly at all 28 rows on the eight common values

```text
active count, active FNV, joint count, nonjoint count,
exposed count, exposed FNV, failure count, failure FNV.     (8)
```

In particular the independent path reproduces all active-family and exposed
body identities, the 25/3 split, and each failure count/FNV in the table.
The frozen cross-consumer checks `(8)` field by field rather than comparing
only the final consequence.

## 4. Typed row consequence

Let `T_4302` be THM-4302's 1,633-row typed union.  Union only the proved fixed
carrier row consequences `(6)`:

```text
T_4303=T_4302 union S_25.
```

The exact new identities are

```text
|T_4303|=1,658,
FNV=43317f1aee06e8bd,
SHA256=bfeb739dcad61dd42bdd9a8b295b6058f3ecee5cc30acd64c469e3b8132393c7,

|U\T_4303|=20,989,
FNV=b0fbaa28440a118f,
SHA256=bbf2dbba58a5492f6e5d136f14940c3bac8b3ddea604b19e3e2b926abb8bad00.
                                                                    (9)
```

The new residual still has maximum endpoint 595, but its complete top layer
is exactly the three-row set `(7)`.

The typed consumer independently re-derives THM-4302's partition, requires
that the primary pair ledger is exactly the derived 28-row maximum layer,
extracts its zero-failure rows, and checks the complete union/residual
partition and every FNV/SHA-256 identity above.

## 5. Scope firewall

- This theorem audits one fixed carrier, namely `(1)`.  It does not assert
  that the three rows in `(7)` cannot be closed by another carrier.
- The 145 bodies are failures of `(1)`, not physical danger witnesses.
- No response minimum, repair set, or inactive exchange is claimed at
  endpoint 595.
- Only the 25 proved row consequences are added to the typed union.  No
  physical objects or carriers are identified across consumers.
- No row below endpoint 595 is tested or claimed.
- No physical-entry construction, arbitrary-pair theorem, semantic-arrival
  theorem, or terminating descent is supplied.
- LRC(14) remains open.

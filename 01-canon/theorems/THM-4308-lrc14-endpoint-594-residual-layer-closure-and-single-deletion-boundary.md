---
id: THM-4308
title: "LRC(14) endpoint-594 residual-layer closure and exact single-deletion boundary"
status: >
  PROVED RELATIVE TO THM-4306/4307 + FINITE-EXACT + DIRECT FULL-UNIVERSE
  RAW AUDIT AND OPTIMIZATION-LEVEL REPRODUCIBILITY PASS. The THM-4307
  3,925-mask carrier closes all 25 rows in the complete endpoint-594 layer of
  the THM-4307 residual. The 22 rows still residual after THM-4306 raise the
  combined typed union to 2,036 and leave 20,611, maximum endpoint 593 on 16
  rows. On the new 25-row carrier layer alone,
  exactly 3,911 of the 3,925 single-mask deletions are safe and 14 are unsafe.
  No simultaneous-deletion, older-target deletion, physical-entry, or LRC(14)
  conclusion follows.
source: root + endpoint594_scout / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4306-lrc14-index-265-recursive-ideal-two-mask-replacement
  - THM-4307-lrc14-endpoint-595-support-threshold-residual-hypergraph-compression
related:
  - THM-4305-lrc14-endpoint-595-pair-tagged-response-exchange
  - THM-4302-lrc14-endpoint-596-response-minimum-four-and-size-preserving-exchange
artifact_root: 05-knowledge/results/lrc14_endpoint594_complete_layer_closure_thm4308
artifact_manifest: 05-knowledge/results/lrc14_endpoint594_complete_layer_closure_thm4308/SHA256SUMS
artifact_manifest_sha256: fdb26a62222d9507e2a4056bcad2d4269d308a923a627360eda539e902e36ba9
primary_scripts:
  - 04-computation/lrc14_endpoint594_complete_layer_closure_thm4308/endpoint594_carrier_audit.cpp
  - 04-computation/lrc14_endpoint594_complete_layer_closure_thm4308/endpoint594_all_singleton_quotient.cpp
  - 04-computation/lrc14_endpoint594_complete_layer_closure_thm4308/endpoint594_singleton_quotient.cpp
  - 04-computation/lrc14_endpoint594_complete_layer_closure_thm4308/typed_endpoint594_consumer.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. The raw verifier
  reconstructs the carrier from canonical component ledgers and scans all
  357,678,750 row-body cases at O3 and O2 with byte-identical ledgers. A
  complete all-mask singleton quotient separately partitions the same body
  universe by zero, one, or at least two joint witnesses and agrees at O3/O2.
  The typed consumer checks the inherited partition and new row consequence.
---

# THM-4308 -- LRC(14) endpoint-594 residual-layer closure and exact single-deletion boundary

**PROVED RELATIVE TO THM-4306/4307 + FINITE-EXACT + DIRECT FULL-UNIVERSE
RAW AUDIT AND OPTIMIZATION-LEVEL REPRODUCIBILITY PASS. LRC(14) REMAINS OPEN.**

## 1. Carrier target, typed target, and inherited carrier

Retain THM-4307's complete endpoint-594 residual layer

```text
K_594^(4307): size=25,
FNV=cce015c81f7121d9,
SHA256=920638d6fb23a8f6492d34cf50e7dc247c2eddfe7ba3f2088c59155e1a56167e.
                                                                    (1)
```

THM-4306 has already typed three of these rows, namely `(173,594)`,
`(381,594)`, and `(383,594)`, by its separate `H_265` deck. Thus its complete
maximum residual layer is the 22-row subset

```text
L_594^(4306): size=22,
FNV=8413f0d2282e4cd6,
SHA256=2a46ac360974ee95b5c468f1f76fb9ddd6b5165fa6e410dd3b6bad02ca93dd54.
                                                                    (2)
```

Retain THM-4307's fixed-pool carrier

```text
C_595^(4307): size=3,925, rank8=3,858, rank9=67,
FNV=6fbd0bffcf0ed78b, all 421 joint masks retained.       (3)
```

This theorem tests `(3)` directly on every row of `(1)`. It does not identify
the THM-4306 rebuilt deck with `C_595^(4307)`.

## 2. Direct raw closure of the full 25-row layer

For every row in `(1)`, the primary program reconstructs `(3)`, recomputes
the exact activity predicate, and scans every labelled nine-body.  Thus it
tests

```text
25*binom(30,9)=357,678,750                               (4)
```

row-body cases.  The complete result is

```text
joint-exposed bodies=46,178,
nonjoint hit incidences=1,757,668,
failures=0,
pair-ledger FNV=996ebeec37e58e98.                        (5)
```

The O3/`NDEBUG` and O2 builds give byte-identical transcripts, pair ledgers,
and empty failure ledgers.  Their common transcript SHA-256 is
`c3fea5ada25f87506c7658a500b31cedb184970a68e4b03f4f25ccb3afc79650`,
and their common pair-ledger SHA-256 is
`22a4dc5dfa98561893fdfbacb1c103867b936e9f9f457b221d0698051ca74086`.
Consequently `C_595^(4307)` closes the full 25-row layer `(1)`. Together with
THM-4307, the same carrier closes its inherited 391-row target and these 25
additional rows, a 416-row union.

## 3. Typed consequence

Only the 22 still-untyped rows `(2)` are unioned with THM-4306; the other
three carrier consequences were already present there. The typed consumer
audits this seam and the inherited partition, obtaining

```text
|T_4306 union L_594^(4306)|=2,036,
FNV=25b8760c2857fc90,
SHA256=0b3481bbca63e0d5e8bc255758e140158d2ebfced3c6e21b55682fe3f1206868,

|U \ (T_4306 union L_594^(4306))|=20,611,
FNV=5cadf1a8e72e8ed6,
SHA256=0192ed5967974b2bf8a41a7df300475fbff5b76f6ca59acbbd5c7294306f3219.
                                                                    (6)
```

The new residual maximum is 593 on exactly 16 rows:

```text
(96,593),  (100,593), (105,593), (147,593),
(186,593), (192,593), (206,593), (210,593),
(220,593), (256,593), (260,593), (294,593),
(332,593), (366,593), (384,593), (520,593),
FNV=5424c07fa724011f,
SHA256=a4bf3d1e9aff29be4fb07d2644912cc098dbd996bc527b4f22dee980f36a6ce6.
                                                                    (7)
```

## 4. Exact all-mask one-deletion boundary on the full layer

For a row-body obligation `o` in `(1)`, let `W_o` be its set of active
carrier witnesses in `(3)`. Since `(5)` proves `W_o` is nonempty, deleting a
single mask `m` fails on `o` exactly when

```text
W_o={m}.                                                   (8)
```

The complete quotient partitions all cases in `(4)` by their number of active
joint witnesses:

```text
zero joint witnesses=46,178,
one joint witness=332,823,
at least two joint witnesses=357,299,749.                 (9)
```

On the zero-joint branch it counts all nonjoint witnesses. On the one-joint
branch it tests whether any nonjoint witness remains. The at-least-two branch
cannot become empty after one deletion. This is a complete decision of `(8)`
for every mask in `C_595^(4307)`.

There are exactly 17 singleton obligations: 3 with a joint witness and 14
with a nonjoint witness. They protect exactly 14 distinct masks, with ordered
mask FNV `0141fea29da37882`; the complete ordered
`(q,r,body,witness)` FNV is `541dc881d5cf3d42`. Hence, on
`K_594^(4307)` alone,

```text
safe single-mask deletions=3,911,
unsafe single-mask deletions=14.                         (10)
```

The O3 and O2 quotient transcripts and both certificate ledgers agree
byte-for-byte. A separate zero-joint implementation reproduces the 14
nonjoint singleton obligations on 11 masks. As a hostile raw control, deleting
protected mask `00b0832c` causes exactly two failures, both at `(96,594)`, on
bodies `054e5001` and `0d4c5001`; the pair-tagged failure FNV is
`ad86682e1019beab`.

## 5. Scope firewall

- The carrier replay covers all 25 rows in THM-4307's endpoint-594 residual
  layer. Only the 22-row subset `(2)` is a new typed consequence over
  THM-4306; the three-row overlap is audited explicitly rather than counted
  twice.
- The exact `3,911/14` split concerns deleting one mask from
  `C_595^(4307)` on the 25-row layer only. It says nothing about simultaneous
  deletions or preservation of THM-4307's older 391-row target.
- Compiler-level O2/O3 agreement is reproducibility of one algorithm, not two
  structurally independent proofs.
- Carrier failures are proof obligations, not physical danger witnesses.
- The carrier acts only on the inherited labelled thirty-speed pool. No
  physical entry, arbitrary-pair theorem, terminating descent, or proof of
  LRC(14) follows.

The cheapest next exact task is the direct raw replay of the 16 rows in `(7)`.
If it fails, the next object is the complete pair-tagged response quotient of
the frozen failures; if it closes, the typed frontier should advance before
another deletion computation.

# THM-4314: endpoint-591 complete layer closure and deletion boundary

**FINITE-EXACT, relative to THM-4313. The statements below concern one fixed
3,925-mask carrier and the labelled 30-speed model. They do not give a
physical entry or prove LRC(14).**

## Complete endpoint-591 closure

Let `C` be THM-4313's final carrier:

```text
|C|=3,925, rank8=3,818, rank9=107,
FNV=a0d08a38c10bdab7, all 421 joint masks retained.
```

The inherited endpoint-591 frontier is the ordered 13-row set

```text
(72,591) (96,591) (100,591) (105,591) (153,591)
(192,591) (210,591) (256,591) (260,591) (294,591)
(366,591) (384,591) (520,591),
```

with ordered FNV `fc332c0697c671c7`. Independent O2 and O3 builds reconstruct
`C` from canonical component ledgers and scan every one of

```text
13*binom(30,9)=185,992,950
```

row-body cases. Their transcripts, pair ledgers, and empty failure ledgers are
byte-identical. The exact totals are

```text
joint-exposed bodies                 28,791
nonjoint hit incidences           1,031,512
failures                                   0
pair-ledger FNV          8191899a3e142a2c
```

The minimum nonjoint hit count in the independent capacity audit is two,
attained on rows `(96,591)`, `(105,591)`, and `(210,591)`. Thus the unchanged
THM-4313 carrier closes the complete endpoint-591 layer.

## Complete deletion classification through size two

For a row-body obligation, call a mask in `C` a witness when it is active on
that row and disjoint from the rank-nine body. A deletion set `D` preserves the
layer exactly when every obligation has a witness in `C - D`.

A second exhaustive O2/O3 program counts all original witnesses, retaining
every obligation of multiplicity at most two. Across the same 185,992,950
obligations it finds

```text
zero-witness obligations       0
one-witness obligations        2
two-witness obligations       26
distinct two-witness pairs    23
```

The two unique witnesses are the rank-eight protected-joint masks

```text
024085d0   for (96,591,051c6208)
20a09640   for (96,591,1d584001).
```

Therefore, for every `D` contained in `C` with `|D| <= 2`:

- a singleton deletion is unsafe exactly for `024085d0` or `20a09640`, so
  2 of 3,925 singletons are unsafe and 3,923 are safe;
- a two-mask deletion is unsafe exactly when it contains one of those two
  critical masks or equals one of the 22 additional pairs flagged in
  `endpoint591_distinct_twohit_witness_pairs.csv`.

There are

```text
C(3925,2)                         = 7,700,850 pairs;
pairs containing a critical mask =     7,847;
additional two-witness pairs     =        22;
unsafe pairs                     =     7,869;
safe pairs                       = 7,692,981.
```

This is complete because deleting at most two masks cannot exhaust an
obligation having at least three distinct original witnesses. No statement is
made for deletion sets of size three or more.

As a protected-joint corollary, retaining all 421 joint masks makes every
single nonjoint deletion safe. Exactly 12 nonjoint pairs are unsafe; these are
the 12 `NN` rows in the distinct-pair ledger. This corollary is not promoted to
arbitrary simultaneous deletion robustness.

## Typed successor

Consuming the 13 certified rows into THM-4313's exact partition gives

```text
typed union: 2,100 rows,
FNV=3b2d991da091a7df,
SHA256=1027325378caa6d4853112fcdb006796c32180b71f82a5a7db8addbec821f01c;

residual: 20,547 rows,
FNV=59ca49a11d140ec5,
SHA256=b53abd545ddd28e95088231f615229a3fbfb0812510f84dc4fdfeda38e76f0c3.
```

The next residual endpoint is 590 on 13 rows,
FNV `44aa8a793d162cf9`, SHA-256
`eefb84e9fec0bcd809830a3c5ed18b0bb2aea1a2f2cb8b4994fafeaf1f5d01ce`:

```text
(72,590) (96,590) (100,590) (105,590) (186,590)
(192,590) (208,590) (210,590) (256,590) (260,590)
(294,590) (416,590) (420,590).
```

No endpoint-590 body replay is included here.

## Scope firewall

- The closure and deletion statements quantify only over the fixed THM-4313
  carrier, the fixed 13 rows, and all rank-nine bodies on 30 labels.
- O2/O3 variants of the capacity, all-witness, and protected-joint scans are
  byte-identical. Normal and optimized typed/verifier paths are also frozen.
- The size-at-most-two deletion theorem is an exact consequence of the
  complete low-multiplicity census, not a heuristic sample.
- There is no globally minimum carrier theorem, endpoint-590 closure,
  arbitrary larger-deletion theorem, physical entry, terminating descent, or
  proof of LRC(14).

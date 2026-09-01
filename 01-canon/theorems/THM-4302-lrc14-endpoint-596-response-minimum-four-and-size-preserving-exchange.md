---
id: THM-4302
title: "LRC(14) endpoint-596 response minimum four and size-preserving exchange"
status: >
  PROVED RELATIVE TO THM-4300 + FINITE-EXACT + COMPLETE RESPONSE-DUAL/COVER
  AUDITS PASS. On the 24 THM-4300 carrier failures at (210,596), the exact
  active rank-eight/rank-nine response-cover minimum is four; rank eight alone
  has exact minimum six. Four rank-nine responses and deletion of 73 masks
  strictly inactive on the complete 363-row endpoint-at-least-596 prefix give
  a 9,019-mask carrier closing that prefix. The typed union has 1,633 rows and
  leaves 21,014, maximum endpoint 595. No physical entry or LRC(14) follows.
source: root / alternating LRC14-JC2 session, 2026-09-01
depends_on:
  - THM-4300-lrc14-size-preserving-response-staircase-and-index-297-ideal
related:
  - THM-4296-lrc14-mixed-rank-deletion-depth-and-recursive-signature-closure
  - THM-4254-fixed-ceiling-band-signed-endpoint-cocycle-cascade
artifact_root: 05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302
artifact_manifest: 05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302/SHA256SUMS
artifact_manifest_sha256: 7c522277ebfc489aad4375b306f4365cf5790bcb5e6c148df0ddc1bee214733b
primary_scripts:
  - 04-computation/lrc14_endpoint596_response_exchange_thm4302/endpoint596_response_exchange_audit.cpp
  - 04-computation/lrc14_endpoint596_response_exchange_thm4302/endpoint596_independent_audit.cpp
  - 04-computation/lrc14_endpoint596_response_exchange_thm4302/typed_union_consumer.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. The full-universe O0/O3
  transcripts agree byte-for-byte. An independent complement-generated
  response audit first reconstructs the 24-body universe by a raw scan of all
  labelled nine-bodies without consuming the failure CSV. The independent raw
  before/after replay and typed row-set consumer agree with the primary path.
---

# THM-4302 -- LRC(14) endpoint-596 response minimum four and size-preserving exchange

**PROVED RELATIVE TO THM-4300 + FINITE-EXACT + COMPLETE
RESPONSE-DUAL/COVER AUDITS PASS. LRC(14) REMAINS OPEN.**

## 1. Inherited universe and statement

Retain THM-4300's labelled thirty-speed pool `P`, threshold `alpha=4/63`,
safe-set notation

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14},
```

and 22,647-row residual universe `U`.  For a residual pair `p=(q,r)`, a mask
`R subset P` is active when

```text
mu(G_((P\R) union {q,r})) >= 4/63.                         (1)
```

A carrier closes `p` when every labelled nine-body `B subset P` is disjoint
from an active carrier mask.  This is a finite sufficient carrier
certificate, not a physical entry.

THM-4300 constructs the augmented mixed carrier

```text
C_+: size 9,088, rank8 9,034, rank9 54,
FNV=55e8588798885ae5.                                      (2)
```

It closes the complete endpoint-`>=597` prefix and eight of the nine rows at
endpoint 596.  At

```text
p_0=(210,596)                                               (3)
```

it has exactly 24 failing labelled nine-bodies.  In the maintained order they
are

```text
B00=0598b001  B01=0708f001  B02=07c81088  B03=07c83008
B04=07cc3000  B05=0fc83000  B06=1518b400  B07=1588b001
B08=158ab000  B09=15a8b000  B10=17883001  B11=1788b000
B12=17893000  B13=17c83000  B14=25007049  B15=25007409
B16=2500f009  B17=2500f088  B18=2502f001  B19=2504d088
B20=2520e088  B21=27801281  B22=27803009  B23=27803081,
FNV=3dbd5b39673070ff.                                       (4)
```

For any active rank-eight or rank-nine mask `R`, define its response on these
obligations by

```text
rho(R)={i in {0,...,23}:R intersect B_i=empty}.             (5)
```

The claims are:

1. the exact cover number of the complete realized mixed response family in
   `(5)` is four;
2. the exact cover number after restricting to rank eight is six;
3. four rank-nine additions and 73 deletions strictly inactive on the complete
   endpoint-`>=596` prefix give a 9,019-mask carrier closing that prefix; and
4. adjoining only this row consequence to THM-4300's typed union gives 1,633
   rows and leaves 21,014 residual rows, with maximum endpoint 595.

## 2. Complete response inventory

The primary audit enumerates the two complete mask universes at `p_0`:

| rank | all masks | active masks | nonempty responders | responder FNV |
|---:|---:|---:|---:|---:|
| 8 | 5,852,925 | 1,267,366 | 9,090 | `2dddbe3405491cdd` |
| 9 | 14,307,150 | 7,100,734 | 138,019 | `2202e93c739926df` |

Responder FNVs are taken in increasing numeric mask order.  The union of the
two ranks realizes exactly 718 distinct nonempty responses, of which 82 are
inclusion-maximal.  Because the bodies in `(4)` fail for `C_+`, the audit also
checks that no nonempty responder is already in `C_+`.

### Mixed lower and upper certificates

Give weight `1/2` to obligations

```text
{B01,B02,B06,B08,B15,B19,B21}                              (6)
```

and zero to the other seventeen.  On every one of the 718 realized responses,
the sum of weights is at most one.  Equivalently, after scaling by two the
complete load histogram is

```text
load 0:255 response types, load 1:355, load 2:108.          (7)
```

Thus `(6)` is a fractional packing of value `7/2`, so every integral response
cover has size at least `ceil(7/2)=4`.

The following four active rank-nine masks attain four.  The response column is
the 24-bit encoding of `(5)`.

| mask | response |
|---:|---:|
| `185f0200` | `00d3c000` |
| `28070a88` | `00002ec3` |
| `02710e02` | `000f4180` |
| `0010c574` | `00e0343c` |

Their bitwise response union is `00ffffff`, and their ordered mask FNV is
`dc0eebaebf688c65`.  Hence the mixed active rank-eight/rank-nine response-cover
minimum is exactly four.

### Rank-eight hostile control

Restrict now to the 220 response types realized by rank-eight masks.  Give
weight `1/2` to

```text
{B00,B01,B06,B08,B15,B19,B20,B22,B23},                    (8)
```

weight one to `B02`, and zero elsewhere.  Every rank-eight response has total
weight at most one.  Its scaled load histogram is

```text
load 0:78 response types, load 1:104, load 2:38.            (9)
```

The dual value is `9/2+1=11/2`, so every integral rank-eight cover has at least
six masks.  The six active rank-eight masks

```text
2010c125:0000303c  12980602:001f4000  183a0280:0041c000
20014a89:00002b70  082244c8:00003c91  08720c08:00a01c82    (10)
```

have response union `00ffffff`.  Therefore the rank-eight-only minimum is
exactly six.  This hostile control shows that rank nine is essential to the
four-mask optimum.

## 3. Complete-prefix inactive exchange

Let

```text
K_596={p in U:endpoint(p)>=596}.                            (11)
```

It has 363 rows.  The primary audit evaluates every one of the

```text
9,088 * 363 = 3,298,944                                    (12)
```

literal mask-row signs for `C_+`.  There are no equality cells.  Exactly 75
masks are inactive at every row of `K_596`; every one is nonjoint and rank
eight, and their increasing-order FNV is `fa143e58f59119f8`.

Let `D` be the lexicographically first 73 of those 75 masks, frozen in
`inputs/delete73.txt`, and let `A` be the four masks in the table above, frozen
in `inputs/additions4.txt`.  Then

```text
|D|=73, FNV(D)=9240b264ab65aa62,
|A|=4,  FNV(A)=dc0eebaebf688c65.                           (13)
```

**Inactive-exchange lemma.**  Suppose a carrier `C union A` closes every row
of a finite row set `K`, and every `d in D subset C` is inactive at every row
of `K`.  Then `(C\D) union A` closes every row of `K`.

Indeed, deleting `D` removes no active mask at any row.  The old active
subfamily and all added active witnesses are unchanged, so every prior body
witness remains.

THM-4300 gives closure of `C_+` on every row of `K_596` except possibly
`p_0`; its exact boundary ledger gives the complete failure set `(4)` at
`p_0`.  The four responses cover `(4)`, so `C_+ union A` closes `p_0`; adding
masks cannot destroy closure at the other 362 rows.  Apply the lemma with the
sign audit above.  The resulting carrier

```text
C_596=(C_+\D) union A,
|C_596|=9,019, rank8=8,961, rank9=58,
FNV=892fef44a9e6b37e                                      (14)
```

closes every row of `K_596`.  All 421 joint masks remain, and two additional
common-inactive rank-eight masks remain unused.

The bookkeeping in `(14)` starts from the augmented 9,088-mask carrier.  It is
not an add-four/delete-four modification of THM-4300's already exchanged
carrier.  Relative to that 9,019-mask carrier, the exact ledger relation is a
five-for-five move: restore `08044436`, add the four masks in `A`, and delete

```text
28101154, 28400c13, 28a002c4, 300086c4, 30418806.          (15)
```

Thus `(14)` genuinely retains the canonical carrier size without conflating
the two deletion ledgers.

## 4. Independent raw replay

The structurally independent audit does not receive THM-4300's failure CSV.
It reconstructs `C_+` from the inherited carrier components and all 76
repairs, then enumerates all

```text
binom(30,9)=14,307,150                                     (16)
```

labelled nine-bodies at `p_0`.  It derives exactly the 24 bodies in `(4)` and
the same FNV.  Its raw summaries are

| carrier | active | joint | nonjoint | exposed bodies | hit range | failures |
|---|---:|---:|---:|---:|---:|---:|
| `C_+` | 3,581 | 265 | 3,316 | 66,310 | -- | 24 |
| `C_+ union A` | 3,585 | 265 | 3,320 | 66,310 | `1..254` | 0 |
| `C_596` | 3,585 | 265 | 3,320 | 66,310 | `1..254` | 0 |

Starting from the raw-derived obligations, this path generates candidate
responders as rank-eight and rank-nine subsets of the 21-label complements,
rather than scanning the two full mask universes.  After deduplication and
the exact activity test it reproduces 9,090 and 138,019 responders, both FNVs,
all 718 response types, and both dual certificates.

It also independently reconstructs the 75-mask common-inactive pool and
checks all `73*363=26,499` deletion signs, with no equalities and sign-ledger
FNV `ee7857ffb11111b2`.  Thus the failure universe, response quotient, and
inactive transfer each have a second exact path.

## 5. Typed row consequence

The independent set consumer starts from THM-4300's maintained partition

```text
|U|=22,647,
|T_4300|=1,624,
|U\T_4300|=21,023.                                        (17)
```

It derives `K_596` directly from the endpoint predicate.  The exact identities
are

```text
|K_596|=363,
FNV=fbf4e6bd7a593649,
SHA256=6fa377dd16aaa12e62fb5e4c6ec36ba30701c9055a1d503b227cb76546f13960,

|T_4300 intersect K_596|=354=|K_597|,
FNV=33b0069ca994b786,
SHA256=e653857c9bfe1ef50e7724cfad05232b3695f88534ca289844a46d914ca52df5.
                                                                    (18)
```

Consequently `K_596` contributes precisely the nine endpoint-596 rows

```text
(96,596), (100,596), (192,596), (210,596), (256,596),
(260,596), (294,596), (306,596), (384,596).                (19)
```

The updated typed union and residual are

```text
|T_4302|=1,633,
FNV=b1c8ecf1dd4a71c5,
SHA256=28084fded429f407188e471183d5645d28eded621967e10e386261bbe52844c0,

|U\T_4302|=21,014,
FNV=7da11cd038486887,
SHA256=2a3ee951deb5b7cfbb4b86aabd4058c8073aae713b42afebabc15e3159deb3b6.
                                                                    (20)
```

The new residual has maximum endpoint 595.  Its complete top layer has 28 rows:

```text
(96,595),  (100,595), (147,595), (186,595), (192,595), (206,595),
(210,595), (220,595), (244,595), (256,595), (260,595), (265,595),
(294,595), (296,595), (306,595), (320,595), (332,595), (338,595),
(346,595), (366,595), (370,595), (372,595), (384,595), (416,595),
(420,595), (440,595), (512,595), (520,595),
FNV=47981ce64825ef2a,
SHA256=c607dab04e4f6849a2226f518771e43b1301d4fc582b47bfa5752c4643c93702.
                                                                    (21)
```

## 6. Scope firewall

- The exact minima four and six concern only the complete active rank-eight
  and rank-nine response families on the fixed 24-body carrier-failure set at
  `(210,596)`.
- The exchange proves a carrier certificate on the inherited labelled
  thirty-speed pool.  It is not a physical-entry construction.
- The 24 bodies are carrier failures before repair, not physical danger
  witnesses.
- Only proved row consequences are unioned.  No carriers, decks, or physical
  objects are identified across consumers.
- No arbitrary-pair theorem, semantic-arrival theorem, or terminating descent
  is supplied.
- LRC(14) remains open.

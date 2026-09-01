---
id: THM-4300
title: "LRC(14) size-preserving response staircase and index-297 ideal"
status: >
  PROVED RELATIVE TO THM-4287 + THM-4295 + THM-4296 + FINITE-EXACT +
  DETACHED LITERAL-WALL AUDITS PASS. A 69-response staircase followed by a
  69-for-69 strictly inactive exchange gives a 9,019-mask carrier closing
  every current row with endpoint at least 597. Separately, the complete
  42-row index-297 singleton ideal has exact common-deck replacement number
  two. The typed consequence union has 1,624 rows and leaves 21,023, maximum
  endpoint 596. No physical entry or LRC(14) follows.
source: root / lrc14-opus-20260831
depends_on:
  - THM-4287-repaired-carrier-endpoint-637-descent
  - THM-4295-lrc14-endpoint-636-minimum-fourteen-exchange-and-recursive-signature-ideals
  - THM-4296-lrc14-mixed-rank-deletion-depth-and-recursive-signature-closure
  - THM-4254-fixed-ceiling-band-signed-endpoint-cocycle-cascade
related:
  - THM-4286-signature-response-nonfactorization-and-two-deck-surgeries
  - THM-4283-endpoint-644-carrier-response-and-signature-fibre-surgery
artifact_root: 05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300
artifact_manifest: 05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300/SHA256SUMS
artifact_manifest_sha256: 16ee3b81d212bd5496de4433a733f43f2a91eca040f9e3f00f60f1aa08142cae
primary_scripts:
  - 04-computation/lrc14_size_preserving_response_staircase_thm4300/endpoint597_exchange_audit.cpp
  - 04-computation/lrc14_size_preserving_response_staircase_thm4300/endpoint597_full_raw_verify.cpp
  - 04-computation/lrc14_size_preserving_response_staircase_thm4300/raw_transcript_consumer.py
  - 04-computation/lrc14_size_preserving_response_staircase_thm4300/index297_detached_audit.cpp
  - 04-computation/lrc14_size_preserving_response_staircase_thm4300/index297_overlap_audit.py
  - 04-computation/lrc14_size_preserving_response_staircase_thm4300/typed_union_consumer.py
  - 04-computation/lrc14_size_preserving_response_staircase_thm4300/typed_union_independent.py
  - 04-computation/lrc14_size_preserving_response_staircase_thm4300/depth_profile.cpp
  - 04-computation/lrc14_size_preserving_response_staircase_thm4300/compressed_probe.cpp
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. O0/O3 literal sign and
  index-297 transcripts agree byte-for-byte. Independent full labelled-body
  replay, augmented/exchanged transcript comparison, rebuilt-deck body scan,
  signature-atlas reconstruction, and two typed set consumers agree.
---

# THM-4300 -- LRC(14) size-preserving response staircase and index-297 ideal

**PROVED RELATIVE TO THM-4287 + THM-4295 + THM-4296 + FINITE-EXACT +
DETACHED LITERAL-WALL AUDITS PASS. LRC(14) REMAINS OPEN.**

## 1. Inherited universe and statement

Retain THM-4254's labelled thirty-speed pool `P`, threshold `alpha=4/63`,
and safe-set notation

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14}.
```

Retain THM-4287's exact residual `U` and THM-4296's ordered joint deck `E`:

```text
|U|=22,647,                  FNV=df5374d4aca67677,
SHA256=14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317,
|E|=421,                     FNV=20d63dd42fe8150e.              (1)
```

For an endpoint pair `p=(q,r)`, call a deletion mask `R subset P` active when

```text
mu(G_((P\R) union {q,r})) >= 4/63.                         (2)
```

A carrier closes `p` when every labelled nine-body `B subset P` is disjoint
from an active mask in the carrier. This is a finite sufficient certificate,
not physical entry.

This theorem proves two distinct consumers. First, 69 exact responses extend
THM-4296's mixed carrier through the complete prefix
`K_597={p in U:endpoint(p)>=597}`; 69 original masks inactive on every prefix
row are then deleted, restoring size 9,019. Second, deleting joint-deck
coordinate 297 and adding two common-active masks gives a separate 422-mask
deck on the complete 42-row singleton ideal `H_297`.

After joining only these row consequences with THM-4296's maintained
1,324-row typed union, the residual is

```text
count=21,023,                  FNV=e93e8089e9dc58c0,
SHA256=ce215cb53a742e5e0d0d4f16e344687da9583ce5dc1417c7ea70399fb0bf70ba,
maximum endpoint=596.                                          (3)
```

## 2. The 69-response endpoint staircase

THM-4296's carrier `C_0` has

```text
|C_0|=9,019, rank8=9,013, rank9=6, FNV=d7f0e06e154e78c2. (4)
```

At endpoint 626 it fails on 995 bodies at `(100,626)` and ten at
`(256,626)`. The latter ten have exact mixed rank-eight/rank-nine response
minimum three, attained here by

```text
100c1485, 210c1096, 0a0a4216.                                (5)
```

A deterministic greedy response cover selects 28 masks for the former 995
bodies. It is irredundant but is **not** asserted minimum. Continuing the
literal failure/response loop gives:

| endpoint | failing blocks before repair | additions | size after layer |
|---:|---|---:|---:|
| 626 | `(100,626):995`, `(256,626):10` | 31 | 9,050 |
| 625 | `(210,625):1` | 1 | 9,051 |
| 624 | `(210,624):1` | 1 | 9,052 |
| 623 | `(210,623):34` | 4 | 9,056 |
| 622 | `(210,622):4`, `(260,622):1` | 2 | 9,058 |
| 621--620 | none | 0 | 9,058 |
| 619 | `(100,619):4`, `(210,619):603`, `(256,619):33` | 20 | 9,078 |
| 618 | none | 0 | 9,078 |
| 617 | `(294,617):1` | 1 | 9,079 |
| 616 | `(210,616)`, `(220,616)`, `(440,616)` | 3 | 9,082 |
| 615 | none | 0 | 9,082 |
| 614 | `(210,614):7` | 1 | 9,083 |
| 613 | none | 0 | 9,083 |
| 612 | `(260,612):1` | 1 | 9,084 |
| 611--603 | none | 0 | 9,084 |
| 602 | `(100,602):5`, `(210,602):2` | 3 | 9,087 |
| 601--600 | none | 0 | 9,087 |
| 599 | `(96,599):7` | 1 | 9,088 |
| 598--597 | none | 0 | 9,088 |

The frozen 76-mask ledger comprises THM-4296's seven inherited post-632
repairs followed by these 69 masks. It gives

```text
repair FNV=64ce5f9d1ec8c4c2,
augmented carrier C_+=9,088 masks, FNV=55e8588798885ae5,
rank8=9,034, rank9=54.                                      (6)
```

A detached scan checks all labelled bodies on all 354 rows of `K_597` and
finds zero failures. At the next layer, eight of nine rows close and
`(210,596)` has exactly 24 failures. No conclusion below 597 is made.

## 3. Strictly inactive exchange restores size 9,019

**Inactive-exchange lemma.** If a carrier `C` closes every row in `K`, and
each `d in D subset C` is inactive at every row of `K`, then `C\D` closes
every row of `K`.

At each row deleting `D` removes no active mask, so the active subfamily, and
therefore every body witness, is unchanged.

The exact sign audit evaluates

```text
9,019 * 354 = 3,192,726                                    (7)
```

literal cells of the original carrier on `K_597`, with zero equalities.
Exactly 76 original masks are inactive throughout; all are nonjoint rank
eight, with FNV `70cf014b73f82b0e`. Delete the lexicographically first 69,
whose FNV is `ce08bd348c1ac6c7`. The lemma gives

```text
C_597=C_+\D,
|C_597|=9,019, FNV=e0fcd06628e1aa37, rank8=8,965, rank9=54. (8)
```

All 421 joint masks remain. O0 and O3 sign transcripts are byte-identical;
their sign ledger FNV is `31605c3bc8dc5311`.

A second program directly replays `C_597`:

```text
354 * binom(30,9) = 5,064,731,100 labelled body-pair tests,
exposed bodies=1,045,532, nonjoint hit incidences=44,773,293,
failures=0, ledger FNV=f9fe3fc694b6aa98.                   (9)
```

Its 354 pair and 40 layer lines agree byte-for-byte with the augmented
transcript through layer 597. This independently checks body coverage rather
than merely invoking the exchange lemma.

## 4. The complete index-297 ideal

Let `I_E(p)` be the inactive coordinate set of the ordered joint deck and

```text
H_297={p in U:I_E(p)={297}}.                               (10)
```

Direct reconstruction from the complete signature atlas gives

```text
|H_297|=42, FNV=211843ee21a19170,
SHA256=ebe58d57e62a60be8325c7618146cb3e18cf9a20b51ee7270615290c109f06e6,
maximum endpoint=550 at (59,550).                          (11)
```

Deleting `E_297=08b204c0` exposes exactly four private bodies:

```text
B0=17087001, B1=2248d208, B2=27446008, B3=3548a008,
FNV=81e5ca5f045bb2cf, support union=374cf209.               (12)
```

Intersecting rank-eight activity over all 42 rows leaves

```text
2,208,823 of 5,852,925 common-active masks,
common-active FNV=dfc36b5aa486b0ea,
response FNV=b1665afb732ece83, ten response classes.         (13)
```

**Response-antichain lemma.** Replacing a finite realized response family by
its inclusion-maximal responses preserves the integral cover number. For a
nonnegative fractional dual, constraints from dominated responses are
redundant. Indeed, replace any dominated response in a cover by a realized
dominator; nonnegative weight on a subset is at most that on its superset.

Here the maximal antichain is exactly

```text
0x5={B0,B2}:    count 2, least mask 089284c0,
0xe={B1,B2,B3}: count 7, least mask 08330481.               (14)
```

No common-active response contains both `B0` and `B3`; hence dual weights
`y0=y3=1`, `y1=y2=0` have load at most one and value two. The two least masks
in `(14)` cover all obligations, so the exact minimum is two.

Deleting coordinate 297 and appending those masks gives

```text
|D_297|=422, FNV=60d75261322593ac.                          (15)
```

All `42*422=17,724` rebuilt mask-row margins are strictly positive, with no
equalities. A fresh scan of all 14,307,150 labelled bodies makes 405,172,243
checks and finds zero uncovered bodies. O0 and O3 reproduce the transcript,
ideal, and margin ledgers byte-for-byte.

## 5. Negative structural controls

Complete `2^21` complement transforms give the following deletion-depth
sequences for the four former `(100,628)` boundary bodies:

| body | at 628 | at 627 | at 626 |
|---|---:|---:|---:|
| `05346408` | 7 | 7 | 8 |
| `15306408` | 8 | 7 | 8 |
| `17581400` | 7 | 6 | 8 |
| `27d01008` | 7 | 6 | 8 |

Endpoint depth is therefore nonmonotone. All four bodies have exact depth
eight at 626 even though the inherited carrier fails them. The endpoint-628
responders `02008327` and `0006e281` have signs `(+,+,-)` at endpoints
`628,627,626`; the mechanism is responder deactivation, not depth above eight.

Colex compression is also false. At `(256,632)`, active mask `0000b28d` has
margin numerator `238551774045984`; the elementary left shift `15 -> 1`
gives `0000328f`, with numerator `-1296044026091136`. Both use the same
grid. Thus the active rank-eight family is not left-compressed.

## 6. Typed proof graph

Let `T` be the inherited 1,324-row typed union, `K=K_597`, and `H=H_297`.
Independent consumers derive `K` from `(1)` and `H` from the signature atlas:

```text
|T|=1,324, |K|=354, |H|=42, |T intersect K|=96,
T intersect H=K intersect H=T intersect K intersect H=empty. (16)
```

Hence `K` adds 258 rows and `H` adds 42. The new union is

```text
count=1,624, FNV=11414a33ab91fef6,
SHA256=ef9102553cd030f67ab1bdb7d6965c3efaf4b0d8aa85daa1092354c9703caf26. (17)
```

Subtracting `(17)` from `(1)` gives `(3)`. Its maximum layer is

```text
(96,596), (100,596), (192,596), (210,596), (256,596),
(260,596), (294,596), (306,596), (384,596).                (18)
```

The primary consumer reads the frozen ideal ledger. A structurally
independent consumer instead derives `H` from the signature bitset, derives
`K` by its endpoint predicate, checks the inherited partition, and reproduces
every count, overlap, FNV, SHA-256, and top row.

## 7. Correction and scope

MISTAKE-534 records a stale discovery diagnostic: a scout printed 437 by
taking the endpoint of terminal lexicographic row `(427,437)`, although its
ideal contained `(59,550)`. The detached audit computes the coordinatewise
maximum 550. No load-bearing object changed.

- The carrier and decks are finite sufficient certificates, not physical
  entry.
- The 24 failures at `(210,596)` are carrier failures, not physical danger.
- The 28-mask endpoint-626 `q=100` cover and later greedy blocks are not
  asserted globally minimum.
- Every deck and carrier remains a separate consumer; only row consequences
  are unioned.
- No uniform depth bound, terminating descent, arbitrary-pair theorem, or
  semantic-arrival theorem is supplied.
- LRC(14) remains open.

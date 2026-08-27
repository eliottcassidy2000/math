---
id: THM-4282
title: "Order-free inactive signatures, exact deck surgery, and carrier descent through endpoint 645"
status: >
  PROVED RELATIVE TO THM-4281 + FINITE-EXACT + DETACHED LITERAL-WALL
  AUDITS PASS. Three complementary rebuilt decks close 586+188 typed rows,
  and an explicit 45-mask augmentation of the inherited carrier closes the
  exact 90-row post-surgery endpoint band through 645. After exact overlap
  removal the theorem contributes 850 post-THM-4281 edges and leaves 23,373,
  with maximum endpoint 644 and exactly seven top rows. No physical entry or
  LRC(14) follows.
source: root/cross-frontier-bridge/2026-08-27
depends_on:
  - THM-4281-rectangle-common-joint-deck-endpoint-670-bridge
related:
  - THM-4276-six-atom-endpoint-671-augmentation-and-one-layer-descent
  - THM-4277-uniform-two-dimensional-outsider-rectangle-common-deck-closure
  - THM-4280-integral-three-channel-fat-contact-observer-and-sharp-five-jet-bound
artifact_root: 05-knowledge/results/lrc14_inactive_signature_deck_surgery_thm4282
artifact_manifest: 05-knowledge/results/lrc14_inactive_signature_deck_surgery_thm4282/SHA256SUMS
artifact_manifest_sha256: 6b7047973f8719f39ddefe90f6df42ae14b29f88be44f89d280a631ee6288182
primary_scripts:
  - 04-computation/lrc14_inactive_signature_deck_surgery_thm4282/final8951_joint_exposure_scan.cpp
  - 04-computation/lrc14_inactive_signature_deck_surgery_thm4282/endpoint650_large_response_greedy.cpp
  - 04-computation/lrc14_inactive_signature_deck_surgery_thm4282/proof_graph_consequence.py
primary_script_sha256:
  - 9986c6d2669737b7661fb2c69739945d4800e9ed0055ded38a25c9e3c9582564
  - b7432a397dbab965f15d8bfb4bb7cc92b823e9911323ec05ca035d6190693d6a
  - 9050fd5cca4261ef09ee4a4f6745c7f7c49e7ac9c3aa1471db7c327e37238ec0
primary_output_sha256:
  - 73dbe44577fd0850e7fec453308f726e9e50bc73f03b89d677cbcf2d161e21db
  - 890f12024337ea9980df53601a400e5795810b721c477088ea9f1bfa6f083706
  - 86502aa6c655151ff0e16b113b2da1716595044f2f628bd0df8094b541081d04
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT, subject to the frozen top-level SHA256SUMS. The index-367
  fibre, the five-mask index-520 surgery, and the delete-seven/add-eight
  index-256 surgery each have exhaustive body scans and detached literal
  controls. The final carrier is reconstructed from canonical dependencies;
  a source-pinned primary scan and a detached scanner independently derive
  and close the exact 90-row band. Normal, optimized, O0, O3, NDEBUG, and
  UBSan controls agree where stated. Exact Python and independent Ruby set
  replays agree on the 850-row union and 23,373-row residual.
---

# THM-4282 -- order-free inactive signatures, exact deck surgery, and carrier descent through endpoint 645

**PROVED RELATIVE TO THM-4281 + FINITE-EXACT + DETACHED LITERAL-WALL AUDITS
PASS. LRC(14) REMAINS OPEN.**

## 1. Statement

Retain THM-4281's labelled pool, threshold, frozen ordered 421-mask joint deck
`E`, and exact post-THM-4281 residual `U_4281`. Thus

```text
|U_4281|=24,223,                  FNV=80ec0687d8c7dba7,
max endpoint=663,                top={(256,663),(366,663),(520,663)}.
```

For a residual pair `p`, let

```text
I_E(p)={i in {0,...,420}: mask E_i is inactive at p}.          (1)
```

Complete exact evaluation of `(1)` gives an order-free inactive-signature
atlas. Three different local modifications of `E` then provide complementary
common-deck certificates:

1. deleting index 367 and appending one fibre-wide responder proves 26 rows;
2. deleting the five inactive indices of `(520,663)` and appending an exact
   minimum of six responders proves 560 post-THM-4281 rows;
3. deleting the seven inactive indices of `(256,663)` and appending eight
   explicit responders proves 188 rows.

The first two families are disjoint, so their union `C` has 586 rows. Starting
from the inherited 8,951-mask final carrier, an explicit nested 45-mask
augmentation closes every row of

```text
K={(q,r) in U_4281\C:r>=645},                         (2)
|K|=90,                       FNV=942995bee7469430.
```

The 188-row third family has overlap 9 with `C` and overlap 5 with `K`; `C`
and `K` are disjoint and the triple overlap is empty. Therefore the exact
proof-graph-new union has

```text
586+188+90-9-5=850                                    (3)
```

rows. Its removal leaves

```text
count=23,373,
FNV=c6ab0ae49ee32273,
SHA256=c3e5bf37887aa57af79cb166fce4a6e933e5daffc26dd8032fdfc52ce31240f3,
maximum endpoint=644,
top={(220,644),(256,644),(258,644),(294,644),
     (366,644),(416,644),(512,644)}.                   (4)
```

The three rebuilt decks and the augmented carrier are separate certificate
nodes. Equation `(3)` is a typed set union, not one globally common deck.

## 2. Exact signature atlas

Testing every mask of `E` at every row of `U_4281` gives

```text
rows                         24,223
tested activation cells     10,197,883
inactive incidences         864,784
equalities                  0
distinct signatures         17,604
signature-weight range      1..421
weight-one rows             5,038 over 215 masks
distinct nonzero columns    421
strict column implications  11
signature FNV               1f991a30499f0ae0
row FNV                     3652b5590b330704
```

The exact top signatures are

```text
(256,663): {75,107,139,374,394,405,417},
(366,663): {367},
(520,663): {57,107,222,275,345}.                       (5)
```

The 5,038 singleton rows force 215 masks into any subfamily that merely
classifies every row of `U_4281` as noncommon for the frozen deck. Those 215
masks cover 24,222 rows. The sole unclassified row is `(70,302)`, with
signature `{78,368}`. Hence there are exactly two minimum classifier bases,
each of size 216: the forced 215 plus index 78 or index 368.

This is a classifier theorem. A 216-mask classifier need not cover any
nine-body and is not a common-activity proof deck.

## 3. Signature-surgery lemma

For a body-covering deck `D` and labelled nine-body `B`, define the exact
response set

```text
R_D(B)={d in D:d intersect B=empty}.                   (6)
```

Let `J subset D` be deleted. Since `D` covers every body, the retained family
`D\J` misses `B` exactly when

```text
R_D(B) subset J.                                      (7)
```

Thus the exposed-obligation family is exactly

```text
X_D(J)={B in binom(P,9):R_D(B) subset J}.              (8)
```

At a fixed pair `p`, an active appended family `A_p` repairs the deletion if
and only if

```text
for every B in X_D(J), some a in A_p has a intersect B=empty. (9)
```

Necessity and sufficiency follow directly by partitioning the possible body
witnesses into retained and appended masks. Consequently the minimum number
of pair-local replacements is the exact set-cover number of `(8)` by the
complete active-response quotient at `p`. For a common certificate on several
pairs, activity must be intersected across those pairs before the same cover
test. This is the mechanism used in all three surgeries below.

## 4. The index-367 singleton fibre

The complete signature fibre of `(366,663)` consists of 26 rows, each with
signature exactly `{367}`. Mask 367 is `0x02188125`. Deleting it exposes
exactly three private bodies

```text
05646408,25065480,35a24080,             FNV=22bd212c1ffec6e2. (10)
```

At `(366,663)`, complete enumeration of all `C(30,8)=5,852,925` masks finds
61 active full responders, so the pair-local replacement minimum is one; the
least is `0x02108325`. Testing those 61 responders across the entire 26-row
fibre leaves exactly two fibre-wide responders,

```text
0x0a188803, 0x0a188a01.                                  (11)
```

Appending the least one in `(11)` after the retained 420 masks gives a
421-mask deck with FNV `87b42cf8a2069177`. It has zero failures on all
`C(30,9)=14,307,150` labelled bodies, and a detached literal-grid audit finds
all `26*421=10,946` activity cells strict. Hence all 26 rows are proved.

The pair-local mask `0x02108325` and the fibre-wide mask `0x0a188803` are
different witnesses for different quantified claims.

## 5. The index-520 five-deletion surgery

At `(520,663)`, delete the five indices in `(5)`, namely

```text
08c0a980,128c8900,009041ac,1aa28002,18868880.           (12)
```

The retained 416 masks expose exactly 53 bodies, with typed response FNV
`1e5d6dfe0b676151`. Complete enumeration at the target gives 2,879,147 active
masks and 7,124 response classes, of which 810 are inclusion-maximal. The
exact minimum replacement number is six. One least-representative witness is

```text
0090492c,018c2114,09a2a040,108c1112,1c81a100,38120016. (13)
```

The lower bound is the frozen integer dual: its 53 weights have numerator sum
57 over denominator 11, while every one of all 7,124 response classes has
weight at most `11/11`. Thus every cover has size at least
`ceil(57/11)=6`, and `(13)` attains it.

Appending `(13)` gives a 422-mask deck with FNV `813801c9bd1676ba` and zero
body failures. Independent endpoint-cocycle and detached literal scans find
the same 145,122 common rows in the post-THM-4277 universe, with zero
equalities. Relative to `U_4281`, it proves 560 new rows after two rows already
closed by the THM-4281 carrier are removed. This family is disjoint from the
26-row index-367 fibre, giving `|C|=586` and

```text
FNV(C)=9769742754535d84,
SHA256(C)=04aed4c107b3244a5e488266c46d1cae3bfffb20cddd831fc2d474c7b8c16a0e.
                                                               (14)
```

The rebuilt 422-mask deck loses 3,503 rows from the old 148,063-row common
family. It is complementary to, not a replacement for, the THM-4281 deck.

## 6. The index-256 explicit surgery

Delete the seven indices in `(5)` for `(256,663)`:

```text
75,107,139,374,394,405,417.                              (15)
```

They expose exactly 71 bodies. Retain the other 414 masks and append, in
order,

```text
0110a550,04871108,10241207,1042d088,
12848902,21141284,30249140,31202206.                     (16)
```

The resulting 422-mask deck has FNV `1c97b54ece61b351`. Exhaustive scanning
finds zero body failures after 405,775,694 short-circuit checks. On `U_4281`,
exact activity first restricts to the 395 rows whose old signature is
contained in `(15)`; testing the eight new masks gives exactly 188 common
rows,

```text
FNV=6588121dbec57bcb,
SHA256=5946ff653c51a74eba09a14430c494074e53b5aba87c3159bd17bafbe9e605d5.
                                                               (17)
```

A self-contained literal-wall engine independently reconstructs the body
decomposition and checks all `188*422=79,336` claimed activity cells, with
zero equalities. Equation `(16)` is an explicit cover witness only. No
minimum-eight claim is made.

## 7. Nested carrier descent

Start with the final THM-4281 carrier of 8,951 masks, FNV
`188f82ab9dd1695a`. Direct joint-exposure scans produce the following exact
descent.

| stage | exact repair action | carrier / FNV | resulting boundary |
|---|---|---|---|
| A | append `084a0a81` at the five failures of `(256,660)` | `8,952 / 9654d18926b1cf5b` | closed through 658 except `(256,657)` |
| B | append `0016580c` at the one failure of `(256,657)` | `8,953 / ad5249bee4b989b9` | closed through 651 |
| C | append exact local covers of sizes `1,4,3,1` at `(294,650),(366,650),(416,650),(512,650)` | `8,962 / 5b9e9a81c1582d06` | only `(256,650)` remains at endpoint 650, with 358 failures |
| D | append an explicit deterministic greedy cover of 31 masks | `8,993 / 05ff3d7ecbaa740c` | endpoint 650 and 649 closed; only `(256,648)` remains through 646 |
| E | append an exact minimum of three masks at the six failures of `(256,648)` | `8,996 / fd899660f14b311c` | every row in `(2)` closed |

The four stage-C covers are pair-local exact minima. Their nine-mask union is

```text
0002409f,
10459084,2009120b,20494209,3008408b,
0084c4c4,0200c325,12650202,
000c9681.                                                (18)
```

At stage D, the complete active response universe has 1,054,426 masks and
39,756 nonempty response candidates on the 358 exact failures. A disjoint
obligation packing gives lower bound 20, while the frozen 31 masks give an
explicit upper bound. No claim is made that 31 is minimum. At stage E the
exact minimum-three witness is

```text
00524a81,0128d00c,04410e81.                              (19)
```

The 45 masks from stages A--E are frozen in order. The source-pinned primary
scanner refuses any pair ledger other than the exact 90 rows in `(2)` and any
augmentation other than that exact 45-mask ledger. It enumerates every
labelled body independently at each pair and finds no exposed body missed by
the active nonjoint carrier. A detached literal-wall scanner reconstructs the
8,951-mask carrier from canonical dependencies, derives the same 90-row band
from the post-surgery residual, and reaches the same zero-failure conclusion.

No row with endpoint at most 644 is claimed by this carrier theorem.

## 8. Exact proof graph

Let `S` be the 188-row family in `(17)`. Exact set replay gives

```text
|C|=586,                    |S|=188,                 |K|=90,
|C intersect S|=9,          |C intersect K|=0,
|S intersect K|=5,          |C intersect S intersect K|=0.      (20)
```

The index-256 surgery therefore contributes 174 rows outside `C union K`.
The ledgers in `(20)` give

```text
C union S union K:
  count=850,
  FNV=8f595510210a5785,
  SHA256=7ad581bccd253e1778b972e8a303207da44534e6b995fa3ba15bd34b2801505b.
                                                               (21)
```

Removing `(21)` from `U_4281` gives exactly `(4)`. The main Python replay and
an independently written Ruby replay agree on every component, intersection,
union, and complement ledger.

## 9. Scope, controls, and non-consequences

- All body universes are the labelled `C(30,9)=14,307,150` universe, and all
  repair universes are the labelled `C(30,8)=5,852,925` universe.
- Every appended repair is rechecked for exact pair activity. A body-cover
  witness, classifier mask, active carrier mask, and common-deck mask are not
  interchangeable types.
- The minimum-one index-367 claim, minimum-six index-520 claim, four small
  stage-C minima, and minimum-three stage-E claim use complete exact response
  quotients. The index-256 eight-mask append and stage-D 31-mask append are
  explicit upper bounds only.
- Nested packet reproduction documents contain their original scratch paths;
  the top-level canonical reproduction routes them to this artifact root and
  supersedes those staging-path assignments.
- The conclusion removes exact fixed-pool residual edges. It supplies no
  exclusive owner, physical entry, semantic arrival, or proof of LRC(14).

`SHA256SUMS` is the authoritative raw-LF artifact ledger.

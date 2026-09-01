---
id: THM-4307
title: "LRC(14) endpoint-595 support-threshold residual-hypergraph compression"
status: >
  PROVED RELATIVE TO THM-4302/4303/4305 + FINITE-EXACT + COMPLETE RESPONSE-DUAL,
  RESIDUAL-HYPERGRAPH, AND FULL-PREFIX RAW AUDITS PASS. The exact mixed
  response minimum on the 145 endpoint-595 failures is ten. A support-threshold
  deletion followed by an exact minimum 37-mask retention repair within the
  fixed D_350 family gives a 3,925-mask carrier closing the complete 391-row
  endpoint-at-least-595 target.
  The typed union has 1,661 rows and leaves 20,986, maximum endpoint 594. No
  physical entry or LRC(14) follows.
source: root / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4302-lrc14-endpoint-596-response-minimum-four-and-size-preserving-exchange
  - THM-4303-lrc14-endpoint-595-twenty-five-row-carrier-closure
  - THM-4305-lrc14-endpoint-595-pair-tagged-response-exchange
related:
  - THM-4300-lrc14-size-preserving-response-staircase-and-index-297-ideal
  - THM-4296-lrc14-mixed-rank-deletion-depth-and-recursive-signature-closure
artifact_root: 05-knowledge/results/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307
artifact_manifest: 05-knowledge/results/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/SHA256SUMS
artifact_manifest_sha256: 639a895f023c29a07802d2a525b94d828940bafab58c393d922c581f4944e7b0
primary_scripts:
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/response_atlas.cpp
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/independent_complement_response_audit.cpp
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/solve_response_cover.py
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/joint_certificate_audit.py
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/activity_support_census.cpp
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/final_full_prefix_raw_audit.cpp
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/residual_retention_atlas.cpp
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/solve_residual_retention.py
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/k350_retention_certificate.cpp
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/independent_repaired_threshold_replay.cpp
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/verify_k350_retention.py
  - 04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307/typed_union_consumer.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. The response atlas scans
  both complete rank universes; an independent complement-generated path
  reproduces it. Exact integer dual/cover pairs certify both set-cover minima.
  The final 3,925-mask carrier is replayed directly on all 5,594,095,650
  row-body cases, while a second implementation reconstructs the support
  threshold, residual hypergraph, and repaired carrier.
---

# THM-4307 -- LRC(14) endpoint-595 support-threshold residual-hypergraph compression

**PROVED RELATIVE TO THM-4302/4303/4305 + FINITE-EXACT + COMPLETE
RESPONSE-DUAL, RESIDUAL-HYPERGRAPH, AND FULL-PREFIX RAW AUDITS PASS. LRC(14)
REMAINS OPEN.**

## 1. Fixed inherited target and pair-tagged obligations

Retain THM-4302's thirty-label pool `P`, threshold `alpha=4/63`, residual
universe `U`, activity predicate, and carrier criterion.  Its carrier is

```text
C_596: size 9,019, rank8 8,961, rank9 58,
FNV=892fef44a9e6b37e, all 421 joint masks retained.         (1)
```

Let `K` be the exact 391-row preservation target formed by the complete
363-row endpoint-`>=596` prefix together with the complete 28-row residual
endpoint-595 layer.  This is not the raw 394-row endpoint band: three
endpoint-595 rows were already typed by older proof nodes.

THM-4303 proves that `(1)` closes 25 of the 28 new rows.  At the remaining
three it has the complete pair-tagged failure ledger

| pair | failing labelled nine-bodies | failure FNV |
|---:|---:|---:|
| `(96,595)` | 116 | `fedacdbff3f31981` |
| `(100,595)` | 13 | `3ac9ac8b4b9ad93f` |
| `(210,595)` | 16 | `a6a226f12c168d3a` |

There are 145 obligations in all, with ordered `(q,r,body)` FNV
`f3d7f95fc38e7b49`.  For an active rank-eight or rank-nine mask `R`, define
its pair-tagged response by

```text
rho(R)={(p,B):R is active at p and R intersect B=empty}.  (2)
```

The pair coordinate in `(2)` is essential: a mask disjoint from a body but
inactive at that body's own row is not a response.

## 2. Complete response atlas

The primary audit scans the complete rank-eight and rank-nine universes:

| rank | all masks | active counts at the three pairs | nonempty responders | responder FNV |
|---:|---:|---:|---:|---:|
| 8 | 5,852,925 | `928827,1116371,1328767` | 51,271 | `cc926b13c719225d` |
| 9 | 14,307,150 | `6076461,6735949,7133375` | 679,004 | `8ff584ab870b72a1` |

Together they realize exactly 14,619 nonempty pair-tagged response types;
2,210 are realized in rank eight.  The complete type-ledger FNV is
`a5f5ebcdeb03ad34`.

A structurally independent audit starts from the 145 body complements rather
than scanning the two full universes.  It enumerates 29,506,050 rank-eight and
42,619,850 rank-nine complement-subset events and then recomputes all 730,275
responder signatures.  It reproduces both responder identities and every one
of the 14,619 response types.

## 3. Exact response-cover minimum ten

The separate mixed minima are

```text
(96,595):7,  (100,595):2,  (210,595):3.                  (3)
```

Each has an explicit cover and exact obligation-weight dual.  Solving the
joint pair-tagged problem saves two masks relative to the sum in `(3)`.

For the joint lower certificate, give the 145 obligations the frozen
nonnegative integer weights in `inputs/response_dual1667.csv`.  Every one of
the complete 14,619 realized responses has load at most

```text
D=1667,                                                     (4)
```

while the total weight is

```text
15030 > 9*1667=15003.                                     (5)
```

Thus nine masks cannot cover the obligations.  The following ten rank-nine
masks attain ten:

```text
02a08325 00b0832c 22c0a124 0228832c 00e12642
1284812c 10923282 10550a81 106a600a 14050236,
ordered FNV=6740cc137170afc5.                              (6)
```

Their response union is all 145 obligations.  Therefore the exact complete
mixed response-cover minimum is ten.  As a hostile control, the exact
rank-eight-only local minima are `17,4,4`; the coupled rank-eight optimizer
value 22 is retained only as a scout because no matching integer lower
certificate is promoted.

## 4. Strict inactivity and its exact obstruction

For a finite row set `S` and carrier `C`, put

```text
f_C(S)=#{m in C:m is inactive at every row of S}.          (7)
```

Each summand `1[S subset I_m]` is decreasing and supermodular, so `f_C` is
decreasing and supermodular; its complementary active-loss function is
increasing and submodular.  More exactly, the supermodular slack at `A,B`
counts masks active somewhere on both exclusive sides but inactive throughout
`A intersect B`.

On the target `K`, `(1)` has exactly two common-inactive masks:

```text
380086a0, 388088c0,
FNV=3b5ca775eedae38b.                                     (8)
```

All `9,019*391=3,526,429` signs are nonzero.  Within the restricted mechanism
that may delete only common-inactive masks, a same-size repair exists iff the
response minimum does not exceed `(7)`.  Equations `(6)` and `(8)` give
`10>2`, so strict-inactive size preservation is impossible.  Its exact
minimum carrier size is

```text
9,019+10-2=9,027.                                         (9)
```

Deleting `(8)` and adding `(6)` attains `(9)` and passes an independent raw
replay of all `28*binom(30,9)=400,600,200` endpoint-595 body-row cases.  This
is a mechanism-specific minimum, not a global carrier minimum.

## 5. Active-redundancy deletion complex

Let

```text
R=C_596 union A_10,                                       (10)
```

where `A_10` is `(6)`.  Thus `|R|=9,029`.  For every pair-tagged body
obligation `o=(p,B)` on `K`, let

```text
W_o={m in R:m is active at p and m intersect B=empty}.    (11)
```

### Deletion-complex lemma

A deletion set `D subset R` preserves carrier closure on `K` iff

```text
W_o is not a subset of D for every obligation o.          (12)
```

Thus safe deletion sets form a hereditary simplicial complex whose minimal
nonfaces are the inclusion-minimal members of the family of witness sets
`{W_o:o an obligation}`.  For a proposed deletion `D_0` and a candidate
retention set `H subset D_0`, restoring `H` repairs `R\D_0` iff `H` hits every
witness set belonging to an obligation exposed by `D_0`.  The minimum repair
is therefore an exact residual hypergraph transversal problem.

The proof is immediate from `(11)`: an obligation fails exactly when every
one of its witnesses was deleted.  This lemma explains why common inactivity
is only the zero-support boundary of a much larger deletion complex.

## 6. Support threshold 350 and exact 37-mask retention

For `m in C_596`, define its row-support size on `K` by

```text
sigma(m)=#{p in K:m is active at p}.                       (13)
```

The complete `9,019*391` sign census gives the support histogram frozen in
`results/activity_support_census.out`.  Let

```text
D_350={m in C_596:m is nonjoint and sigma(m)<=350}.        (14)
```

Then

```text
|D_350|=5,141,
FNV=03921cf597ee9863,                                     (15)
```

and all 421 joint masks are retained.  The direct full-prefix replay of
`R\D_350` finds exactly 84 failed pair-tagged bodies on 21 rows, with ordered
obligation FNV `3a4207b0cb910c10`.  This is the hostile boundary, not a
successful carrier.

Restrict `(11)` to these 84 obligations and the deleted masks in `(14).` The
complete residual quotient has

```text
70 responding deleted masks, FNV=b18cf54b71b4a093,
53 response types, type FNV=578ebd478ecee5a.              (16)
```

The 37 indices in `inputs/retention_packing37_indices.txt` form an integral
packing: every one of the 53 response types contains at most one packed
obligation.  Hence every retention cover has at least 37 masks.  The 37 masks
in `inputs/retain37.txt` cover all 84 obligations.  Their ordered mask FNV is
`531da886ae2e455a`, and their mask/response FNV is `061eaf6e2f6fbc8f`.
Consequently the exact residual retention minimum is 37.

Let

```text
D_final=D_350\H_37,
|D_final|=5,104,
FNV=ff4c932f9a7adac8.                                    (17)
```

By the deletion-complex lemma, the carrier

```text
C_595=R\D_final,
|C_595|=3,925,
rank8=3,858, rank9=67,
FNV=6fbd0bffcf0ed78b                                      (18)
```

closes `K`.  Moreover, 3,925 is the exact minimum among carriers of the fixed
form `R\(D_350\H)` with `H subset D_350`; this is not a global carrier lower
bound.

The primary raw program reconstructs `(18)` from the frozen component ledgers
and tests all

```text
391*binom(30,9)=5,594,095,650                             (19)
```

row-body cases simultaneously.  It finds zero failures, 1,418,344
joint-exposed bodies, 47,375,188 nonjoint hit incidences, and pair-ledger FNV
`bb28f7e567c4a4b0`.  A second implementation independently reconstructs the
support threshold, 84 residual witness edges, cover/packing certificates, and
final carrier identity.  It then tests all 5,588,707,348 body-row cases whose
witness sets can change under `D_final`, finds zero failures, and has row-ledger
FNV `52986f3355da0780`; every untested body retains its complete old witness
set by construction.

As a hostile stopping control, the analogous nonjoint support-`<=360`
deletion leaves 391 failed bodies on 36 rows before repair.  No claim is made
that threshold 350, its retention family, or `(18)` is globally optimal.

## 7. Typed row consequence

Carrier `(18)` closes the 363-row endpoint-`>=596` prefix and all 28 rows of
the residual endpoint-595 layer.  Unioning only these row consequences with
THM-4302 gives

```text
|T_4307|=1,661,
FNV=5bdd2ebf09e9404a,
SHA256=de00493a80ca88eb4ed802be00fce19967f0978508439bd07afc7393708a4b62,

|U\T_4307|=20,986,
FNV=606bf18913a49a14,
SHA256=67561f7f0c5c3a32155811e9978b42b2393c10ea8964387eb712cea6a6683f50.
                                                                    (20)
```

The residual maximum drops to endpoint 594 on exactly 25 rows, with FNV
`cce015c81f7121d9` and SHA-256
`920638d6fb23a8f6492d34cf50e7dc247c2eddfe7ba3f2088c59155e1a56167e`.

## 8. Scope firewall

- The minimum ten concerns the complete rank-eight/rank-nine response family
  on the fixed 145 carrier failures only.
- The minimum 37 and carrier-size minimum 3,925 are relative to the fixed
  support-threshold deletion family `(14)`, not all possible carriers.
- The 145 and 84 body ledgers are carrier failures, not physical danger
  witnesses.
- The carrier acts on the inherited labelled thirty-speed pool and supplies
  finite row certificates only.  It is not a physical-entry construction.
- Only proved row consequences are unioned; carriers and decks from distinct
  proof nodes are never identified.
- No arbitrary-pair theorem, semantic-arrival theorem, terminating descent,
  or proof of LRC(14) follows.

---
id: THM-4305
title: "LRC(14) endpoint-595 response minimum five and inactive-exchange obstruction"
status: >
  PROVED RELATIVE TO THM-4303 + FINITE-EXACT + INDEPENDENT RESPONSE/CARRIER
  AUDITS PASS. On the 145 pair-tagged failures of THM-4303's fixed C596
  carrier, the exact arbitrary-rank response-cover minimum is five. The exact
  rank-eight/rank-nine minimum is ten, and rank eight alone has minimum 22.
  Five masks close all three endpoint-595 rows and lower the typed residual
  maximum to 594, but the complete protected 391-row carrier audit finds only
  two strictly common-inactive masks. Thus the repaired carrier has size
  9,024 (or 9,022 after the safe two deletions); a size-9,019 continuation
  requires three redundancy-aware deletions. No physical entry or LRC(14)
  follows.
source: root / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4303-lrc14-endpoint-595-twenty-five-row-carrier-closure
  - THM-4296-lrc14-mixed-rank-deletion-depth-and-recursive-signature-closure
related:
  - THM-4302-lrc14-endpoint-596-response-minimum-four-and-size-preserving-exchange
artifact_root: 05-knowledge/results/lrc14_endpoint595_pair_tagged_response_exchange_thm4305
artifact_manifest: 05-knowledge/results/lrc14_endpoint595_pair_tagged_response_exchange_thm4305/SHA256SUMS
artifact_manifest_sha256: a7d4bb175ee4fcfd2b26feb14b5a9814143fe5dbdec6294ef1de577e09d8fd5a
primary_scripts:
  - 04-computation/lrc14_endpoint595_pair_tagged_response_exchange_thm4305/pair_tagged_atlas.cpp
  - 04-computation/lrc14_endpoint595_pair_tagged_response_exchange_thm4305/exact_cover.py
  - 04-computation/lrc14_endpoint595_pair_tagged_response_exchange_thm4305/rank_free_exact_cover.py
  - 04-computation/lrc14_endpoint595_pair_tagged_response_exchange_thm4305/rank_free_cover_independent.cpp
  - 04-computation/lrc14_endpoint595_pair_tagged_response_exchange_thm4305/protected391_activity.cpp
  - 04-computation/lrc14_endpoint595_pair_tagged_response_exchange_thm4305/typed_union_consumer.py
audit_scripts:
  - 04-computation/lrc14_endpoint595_pair_tagged_response_exchange_thm4305/rank_free_two_mask_intersection.cpp
  - 04-computation/lrc14_endpoint595_pair_tagged_response_exchange_thm4305/rank_free_two_mask_nextclosure.cpp
  - 04-computation/lrc14_endpoint595_pair_tagged_response_exchange_thm4305/protected391_independent.cpp
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. O0/O3 complete rank-eight/
  rank-nine atlases agree byte-for-byte. The arbitrary-rank lower bound has
  exact full-class CP-SAT replays under materially different configurations;
  an independent C++ literal-wall replay verifies the five-mask upper
  witness. NextClosure and intersection-closure programs independently give
  the same two-mask structural obstruction. Two independent carrier/grid
  implementations reproduce all protected-family signs and the exact
  common-inactive pair.
---

# THM-4305 -- LRC(14) endpoint-595 response minimum five and inactive-exchange obstruction

**PROVED RELATIVE TO THM-4303 + FINITE-EXACT + INDEPENDENT
RESPONSE/CARRIER AUDITS PASS. LRC(14) REMAINS OPEN.**

## 1. Inherited pair-tagged response problem

Retain THM-4296's rank-free deletion lemma and THM-4303's labelled
thirty-speed pool `P`, threshold `alpha=4/63`, literal pair grids, and fixed
carrier

```text
C596: size 9,019, rank8 8,961, rank9 58,
FNV=892fef44a9e6b37e.                                      (1)
```

For a pair `p=(q,r)` and mask `R subset P`, say that `R` is active at `p`
when

```text
mu(G_((P\R) union {q,r})) >= 4/63.                         (2)
```

THM-4303 proves that `(1)` closes the complete endpoint-`>=596` prefix and
25 of the 28 residual endpoint-595 rows. Its only failing rows and labelled
nine-body counts are

```text
(96,595):116,  (100,595):13,  (210,595):16.                (3)
```

Their 145 ordered `(q,r,B)` obligations have FNV
`f3d7f95fc38e7b49`. For any mask `R`, its pair-tagged response is

```text
rho(R)={(q,r,B) in (3):R is active at (q,r) and R intersect B=empty}.
                                                                    (4)
```

The first claim is that the exact cover number of the complete response
family `(4)`, with no rank restriction, is five. The second claim records the
sharper tariffs inside ranks eight and nine. The third claim determines the
strictly common-inactive deletion capacity on the already protected rows.
All are fixed-pool carrier statements, not physical entry.

## 2. Complete rank-eight/rank-nine response atlas

The complete atlas scans every mask in both universes, applies each pair's
own literal grid, and keeps the entire 145-bit tagged response. There are no
equality cells. The inventory is:

| rank | masks | global responders | responder FNV | response types |
|---:|---:|---:|---:|---:|
| 8 | 5,852,925 | 51,271 | `cc926b13c719225d` | 2,210 |
| 9 | 14,307,150 | 679,004 | `8ff584ab870b72a1` | -- |
| mixed | -- | -- | -- | 14,619 |

The pairwise responder counts/FNVs are

| rank | `(96,595)` | `(100,595)` | `(210,595)` |
|---:|---:|---:|---:|
| 8 | 37,261 / `acfe62887fe8b58b` | 7,661 / `9dc01c8a42f7bca8` | 8,793 / `678e88f84149a31b` |
| 9 | 501,371 / `a346172c10dd8be4` | 106,594 / `ba0cd6c95ddc561d` | 130,336 / `c08228f4a7590666` |

No responder is already in `(1)`, and neither rank contains a responder for
an entire family in `(3)`.

### Mixed rank-eight/rank-nine minimum ten

The ten rank-nine masks

```text
022083a5 00b0832c 22c0a124 0228832c 10e12602
1284812c 10923282 10550a81 106a600a 14050236              (5)
```

cover all 145 obligations. Their ordered FNV is `76e1eaf083842155`.

For the lower bound, the frozen integral dual assigns nonnegative weights
with denominator 10,000 to 50 obligations. Its numerator sum is 90,128, while
the load of every one of the 14,619 realized response types is at most 9,999.
Thus every cover has size `k` satisfying

```text
90,128 <= 10,000 k,
```

so `k>=10`. Together with `(5)`, the mixed rank-eight/rank-nine minimum is
exactly ten.

### Rank-eight-only minimum twenty-two

The frozen 22-mask rank-eight cover has ordered FNV `a9732e5fd06361c8`.
After removing dominated responses, 430 inclusion-maximal rank-eight types
remain. A deterministic one-worker CP-SAT feasibility audit of their complete
145-row incidence matrix, with `sum x<=21`, returns `INFEASIBLE`; the frozen
run uses OR-Tools 9.15.6755, subsolver `max_lp`, seed zero, 496 branches, and
zero conflicts. Hence the rank-eight-only minimum is exactly 22.

## 3. Arbitrary-rank exact minimum five

Rank-free activity must retain every literal-wall failure class
(MISTAKE-535). The three full geometries have respectively

```text
pair              grid             cells  distinct failure classes
(96,595)     36482318832960          8515          2453
(100,595)    91205797082400          8523          2519
(210,595)    18241159416480          8743          2543.   (6)
```

The following five masks attain a complete cover. All five are active on all
three rows; ticks are the pair-specific numerators
`63*mass-4*grid` and are not compared across grids.

| mask | rank | tagged counts `96/100/210` | ticks `96/100/210` |
|---:|---:|---:|---:|
| `00612a76` | 11 | `42/1/5` | `1234151772840 / 13335326729046 / 5266211229972` |
| `00a183f6` | 12 | `67/0/3` | `460514766222 / 38150010965556 / 10426181645376` |
| `024d8b32` | 12 | `17/0/3` | `2180791687914 / 33821374934088 / 5163277742154` |
| `0280a1ae` | 10 | `44/1/0` | `1015644756798 / 36013928708016 / 4729826045196` |
| `10110bf6` | 12 | `39/13/9` | `6001475952780 / 31100972723226 / 8198256331188` |

Their response union is `116/13/16`; there are 244 hit incidences. The ordered
mask FNV is `ebea2511eb7fa46f`, and the frozen pair-tagged response ledger FNV
is `8e13cc00c5bca917`.

For the lower bound, the exact finite model gives each of four mask slots 30
Boolean support variables. For every full failure class `F`, a Boolean is
constrained to be exactly the conjunction `F subset R`. Its weighted sum is
the exact literal mass in `(2)`, and the activity Boolean is constrained in
both directions by `63*mass >= 4*grid`. For every obligation, disjointness is
the exact conjunction of the nine negated support variables, and a hit is the
conjunction of activity and disjointness. Every one of the 145 obligations is
required to have a hit. Numeric ordering of the four masks is the only
symmetry breaker and loses no solution.

The unrestricted four-mask model is infeasible. The frozen one-worker hostile
run disables symmetry and still returns `INFEASIBLE`; a separate pure-
feasibility configuration also returns `INFEASIBLE`. Positive controls find a
four-mask cover of the 116-body `(96,595)` family and a four-mask joint cover
of the `(100,595)` and `(210,595)` families, so the contradiction is genuinely
cross-family rather than a vacuous activity encoding. Consequently the exact
arbitrary-rank response-cover minimum is five.

As a structural lower-rank control, two independent formal-concept programs
also prove that two arbitrary-rank masks cannot cover all 145 obligations.
The full no-cooccurrence relation has 548 concepts. The maximal whole-family
complements at `(96,595)` and `(210,595)` have negative ticks
`-57328586712180` and `-10387937227122`; monotonicity forces both masks in any
putative pair to be active on both rows, but none of the 548 maximal concepts
has that signature.

## 4. Carrier repair and exact inactive-exchange capacity

Let `A5` be the five masks in the table above. Since additions cannot destroy
an existing response, and `A5` covers every failure in `(3)`,

```text
C595+=(C596 union A5),
|C595+|=9,024,
ranks 8/9/10/11/12 = 8,961/58/1/1/3.                       (7)
```

It closes the complete 391-row family

```text
K*={rows of the inherited universe with endpoint >=596}
   union {the complete 28-row residual endpoint-595 layer}.             (8)
```

The set `(8)` consists of 363+28 rows and has canonical pair FNV
`c732a1532c12b9f6`.

The complete audit of all

```text
9,019 * 391 = 3,526,429                                   (9)
```

carrier-row signs has no equality cells. Exactly two masks in `(1)` are
strictly inactive on every row of `(8)`:

```text
380086a0, 388088c0,                                       (10)
```

both nonjoint rank eight, ordered FNV `3b5ca775eedae38b`. Every other carrier
mask is active somewhere in `(8)`. Thus deleting `(10)` from `(7)` preserves
all 391 closures and gives a 9,022-mask carrier.

This capacity is exact only for the simple common-inactive exchange. Since
Section 3 proves that every additive repair of the inherited 145 obligations
needs at least five masks, the two-mask slack in `(10)` cannot restore the
9,019 size. Any size-9,019 continuation through endpoint 595 must justify at
least three deletions that are active somewhere in `(8)` by a finer
row-dependent redundancy or replacement argument.

## 5. Typed row consequence

THM-4303's typed union has 1,658 rows and leaves 20,989 residual rows, whose
maximum endpoint is 595 on precisely `(3)`. Unioning the three newly proved
row consequences gives

```text
|T4305|=1,661,
FNV=5bdd2ebf09e9404a,
SHA256=de00493a80ca88eb4ed802be00fce19967f0978508439bd07afc7393708a4b62,

|U\T4305|=20,986,
FNV=606bf18913a49a14,
SHA256=67561f7f0c5c3a32155811e9978b42b2393c10ea8964387eb712cea6a6683f50.
                                                                    (11)
```

The residual maximum drops to 594. Its complete 25-row top layer has FNV
`cce015c81f7121d9` and SHA-256
`920638d6fb23a8f6492d34cf50e7dc247c2eddfe7ba3f2088c59155e1a56167e`.

## 6. Independent audits and validity boundary

- O0 and O3 complete rank-eight/rank-nine atlas runs and their 2,023,282-byte
  response CSVs agree byte-for-byte after normalizing the output path.
- A direct C++ full-wall replay, independent of the CP-SAT geometry code,
  reproduces every tick, all 145 hits, both five-mask FNVs, and the exact rank
  pattern in `(7)`.
- The arbitrary-rank four-mask infeasibility survives one-worker/no-symmetry
  and pure-feasibility configurations. The positive subfamily controls are
  frozen beside those outputs.
- NextClosure and iterative row-intersection enumerations independently
  reproduce the two-mask formal-concept obstruction.
- The primary and structurally independent protected-family programs
  reconstruct `(1)` and `(8)` separately and reproduce the Boolean sign FNV
  `855c6f82deef952f`, exact-margin FNV `eee2c2c38d5b2e55`, row-activity FNV
  `61549efd32b7ed3b`, zero equalities, and precisely `(10)`.
- The typed consumer derives the complete partitions in `(11)` from the
  inherited universe and the independently replayed five-mask response; it
  does not infer rows below endpoint 595.

## 7. Scope firewall

- The exact minima are relative to the frozen 145 pair-tagged failures and
  the literal activity predicate `(2)`. They are not minima for all possible
  carriers or physical configurations.
- `C595+` is a sufficient fixed-pool carrier certificate. It is not a
  physical counterexample entry, semantic-arrival map, or arbitrary-pair
  theorem.
- The common-inactive result concerns exactly `(8)`. It gives no deletion
  permission at lower rows, and it does not prove that active masks cannot be
  removed after a response-redundancy audit.
- The three additional deletions needed to recover size 9,019 are open.
- The typed union records only proved row consequences. Finite closure,
  physical entry, and LRC(14) remain open.

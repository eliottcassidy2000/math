---
id: THM-2896
title: "Seven-body adaptive six-cover hitting-gate atlas"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.  Every one
  of the 3432 seven-body roots has a globally sealed adaptive six-cover
  hitting gate of least size at most 25.  The unique K=25 root is
  (1,8,10,11,12,13,14); all others have K<=23.  This is a hitting reduction
  only, not an apex-carrier closure or LRC(14).
source: codex/lrc-rank-impossible-overlap-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
related:
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2892-eight-body-five-slot-heavy-triangle-closure
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
---

# THM-2896 -- seven-body adaptive six-cover hitting-gate atlas

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.**

## 1. Statement

For every

```text
E in C({1,...,14},7),
```

let `G_E` be its gap-`1/14` good set, put `h_E=|G_E|`, and define

```text
c_E(w)=|G_E intersect D_w|,                  w>=15.       (1)
```

Rank the external speeds by decreasing `c_E(w)`, breaking ties by increasing
speed, and write the resulting values as

```text
q_1(E)>=q_2(E)>=... .                                    (2)
```

Then there is a least integer `K(E)<=25` such that

```text
q_(K+1)(E)+...+q_(K+6)(E)<h_E.                           (3)
```

Consequently every six-speed set whose danger combs cover `G_E` meets the
first `K(E)` ranked external speeds.

The maximum is sharp inside this universe:

```text
K(1,8,10,11,12,13,14)=25,                               (4)
```

and this is the unique root with `K=25`.  Every other root has `K<=23`;
no root has `K=24`.

## 2. Hitting proof

Suppose a six-set `Q` avoids the first `K` ranked speeds.  Subadditivity and
the definition of the ranking give

```text
|G_E intersect union_(w in Q) D_w|
  <=sum_(w in Q)c_E(w)
  <=q_(K+1)+...+q_(K+6)
  <h_E.                                                  (5)
```

Thus `Q` cannot cover `G_E`, proving the hitting conclusion.  The strict
inequality in `(3)` is load-bearing; equality is never promoted.

## 3. Global top-forty seal

The verifier first scans every speed `15<=w<=1600`.  Let `q_40^0` be the
fortieth value in that finite scan, and let `r_E` be the number of interval
components of `G_E`.  On all `3432` roots,

```text
q_40^0>h_E/7.
```

THM-735(ii) gives the strict discrepancy estimate

```text
c_E(w)<h_E/7+(99/70)r_E/(7w).                           (6)
```

Put

```text
T_E=(99/70)r_E/[7(q_40^0-h_E/7)],
W_E=max(1601,ceil(T_E)).                                (7)
```

The verifier scans the missing finite head through `W_E-1`.  For every
unscanned `w>=W_E`, `(6)` and `(7)` give

```text
c_E(w)<h_E/7+(99/70)r_E/(7W_E)<=q_40^0.                 (8)
```

Hence the retained top forty are global.  Equality in the majorant at the
last step of `(8)` is safe because the actual discrepancy bound is strict.

The largest exact threshold is

```text
327134808/148661
```

at `E=(1,3,8,10,12,13,14)`, so the largest tail start is `2201`.  The
largest speed retained anywhere in a top-forty ledger is `574`.

The original top-thirty proposal fails honestly: for the root in `(4)`, no
`K<=24` can be read from thirty ranks.  Widening to forty exposes the strict
rank-`26` through rank-`31` margin

```text
5703505/4933292364>0,
```

and no assertion is relaxed.

## 4. Exact atlas

The global least-`K` distribution is

```text
K:count =
2:1, 3:3, 4:7, 5:33, 6:60, 7:140, 8:207, 9:285,
10:375, 11:436, 12:429, 13:385, 14:333, 15:245,
16:188, 17:112, 18:76, 19:55, 20:30, 21:20,
22:9, 23:2, 25:1.                                      (9)
```

By stratum:

```text
E subset {1,...,12}:                       792 roots, max K=21;
exactly one of {13,14}:                   1848 roots, max K=23;
both 13 and 14:                            792 roots, max K=25. (10)
```

The fixed `K=12` gate passes only

```text
485+1083+408=1976
```

roots in the three strata and fails on `1456`; adaptivity is therefore
structural.  The smallest positive adaptive margin is

```text
1669/6524317800
```

at `E=(2,3,8,9,11,12,13)`, with least `K=10`.

The computation performs four distinct vector/scalar controls on every
root, for `13728` controls total.  It locks the complete top-forty ledger,
all stratum distributions and extrema, and the three stratum digests.

## 5. Verification

```text
04-computation/lrc14_j6_all_root_adaptive_gate_atlas_codex_20260729.py
SHA-256 fc36f26d4c8da5b005465696b954eec700c080376eef9ee5ba74a7111def99d7

05-knowledge/results/lrc14_j6_all_root_adaptive_gate_atlas_codex_20260729.out
SHA-256 3081b93a870faacb31d205e43f7ca87872d7a9f196f4774d8740ced6a314d80b
```

Ordinary and optimized replays are byte-identical.  The complete ledger
digest is

```text
a76639c23fc613e664c2f4f35a492c2658b611ec8ec31388728d2883d11e4517.
```

The proof, tail boundary, least-`K` search, and scope were independently
audited.

## 6. Scope

This theorem gives a finite first-apex gate for every root in the exact
seven-body/six-slot rung left by THM-2892.  It does not prove that any
selected apex carrier is terminal.  THM-2893's ranked suffix refinement and
the scalar-bootstrap battery reduce that downstream workload.  THM-2897
identifies the next cheap terminal test as the suffix singleton cap plus
twice the suffix pair-union cap, but does not assert that it passes.  Seed
branches still require nonscalar cap/flag proofs.  Thus this theorem does
not prove the seven-body rung or LRC(14).

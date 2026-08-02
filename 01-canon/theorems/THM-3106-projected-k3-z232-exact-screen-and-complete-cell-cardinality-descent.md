---
id: THM-3106
title: "Projected k3 z232 exact screen and complete-cell cardinality descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  All three rows in
  the pinned projected k=3, z1=232 necessary layer are empty.  The projected
  ledger is 374322 and its cap is z1<=231.  The nine z1=231 rows remain
  occupied.  No LRC(14) claim is made.
source: codex-thm3094-hostile-audit-2026-08-02
audit: >
  An independent audit replayed normal and optimized modes byte-for-byte,
  rederived the scalar terminal and all implication directions, and rebuilt
  the 37,702-cell carrier directly on the full Z/LZ grid without inherited
  ranges or serialized low/high tables.  Every one of the 24 projected
  supports exceeds ceil(d/7), with minimum actual slack six; hashes,
  dependencies, the disjoint ledger, endpoint typing, and docs all pass.
depends_on:
  - THM-3102-projected-k3-z233-exact-screen-and-complete-cell-cardinality-descent
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
script: 04-computation/lrc14_j7_k3_z232_exact_screen_complete_cell_cardinality_descent_thm3106.py
output: 05-knowledge/results/lrc14_j7_k3_z232_exact_screen_complete_cell_cardinality_descent_thm3106.out
script_sha256: 9c38e808e22c9ac376217b9c76da69f198a14a4060cd4af4bf3e10b2c6a604f6
output_sha256: 286bb3e31e1ef28af5640a8ebe5e0df5e57b3e7253aa4c1e5ac0483beb9d0e63
semantic_sha256: a14adcd1e52323baf2b791f55d5846c4fbd422e3441fa2411138bdd98acca3d3
hash_basis: LF-normalized bytes
---

# THM-3106 -- projected k3 z232 exact screen and complete-cell cardinality descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Candidate statement

In the pinned THM-2941 projected `k=3` necessary atlas, the complete
`z_1=232` layer is

```text
3 rows = 3 wall + 0 order:                                (1)

E_A=(1,2,6,9,12,14),    L=3528,   H=348;
E_B=(1,5,6,9,12,14),    L=17640,  H=1738;
E_C=(1,9,10,11,12,14),  L=194040, H=19111.
```

Here `H=floor(13L/132)+1` is the inherited high floor.  A fresh evaluation
through the promoted THM-3078 screen gives the exact partition

```text
147 states = 58 crude + 65 exact-status + 24 residual.     (2)
```

The `24` residual states all belong to `E_C`.  Their relaxed two-high gap is
strictly positive.  Since every row in `(1)` is a wall row, the inherited
wall gate requires at least one high suffix label; therefore the residual
sector has exactly one high label.  The complete one-high enumeration has
`24` cases, and complete-cell translated-band cardinality closes every one.
Thus

```text
all three projected k=3, z_1=232 necessary rows are empty. (3)
```

Composition with THM-3102 gives the disjoint layer subtraction

```text
374325-3=374322,
projected k=3 cap: z_1<=231.                               (4)
```

The next layer is occupied by `9=5+4` wall/order rows.  It is not screened or
closed here.

## 2. Pinned universe and the `147`-state screen

The companion pins the promoted source, output, and semantic hashes of
THM-3078 and THM-3102.  It reparses all `6,060` structured rows in the
THM-2941 atlas, whose LF-normalized hash is

```text
cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda. (5)
```

The neighboring census is recomputed as

```text
z_1=233: 62 = 45 wall + 17 order;
z_1=232:  3 =  3 wall +  0 order;
z_1=231:  9 =  5 wall +  4 order.                          (6)
```

The exact ordered tuple of the three rows in `(1)` has digest

```text
9bda2b6c4582541b8303156c2c9a47bf8b44d6d0ca9ceeb0b30d78860eebb796. (7)
```

The rowwise screen is:

| body | states | crude | status | residual | Farkas route | stage SHA-256 |
|---|---:|---:|---:|---:|---|---|
| `E_A` | `2` | `0` | `2` | `0` | `2` legacy | `113730b2446c1be0da689117f6f408374520700cb0638e2514b89249f920a7d1` |
| `E_B` | `2` | `2` | `0` | `0` | none | `6ccc7a3c1f50661717767cfa5f0bd575af1a567bd7d8695ff89756197cf7b1ca` |
| `E_C` | `143` | `56` | `63` | `24` | `63` legacy | `883c78fa8e0ef3b5259a7f33934cb30de90acf4c6c591eff329f34099285f4a6` |

Every one of the `65` status exclusions passes the promoted legacy
full-table exact Farkas verifier.  THM-3078's narrowly repaired zero-good
direct branch is never invoked:

```text
direct Farkas certificates:  0;
legacy Farkas certificates: 65.                            (8)
```

The complete screen-record digest is

```text
3ca9c8569f05d8adf0ff04586347700597d7875794525452b13031e15acd909a. (9)
```

The direction of this screen is essential.  Its ray quotient relaxes global
suffix order while preserving the fixed first label `232`, distinct
maximizing labels, their denominators, and the wall at-least-one-high gate.
Thus every actual projected packet maps into the enumerated state set.
Crude upper exclusions and exact Farkas contradictions close that enlarged
set; no converse implication from the quotient is used.

## 3. The exact `24`-mask boundary and high-count direction

The complete residual bank is supported on `E_C` and has SHA-256

```text
588fa690d15f18a29083f40f422507462f895c8880df0fb75c6ea8c394189a66. (10)
```

It is exactly

```text
(d,9702,10780,24255),

d in {8,24,40,56,72,88,120,168,264,280,360,440,
      504,616,792,840,1320,1848,2520,3080,3960,5544,9240},

and the final mask (9702,10780,24255,27720).                (11)
```

Let `u_d` be the rigorous high-ray supremum for denominator `d`.  The
terminal grants repeated labels and gives both relevant suffix slots their
independent suprema.  Its duplicate-permitting two-high gap is therefore an
upper relaxation of every physical assignment with at least two high labels.
The exact minimum is

```text
918906100793/174485171984760 > 0,                          (12)
```

with denominator-mask witness `(9702,10780,24255,27720)`.  Hence even this
enlarged state set contains no assignment with two or more high labels.

The wall grammar in `(1)` independently forces at least one high suffix
label.  Equations `(1)` and `(12)` therefore leave exactly one high label.
The one-high enumerator retains every finite low label, including negative
scalar contributions, enforces distinct literal labels, and grants the high
slot its ray supremum.  Its `24` cases are a superset of actual literal
one-high assignments.  They have one fixed low-label set,

```text
(234,260),                                                (13)
```

and exact terminal profile

```text
24 cases = 23 coarse-cardinality + 1 exact-cardinality;
0 max-gap; 0 failures; 0 unit checks;
minimum terminal certificate slack: 1.                    (14)
```

The complete terminal-record and per-case digests are

```text
56cd0b886afdb39ca3423c0cc5240f06298075138f6336ab31e65fbd65e2907c,
f44f893f55db7160d032758170d876417de1fdd2ac1e69268fcea0eda295b495. (15)
```

There are also `23` zero-high scalar passes.  They are hostile controls, not
closures: the inherited wall gate excludes them because it requires a high
suffix label.  They are not included among the `24` one-high certificates.

## 4. Direct full-grid carrier and complete-cell implication

Let `L=194040`.  For the fixed labels

```text
E_C union {232,234,260},                                  (16)
```

the companion directly scans every cell `j in Z/LZ`.  Because every fixed
label is below `L`, the complete closed cell represented by `j` is safe from
the strict-open comb of a label `y` exactly when, for
`r=jy mod L`,

```text
14r>=L,                    14(r+y)<=13L.                  (17)
```

The resulting fixed-safe carrier `J` has

```text
|J|=37702,
SHA-256(J)=656bb385012dd08368a49ef6d365dbaa34362d6e6707313afbfb20128e4945fa.
                                                               (18)
```

This full-grid implementation is compared cell-for-cell with both the
inherited scalar predicate and the vectorized body-range implementation;
all three sets are identical.

For a possible high denominator `d|L`, project the carrier modulo `d`:

```text
S_d={j mod d:j in J},                 kappa(d)=ceil(d/7). (19)
```

One residue contains at most `L/d` cells, giving the coarse lower bound

```text
|S_d|>=ceil(|J|/(L/d)).                                   (20)
```

For `d=8`, `(20)` equals `2=kappa(8)`, so the support is computed exactly:
`|S_8|=8`, with actual slack `6`.  For the other `23` cases `(20)` already
exceeds `kappa(d)`; its smallest certificate slack is `1`, at `d=24`.
The companion nevertheless computes every exact support as a hostile audit:

| high denominator `d` | exact `|S_d|` | `kappa(d)` | terminal method |
|---:|---:|---:|---|
| `8` | `8` | `2` | exact |
| `24` | `24` | `4` | coarse |
| `40` | `40` | `6` | coarse |
| `56` | `56` | `8` | coarse |
| `72` | `72` | `11` | coarse |
| `88` | `88` | `13` | coarse |
| `120` | `120` | `18` | coarse |
| `168` | `168` | `24` | coarse |
| `264` | `264` | `38` | coarse |
| `280` | `280` | `40` | coarse |
| `360` | `360` | `52` | coarse |
| `440` | `440` | `63` | coarse |
| `504` | `504` | `72` | coarse |
| `616` | `616` | `88` | coarse |
| `792` | `792` | `114` | coarse |
| `840` | `840` | `120` | coarse |
| `1320` | `1320` | `189` | coarse |
| `1848` | `1848` | `264` | coarse |
| `2520` | `2520` | `360` | coarse |
| `3080` | `3080` | `440` | coarse |
| `3960` | `3960` | `566` | coarse |
| `5544` | `5544` | `792` | coarse |
| `9240` | `9240` | `1320` | coarse |
| `27720` | `21846` | `3960` | coarse |

Thus the minimum actual support slack is `6`, and every one of the `24`
supports satisfies

```text
|S_d|>kappa(d).                                           (21)
```

The complete support-record digest, including all `24` support hashes, is

```text
828f0631b5743d0ecbf8f767a5f187bc373cb7aec20bbdf4850f3b81fc801772. (22)
```

THM-2984 proves that an arbitrarily translated strict high-danger band has
capacity `kappa(d)=ceil(d/7)`.  This is the translated count, not the smaller
centered quantity `beta(d)`.  Equation `(21)` therefore supplies, at every
local coordinate and for every high label on its exact denominator ray, a
fixed-safe complete cell missed by the high comb.  Hence

```text
P_(E_C,Z)=T.                                              (23)
```

The direction from THM-2941 `(25g)` used here is

```text
literal completion  ==>  P_(E_C,Z) subset U_A.            (24)
```

It is not the converse and is not inferred from an intermediate scalar
mass.  THM-1166 bounds the union for the three distinct aligned labels by

```text
mu(U_A)<=36/91<1.                                         (25)
```

Equations `(23)--(25)` contradict one another and close every case in
`(14)`.

The weak inequalities in `(17)` are load-bearing.  Their least endpoint
slack is exactly zero at cell `14850`, body label `14`, left endpoint.
Replacing them by strict inequalities would discard a genuine safe cell.

## 5. Exact evidence

The principal digests are:

| object | SHA-256 |
|---|---|
| atlas | `cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda` |
| `z_1=232` row order | `9bda2b6c4582541b8303156c2c9a47bf8b44d6d0ca9ceeb0b30d78860eebb796` |
| screen record | `3ca9c8569f05d8adf0ff04586347700597d7875794525452b13031e15acd909a` |
| residual bank | `588fa690d15f18a29083f40f422507462f895c8880df0fb75c6ea8c394189a66` |
| terminal record | `56cd0b886afdb39ca3423c0cc5240f06298075138f6336ab31e65fbd65e2907c` |
| support record | `828f0631b5743d0ecbf8f767a5f187bc373cb7aec20bbdf4850f3b81fc801772` |
| direct carrier | `656bb385012dd08368a49ef6d365dbaa34362d6e6707313afbfb20128e4945fa` |
| `z_1=231` row order | `ff12b3d3c41a10aa38da0a1483bfb506b4c6a9a771c2325fb4df1df867088f5e` |

Run

```text
python 04-computation/lrc14_j7_k3_z232_exact_screen_complete_cell_cardinality_descent_thm3106.py --processes 3
python -O 04-computation/lrc14_j7_k3_z232_exact_screen_complete_cell_cardinality_descent_thm3106.py --processes 3 --output <optimized-output>
```

The companion contains no truth-bearing Python `assert`.  Fresh normal and
optimized runs pass every frozen gate; their transcripts are LF-byte-identical
to the stored output and end in `all_exact_controls=PASS`.  Their hashes are

```text
script:   9c38e808e22c9ac376217b9c76da69f198a14a4060cd4af4bf3e10b2c6a604f6
output:   286bb3e31e1ef28af5640a8ebe5e0df5e57b3e7253aa4c1e5ac0483beb9d0e63
semantic: a14adcd1e52323baf2b791f55d5846c4fbd422e3441fa2411138bdd98acca3d3
```

## 6. Exact handoff and scope

The `z_1=231` census in `(6)` consists exactly of:

| body | `L` | `H` | branch |
|---|---:|---:|---|
| `(1,2,3,8,10,12)` | `1680` | `166` | order |
| `(1,2,4,8,10,12)` | `1680` | `166` | order |
| `(1,2,4,9,10,12)` | `2520` | `249` | wall |
| `(1,2,5,9,10,12)` | `2520` | `249` | wall |
| `(1,3,4,8,10,12)` | `1680` | `166` | order |
| `(1,4,5,9,10,12)` | `2520` | `249` | wall |
| `(1,4,8,9,10,12)` | `5040` | `497` | wall |
| `(1,5,9,11,12,13)` | `360360` | `35491` | wall |
| `(2,3,4,8,10,12)` | `1680` | `166` | order |

This table is only an occupied handoff.  None of its nine rows has been
screened by the companion, and no emptiness claim is made at `z_1=231`.

**Scope.**  This theorem acts only in the pinned projected
`k=3` necessary atlas.  It does not classify physical covers outside that
projection, say anything about arbitrary `k<=1` packets or the final rung,
or prove LRC(14).

QED.

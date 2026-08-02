---
id: THM-3109
title: "Projected k3 z231 exact screen and complete-cell cardinality descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  All nine rows
  in the pinned projected k=3, z1=231 necessary layer are empty.  The
  projected ledger is 374313 and its cap is z1<=230.  The fifty z1=230
  rows remain occupied.  No LRC(14) claim is made.
source: root/lrc14-projected-k3-z231-2026-08-02
audit: >
  An independent referee rebuilt all nine screens and the sole terminal,
  scanned every one of the 360360 complete cells without serialized carrier
  data, and checked both projected supports.  A first immutable audit caught
  and repaired a reversed prose inclusion: the weak complete cells form an
  inner carrier in the strict-open safe set.  Repaired normal, optimized,
  and stored transcripts are byte-identical; dependencies, hashes, semantic
  digest, directions, ledger, next-layer census, docs, and diff checks pass.
  A later MISTAKE-331/333 evidence repair retained every exact certificate
  check but repinned the screen/semantic records to the canonical
  nineteen-field instance/result rows and basis-invariant branch counts.
depends_on:
  - THM-3106-projected-k3-z232-exact-screen-and-complete-cell-cardinality-descent
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
script: 04-computation/lrc14_j7_k3_z231_exact_screen_complete_cell_cardinality_descent_thm3109.py
output: 05-knowledge/results/lrc14_j7_k3_z231_exact_screen_complete_cell_cardinality_descent_thm3109.out
script_sha256: 1f74f2b2368c04f514f2c388b54c70a9ee66c9387fbc437093884b807b3eb23c
output_sha256: e6d56fa6a419ffe229b8334090e02d98c9d2cdf5f2fa5e24baedd5f722dadf70
semantic_sha256: 5be5c2fc680d6873600e77227f51264d74a7cd353652795e5ef74215e0fda843
hash_basis: LF-normalized bytes
---

# THM-3109 -- projected k3 z231 exact screen and complete-cell cardinality descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The replay evidence was later repinned under MISTAKE-331/333.  Exact rational
dual verification is unchanged; persisted hashes exclude the optimizer's
noncanonical basis choice and contradiction scale.

## 1. Exact layer statement

In the pinned THM-2941 projected `k=3` necessary atlas, the complete
`z_1=231` layer consists of

```text
9 rows = 5 wall + 4 order.                                (1)
```

Their ordered tuple has SHA-256

```text
ff12b3d3c41a10aa38da0a1483bfb506b4c6a9a771c2325fb4df1df867088f5e. (2)
```

A fresh THM-3078 screen gives

```text
226 states = 127 crude + 97 exact-status + 2 residual.    (3)
```

Eight rows are empty already at this stage.  The two residual masks lie on
the single wall body

```text
E=(1,5,9,11,12,13),       L=360360,       H=35491.       (4)
```

Their duplicate-permitting two-high gap is positive, so the wall gate
forces exactly one high suffix.  The two resulting one-high cases have one
fixed low-label set `(243,351)`.  A direct scan of all `360360` complete
cells gives a `67700`-cell weak safe carrier; its two projected supports
strictly exceed the translated high-band capacities.  THM-2941 completion
and THM-1166 then contradict literal completion.  Therefore

```text
all nine projected k=3, z_1=231 necessary rows are empty. (5)
```

The disjoint ledger consequence is

```text
374322-9=374313,
projected k=3 cap: z_1<=230.                              (6)
```

The next layer has `50=48+2` wall/order rows and is not screened here.

## 2. The complete screen

The companion pins the source, output, and semantic hashes of THM-3106,
then reparses all `6060` structured rows of the THM-2941 atlas.  The exact
rowwise result is:

| body | branch | states | crude | status | residual |
|---|---|---:|---:|---:|---:|
| `(1,2,3,8,10,12)` | order | 3 | 0 | 3 | 0 |
| `(1,2,4,8,10,12)` | order | 1 | 0 | 1 | 0 |
| `(1,2,4,9,10,12)` | wall | 1 | 0 | 1 | 0 |
| `(1,2,5,9,10,12)` | wall | 1 | 1 | 0 | 0 |
| `(1,3,4,8,10,12)` | order | 19 | 0 | 19 | 0 |
| `(1,4,5,9,10,12)` | wall | 3 | 2 | 1 | 0 |
| `(1,4,8,9,10,12)` | wall | 6 | 0 | 6 | 0 |
| `(1,5,9,11,12,13)` | wall | 191 | 124 | 65 | 2 |
| `(2,3,4,8,10,12)` | order | 1 | 0 | 1 | 0 |

All `97` status exclusions pass the inherited exact legacy Farkas verifier;
the repaired direct branch is not invoked.  The complete screen-record
digest is

```text
8fa535bb4e987a6ee21657503b769c202d743d99577cf8e05019a5ec93f525f1. (7)
```

The ray quotient enlarges the physical assignment universe while retaining
the fixed first label, denominator rays, distinct maximizing labels, and
wall/order grammar.  Crude estimates and exact Farkas contradictions are
therefore applied only to a superset of actual projected packets.  No
converse from the quotient is used.

## 3. The two-mask terminal

The complete residual bank is

```text
(1560,3080,10296,40040),
(1560,3080,40040,51480),                                (8)
```

with digest

```text
89963fe483156389680ed3a1bf1c6feb92fdc00f72ce6cffaccb3d6ca1e88cc8. (9)
```

Granting repeated labels and independent high-ray suprema makes the
two-high test an upper relaxation.  Its exact gap is

```text
53144999/13513175949 > 0.                               (10)
```

Thus no assignment has two high suffixes.  Since `(4)` is a wall row, at
least one high suffix is mandatory, leaving exactly one.  Two zero-high
scalar passes are retained as hostiles but excluded by that wall gate.

The exact one-high enumeration has two cases, both with literal lows

```text
(3080,351),                 (40040,243).                (11)
```

The high denominators are `10296` and `51480`; all negative low
contributions remain present.  Both cases pass the inherited coarse
cardinality gate, whose minimum certificate slack is `464`.  The complete
case digest is

```text
d75fcc91dfa1ed6c2e1240c609376f74e077457d61caac16da36a36d9ff89b91. (12)
```

## 4. Independent complete-cell reconstruction

For every `j in Z/360360Z`, the companion tests the fixed labels

```text
E union {231,243,351}.                                   (13)
```

Writing `r=jy mod L`, the complete closed cell is weakly safe from the
strict-open comb of label `y` exactly when

```text
14r>=L,                      14(r+y)<=13L.               (14)
```

The direct full-grid set has

```text
|J|=67700,
SHA-256(J)=02c3c649a0cdd6808479821077de9fe142fde6efc7fc6e3daa5380276d403b93.
                                                               (15)
```

It agrees cell-for-cell with both inherited scalar and vector
implementations.  The weak inequalities are load-bearing: the minimum
endpoint slack is zero at cell `35100`, body label `11`, left endpoint.
Because the danger comb is strict-open, every closed cell selected by
`(14)` is wholly contained in the safe set.  Thus `J` is an inner
complete-cell carrier; it need not contain safe points from partial boundary
cells.

For a high denominator `d`, project `J` modulo `d` and compare with the
translated high-band capacity `kappa(d)=ceil(d/7)`:

| `d` | exact support | `kappa(d)` | actual slack | support SHA-256 |
|---:|---:|---:|---:|---|
| `10296` | `10296` | `1471` | `8825` | `74673614815c2bcf6fbdb23d98e5e20eb4edbf37d3c40243eb3d3b199a39bf37` |
| `51480` | `39624` | `7355` | `32269` | `e28a6df29f568e4b3835555cfe11ab504413cea3202a19e7f5f575b86eb2fb5f` |

In both cases the support is larger than any translated strict high-danger
band.  Hence every translated high band misses at least one whole
fixed-safe cell, and THM-2984 gives `P_(E,Z)=T`.  The implication from
THM-2941 used here is only

```text
literal completion  ==>  P_(E,Z) subset U_A.            (16)
```

THM-1166 supplies `mu(U_A)<=36/91<1`, contradicting `(16)` and closing both
one-high cases.

## 5. Exact evidence and scope

Run

```text
python 04-computation/lrc14_j7_k3_z231_exact_screen_complete_cell_cardinality_descent_thm3109.py --processes 3
python -O 04-computation/lrc14_j7_k3_z231_exact_screen_complete_cell_cardinality_descent_thm3109.py --processes 3 --output <optimized-output>
```

The companion contains no truth-bearing Python `assert`.  Fresh normal and
optimized runs recompute the whole screen and carrier, pass every frozen
gate, and byte-match the stored transcript.  Their LF hashes are

```text
script:   1f74f2b2368c04f514f2c388b54c70a9ee66c9387fbc437093884b807b3eb23c
output:   e6d56fa6a419ffe229b8334090e02d98c9d2cdf5f2fa5e24baedd5f722dadf70
semantic: 5be5c2fc680d6873600e77227f51264d74a7cd353652795e5ef74215e0fda843
```

The `z_1=230` layer has `50=48` wall `+2` order rows, ordered digest

```text
f5129616fe3d7f5884350abe7e77bb2684798e213600e266cd456bc9c582f058. (17)
```

This is an occupied, unscouted handoff only.  The theorem acts solely in the
pinned projected `k=3` necessary atlas.  It does not classify physical
covers outside that projection, treat arbitrary `k<=1` packets or the final
rung, or prove LRC(14).

**QED.**

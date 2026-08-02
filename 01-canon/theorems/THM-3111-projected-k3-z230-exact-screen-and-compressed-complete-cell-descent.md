---
id: THM-3111
title: "Projected k3 z230 exact screen and compressed complete-cell descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  All fifty rows
  in the pinned projected k=3, z1=230 necessary layer are empty.  The
  projected ledger is 374263 and its cap is z1<=229.  No LRC(14) claim is
  made.
source: root/lrc14-projected-k3-z230-2026-08-02
audit: >
  Independent checkpoints reconstructed all fifty screens, including a
  separate replay of the three largest rulers, then rebuilt all six
  terminals and eight direct/scalar/vector carriers.  Hostile audit caught
  and repaired an unpromoted raw-case versus case-certificate digest type
  error.  Repaired normal, optimized, and stored transcripts are
  byte-identical; dependencies, hashes, directions, ledger, handoff, docs,
  and diff checks pass.  A later MISTAKE-331/333 evidence repair retained all
  exact dual checks while repinning the screen and semantic records to the
  canonical nineteen-field instance/result rows and basis-invariant branch
  counts.
depends_on:
  - THM-3109-projected-k3-z231-exact-screen-and-complete-cell-cardinality-descent
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
script: 04-computation/lrc14_j7_k3_z230_exact_screen_compressed_complete_cell_descent_thm3111.py
output: 05-knowledge/results/lrc14_j7_k3_z230_exact_screen_compressed_complete_cell_descent_thm3111.out
script_sha256: 42323171481deba2371eed9947b2079976cb367dac340cf58b8f1f0c0afb5082
output_sha256: 54a82d696c592162bbe3f98a3dd34e092967a0a6320e9931adb2866570cc5813
semantic_sha256: 4ff290e285dbb748dac71e1b885ce220dbfe04ec6f236ea97a5526bc27baa497
hash_basis: LF-normalized bytes
---

# THM-3111 -- projected k3 z230 exact screen and compressed complete-cell descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The replay evidence was later repinned under MISTAKE-331/333.  Every exact
rational dual is still rebuilt and checked, but persisted hashes no longer
depend on a solver-selected basis or contradiction normalization.

## 1. Exact layer statement

In the pinned THM-2941 projected `k=3` necessary atlas, the complete
`z_1=230` layer has

```text
50 rows = 48 wall + 2 order,                              (1)
```

with ordered-tuple digest

```text
f5129616fe3d7f5884350abe7e77bb2684798e213600e266cd456bc9c582f058. (2)
```

A fresh THM-3078 screen gives

```text
4156 states = 2437 crude + 1648 exact-status + 71 residual. (3)
```

All `1648` status exclusions pass the inherited exact verifier.  The two
order rows contribute `8=8+0+0` states and are already empty.  The
seventy-one residual masks occur on exactly six wall bodies.  Every one of
their duplicate-permitting two-high gaps is positive.  The wall gate then
forces exactly one high suffix, and the resulting seventy-one cases factor
through only eight fixed `(body, low-label pair)` carriers.  Sixty-eight
cases close by a coarse complete-cell lower bound and three tied cases by
exact projected support.  Therefore

```text
all fifty projected k=3, z_1=230 necessary rows are empty. (4)
```

The disjoint ledger consequence is

```text
374313-50=374263,
projected k=3 cap: z_1<=229.                              (5)
```

The next layer has `43=36+7` wall/order rows and ordered digest

```text
7449dd7ad70cf3c76c32edb2cc509e29989ac008c2e9a9968bceaabf3e2a2587. (6)
```

## 2. The complete screen and the six wall residues

The companion pins THM-3109, reparses all `6060` atlas rows, and recomputes
all fifty tasks.  Its complete screen-record digest is

```text
3c7681f663bbb9bbe6f0483918474df76c935dc6a9a6006499848b9898d57477. (7)
```

The residual profile is:

| body | `L` | masks | zero-high hostiles | one-high cases | low-pair carriers |
|---|---:|---:|---:|---:|---:|
| `(1,7,9,10,11,12)` | `194040` | `19` | `18` | `19` | `1` |
| `(1,7,10,11,12,13)` | `840840` | `12` | `12` | `12` | `1` |
| `(1,8,10,11,12,14)` | `129360` | `4` | `4` | `4` | `1` |
| `(1,9,10,11,12,14)` | `194040` | `34` | `34` | `34` | `3` |
| `(2,9,10,11,12,14)` | `194040` | `1` | `0` | `1` | `1` |
| `(3,9,10,11,12,14)` | `194040` | `1` | `0` | `1` | `1` |

The ray quotient enlarges the physical assignment universe while retaining
the fixed first label, denominator rays, distinct maximizing labels, and
wall/order grammar.  Crude and exact-status contradictions therefore close
only a superset of actual projected packets.  No converse from the quotient
is used.

## 3. The exact-one-high gate

Granting repeated labels and independent high-ray suprema makes the
two-high calculation an upper relaxation.  Its six exact gaps are

| body | duplicate-two-high gap |
|---|---:|
| `(1,7,9,10,11,12)` | `113487288713/32897729721780` |
| `(1,7,10,11,12,13)` | `5173786082623/1879348592763460` |
| `(1,8,10,11,12,14)` | `1349719633/341778176376` |
| `(1,9,10,11,12,14)` | `8059115839/2508612085980` |
| `(2,9,10,11,12,14)` | `14389939/3132625860` |
| `(3,9,10,11,12,14)` | `49827592/15730119405` |

All are positive, so no residual assignment can have two high suffixes.
The companion explicitly checks that all six residual bodies are wall rows;
each therefore has at least one high suffix.  Exactly one high suffix remains.

The `68` zero-high scalar passes are retained as hostile witnesses and are
excluded solely by this wall premise.  Exhaustive one-high enumeration gives
exactly `71` cases, keeps every negative low contribution, and still uses a
high-ray supremum, so it remains an upper relaxation rather than a converse.
The aggregate six-tuple terminal-record digest is

```text
90b1768cf0165683ceb086359180563081a0ff89a3876c88a47c7240ca036e37. (8)
```

## 4. Eight complete-cell carriers

For a wall body `E`, low-label pair `B`, and ruler `L`, let `J(E,B)` be the
set of complete cells `j in Z/LZ` for which every label

```text
E union {230} union B                                      (9)
```

satisfies, with `r=jy mod L`,

```text
14r>=L,                         14(r+y)<=13L.              (10)
```

The weak inequalities are deliberate.  The comb dangers are strict-open,
so every closed complete cell selected by `(10)` is wholly contained in the
fixed safe set.  Thus `J(E,B)` is an inner carrier; it need not contain safe
points lying in partial boundary cells.

The seventy-one one-high cases group into the following eight carriers:

| body | lows `B` | `|J(E,B)|` | minimum weak endpoint |
|---|---|---:|---|
| `(1,7,9,10,11,12)` | `(234,260)` | `37486` | `(0,13860,1,'L')` |
| `(1,7,10,11,12,13)` | `(312,410)` | `155324` | `(0,102101,230,'R')` |
| `(1,8,10,11,12,14)` | `(260,312)` | `24894` | `(0,9900,14,'L')` |
| `(1,9,10,11,12,14)` | `(260,312)` | `37772` | `(0,14850,14,'L')` |
| `(1,9,10,11,12,14)` | `(234,260)` | `38078` | `(0,14850,14,'L')` |
| `(1,9,10,11,12,14)` | `(243,260)` | `37950` | `(0,14850,14,'L')` |
| `(2,9,10,11,12,14)` | `(234,260)` | `38322` | `(0,6930,2,'L')` |
| `(3,9,10,11,12,14)` | `(243,260)` | `41740` | `(0,4620,3,'L')` |

Every direct full-grid set agrees cell-for-cell with the inherited scalar and
vector implementations.  The carrier-record digest is

```text
a8e749413d276397dd8ba521a62bad219e82639d7c53669030c41ec166ebeb5d. (11)
```

Grouping is only by `(body,B)`: it reuses a common carrier without
identifying distinct high cases or mixing bodies.

## 5. Coarse and exact translated support

For a high denominator `d`, a translated strict high-danger band occupies at
most

```text
kappa(d)=ceil(d/7)                                        (12)
```

classes modulo `d`.  This is the translated capacity; the centered `beta`
forbidden by MISTAKE-334 is not used.  If

```text
ceil(|J|/(L/d)) > kappa(d),                               (13)
```

the coarse pigeonhole lower bound already leaves a whole fixed-safe cell
outside every translated high band.  This closes `68` cases; the minimum
coarse lower-bound margin is `2`.

In precisely three equality cases `(13)` is inconclusive.  Direct support
projection gives

| body | `d` | exact support | `kappa(d)` | slack |
|---|---:|---:|---:|---:|
| `(1,7,9,10,11,12)` | `8` | `8` | `2` | `6` |
| `(2,9,10,11,12,14)` | `8` | `8` | `2` | `6` |
| `(3,9,10,11,12,14)` | `2` | `2` | `1` | `1` |

Thus all three also leave a whole fixed-safe cell outside every translated
high band.  The full seventy-one-case direct-record digest is

```text
c2640270e1af368ea4c40856ef9b8254d16214bea84ea2be6197f75f6d17e761. (14)
```

THM-2984 then gives `P_(E,Z)=T`.  The only THM-2941 direction used is

```text
literal completion  ==>  P_(E,Z) subset U_A.             (15)
```

Together `P_(E,Z)=T` and `(15)` give `T subset U_A`, contradicting
THM-1166's `mu(U_A)<=36/91<1`.  This closes every one-high case and proves
`(4)`.

## 6. Boundary, evidence, and scope

The carrier compression begins only after the positive two-high gap and wall
gate have forced exactly one high label.  It forgets joint two-high
compatibility.  If a later layer has a nonpositive duplicate-two-high gap,
these eight-carrier arguments do not apply without a common two-high carrier
or another compatibility sidecar.

The companion contains no truth-bearing Python `assert`.  Its semantic digest
is

```text
4ff290e285dbb748dac71e1b885ce220dbfe04ec6f236ea97a5526bc27baa497. (16)
```

An independent reconstruction reproduces all fifty ordered screen tuples;
a separate replay of the three largest rulers has selected-screen digest
`78937062dc5092adde8cb4ae53db6c7ba2879c57cb7da2916923f634fe082a35`.
The hostile audit also caught an unpromoted evidence-type error: the first
candidate compared raw `sha256(repr(cases))` values with inherited row-17
case-plus-certificate digests.  The repaired companion pins and reports the
two objects separately, requires every residual body to be a wall, and
distinguishes certificate, coarse-lower-bound, and exact-support slack.

Run

```text
python 04-computation/lrc14_j7_k3_z230_exact_screen_compressed_complete_cell_descent_thm3111.py --processes 12
python -O 04-computation/lrc14_j7_k3_z230_exact_screen_compressed_complete_cell_descent_thm3111.py --processes 12 --output <optimized-output>
```

The companion contains no truth-bearing Python `assert`.  Fresh repaired
normal and optimized runs recompute all fifty screens, all six terminals,
and all eight carriers and byte-match the stored transcript.  Their LF hashes
are

```text
script:   42323171481deba2371eed9947b2079976cb367dac340cf58b8f1f0c0afb5082
output:   54a82d696c592162bbe3f98a3dd34e092967a0a6320e9931adb2866570cc5813
semantic: 4ff290e285dbb748dac71e1b885ce220dbfe04ec6f236ea97a5526bc27baa497
```

This result acts solely in the pinned projected `k=3` necessary atlas.  It
does not classify physical covers outside that projection, treat arbitrary
`k<=1` packets or the final rung, or prove LRC(14).

**QED.**

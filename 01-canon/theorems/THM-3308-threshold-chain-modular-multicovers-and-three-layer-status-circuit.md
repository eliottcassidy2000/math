---
id: THM-3308
title: "Threshold-chain modular multicovers and the exact three-layer status circuit"
status: >
  PROVED + VERIFIED-EXACT, in the declared finite projected atlas.  All `271`
  exact common-status eliminations in the second intrinsic-ruler cost prefix
  at projected `k=3,z1=216` admit deterministic zero-constant nonnegative
  integral modular multicover certificates.  Their minimum binary-layer
  support taxonomy is `260` one-layer, `10` two-layer, and `1` three-layer.
  The eleven states formerly presented as a weighted exact-Farkas residue are
  exactly the multi-layer states.  Exact feasible common tables for all `59`
  addressed smaller tail supports prove that these eleven support minima hold
  in the full rational dual cone, not only in the integral modular template.
  This recertifies existing projected status kills and changes no row ledger;
  it is not a physical-entry, arbitrary-`k`, rung, or LRC(14) theorem.
audit: >
  The exact companion pins the independent cost-prefix audit and reconstructs
  all `271` deterministic instances from its older sources, without using or
  hashing the inherited solver-selected dual vector.  It exhausts bounded
  integral modular covers, checks every pointwise Boolean inequality on all
  sixteen patterns, proves coordinate minimality and primitive normalization,
  and constructs lexicographically deterministic exact rational feasible
  tables for every proper tail support of each distinct multi-layer instance.
  Row `64` is a hostile control: all `57` actual status kills have one-layer
  covers, while the first of eight status survivors has an explicit feasible
  common table and no multicover contradiction.  Normal and optimized modes
  byte-match the stored output; the source has zero assertion nodes and zero
  floating literals.
source: root/creative-synthesis-recover/2026-08-03
depends_on:
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
  - THM-3264-projected-k3-z216-low-cost-gcd8-seventeen-row-terminal-descent
  - THM-3281-projected-k3-z216-three-natural-wall-family-screen-descent
related:
  - THM-3282-complete-boolean-seam-decoder-width-spectrum-and-maximal-cut-class
script: 04-computation/lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.py
output: 05-knowledge/results/lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.out
script_sha256: e7462eabab773133688079141708d2377742e2d6f9a9480089666956f472020c
output_sha256: f83d7ca1a3e2557330ec7ed655dc756cef21d6abdb52326d3623563c2abc407d
semantic_sha256: bd7ac011d960afe7c4566ac9c17e9df3d8c2cf710e6023a715cb6a6d1ce1408d
hash_basis: LF-normalized bytes
---

# THM-3308 -- threshold-chain modular multicovers and the exact three-layer status circuit

**PROVED + VERIFIED-EXACT in the declared finite projected atlas.**

The theorem replaces `271` generic LP dual verdicts by a deterministic
certificate compiler.  Its intrinsic complexity coordinate is the number of
nested load layers needed by a support-minimal contradiction.

## 1. Common-table relaxation

Fix one of the `271` exact-status instances in the second audited intrinsic-
ruler cost prefix at projected `k=3,z1=216`.  It consists of

```text
q                         total status mass,
m_i, 0<=i<4              four prescribed one-marginals,
c(P), P subseteq [4]     capacity of each Boolean status pattern,
H_t                       required target mass at load at least t.           (1)
```

A common status table is a nonnegative rational vector `(x_P)` satisfying

```text
sum_P x_P=q,
sum_(P contains i) x_P=m_i,                                  (2)
sum_(P:c(P)>=t) x_P >= H_t                                   (3)
```

for every relevant threshold `t`.  The word *common* is essential: all tail
constraints use one joint table.  Separate optimizing tables at different
thresholds do not certify `(2)--(3)`.

## 2. Threshold-chain multicover lemma

For a finite set `T` of thresholds define the layer-count function

```text
f_T(P)=sum_(t in T) 1[c(P)>=t].                              (4)
```

Suppose `w=(w_0,w_1,w_2,w_3)` is a nonnegative integral vector satisfying the
pointwise modular majorization

```text
f_T(P) <= sum_(i in P) w_i             for all 16 patterns P. (5)
```

Multiplying `(5)` by `x_P`, summing, and applying `(2)--(3)` gives

```text
sum_(t in T) H_t <= sum_i w_i m_i.                            (6)
```

Consequently

```text
Delta(T,w)=sum_(t in T)H_t-sum_i w_i m_i > 0                 (7)
```

is an exact Farkas contradiction.  It has coefficient one on each selected
tail row, zero constant term, and nonnegative integral coordinate weights.

If `|T|=k`, then `0<=f_T<=k`.  Any coordinate weight larger than `k` may be
truncated to `k` without invalidating `(5)`.  Hence the bounded search

```text
T subseteq relevant thresholds,       w in {0,...,|T|}^4      (8)
```

is exhaustive for this integral modular template.

## 3. Complete finite taxonomy

Applying `(8)` to all `271` deterministic instances gives

```text
minimum layer support 1: 260 instances,
minimum layer support 2:  10 instances,
minimum layer support 3:   1 instance.                        (9)
```

The preceding solver-free taxonomy refines exactly as follows:

| earlier branch | one layer | two layers | three layers |
|---|---:|---:|---:|
| coordinate union | 248 | 0 | 0 |
| zero-marginal-reduced union | 1 | 0 | 0 |
| two-fan | 11 | 0 | 0 |
| weighted exact-Farkas residue | 0 | 10 | 1 |

Thus the old eleven-state residue is not a separate LP species.  It is exactly
the set of states for which one load threshold is insufficient.

There are eight distinct deterministic instances among the eleven addresses.
Their support-minimal multicovers are:

| atlas row / divisors | multiplicity | `T` | `w` | `Delta` |
|---|---:|---|---|---:|
| `134 / (385,1386,1980,3080)` | 1 | `(3,4)` | `(2,1,1,2)` | 66 |
| `121 / (45,280,455,3640)` | 1 | `(3,4)` | `(2,1,2,1)` | 194 |
| `17 / (49,63,196,252)` | 1 | `(1,2)` | `(2,2,1,1)` | 18 |
| `27 / (49,126,196,252 or 1764)` | 2 | `(3,5)` | `(2,1,1,0)` | 6 |
| `56 / (49,126,196,252)` | 1 | `(2,5)` | `(2,1,1,0)` | 12 |
| `138 / (18,49,196,882)` | 1 | `(1,4,5)` | `(1,3,1,1)` | 11 |
| `138 / (49,63,196,252 or 1764)` | 2 | `(1,2)` | `(2,2,1,1)` | 18 |
| `298 / (49,63,196,252 or 1764)` | 2 | `(1,2)` | `(2,2,1,1)` | 18 |

Every displayed cover is pointwise checked on all sixteen patterns,
coordinatewise minimal, and primitive.

## 4. Full-dual support minimality

Minimum support in the integral template does not by itself rule out a
different rational Farkas vector on fewer tail rows.  The exact companion
therefore audits the primal subsystem behind every proper tail support.

For each of the eight distinct multi-layer instances and every threshold set
of cardinality below the claimed minimum, it constructs a nonnegative
rational table satisfying `(2)` and all selected inequalities `(3)`.  The
construction enumerates exact vertices of the bounded sixteen-column
polytope: the five equality rows in `(2)`, every possible independent subset
of active tail equalities, and the necessary zero coordinates.  All tables
are rechecked directly.

The resulting census is

```text
59 addressed proper tail supports,
47 distinct proper tail supports after repeated-instance identification,
47/47 with an exact feasible common table.                   (10)
```

Primal feasibility rules out every infeasibility certificate supported on
that proper subsystem, independent of dual normalization, signs, constant
term, or rational tail weights.  Therefore the ten support-two claims and the
one support-three claim in `(9)` are minimum-support statements in the full
rational dual cone.

For the exceptional instance

```text
row 138, divisors (18,49,196,882),
q=294, m=(147,42,84,126),
load histogram ((0,42),(1,36),(2,6),(3,36),(4,106),(5,56),(6,12)),          (11)
```

every singleton and pair among its six nontrivial tail rows has an exact
feasible common table.  Yet `T=(1,4,5)` and `w=(1,3,1,1)` give

```text
H_1+H_4+H_5 = 252+174+68 = 494,
w.m = 147+126+84+126 = 483,
Delta=11.                                                     (12)
```

This proves that nested tail events do not force every common-table
obstruction to have a two-layer certificate.

## 5. Submodularity boundary

The layer-count functions selected in the eleven multi-layer addresses are
submodular in eight cases and non-submodular in three.  Exact violations are
supplied by the two row-`27` copies and row `56`; for example row `27` has

```text
f(2)+f(8)-f(10)-f(0)=-1.                                    (13)
```

None of the eleven functions is supermodular.  Since the modular majorants
remain valid in all eleven cases, the proved mechanism is

```text
nested capacity tails -> layer-count function -> modular majorant,          (14)
```

not a polymatroid or submodular-rank reduction.

## 6. Hostile survivor and scope

The older row `64` status stage has

```text
115 states = 50 crude + 57 exact status + 8 residual.         (15)
```

All `57` genuine status eliminations have one-layer modular covers.  For the
first residual state, with divisors `(80,2156,4312,5390)` and `q=43120`, the
exact sparse table

```text
x_0=28046,
x_3=3234, x_5=3234, x_6=2446, x_8=5680, x_14=480             (16)
```

satisfies the total, all marginals, and the required tail.  The compiler
correctly finds no contradiction.  THM-3264 closes that row later by terminal
machinery; `(16)` is a negative control only for the common-status stage.

The theorem recertifies `271` already selected eliminations inside a necessary
projected `k=3,z1=216` atlas.  It changes no ledger count and restores none of
the quotient losses: physical entry, endpoint, owner, phase, current,
arbitrary `k<=1`, the rung, and LRC(14) all remain open.

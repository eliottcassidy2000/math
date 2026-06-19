---
id: HYP-2643
title: LRC(14) fold multiplicity transport - the k=9 near-AP defect is a tiny target shift
status: OPEN; exact scout complete
source: codex-2026-06-19-S29
depends_on:
  - HYP-2642
  - HYP-2640
  - HYP-2638
  - HYP-2637
  - HYP-2639
related:
  - HYP-2641
  - HYP-2118
  - HYP-2122
  - THM-400
  - THM-531
  - OPEN-Q-108
---

# HYP-2643 - LRC(14) Fold Multiplicity Transport

## Claim

The fold multiplicity idea should be stored as a target profile, not as a
scalar count.

For a cluster `E` with `0 in E`, define the nontrivial visible-fold profile

```text
F_E(c) = #{0<a<b in E : a+b=c in E}.
```

The AP-relative transport coordinates are

```text
lost_1(E)   = sum_c max(0, F_AP(c)-F_E(c))/c
gained_1(E) = sum_c max(0, F_E(c)-F_AP(c))/c
net_1(E)    = lost_1(E)-gained_1(E).
```

Raw fold count misses the k=9 boundary: AP9 and the tight non-AP
`E=(0,1,2,3,4,5,6,7,9)` both have `12` nontrivial folds.  The profile sees the
defect exactly:

```text
AP9:      F(8)=3
nearAP9: F(8)=0, F(9)=3
net_1 = 3/8 - 3/9 = 1/24.
```

So the binding non-AP is not a loss of fold count; it is the smallest possible
outward transport of the top fold target.

This complements the newly landed HYP-2642 wall-transfer certificate for the
same endpoint defect.  HYP-2642 measures the exact `L_y` loss by common-wall
missed-sector transfer; HYP-2643 isolates the fold-profile motion behind that
loss before the sector functional is applied.

## Evidence

Script:
`04-computation/lrc14_fold_multiplicity_transport_codex_s29.py`

Stored output:
`05-knowledge/results/lrc14_fold_multiplicity_transport_codex_s29.out`

Single-end defect exact ledgers:

| k | tight row | `L_y` drop from AP | `p0` drop from AP | exact fold transport |
|---:|---|---:|---:|---:|
| 8 | `(0,1,2,3,4,5,6,8)` | `0.120850340` | `0.118877551` | `3/7 - 2/8 = 5/28` |
| 9 | `(0,1,2,3,4,5,6,7,9)` | `0.005587050` | `0.014455782` | `3/8 - 3/9 = 1/24` |
| 10 | `(0,1,2,3,4,5,6,7,8,10)` | `0.007665659` | `0.010600907` | `4/9 - 3/10 = 13/90` |

The exact formula for the `s=1` end-defect
`E_k^+={0,...,k-2,k}` is

```text
net_1(E_k^+) = floor((k-2)/2)/(k-1) - (floor((k-1)/2)-1)/k.
```

For k=9 this is unusually tiny because the transported multiplicity is the
same on both sides: three folds move from target `8` to target `9`.

Bounded-bank evidence:

- In the critical k=9 bank with `max(E)<=13`, the unique top non-AP is also
  the unique row in the tiny `net_1 in (0,0.05)` bucket.
- The next fold-transport bucket starts at `net_1=0.175` and its best row has
  `L_y=0.469349962`, far below the binding row `0.487288990`.
- k=10 shows the same local behavior in the default bank: the top non-AP is
  the smallest positive transport row.
- k=8 is the caveat: its cap is looser, and the top non-AP is a left-hole
  clipped AP `(0,2,3,4,5,6,7,8)` with larger `net_1=0.842857143`. Thus fold
  transport is not a universal scalar order; it is the near-AP boundary
  coordinate for the razor-thin k=9 piece.

## Interpretation

HYP-2640 said raw relation rank is a switch, not a ruler.  HYP-2643 makes the
next refinement: fold multiplicity count is also too coarse.  The useful object
is the fold target measure.  In the k=9 binding row, the profile is nearly AP;
only the top target moves one denominator outward.

This suggests a proof route for the finite near-AP boundary:

```text
non-AP near the AP envelope
  -> positive AP-relative fold transport;
critical k=9 boundary
  -> smallest possible transport is exactly 1/24;
sector functional drop from this transport
  -> exceeds the AP-to-cap slack after finite cap accounting.
```

The last arrow is not proved.  The computation isolates the correct local
ledger for it.

## Proof Route

1. Prove a clipped-AP fold transport lemma: among non-AP k=9 rows in the
   high-`L_y` envelope, the minimum positive `net_1` is `1/24`, uniquely
   realized by `(0,1,2,3,4,5,6,7,9)`.
2. Convert `net_1=1/24` plus the exact target move `3/8 -> 3/9` into a sector
   distribution loss, preferably directly at `p0` or the `L_y` dual.
3. Use HYP-2638's finite small-excess certificate to keep the remaining
   near-AP candidates bounded.
4. Send larger fold-transport defects to the safer envelope or to HYP-2639's
   signed shell machinery.

## Tournament Analysis

Pairwise observable: proof usefulness after rank saturation.  Switch/gauge:
lexicographic priority declared by the exact near-AP separation.

Hamiltonian path:

```text
fold_transport_profile
> reciprocal_fold_mass
> fold_target_centroid
> visible_fold_count
> pair_sum_energy
> raw_relation_rank
> raw_runner_vertices
```

Fingerprint: score histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, no directed
3-cycles, one Hamiltonian path.

Assumption challenge: the vertex set is not runners, arcs, or raw folds.  The
vertices are proof quotients.  The quotient preserves exactly what scalar fold
count destroys: where the fold multiplicity sits.

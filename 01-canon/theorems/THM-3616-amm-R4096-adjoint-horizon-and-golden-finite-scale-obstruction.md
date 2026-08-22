---
id: THM-3616
title: "AMM R=4096 adjoint horizon and golden finite-scale obstruction"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  At the pre-registered hostile epoch R=4096, the archived failing
  Rule-A offset D0=88 dies at row 1014.  Its positive truncated-adjoint wall
  changes sign between cuts 382 and 383, so every admissible continuation
  surviving that fatal inequality must differ from Rule A by row 382.  The
  exact ratio q/R=191/2048 lies below (3-sqrt(5))/8 by
  (-577+256sqrt(5))/2048; its error is larger than at R=2048, refuting the
  proposed monotone finite-scale golden sharpening.  No eventual limit,
  alternative feasibility, uniform extractor, or AMM constant is claimed.
source: kps-s188 + agent Anscombe / THM-3601 continuation, 2026-08-21
audit: >
  ACCEPT.  An independent full-state reconstruction recovered the 4096-degree
  profile and fatal row, rederived the top-distance adjoint recurrence, and
  directly cancelled every free coefficient at cuts 382 and 383.  It obtained
  the complete negative-cut interval, boundary integers, digest, bit lengths,
  active-cell count, and exact radical gaps, and checked the departure-row
  indexing and scope.  Normal, optimized, and stored transcripts agree.
depends_on:
  - THM-3601-amm-r2048-terminal-failure-adjoint-horizon
related:
  - THM-3597-amm-dyadic-rule-a-transition-and-adjoint-horizons-through-R1024
  - THM-3588-amm-r512-truncated-adjoint-pascal-repair-horizons
script: 04-computation/amm_binary_sturmian_R4096_adjoint_horizon_thm3616.py
output: 05-knowledge/results/amm_binary_sturmian_R4096_adjoint_horizon_thm3616.out
script_sha256: c6753be05870700747aadb3f43e41a24e80f84f37b79e07fe70960f71a9e864a
output_sha256: 519e52f14ad8dc8cb42a182d769b00847eed1a19783958184adf36f87c8d679f
semantic_sha256: 4cd0b160559ea7b6820b90f7ccb204d76206d3271517b6e6afa06e5dce6a5d9a
hash_basis: LF-normalized bytes
---

# THM-3616 -- AMM R=4096 adjoint horizon and golden finite-scale obstruction

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  This executes the hostile test pre-registered in THM-3601.  The
test confirms the next exact adjoint horizon but rejects the simplest reading
of the apparent golden scaling.

## 1. Exact fatal row and adjoint wall

For `0<=i<4096`, use the inherited exact Fibonacci--Lucas comparison to set

```text
d_i(D0)=floor(log_5(phi^(2(4096+i))))+D0.                (1)
```

On the archived last-failing Rule-A offset `D0=88`, the exact trace has

```text
first death j=1014,
(d_0,d_1014,d_4095)=(2537,3143,4986).                   (2)
```

Use the same halved junk states, Lucas bounds, and Pascal transitions as
THM-3601.  If an admissible alternative agrees with this Rule-A trace through
row `s-1`, transposed propagation of the fatal-top evaluation gives a
nonnegative integer multiplier certificate and a necessary inequality

```text
0 <= B_s.                                                (3)
```

The complete exact sweep proves

```text
B_382>0>B_383,                                           (4)
```

so agreement through row 382 would force the false inequality `0<=B_383`.
Every admissible continuation that avoids this fatal inequality must depart
from Rule A at or before row

```text
q=382.                                                   (5)
```

The two boundary integers have signs and bit lengths

```text
B_382: positive, 2452 bits,
B_383: negative, 2448 bits,                              (6)
```

with compact boundary digest

```text
b92b6be9ea381a4a1789f02765a561bd93b115235351aa60432ad9e50d6cc728, (7)
```

and the first negative cut has `199396` active multiplier cells.

## 2. Compressed exact recurrence

The computation is not a dense matrix power.  Starting from the fatal row,
encode the top-distance multipliers by a polynomial `M_i(z)`.  With

```text
M_1013(z)=1,                                             (8)
```

the exact backward step is

```text
flat degree step:    M_(i-1)=(1+z)M_i,
rising degree step:  M_(i-1)=Pi_(>=0)[z^(-1)(1+z)^2 M_i]. (9)
```

Here `Pi_(>=0)` discards negative powers only.  The initial 1366 top-distance
cells cover the entire hostile horizon because `1365>1014`.  A shrinking
forward pass pairs the resulting nonnegative multipliers with the exact clamp
charges, and one suffix sum produces every cut value.  As an independent
indexing control, the same implementation reproduces all THM-3601 `R=2048`
boundary bits, digest, and active metadata before computing `(2)--(7)`.

## 3. The golden finite-scale hypothesis fails its hostile test

The three distinct normalized coordinates at `R=4096` are

```text
j/R = 507/2048 = 0.24755859375,
q/j = 191/507  = 0.3767258382...,
q/R = 191/2048 = 0.09326171875.                         (10)
```

THM-3601 nominated

```text
theta_gold=phi^(-2)/4=(3-sqrt(5))/8
          =0.09549150281252628795...                    (11)
```

as a possible continuum adjoint zero.  At the new scale the discrepancy is
exactly

```text
q/R-theta_gold=(-577+256sqrt(5))/2048
              =-0.00222978406252628795... <0.           (12)
```

The sign is rigorous because

```text
577^2 > 5*256^2.                                        (13)
```

At `R=2048` the corresponding discrepancy was

```text
(-573+256sqrt(5))/2048<0.                               (14)
```

Hence the absolute-error increase is not merely numerical:

```text
|error_4096|-|error_2048|=4/2048=1/512>0.               (15)
```

The magnitudes are approximately `0.0002766591` and `0.0022297841`.  Thus the
pre-registered next epoch moves exactly `1/512` farther from `(11)`.  This
refutes monotone finite-scale sharpening toward `theta_gold` and the simplest
identification of the observed ratios with a scale-free golden coordinate.

It does not refute eventual convergence along a subsequence or after adding a
missing coordinate.  The exact likely sidecars are the normalized offset
`D0/R`, the local Sturmian phase of the degree word, and a clamp-saturation
statistic; none is yet proved to control the horizon.

## 4. Scope and reproduction

Reproduce with

```bash
python3 04-computation/amm_binary_sturmian_R4096_adjoint_horizon_thm3616.py
python3 -O 04-computation/amm_binary_sturmian_R4096_adjoint_horizon_thm3616.py
```

This theorem certifies the archived failing trace `R=4096,D0=88`.  It does
not expensively replay the adjacent archived survivor `D0=89`.  As in
THM-3601, the adjoint wall constrains only continuations agreeing with the
failed Rule-A trace up to a cut: it does not construct a repair, prove the
full entry polytope infeasible, or establish entry into the superblock basin.
It proves no alternative feasibility, asymptotic horizon law, uniform
extractor, value of the AMM constant `C*`, or improvement to its bound.

The independent audit reconstructed the degree word and fatal row with a
full-state, feed-aware replay, then derived the compressed recurrence anew.
Direct coefficient ledgers at both boundary cuts leave zero residue in every
free coordinate.  It found negative cuts exactly `383,...,1013`, reproduced
the boundary digest and `199396` active cells, and checked that cut `383`
assumes agreement through row `382`.  Finally it verified
`577^2-5*256^2=5249` and `573^2-5*256^2=649`, so the two signed errors and
their absolute difference `1/512` are exact.

**QED.**

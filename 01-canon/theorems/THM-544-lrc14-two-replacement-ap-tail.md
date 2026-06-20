---
id: THM-544
title: LRC14 two-replacement AP-tail theorem - every three-hole two-tail AP row is at least the AP one-hole second value
status: PROVED (two-comb and one-comb rational cutoffs plus finite exact interval certificate)
source: codex-2026-06-19-S35
depends_on:
  - THM-541
  - THM-543
  - HYP-2654
related:
  - HYP-2651
  - HYP-2653
  - HYP-2650
  - HYP-2652
  - HYP-2569
external: Lonely Runner Conjecture, n=14
---

# THM-544 - Two-Replacement AP-Tail Theorem

Let

```text
C_{a,b,c,r,s} = ({1,2,...,13} \ {a,b,c}) union {r,s},
```

where `1 <= a < b < c <= 13` and `14 <= r < s` are integers.  Define

```text
G_{a,b,c,r,s} = {t in [0,1): ||d t|| > 1/14 for every d in C_{a,b,c,r,s}}.
```

Then

```text
meas(G_{a,b,c,r,s}) >= 426/35035.
```

Thus no two-replacement AP-tail row lies below the AP one-hole second value.

The exact minimum found in the finite residue is

```text
holes=(4,6,10), tails=(14,15),
C=(1,2,3,5,7,8,9,11,12,13,14,15),
meas(G_C)=14249/252252.
```

Its exact gap above the AP one-hole second value is

```text
14249/252252 - 426/35035 = 1141/25740.
```

## Certificate

Script:

```text
04-computation/lrc14_two_replacement_ap_tail_theorem_codex_s35.py
```

Stored output:

```text
05-knowledge/results/lrc14_two_replacement_ap_tail_theorem_codex_s35.out
```

The script uses exact rational interval arithmetic and prints exact fractions
only.

## Two-Comb Cutoff

For a three-hole base

```text
B_{a,b,c} = {1,...,13} \ {a,b,c},
```

let `G_B` have measure `M_B` and `c_B` interval components.  For two replacement
danger combs `D_r` and `D_s`, the union bound and exact one-comb component
estimate give

```text
meas(G_B \ (D_r union D_s))
  >= (5/7)M_B - 2c_B/(7r) - 2c_B/(7s).
```

Every three-hole base has positive exact two-tail slack:

```text
5M_B - 7*(426/35035) > 0.
```

Therefore all pairs with both tails at least

```text
R_B = ceil(4c_B/(5M_B - 7*(426/35035)))
```

are certified above the threshold.  The largest two-comb cutoff is

```text
holes=(4,5,6), M=79/1155, c=16, 5M-7Q=551/2145, R=250.
```

## Fixed-Tail Cutoff

For each remaining smaller tail `14 <= r < R_B`, the script computes the exact
safe set after adding `r`.  If this row has measure `M_{B,r}` and component
count `c_{B,r}`, then for the second tail

```text
meas(G_{B,r} \ D_s) >= (6/7)M_{B,r} - 2c_{B,r}/(7s).
```

Every fixed-tail row has positive exact slack

```text
6M_{B,r} - 7*(426/35035) > 0.
```

The largest active fixed-tail cutoff is

```text
holes=(4,5,6), r=14,
M=79/1155, c=16, 6M-7Q=148/455, S=99.
```

## Finite Residue

The certificate covers:

```text
three-hole bases: 286
fixed smaller-tail rows checked for cutoff: 24824
finite two-tail pairs checked exactly: 129991
rows below 426/35035: 0
```

The exact finite minimum is the row displayed above:

```text
holes=(4,6,10), tails=(14,15),
safe=14249/252252,
gap_above_Q=1141/25740.
```

The four old drop-6 mouth intervals survive in that minimum with total old
survivor mass `7/858`, but the theorem does not need mouth retention: every
row in the layer is already above the AP-second threshold.

## Proof Use

THM-543 showed that the one-replacement AP-tail layer has exactly one
below-second row, the mouth-retaining `(6,10)->20` tail.  THM-544 shows that the
next AP-tail layer, with three AP holes and two replacement tails, has no
below-second rows at all.

This narrows HYP-2654 further.  The bounded AP-tail obstruction is not in the
first two replacement layers.  Any remaining below-second bounded AP-tail row
must start at three or more replacement tails, or else the proof must leave the
AP-tail model and use state-word damage, additive/Freiman structure, or the
far-element discrepancy branch.

## Tournament Analysis

Vertices: two-tail proof gates.

Pairwise observable: exact eliminations before the threshold `Q=426/35035`.

Switch/gauge: two-comb cutoff, then fixed-tail one-comb cutoff, then exact
finite scan.

Hamiltonian path:

```text
base_two_tail_slack > fixed_tail_slack > finite_scan > best_gap > mouth_survivor
```

Fingerprint: transitive proof-obligation order; directed `3`-cycles `0`.

Assumption challenge: the useful tournament vertices are proof gates, not
runners or arcs.  The quotient preserves the exact predicate
`meas(G_C) >= 426/35035` for the whole two-replacement layer.

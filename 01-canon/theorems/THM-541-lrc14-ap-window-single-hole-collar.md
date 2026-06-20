---
id: THM-541
title: LRC14 AP-window single-hole collar - among {1,...,13} with one deletion, the drop-6 core uniquely minimizes the level-1/14 fixed-observer safe measure
status: PROVED (finite exact addressed wall certificate; strict rational comparisons)
source: codex-2026-06-19-S34
depends_on:
  - THM-523
  - HYP-2651
related:
  - HYP-2650
  - HYP-2652
  - HYP-2648
  - THM-526
external: Lonely Runner Conjecture, n=14
---

# THM-541 - AP-window single-hole collar

For `e in {1,...,13}`, let

```text
C_e = {1,2,...,13} \ {e}
```

and define the fixed-observer level-`1/14` safe set

```text
G_e = {t in [0,1): ||c t|| > 1/14 for every c in C_e}.
```

Then

```text
meas(G_e) >= 7/858
```

for every `e`, with equality if and only if `e=6`.  The next value is

```text
426/35035
```

at `e=12`, separated from the minimum by

```text
426/35035 - 7/858 = 841/210210.
```

## Certificate

Script:

```text
04-computation/lrc14_ap_window_single_hole_certificate_codex_s34.py
```

Stored output:

```text
05-knowledge/results/lrc14_ap_window_single_hole_certificate_codex_s34.out
```

The script uses exact rational arithmetic.  For a speed `v`, the danger walls at
target `1/14` are

```text
L(v,a) = (14a - 1)/(14v),
R(v,a) = (14a + 1)/(14v),
```

with the wrapped left wall represented periodically as `L(v,v)`.  Each safe
component is an open gap from a right wall to a left wall,

```text
R(v,a) -> L(w,b),
```

and has length

```text
((14b - 1)v - (14a + 1)w)/(14vw).
```

The exact table is:

```text
rank | drop | safe measure | delta from 7/858
   1 |    6 |        7/858 |                0
   2 |   12 |    426/35035 |       841/210210
   3 |   10 |   1520/63063 |      2011/126126
   4 |    4 |      97/4004 |        193/12012
   5 |    2 |       11/364 |        265/12012
   6 |   13 |  6617/194040 |    65441/2522520
   7 |    5 |     103/2860 |         239/8580
   8 |    3 |       11/273 |         193/6006
   9 |    8 |    950/21021 |        519/14014
  10 |    9 |   3273/56056 |      8447/168168
  11 |   11 |   4733/76440 |     45203/840840
  12 |    1 |         1/14 |         190/3003
  13 |    7 |   4319/51480 |       3899/51480
```

The minimizing core is

```text
C_6 = (1,2,3,4,5,7,8,9,10,11,12,13),
```

with four addressed components:

```text
[29/182, 9/56]     len=1/728  R(13,2) -> L(12,2), det=3
[29/168, 27/154]   len=5/1848 R(12,2) -> L(11,2), det=5
[127/154,139/168]  len=5/1848 R(11,9) -> L(12,10), det=5
[47/56, 153/182]   len=1/728  R(12,10) -> L(13,11), det=3
```

Thus

```text
meas(G_6) = 2*(1/728) + 2*(5/1848) = 7/858.
```

The next competitor `e=12` has:

```text
[15/182, 13/154]   len=2/1001
[29/70, 41/98]     len=1/245
[57/98, 41/70]     len=1/245
[141/154, 167/182] len=2/1001
```

with total `426/35035`.

## Proof

For a fixed core `C_e`, form the finite union of closed danger arcs

```text
D_e = union_{c in C_e} {t : ||c t|| <= 1/14}.
```

The endpoints of `D_e` are exactly the wall points `L(v,a)` and `R(v,a)` above.
After sorting and merging the danger intervals, the components of `[0,1)\D_e`
are precisely the gaps from a right wall to a left wall that remain uncovered.
Endpoint inclusion is irrelevant for measure, so this gives the measure of the
strict set `G_e`.

The companion certificate carries out that merge for all `13` possible drops in
`Q`.  It verifies each component endpoint has the advertised owner, each
component midpoint is outside every danger arc, and the component list equals
the complement of the merged danger union.  Summing the listed rational lengths
gives the table above.  Every entry except `e=6` is strictly larger than
`7/858`, and the smallest strict difference is `841/210210`.  Hence the
minimum is unique and equal to `7/858`.  QED.

## Proof Use

This closes the first proof obligation isolated in HYP-2651:

```text
For C=[1,13]\{e}, meas(G_C) is minimized uniquely at e=6.
```

It does not prove the full OPEN-Q-108 uniform fattening lemma.  The next sharp
obligation is now the near-collar mouth-retention theorem:

```text
If a positive primitive 12-core C has meas(G_C) < 426/35035,
then C retains the drop-6 mouth template.
```

The theorem also explains why HYP-2650/HYP-2652's addressed-wall lesson is not
cosmetic.  The scalar measure is the final valuation, but the proof carrier is
the signed boundary address `R(v,a) -> L(w,b)`.  The drop-6 collar is the small
determinant pattern `[3,5,5,3]` between neighboring high-speed walls
`13,12,11,12,13`.

HYP-2654 records the first correction to the naive wording: the template is not
exact-row rigidity.  The AP-tail row `(1,2,3,4,5,7,8,9,11,12,13,20)` lies below
`426/35035`, but it preserves the four drop-6 mouth intervals and adds new
mouth mass `1/980`.

## Tournament Analysis

Vertices: the deleted AP-window positions `e=1..13`.

Pairwise observable: lower exact safe measure wins.

Switch/gauge: retain `R/L` boundary-owner addresses before scalarizing to
measure.

Hamiltonian path:

```text
drop-6 > drop-12 > drop-10 > drop-4 > drop-2 > drop-13 > drop-5
> drop-3 > drop-8 > drop-9 > drop-11 > drop-1 > drop-7
```

Fingerprint: strict total order, score histogram
`{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1}`,
and no directed `3`-cycles.

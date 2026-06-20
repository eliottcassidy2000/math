---
id: THM-543
title: LRC14 one-replacement AP-tail theorem - the full two-hole one-tail AP layer has only the (6,10)->20 below-second row
status: PROVED (periodic-comb rational cutoff plus finite exact interval certificate)
source: codex-2026-06-19-S35
depends_on:
  - THM-541
  - THM-542
  - HYP-2654
related:
  - THM-544
  - HYP-2651
  - HYP-2655
  - HYP-2653
  - HYP-2650
  - HYP-2652
  - HYP-2569
external: Lonely Runner Conjecture, n=14
---

# THM-543 - One-Replacement AP-Tail Theorem

Let

```text
C_{a,b,r} = ({1,2,...,13} \ {a,b}) union {r},
```

where `1 <= a < b <= 13` and `r >= 14` is an integer.  Define

```text
G_{a,b,r} = {t in [0,1): ||c t|| > 1/14 for every c in C_{a,b,r}}.
```

Then

```text
meas(G_{a,b,r}) < 426/35035
```

if and only if

```text
(a,b,r) = (6,10,20).
```

In that exceptional case

```text
C_{6,10,20} = (1,2,3,4,5,7,8,9,11,12,13,20)
```

and

```text
meas(G_{6,10,20}) = 3859/420420 = 7/858 + 1/980.
```

The four old drop-6 mouth intervals from THM-541 survive intact:

```text
[29/182, 9/56]
[29/168, 27/154]
[127/154, 139/168]
[47/56, 153/182].
```

## Certificate

Script:

```text
04-computation/lrc14_one_replacement_ap_tail_theorem_codex_s35.py
```

Stored output:

```text
05-knowledge/results/lrc14_one_replacement_ap_tail_theorem_codex_s35.out
```

The script uses exact rational interval arithmetic throughout.  The proof
comparisons and printed certificate use `Fraction` values, not floating-point
approximations.

## Periodic-Comb Cutoff

For a two-hole base

```text
B_{a,b} = {1,...,13} \ {a,b},
```

let `G_{a,b}` be the level-`1/14` safe set, with measure `M_{a,b}` and
`c_{a,b}` interval components.  Let `D_r` be the danger comb for speed `r`.
For every interval component of `G_{a,b}`, full `1/r` periods remove exactly a
`1/7` share, while the two endpoint partial periods remove at most two danger
teeth of total length `2/(7r)`.  Therefore

```text
meas(G_{a,b} \ D_r) >= (6/7)M_{a,b} - 2c_{a,b}/(7r).
```

For every pair `1 <= a < b <= 13`, the exact denominator slack

```text
6M_{a,b} - 7*(426/35035)
```

is positive.  Thus all

```text
r >= R_{a,b} = ceil(2c_{a,b}/(6M_{a,b}-7*(426/35035)))
```

are certified at or above `426/35035`.

The weakest slack and largest cutoff occur at the same resonant base:

```text
holes=(6,10), M=313/9702, c=8,
6M-7Q=11399/105105, R=148.
```

The full exact cutoff table for all `78` pairs is in the stored output.

## Finite Residue

After the periodic-comb cutoff, only `3277` exact rows with `14 <= r < R_{a,b}`
remain.  The finite certificate finds a single row below `426/35035`:

```text
holes=(6,10), r=20,
safe=3859/420420,
old_survivor=7/858,
new_mass=1/980.
```

All assertions are made in the script as exact `Fraction` equalities and
inequalities.

## Proof Use

THM-542 closed the one-tail subcase where one removed AP speed was already `6`.
THM-543 closes the entire one-replacement AP-tail layer.  No two-hole AP base
without the central `6` hole can fall below the AP one-hole second value after a
single replacement tail, and the only below-second central-hole row is the
already-known `(6,10)->20` row.

This sharpens HYP-2654: in the one-replacement AP-tail layer, below-second means
old drop-6 mouth retention exactly.  The remaining near-collar problem is no
longer a one-tail problem; it must enter multi-hole/multi-tail templates,
state-word damage, or the far-element plateau/discrepancy branch.

THM-544 closes the next layer: every three-hole/two-replacement AP-tail row is
already at least `426/35035`, with exact finite minimum
`50189/3223220=426/35035+1571/460460`.

The concurrent HYP-2653/HYP-2655 exact decorrelation work shows why this split
matters: a small uniform far-tail constant fails on resonant multiscale cores,
so bounded AP-tail layers should be certified by exact rational cutoffs, while
genuinely wide branches need the joint plateau/Delta recursion of HYP-2655.

## Tournament Analysis

Vertices: proof-obligation gates for the one-replacement AP-tail layer.

Pairwise observable: exact row eliminations before the AP-second threshold.

Switch/gauge: use the periodic-comb denominator cutoff before finite wall scan.

Hamiltonian path:

```text
denominator_slack > comb_cutoff > finite_exact_scan > mouth_survivor > relation_rank
```

Fingerprint: transitive proof-obligation order; directed `3`-cycles `0`.

Assumption challenge: tournament vertices do not need to be runners or arcs.  In
this proof the useful vertices are gates in the proof quotient, while the
preserved LRC predicate is exact safe measure below `426/35035`.

---
id: HYP-2654
title: LRC14 near-collar AP-tail template - below the AP one-hole second value, retain the drop-6 mouth geometry
status: OPEN; structured AP-tail evidence
source: codex-2026-06-19-S34
depends_on:
  - THM-541
  - THM-542
  - THM-543
  - THM-544
  - HYP-2651
related:
  - HYP-2659
  - HYP-2655
  - HYP-2653
  - HYP-2650
  - HYP-2652
  - HYP-2648
  - HYP-2569
  - OPEN-Q-108
---

# HYP-2654 - Near-Collar AP-Tail Template

## Claim

The next proof obligation after THM-541 is not:

```text
meas(G_C) < 426/35035 forces C = {1,...,13}\{6}.
```

That exact-row statement is already too strong in the AP-tail pocket.  The
correct target is a template statement:

```text
meas(G_C) < 426/35035 forces the drop-6 mouth geometry to survive.
```

Concretely, the four drop-6 safe components

```text
[29/182, 9/56]
[29/168, 27/154]
[127/154, 139/168]
[47/56, 153/182]
```

should either survive as a full sub-reservoir of `G_C`, or the row should pay at
least the AP one-hole second value `426/35035`.

## Evidence

Script:

```text
04-computation/lrc14_near_collar_tail_template_codex_s34.py
```

Stored output:

```text
05-knowledge/results/lrc14_near_collar_tail_template_codex_s34.out
```

The scout enumerates structured AP-tail cores

```text
({1,...,13} \ holes) union replacements,
```

with `len(replacements)=len(holes)-1`, `max(replacement)<=35`, and up to three
AP holes.  It scanned `67,795` primitive rows.

Rows below the AP one-hole second value:

```text
rank | safe          | holes   | repl | old drop-6 survivor | new mass
   1 | 7/858         | (6,)    | ()   | 7/858               | 0
   2 | 3859/420420   | (6,10)  | (20) | 7/858               | 1/980
```

Thus every below-threshold AP-tail row in this pocket:

1. has `6` among the AP holes;
2. leaves the old drop-6 mouth mass undamaged (`old_survivor=7/858`);
3. is the drop-6 collar plus optional new mouth mass.

The best row below the cutoff but not equal to the exact drop-6 core is

```text
(1,2,3,4,5,7,8,9,11,12,13,20)
```

with holes `(6,10)`, replacement `20`, and

```text
meas(G_C) = 3859/420420 = 7/858 + 1/980.
```

The first row at the AP one-hole second value is still the endpoint drop:

```text
(1,2,3,4,5,6,7,8,9,10,11,13)
```

with `426/35035`, and it has no old drop-6 survivor mass.

## First Two AP-Tail Replacement Layers Closed

THM-542 first proved the infinite one-tail AP-drop-6 subcase exactly.  For

```text
C_{h,r} = ({1,...,13} \ {6,h}) union {r},  h != 6, r >= 14,
```

the only row below `426/35035` is `(h,r)=(10,20)`, with

```text
meas(G_C) = 3859/420420 = 7/858 + 1/980.
```

The proof uses a rational periodic-comb cutoff:

```text
meas(G_h \ D_r) >= (6/7)meas(G_h) - 2c_h/(7r),
```

then checks only `863` exact finite rows below the resulting cutoffs.  This
turns the first AP-tail mouth-retention statement from evidence into a theorem.

THM-543 then closes the full one-replacement AP-tail layer.  For every

```text
C_{a,b,r} = ({1,...,13} \ {a,b}) union {r},
1 <= a < b <= 13, r >= 14,
```

the exact inequality

```text
meas(G_C) < 426/35035
```

holds if and only if `(a,b,r)=(6,10,20)`.  The certificate checks only `3277`
finite rows after exact rational comb cutoffs over all `78` two-hole bases; the
weakest slack and largest cutoff are both the resonant base `(6,10)`:

```text
M=313/9702, c=8, 6M-7Q=11399/105105, R=148.
```

Thus the one-replacement AP-tail layer cannot create a below-second row without
retaining the four old drop-6 mouth intervals exactly.

THM-544 closes the next AP-tail layer.  For every

```text
C_{a,b,c,r,s} = ({1,...,13} \ {a,b,c}) union {r,s},
1 <= a < b < c <= 13, 14 <= r < s,
```

one has

```text
meas(G_C) >= 426/35035.
```

The exact proof uses a two-comb cutoff over all `286` three-hole bases, then a
fixed-tail one-comb cutoff over `24,824` exact smaller-tail rows, leaving
`400,090` finite two-tail pairs.  None is below threshold.  The finite minimum
is

```text
holes=(4,6,10), tails=(20,46),
meas(G_C)=50189/3223220,
old_survivor=1/364,
meas(G_C)-426/35035=1571/460460.
```

So the first layer has one harmless mouth-retaining below-second exception, and
the second replacement layer has no below-second exceptions at all.

## Refined Proof Target

The HYP-2651 near-collar target should be replaced by:

```text
If C is a primitive positive 12-core and meas(G_C) < 426/35035,
then G_C contains the four drop-6 mouth intervals, up to scale-preserving
fixed-observer address equivalence.
```

An even more proof-facing split is:

```text
either old drop-6 mouth damage > 0 or no central drop-6 hole
    => meas(G_C) >= 426/35035,
else
    meas(G_C) >= 7/858 from the surviving mouth reservoir.
```

This is exactly the addressed-wall lesson of HYP-2650/HYP-2652.  The scalar
cutoff does not identify one row; it identifies a boundary-owner template.

## Tournament Analysis

Vertices are proof-obligation observables for the AP-tail near-collar rows:

```text
drop6_mouth_retained
old_mouth_damage
new_mouth_mass
hole_count
sumset_excess
raw_tail_size
```

Pairwise observable: which observable separates rows below `426/35035` from the
rest.

Switch/gauge: classify by retained drop-6 mouth mass before scalar safe measure.

Hamiltonian path:

```text
drop6_mouth_retained
> old_mouth_damage
> new_mouth_mass
> hole_count
> sumset_excess
> raw_tail_size
```

Fingerprint in the stored pocket: below-threshold rows with undamaged old mouths
`2/2`; below-threshold rows with hole `6`: `2/2`; directed priority tournament
transitive.

## Scope

This is not a proof of the full near-collar theorem.  The broad exact scan
through `[1,22]` was attempted during this session but timed out before
completion, so HYP-2654 is structured AP-tail evidence, not a bounded-box
theorem.  Its value is that it corrects the next target and now removes the
first two replacement AP-tail layers: prove mouth-retention rigidity, not
exact-row rigidity, for deeper multi-tail/state-word damage and
far/discrepancy branches.

Namespace note: this packet was renumbered after origin/main landed KPS HYP-2653
for the far-element decorrelation rate.  The later HYP-2655 update refutes the
single small uniform-constant dovetail and re-localizes the wide branch to a
joint plateau/Delta recursion.  These results are complementary: HYP-2655
supports the wide/far branch, while HYP-2654 isolates the bounded AP-tail
mouth-retention branch.  THM-543 and THM-544 follow the complementary strategy
of certifying bounded AP-tail cores by exact rational finite cutoffs.

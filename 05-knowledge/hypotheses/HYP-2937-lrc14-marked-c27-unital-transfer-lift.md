---
id: HYP-2937
title: Marked C27 shell transfers lift into q=3 unital blocks only after AP/Goddyn-Wong labels are fixed
status: COMPUTATIONAL SCOUT / positive incidence carrier with cyclic-carry guardrail
source: codex-2026-06-24
tags: [lrc14, c27, unital, goddyn-wong, ap, marked-transfer, pair-incidence]
related:
  - THM-407
  - THM-417
  - HYP-2891
  - HYP-2894
  - HYP-2920
results:
  - 04-computation/lrc14_hyp2937_c27_unital_lift_codex.py
  - 05-knowledge/results/lrc14_hyp2937_c27_unital_lift_codex.out
---

# HYP-2937: marked C27 transfers into q=3 unital blocks

This scout tests the prompt:

```text
lift marked C27 transfers into q=3 unital 4-point blocks after AP/Goddyn-Wong
labels are attached.
```

The computation is in
[lrc14_hyp2937_c27_unital_lift_codex.py](/home/bigo/math/04-computation/lrc14_hyp2937_c27_unital_lift_codex.py:1),
with output at
[lrc14_hyp2937_c27_unital_lift_codex.out](/home/bigo/math/05-knowledge/results/lrc14_hyp2937_c27_unital_lift_codex.out:1).

## Construction

The script builds the classical Hermitian unital of order `q=3`:

```text
points = q^3+1 = 28
blocks = 63
block size = q+1 = 4
lambda = 1
replication = 9
```

It verifies the design exactly:

```text
block_size_hist = {4: 63}
point_rep_hist  = {9: 28}
pair_rep_hist   = {1: 378}
```

The 27 affine points are labelled by residues in `C27` via base-3 digits:

```text
r = d0 + 3*d1 + 9*d2  ->  (x,y) in GF(9)^2,  y^3+y = x^4.
```

The last point is labelled `inf`.

This is the honest lift:

```text
C27 residue labels -> F3^3 affine chart of the unital.
```

It is not a cyclic-action lift.  That distinction is the main guardrail.

## AP/Goddyn-Wong labels

At `n=14`, the pair-sum modulus is:

```text
C = 2n-1 = 27.
```

The 13 antipodal shell representatives are `1,...,13`.  THM-407 collapses them
to three `3`-adic strata:

```text
gcd(shell,27)=1  -> 9 shells
gcd(shell,27)=3  -> 3 shells
gcd(shell,27)=9  -> 1 shell
```

The AP row has all shells filled once.  The Goddyn-Wong row replaces `12` by
`24`, which is `-3 mod 27`.  So the AP-to-GW marked transfer is:

```text
shell 12 -> shell 3
hole shell12 -> collision shell3.
```

This is also the shell-doubling move:

```text
12 -> fold(2*12 mod 27) = fold(24) = 3.
```

## Transfer Lift

For each shell doubling transfer:

```text
a -> fold(2a mod 27),
```

the script asks which q=3 unital blocks carry the source-target point-pairs.
Because the unital has `lambda=1`, every point-pair lands in a unique block.

The shell transfer table includes:

```text
3  -> 6   gcd3 -> gcd3   4 pairs, 4 carrier blocks
6  -> 12  gcd3 -> gcd3   4 pairs, 4 carrier blocks
12 -> 3   gcd3 -> gcd3   4 pairs, 4 carrier blocks  (GW petal)
9  -> 9   gcd9 -> gcd9   2 pairs, 1 carrier block
```

The GW transfer is carried by four unital blocks:

```text
{3, 5, 15, 17}
{3, 12, 21, inf}
{6, 15, 24, inf}
{12, 14, 24, 26}
```

These are the four blocks containing the four cross-pairs between:

```text
hole shell12      = {12,15}
collision shell3  = {3,24}.
```

## Interpretation

Positive result:

```text
The q=3 unital is a real pair-incidence carrier for marked C27 transfers.
```

It gives a canonical `lambda=1` way to ask which four-point packet contains a
given AP/Goddyn-Wong transfer pair.

Negative guardrail:

```text
The natural unital chart is F3^3, not cyclic C27.
```

The base-3 labelling erases cyclic carry.  Therefore the block statistics are
label-dependent unless a stronger `C27`-compatible marking is proved.  This
agrees with HYP-2894: the unital should not be promoted into a standalone LRC
invariant.

## Proof Use

The live use is:

```text
first attach AP/Goddyn-Wong shell labels,
then use q=3 unital blocks as a pair-incidence ledger.
```

The unital can organize marked transfers after the LRC meaning has been fixed.
It cannot by itself decide tightness, because the tightness data lives in the
marked Farey node and in cyclic/magnitude information that the `F3^3` unital
chart does not preserve.

Candidate next computation:

```text
For each residual row near the AP/GW frontier, push its C27 shell defects into
the four unital carrier blocks for the relevant transfer, then test whether
the signed leakage concentrates in the two inf-blocks or the two finite blocks.
```


---
id: HYP-2942
title: LRC14 C27 unital block-lift obstruction and calibrated pair-completion
status: PROOF-INTERFACE / global no-go plus branch-local calibrated positive lift; not a proof
source: codex-2026-06-24-S140
related:
  - HYP-2941
  - HYP-2940
  - HYP-2939
  - HYP-2938
  - HYP-2937
  - HYP-2936
  - HYP-2935
  - HYP-2932
  - HYP-2894
  - HYP-2892
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_c27_unital_block_lift_codex_s140.py
  - 05-knowledge/results/lrc14_c27_unital_block_lift_codex_s140.out
---

# HYP-2942: The q=3 unital is a branch-local C27 chart and AP/GW pair-completion forum

S140 tests the prompt:

```text
lift HYP-2937 marked C27 transfers into q=3 unital 4-point blocks after
AP/Goddyn-Wong labels are attached.
```

The synthesized answer has two layers:

```text
global raw lift of all low-frontier transfers: no
branch-local lift and calibrated AP/GW pair-completion: yes
```

The negative result is a useful obstruction, not a failure of the analogy.  The
positive result keeps the unital as a precise relational forum for transfer
packets.

## Computation

The script
`04-computation/lrc14_c27_unital_block_lift_codex_s140.py` stores output at
`05-knowledge/results/lrc14_c27_unital_block_lift_codex_s140.out`.

It constructs the Hermitian unital over

```text
GF(9)=F3[w]/(w^2+1)
```

using

```text
x^4 + y^4 + z^4 = 0
```

in `PG(2,9)`.  The sanity checks are:

```text
points = 28
non-tangent 4-point blocks = 63
line intersection histogram = {1:28, 4:63}
point replication = 9
every pair of points lies in exactly one block
block intersection histogram = {0:945, 1:1008}
```

Thus two distinct q=3 unital blocks never share two points.

## Raw C27 Residue-Pair Lift

The first lift sends a marked transfer to the antipodal C27 residue packet:

```text
H[a] -> D[d]  maps to  {a, 27-a, d, 27-d}.
```

The low-frontier packets become:

```text
GW H12->D3          {3, 12, 15, 24}
K33 H12->D9         {9, 12, 15, 18}
petal H10->D7       {7, 10, 17, 20}
petal H13->D1       {1, 13, 14, 26}
```

The global low-frontier raw lift fails because

```text
GW block  cap  K33 block = {12,15}.
```

Since the q=3 unital is a `2-(28,4,1)` design, the same pair cannot lie in two
different blocks.  A single raw unital chart cannot contain both the GW
`H12->D3` block and the K33 `H12->D9` near-miss block.

## AP/GW Global-Label No-Go

The tempting labelled model

```text
H[a] -> D[d]  maps to  {AP, GW, H_a, D_d}
```

also fails globally.  Every transfer block contains `{AP,GW}`, so more than one
such block repeats a pair.  AP/GW cannot be universal vertices repeated inside
every transfer block.  They must be a calibrated anchor, branch colors, or
external metadata.

## Positive Branch-Local Lift

The pair-repeat obstruction kills a single atlas, not local charts.  The
following desired block systems have no repeated C27 pair and embed into the
q=3 unital incidence design:

```text
tight chart:
  GW H12->D3
  petal H10->D7
  petal H13->D1

K33 chart:
  K33 H12->D9
  petal H10->D7
  petal H13->D1
```

The S138 two-hole frontier rows lift as two-block splices:

```text
drop(10,12)->add(20,24)
  = petal H10->D7  plus  GW H12->D3

drop(10,12)->add(20,36)
  = petal H10->D7  plus  K33 H12->D9
```

This matches the HYP-2940 reading: the two-hole rows are splices of the unit
petal with one of the two `12`-branch packets.

## Calibrated AP/GW Pair-Completion Lift

The second positive layer attaches the 28 labels

```text
AP, GW, H1..H13, D1..D13
```

to one Hermitian unital and calibrates one block as

```text
{AP, GW, H12, D3}.
```

Then the Goddyn-Wong transfer `H12 -> D3` is literally the AP/GW anchor block.
The other transfer pairs still have unique q=3 completion blocks.  The script
records one deterministic calibration:

```text
GW H12->D3:
  {AP, D3, GW, H12}

near-miss H12->D9:
  {D10, D9, H12, H9}

petal H10->D7:
  {D7, D9, H10, H13}

petal H13->D1:
  {D1, D11, D12, H13}
```

The two S138 splices are separated by incidence:

```text
10+12 -> 20+24:
  petal block plus AP/GW block, disjoint, union size 8

10+12 -> 20+36:
  petal block plus near/K33 block, sharing D9, union size 7
```

The linked packet is:

```text
AP/GW --H12-- near/K33 --D9-- petal10.
```

## Interpretation

The proof-use distinction is now sharp:

```text
raw residue lift: detects forbidden global superposition of the two H12 branches
calibrated lift: supplies unique 4-point completions for labelled transfer pairs
```

This gives a concrete packet grammar:

```text
low-gap residual
  -> AP/GW anchor block
  -> unit petal block
  -> AP/GW--near/K33--petal chain
```

The first route is the known tight floor, the second routes to the C27
petal/two-swap discharge, and the third is the candidate HYP-2908 / THM-572
state-lift packet.

## Guardrail

This is not a canonical AP8 pair-slot design and not scalar evidence for
LRC14.  It is a branch-local, AP/GW-calibrated, pair-unique incidence forum.
Any future proof that tries to put both `12` branches into one unital object
must either split the H12 pair by an added branch label or explicitly move to a
multi-chart atlas.

## Tournament Analysis

S140 uses proof obligations and lifted packets as tournament vertices, not
runners.  The carrier observable is:

```text
theorem-scale retention,
C27 predicate retention,
lambda=1 incidence,
finite checkability,
anti-scalar guard.
```

The carrier tournament is transitive:

```text
exact M/Farey branch
> C27 marked transfer
> unital pair-repeat obstruction
> branch-local q3 unital chart
> calibrated AP/GW pair-completion
> S138 two-block splice path
> global AP/GW vertex model
> raw unital analogy
```

The calibrated packet tournament uses:

```text
AP/GW visibility,
unit/nonunit content,
component depth,
K33 flag,
smaller unital union as pressure/tie gauge.
```

It is also transitive in the recorded calibration, with the linked and disjoint
two-block splices appearing at the top of the packet comparison.

## POKE Forum Contribution

Treat q=3 unital blocks as a branch-local pair-completion forum for C27
transfer packets.  The most actionable new rule is:

```text
do not merge GW H12->D3 and K33 H12->D9 in one raw C27 unital chart unless
the H12 pair is split by additional structure.
```

The complementary constructive rule is:

```text
after AP/GW calibration, inspect the unique completion block for each marked
transfer pair and compare block intersections before attempting a state lift.
```

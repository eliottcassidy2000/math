---
id: HYP-2942
title: LRC14 C27 unital block-lift obstruction and branch-local charts
status: PROOF-INTERFACE / global no-go plus branch-local positive lift; not a proof
source: codex-2026-06-24-S140
related:
  - HYP-2941
  - HYP-2940
  - HYP-2939
  - HYP-2938
  - HYP-2937
  - HYP-2936
  - HYP-2894
  - HYP-2892
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_c27_unital_block_lift_codex_s140.py
  - 05-knowledge/results/lrc14_c27_unital_block_lift_codex_s140.out
---

# HYP-2942: The q=3 unital is a branch-local C27 chart, not a global atlas

S140 tests the concrete prompt:

```text
Can HYP-2937 marked C27 transfers lift into q=3 unital 4-point blocks after
AP/Goddyn-Wong labels are attached?
```

The answer is split:

```text
global lift of all low-frontier transfers: no
branch-local lift and S138 two-block splice lift: yes
```

The useful obstruction is not numerical.  It is the `lambda=1` pair uniqueness
of the `2-(28,4,1)` q=3 unital.

## Computation

The script
`04-computation/lrc14_c27_unital_block_lift_codex_s140.py` stores output at
`05-knowledge/results/lrc14_c27_unital_block_lift_codex_s140.out`.

It constructs the Hermitian q=3 unital through the S105 code and verifies:

```text
points = 28
blocks = 63
block size = 4
block intersection histogram = {0: 945, 1: 1008}
```

Thus two distinct unital blocks never share two points.

## Raw C27 Residue-Pair Lift

The most natural C27 transfer-to-block map is:

```text
H[a] -> D[d]  maps to  {a, 27-a, d, 27-d}.
```

For the low-frontier transfer packets this gives:

```text
GW H12->D3          {3, 12, 15, 24}
K33 H12->D9         {9, 12, 15, 18}
petal H10->D7       {7, 10, 17, 20}
petal H13->D1       {1, 13, 14, 26}
```

The global low-frontier lift fails immediately:

```text
GW block  cap  K33 block = {12,15}.
```

Since a q=3 unital is a `2-(28,4,1)` design, the same pair cannot lie in two
different blocks.  Therefore one q=3 unital chart cannot simultaneously contain
both the GW `H12->D3` block and the K33 near-miss `H12->D9` block under this
raw residue-pair interpretation.

## AP/GW Global-Label Model

The tempting labelled model

```text
H[a] -> D[d]  maps to  {AP, GW, H_a, D_d}
```

is even worse.  Every transfer block contains the pair `{AP,GW}`, so more than
one transfer block already violates pair uniqueness.  Thus AP/GW labels cannot
be universal unital points repeated in every block.  They must be chart labels,
branch colors, or external metadata.

## Positive Branch-Local Lift

The obstruction only kills a single global atlas.  It does not kill local
branch charts.

The following desired block systems have no repeated C27 pair and embed into
the q=3 unital incidence design:

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

The S138 genuine two-hole frontier rows also lift, but not as new single
blocks.  They lift as two-block splices:

```text
drop(10,12)->add(20,24)
  = petal H10->D7  plus  GW H12->D3

drop(10,12)->add(20,36)
  = petal H10->D7  plus  K33 H12->D9
```

This exactly matches HYP-2940: the two-hole rows are splices of the unit petal
with one of the two `12`-branch packets.

## Tournament Analysis

Vertices:

```text
exact M/Farey branch,
C27 marked transfer,
unital pair-repeat obstruction,
branch-local q3 unital chart,
S138 two-block splice path,
global AP/GW vertex model,
raw unital analogy.
```

Pair observable, ordered lexicographically:

```text
theorem-scale retention,
C27 predicate retention,
lambda=1 incidence,
finite checkability,
anti-scalar guard.
```

The lift-obligation tournament is transitive:

```text
exact M/Farey branch
> C27 marked transfer
> unital pair-repeat obstruction
> branch-local q3 unital chart
> S138 two-block splice path
> global AP/GW vertex model
> raw unital analogy.
```

Fingerprint:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
c3 = 0
Hamiltonian path count = 1
```

## Proof-Route Consequence

The q=3 unital should now be used as a branch-local pair-unique chart or a
two-block splice grammar, not as a universal design atlas for all C27 marked
transfers.

The pair-repeat obstruction is itself useful:

```text
the shared H12 hole pair separates the tight GW branch from the K33 near-miss
branch.
```

If a future proof tries to put both branches into one unital object, it must
either split the H12 pair by an additional branch label or use multiple unital
charts.  Either move changes the preserved unit and must be declared explicitly.

## Honest Status

This is not an LRC14 proof.  It is a structural filter.  It says:

```text
unital blocks can organize one C27 branch at a time;
they cannot globally unify the GW and K33 12-branch transfers without
violating lambda=1 pair uniqueness.
```

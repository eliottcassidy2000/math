---
id: HYP-2944
title: LRC14 affine-depth unital chain
status: PROOF-INTERFACE / calibrated depth-14 packet grammar; not a proof
source: codex-2026-06-24-S142
related:
  - HYP-2943
  - HYP-2942
  - HYP-2941
  - HYP-2940
  - HYP-2939
  - HYP-2938
  - HYP-2937
  - HYP-2935
  - HYP-2932
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_affine_depth_unital_chain_codex_s142.py
  - 05-knowledge/results/lrc14_affine_depth_unital_chain_codex_s142.out
---

# HYP-2944: LRC14 Affine-Depth Unital Chain

HYP-2941 turned the prompt's triangular/perfect-number equality into an
affine-depth guardrail.  HYP-2942 supplied a calibrated q=3 unital path:

```text
AP/GW --H12-- near/K33 --D9-- petal10.
```

HYP-2944 composes the HYP-2942 transfer components as words in

```text
a(x)=x/2,    b(x)=x+1.
```

The script `04-computation/lrc14_affine_depth_unital_chain_codex_s142.py`
stores output at
`05-knowledge/results/lrc14_affine_depth_unital_chain_codex_s142.out`.

## Component Depth Rule

For a marked C27 transfer `H[h]->D[d]`, set

```text
depth = 1 + depth_gcd(h) + depth_gcd(d),
```

where the C27 strata have depths

```text
unit -> 0,   gcd3 -> 1,   gcd9 -> 2.
```

For the calibrated HYP-2942 components:

```text
GW 12->24          depth 3, block {AP,D3,GW,H12}
near-miss 12->36   depth 4, block {D10,D9,H12,H9}
petal 10->20       depth 1, block {D7,D9,H10,H13}
petal 13->26       depth 1, block {D1,D11,D12,H13}
```

The rule is not asserted as canonical.  It is a compact packet grammar that
retains the unit/nonunit C27 stratum and the calibrated unital block path.

## Two-Block Splice Strips

Writing a component of depth `d` as `b a^d`, the S138 two-block rows become:

```text
P10 then GW:
  component depths [1,3]
  suffix depths [4,3]
  depth sum 7
  beta = 3/16
  block union size 8, disjoint

P10 then near/K33:
  component depths [1,4]
  suffix depths [5,4]
  depth sum 9
  beta = 3/32
  block union size 7, shared D9
```

Reversing either order changes the suffix-depth sum.  Thus the affine packet is
noncommutative, matching the branch-local q=3 unital guardrail.

## Linked Depth-14 Chain

The calibrated HYP-2942 linked path uses the component order

```text
GW -> near/K33 -> petal10.
```

The component depths are

```text
[3,4,1].
```

Composed as affine words, this gives

```text
word = b a^3 b a^4 b a
suffix depths = [8,5,1]
depth sum = 14
beta = 137/256
```

The same multiset of component depths has six possible orders.  Their suffix
depth sums are:

```text
(1,3,4) -> 19
(1,4,3) -> 18
(3,1,4) -> 17
(3,4,1) -> 14
(4,1,3) -> 15
(4,3,1) -> 13
```

So the HYP-2942 linked order is the unique order producing `14`.  This is the
new LRC14 signal: `14` appears as an order-sensitive affine-depth invariant on
the calibrated unital chain, not as a raw scalar equality.

## Tournament Analysis

Tournament vertices are proof carriers, not runners.  The pairwise observable is

```text
(M/Farey retention, unital incidence retention, affine order retention,
 finite checkability, state-lift fit, anti-scalar guard).
```

The conservative carrier tournament is transitive:

```text
exact M/Farey branch
> C27 unital calibrated chain
> affine depth-14 signature
> two-block splice strips
> Kpq/K33 owner packet
> triangular/perfect product lane
> raw scalar equality.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
c3=0
hp=1
```

## Proof Target

The refined local target is:

```text
In the second-gap frontier, any residual that reaches the calibrated
AP/GW--H12--near/K33--D9--petal10 path should carry the depth-14 affine
signature.  Unit-only splices stay in lower triangular strips and should
discharge by C27 petal/two-swap rigidity; depth-14 linked packets should feed
the HYP-2908 / THM-572 forbidden-H state lift.
```

This is not a proof of LRC14.  It is a concrete packet grammar that turns the
triangular/perfect-number prompt into a branch-local, order-sensitive invariant
on known low-frontier transfer packets.

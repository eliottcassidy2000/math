---
id: HYP-2903
status: SHARPENED PROOF OBLIGATION / guardrail
source: codex-2026-06-22
tags: [lrc14, bonferroni, venn, far-packets, activation-depth, open-q-108]
depends_on:
  - HYP-2901
  - HYP-2902
  - HYP-2900
  - THM-548
  - THM-557
  - THM-563
  - THM-564
related:
  - OPEN-Q-108
  - HYP-2788
results:
  - 04-computation/lrc_depth_corrected_bonferroni_codex.py
  - 05-knowledge/results/lrc_depth_corrected_bonferroni_codex.out
---

# HYP-2903: the wide Bonferroni target must be corrected by activation depth

The corrected Venn/Legendre picture suggests a third-order far-packet target

```text
p0(B union F) <= p0(B) + T_1 + T_2 + T_3,
```

where `T_r` is the total `r`-far Newton packet over the far runner set `F`.
This is the right target only after an activation-depth split, in rows whose
first live packet occurs by depth two.

It is false as a bare statement.  Exact rational integration gives the
counterexample

```text
B={0,1,2,3}
F={16,19,22,25,28}

T_1 = 0
T_2 = 0
T_3 = 638291/6144600
T_4 = 17921/6144600
T_5 = -1/931
T_{>=4} = 11321/6144600 > 0.
```

So third order is a lower truncation in this row:

```text
p0(B union F) - (p0(B)+T_1+T_2+T_3) = 11321/6144600.
```

The reason is structural, not numerical.  This base is not edge-active:
the first nonzero packet is the triple packet.  Once the leading packet shifts
from `T_2` to `T_3`, the first upper truncation shifts from order `3` to order
`4`.

Incoming KPS S31u supplies the complementary high-depth slack failure:

```text
B={0,3,9}
F={27,28,42,43,47}

T_1=T_2=T_3=0
T_4=453955/10695132
T_5=-429/56588
p0=186437/5347566.
```

So Bonferroni-3 also fails when the first live packet is at depth `4`, but
this row is not cap-threatening; `p0` is far below the LRC14 cap.  This matches
the KPS S31u readout: Bonferroni-3 is a binding-leg handle, not a universal
wide closure.

## Corrected sharp obligation

Define the activation depth

```text
r0(B,F) = min { r : T_r(B,F) != 0 }.
```

The wide branch should be split by activation depth:

```text
corner/edge rows      r0 <= 2: prove T_{>=4} <= 0, hence order-3 upper bound.
triple-active rows    r0 = 3: prove T_{>=5} <= 0, hence order-4 upper bound.
higher-depth rows     r0 >= 4: prove a slack bound p0 << cap or route to the
                      single-block/decorrelated / finite residual atlas.
```

This preserves the useful S31t insight while removing the overclaim.  The
Legendre three-corner Venn supplies the address labels, but the applicable
Bonferroni truncation level is determined by how many far runners are needed
before the base can cover all six sectors.

## Proof relevance

For the current LRC14 proof order, HYP-2903 changes the sharp obligation from

```text
prove a universal Bonferroni-3 tail sign
```

to

```text
prove activation-depth Bonferroni signs, then close each activation class
with the already separated doublet/triple/single-block estimates.
```

The practical target remains strongest in the edge-active rows: bound
`T_2+T_3` by the cap and prove the `r>=4` tail is nonpositive.  Corner-active
rows with `T_1 != 0` can still use the same order-3 upper target when their
post-third tail is nonpositive.  But
triple-active rows must not be forced through that same third-order truncation.

Post-rebase S46 sharpens the role of this obligation: if the Node-3
equidistribution/induction skeleton is made effective, then the remaining
irreducible work is precisely this bounded Node-2 cap side.  So HYP-2903 should
be read as a Node-2 scope correction, not as a replacement for the analytic
Node-3 input.

## Assumption challenge

Tournament/tiling vertices cannot be raw far runners alone.  The preserved
predicate is the first live Newton layer of the base/far interaction.  The
quotient destroys runner identity but preserves activation depth, Venn layer,
and the sign of the remaining tail.

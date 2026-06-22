---
id: HYP-2903
status: SHARPENED PROOF OBLIGATION / exact counterexamples
source: codex-2026-06-22; extended by codex-2026-06-22-S115
tags: [lrc14, bonferroni, venn, far-packets, activation-depth, missing-depth, newton-packets, open-q-108]
depends_on:
  - HYP-2901
  - HYP-2902
  - HYP-2900
  - THM-548
  - THM-557
  - THM-563
  - THM-564
related:
  - HYP-2708
  - HYP-2788
  - HYP-2889
  - HYP-2890
  - OPEN-Q-108
results:
  - 04-computation/lrc_depth_corrected_bonferroni_codex.py
  - 05-knowledge/results/lrc_depth_corrected_bonferroni_codex.out
  - 04-computation/lrc_bonferroni_depth_guard_codex_s115.py
  - 05-knowledge/results/lrc_bonferroni_depth_guard_codex_s115.out
---

# HYP-2903: the wide Bonferroni target must be corrected by missing-depth parity

The corrected Venn/Legendre picture suggests a third-order far-packet target

```text
p0(B union F) <= p0(B) + T_1 + T_2 + T_3,
```

where `T_r` is the total `r`-far Newton packet over the far runner set `F`.
This was the natural first target in the binding edge/corner regime, but it is
not a universal wide theorem.

Exact rational integration first gave the activation-depth counterexample

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

The reason is structural.  This base is not edge-active: the first nonzero
packet is the triple packet.  Once the leading packet shifts from `T_2` to
`T_3`, the first upper truncation shifts from order `3` to order `4`.

Incoming KPS S31u supplies the complementary high-depth slack failure:

```text
B={0,3,9}
F={27,28,42,43,47}

T_1=T_2=T_3=0
T_4=453955/10695132
T_5=-429/56588
p0=186437/5347566.
```

Bonferroni-3 also fails when the first live packet is at depth `4`, but this
row is not cap-threatening; `p0` is far below the LRC14 cap.  This matches the
KPS S31u readout: Bonferroni-3 is a binding-leg handle, not a universal wide
closure.

The S115 computation shows that activation depth is still not the full
invariant.  A second exact counterexample is already edge-active:

```text
B={0,1,7,10,13}
F={15,23,24,31}

T_0 = 0
T_1 = 0
T_2 = 590813/13625430
T_3 = 631/42780
T_4 = 307/598920
T_{>=4} = 307/598920 > 0.
```

Thus `r0=2` does not by itself imply the Bonferroni-3 upper bound.  The
surviving object is the missing-depth ledger.

## Exact pointwise object

Fix a slow time `x`.  Let `M(x)` be the set of inner sectors not already hit by
the base `B`, and put `d(x)=|M(x)|`.  If far runner `i` has sector color
`c_i(x)`, then the pointwise full-cover predicate contributed by the far
runners is

```text
prod_{a in M(x)} OR_{i : c_i(x)=a} z_i.
```

Therefore the integrated Newton packet is

```text
T_r = sum_d (-1)^(r+d) A_{d,r},
```

where `A_{d,r}` is the measure-weighted count of `r`-subsets of far runners
whose colors all lie in the `d` missing sectors and cover those `d` sectors.
Equivalently, at each `x`,

```text
T_r(x) = (-1)^(r+d(x)) C_{d(x),r}(x).
```

The high tail after order three is not controlled by Venn containment alone:

```text
Tail_3 = sum_{r>=4} T_r
       = sum_{d,r>=4} (-1)^(r+d) A_{d,r}.
```

The KPS binding-style row

```text
B={0,1,2,3,4}
F={20,21,22,23,24}
```

has negative tail `-13141/212520`, because its depth-3 negative contribution
dominates the positive depth-2 and depth-4 contributions.  The k=8
counterexample

```text
B={0,1,2,3}
F={16,28,29,32}
```

has positive tail `19/68208`, because the positive depth-4 high packet exceeds
the negative depth-3 high packet:

```text
d=3 tail = -163/45472
d=4 tail =  527/136416.
```

This is the concrete proof correction: the tail sign is a depth-parity
inequality, not a Venn-containment theorem and not an activation-depth theorem.

## Corrected sharp obligation

Define the activation depth

```text
r0(B,F) = min { r : T_r(B,F) != 0 }.
```

Use it as a routing variable, not as the final sign invariant:

```text
corner/edge rows      r0 <= 2: inspect the depth-parity tail before using
                      order 3 as an upper bound.
triple-active rows    r0 = 3: shift the candidate upper target to order 4
                      and still inspect the depth-parity tail.
higher-depth rows     r0 >= 4: prove a slack bound p0 << cap or route to the
                      single-block/decorrelated / finite residual atlas.
```

Inside each activation class, prove the actual depth-parity guard

```text
sum_{r>=4, r+d even} A_{d,r}
  <= sum_{r>=4, r+d odd} A_{d,r} + cap-slack(B,F).
```

Rows with too much even-depth high mass must be discharged by the existing
single-block/decorrelated extremality, by KPS S31u's spread-far slack branch,
or by a finite residual atlas.  This preserves the useful S31t insight while
removing the overclaim.  The Legendre three-corner Venn supplies the address
labels, but the applicable Bonferroni truncation and remaining sign are
determined by the depth profile of the base/far interaction.

## Proof relevance

For the current LRC14 proof order, HYP-2903 changes the sharp obligation from

```text
prove a universal Bonferroni-3 tail sign from Venn containment
```

to

```text
prove the depth-parity high-tail inequality in the binding leg, then close
positive-tail activation/depth classes with cap slack, doublet/triple
estimates, spread-far slack, or decorrelation.
```

Post-rebase S46 sharpens the role of this obligation: if the Node-3
equidistribution/induction skeleton is made effective, then the remaining
irreducible work is precisely the bounded/binding Node-2 cap side.  So
HYP-2903 should be read as a Node-2 scope correction, not as a replacement for
the analytic Node-3 input.

## Assumption challenge

Tournament/tiling vertices cannot be raw far runners alone.  The preserved
predicate is the depth-labelled Newton packet of the base/far interaction.
The quotient destroys runner identity but preserves missing-sector depth,
activation depth, Venn layer, cap slack, and the sign of the remaining tail.

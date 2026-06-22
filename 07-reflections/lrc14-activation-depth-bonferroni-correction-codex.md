---
source: codex-2026-06-22
status: reflection / proof-obligation sharpening
tags: [lrc14, bonferroni, venn, activation-depth, open-q-108]
---

# The Bonferroni truncation level is an activation-depth address

The S31t/HYP-2901 Venn readout was directionally right: the far-runner
expansion is a Newton/Mobius expansion, and the doublet/triple layers are the
right named objects.  The overclaim was treating order three as the same
truncation level for every wide row.

Exact rational integration gives the guardrail:

```text
B={0,1,2,3}
F={16,19,22,25,28}

T1=T2=0
T3=638291/6144600
T4=17921/6144600
T5=-1/931
T_{>=4}=11321/6144600 > 0.
```

So order three underestimates `p0`.  This is not a numerical quirk; the first
live packet is the triple.  If the first live packet is depth three, the first
upper truncation moves to order four.

The sharper obligation is therefore:

```text
r0=min{r:T_r != 0}

r0<=2: prove T_{>=4}<=0, close with order 3.
r0=3:  prove T_{>=5}<=0, close with order 4.
r0>=4: prove p0<<cap or route to single-block/decorrelated extremality.
```

This adds one missing address to the Legendre/Venn sheaf: activation depth.
It preserves the LRC predicate better than raw runner vertices or a scalar
Bonferroni order.  The practical proof target is still the edge-active branch,
where `T2+T3` should be compared to the cap and the `r>=4` tail should help.
But triple-active rows must not be pushed through a false third-order upper
bound.

For the next proof pass, the useful split is:

1. Prove an activation-depth classifier for the wide families currently routed
   to S31t.
2. In the `r0<=2` class, prove signed/decreasing `T_{>=4}<=0`.
3. In the `r0=3` class, prove `T_{>=5}<=0` and compare the fourth-order
   truncation to the cap or route it to THM-557-style single-block estimates.
4. In the `r0>=4` class, use KPS S31u's scope correction: these are spread-far
   slack rows unless a separate finite atlas says otherwise.

Incoming mac-mini S46 reinforces the placement of this correction.  If its
Node-3 equidistribution/induction skeleton is made effective, the LRC14 burden
falls back onto bounded Node-2 cap statements.  Activation depth is therefore
the address needed inside Node 2; it is not a substitute for the Node-3 Weyl
input.

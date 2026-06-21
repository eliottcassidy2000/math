---
id: HYP-2790
title: THM-563 bounded-base period maxima should be explainable by Boolean/type cut slack
status: OPEN; claimed synthesis target, computation underway
source: codex-2026-06-21
depends_on:
  - THM-563
  - HYP-2788
  - HYP-2751
  - HYP-2752
  - HYP-2789
related:
  - HYP-2780
  - HYP-2781
  - HYP-2770
  - HYP-2744
  - THM-534
  - OPEN-Q-108
---

# HYP-2790: Period-Max / Boolean-Type Bridge

## Claim Under Test

`THM-563` turns the single-far signed deviation into a finite endpoint-period
maximum.  `HYP-2788` says the wide near-cap branch routes to
single-perturbation bounded scaffolds.  The remaining finite object is therefore
a list of bounded bases `B` and the inequality

```text
periodmax(B) < 15 * (cap_k - Plat(B)).
```

This hypothesis tests whether the same Boolean/type slack discovered in
`HYP-2751` and the containment-cut refinements of `HYP-2752` explain why the
dangerous bounded bases satisfy the period-max inequality.

## Planned Finite Ledger

For each bounded base `B` with `0 in B`, `max(B)<=14`, and final row size
`k=|B|+1`, compute:

```text
Plat(B) = p0(B) + p1(B)/7,
margin(B) = cap_k - Plat(B),
periodmax(B) = max_w w*Delta_w(B),
ratio(B) = periodmax(B) / margin(B),
```

then compare `ratio(B)` with Boolean/type features of the missed-sector law of
`B`, especially:

```text
T1=(1,(1)), T2sep=(2,(1,1)), T2adj=(2,(2)),
Phi_low = 21*T1 + 57*T2sep + 2*T2adj.
```

The first goal is not to prove the inequality.  It is to determine whether the
top period-max pressure rows form a small family already visible to the
Boolean/type quotient, sorted-cell leak quotient, or affine gap-word quotient.

## Why This Is Not A Duplicate

`HYP-2751` is a full k=8 bounded-row atom-run-type certificate.  `HYP-2752`
shows fixed-k low-level Mobius/containment cuts exist but are non-uniform and
not structurally free.  `HYP-2790` asks a different question: after the wide
branch has been reduced to single-far bounded scaffolds, can the finite
period-max pressure itself be ranked or certified by the Boolean/type slack?

If yes, this gives a route from:

```text
wide near-cap -> single-perturbation bounded scaffold
              -> THM-563 endpoint-period max
              -> Boolean/type finite ledger
```

instead of a separate ad hoc period scan for every bounded base.

## Assumption Challenge

The vertices for this test are bounded bases/scaffolds, not final rows or
runner arcs.  The quotient preserves missed-sector cyclic shape and endpoint
period pressure; it destroys speed ownership inside wall atoms and the full
Fourier phase distribution.  If the test fails, the likely missing coordinate
is the endpoint-period numerator itself, not another scalar energy statistic.

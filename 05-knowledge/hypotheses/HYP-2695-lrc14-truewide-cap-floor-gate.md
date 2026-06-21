---
id: HYP-2695
title: LRC14 true-wide Bonferroni4 cap-floor gate
status: OPEN
source: codex-2026-06-20-S60
tangent: T930
depends_on:
  - THM-535
  - THM-556
  - HYP-2693
  - HYP-2675
related:
  - HYP-2701
  - HYP-2694
  - HYP-2691
  - HYP-2692
  - HYP-2683
  - HYP-2684
  - HYP-2676
  - HYP-2677
  - THM-547
  - THM-548
  - OPEN-Q-108
---

# HYP-2695 - LRC14 True-Wide Bonferroni4 Cap-Floor Gate

## Claim

HYP-2693 should be sharpened by splitting the exact cap into a universal floor
and a small dividend:

```text
floor_k    = (k - 6) / 7,
dividend_k = cap_k - floor_k.
```

THM-535 proves `cap_k >= floor_k`.  Therefore, on any row where the
Bonferroni4 expression satisfies `U4(E)<=floor_k`, the HYP-2693 cap gate
follows without using exact cap minimizers.

The proposed sharp target is:

```text
second_largest(E)>14, |E|=k>=9
    => U4(E) = p0(E)+p5(E)+5p6(E) <= (k-6)/7.
```

The row `k=8` is genuinely exceptional.  True-wide rows can exceed the
subadditive floor `2/7`, but the excess appears finite-template and is covered
by the exact `k=8` cap dividend `cap_8-floor_8=563/5880`.

Thus the proof split becomes:

```text
k>=9 true-wide  -> prove Bonferroni4 floor gate
k=8 true-wide   -> finite cap-dividend templates
boundary/AP     -> HYP-2691 transfer-DP / low-state template branch
```

This is complementary to incoming HYP-2694's single-block wide-cover
extremizer claim.  HYP-2694 attacks the decorrelated wide-cover partition
function and identifies the coherent block as the comfortable-margin extremal;
HYP-2695 is the final-row Bonferroni currency split after the true-wide row has
been routed to a high-tail or finite-template branch.

## Evidence

The exact scout `04-computation/lrc14_truewide_cap_floor_gate_codex_s60.py`
stores its run in
`05-knowledge/results/lrc14_truewide_cap_floor_gate_codex_s60.out`.
It asserts THM-556's identity `U4=p0+p5+5p6` on every audited row and compares
`U4` with both `floor_k` and `cap_k`.

Cap-currency table:

```text
k=8:  floor=2/7, cap=2243/5880, dividend=563/5880
k=9:  floor=3/7, cap=1979/4004, dividend=263/4004
k=10: floor=4/7, cap=55/91, dividend=3/91
k=11: floor=5/7, cap=66/91, dividend=1/91
k=12: floor=6/7, cap=6/7, dividend=0
```

Exact staged boxes:

```text
k=8,  B18: truewide=16359, floor violations=3, cap violations=0
k=9,  B18: truewide=27020, floor violations=0, cap violations=0
k=10, B16: truewide=3432,  floor violations=0, cap violations=0
k=11, B16: truewide=3003,  floor violations=0, cap violations=0
k=12, B16: truewide=2002,  floor violations=0, cap violations=0
```

The `k>=9` tight rows in the audit are:

```text
k=9, E=(0,4,6,8,10,12,14,15,16):
  U4=391/980, floor-U4=29/980, cap-U4=3338/35035.

k=10, E=(0,2,4,6,8,10,12,14,15,16):
  U4=307/588, floor-U4=29/588, cap-U4=629/7644.

k=11, E=(0,2,4,6,8,9,10,12,14,15,16):
  U4=73/126, floor-U4=17/126.

k=12, E=(0,4,6,8,9,10,11,12,13,14,15,16):
  U4=274891/420420, floor-U4=85469/420420.
```

The `k=8` exception is exact and small:

```text
E=(0,3,6,9,12,14,15,18):
  p0=391/1764, tail45=8/105,
  U4=2627/8820,
  floor-U4=-107/8820,
  cap-U4=295/3528.
```

Only three true-wide `k=8` rows in the `B18` box exceed the floor; none exceed
cap.  The worst floor debt `107/8820` is much smaller than the cap dividend
`563/5880`.

Boundary/AP rows confirm why the floor-gate statement must keep the true-wide
stratum.  AP9 and the doubled AP boundary row have
`U4=1621/2940`, so they fail both floor and cap for `k=9`; the boundary leader
`(0,2,4,6,8,10,12,14,15)` fails the floor by `39/1960` but still passes cap.

## S64 Survival-Currency Refinement

HYP-2701 changes the coordinate for this gate.  Since `sum_t p_t=1`,

```text
floor_k - U4(E)
  = p1+p2+p3+p4 - 4p6 - (13-k)/7.
```

Thus the `U4` floor gate is exactly a lower bound on survival middle mass after
charging the fully-missed tail.  The exact cap gate uses the same left side with
right side `1-cap_k`.  The S64 scan extends the true-wide boxes to k=10,11,12
at B17 and finds the same qualitative pattern: no cap failures, no `k>=9`
floor failures, and the only audited floor failures are three cap-safe k=8
rows.  Its far-count ledger sharpens this hypothesis' proof route: prove a
two-far survival-currency lemma first, then prove a separate `far_count>=3`
margin lemma, while keeping k=8 on the finite dividend branch.

The HYP-2701 two-far addendum further rewrites the first branch as a THM-548
boundary-value problem: the decorrelated two-far death-chain currency is already
floor-safe for every bounded core scanned, so the open inequality is exactly
positive boundary margin versus negative signed two-far deviation.

## Proof Route

1. **Floor gate for `k>=9` true-wide rows.**  Use THM-556 to work only with
   `p0+p5+5p6`.  Prove that two-or-more-far true-wide structure suppresses the
   high missed-sector tail enough to fit below `(k-6)/7`.
2. **State-mass/decorrelation split.**  HYP-2683 suggests that high-risk
   true-wide rows have concentrated missed-state support and lower entropy,
   while HYP-2684 supplies the BV/Fourier nonresonant error form.  The floor
   gate should be proved by combining these: finite low-state packets plus
   decorrelated high-state rows.
3. **Finite resonant routers.**  HYP-2676/HYP-2677 show that same-sign packet
   models and tournament sign quotients locate finite Ruzsa/Freiman pockets,
   but raw additive energy is too coarse.  The resonant high-tail exceptions
   should retain packet magnitudes, sector-state profiles, divisibility, and
   low-height relation type before scalarizing.
4. **Cap-dividend row.**  Do not force `k=8` through the floor proof.  Its
   anomalous cap minimizer in THM-535 is real proof currency; the exact
   dividend `563/5880` should close a finite template ledger.

## Assumption Challenge

This hypothesis uses proof obligations as vertices: floor gate, cap dividend,
finite boundary template, state-mass deficit, decorrelated plateau, and
packet/Ruzsa resonant router.  It does not treat runners or arcs as the default
vertex set.  The quotient preserves the LRC implication because
`p0<=U4<=floor_k<=cap_k` closes the sector cap predicate for `k>=9`; it
destroys sector labels and phase, so suspected near-equality rows must be
lifted back to addressed state packets before any theorem is claimed.

## Status

No LRC14 proof is claimed.  The new sharp proof obligation is the true-wide
floor gate for `k>=9`, plus a finite cap-dividend treatment for true-wide
`k=8`.  The exact audit makes the old HYP-2693 target more structured:
most true-wide rows should not spend any exact cap dividend at all.

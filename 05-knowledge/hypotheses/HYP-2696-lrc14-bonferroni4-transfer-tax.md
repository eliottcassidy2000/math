---
id: HYP-2696
title: LRC14 Bonferroni4 transfer-tax router
status: OPEN
source: codex-2026-06-20-S62
tangent: T932
depends_on:
  - THM-555
  - THM-556
  - THM-558
  - HYP-2691
  - HYP-2693
related:
  - HYP-2675
  - HYP-2676
  - HYP-2677
  - HYP-2680
  - HYP-2681
  - HYP-2682
  - HYP-2683
  - HYP-2684
  - HYP-2694
  - HYP-2695
  - HYP-2692
  - THM-557
  - OPEN-Q-108
---

# HYP-2696 - LRC14 Bonferroni4 Transfer-Tax Router

## Claim

THM-558 should be used as the local accounting law between the sector-state
transfer DP of HYP-2691 and the true-wide Bonferroni4 gate of HYP-2693.

For a prefix insertion `P -> P union {e}`,

```text
Delta U4 = mass(1 -> 0) - mass(5 -> 4) - 4 mass(6 -> 5),
U4 = p0 + p5 + 5p6.
```

Therefore a proof of the HYP-2693 gate should not try to control all missed
states uniformly.  It should classify or bound the possible unpaid
`1 -> 0` closures after subtracting the high-tail tax from `5 -> 4` and
`6 -> 5` transitions.

The proposed router is:

```text
unpaid closure packets
  -> finite AP/dyadic/cube-root/Ruzsa templates when structured,
  -> Weyl/BV/decorrelation bounds when true-wide and high-state.
```

The incoming THM-557/HYP-2694 single-block theorem supplies the complementary
wide-cover compression target: coherent far blocks have a proved decorrelated
extremizer and explicit diagonal-freeze error.  HYP-2696 is the local
final-row ledger to apply after that compression, or to diagnose the residual
unpaid closure packets that compression does not yet cover.

## Evidence

The script `04-computation/lrc14_bonferroni4_transfer_tax_codex_s62.py`
stores its run in
`05-knowledge/results/lrc14_bonferroni4_transfer_tax_codex_s62.out`.
It imports the HYP-2691 row bank and verifies THM-558 by assertion for every
audited step and insertion schedule.

The global best-schedule summary shows the split already visible in named
rows:

```text
AP9:
  U4=1621/2940, cap9-U4=-6002/105105, tail45=53/392.

one-gap-top9:
  U4=1531/2940, cap9-U4=-5569/210210, tail45=5/42.

direct-risk true-wide leader (0,4,6,8,10,12,14,15,16):
  U4=391/980, cap9-U4=3338/35035, tail45=1/14.

boundary leader (0,2,4,6,8,10,12,14,15):
  U4=879/1960, cap9-U4=12833/280280, tail45=113/1470.

dyadic block:
  U4=103/336, cap9-U4=9019/48048, tail45=1/21.

three-cluster true-wide row:
  U4=849937/3176880,
  cap9-U4=102996349/454293840,
  tail45=23407/529480.
```

The AP and one-gap rows fail the final `U4<=cap` gate because a finite
template high tail survives.  The direct true-wide, dyadic, doubled-odd,
three-cluster, and AP-triple rows have enough transfer tax that the final
`U4` remains comfortably below cap.

The exact run also confirms the order-dependent proof burden.  For example,
the dyadic block is best routed by the `dyadic-tower` schedule, while AP-like
and AP-triple rows prefer increasing order in the transfer-tax key.  The
Tournament Analysis vertices are insertion schedules/proof states, not
runners; the observed schedule tournaments are transitive in the audited bank.

## Proof Route

The useful inequality is not just `p0<=U4`.  The transfer identity gives a
ledger:

```text
final U4 - 5 = total_close1 - total_tax,
total_tax = total_hit5 + 4 total_hit6,
final_tail45 = 5 - total_tax.
```

Thus a row fails the Bonferroni4 gate only if its accumulated one-missed
closures are not paid down by enough high-tail tax.  This suggests three
sublemmas.

1. **Finite unpaid templates.** If many `1 -> 0` closures happen with small
   transfer tax, then the prefix state is low-growth/low-support and should
   land in the existing AP-prefix, dyadic, cube-root/AP-triple, or Ruzsa
   packet atlases.
2. **True-wide high-tail tax.** In the true-wide branch, prove that closures
   either occur after enough `5 -> 4` and `6 -> 5` tax has been paid, or are
   small by HYP-2684 Weyl/BV decorrelation.
3. **Phase-lift exceptions.** If a true-wide row has unexpectedly large
   unpaid closure mass, lift it from missed-counts back to the labelled
   sector-state word and phase atlases HYP-2681/HYP-2682 before scalarizing.

This is compatible with HYP-2692: the apex-divisor arithmetic organizes the
resonance classes, while the effective archimedean quantity is the signed
leading residual.  HYP-2696 gives the same idea as a finite transfer ledger
for the Bonferroni4 certificate.

It also dovetails with HYP-2694/THM-557: the coherent-block theorem pushes
the true-wide wide-cover branch toward single-block extremizers, while
THM-558 measures the exact signed Bonferroni pressure when those compressed
rows are inserted one runner at a time.

## Assumption Challenge

The vertices in this hypothesis are not runners, arcs, or raw sector labels.
They are transition types, missed-state counts, insertion schedules, and proof
obligations.  This quotient preserves the LRC-relevant upper bound
`p0(E)<=U4(E)` and the exact one-step changes in `U4`.  It destroys the
cyclic sector labels and phase information, so any near-equality or exception
must be lifted back to the labelled sector-state DP.

## Status

THM-558 is proved.  HYP-2696 is the open proof program: use the transfer-tax
ledger to turn HYP-2693's final true-wide gate into a classification of
unpaid one-missed closure packets.

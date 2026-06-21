---
id: HYP-2719
title: LRC14 Weyl error is odd-support-envelope dominated, not a signed odd cone
status: OPEN; exact scout evidence
source: codex-2026-06-21-S69
depends_on:
  - HYP-2718
  - THM-561
  - HYP-2717
  - HYP-2716
  - HYP-2715
related:
  - HYP-2698
  - THM-558
  - HYP-2676
  - OPEN-Q-108
---

# HYP-2719 - Odd-Support Weyl Error Envelope

## Claim

The HYP-2718 origin-atom/Weyl error should be controlled by an odd-support
envelope in the factorial-moment basis, but not by a naive signed odd cone.

For the actual line law minus the shared-slow-`x` independent carrier product
law, write

```text
W_j = Delta E[binom(T,j)]
    = sum_{|R|=j} (z_product(R)-z_actual(R)),
```

where `T` is the number of missed inner sectors.  THM-561 gives

```text
Q_0 = ProductCover - ActualCover
    = sum_j (-1)^j W_j
    = (sum_{j even} W_j) - (sum_{j odd} W_j).
```

Thus "odd support dominated" must be made precise.  The exact scout supports
this version:

```text
OddAbs = sum_{j odd} |W_j|
EvenAbs = sum_{j even} |W_j|
```

often has `OddAbs >= EvenAbs`, and aggregate odd support has the larger
unsigned envelope.  But the signed origin atom is usually produced by
cancellation between two negative signed aggregates, so signed odd dominance is
false in most tested rows.  The useful theorem should therefore be:

```text
high-height/nonresonant tail: odd-support L1 envelope controls available mass;
low-height signed-even-led packets: finite HYP-2717/HYP-2714 ledger;
origin atom Q_0: bounded only after this split.
```

## Exact Scout

Script:

```text
04-computation/lrc14_odd_support_weyl_error_codex_s69.py
```

Stored output:

```text
05-knowledge/results/lrc14_odd_support_weyl_error_codex_s69.out
```

The scout evaluates two-carrier rows, the existing HYP-2715 split bank, and
the S68 relation-height examples.  It prints exact `W_j`, atom deltas `Q_t`,
odd/even signed support, odd/even `L1` support, and a tournament over
odd-support pressure.

Aggregate signal over the ten-row bank:

```text
aggregate even_L1 = 1.041714308...
aggregate odd_L1  = 1.169709905...
odd_L1_share      = 0.528939630...
flag_counts       = {(signed_odd_dominates=False, sign_from_odd=False): 8,
                     (signed_odd_dominates=True,  sign_from_odd=True): 2}
```

So odd support dominates the unsigned envelope in aggregate, while signed odd
dominance occurs only on two positive-`Q_0` arithmetic phases.  The largest
odd-share row in the scout is the five-2-block row:

```text
odd_L1_share = 0.632953886...
Q_0          = -2447628624709/93106921650624.
```

The signed aggregates there are both negative:

```text
even_support = -0.128046613...
odd_support  = -0.101758250...
Q_0          = even_support - odd_support < 0.
```

This is the main correction: odd support can dominate the available fluctuation
without dominating the final signed origin atom.

## Proof Route

1. Use THM-561 to keep the factorial support layers `W_j` until the last step.
2. Prove an odd-support `L1` envelope for the high-height/high-denominator
   HYP-2717 carrier tail.
3. Identify signed-even-led low-height packets by relation vector, phase, and
   generated miss-zeta word; route them to the finite HYP-2714 ledger.
4. Only then evaluate the origin atom `Q_0=Even-Odd` and compare it with
   `cap_k-ProductCover`.

This keeps the user's odd-support intuition, but places it in the right
currency: an envelope before scalarization, not a pointwise signed dominance
claim.

## Tournament Analysis

Vertices are tested split rows.  Pairwise observable:

```text
larger OddAbs/(OddAbs+EvenAbs),
then signed odd dominance,
then larger cap-risk ratio |Q_0|/(cap-product).
```

Switch/gauge: pass from residual masks to factorial-support parity before
origin-atom scalarization.  The scout tournament is transitive:

```text
five 2-blocks
> two 4-blocks wider gap
> two 4-blocks high relation phase
> two 4-blocks moderate gap
> two positive-Q0 4-block phases
> two 4-blocks ratio 2:1
> 3+3+2 split
> seven singleton carriers
> 5+3 split.
```

Fingerprint: score histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}`,
zero directed 3-cycles.

## Status

This is not an LRC(14) proof.  It corrects the phrase "odd support dominated"
into a proof-shaped split: odd-support envelope for the tail, finite ledger for
signed-even-led packets, then the HYP-2718 origin-atom bound.


---
id: HYP-2693
title: LRC14 true-wide Bonferroni4 high-tail gate
status: OPEN
source: codex-2026-06-20-S59
tangent: T929
depends_on:
  - THM-556
  - THM-555
  - HYP-2675
  - HYP-2683
  - HYP-2684
  - HYP-2691
related:
  - HYP-2692
  - THM-546
  - THM-548
  - HYP-2676
  - HYP-2677
  - HYP-2679
  - HYP-2680
  - HYP-2681
  - HYP-2682
  - OPEN-Q-108
---

# HYP-2693 - LRC14 True-Wide Bonferroni4 Gate

## Claim

The HYP-2675 LRC14 wide proof obligation should split at the row's second
largest speed.

For a row `E` in the seven-sector model, let `p_t(E)` be the exact measure of
times at which the row misses exactly `t` of the six inner sectors.  Let
`cap_k` denote the HYP-2675 sector cap for `k=|E|`.

The proposed true-wide gate is

```text
second_largest(E) > 14 and span(E)>14
    => p_0(E) + p_5(E) + 5 p_6(E) <= cap_k.
```

By THM-556, the left side is exactly the level-4 Bonferroni upper expression

```text
U4(E) = 1 - S_1(E) + S_2(E) - S_3(E) + S_4(E).
```

Thus the true-wide branch can be discharged by proving a high-tail bound:
even if `p0` is not tracked directly, the extra mass at five or six missed
inner sectors is too small to consume the cap margin.

Rows with `second_largest(E)<=14` but `span(E)>14` are boundary/AP-template
rows.  They are not expected to satisfy the level-4 gate uniformly; they route
to HYP-2691's finite low-state templates, including AP-prefix append, dyadic,
cube-root/AP-triple, and Ruzsa/Freiman packet ledgers.

## Evidence

The script `04-computation/lrc14_truewide_bonferroni4_gate_codex_s59.py`
stores its run in
`05-knowledge/results/lrc14_truewide_bonferroni4_gate_codex_s59.out`.  It
asserts THM-556 on every audited row and scans exact primitive boxes:

```text
k=8,  bound=18: truewide=16359, violations=0
k=9,  bound=18: truewide=27020, violations=0
k=10, bound=16: truewide=3432,  violations=0
k=11, bound=15: truewide=0,     violations=0
k=12, bound=15: truewide=0,     violations=0
```

The named AP and boundary rows show why level 4 is not a universal proof:

```text
AP9 and doubled AP boundary:
  p0=2447/5880, tail45=p5+5p6=53/392,
  U4=1621/2940, cap9-U4=-6002/105105.

boundary leader (0,2,4,6,8,10,12,14,15):
  p0=437/1176, tail45=113/1470,
  U4=879/1960, cap9-U4=12833/280280.
```

The true-wide leaders in the scanned boxes remain below cap after adding the
exact high tail:

```text
k=8 tightest true-wide row (0,3,6,9,12,14,15,18):
  p0=391/1764, tail45=8/105,
  U4=2627/8820, cap8-U4=295/3528.

k=9 direct true-wide leader (0,4,6,8,10,12,14,15,16):
  p0=321/980, tail45=1/14,
  U4=391/980, cap9-U4=3338/35035.

k=10 tightest true-wide row (0,2,4,6,8,10,12,14,15,16):
  p0=265/588, tail45=1/14,
  U4=307/588, cap10-U4=629/7644.
```

The exact equality of AP9 and the doubled AP boundary row confirms that
integer dilation can preserve the missed-sector profile while crossing the
span threshold.  This is the same finite-template phenomenon seen in
HYP-2691's AP-prefix transfer audit.

Incoming HYP-2692 gives the adjacent apex-divisor / leading-order-residual
view: the arithmetic `Q(sqrt(-7))` structure indexes resonances, but the
effective bound should be on the leading summed residual `R_{s0}`.  HYP-2693 is
compatible with that redirect.  It supplies a final-row upper target after the
leading residual and finite-template branches have been separated: the
true-wide row must still pass the exact high-tail inequality
`p0+p5+5p6<=cap`.

## Proof Route

THM-556 turns the level-4 sieve into a two-part target:

```text
p0(E) + high_tail(E) <= cap_k,
where high_tail(E)=p5(E)+5p6(E).
```

This suggests the following branch proof for HYP-2675.

1. **True-wide high-state branch.** Prove the gate above using sector moment
   constraints, missed-state entropy/support constraints from HYP-2683, and
   Weyl/BV decorrelation from HYP-2684.  The target is not to bound every
   inclusion-exclusion term separately; it is to suppress the five-six missed
   tail after the fourth Bonferroni cancellation has killed `p1..p4`.
2. **Boundary/AP low-state branch.** When the second-largest speed is at most
   `14`, route the row to HYP-2691's finite address templates.  The failing
   AP and doubled AP rows are expected here, and their tail is structural
   rather than a high-state decorrelation error.
3. **Resonant exceptions.** If a true-wide row has unusually large
   `high_tail`, classify its low-height relation packet through the existing
   AP-triple/cube-root and Ruzsa/Freiman atlases before applying the high-state
   bound to the residual.

## Assumption Challenge

This hypothesis uses vertices given by missed-sector cardinalities and row
strata, not by runners, arcs, or individual sectors.  The quotient preserves a
rigorous upper bound for the LRC predicate `p0(E)<=cap_k`, because
`p0<=p0+p5+5p6=U4`.  It destroys sector labels and cyclic packet phase, so any
near-equality or failure case must be lifted back to the sector-state DP and
phase atlases before being treated as a proof.

## Status

No LRC14 proof is claimed.  The new sharp target is the exact true-wide gate

```text
second_largest(E)>14 => p0(E)+p5(E)+5p6(E)<=cap_|E|.
```

The complementary boundary/AP branch remains the HYP-2691 finite-template
program.

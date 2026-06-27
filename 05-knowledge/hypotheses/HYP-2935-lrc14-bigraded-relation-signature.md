---
id: HYP-2935
title: LRC14 bigraded summand/multiplicand relation signature
status: PROOF-INTERFACE / typed relation-channel refinement; not a proof of LRC14
source: codex-2026-06-23-S134
related:
  - HYP-2934
  - HYP-2933
  - HYP-2932
  - HYP-2931
  - HYP-2930
  - HYP-2639
  - HYP-2640
  - HYP-2637
  - HYP-2083
  - HYP-2161
  - HYP-2908
  - THM-572
---

# HYP-2935: Bigraded summand/multiplicand relation signature

S134 reconnects the last Farey-product thread with the older summand graph and
multiplicand graph work.  The key rule is:

```text
addition creates relation shells; multiplication tests their visibility.
```

This does not replace the S130-S133 Farey split.  It refines it by recording,
inside each branch, whether the additive relation is observer-visible, hidden,
balanced, or multiplicand-clearable.

## Computation

The script
`04-computation/lrc14_bigraded_relation_signature_codex_s134.py` stores output
at
`05-knowledge/results/lrc14_bigraded_relation_signature_codex_s134.out`.

For each speed row `S`, it records:

```text
M(S)=p/q and Farey branch
pair-sum support and excess over the AP minimum
visible folds a+b=c in S
balanced pair-sum collisions a+b=c+d split by visible/hidden shell
exact-denominator divisibility blockers
C=27 antipodal shell profile
K_{p,q} rank: star / two-block / K33-wall
```

The selected-row signatures show why scalar additive energy is not enough:

```text
row                 M      branch                 sum+ fold visC hidC  K
AP                  1/14   tight-floor              0   36   55   70   star
GW 12->24           1/14   tight-floor             11   31   40   50   star
petal 10->20        2/27   C27-petal-two-block      7   32   39   51   two-block
petal 13->26        2/27   C27-petal-two-block     10   30   40   55   two-block
near-miss 12->36    3/41   K33-unit-excess         11   30   40   50   K33-wall
shiftAP +10         5/16   nonunit-excess           0    2    0  125   K33-wall
```

The shifted AP calibration is the HYP-2639 warning in miniature.  It has the
same AP-style sumset excess as AP, but its observer-visible folds collapse from
`36` to `2` and its hidden balanced payload jumps to `125`; it is therefore
very safe despite identical raw additive shape.

## Bank-Level Readout

On the S130 `749`-row bank, exact Farey branch remains the best primary split.
The typed relation signature then tells which old graph channel to keep inside
the split:

```text
branch                  rows  avg sum+  avg fold  avg hidden collision  K-ranks
C27-petal-two-block        2      17/2        31                 53     two-block
K33-unit-excess            1        11        30                 50     K33-wall
q-parent-star             54    248/27     275/9            1510/27     star
tight-floor                2      11/2      67/2                 60     star
```

The C27 `p=2` branch carries a two-block multiplicand tag plus a typed
antipodal shell.  The `p>=3` branch carries three-owner incidence.  Raw
sumset, raw energy, or raw product alone does not know this routing.

## Tournament Analysis

The tournament vertices are typed relation channels:

```text
q/Farey branch
C27 typed shell
multiplicand clearance
Kpq incidence
visible folds
hidden balanced shells
raw sumset/energy
raw runner vertices
```

The pair observable is a role score:

```text
theorem scale, branch separation, sign visibility,
multiplicand clearance, old-repo maturity, scalar-decoy resistance.
```

The resulting tournament is transitive:

```text
q/Farey branch
> C27 typed shell
> multiplicand clearance
> Kpq incidence
> visible folds
> hidden balanced shells
> raw sumset/energy
> raw runner vertices
```

An old-repo-maturity-only gauge would flip `5` channel pairs.  That is the
guardrail: the old summand/multiplicand tools are mature and valuable, but they
must be inserted after the binding-scale and branch labels, not used as a raw
scalar replacement.

## Proof Target

Refine the current proof interface to:

```text
1. Keep exact q/Farey branch as theorem scale.
2. Inside p=2, use C=27 typed shell plus multiplicand clearance.
3. Inside p>=3, use Kpq/K33 incidence plus owner packets.
4. Inside relation-rich residuals, split visible folds from hidden balanced
   shells before applying additive-energy or Freiman bounds.
```

The next useful lemma is:

```text
Every remaining non-AP/GW q=14 atom either has a C27 typed-shell defect
handled by petal/lift rigidity, or has a three-owner K33 incidence packet
with enough sign-visible relation mass to feed the HYP-2908 tournament-state
lift.
```

This is not a proof.  It is the typed address that prevents the proof from
collapsing AP, shifted AP, Goddyn-Wong, `2/27` petals, and `3/41` near-misses
into the same scalar additive-energy bucket.

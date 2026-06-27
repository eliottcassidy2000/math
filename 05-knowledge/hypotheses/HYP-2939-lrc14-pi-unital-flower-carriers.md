---
id: HYP-2939
title: LRC14 pi/unital flower carrier synthesis
status: PROOF-INTERFACE / symbolic carrier guardrail; not a proof of LRC14
source: codex-2026-06-23-S137
related:
  - HYP-2938
  - HYP-2937
  - HYP-2936
  - HYP-2935
  - HYP-2934
  - HYP-2933
  - HYP-2932
  - HYP-2931
  - HYP-2930
  - HYP-2908
  - HYP-2823
  - HYP-2887
  - HYP-2891
  - THM-572
---

# HYP-2939: LRC14 Pi/Unital Flower Carrier Synthesis

The pi and unital prompts add useful carriers only if the quotient is named
before it is used.  They do not replace `q`, exact `M(S)`, or the current
C27/K33 branch labels.

This entry was rebased over the two concurrent S137 pi/unital guardrails
registered as `HYP-2938-lrc14-pi-unital-flower-alias-guardrail.md` and
`HYP-2938-pi-unital-flower-unit-guardrail.md`.  It keeps their normalization
warnings and adds the carrier-tournament proof interface.

## Pi Approximants

The script `04-computation/lrc14_pi_unital_flower_carriers_codex_s137.py`
stores output at
`05-knowledge/results/lrc14_pi_unital_flower_carriers_codex_s137.out`.

It computes:

```text
22/7 - pi          =  0.001264489267349619...
cuberoot(31) - pi  = -0.000212001198400234...
pi^3 - 31          =  0.006276680299820175...
```

Thus `cuberoot(31)` is numerically closer to `pi` than `22/7`, but `22/7`
is the one that explains the flower-family quotient.

## Flower Quotient

If each petal is placed after a turn of `1/pi` radians, then

```text
1/pi ~= 7/22.
```

Modulo one radian, this is a step `+7` on `Z/22`.  Since `gcd(7,22)=1`, the
rational model visits all `22` radian-residue families.  The exact turn has
only a small drift:

```text
22/pi - 7 = 0.002817496043394773...
```

The guardrail is that this is not a full-circle period.  The full-circle
rotation number is

```text
(1/pi)/(2pi) = 1/(2pi^2) = 0.050660591821168885...
```

The best full-circle rational with denominator at most `22` is `1/20`, and
substituting `pi=22/7` gives `49/968`, not a denominator-22 closure.  Therefore
the "22 families" claim is a radian-address quotient, not a literal closed
flower on the circle.

## Unital Distinctions

The geometric-combinatorics unital with `q=3` has parameters

```text
2-(28,4,1),   blocks=63,   point replication=9.
```

These numbers align with the LRC14 shell package:

```text
28 = 2*14,
27 = 28-1 = C = 2*14-1,
4  = q+1,
63 = 7*9,
31 = 27+4.
```

So `cuberoot(31)` is best read as a cubic-shell mnemonic:

```text
31 = q^3 + (q+1) = 27+4     for q=3.
```

This is separate from algebraic "unital" usage.  The algebraic lesson is a
unit-preserving quotient rule: a map from exact LRC data to a carrier should
send the identity/floor object to the carrier identity, rather than erasing
AP/GW, unit-visible C27 holes, or exact Farey branch information.  Unit groups
are a third notion: the `+7` step on `Z/22` and the unit/nonunit strata in the
`C=27` shell quotient are unit-group facts, not geometric unitals.

## Tournament Analysis

Tournament vertices are proof carriers, not runners.  The pairwise observable is

```text
(branch retention, typed visibility, unit preservation, pair coverage,
 state-lift fit, scalar-decoy resistance).
```

The resulting tournament is transitive:

```text
exact_M_Farey_branch
> marked_C27_shell_transfer
> bigraded_relation_signature
> Kpq_K33_incidence
> geometric_unital_28_4_1
> algebraic_unital_maps
> octa_Clebsch_packet
> PZ_degree4_gateway
> flower_22_radian_quotient
> cuberoot31_cubic_shell.
```

The ranking keeps the proof order intact.  Pi approximants are suggestive
carriers; exact Farey branch and C27 shell transfer remain the theorem-facing
data.

## Proof Target

The useful new target is a finite pair-coverage interface:

```text
After exact M/Farey branch and C27 transfer are attached, try to build a
q=3 unital-style pair-unique packet for the p>=3 K33 owner branch, and test
whether it can feed the HYP-2908 / THM-572 forbidden-H7 state lift.
```

This is deliberately weaker than an LRC14 proof.  It says that the q=3
geometric unital package may be a finite pair-coverage carrier for the
three-owner K33 branch, while `22/7` and `cuberoot(31)` remain labelled
quotient signals rather than scalar theorems.

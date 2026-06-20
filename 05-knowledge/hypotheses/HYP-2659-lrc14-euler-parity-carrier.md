---
id: HYP-2659
title: LRC14 Euler parity carrier - odd-shell carry plus mouth survival before scalar measure
status: OPEN; exact carrier scout
source: codex-2026-06-19-S36
depends_on:
  - THM-541
  - THM-543
  - THM-544
  - HYP-2651
  - HYP-2654
  - HYP-2656
  - HYP-2658
related:
  - HYP-2655
  - HYP-2648
  - HYP-2646
  - HYP-2630
  - HYP-2628
  - OPEN-Q-108
---

# HYP-2659 - LRC14 Euler Parity Carrier

## Claim

The Euler/Glaisher identity

```text
prod_{n>=1} (1+x^n) = prod_{m odd} 1/(1-x^m)
```

and the labelled suspension bijection between arbitrary graphs on `n-1`
vertices and even graphs on `n` vertices point to the same proof carrier:
free binary choices should be moved into a constrained parity-completed object
before scalarizing to measure.

For LRC14 AP-tail cores, the useful constrained object is the odd-shell carry
profile

```text
speed 2^a*m, m odd  ->  add 2^a units to odd shell m.
```

The working route is:

```text
odd-shell carry profile + exact drop-6 mouth survivor
    -> AP-tail/state-template routing
    -> scalar safe measure only at the final gate.
```

This does not prove LRC(14).  It refines the bounded near-collar proof
obligation by explaining why HYP-2656's odd base is rigid, why dyadic tails are
not independent strangers, and why the one-tail and two-tail layers behave
differently.

## Evidence

Script:

```text
04-computation/lrc14_euler_parity_carrier_codex_s36.py
```

Stored output:

```text
05-knowledge/results/lrc14_euler_parity_carrier_codex_s36.out
```

The script verifies, using exact integer or rational arithmetic:

1. `prod(1+x^n)=prod_odd 1/(1-x^m)` through `x^80`;
2. the binary carry bijection between distinct partitions and odd partitions
   for sample totals `14`, `20`, and `30`;
3. the labelled graph/even-graph suspension count `2^C(n-1,2)` for `n=3..7`;
4. exact LRC14 AP-collar safe measures, old-mouth survivors, and odd-shell
   carry profiles.

The key exact rows are:

```text
THM-541 drop-6 collar
  core=(1,2,3,4,5,7,8,9,10,11,12,13)
  safe=7/858
  old_survivor=7/858
  odd_carry_profile={1:15, 3:5, 5:3, 7:1, 9:1, 11:1, 13:1}

THM-543 one-tail exception
  core=(1,2,3,4,5,7,8,9,11,12,13,20)
  safe=3859/420420 = 7/858 + 1/980
  old_survivor=7/858
  delta_vs_drop6={5:+2}

THM-544 corrected two-tail minimum
  core=(1,2,3,5,7,8,9,11,12,13,20,46)
  safe=50189/3223220 = 426/35035 + 1571/460460
  old_survivor=1/364
  delta_vs_drop6={1:-4, 5:+2, 23:+2}
```

Reading: the unique one-tail below-second exception is pure carry inside an
already-present odd shell (`5`) and retains every old drop-6 mouth interval.
The corrected two-tail minimum introduces a new odd shell (`23`) and damages
the old mouth, but the exact scalar measure has already crossed the
AP-second threshold.

The AP second row is different again:

```text
drop 12 instead of drop 6
delta_vs_drop6={3:-2}
old_survivor=0
safe=426/35035.
```

So the second value is not merely "small scalar measure"; it is the first gate
where either old-mouth survival disappears or new-shell/mouth-damage debt has
already been paid.

## Relationship to HYP-2656

HYP-2656 shows on the single-deletion layer that the odd sub-core is constant
and the even speeds are dyadic refinements.  HYP-2659 extends that reading to
the AP-tail route: track dyadic carries by odd shell, then track which carries
preserve or damage the addressed drop-6 mouth components.

The slogan is:

```text
Glaisher carry accounts for speed content;
mouth survivor accounts for fixed-observer geometry.
```

Neither carrier alone is enough.  Carry without geometry cannot distinguish
harmless shell reinforcement from mouth damage; geometry without carry misses
why the near-collar rows are dyadic refinements rather than generic tails.

## Next Proof Obligation

Try to prove a parity-carrier rigidity lemma for bounded AP-tail/state-word
rows:

```text
If an LRC14 12-core has meas(G_C) < 426/35035,
then either
  (a) its odd-shell carry differs from the drop-6 profile only inside
      already-present shells and the four drop-6 mouths survive, or
  (b) the row pays at least 426/35035.
```

The exact one-tail and two-tail theorems are the first two layers of this
statement.  The remaining risk is deeper multi-tail/state-word damage or a
handoff to the HYP-2655 joint plateau/Delta recursion for genuinely wide rows.

## Tournament Analysis

Vertices are carrier lenses, not runners:

```text
odd_shell_carry
mouth_survivor
apex_parity_completion
state_word_damage
raw_speed_set
```

Pairwise observable: which lens preserves the predicate
`meas(G_C) < 426/35035` before scalarization.

Switch/gauge: classify by odd-shell carry and old-mouth survivor before using
total safe measure.

Hamiltonian path:

```text
odd_shell_carry > mouth_survivor > state_word_damage
> apex_parity_completion > raw_speed_set
```

Fingerprint: transitive proof-lens order in this scout; directed `3`-cycles
`0`.

Assumption challenge: tournament vertices do not have to be runners or arcs.
Here the relevant vertices are proof carriers.  The quotient preserves dyadic
speed content and old-mouth survival, but it destroys most fine endpoint
ownership away from the drop-6 mouth intervals.

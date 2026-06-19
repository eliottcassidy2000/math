---
id: HYP-2634
title: LRC(14) two-large lift-opposition atlas - bounded reciprocal signs are controlled by integer lift shape, not finite packet weight alone
status: PARTIALLY-CONFIRMED
source: codex-2026-06-19-S25
depends_on:
  - HYP-2633
  - HYP-2632
  - HYP-2614
  - THM-538
related:
  - HYP-2619
  - OPEN-Q-108
---

# HYP-2634 - LRC(14) Two-Large Lift-Opposition Atlas

## Claim

HYP-2633 found the first concrete lift-opposition warning: two `4+2` QR
finite packets with the same HYP-2632 weight `-25U` can have opposite bounded
reciprocal signs after integer lifting.  The example is:

```text
(1,2,8,9,15,22)   packet (1,1,1,1,2,2)   lift at H=16 positive
(1,4,8,11,15,22)  packet (1,1,1,1,4,4)   lift at H=16 negative
```

The next proof obligation is to decide whether this is an isolated accident or
a stable lift-shape phenomenon.  The proposed atlas should enumerate families
of integer representatives for the same finite two-large packet, compare their
bounded reciprocal signs, and record the low-height relation features that
drive sign opposition.

The working target is:

```text
finite packet class
+ integer lift profile
+ low-height relation-shell signature
=> predicted bounded reciprocal sign / cancellation pressure.
```

This is still a guardrail, not a proof of LRC(14).  Its purpose is to identify
the right discrepancy statistic for the HYP-2633 residue-lift
equidistribution / Abel-summation lemma.

## Computation

Script:

- `04-computation/lrc14_two_large_lift_opposition_codex_s25.py`
- output: `05-knowledge/results/lrc14_two_large_lift_opposition_codex_s25.out`

The script studies the seed family

```text
S_a = (1, a, 8, a+7, 15, 22),  a=2,3,4,5,6.
```

All rows are `4+2` finite packets.  HYP-2632 sees only the QR/NQR class:

```text
a=2,4      finite S/U = -25
a=3,5,6    finite S/U = -18
```

The bounded reciprocal lift through `H=16` is different:

```text
a  chi7  finite S/U        H8              H12              H16        sign
2   +1      -25      +0.0024787009  +0.0026649838  +0.0026761437   +
3   -1      -18      +0.0019992386  +0.0024456948  +0.0024241380   +
4   +1      -25      +0.0000033498  -0.0002201828  -0.0001305137   -
5   -1      -18      +0.0014682845  +0.0022436214  +0.0023090364   +
6   -1      -18      +0.0022698457  +0.0024658221  +0.0024711467   +
```

Only `a=4` flips negative in this seed family.

## Mechanism

The opposition is controlled by exact integer relation motifs, not by the
finite packet weight alone.  For `S_a=(1,a,8,a+7,15,22)`, the following
relation vectors have defects `sum_i S_a[i] n_i = c*a+d`:

```text
motif                  n                         defect    exact a    term z
universal_h2_positive  ( 1, 1,-1,-1,-2, 2)       0         2..6       +7.403726469e-05
a2_h2_negative         (-1, 2,-1,-1,-2, 2)       a-2       2          -3.135619413e-05
a4_h3_negative         (-1, 3,-1,-1, 2,-1)       2a-8      4          -7.233053203e-05
a4_h4_negative         (-4, 4,-1, 3,-1,-1)       7a-28     4          -3.655305842e-05
```

Thus every `S_a` has the universal positive height-2 relation, but the `a=4`
lift alone has the larger negative height-3/4 resonances.  The shell-group
contrast shows the cancellation:

```text
group       a=2 shell sum     a=2 cumulative    a=4 shell sum     a=4 cumulative
h=2        +0.0004290642     +0.0004290642     +0.0004465405     +0.0004465405
h=3..4     +0.0009962278     +0.0014252920     -0.0004119506     +0.0000345899
h=5..6     +0.0010907115     +0.0025160035     +0.0000631109     +0.0000977008
h=8..12    +0.0001489803     +0.0026649838     -0.0003178836     -0.0002201828
h=13..16   +0.0000111599     +0.0026761437     +0.0000896690     -0.0001305137
```

This makes the hidden integer sequence explicit: the dangerous lift is detected
by low-height relation defects such as `2a-8` and `7a-28`, not by the
Legendre class of `a` alone.

## Offset Scan

The script also scans consecutive-ladder variants where the residue-1 core is
four consecutive lifts and the repeated QR pair is two consecutive lifts.  At
height `H=12`, sign opposition is localized:

```text
core_start  pair_start       a=2 H12       a=4 H12    opposite?
0           0             +0.002664984   -0.000220183    True
0           1             +0.001283113   +0.002632105    False
0           2             +0.002282046   +0.001892945    False
1           0             +0.001684260   +0.001670322    False
1           1             +0.002227558   +0.002971846    False
1           2             +0.002313445   +0.002254372    False
```

Moving the repeated pair one period up removes the `a=4` low-height resonance.

## Tournament Analysis

Candidate vertices included runners, residues, integer lift offsets,
relation-shell events, boundary faces, additive-frequency shells, and proof
obligations.

Chosen Hamiltonian path:

```text
low_height_relation_motifs
> integer_lift_offsets
> residue_lift_discrepancy
> additive_frequency_shells
> finite_packet_weight
> boundary_face_mass
> raw_support_values
> raw_runner_vertices
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles = 0
SCC_sizes = [1,1,1,1,1,1,1,1]
hamiltonian_paths = 1
```

The quotient preserves the bounded signed reciprocal-tail predicate.  It
destroys raw witness times and finite packet labels.

## Status

Partially confirmed.  The first opposition pattern is now explained by exact
low-height integer relation defects.  The next step is to broaden the relation
motif sieve beyond the start-aligned seed family and turn it into the
discrepancy statistic used in HYP-2633's Abel-summation lemma.

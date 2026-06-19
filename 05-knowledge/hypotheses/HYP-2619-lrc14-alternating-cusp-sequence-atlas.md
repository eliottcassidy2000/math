---
id: HYP-2619
title: LRC(14) alternating cusp sequence atlas
status: OPEN
source: codex-2026-06-19-S15
depends_on:
  - THM-538
  - HYP-2608
  - HYP-2614
  - HYP-2615
  - HYP-2617
related:
  - HYP-2618
  - HYP-2616
  - HYP-2612
  - HYP-2613
  - MISTAKE-078
  - OPEN-Q-108
---

# HYP-2619 - LRC(14) Alternating Cusp Sequence Atlas

## Claim

The "large absolute mass but tiny signed mass" phenomenon in the LRC(14)
support-6 residual is not a numerical accident.  It is an alternating-series
structure visible in three compatible sequence layers:

1. conjugation-paired residue coefficients,
2. shell-by-shell signed reciprocal sums on exact relation hyperplanes,
3. projective coimage fibers after reducing speed residues mod `7`.

The proof target should therefore be a two-stage signed theorem:

```text
raw cusp absolute mass
  -> signed shell variation
  -> alternating shell net / coimage fiber sum,
```

after finite low-height wall deletion.  A proof that treats all signs as
positive attacks the wrong object; it is the analogue of replacing a
conditionally convergent alternating series by a divergent all-positive one.

This does not prove LRC(14).  It gives the next analytic shape: finite wall
ledger plus class-by-class Dedekind/cotangent summation over non-null coimage
fibers, with explicit warning that monotone decay of the largest coimage fiber
is false beyond `d=13`.

## Computation

Script:

- `04-computation/lrc14_alternating_cusp_sequences_codex_s15.py`
- output: `05-knowledge/results/lrc14_alternating_cusp_sequences_codex_s15.out`

The script reuses HYP-2614's exact residue coefficient `C_d(r)` and HYP-2617's
projective coimage classes.  Individual residue coefficients are complex, so
the first signed real atoms are conjugate pairs:

```text
r <-> -r mod 7.
```

This is the finite residue version of the CM/conjugation quotient: the signed
object appears only after pairing.

## Residue Sign-Balance Sequence

For `d=6..16`, the conjugation-paired residue totals have comparable positive
and negative piles:

```text
d      +pairs  -pairs        net             abs total       abs/net
6      12088   11240   ~0                  111.806          inf
7      12088   11240   ~0                   65.50005        inf
8      12088   11240   ~0                   38.32543        inf
9      11458   11870   ~0                   22.69421        inf
10     11683   11645   ~0                   14.31111        inf
11     10493   12835    0.01365212           9.92566        727.04172
12     10103   13225    0.04095636           7.158129       174.77455
13     10085   13243    0.0741115            5.297627        71.481848
14     10565   12763    0.1053163            4.004875        38.027096
15     10565   12763    0.1297946            3.121947        24.052977
16     10805   12523    0.1457951            2.458847        16.865088
```

The max paired coefficient sequence is:

```text
0.286199688, 0.143099844, 0.0749570612, 0.0545142264,
0.0385215324, 0.0269789794, 0.0189053955, 0.0133196091,
0.00945915589, 0.00739739336, 0.0057554252
```

Readout: the exact zero net through `d=10`, followed by a tiny signed drift
from `d=11`, is the sequence-level version of the user's alternating-series
diagnosis.  The all-positive pile is large and smooth; the signed object is
small, structured, and quotient-dependent.

## Shell Alternation Ladders

For the five named S12/S14 support rows, the raw absolute mass splits into
cusp-to-shell cancellation and shell-to-net cancellation:

```text
case                 raw/variation   variation/net   raw/net
AP core              29.0667         1               29.0667
resonant 21          103.843         2.0898461       217.015
dissociated 211      45.2518         1.287521        58.2626
k=9 wide 68          53.8882         20.760695       1118.7569
k=10 wall 22         44.8267         9.9857761       447.6290
```

The AP core has no shell sign flips; it needs only cusp-to-signed-shell
collapse.  The k=9 wide support and k=10 wall support are different: they need
the second alternating-shell summation lemma as well.  This distinguishes the
proof obligations instead of throwing all supports into one absolute envelope.

Shell sign words:

```text
AP core              ++++++++++++++++++++++++
resonant 21          -+---++----+-+-+-+-+--+
dissociated 211      ---++++++---+--------
k=9 wide 68          +--++++--+-+++++------+
k=10 wall 22         ++-+++++-+---+++--+--
```

## Boundary-Face Sequences

Final-shell relation counts are large:

```text
5839636, 1991422, 378826, 854792, 1113506
```

The final-shell abs/signed ratios are:

```text
211.94099, 9220.9535, 13.440976, 34.887488, 116.21417
```

One-face boundary relation counts:

```text
5291542, 1825394, 353652, 783346, 996976
```

One-face abs/signed ratios:

```text
202.38572, 3484.6403, 13.323897, 34.803376, 112.73209
```

Readout: integer cusp counts are pre-quotient mass.  The proof-relevant
fraction is the signed collapse ratio, not the absolute number of boundary
relations.

## Coimage Extension

Extending the HYP-2617 projective coimage atlas to `d=16` gives:

```text
max |S_d| sequence:
2.18524924, 1.40480308, 0.825043079, 0.439598243, 0.278503234,
0.264266024, 0.221554394, 0.210632698, 0.259589278, 0.291590164,
0.307196067

coimage-null count sequence:
142, 113, 80, 43, 34, 3, 3, 3, 3, 3, 3

coimage <0.01 count sequence:
142, 113, 80, 77, 43, 35, 40, 25, 37, 18, 19
```

Important correction: the largest coimage fiber does not keep decreasing after
`d=13`; it rebounds at `d=14,15,16`, with the all-ones class becoming dominant.
Thus the final analytic theorem should not rely on monotone max-fiber decay.
It should be class-by-class and sign-aware.

## Tournament Analysis

The script explicitly rejects runner vertices as the default.  It considered
runners, residue tuples, shell heights, boundary faces, Fourier modes, coimage
classes, and proof obligations.  The chosen vertices are proof-obligation
sequence families:

```text
raw_absolute_volume
coefficient_sign_inventory
shell_alternation_ladder
boundary_cusp_faces
projective_coimage_nullity
height_one_wall_ledger
dedekind_tail_bound
```

The tournament is transitive, with Hamiltonian proof path:

```text
dedekind_tail_bound
> projective_coimage_nullity
> shell_alternation_ladder
> height_one_wall_ledger
> boundary_cusp_faces
> coefficient_sign_inventory
> raw_absolute_volume
```

This quotient preserves the signed-tail predicate and the recursive `d/k`
sequence data while discarding witness-time geometry.

## Status

LRC(14) remains open.  HYP-2619 improves the route by separating the support-6
residual into:

1. finite low-height wall rows,
2. cusp-to-signed-shell collapse,
3. alternating shell summation,
4. finite projective coimage fiber summation.

The next theorem to try is a signed reciprocal-tail estimate over each non-null
projective coimage class after the HYP-2616-style wall ledger is removed.  This
is the support-six analogue of HYP-2618's packet-address lesson: retain the
finite address first, then evaluate the signed compatible mass.

---
id: HYP-3041
title: LRC14 AP-tail q13 puncture and reciprocal fixed-point atlas
status: PROVED family lemma / proof-interface atlas; not a proof of LRC14
source: codex-2026-06-26-S204
tangent: T1122
related:
  - HYP-3033
  - HYP-3032
  - HYP-3031
  - HYP-3029
  - HYP-3028
  - HYP-3027
  - HYP-3024
  - HYP-3023
  - HYP-3017
  - HYP-2963
  - OPEN-Q-108
---

# HYP-3041: AP-Tail Puncture Atlas

## Claim

For the AP-tail family

```text
S_m = {1,2,...,12,m},
```

the LRC14 threshold is discharged by two explicit clocks:

```text
m not divisible by 13:  t = 1/13
m = 13s, s >= 2:        t = s/(13s+1)
m = 13:                 AP boundary atom
```

Thus every `S_m` with `m != 13` has a strict lonely witness above `1/14`.

The two coarse largest-stalk residuals from HYP-3029 are exactly instances
where the quotient retained the same mod-14 owner strip and coarse local
geometry, but forgot the `13 | m` puncture clock.

## Proof

Let `S_m={1,...,12,m}`.

If `13` does not divide `m`, take `t=1/13`.  For every `j=1,...,12`,
`j/13` is nonzero modulo `1` and has distance at least `1/13` from an integer.
The tail `m/13` is also nonzero modulo `1`, so it too has distance at least
`1/13`.  Hence

```text
min_{v in S_m} ||v t|| >= 1/13 > 1/14.
```

If `m=13s` with `s>=2`, take

```text
t = s/(13s+1).
```

Runner `1` has distance `t`.  The tail runner satisfies

```text
m t = 13s * s/(13s+1) = s - s/(13s+1),
```

so its distance to an integer is also `t`.

For `2<=j<=12`, the point `j t` lies between `0` and `1`, and its distance is
either `j t` or `1-j t`.  The first is at least `t`; the second is smallest
at `j=12` and equals

```text
1 - 12s/(13s+1) = (s+1)/(13s+1) > s/(13s+1) = t.
```

Therefore every AP-core runner has distance at least `t`, while runners `1`
and `m` bind at height `t`.  Finally,

```text
t > 1/14  <=>  14s > 13s+1  <=>  s > 1.
```

The lone non-strict multiple is `s=1`, i.e. `m=13`, the AP equality atom.

## Computation

Script:

```text
04-computation/lrc14_ap_tail_puncture_atlas_codex_s204.py
```

Stored output:

```text
05-knowledge/results/lrc14_ap_tail_puncture_atlas_codex_s204.out
```

The script checks `m=13..400` and finds:

```text
counts={'ap_boundary': 1, 'q13_puncture': 358, 'fixed_point_tail': 29}
non_strict_except_AP=0
min_strict_margin=1/378
```

The minimum strict margin occurs at `m=26`, where `t=2/27`.

## HYP-3029 Residual Teeth

HYP-3029 left exactly two mixed-route fibers after the coarse largest-component
stalk on the hard AP/GW automatic word `MFCMMCCFFFCCC`.  HYP-3033 names the
neighboring stored-ledger problem as residual certificate scheduling; this file
discharges the two AP-tail teeth by family formula.

The S204 atlas identifies them as:

```text
13->104: m=104, mod14=6, mod13=0, fixed_point_tail, t=8/105
13->118: m=118, mod14=6, mod13=1, q13_puncture,    t=1/13

13->117: m=117, mod14=5, mod13=0, fixed_point_tail, t=9/118
13->159: m=159, mod14=5, mod13=3, q13_puncture,    t=1/13
```

So the collision is not a new F7 sector.  The coarse stalk remembered the
mod-14 owner strip:

```text
('1+') -> ('1+','6-') -> ('6-')
('1+') -> ('1+','5-') -> ('5-')
```

but forgot whether `m` kills the `q=13` witness.  Reattaching the
`13|m` bit, or the explicit AP-tail certificate, makes the target-fiber
coarse-stalk ladder route-pure:

```text
residue_terminal              mixed_route=27
coarse_stalk                  mixed_route=2
coarse_plus_q13_bit           mixed_route=0
coarse_plus_tail_certificate  mixed_route=0
coarse_plus_peak_height       mixed_route=0
exact_stalk                   mixed_route=0
magnitude_cocycle             mixed_route=0
```

## HYP-3031 Repair Class

In HYP-3031 language, these residuals are `nested_refinement` with an
`owner_strip` witness:

```text
mod-14 owner strip      row/column shadow
q=13 puncture bit       missing mixed coordinate
fixed-point tail clock  reciprocal owner-strip repair
AP-tail formula         nested family descent
```

The quotient may forget exact stalk height only after the AP-tail formula has
been proved; otherwise it erases the mixed coordinate that distinguishes a
direct `q=13` proof from a fixed-point covering proof.

## Tournament Analysis

Vertices are repair clocks and proof carriers, not runners:

```text
residue_terminal
coarse_stalk
coarse_plus_q13_bit
coarse_plus_tail_certificate
coarse_plus_peak_height
exact_stalk
magnitude_cocycle
```

Pairwise observable:

```text
route_purity,
max_mixed,
exact_period_retention,
owner_strip,
local_geometry,
family_proof_scope,
proof_cost
```

The tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

The first Hamiltonian path is:

```text
coarse_plus_tail_certificate
> coarse_plus_q13_bit
> coarse_plus_peak_height
> exact_stalk
> magnitude_cocycle
> coarse_stalk
> residue_terminal
```

## Next Use

Search the remaining automatic-fiber residuals for rows that are AP-core plus
one or two tail punctures.  Before invoking heavy Fejer, topology, or THM-572
machinery, try to prove a reciprocal fixed-point certificate of the form

```text
t = source_multiplier / (tail + source_shift).
```

This turns HYP-3031's abstract `zeta_repair_class` into a practical rule:
when a quotient remembers the owner strip but forgets a prime clock, attach
the missing exact-period puncture or prove the reciprocal fixed point.

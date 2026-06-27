---
id: HYP-3088
title: LRC14 Conjecture 7.1 bridge is refuted in raw time and repaired in normalized ruler coordinates
status: REFUTED-AND-REPAIRED / proof-interface target
source: codex-2026-06-27-S255
related:
  - THM-574
  - THM-573
  - THM-566
  - THM-565
  - THM-530
  - HYP-3087
  - HYP-2758
  - HYP-2072
  - OPEN-Q-108
external:
  - arXiv:2604.23906
---

# HYP-3088: Conjecture 7.1 Refuted, Normalized-Arc Route Repaired

## Correction

The earlier HYP-3088 formulation identified the paper's Conjecture 7.1 with a
uniform largest-lonely-arc floor in the original time circle. That identification
is false.

THM-574 proves that Conjecture 7.1 is false for `k=13` as stated. For every
proposed denominator bound, the divisor-loaded rows

```text
S_B = {1,2,...,11,13,84*lcm(1,...,B)}
```

are primitive and non-tight, but have no witness in `(1/d)Z` for any `d <= B`.
The direct safe components necessarily shrink like `O(1/N_B)` because the
loaded apex alone has safe gaps of length `6/(7N_B)`.

So the paper's witness-time conjecture is not the project witness route
verbatim.

## Repaired bridge

The useful bridge is instead normalized:

```text
paper k=13, k+1=14=2*7
  -> prime-field Proposition 4.1 no longer applies
  -> paper fallback wants c=2 and c=7 lifts
  -> project has dyadic descent + THM-573 level-7 sieve
  -> residual is primitive covering rows with <=6 multiples of 7
  -> prove a normalized slow/ruler-coordinate arc floor there
  -> use THM-565 finite-ruler sampling to return to real witnesses.
```

In this coordinate, divisor loading is not a refutation. It is exactly the
reason to avoid raw time arcs: the large speed makes direct components tiny,
while the slow coordinate keeps the finite-ruler object at bounded scale.

## Incoming S61 integration

The concurrent mac-mini S61 checkpoint confirms both sides of this repair.
Its exact `I(k,p,1)` apex bridge reports:

```text
no (1/7)-witness  <=>  some coordinate is 0 mod 7
one c=2 lift:
no (1/14)-witness <=>  some coordinate is 0 mod 14
```

So the paper's first improper-set bottleneck really does descend to the
project's covering condition: odd multiples of `7` are rescued by the
`1/14` ansatz, while multiples of `14` survive as covering rows needing a
finer/continuous witness.

S61 also measures the direct lonely set for `{1,...,12,14V}`. Direct arcs are
bounded below only in the finite `V <= V*` range and then get shredded by the
apex at scale `1/V`. This is the same failure mode THM-574 exposes with
divisor-loaded rows. Read the S61 "crossover near `V*`" as a finite-ruler
threshold and a cue to peel the apex, not as a rescue of the raw absolute
denominator Conjecture 7.1.

## Two live proof obligations

**O1: normalized arc floor in the THM-573 residual.** After THM-573, every row
with at least seven multiples of `7` is closed. The live residual has at most
six multiples of `7`, including the divisor-loaded rows that refute raw
Conjecture 7.1. The theorem target is:

```text
S = P union {V-e : e in E}, with V large and residual <=6 mod-7 multiples
  -> G(P,E) has measure c >= c0 and component count m <= m0
  -> THM-565 gives a finite-ruler witness for V > m0/c0
  -> finite complement closes by exact check / packet ledger.
```

This keeps THM-530's `m_P=14249/252252` and THM-565's arc-count technology,
but it must be stated in the normalized slow coordinate, not as a direct
Conjecture 7.1 denominator claim.

**O2: initial improper-set sieve `I(k,p,1)`.** The paper names computing
`I(k,p,1)` as the main k=13 bottleneck. HYP-2072 already gives the matching
algorithmic shape for `k+1=14=2*7`: a mod-7 oracle settles most residues, an
apex-stuck peel removes most of the residual, and the genuine hard ratio-cover
part is small. This should be integrated with the normalized-arc certificate:
the sieve should retain the coordinate that certifies normalized witness mass,
not ask for a raw absolute denominator.

## Back-and-forth readout

Failure of O1 in raw time explains how to formulate O2 correctly. The
divisor-loaded family sits inside the THM-573 residual (`mult7_count=2` in the
scout), so any `I(k,p,1)` pruning rule based on small absolute denominators is
guaranteed to lose. Conversely, O2 suggests what O1 must preserve: the
`c=2,7` lift residue word and the normalized finite-ruler arc, because those
survive divisor loading.

## Tournament Analysis

Vertices are proof obligations / coordinate carriers, not runners. Alternate
vertex sets considered and rejected as primary vertices: runners, raw arcs,
fixed circle sections, residues mod `14`, residues mod `7`, denominators,
danger endpoints, lift factors, Fourier modes, and packet rows.

Pairwise observable: which carrier retains the LRC predicate under legal
scaling/descent, survives divisor loading, reduces the paper's `I(k,p,1)` or
`c=14` bottleneck, and avoids a false denominator assumption.

The S255 scout ranks:

```text
normalized_ruler_arc_floor
  > level7_c7_lift_sieve_THM573
  > two_tier_CRT_I_k_p_1_sieve
  > THM566_divisor_loaded_obstruction
  > direct_t_largest_arc_floor
  > prime_field_Prop4_1_shortcut
  > raw_Conjecture_7_1_denominator_floor
```

## Evidence

`04-computation/lrc14_conj71_refutation_normalized_arc_codex_s255.py` records:

- exact non-tight witnesses for `B=6,10,14,26,41,67`;
- no denominator hits up to `B` in those rows;
- direct component shrinkage: `B=6` has largest direct component `1/5880`;
- the tight AP row has no positive components, while `{1..11,13,84}` has
  positive direct components but already shows apex-scale shrinkage.

The hypothesis is therefore no longer "Conjecture 7.1 implies LRC14." It is:
prove the normalized version that the project already needs, and use the
paper's `c=2,7` lift language to sharpen the residual ledger.

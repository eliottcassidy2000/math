---
id: HYP-3452
title: LRC14 AP84 tail component phase
status: EVIDENCE / exact AP-tail component-cover phase audit; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3433 endpoint-spine law and HYP-3450/HYP-3451 component-cover graph router
tangent: T1412
technique: LTI-412
tournament_technique: LTT-312
script: 04-computation/lrc14_ap84_tail_component_phase_codex_20260629.py
result: 05-knowledge/results/lrc14_ap84_tail_component_phase_codex_20260629.out
reflection: 07-reflections/lrc14-ap84-tail-component-phase-codex-20260629.md
related:
  - HYP-3457
  - HYP-3456
  - HYP-3454
  - HYP-3453
  - HYP-3451
  - HYP-3450
  - HYP-3436
  - HYP-3435
  - HYP-3434
  - HYP-3433
  - HYP-3431
  - HYP-3429
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3422
  - HYP-3417
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3452: LRC14 AP84 Tail Component Phase

## Claim

HYP-3433 found that the canonical tail

```text
S_m = {1,2,...,11,13,84m}
```

has an eventual harmonic endpoint-spine law from `m >= 5`.  HYP-3450 and
HYP-3451 then showed that AP-with-`84m` rows are the component-cover graph
danger family: connected dead-cover projections, high dead fraction, and few
low-rank escapes.

HYP-3452 stitches those two observations together.  In the checked range
`m=1..70`, the component-cover graph has the same phase transition:

```text
m >= 5:
  best component is rank 1
  endpoint labels are E:84m/E:84m
  best interval is [(14ceil(48m/7)+1)/(588m),
                    (14ceil(48m/7)+13)/(588m)]
  best length is 1/(49m)
```

The transient is finite and sharper than before:

```text
m=1,2: max paired dead-cover rank = 3
m>=3:  max paired dead-cover rank <= 2
m=1..4: best escape is mixed E:84m/B1:5
m>=5:  best escape is pure E:84m/E:84m
```

Thus the AP-tail graph proof can split into the HYP-3457 finite packet for
`m=1..4` plus a rank-one harmonic tail for `m>=5`.

## Exact Readout

Script:

```text
04-computation/lrc14_ap84_tail_component_phase_codex_20260629.py
```

Stored result:

```text
05-knowledge/results/lrc14_ap84_tail_component_phase_codex_20260629.out
```

Readout:

```text
canonical_family=S_m={1,2,...,11,13,84m}, checked_m=1..70
rank_one_EE_endpoint_phase_start=5
hyp3433_interval_failures_after_phase=[]
dead_pair_rank_le_2_start=3
connected_dead_projection_start=1
low_rank_escape_equals_alive_failures=[]
max_dead_fraction=361/373 at m=35
```

The dead projection is connected throughout the checked range.  This keeps the
HYP-3451 graph theorem target alive but makes the base family cleaner: the
largest local cover rank disappears by `m=3`, and the rank-one endpoint-spine
certificate begins by `m=5`.

## Escape Count Clock

The number of escape components is not constant.  In the checked range it
obeys a period-`35` Beatty correction:

```text
escapes(m) = 2 * (floor(12m/35) + d[m mod 35])
```

with residues listed as `m=1..35`, final entry representing residue `0`:

```text
d = [2,2,1,1,1,1,1,2,1,1,2,1,1,1,1,2,2,1,1,1,1,2,2,1,1,2,1,1,2,1,2,2,1,1,0].
```

This explains why raw dead fraction is a poor proof guide.  HYP-3451's smaller
bank saw the maximum dead fraction at `m=5`; extending the AP-tail family to
`m=70` moves the checked maximum to `m=35`, because residue `0 mod 35` minimizes
the escape correction.  The theorem-relevant event is not that scalar peak.
It is the `m=5` transition to a rank-one `E:84m` endpoint component.

## Candidate Lemma

For every `m>=5`, prove directly that

```text
I_m = [(14ceil(48m/7)+1)/(588m),
       (14ceil(48m/7)+13)/(588m)]
```

is an `E_safe` component for the tail half-speed `42m` and survives the odd
two-branch bad cover.  Equivalently, prove that the HYP-3433 endpoint interval
is not only an endpoint-spine winner but also the HYP-3450 component-cover
escape.

HYP-3457 proves the small transients `m=1..4` separately, where the best
component is still mixed against `B1:5`.  The paired-cover graph side of the
proof should use:

```text
dead projection connected for all checked m
max paired dead-cover rank <= 2 for m>=3
mod-35 escape count clock for the remaining component boundary count
```

HYP-3456 supplies the all-`m` floor count over the HYP-3431 low corridors cut
by the moving tail period.  The denominator `35` is the residue clock produced
by the low corridor endpoints against `42m`; it is now a finite floor identity
rather than an analytic density.

## Proof Use

This packet improves the HYP-3451 next hook.  Instead of trying to prove the
AP-with-`84` base case and then vaguely "lift to AP-with-`84m`", the proof can
aim for a three-part AP-tail theorem:

```text
finite transient: HYP-3457 packet for m=1..4
rank-one harmonic tail: m>=5
residue boundary count: mod-35 Beatty correction
```

That gives a more legal quotient.  It retains endpoint labels, component
escape status, paired-cover rank, and dead-projection connectivity.  It does
not rely on raw dead fraction, which is now known to peak at a residue artifact
inside the checked AP-tail range.

## Tournament Analysis

Tournament vertices are proof carriers for the AP-tail phase, not runners.

```text
score_hist={20:1,23:1,49:1,51:1,55:1,56:1,58:1}
directed_3cycles=0
hamiltonian_path=
  m5_rank_one_endpoint_phase
  -> H3433_tail_address_certificate
  -> connected_dead_projection_family
  -> mod35_beatty_escape_clock
  -> m3_dead_pair_rank_drop
  -> dead_fraction_peak_scalar
  -> raw_component_count_slope
```

## Assumption Challenge

Alternate vertices considered:

```text
runners, m-values, E_safe components, dead components, branch blockers,
endpoint labels, tail gaps, and residues mod 35.
```

The chosen quotient preserves the branch-union escape predicate, endpoint
rank/labels, dead-cover rank, and dead-projection connectivity.  It destroys
full wall geometry away from the AP-tail family, so it is a base-case and
induction carrier, not a global proof of LRC14 by itself.

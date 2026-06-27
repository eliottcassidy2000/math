---
id: HYP-3018
title: LRC14 active-bottleneck normal-fan carrier
status: EVIDENCE / exact finite scout and proof-interface carrier; not a proof
source: codex-2026-06-26-S182
tangent: T1101
script: 04-computation/lrc14_active_bottleneck_normal_fan_codex_s182.py
result: 05-knowledge/results/lrc14_active_bottleneck_normal_fan_codex_s182.out
related:
  - HYP-3017
  - HYP-3016
  - HYP-3015
  - HYP-3014
  - HYP-3013
  - HYP-3012
  - HYP-3011
  - HYP-3010
  - HYP-3009
  - HYP-3008
  - HYP-2996
  - HYP-2975
  - HYP-2974
  - HYP-2963
  - THM-572
  - LTI-166
  - LTI-165
  - LTI-164
  - OPEN-Q-108
---

# HYP-3018: LRC14 Active-Bottleneck Normal-Fan Carrier

## Claim

A better proof carrier than a raw barcode is the active bottleneck normal fan
of the lonely profile

```text
m_S(t)=min_v ||v t||.
```

For each strict safe bar, keep the left endpoint owner set, peak bottleneck
owner set, right endpoint owner set, support speeds, residue sums mod `14`, and
peak time/height.  These records expose the local inequalities that a Fejer,
endpoint-owner, tope/cocircuit, or state-lift certificate must discharge.

The carrier preserves the LRC predicate because a positive bar is still a
strict lonely interval.  It retains more proof data than HYP-3015's barcode and
repairs the HYP-3016/HYP-3017 warning that automatic sidecars can mix boundary
and open rows inside the same quotient fiber.

## Computation

Script:

```text
04-computation/lrc14_active_bottleneck_normal_fan_codex_s182.py
```

Stored output:

```text
05-knowledge/results/lrc14_active_bottleneck_normal_fan_codex_s182.out
```

The script uses exact Fraction arithmetic on triangular-wave cells.  It
computes global maximin peaks, strict threshold bars at `1/14`, and owner signs
at bar endpoints and peaks.

Named-row readout:

```text
AP13 and GW 12->24:
  zero bars, M=1/14, six global boundary peaks
  peak owner pairs include {1+,13-}, {5+,9-}, {3+,11-}
  every displayed boundary peak pair has residue sum 0 mod 14

K33 12->36:
  two bars, safe_mass=1/1260
  peak support (5,36), peak owner residue sum 13

petal 10->20 and P10+GW:
  two bars with peak support (7,20)

petal 13->26:
  two bars with peak support (1,26)

covering 12->84:
  eight bars; top support (5,84), lowest-persistence support (13,84)

fibbinary first13:
  38 bars; peak owner count histogram {2:28, 3:6, 4:4}

Moser-de Bruijn first13:
  64 bars; peak owner count histogram {2:56, 3:6, 4:2}
```

Aggregate exact fingerprints over the tested rows:

```text
total_positive_bars=138
peak_owner_count_hist={2:114, 3:18, 4:6}
endpoint_owner_count_hist={1:276}
bar_support_size_hist={2:86, 3:42, 4:8, 5:2}
bars_with_endpoint_zero_sum=8
```

The low-persistence discharge queue starts with Moser and covering bars:

```text
Moser support (17,80), persistence 1/1358
Moser support (16,81), persistence 1/1358
covering support (13,84), persistence 1/1358
fibbinary support (20,21), persistence 1/574
fibbinary support (9,32), persistence 1/574
```

## Mixed Automaton-Fiber Repair

HYP-3016/HYP-3017 show that automaton and residue-language fibers can mix
AP/Goddyn-Wong boundary atoms with open rows.  The normal fan separates those
inside the same coarse fiber:

```text
AP13:
  zero_bar boundary_peak_owners={1+,13-}, {5+,9-}, ...

mixed_AP_fiber_12_to_26:
  open_bar peak=5/12 support=(5,7)

mixed_AP_fiber_12_to_96:
  open_bar peak=42/101 support=(5,96)

GW 12->24:
  zero_bar boundary_peak_owners={1+,13-}, {5+,9-}, ...

mixed_GW_fiber_12_to_38 and 12_to_52:
  open_bar peak=5/12 support=(5,7)
```

So the better carrier is not "automaton word plus M"; it is the active local
support that tells which inequalities form the witness or boundary atom.

## Tournament Analysis

Vertices are proof carriers and local certificate schemas, not runners:

```text
active_bottleneck_normal_fan
lonely_profile_persistence_barcode
tope_cocircuit_wall_language
fejer_interval_certificate
exact_safe_component_front
automaton_fiber_magnitude_cocycle
raw_maximin_scalar
raw_automaton_word
```

Pairwise observable:

```text
(predicate, local_bottleneck_support, endpoint_geometry,
 quotient_purity_repair, stability, handoff_quality, anti_scalar)
```

Switch/gauge: majority comparison of the observable vector, with the listed
carrier order as the tie Hamiltonian path.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
score_order=active_bottleneck_normal_fan > lonely_profile_persistence_barcode > tope_cocircuit_wall_language > fejer_interval_certificate > exact_safe_component_front > automaton_fiber_magnitude_cocycle > raw_maximin_scalar > raw_automaton_word
```

## Packet Fields

Add a sidecar group to HYP-2963-style packet records:

```text
normal_fan_left_owner_set
normal_fan_peak_owner_set
normal_fan_right_owner_set
normal_fan_support_speeds
normal_fan_support_size
normal_fan_owner_residue_sums
normal_fan_endpoint_zero_sum_flag
normal_fan_peak_time
normal_fan_peak_height
normal_fan_local_equation_count
```

These fields should travel with exact `M`, barcode bars, automaton sidecars,
endpoint-owner labels, Fejer interval centers, and state-lift route labels.

## Assumption Challenge

Alternate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid
circuits, persistence bars, active bottleneck sets, and proof obligations.

The chosen quotient preserves the LRC predicate "there is an open safe bar
above `1/14`" and retains the local active constraints that make the witness
checkable.  It destroys global runner identity only after the active support,
endpoint owners, exact peak, and route labels are retained.  A future quotient
may collapse normal-fan classes only after proving boundary/open status is
constant on the collapsed class or emitting a named residual.

## Next Pulls

1. Add the normal-fan sidecar to the full HYP-2963 packet bank.
2. Test whether each automatic mixed fiber becomes pure after adding peak
   support and endpoint residue sums.
3. Use lowest-persistence support pairs as the priority list for interval
   Fejer and endpoint-owner certificates.
4. Compare peak support pairs such as `(5,36)`, `(7,20)`, `(1,26)`,
   `(13,84)`, `(16,81)`, and `(17,80)` against C27/K33/state-lift labels.

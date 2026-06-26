---
id: HYP-3015
title: LRC14 lonely-profile persistence barcode carrier
status: EVIDENCE / exact finite scout and proof-interface carrier; not a proof
source: codex-2026-06-26-S179
tangent: T1099
script: 04-computation/lrc14_lonely_profile_persistence_barcode_codex_s179.py
result: 05-knowledge/results/lrc14_lonely_profile_persistence_barcode_codex_s179.out
related:
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
  - LTI-164
  - LTI-163
  - OPEN-Q-108
---

# HYP-3015: LRC14 Lonely-Profile Persistence Barcode Carrier

## Claim

The next useful creative angle is to stop treating exact safe components as a
flat yes/no or mass scalar.  For a row `S`, keep the full threshold-superlevel
profile

```text
m_S(t) = min_{v in S} ||v t||
```

and record each connected component of `{t : m_S(t) > 1/14}` as a persistence
bar with exact length, peak time, peak height, and persistence margin above
`1/14`.

This is not a proof of LRC14.  It is a proof-interface carrier: a primitive
row is certified open only when it has at least one positive bar, while AP and
Goddyn-Wong remain zero-bar boundary atoms.  The barcode keeps more information
than raw `M`: it records topology, local robustness, and where a certificate
can anchor.

## Computation

Script:

```text
04-computation/lrc14_lonely_profile_persistence_barcode_codex_s179.py
```

Stored output:

```text
05-knowledge/results/lrc14_lonely_profile_persistence_barcode_codex_s179.out
```

The script cuts the time circle at all triangular-wave breakpoints `k/(2v)`.
On each cell, every function `||v t||` is affine, so the maximum of
`min_v ||v t||` and the exact threshold-safe intervals are computed with
Fraction arithmetic.  No grid sampling is used.

Named-row readout:

```text
AP13:               M=1/14, bar_count=0, safe_mass=0
GW 12->24:          M=1/14, bar_count=0, safe_mass=0
K33 12->36:         M=3/41, bar_count=2, safe_mass=1/1260
petal 10->20:       M=2/27, bar_count=2, safe_mass=1/980
petal 13->26:       M=2/27, bar_count=2, safe_mass=1/182
covering 12->84:    M=7/89, bar_count=8, safe_mass=563/105105
fibbinary first13:  M=3/25, bar_count=38, safe_mass=66077/399840
Moser first13:      M=1/6,  bar_count=64, safe_mass=4264747/40348854
```

The smallest positive persistence among the tested rows occurs in the covering
and Moser bar families:

```text
covering 12->84:       min persistence = 1/1358
Moser-de Bruijn row:   min persistence = 1/1358
K33 12->36:            min persistence = 1/574
fibbinary row:         min persistence = 1/574
```

This separates three proof currencies that raw safe mass mixes:

```text
bar existence     proves strict openness
bar height        gives threshold stability
bar length/count  gives topology and certificate anchoring budget
```

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
lonely_profile_persistence_barcode
exact_safe_component_front
fejer_interval_certificate
ramanujan_exact_period_projector
automaton_gap_state
divisor_abundancy_side_channel
raw_maximin_scalar
raw_sequence_name
```

Pairwise observable:

```text
(predicate, topology, exact_height, endpoint_geometry,
 stability, certificate_handoff, anti_scalar)
```

Switch/gauge: majority comparison of the retained-coordinate vector, with the
listed carrier order as the tie Hamiltonian path.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

The transitive order places the persistence barcode first because it keeps the
LRC predicate, topology, exact height, stability margin, and certificate handoff
interface.  Raw `M` is useful, but it forgets whether the witness is a wide
robust component or a thin near-boundary sliver.

## Packet Fields

Add a sidecar group to HYP-2963-style packet records:

```text
lonely_profile_bar_count
lonely_profile_total_length
lonely_profile_longest_bar
lonely_profile_peak_height
lonely_profile_peak_time
lonely_profile_persistence_margin
lonely_profile_component_signature
```

The field group should travel with existing exact `M`, qdiv, endpoint-owner,
Haar/Baire, Fejer, Ramanujan, automaton, and divisor labels.

## Assumption Challenge

Alternate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid
circuits, persistence bars, and proof obligations.

The chosen quotient preserves the LRC predicate "there is an open safe bar
above `1/14`."  It destroys raw runner identity only after exact scale, safe
topology, and peak geometry are retained.  A future quotient may collapse bars
only if it proves the LRC predicate is constant on the collapsed barcode class
or emits a named residual such as K33/THM-572/F7.

## Next Pulls

1. Run the barcode sidecar over the full HYP-2963 packet bank and group rows by
   `(route, bar_count, longest_bar, min_persistence)`.
2. Identify the lowest-persistence positive families; these are likely the
   finite rows that should receive interval Fejer or endpoint-owner certificates
   first.
3. Compare barcode classes against HYP-3012 induced automaton tournament
   classes and HYP-3013 divisor/abundancy fields to find mixed fibers.
4. Convert the strongest bars into rational interval centers for HYP-2981 /
   HYP-2974 Fejer certificate manifests.

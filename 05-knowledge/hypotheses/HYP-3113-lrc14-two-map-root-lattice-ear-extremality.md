---
id: HYP-3113
title: LRC14 two-map root-lattice-ear extremality synthesis
status: SYNTHESIS + exact map scout; not a proof
source: codex-2026-06-27-S265
script: 04-computation/lrc14_two_map_root_lattice_ear_extremality_codex_s265.py
result: 05-knowledge/results/lrc14_two_map_root_lattice_ear_extremality_codex_s265.out
tangent: T1189
technique: LTI-250
tournament_technique: LTT-148
related:
  - HYP-3112
  - HYP-3111
  - HYP-3110
  - HYP-3109
  - HYP-3108
  - HYP-3107
  - HYP-3106
  - HYP-3105
  - HYP-3104
  - HYP-3103
  - HYP-3102
  - HYP-3101
  - HYP-3098
  - HYP-3097
  - HYP-3096
  - HYP-3095
  - HYP-3062
  - HYP-2879
  - HYP-2874
  - THM-577
  - OPEN-Q-108
---

# HYP-3113: LRC14 Two-Map Root-Lattice-Ear Extremality Synthesis

## Claim

The LRC14 extremality frontier should be read through two coupled maps:

```text
root-curve map:
  single value -> full PGF coefficients -> root locus
  -> Lee-Yang zero-free region / discriminant break
  -> packet-keyed extremality certificate

memory-lattice-ear map:
  raw runner set -> Savitch-style midpoint sidecar recursion
  -> relation-lattice Bravais shape and successive minima
  -> strong/odd/nested ear certificate grammar
  -> first-obstruction cocycle / packet sheaf legal exit
```

The recurring lesson is the same in both maps: a single value is not the
object.  `p0`, cap value, `H`, covolume, shortest vector, root count, or local
exchange optimality can all be correct shadows while the proof-relevant
coordinate is hidden.  The proof carrier is the whole curve, whole lattice
shape, or whole certificate decomposition, plus a declared legal-forgetting
sidecar.

## Prompt Translation

Savitch's theorem is used here as a memory-compression model: reachability can
be certified by recursive midpoint witnesses, but only if each midpoint keeps
the sidecar that the next operation needs.  The LRC transfer is not complexity
theory; it is a proof-interface warning.  A bounded recursive witness route
must record which coordinate is destroyed at each split.

Bravais lattices are used as a relation-lattice shape warning.  A lattice is
not determined by a shortest vector or covolume alone; its Gram shape,
automorphism chamber, anisotropy, successive minima, and low-height walls can
carry the extremal obstruction.  This extends HYP-3062: add a Bravais shape
class before using lattice pressure as an LRC proof carrier.

Lee-Yang extremality is used as the root-curve target suggested by HYP-3103:
do not measure only `G_N(0)=p0`; measure the full miss-count PGF and its zeros.
The bold route is a zero-free-region certificate near the LRC evaluation, with
the k=8/k=9 extremal break treated as a discriminant/root-collision event.

The `exp(-lambda S^4 - b S^2) dS` cue supplies the quartic stabilizer signal:
positive quartic confinement says the fourth cumulant can stabilize a bad
quadratic direction.  In LRC language, `quartic_cumulant_stabilizer` should be
measured against boundary mass, gK8/S2 moment pressure, relation-lattice
anisotropy, and far-block curvature tails.

Ear decomposition supplies the certificate grammar.  Strongly connected
digraphs have ear decompositions; factor-critical graphs have odd ear
decompositions; 2-vertex-connected series-parallel graphs have nested ear
decompositions.  The LRC transfer is to proof-obligation graphs: when a root
collision or observer gluing failure appears, ask whether the proof graph has a
strong ear, an odd parity ear, a nested refinement ear, or a crossing ear that
must route to K33/THM-572/F7 debt.

## Evidence From The Scout

Added exact synthesis scout, rebased after the HYP-3108/HYP-3109 exact
Lee-Yang root-curve work, HYP-3112's one-runner ear-payload ledger, and the
HYP-3111 carrier-sidecar atlas:

- `04-computation/lrc14_two_map_root_lattice_ear_extremality_codex_s265.py`
- `05-knowledge/results/lrc14_two_map_root_lattice_ear_extremality_codex_s265.out`

Incoming S262/HYP-3108 sharpens the Bravais side of this lane.  Its bounded
`{0}+7` scan found that high `p0` tracks no-real-root PGF strata, large
nearest-root radius, and high residue entropy, while the largest mod-7
Bravais peak is negatively correlated with `p0`.  Therefore the S265
`Bravais_relation_shape_class` should not be reduced to covolume, shortest
vector, or Bragg peak.  It should record reciprocal flatness, entropy, and
strict-descent trap data as sidecars before the root-lattice-ear map is used
as a proof carrier.

Incoming S264/HYP-3111 sharpens the certificate side.  Its Minkowski q-body
readout gives a real volume threshold only after the body and predicate are
declared, and its proof-circuit ledger makes root/ear exits multi-sidecar
certificates rather than one-step shortcuts.  Thus S265 should treat
`Bravais_shape_tensor`, `Savitch_midpoint_certificate`, and
`strong_ear_spine` as jointly legal only when they route through an
observer-gluing certificate or finite-address packet.

Incoming S262b/HYP-3112 sharpens the ear side.  The one-runner payload
identity `q_full[t]=q_base[t]-A_t+A_{t+1}` says root motion is controlled by a
retained extension payload, not by the scalar value alone.  S265 therefore
uses ear payloads as one of the proof-carrier signals inside the two-map
portfolio rather than as a competing namespace.

Map A, the root-curve extremality tournament, uses vertices as signals rather
than runners.  Its pairwise observable is majority retention over:

```text
predicate, curve, zeros, cumulants, extremality, transfer, computable, proof_exit
```

Readout:

```text
score_hist={0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 7:3, 9:1}
directed_3cycles=1
Hamiltonian_path_count=3
top path:
root_curve_zero_locus
> Lee_Yang_zero_free_region
> PGF_discriminant_break
> tournament_Iomega_root_spectrum
> miss_count_PGF_root_stratum
> phi4_quartic_cumulant_stabilizer
> full_PGF_coefficients
> fugacity_rank_curve
> single_value_p0
> raw_scalar_rank
```

The nontrivial SCC among `Lee_Yang_zero_free_region`,
`PGF_discriminant_break`, and `tournament_Iomega_root_spectrum` is useful: the
zero-free region, root-collision boundary, and tournament hard-core-gas roots
are not linearly ordered.  They form a feedback triad and should be compared
as a portfolio.

Map B, the memory-lattice-ear certificate tournament, uses proof carriers as
vertices.  Its axes are:

```text
predicate, compression, lattice, ears, parity, quotient, computable, proof_exit
```

Readout:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
Hamiltonian_path_count=1
top path:
packet_sheaf_legal_exit
> first_obstruction_cocycle
> Savitch_midpoint_certificate
> Bravais_shape_tensor
> strong_ear_spine
> successive_minima_covolume
> odd_ear_parity_spine
> nested_ear_series_parallel_spine
> relation_lattice_basis
> raw_runner_vertices
```

This contrast matters.  The analytic/root map is cyclic and portfolio-like.
The memory/lattice/ear map is a stricter certificate ladder.  That suggests
the workflow:

1. use the root map to detect which extremal mechanism is active;
2. use the certificate map to decide what sidecar must be kept before a proof
   quotient is legal.

## Two Comprehensive Maps

### Map 1: Root-Curve Extremality

```text
raw_scalar_rank
  -> single_value_p0
  -> full_PGF_coefficients
  -> miss_count_PGF_root_stratum
  -> PGF_discriminant_break
  -> Lee_Yang_zero_free_region
  -> root_curve_zero_locus
```

Side bridges:

```text
phi4_quartic_cumulant_stabilizer
  -> boundary_mass_charge / gK8 moment pressure / relation-lattice anisotropy

tournament_Iomega_root_spectrum
  -> parallel hard-core gas comparison, not a direct LRC value formula

fugacity_rank_curve
  -> how row rankings change as z moves away from the LRC evaluation
```

The root map preserves the whole PGF curve and root geometry.  It destroys
direct interval topology unless endpoint owners, component/barcode data, and
packet-family labels are attached.

### Map 2: Memory-Lattice-Ear Certificates

```text
raw_runner_vertices
  -> relation_lattice_basis
  -> successive_minima_covolume
  -> Bravais_shape_tensor
  -> Savitch_midpoint_certificate
  -> strong_ear_spine / odd_ear_parity_spine / nested_ear_series_parallel_spine
  -> first_obstruction_cocycle
  -> packet_sheaf_legal_exit
```

Side bridges:

```text
strong ear:
  proof-obligation graph is reachable/strong after extensions

odd ear:
  parity-critical or factor-critical style residual, possible AP/GW or F7 debt

nested ear:
  nested-refinement O2 route; crossing/non-nested ear suggests K33/O3 handoff
```

This map preserves legal proof transport.  It destroys the comfort of a single
runner set or scalar rank, deliberately.

## Coupled Readout

The map alignment table is the actual contribution:

```text
single_value_p0 -> root_curve_zero_locus
  matches
raw_runner_vertices -> packet_sheaf_legal_exit
```

Both say the scalar/rootless object is only a shadow.

```text
PGF_discriminant_break
  matches
strong_ear_spine / nested_ear_series_parallel_spine
```

A root collision should name the structural ear where the proof graph changes.

```text
Lee_Yang_zero_free_region
  matches
Savitch_midpoint_certificate
```

Zero exclusion may be certifiable recursively, but every midpoint split must
retain the sidecar that makes the next quotient legal.

```text
phi4_quartic_cumulant_stabilizer
  matches
Bravais_shape_tensor / successive_minima_covolume
```

Quartic stabilization is a fourth-moment face of relation-lattice shape, not a
free scalar tail estimate.

```text
miss_count_PGF_root_stratum
  matches
first_obstruction_cocycle
```

A near-real root or real-root jump should be treated as a first-obstruction
alarm.

## Bold Hypotheses

1. **Lee-Yang confinement route.**  For the bounded cap banks behind THM-577,
   the extremal AP/consecutive rows may be characterized by a zero-free region
   of `G_N(z)` near the LRC evaluation, not just by `G_N(0)`.

2. **Discriminant break principle.**  The k=8/k=9 optimizer break is a
   discriminant event of the miss-count PGF: a pair of complex roots collides
   with the real axis or crosses a certified root-argument boundary.

3. **Quartic stabilizer principle.**  Rows that fool the quadratic/pair-Pascal
   shadow are separated by a fourth cumulant tied to gK8/S2 and relation-lattice
   anisotropy.  The phi4 analogy says to measure the stabilizing `S^4` term,
   not to trust the quadratic well.

4. **Bravais shape beats covolume.**  The relation-lattice side of LRC14 will
   not close from rank/covolume/short vector alone.  The missing feature is a
   Bravais-style shape class: Gram anisotropy, automorphism chamber,
   successive-minima ratios, and low-height wall orientation.

5. **Savitch midpoint proof compression.**  The observer-gluing frontier may
   admit a recursively checked midpoint certificate: split a proof graph into
   two smaller chart-reachability obligations, but emit a first-obstruction
   sidecar whenever the midpoint forgets endpoint owner, root stratum, lattice
   shape, or ear parity.

6. **Ear type predicts terminal route.**  Strong ears are ordinary gluing;
   odd ears are parity/AP-GW/F7 alarms; nested ears are O2 covering discharge;
   non-nested crossing ears are O3/K33/THM-572 handoff signals.

7. **Root-lattice-ear resonance portfolio.**  The most useful next invariant
   is not any one of PGF roots, Bravais shape, or ears.  It is the portfolio
   recording whether a root-locus event, lattice-shape anisotropy, and
   proof-graph ear type co-occur on the same packet.

## New Signals To Measure

```text
PGF_zero_locus_signature
nearest_zero_to_LRC_evaluation
Lee_Yang_confinement_margin
PGF_discriminant_break_index
quartic_cumulant_stabilizer
root_curve_vs_Iomega_arc_distance
Savitch_midpoint_sidecar_depth
Bravais_relation_shape_class
successive_minima_anisotropy
ear_certificate_type
odd_ear_parity_debt
nested_ear_crossing_defect
root_lattice_ear_resonance_portfolio
```

Add these beside HYP-3104's existing maximizer signals, especially
`miss_count_PGF_root_stratum`, `boundary_mass_charge`,
`decorrelation_inward_flux`, `exchange_trap_index`, and
`first_live_interaction_order`.

## Tournament Analysis

Assumption challenge: the vertices need not be runners or arcs.  This session
explicitly considered runners, gaps, fixed circle sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes, PGF roots, relation
lattices, Bravais cells, strong ears, odd ears, nested ears, obstruction
syndromes, and proof obligations.

Chosen vertices:

- Map A: signal families on the PGF/root/fugacity side.
- Map B: proof carriers on the memory/lattice/ear side.

Pairwise observables and gauges are declared in the scout.  The preserved LRC
predicate is not the original row identity; it is the proof obligation plus
the named sidecar needed to make forgetting legal.  The destroyed coordinates
are recorded in each vertex ledger.

Tie Hamiltonian paths:

```text
root_curve_zero_locus
> Lee_Yang_zero_free_region
> PGF_discriminant_break
> miss_count_PGF_root_stratum
> phi4_quartic_cumulant_stabilizer
> tournament_Iomega_root_spectrum
> full_PGF_coefficients
> fugacity_rank_curve
> single_value_p0
> raw_scalar_rank
```

```text
packet_sheaf_legal_exit
> first_obstruction_cocycle
> Savitch_midpoint_certificate
> Bravais_shape_tensor
> successive_minima_covolume
> relation_lattice_basis
> strong_ear_spine
> odd_ear_parity_spine
> nested_ear_series_parallel_spine
> raw_runner_vertices
```

## Next Tests

1. Join HYP-3103's miss-count PGF root data to HYP-3104's exchange-trap and
   inclusion-exclusion anatomy rows; record `PGF_discriminant_break_index`.
2. Compute `nearest_zero_to_LRC_evaluation` and `Lee_Yang_confinement_margin`
   across the bounded cap bank.
3. Add `quartic_cumulant_stabilizer` to the gK8/S2/Clebsch moment rows and
   compare against boundary/middle mass transfer.
4. Add `Bravais_relation_shape_class` and `successive_minima_anisotropy` to
   relation-lattice packets already carrying HYP-3062 fields.
5. Build a proof-obligation graph from observer-gluing rows and classify its
   ear type: strong, odd, nested, or crossing.
6. Test whether root collisions, Bravais anisotropy jumps, and crossing ears
   occur on the same residual rows; if yes, promote
   `root_lattice_ear_resonance_portfolio` into the HYP-2963 packet ledger.

## Dependencies

Builds on HYP-3104's maximizer signal atlas, HYP-3103's miss-count PGF zeros
and perspective groupoid notes, HYP-3102 first-obstruction cocycles,
HYP-3101 component-bound topology, HYP-3098 observer gluing, HYP-3097
Pascal/equivalence shadows, HYP-3096 polynomial-method witness route,
HYP-3062 Roth-Minkowski lattice fence, HYP-2879 strong-ear atom calculus, the
tournament root-spectrum notes, and THM-577's cap value closure.

---
id: HYP-3073
status: SYNTHESIS / finite proof-interface scout; not a proof of LRC14
source: codex-2026-06-26-S239
tags: [lrc14, signed-polymers, dirichlet-energy, schur-complement, cross-disciplinary, tournament-analysis]
script: 04-computation/lrc14_polymer_dirichlet_bridge_codex_s239.py
result: 05-knowledge/results/lrc14_polymer_dirichlet_bridge_codex_s239.out
related:
  - HYP-3072
  - HYP-3071
  - HYP-3070
  - HYP-3069
  - HYP-3066
  - HYP-3037
  - HYP-2645
  - HYP-2632
  - HYP-2540
  - OPEN-Q-108
---

# HYP-3073: Renormalized Polymer / Dirichlet Bridge for LRC14

This pass was renumbered from the local HYP-3072 stub after the rebase revealed
that incoming S238 had already claimed HYP-3072 for the cross-carrier pullback
portfolio.  HYP-3073 should therefore read HYP-3072 as upstream signal: the
portfolio chooses promising proof carriers, while this note asks whether two
older non-median carriers can be made theorem-facing.

The synthesis deliberately avoids another median-center continuation.

1. **Renormalized signed polymers.**  The old absolute Mayer/polymer gas failed
   for the LRC14 relation lattice: the raw absolute activity is not summable in
   the useful form.  HYP-3073 does not revive that argument.  It retypes
   polymers by packet class and signed activity: AP core, near-AP wall,
   odd-sector/AP-like current, repeated-residue character, and wide/Sidon
   finite-cell packets become different species.
2. **Dirichlet/Schur sidecar energy.**  The residual capacitor and cycle-class
   work can be read as an electrical network.  Sidecars are boundary
   conditions, not decorations.  A legal quotient must preserve boundary energy
   to a named discharge, or explicitly route the lost current to AP/GW,
   owner-strip descent, a known cycle span, or named THM-572/F7 debt.

## Scout Evidence

The S239 script cross-reads HYP-2540 signed-polymer/Riesz work, HYP-2645
Poisson-theta finite-cell data, HYP-2632 repeated-residue character channels,
HYP-3037 residual capacitor cuts, HYP-3071 cycle-class observability, and the
incoming HYP-3072 cross-carrier pullback portfolio.

From the HYP-2645 relation-packet table:

```text
packet        R6_box2  corr      corr/R6      density_vs_AP  coherence
AP_core          982   0.302731  0.00030828      1.000000    0.009661
near_AP          924   0.183854  0.00019898      0.940937    0.006048
odd_AP           546   0.213476  0.00039098      0.556008    0.009136
primes_ish       414   0.009649  0.00002331      0.421589    0.000474
squares_sidon    156   0.000536  0.00000344      0.158859    0.000043
```

Raw relation density is therefore not theorem currency.  The odd-sector packet
has fewer R6 relations than the near-AP packet but larger signed correction.
Wide/Sidon packets have tiny signed activity and should be routed through
finite-cell/Poisson-theta control, not through an absolute Mayer envelope.
AP-like packets are the slow sector and need bounded-core, Riesz-product, or
cross-carrier pullback certificates.

The Dirichlet toy network uses S237 first-tooth counts and the HYP-3071 cycle
template rank:

```text
raw_route_effective_conductance = 1/2
legal_sidecar_effective_conductance = 9
conductance_ratio_legal/raw = 18
raw_min_cut = 1
legal_min_cut_to_discharge = 27
legal_min_cut_to_phantom_f7 = 1
```

Interpretation: collapsing to raw route labels leaves a unit-capacity fake cut.
Keeping topology/stalk/cycle sidecars Schur-complements to a positive energy
certificate.  The phantom F7 coordinate remains a one-unit boundary exit,
matching HYP-3071's "outside the known cycle image" reading.

## Proposed Packet Fields

```text
signed_polymer_packet_type
signed_activity_budget
finite_cell_route
renormalized_activity_exit
dirichlet_sidecar_graph_id
dirichlet_boundary_potential
schur_complement_conductance
sidecar_energy_exit
phantom_f7_boundary_atom
```

## Tournament Analysis

Vertex candidates considered: runners, route labels, packet rows,
relation-lattice supports, polymers, residual currents, sidecar boundary
conditions, cycle-class atoms, cross-carrier rows, and proof obligations.

Chosen vertices: proof carriers and renormalization/energy obligations.

Pairwise observable: retained payload, positive test measure, signed
cancellation, boundary energy, finite-cell exactness, discharge specificity,
and cross-carrier reuse.  The switch orients toward lexicographically stronger
payload.

The scout tournament is transitive:

```text
score_hist = {0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
tie_path =
  renormalized_signed_polymer >
  dirichlet_schur_certificate >
  cross_carrier_pullback_portfolio >
  cycle_class_observability >
  riesz_positive_test_measure >
  residual_capacitor_min_cut >
  poisson_finite_cell >
  repeated_residue_character >
  absolute_mayer_shadow >
  raw_route_scalar
```

## Assumption Challenge

This pass does not assume tournament vertices are runners, route labels, or
median centers.  The useful vertices are proof carriers, typed polymers, and
sidecar-energy obligations.  The preserved LRC predicate is that a
sidecar-completed residual cannot cover/discharge illegally.  The destroyed
data are individual runner labels, raw relation counts, and raw route labels.
The quotient is legal only when signed activity is renormalized by packet type
or boundary current is preserved by a Dirichlet/Schur certificate.

## Theorem Targets

1. **Polymer target:** define packet activities by
   `(packet_type, signed_activity, finite_cell_route)`, then prove wide/Sidon
   and repeated-residue activities are summable after AP-like cores are removed
   or certified.
2. **Dirichlet target:** build the actual HYP-2963 residual sidecar graph and
   show every admissible Schur complement preserves positive conductance to
   named exits or isolates phantom F7 as a concrete boundary atom.

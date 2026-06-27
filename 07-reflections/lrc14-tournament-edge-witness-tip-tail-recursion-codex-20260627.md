# LRC14 tournament edge-witness tip-tail recursion

The useful shift is:

```text
edge != orientation
edge = information channel between two proof states
```

A tournament arc `u -> v` is not yet a witness.  It becomes a witness only
when it says what information is lost at the tail, what information is demanded
at the tip, and why those two operations commute for the LRC predicate.

## Core Definition

For a proof quotient `q` and directed edge `e=(tail -> tip)`, define:

```text
EdgeWitness_q(e) =
  (Tail_q(e), Tip_q(e), Orbit_q(e), Comm_q(e), Exit_q(e)).
```

`Tail_q(e)` is the information needed to delete or condition on the tail.
`Tip_q(e)` is the information needed to extend or observe at the tip.
`Orbit_q(e)` is the observer-cut payload orbit modulo the visible symmetries.
`Comm_q(e)` measures the failure of tail deletion and tip extension to commute.
`Exit_q(e)` is positive `Phi/P`, finite address, observer gluing, coordinate
resurrection, Node-3 decorrelation, or named residual debt.

This is just HYP-3054/HYP-3056 with edge language sharpened.  HYP-3049 gives
the finite toy model: the missing object at the first A000568 failure is not
more node memory, but the old-root/new-observer directed edge payload.

The HYP-3125/HYP-3126/HYP-3127/HYP-3128/HYP-3129/HYP-3130/HYP-3131/HYP-3132/HYP-3133/HYP-3134/HYP-3135/HYP-3136/HYP-3137
floor-decoupling lanes add a concrete analytic instance: the tip can be a
single-far Lee-Yang factor, the tail can be the bounded core left after
deletion or averaging, and the witness is the commutator between
contract-tip-first and decouple-wide-tip-first.  HYP-3138 adds the even-fold
k=8 warning: reflection-fold compression is useful only when an odd-coordinate
resurrection sidecar says which coordinate can be recovered downstream.
HYP-3139 converts that warning into block-page data: inner shell, center
coupling, boundary leakage, and phi4 sign have to be retained.  HYP-3140 adds
the fiber-PGF version of the same discipline for `Rprime`: keep the sheet PGF,
Q-masked PGF, conditional first moment, and SPEC/Fourier sidecar.  The
raw zero-free assertion is
not enough; HYP-3128 shows the overcrowded tail obstruction, HYP-3129 supplies
the SPEC certificate for the surviving floor, HYP-3130 adds the
Gaussian/Beurling-Selberg minorant tail certificate, and HYP-3131 tests the
far-zero-push monotonicity route.  HYP-3132 says the bounded-core residue is a
solvable De Moivre/phi4 biquadratic, while HYP-3133 adds the A000568
extension-shadow layer between sector word and paired child deck.  HYP-3134
then identifies the A000568 layer as a global-consistency quotient between the
lower sector envelope and upper paired-child envelope, while HYP-3135 says the
De Moivre/resolvent clue is a middle elementary-symmetric pair/triple payload.
The edge has to name the contraction order word, factor id, zero-free region,
rate budget, SPEC certificate, minorant tail certificate, far-zero-push status,
bounded-core resolvent status, A000568 envelope/gluing status, resolvent
middle payload, integrated-floor child/gluing status, GF coefficient layer,
root locus, log-derivative/cumulant status, obstruction status, and
bounded-core floor exit.

## Recursive Tip And Tail

Tail recursion asks:

```text
If I delete this source endpoint, which coordinate stops being reconstructible?
```

Tip recursion asks:

```text
If I extend at this target endpoint, which observer payload is newly demanded?
```

The edge witness is the square between them.  If the square commutes, the edge
is a legal proof transport.  If it does not, the commutator is the next missing
sidecar.

That gives a clean place for the current threads:

- HYP-3049: old-root/new-observer edge role and cross-sector orientation.
- HYP-3053: rectangle/hourglass line-cycle residues.
- HYP-3054/HYP-3056: observer-cut orbit and discharge mode.
- HYP-3118: coordinate-resurrection cover for the commutator.
- HYP-3116/HYP-3117: missing-input vector decreases along a legal edge.
- HYP-3111/HYP-3115: Minkowski relation pressure, circuit complexity, Ising
  walls, and de Moivre quintic residuals become edge payloads, not standalone
  scalar criteria.
- Incoming phi4 HYP-3122: the cap dip and `kappa4` sign are quartic edge
  payloads, not rowwise scalar corrections.
- HYP-3120: finite-address/Phi receiver for edge transport.
- HYP-3121: `R-safe -> Q-lonely` is an event-edge whose witness is the
  decorrelation floor and low-frequency resonance payload.

## Creative Hypotheses

The edge witness may be the right atomic unit for the proof.  Vertices are
static; edges are where a proof quotient is asked to survive an operation.

Possible definitions:

```text
edge_information_gain(e) =
  rank of predicate-relevant coordinates retained by EdgeWitness(e)
  minus rank of raw orientation u->v.

edge_commutator_debt(e) =
  minimal repair-cover rank needed to make tail deletion and tip extension
  commute.

edge_entropy_drop(e) =
  number of coarse-fiber packet pairs separated by the edge payload.

edge_terminality(e) =
  first exit in {Phi/P, finite_address, observer_gluing,
                 coordinate_resurrection, decorrelation_floor,
                 named_residual_debt}.
```

For the Minkowski/circuit/Ising/de Moivre lane, the stronger packet is:

```text
MCIQEdge(e) =
  (minkowski_relation_wall_class,
   successive_minima_proxy,
   proof_circuit_size_depth_fanin,
   circuit_uniformity_guard,
   ising_domain_wall_id,
   domain_wall_legal_exit,
   phi4_quartic_cumulant_delta,
   phi4_lambda_sign,
   de_moivre_quintic_residual_delta,
   de_moivre_branch_orbit_word).
```

This reads HYP-3115's `10084` one-swap Ising domain-wall edges as witnesses:
tail side = the row before the swap, its root stratum, relation wall, and
proof-circuit inputs; tip side = the row after the swap, its collision ear,
branch orbit, and terminal exit.  A wall edge that cannot say whether it is a
legal root-collision ear, observer-gluing exit, finite-address exit, or
forbidden wall is not proof evidence yet.

Wild guess: the recurring structures all express edge commutators:

```text
cross-sector orientation
endpoint-owner transfer
Lee-Yang ear payload
Minkowski low-height wall owner
Ising domain-wall legality
phi4 quartic cumulant sign
de Moivre fifth-root branch orbit
De Moivre biquadratic resolvent status
hidden fold/virtual sum
magnitude cocycle
first obstruction cochain
rectangle/hourglass residue
low-frequency decorrelation correction
minorant-tail certificate
far-zero-push monotonicity
A000568 extension shadow
A000568 envelope/global-consistency gluing status
resolvent pair/triple middle payload
generating-function coefficient layer
PGF root-locus status
log-derivative/cumulant payload
```

They may be different projections of the same "edge remembers what a vertex
quotient forgot" principle.

## What To Measure

The next packet rows should add:

```text
edge_witness_packet_id
edge_tail_payload_word
edge_tip_payload_word
tail_delete_recursion_depth
tip_extend_recursion_depth
tip_tail_commutator_defect
edge_cut_payload_orbit_id
old_new_endpoint_role
cross_sector_orientation_word
edge_gf_carrier_type
edge_coefficient_payload_layer
edge_pgf_root_locus_status
edge_log_derivative_cumulant_status
edge_k8_reflection_fold_adjoint_status
edge_odd_coordinate_resurrection_status
edge_a000568_extension_shadow
edge_a000568_envelope_position
edge_global_consistency_class
edge_child_gluing_status
edge_information_gain_rank
edge_predicate_delta
edge_coordinate_resurrection_cover
edge_missing_input_delta
edge_proof_circuit_size_depth_fanin
edge_circuit_uniformity_guard
edge_circuit_uniformity_guard_status
edge_phi_p_activation_delta
edge_minkowski_relation_wall_class
edge_minkowski_covolume_threshold_status
edge_successive_minima_proxy
edge_ising_domain_wall_id
edge_ising_partition_zero_locus_status
edge_domain_wall_legal_exit
edge_ear_payload_vector
edge_phi4_quartic_cumulant_delta
edge_phi4_lambda_sign
edge_de_moivre_quintic_residual_delta
edge_de_moivre_auxiliary_quadratic_status
edge_de_moivre_biquadratic_resolvent_status
edge_resolvent_middle_payload_status
edge_resolvent_pair_triple_layer
edge_signed_spec_resolvent_packet_status
edge_de_moivre_branch_orbit_word
edge_rectangle_hourglass_residue
edge_decorrelation_floor_status
edge_asano_contraction_order_word
edge_single_far_factor_id
edge_zero_free_region_status
edge_wide_decoupling_rate_bound
edge_spec_certificate_status
edge_lee_yang_obstruction_status
edge_minorant_tail_certificate_status
edge_far_zero_push_status
edge_bounded_core_floor_exit
edge_terminal_exit_or_debt
```

The most important first target is HYP-3115's Ising domain-wall graph.  It has
`10084` root-stratum wall edges.  HYP-3141 says each wall edge should be
classified by tail payload, tip payload, and commutator defect before we use
it as evidence for extremality.

## Tournament Readout

The scout ranks edge-witness carriers, not runners or raw arcs:

```text
recursive_tip_tail_edge_witness
> asano_tip_contraction_edge
> observer_cut_payload_orbit
> finite_address_phi_edge
> ising_domain_wall_edge
> coordinate_resurrection_edge
> directed_edge_perspective
> endpoint_owner_transfer_edge
> phi4_quartic_stabilizer_edge
> lee_yang_ear_motion_edge
> de_moivre_stationary_quintic_edge
> proof_circuit_missing_input_edge
> minkowski_relation_wall_edge
> node3_decorrelation_edge
> source_sink_apex_arc
> rectangle_hourglass_flow_edge
> raw_tournament_arc
```

The raw arc is last.  That is the main guardrail.  Tournament edges are
valuable as witnesses only when they transport proof information.

Related: HYP-3141, HYP-3140, HYP-3139, HYP-3138, HYP-3137, HYP-3136, HYP-3135, HYP-3134, HYP-3133, HYP-3132, HYP-3131, HYP-3130, HYP-3129, HYP-3128, HYP-3127, HYP-3126, HYP-3125, HYP-3124, HYP-3122, HYP-3121, HYP-3120, HYP-3119,
HYP-3118, HYP-3117, HYP-3116, HYP-3115, HYP-3112, HYP-3056, HYP-3054,
HYP-3053, HYP-3049, HYP-3111, HYP-3110, HYP-3103, HYP-3062, HYP-2008,
OPEN-Q-108.

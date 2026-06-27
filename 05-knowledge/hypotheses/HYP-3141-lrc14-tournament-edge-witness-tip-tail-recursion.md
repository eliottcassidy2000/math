---
id: HYP-3141
title: LRC14 tournament edge-witness tip-tail recursion
status: SYNTHESIS / executable edge-witness scout and proof-packet proposal; not a proof
source: codex-2026-06-27
tangent: T1206
technique: LTI-267
tournament_technique: LTT-165
script: 04-computation/lrc14_tournament_edge_witness_tip_tail_codex_20260627.py
result: 05-knowledge/results/lrc14_tournament_edge_witness_tip_tail_codex_20260627.out
reflection: 07-reflections/lrc14-tournament-edge-witness-tip-tail-recursion-codex-20260627.md
related:
  - HYP-3140
  - HYP-3139
  - HYP-3138
  - HYP-3137
  - HYP-3136
  - HYP-3135
  - HYP-3134
  - HYP-3133
  - HYP-3132
  - HYP-3131
  - HYP-3130
  - HYP-3129
  - HYP-3128
  - HYP-3127
  - HYP-3126
  - HYP-3125
  - HYP-3124
  - HYP-3123
  - HYP-3122
  - HYP-3121
  - HYP-3120
  - HYP-3119
  - HYP-3118
  - HYP-3117
  - HYP-3116
  - HYP-3115
  - HYP-3113
  - HYP-3111
  - HYP-3110
  - HYP-3112
  - HYP-3109
  - HYP-3103
  - HYP-3062
  - HYP-3056
  - HYP-3054
  - HYP-3053
  - HYP-3049
  - HYP-3045
  - HYP-2008
  - HYP-2112
  - HYP-2108
  - THM-571
  - HYP-2968
  - OPEN-Q-108
---

# HYP-3141: LRC14 Tournament Edge-Witness Tip-Tail Recursion

## Claim

Tournament edges should be treated as information-bearing proof witnesses, not
only as binary orientations.  A directed edge

```text
e = (tail -> tip)
```

is theorem-facing only after it carries:

```text
edge_witness(e) =
  (tail_delete_payload,
   tip_extend_payload,
   observer_cut_payload_orbit,
   tip_tail_commutator_defect,
   coordinate_resurrection_cover,
   predicate_delta,
   terminal_exit_or_named_debt)
```

The recursive rule is:

```text
delete tail, then extend tip
  must agree with
extend tip, then delete tail
```

up to a named discharge: zero, reconstruction, dual annihilation, familywise
descent, AP/GW boundary stop, finite-address/observer-gluing exit, positive
`Phi/P`, Node-3 decorrelation floor, or named residual debt.

Read abstractly, the edge is an information channel.  The raw orientation
`u -> v` is telemetry.  The witness is the transported sidecar that lets the
LRC predicate survive recursive deletion and extension.

## Why This Extends The Current Frontier

HYP-3049 already showed that the first A000568/rooted-perspective failure is
fixed by an ordered-pair edge payload: old-root/new-observer role, incident
word, sector deck, and cross-sector orientation.  HYP-3054/HYP-3056 abstracted
that as observer-extension/cut payload orbits.  HYP-3118 turned destroyed
coordinates into repair covers and adjoint sections.  HYP-3120 made those
fields part of a finite packet receiver.  HYP-3121 then added the covering
case as a decorrelation edge between `R-safe` and `Q-lonely` lift events.

HYP-3141 extends HYP-3124's exact edge-recursion census: the proof edge is the
smallest unit that can
transport missing information from one proof state to the next.

The HYP-3111/HYP-3115 Minkowski/circuit/Ising/de Moivre lane and the incoming
phi4 quartic-stabilizer lane become edge-level in this language.  A one-swap
Ising wall is not evidence until its edge witness says which relation wall,
proof-circuit input, root-collision ear, phi4 cumulant sign, de Moivre branch
orbit, and terminal exit travel from tail to tip.  Raw Minkowski volume, a
fitted circuit threshold, an Ising energy, a quartic dip, or a quintic
residual is only telemetry unless it is attached to that edge packet.

## Recursive Definitions

For a proof quotient `q:X -> Y` and an edge `e=(a->b)`, define:

```text
Tail_q(e) =
  payload needed to delete, condition on, or contract the tail a;

Tip_q(e) =
  payload needed to add, observe, or extend at the tip b;

Orbit_q(e) =
  observer-cut payload orbit modulo the visible automorphism group;

Comm_q(e) =
  Extend_tip(Delete_tail(x)) - Delete_tail(Extend_tip(x)).
```

The edge is safe for the LRC predicate only if `Comm_q(e)` is discharged.  If
not, `Comm_q(e)` is not an error term; it is the first missing proof
coordinate.

Tail recursion asks what can be forgotten when the source endpoint is removed.
Tip recursion asks what must be added when the target endpoint is observed.
The edge witness is the compatibility square between those two recursions.

## Two Maps To Keep Separate

HYP-3124 supplies the exact combinatorial map:

```text
raw edge
  -> outside four-sector deck
  -> paired tail/tip endpoint-deletion children
  -> recursive child edge decks.
```

HYP-3141 supplies the proof-information overlay:

```text
recursive edge deck
  -> carried proof sidecar
  -> tip-tail commutator
  -> terminal discharge or named debt.
```

The gold is in keeping both maps.  The first map says which local edge
quotients are distinguishable.  The second map says which distinctions matter
to the LRC predicate.

Define the depth-`k` bideck:

```text
BiDeck_k(e) =
  (TailDeck_k(e), TipDeck_k(e), CrossDeck_k(e), ExitDeck_k(e)).
```

`TailDeck_k` recursively deletes or contracts from the source side.
`TipDeck_k` recursively extends or observes from the target side.  `CrossDeck`
records the observer-cut orbit, sector orientation, and sidecar transport
between them.  `ExitDeck` records the first legal stop.  The speculative
theorem target is:

```text
LRC predicate is constant on an edge quotient
  iff
BiDeck_k(e) has zero unresolved commutator
for a bounded k or names a standard residual sector.
```

Wild definition: treat a tournament edge as a proof boundary operator,

```text
partial_W(e) = TipPayload(e) - TailPayload(e).
```

Then `Comm(e)` is the obstruction to `partial_W^2=0` after a quotient.  In
this language, Lee-Yang ears, Ising walls, Minkowski relation walls, phi4
quartic stabilizers, de Moivre branch jumps, and circuit missing-input gates
are not analogies.  They are candidate bases for the same edge-boundary
information complex.

## Asano, SPEC, And Wide-Decoupling Integration

Incoming HYP-3125/HYP-3126/HYP-3127/HYP-3128/HYP-3129/HYP-3130/HYP-3131/HYP-3132/HYP-3133/HYP-3134/HYP-3135/HYP-3136/HYP-3137/HYP-3138/HYP-3139/HYP-3140 add a stronger analytic
interpretation of the same tip/tail rule.  HYP-3125 supplies the edge-floor
packet ledger, HYP-3127 promotes the multi-far floor to an Asano contraction
of single-far factors, HYP-3128 proves the Q-block Lee-Yang factor while
exposing the overcrowded tail obstruction, and HYP-3129 supplies the elementary
SPEC certificate for the surviving quasi-independence floor.  HYP-3130 adds
the Gaussian/Beurling-Selberg minorant and uniform-tail certificate, while
HYP-3131 says far elements push the PGF zeros outward and reduce the multi-far
work to a bounded-core Lee-Yang problem.  HYP-3132 identifies the k=8 bounded
core as a solvable De Moivre/phi4 biquadratic resolvent, and HYP-3133 inserts
the A000568 one-extra-vertex extension shadow between sector words and paired
child decks.  HYP-3134 sharpens that shadow into an edge-envelope quotient:
raw four-sector data is the lower envelope, paired tail/tip child data is the
safe upper envelope, and A000568 is the named global-consistency class between
them.  HYP-3135 adds the algebraic analogue: the De Moivre/resolvent lesson is
to keep the middle elementary-symmetric pair/triple payload before collapsing
to a scalar residual.  HYP-3136 packages the assembled factorization
`L(S)=Rprime*meas(R-safe)*meas(Q-lonely)` and says which child carries the
Q/apex floor, bounded-core tail, signed-SPEC coupling, and legal quotient
payloads.  HYP-3137 supplies the generating-function payload atlas, so an edge
should record whether the useful GF data is coefficient layer, root locus,
log-derivative/cumulant, quotient legality, signed tail certificate, or
terminal exit; the completed atlas says the smallest full-cover packet is
signed SPEC resonance series plus A000568 cycle-index quotient plus
miss-count PGF root locus.  HYP-3138 sharpens the k=8 bounded-core tail into
a reflection-fold coordinate-resurrection test: the even fold is useful only
with an odd-coordinate resurrection table or named finite-address /
observer-gluing debt.  HYP-3139 then makes the reflection block concrete at
matrix level: the edge witness should remember the inner `2x2` shell page,
center coupling, boundary-sector leakage, and HYP-3122 phi4 sign before
calling the bounded-core dip discharged.  HYP-3140 supplies the fiber-PGF
version of the remaining `Rprime` floor: the edge witness should remember the
14-sheet PGF, Q-masked PGF, conditional first moment, and SPEC/Fourier
sidecar before scalarizing to an `Rprime` value.  HYP-3126 remains the
wide-rate limit: a sufficiently wide speed decorrelates from the bounded core
at an elementary `O(1/w)` rate and contributes its mean factor `6/7`.

For edge-witness language, this suggests a new carrier:

```text
asano_tip_contraction_edge(e) =
  (single_far_factor,
   contraction_order_word,
   zero_free_region_certificate,
   spec_or_elementary_decoupling_rate,
   lee_yang_obstruction_status,
   minorant_tail_certificate,
   far_zero_push_status,
   a000568_global_consistency_class,
   resolvent_middle_payload_status,
   bounded_core_floor_exit).
```

Tail recursion is the bounded core left after deleting or averaging a wide
tip.  Tip recursion is the insertion of one far factor before Asano
contraction.  The commutator is the defect between "contract tips first" and
"decouple a wide tip first".  A legal edge either proves these two orders
agree within the HYP-3126/HYP-3129/HYP-3130 rate budget, uses HYP-3131's
far-zero-push monotonicity, or routes the defect to the HYP-3128 bounded-core /
Lee-Yang obstruction.  It may quotient away paired child data only after the
HYP-3134 global-consistency/gluing status is named, and it may compress the
bounded-core algebra only after the HYP-3135 middle resolvent packet is
recorded.  This gives a concrete
proof-facing use for the
abstract boundary operator above:

```text
partial_W(far_tip_edge)
  = contracted_multi_far_floor
    - (6/7) * bounded_core_floor.
```

New signal to add beside the existing packet fields:

```text
edge_asano_contraction_order_word
edge_single_far_factor_id
edge_zero_free_region_status
edge_gf_carrier_type
edge_coefficient_payload_layer
edge_pgf_root_locus_status
edge_log_derivative_cumulant_status
edge_wide_decoupling_rate_bound
edge_spec_certificate_status
edge_lee_yang_obstruction_status
edge_minorant_tail_certificate_status
edge_far_zero_push_status
edge_a000568_envelope_position
edge_global_consistency_class
edge_child_gluing_status
edge_resolvent_middle_payload_status
edge_resolvent_pair_triple_layer
edge_signed_spec_resolvent_packet_status
edge_k8_reflection_fold_adjoint_status
edge_odd_coordinate_resurrection_status
edge_reflection_core_block_status
edge_inner_shell_bound_status
edge_center_boundary_leakage_status
edge_fiber_pgf_word
edge_q_masked_fiber_pgf_word
edge_conditional_first_moment_floor_status
edge_spec_resonance_lattice_status
edge_bounded_core_floor_exit
```

## Scout Results

Executable scout:

```text
04-computation/lrc14_tournament_edge_witness_tip_tail_codex_20260627.py
05-knowledge/results/lrc14_tournament_edge_witness_tip_tail_codex_20260627.out
```

The scout ranks seventeen edge-witness carrier types.  Top readout:

```text
recursive_tip_tail_edge_witness score=47
asano_tip_contraction_edge      score=46
observer_cut_payload_orbit      score=34
finite_address_phi_edge         score=33
ising_domain_wall_edge          score=33
coordinate_resurrection_edge    score=32
directed_edge_perspective       score=31
endpoint_owner_transfer_edge    score=30
phi4_quartic_stabilizer_edge    score=30
de_moivre_stationary_quintic_edge score=29
lee_yang_ear_motion_edge        score=29
proof_circuit_missing_input_edge score=27
minkowski_relation_wall_edge    score=24
node3_decorrelation_edge        score=24
source_sink_apex_arc            score=24
rectangle_hourglass_flow_edge   score=24
raw_tournament_arc              score=-8
```

Tournament fingerprint:

```text
vertices = edge-witness carrier types, not runners or raw arcs
score_hist = {-8:1,24:4,27:1,29:2,30:2,31:1,32:1,33:2,34:1,46:1,47:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
priority_path =
  recursive_tip_tail_edge_witness
  -> asano_tip_contraction_edge
  -> observer_cut_payload_orbit
  -> finite_address_phi_edge
  -> ising_domain_wall_edge
  -> coordinate_resurrection_edge
  -> directed_edge_perspective
  -> endpoint_owner_transfer_edge
  -> phi4_quartic_stabilizer_edge
  -> lee_yang_ear_motion_edge
  -> de_moivre_stationary_quintic_edge
  -> proof_circuit_missing_input_edge
  -> minkowski_relation_wall_edge
  -> node3_decorrelation_edge
  -> source_sink_apex_arc
  -> rectangle_hourglass_flow_edge
  -> raw_tournament_arc
```

The negative control is important: a raw tournament arc ranks last.  The edge
becomes a witness only when it carries tip/tail payload and a commutator
discharge.

## New Signals

Add these fields to HYP-2963/HYP-3098/HYP-3107/HYP-3112/HYP-3115 packet rows:

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
edge_k8_reflection_fold_adjoint_status
edge_odd_coordinate_resurrection_status
edge_reflection_core_block_status
edge_inner_shell_bound_status
edge_center_boundary_leakage_status
edge_fiber_pgf_word
edge_q_masked_fiber_pgf_word
edge_conditional_first_moment_floor_status
edge_spec_resonance_lattice_status
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

## Hypotheses To Test

1. **Edge-witness theorem.**  Every proof-critical LRC14 edge either has a
   discharged tip-tail commutator or emits a named missing coordinate.
2. **Endpoint-role coloop.**  Old/new endpoint role behaves like a coloop in
   the current repair matroid: quotienting it away is legal only with an
   endpoint-owner, finite-address, or observer-cut sidecar.
3. **Lee-Yang ears are edge witnesses.**  The one-runner ear payload `A_t` is
   exactly the tip-extension component of the edge witness for root motion.
4. **Node-3 as event-edge witness.**  HYP-3121's decorrelation floor is an
   edge witness between two event channels, not a Bonferroni scalar.  Its
   tail is `R-safe`, its tip is `Q-lonely`, and its payload is the retained
   low-frequency resonance correction.
5. **Source-sink apex as boundary edge.**  HYP-2008's source-sink arc is the
   minimal edge witness for the transitive/wrap split: it selects the easy
   semicircle branch or the AP-hard back-arc branch.
6. **MCIQ wall packet.**  HYP-3115's `10084` one-swap Ising domain-wall edges
   should be converted into packets
   `(relation_wall, circuit_missing_input, root_collision_ear,
   phi4_cumulant_sign, de_moivre_branch_orbit, legal_exit)`.  A wall that cannot fill those fields
   is phase telemetry, not proof evidence.
7. **Minkowski as an edge sidecar.**  Relation-rich edges are promising only
   when the convex body, covolume/successive-minima proxy, and low-height wall
   deletion status are retained.  Raw short-relation pressure is legal only
   after the edge names why the wall was deleted, reconstructed, or left as
   residual debt.
8. **phi4 as quartic edge stabilizer.**  The cap dip should be measured as an
   edge payload `edge_phi4_quartic_cumulant_delta` with `edge_phi4_lambda_sign`,
   not as a row scalar.  A legal edge either turns the negative `kappa4` /
   positive `lambda` sign into a uniform dip discharge or names the
   phi4/Lee-Yang residual debt.
9. **de Moivre as branch-coordinate alarm.**  The translated `G'(z)` quintic
   residual should be measured across root-collision edges, not optimized
   rowwise.  Its useful output is a solvable local branch model or a named
   fifth-root branch debt.
10. **k=8 reflection fold is an adjoint, not erasure.**  The bounded-core
    biquadratic fold may be used only when
    `edge_k8_reflection_fold_adjoint_status` and
    `edge_odd_coordinate_resurrection_status` certify the odd leakage table or
    name finite-address / observer-gluing debt.
11. **Generating-function payload is edge-local.**  Full PGF/root-locus
    information is retained only as coefficient layer, root locus,
    log-derivative/cumulant, quotient legality, signed tail certificate, or
    terminal-exit fields on the edge witness; a single PGF scalar or root
    count is telemetry.
12. **A000568 envelope as legal forgetting test.**  A sector quotient may
    forget paired tail/tip children only after `edge_a000568_envelope_position`,
    `edge_global_consistency_class`, and `edge_child_gluing_status` are named.
    Otherwise A000568 is a middle count, not a proof quotient.
13. **Resolvent middle payload as branch tax.**  A De Moivre/phi4 edge may
    compress to a scalar residual only after `edge_resolvent_pair_triple_layer`
    and `edge_signed_spec_resolvent_packet_status` record the middle
    elementary-symmetric payload.
14. **Information rank order.**  The useful edge ranking should follow
   information retained for the next proof operation, not tournament novelty
   or scalar edge count.

## Assumption Challenge

Candidate vertices considered:

```text
runners
raw tournament vertices
raw arcs
directed edges
ordered pairs
tip/tail endpoint roles
observer-cut payload orbits
sector-deck matrices
source-sink apex arcs
Lee-Yang ear edges
Asano contraction tip edges
phi4 quartic-stabilizer edges
Minkowski relation-wall edges
Ising domain-wall edges
stationary de Moivre quintic branches
proof-circuit gates
finite-address proof edges
coordinate-resurrection sections
Node-3 decorrelation event pairs
proof obligations
A000568 extension shadows
bounded-core biquadratic resolvents
k8 reflection-fold adjoints
generating-function payload carriers
```

The selected vertices are edge-witness carrier types.  This quotient preserves
progress toward LRC14 through exact `Phi/P`, finite address, observer gluing,
coordinate resurrection, quasi-independence, or named residual debt.  It
destroys raw labels, old/new endpoint role, incident words, cross-sector
orientation, root-motion payload, line-cycle defects, relation-wall owners,
Ising wall legality, phi4 cumulant signs, de Moivre branch data, circuit
uniformity, Asano contraction order, single-far factor identity, zero-free
certificates, wide-decoupling error budgets, SPEC certificates, Lee-Yang
obstruction status, minorant tail certificates, far-zero-push status, and
low-frequency resonance corrections, A000568 extension shadows, full
generating-function coefficient/root-locus/log-derivative payloads, k=8
reflection-fold odd-coordinate resurrection data, and bounded-core
biquadratic resolvent status when only the arc sign remains.

## Next Work

1. Attach edge-witness fields to HYP-3115 one-swap Ising domain-wall edges.
2. Add `edge_cut_payload_orbit_id`, `old_new_endpoint_role`, and
   `cross_sector_orientation_word` to HYP-3098 observer-gluing rows.
3. For HYP-3112 ear payload rows, record the edge witness
   `(base row -> added runner)` with tail/tip recursion depths.
4. For HYP-3121 covering/decorrelation rows, represent `R-safe -> Q-lonely`
   as an event-edge witness and store the low-frequency resonance payload.
5. For HYP-3125/HYP-3127 floor/Asano rows, record
   `edge_asano_contraction_order_word`,
   `edge_single_far_factor_id`,
   `edge_zero_free_region_status`,
   `edge_wide_decoupling_rate_bound`,
   `edge_spec_certificate_status`,
   `edge_lee_yang_obstruction_status`,
   `edge_minorant_tail_certificate_status`,
   `edge_far_zero_push_status`, and
   `edge_bounded_core_floor_exit`, using HYP-3126/HYP-3129/HYP-3130 for
   rate/SPEC/minorant bounds, HYP-3131 for far-zero-push monotonicity, and
   HYP-3128 for the forbidden tail-contraction obstruction.  Then compare the
   contract-tip-first and decouple-wide-tip-first paths as the edge
   commutator.
6. Add the HYP-3111/HYP-3115 fields
   `edge_minkowski_relation_wall_class`,
   `edge_minkowski_covolume_threshold_status`,
   `edge_successive_minima_proxy`,
   `edge_proof_circuit_size_depth_fanin`,
   `edge_circuit_uniformity_guard`,
   `edge_circuit_uniformity_guard_status`,
   `edge_phi4_quartic_cumulant_delta`,
   `edge_phi4_lambda_sign`,
   `edge_ising_partition_zero_locus_status`,
   `edge_de_moivre_quintic_residual_delta`,
   `edge_de_moivre_auxiliary_quadratic_status`,
   `edge_de_moivre_biquadratic_resolvent_status`, and
   `edge_de_moivre_branch_orbit_word` before using a wall edge as proof
   evidence.
7. Add `edge_a000568_extension_shadow` between sector words and paired
   child decks when importing HYP-3133 into HYP-3125/HYP-3129 edge-floor rows.
8. Join HYP-3134's envelope quotient to those rows with
   `edge_a000568_envelope_position`, `edge_global_consistency_class`, and
   `edge_child_gluing_status`, so the quotient records what paired local data
   was legally forgotten.
9. Join HYP-3135's resolvent packet with
   `edge_resolvent_middle_payload_status`, `edge_resolvent_pair_triple_layer`,
   and `edge_signed_spec_resolvent_packet_status` before using De Moivre or
   phi4 scalar residues as proof evidence.
10. Join HYP-3138's k=8 reflection fold with
   `edge_k8_reflection_fold_adjoint_status` and
   `edge_odd_coordinate_resurrection_status` before using the biquadratic
   fold as a legal bounded-core quotient.
11. Join HYP-3139's reflection-block scout with
   `edge_reflection_core_block_status`, `edge_inner_shell_bound_status`, and
   `edge_center_boundary_leakage_status` before using a k=8 block page as
   proof evidence.
12. Join HYP-3140's fiber-PGF certificate with `edge_fiber_pgf_word`,
   `edge_q_masked_fiber_pgf_word`,
   `edge_conditional_first_moment_floor_status`, and
   `edge_spec_resonance_lattice_status` before scalarizing the Rprime floor.
13. Join HYP-3137's generating-function atlas with
   `edge_gf_carrier_type`, `edge_coefficient_payload_layer`,
   `edge_pgf_root_locus_status`, and
   `edge_log_derivative_cumulant_status` before using a PGF scalar, root
   count, determinant ratio, or ordinary sequence value as proof evidence.
14. Reject any future tournament-edge proof that cannot fill
   `tip_tail_commutator_defect` and `edge_terminal_exit_or_debt`.

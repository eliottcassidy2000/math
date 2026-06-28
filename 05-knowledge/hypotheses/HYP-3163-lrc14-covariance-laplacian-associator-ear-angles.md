# HYP-3163: LRC14 covariance Laplacian and associator-ear proof angles

Status: SYNTHESIS / proof-angle proposal; not a proof.

Source: codex-2026-06-28.  This continues HYP-3160's k=8
variance/total-covariance target, HYP-3153's executed
Lee-Yang/Worpitzky/quartic compression packet, HYP-3154's Joukowski/De Moivre
bridge, HYP-3161's exhaustive covariance-max and ferromagnetic-transition
verification, HYP-3162's cyclotomic-ideal / first-cubic-apex calibration, and
HYP-3199's n=4 cover-versus-section warning.  Post-rebase HYP-3200 adds the
exact anchored-bank cumulant audit, and HYP-3201 adds the conditional-entropy
law-defect audit for quotient legality.  KPS-S31al adds a separate
HYP-3201-labeled Toeplitz/Szego and random-current packet; this file
disambiguates that lane as `kps-S31al` because the HYP-3201 label is also used
by the S279 law-defect detail file.

## Claim

After HYP-3161, the live question is no longer whether consecutive maximizes
total empty-sector covariance in the bounded k=8 bank.  It does, exhaustively:
`0` beaters among `3432` bounded clusters.  The remaining proof target is to
explain that finite fact by a transferable certificate.  Split the target into
two deliberately different proof currencies:

1. A degree-two covariance route: prove that the consecutive/AP packet
   maximizes total empty-sector covariance by turning the covariance matrix
   into a Laplacian, Monge, or conditionally positive type certificate.
2. A degree-three odd route: treat the Worpitzky `-9S3` residue as an
   associator or 3-cocycle obstruction, then control it by odd-ear recursion
   or by charging it to even covariance slack.

This is not meant to replace the Lee-Yang root curve.  It says the root curve
and the Joukowski/De Moivre real-rootedness-defect bridge should emit two
sidecars: an even ferromagnetic susceptibility certificate and an odd
associator/ear certificate.

Post-rebase HYP-3153 supplies the finite packet this proposal should consume:
`q0=q6*R^6`, k=8/9/10 `L_y <= cap` margins, the exact k=8 split
`10S0 - 10S1 + 10S2 - 9S3 + 6S4`, Worpitzky `-1/3` edge mode, the
biquadratic `u^4-5u^2+4`, Newton defects, and an odd-ear witness grammar.  In
this HYP-3163 packet, those are inputs to the proof-mechanism search.

Final rebase adds HYP-3162 as a calibration sidecar rather than a competing
route: the 7th-cyclotomic ideal gives the de Moivre angles in the cubic real
subfield of `Q(zeta_7)`, while `n=14` is the first `2q` apex where
`phi(q)/2` reaches `3`.  The scout below should therefore also record whether
the covariance/associator split preserves the cubic cyclotomic defect instead
of merely fitting rational `L_y` margins.

Second rebase adds HYP-3200 as a validation sidecar: in primitive normal form,
consecutive is rank `0/3431` for `Sigma kappa_2`; the all-bank `1/3432` rank
is explained by the nonprimitive dilation twin; no bounded row has
`Sigma kappa_3/S3 = 1/7` exactly; and `kappa_4` is a phi4 stabilizer signal,
not an extremal theorem.  So the scout should target primitive-normal-form
covariance extremality and should keep the odd `1/7`-looking ratio as a
discarded scalar compression, not as proof currency.

Third rebase adds HYP-3201 as the quotient-legality meter.  Ordinary row
entropy remains the wrong extremal scalar, but conditional entropy
`H(target | quotient)` is exactly the right way to detect a lossy compression:
zero means the target function factors through the quotient; positive residual
means a typed sidecar is still owed.  The HYP-3163 scout should therefore emit
law-defect entropy only as a sidecar audit on proposed PSD/Monge/fold/ear
quotients, never as the objective being maximized.

Fourth rebase adds the KPS-S31al Toeplitz/random-current packet.  Its
Caratheodory-Toeplitz readout is the most proof-facing upgrade to Angle A so
far: build the Toeplitz moment matrix `T[j,k]=q_{|j-k|}` from the miss-count
law and test the Caratheodory PSD margin.  KPS reports that consecutive
maximizes `lambda_min(T)` over all `3432` bounded k=8 clusters.  This turns
the Lee-Yang/de Moivre circle into a Szego/Caratheodory-Fejer route.  The
companion ferromagnetic route needs a random-current or coupling-manifold
partial order, because naive speed-path moves toward consecutive are not
monotone in `Sigma kappa_2`.

## Angle A: covariance as a Laplacian/Gram certificate

HYP-3160 says entropy is the wrong extremal language and that the live even
target is

```text
Sigma kappa_2 = sum_{i<j} Cov(X_i, X_j),
X_i = indicator(sector i is empty).
```

HYP-3161 strengthens this from a scan target to an exact finite datum:
`consec_8` maximizes `Sigma kappa_2` over all `3432` bounded k=8 clusters, and
the consecutive rows undergo a ferromagnetic sign transition at `k=5 -> 6`.
The certificate should therefore explain a ferromagnetic ground state, not
rediscover the maximum.

The new attempt is to treat the covariance kernel as proof geometry rather
than as a scalar.  If `K_E(i,j)=Cov_E(X_i,X_j)` is the sector co-emptiness
kernel for a row `E`, then the total covariance is a quadratic form.  The
consecutive row should be attacked by one of the following equivalent-looking
certificates:

```text
PSD certificate:        K_consec - K_E >= 0 on the legal sidecar subspace
Laplacian certificate:  <f, L_E f> >= <f, L_consec f> for all legal f
Monge certificate:      K(d1)+K(d2) >= K(d3)+K(d4) under cyclic uncrossing
Schur certificate:      distance-profile majorization maximizes 1^T K 1 at AP
```

The bold guess is that the exhaustively verified `consec_8` maximum is a disguised
rearrangement theorem: the AP/consecutive distance profile is the unique
ferromagnetic ground state for the pair kernel, and every non-AP row either
has a PSD slack vector or emits named finite-address/observer-gluing debt.

New measurable signals:

```text
covariance_kernel_distance_profile
covariance_matrix_spectrum
PSD_dual_slack_vector
legal_zero_sum_subspace
laplacian_cut_energy
M_matrix_or_stieltjes_status
monge_four_point_defect
cyclic_uncrossing_slack
conditional_negative_type_status
schur_convex_rearrangement_status
ferromagnetic_susceptibility_gap
finite_sidecar_exception_id
cyclotomic_cubic_defect
de_moivre_angle_slack
FM_AFM_bridge_status
primitive_normal_form_status
dilation_twin_exception_id
one_seventh_claim_status
kappa4_stabilizer_sidecar_status
law_defect_entropy_bits
target_function_factor_through_status
sidecar_zero_defect_status
action_determinism_status
toeplitz_lambda_min_margin
caratheodory_psd_margin
szego_fejer_route_status
random_current_coupling_order_status
speed_path_nonmonotonicity_count
```

Immediate proof target: for the HYP-3161/HYP-3200 bounded k=8 bank in
primitive normal form and the HYP-3153 packet columns, compute the legal
sector-distance kernel and test whether every competitor's loss is explained
by a PSD/Monge/Laplacian slack certificate, with the all-bank dilation twin
and any raw PSD failures routed to already named sidecars:
sector-6 boundary leakage,
finite-address observer data, odd-coordinate resurrection, or n=4 canary
compression debt.  If the reduced kernel is Monge or conditionally positive,
the covariance maximum can plausibly be proved without touching the full
degree-six PGF.

## Angle B: Worpitzky residue as associator/odd-ear cocycle

The odd face should not be forced into a variance proof.  HYP-3147 and
HYP-3151 show the small local model: the n=3 two-class kernel has nontrivial
mode `-1/3`, the transitive fiber has a minority edge, and the Worpitzky row
`1,4,1` is lost by scalar score-class compression.  HYP-3160 then identifies
the k=8 odd residue as the non-associative part.

Define an associator signal using empty-sector indicators:

```text
A(i,j,k) =
  E[X_i X_j X_k]
  - E[X_i X_j]E[X_k]
  - E[X_i X_k]E[X_j]
  - E[X_j X_k]E[X_i]
  + 2 E[X_i]E[X_j]E[X_k].
```

This is just a normalized third cumulant candidate.  The proposed proof move
is not that `A` vanishes.  It should not vanish at the hard node.  The move is
to make its boundary visible:

```text
pair covariance = 1-chain / commutative face
triple associator = 2-chain / odd face
Worpitzky minority-edge mode = local boundary defect
odd ear = recursive generator for nonzero associator debt
```

Ear decompositions supply the recursion language: strong connectivity has
ears, factor-critical graphs have odd ears, and nested ears mark
series-parallel tameness.  The conjectural LRC14 use is:

```text
Every nonconsecutive k=8 packet admits either
  (a) an ear-shortening move that lowers odd associator debt without lowering
      even covariance slack, or
  (b) a finite obstruction packet already covered by the n=3/K3 minority-edge,
      n=4 canary/filler, reflection-fold, or observer-gluing sidecars.
```

New measurable signals:

```text
associator_triple_cocycle
associator_boundary_vector
odd_ear_surplus
ear_shortening_move_status
worpitzky_boundary_mode
k3_minority_edge_lift_count
n4_canary_associator_source
even_covariance_odd_associator_exchange
pentagon_associator_defect
terminal_odd_debt_or_discharge
```

This gives a concrete remaining target for the odd part: prove an inequality
of the form

```text
odd_associator_excess(E)
  <= even_covariance_slack(consec, E) + named_finite_sidecar_debt(E).
```

That would let the covariance Laplacian route close the even part while the
odd-ear route closes, bounds, or explicitly residualizes the Worpitzky face.

## Tournament Analysis

Vertex choice is the proof-carrier set, not runners and not raw arcs:

```text
covariance_laplacian_psd_certificate
monge_rearrangement_certificate
sector_distance_kernel_profile
finite_sidecar_exception_packet
associator_triple_cocycle
odd_ear_shortening_certificate
k3_minority_edge_lift
n4_canary_filler_source
lee_yang_root_curve_ledger
raw_scalar_bimodality_value
```

Pairwise observable: for two carriers `A,B`, compare which one preserves the
HYP-3160 target predicate, which one names destroyed coordinates, and which
one can discharge or reduce a proof obligation without scalarizing the root
curve.

Switch/gauge: `A -> B` if `A` retains all data retained by `B` for the current
target and also exposes a certificate or sidecar that `B` forgets.  Ties use
the Hamiltonian path

```text
covariance_laplacian_psd_certificate
-> monge_rearrangement_certificate
-> sector_distance_kernel_profile
-> associator_triple_cocycle
-> odd_ear_shortening_certificate
-> finite_sidecar_exception_packet
-> k3_minority_edge_lift
-> n4_canary_filler_source
-> lee_yang_root_curve_ledger
-> raw_scalar_bimodality_value
```

No executable tournament fingerprint is claimed yet.  The next computation
should build this proof-carrier tournament and report score histogram, SCCs,
directed cycles, edge flips, and Hamiltonian-path count.

## Assumption Challenge

Alternate vertex sets considered: runners, gaps, fixed circle sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
matroid circuits, sector pairs/triples, ear generators, and proof obligations.

Chosen vertices: proof obligations plus sector-pair and sector-triple
certificates.  This preserves the LRC predicate "consecutive maximizes the
k=8 covariance/bimodality packet after legal sidecars" and the named odd
residue.  It destroys exact timing, root angles, and representative identity
unless those are supplied by the Lee-Yang root-curve, finite-address, n=4
canary, or K3 minority-edge sidecars.

Challenged assumption: the useful tournament vertices do not have to be
runners, arcs, or score classes.  At this stage the best vertices are the
functions and proof obligations themselves.

## Next Pull

Build a small k=8 scout that emits the two ledgers side by side:

```text
row_id
q_vector
L_y
Sigma_kappa_2
covariance_kernel_distance_profile
monge_four_point_defect
PSD_dual_slack_vector
associator_triple_cocycle
odd_ear_surplus
even_covariance_odd_associator_exchange
cyclotomic_cubic_defect
de_moivre_angle_slack
FM_AFM_bridge_status
primitive_normal_form_status
dilation_twin_exception_id
one_seventh_claim_status
kappa4_stabilizer_sidecar_status
law_defect_entropy_bits
target_function_factor_through_status
sidecar_zero_defect_status
action_determinism_status
toeplitz_lambda_min_margin
caratheodory_psd_margin
szego_fejer_route_status
random_current_coupling_order_status
speed_path_nonmonotonicity_count
sidecar_exception_id
terminal_exit_or_named_debt
```

First target bank: the primitive bounded k=8 rows already used by HYP-3138,
HYP-3139, HYP-3142, HYP-3160, HYP-3161, and HYP-3200, joined with HYP-3153's
exact Lee-Yang/Worpitzky/quartic packet columns and HYP-3162's cyclotomic
angle defect columns.  The desired readout is not just whether consecutive
wins, but whether every primitive loss decomposes into even PSD/Monge slack
plus odd associator/ear debt while preserving the `q=7` cubic-apex signal and
quarantining the nonprimitive dilation twin.  Any quotient used in this scout
must also report whether its target function has zero HYP-3201 law-defect
entropy or a named typed sidecar.  The Toeplitz page should run in parallel
as the current best spectral proof route: compare `lambda_min(T)` and its
Caratheodory slack against the covariance/Monge/associator decomposition.

## Pointers

HYP-3201, HYP-3200, HYP-3162, HYP-3161, HYP-3160, HYP-3154, HYP-3153, HYP-3152,
HYP-3151, HYP-3150, HYP-3199, HYP-3149, HYP-3148, HYP-3147, HYP-3146, HYP-3142, HYP-3139,
HYP-3138, HYP-3137, HYP-3132, HYP-3124, HYP-3118, THM-577, OPEN-Q-108,
T1219, LTI-280, LTT-178.

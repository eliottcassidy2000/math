---
id: HYP-3081
title: Robbins-Fermat-Catalan branch tournaments and strong proof orientation
status: SYNTHESIS / proof-network guardrail and branch-kernel carrier; not a proof
source: codex-2026-06-26-S249
tangent: T1165
technique: LTI-230
tournament_technique: LTT-128
related:
  - HYP-3079
  - HYP-3078
  - HYP-3077
  - HYP-3076
  - HYP-3075
  - HYP-3074
  - HYP-3073
  - HYP-3072
  - HYP-3071
  - HYP-3070
  - HYP-3069
  - HYP-3068
  - HYP-3067
  - HYP-3066
  - HYP-3058
  - HYP-3057
  - HYP-3056
  - HYP-3054
  - HYP-3053
  - HYP-3052
  - HYP-3051
  - HYP-3050
  - HYP-3049
  - HYP-3048
  - HYP-3047
  - HYP-3009
  - HYP-2997
  - HYP-2995
  - HYP-2990
  - HYP-2981
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3081: Robbins-Fermat-Catalan Branch Tournaments And Strong Proof Orientation

## Claim

Robbins' theorem should be promoted from a metaphor to a proof-network test:
a connected graph admits a strong orientation exactly when it has no bridges.
The LRC14 translation is that a proof quotient is unsafe whenever it creates a
load-bearing bridge whose forgotten coordinate is needed by the next
certificate step.

HYP-3056 makes the missing coordinate an observer-cut orbit.  HYP-3057 says
small tournament numbers must carry value-origin tags before they become proof
objects.  HYP-3058 says Fermat-Catalan reciprocal curvature is a sidecar
pressure, not a proof shortcut.  HYP-3067 and HYP-3068 add median-center,
owner-object, and root-object sidecars before route compatibility is trusted;
HYP-3069 closes named route/certificate states under Boolean median centers.
HYP-3070 adds raw-vs-legal route-triple center control, HYP-3071 turns
the residual cycle-class/F7 coordinate into an observability matrix,
HYP-3072 adds a cross-carrier pullback resonance stub, HYP-3073 reserves
the polymer/Dirichlet-energy bridge rather than an anonymous bucket, HYP-3074
supplies the route-state closure median interface that branch contractions
must preserve, HYP-3075/HYP-3076 mark arithmetic coincidences as retained
sidecars rather than scalar proof shortcuts, HYP-3078 makes q-Pochhammer and
modular-cusp coincidences legal only after finite polar debt is named,
HYP-3079 gives the Lean-facing finite-tail and sixth-power padding interface,
and HYP-3077 shows the completed state hull can be legal while still producing
nonterminal scheduler centers.
The S248 addendum adds the next orientation layer:

```text
proof-safe quotient = bridgeless proof graph + oriented branch corridors
                     + small tournament kernels + finite no-lift guards
```

Once the proof graph is strongly oriented, it should look like branches
connecting small tournament kernels.  The kernels are the local binary
switchboards: AP/GW boundary atoms, K33 state-lift exits, C27 petal/two-block
exits, q=23 Haar squares, A000568 edge-sector decks, residual capacitors,
diagonal-layer rectangle/hourglass cells, and automaton or power-lift guards.
The branches are the proof corridors that carry certificates between them.

Fermat-Catalan enters as the branch-terminal discipline.  Perfect-power and
power-sum coincidences are not proof shortcuts; they are no-lift guards.  A
branch containing a power lift must either strictly decrease a labelled
complexity, land in a finite named exception packet, or carry cyclotomic and
p-adic sidecars proving that the lift cannot proliferate.

## Branch-Tournament Orientation

For a quotient `q:X->Y`, define its proof graph:

```text
Gamma_q = (V_q, E_q)
```

where vertices are proof obligations or local kernels, not runners by default.
Useful vertex choices include:

```text
observer-cut payload orbits
endpoint-owner boundary atoms
Haar zeta squares
residual capacitor plates
K33 state-lift obligations
C27 petal transfer packets
automaton states
Fermat-Catalan exponent guards
diagonal-layer rectangle/hourglass residues
matrix observability columns
Fejer/Ramanujan certificate obligations
```

An undirected edge means that a proof move needs information from one vertex
to validate the other.  An orientation `u -> v` means that the certificate at
`u` can be pushed forward to `v` while preserving the LRC predicate it claims
to preserve.  The orientation is strong if every vertex can reach every other
vertex by certificate-preserving moves.

A **small tournament kernel** is a finite tournament `K` attached to a vertex
or interface of `Gamma_q`.  Its vertices are the possible local exits or
payload choices, and its edge gauge is the declared pairwise proof preference:
which exit preserves more LRC predicate, destroys less sidecar data, or has an
earlier certified discharge.  Its count or isomorphism class is not admissible
unless the HYP-3057 `value_origin_type` has also been declared.

A **branch corridor** is a path between two kernels whose internal vertices are
routine translations: automaton transitions, exact-period refinements,
Haar/Fejer pushes, endpoint-owner transports, sidecar reconstructions, or
state-lift handoffs.  Contracting each kernel to a node and each corridor to
an edge gives the skeleton `Skel(Gamma_q)`.

The proposed admissibility test is:

```text
q is branch-orientation safe only if Skel(Gamma_q) has no naked bridge.
```

A bridge is **naked** when deleting it disconnects the proof graph and neither
side carries a reconstruction, dual annihilator, exact/coboundary discharge,
family descent, AP/GW boundary stop, Fermat-Catalan finite exception, nor
named THM-572/F7 debt.

## Ear-Decomposition Reading

Robbins' proof uses ear decompositions: a bridgeless graph can be built by
starting from a cycle and adding ears between old vertices.  The LRC analogue
is a proof sheaf built from a seed tournament kernel plus oriented ears:

```text
seed kernel
  + oriented branch between two existing kernels
  + oriented branch between an old kernel and a new finite kernel
  + finite residual ear named as debt
```

Every legal ear must satisfy:

1. its endpoint kernels have declared tournament isomorphism classes;
2. its internal branch states preserve the claimed LRC predicate;
3. its reverse verification path is present, reconstructed, dual-certified, or
   named as debt;
4. any power-lift or lacunary transition has a Fermat-Catalan/gap guard; and
5. contracting the ear does not merge AP/GW boundary atoms with positive-open
   packets unless a sidecar separates them.

This suggests a concrete proof strategy.  Instead of asking for one scalar
that separates all rows, enumerate the achievable small tournament kernels and
show every non-AP/GW branch can be oriented into a strict-open or certificate
kernel.  AP/GW are then the only closed boundary seed cycle; all other ears
must attach with a positive exit, a K33/C27/q23/Fejer route, or named debt.

## Fermat-Catalan No-Lift Rule

HYP-3009 already records that the unit-excess lane

```text
M = p/(14p-1)
```

has perfect-power stress points such as:

```text
p=2: q=27=3^3
p=4: p=4=2^2
p=8: p=8=2^3
```

The branch-orientation interpretation is:

```text
power branch = branch corridor with exponent vector debt
```

Such a branch is legal only with fields:

```text
power_lift_guard
fermat_catalan_residual
exponent_vector
cyclotomic_factor_label
p_adic_valuation_label
finite_exception_id
branch_complexity_drop
```

The Fermat-Catalan conjectural shape says that supercritical exponent
coincidences should be finite after common factors are controlled.  In this
repo we cannot use that as a black box proof.  We can use it as a necessary
condition: no LRC branch is allowed to hide an infinite family of power-lift
coincidences behind a scalar automatic word or denominator value.

## Necessary Conditions

1. **No naked proof bridge.**  If a quotient makes one dependency edge the
   only route between two certificate regions, the missing coordinate must be
   retained or discharged.
2. **Kernel endpoint declaration.**  Every branch endpoint must name a small
   tournament kernel: its vertex set, pairwise observable, gauge, tie path, and
   isomorphism class or achievability set.
3. **Two-way certificate transport.**  A forward certificate push must have a
   reverse verification path: reconstruction, dual annihilation, exactness,
   descent, boundary stop, or residual debt.
4. **No-lift power discipline.**  Perfect-power, power-sum, lacunary, or
   automatic coincidences require exponent, p-adic, cyclotomic, and finite
   exception sidecars before they can be used as branch contractions.
5. **Controlled forgetting commutes with orientation.**  Passing to a quotient
   must preserve strong orientation of the proof graph after all discharged
   edges are contracted.  If orientation fails, the quotient has created a
   Robbins bridge.
6. **Achievable-kernel subset.**  For each small `n`, record the subset
   `Ach_n(q)` of tournament isomorphism classes realizable as branch kernels.
   A proof by tournament isomorphism classes needs this subset, not only the
   full A000568 count.
7. **Bridge-to-residual naming.**  If a bridge cannot be oriented into a
   reversible branch, it must become a named residual: AP/GW boundary, C27,
   K33/THM-572, q=23 capacitor, Fejer/Ramanujan handoff, or new F7.

## Tournament Analysis

The vertices for this pass are branch carriers, not runners:

```text
observer_cut_payload_orbit
Robbins_no_bridge_assembly
small_tournament_kernel
endpoint_owner_closed_H1
residual_capacitor_cut
K33_state_lift_branch
C27_petal_transfer_branch
q23_Haar_square_branch
diagonal_rectangle_hourglass_branch
Fermat_Catalan_no_lift_guard
automaton_gap_branch
Fejer_Ramanujan_certificate_branch
raw_scalar_shadow
```

Pairwise observable:

```text
preserves_boundary_open_status
preserves_route_handoff
preserves_owner_current
retains_or_reconstructs_observer_cut_payload
keeps_strong_orientation_after_quotient
blocks_infinite_power_lifts
has_named_residual_exit
```

Switch:

```text
A -> B
```

when `A` preserves every proof predicate preserved by `B`, destroys fewer
observer-cut or branch-orientation coordinates, or has an earlier certified
discharge at equal preservation.

One Hamiltonian tie path is:

```text
observer_cut_payload_orbit
> Robbins_no_bridge_assembly
> small_tournament_kernel
> endpoint_owner_closed_H1
> residual_capacitor_cut
> K33_state_lift_branch
> C27_petal_transfer_branch
> q23_Haar_square_branch
> diagonal_rectangle_hourglass_branch
> Fermat_Catalan_no_lift_guard
> automaton_gap_branch
> Fejer_Ramanujan_certificate_branch
> raw_scalar_shadow
```

Assumption challenge: the vertices could instead be gaps, fixed circle
sections, section boundaries, wall crossings, residues, cover arcs, Fourier
modes, matroid circuits, branch ears, endpoint interfaces, or proof
obligations.  The branch-kernel quotient preserves boundary/open status, route
handoff, owner current, observer-cut payload, and named residual exits.  It
destroys raw runner labels, raw phase location, and scalar magnitude unless
those are reattached as sidecars.

## Next Pull

Build a branch-kernel ledger over HYP-2963 coarse fibers:

```text
base_quotient
proof_graph_vertex_set
branch_id
endpoint_kernel_left
endpoint_kernel_right
kernel_iso_class_left
kernel_iso_class_right
achievable_kernel_subset
oriented_forward_certificate
reverse_verification_mode
bridge_status
observer_cut_payload_orbit_id
power_lift_guard
finite_exception_id
residual_debt_name
```

The first concrete target is to take the S220 observer-cut orbit ledger and
ask whether each mixed route/status fiber becomes bridgeless after adding the
known sidecars.  Any remaining bridge is a proof obligation, not an analogy.

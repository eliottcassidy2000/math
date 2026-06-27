---
id: HYP-3083
title: LRC14 Hurwitz finite-address branch-closure spine
status: SYNTHESIS / finite-address covering and branch-closure remaining-obligation map; not a proof
source: codex-2026-06-27-S252
updated: codex-2026-06-27-S252b+S253
tangent: T1167
technique: LTI-232
tournament_technique: LTT-130
related:
  - HYP-2900
  - HYP-2906
  - HYP-2909
  - HYP-2996
  - HYP-2997
  - HYP-2963
  - HYP-3082
  - HYP-3081
  - HYP-3080
  - HYP-3079
  - HYP-3078
  - HYP-3077
  - HYP-3075
  - THM-523
  - THM-525
  - THM-534
  - THM-566
  - THM-570
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-3083: LRC14 Hurwitz Finite-Address Branch-Closure Spine

This file completes the S252/S252b/S253 synthesis lane for the user's request
to integrate incoming work, past LRC14 proof work, Hurwitz's theorem, and the
q-Pochhammer/full-modular-group cusp cue into a sharper account of what remains
in the LRC14 proof.

The claim is not that a Hurwitz theorem or modular q-series directly proves
LRC14.  The claim is that all Hurwitz-shaped lessons point to the same
proof-interface rule:

```text
an infinite tail, extremal scalar, or branch contraction is legal only after
the finite address that generates it has been retained or discharged
```

The relevant Hurwitz readings are:

- Diophantine Hurwitz/Markov: a sharp approximation constant is explained by a
  finite continued-fraction or Markov address, not by the scalar alone.
- Complex-analytic Hurwitz: zeros/poles cannot appear silently under locally
  uniform holomorphic limits; a q-Pochhammer limit must retain zero/pole
  divisor data.
- Hurwitz/Vieta surfaces: an infinite Diophantine orbit is legal only after
  its finite seed and mutation word are retained.

For the full modular group, a modular function must be invariant under the
modular transformations and meromorphic at the cusp `i infinity`.  With
`q=exp(2*pi*i*z)`, this means its Laurent q-expansion has only finitely many
negative-power terms.  In LRC14 language, a quotient may use an infinite
product, partition, periodicity, or mutation tail only after its polar debt is
finite and every polar term has a named packet exit.

Incoming work supplies adjacent pieces:

- HYP-3075: Hurwitz/Markov/Pell/cannonball and modular-cusp finite-address
  sidecars.
- HYP-3077: legal route-state median hulls can exist while remaining
  nonterminal scheduler centers.
- HYP-3078/HYP-3079: q-series tails are legal only after finite principal
  part, transformation, multiplier, and zero/pole persistence data are named.
- HYP-3080: sixth-power equalities are proof data only with lane tuple,
  primitive gcd, shared-term filter, residue words, and collision exit.
- HYP-3081/HYP-3082: proof quotients must not create naked Robbins bridges;
  the current HYP-2963 packet-bank protected graph has no naked bridges after
  sidecars are attached.
- HYP-2996/HYP-2963: the bounded packet bank already has exact `M`,
  q-threshold, route, residual-section, and Haar-grid exits for 21,913 rows,
  with zero anonymous hard non-AP/GW remainder.

The synthesis target is:

```text
global primitive row
  -> finite-address packet bank
  -> protected branch graph
  -> terminal discharge or named residual debt
  -> formal witness readout M >= 1/14
```

## Covering-Bound Priority

The incoming S59 redirect is correct in its main prioritization: the proof
should not spend its critical path on the AP/Goddyn-Wong census.  For the
non-strict LRC14 theorem, `t=1/14` already handles every row with no multiple
of `14`, because every nonzero residue mod `14` has distance at least `1/14`
from an integer at that time.  THM-523 strengthens this to the full
small-divisor covering condition for any strict counterexample.  The proof is
therefore a covering bound.

Past work sharpens that slogan.  After THM-571, the covering bound is no
longer the whole "some multiple of 14" world: rows with at least seven
multiples of `14` are discharged by the apex-majority gamma descent.  After
HYP-2906, a single top speed satisfying `v_max > 13 v_second` is discharged by
a connected-safe-arc tooth comparison, and HYP-2900 supplies the larger Node-3
equidistribution/induction picture.  THM-566 rules out a naive uniform
bounded-denominator witness.

Thus the live residual is:

```text
primitive covering 13-row
  with 1..6 multiples of 14,
  no one-large-speed ratio peel,
  no small-q witness,
  no apex-majority gamma descent,
  and no legal quotient that forgets the finite address of its obstruction.
```

That is much smaller, but still open.

## Integrated Proof Spine

The improved argument should be stated as a finite-address exponential covering
spine:

```text
S primitive, |S|=13
  -> q-witness gate
  -> covering row with a multiple of 14
  -> apex-majority gamma descent or low-apex residue
  -> one-large-speed interval peeler or top-balanced residue
  -> finite-address packet P(S)
  -> protected branch graph
  -> terminal discharge or named residual debt
  -> formal witness readout M(S) >= 1/14
```

The word "exponential" keeps the useful part of the S59 hyperoperation
reframe: multiples of `14` are periodic on the `1/14` lift, and the gamma trick
uses that periodicity.  Hurwitz/Markov/q-cusp adds the guardrail: periodicity,
q-series tails, route medians, power equalities, and scalar extrema are legal
only after the finite address that made them special is retained.

The finite address should contain at least:

```text
P(S) = (exact scale, multiple-of-14 profile, owner/topology, route,
        primitive period deck, safe/open status, analytic certificate,
        endpoint shell, cusp/principal-part debt, relation/collision sidecars,
        branch kernel, destroyed-coordinate exit, formal status).
```

The missing theorem is not "find one scalar separator."  It is:

```text
Every primitive low-apex, top-balanced covering residual admits P(S), and
P(S) preserves the LRC predicate M(S) >= 1/14 until a terminal witness,
boundary equality, contradiction, or named residual debt is reached.
```

## Theorem-Facing Gates

1. **Witness gate.**  If a row omits a multiple of some `q in {2,...,14}`,
   THM-523 gives the constructive `1/q` witness.  AP/Goddyn-Wong equality
   atoms live on this easy side; they explain the tight floor but are not the
   hard covering proof.

2. **Apex-majority gate.**  If at least seven speeds are multiples of `14`,
   THM-571 gives `M(S)>1/14` using the `14`-phase sieve and, when necessary,
   the `7`-phase gamma descent.  This removes the old many-14-multiples
   residual.

3. **Scale-separation gate.**  If the top speed is larger than `13` times the
   second speed, HYP-2906 peels it by a connected seed-safe interval.  HYP-2900
   is the wider analytic Node-3 version: large committed speeds remove the
   expected `1/7` of a positive seed-safe set, pending effective
   Erdos-Turan/Weyl bounds.

4. **Finite-address packet gate.**  The remaining low-apex, top-balanced
   covering residue must emit a packet retaining exact `M`, q-threshold,
   binding denominators, safe topology, endpoint owners, route labels, C27/K33
   status, Ramanujan/q-cusp sidecars, sixth-power/relation sidecars, and
   destroyed-coordinate exits.  This is HYP-2963 plus the
   HYP-3075/HYP-3078/HYP-3079/HYP-3080 finite-address discipline.

5. **Protected branch graph.**  HYP-3082 changes the next theorem obligation
   from "find one scalar" to "prove every primitive residual reaches the
   protected graph and then discharge each protected exit."  A naked bridge is
   a missing sidecar/debt, not a usable theorem shortcut.

6. **Terminal exits.**  The exits are direct q-witness, strict safe component
   or covering-moment discharge, C27/petal discharge, AP/GW boundary equality,
   K33/THM-572 state-lift contradiction, or named formal residual debt.
   THM-572 is already Lean-formalized as a closure theorem once the missing
   LRC-to-tournament state lift is constructed.

## q-Pochhammer / Modular-Cusp Merge

HYP-3078 and HYP-3079 are the q-expansion discipline for the covering branch.
The exact scouts verify:

```text
(q;q)_infty = 1 - q - q^2 + q^5 + q^7 - q^12 - ...
1/(q;q)_infty = 1,1,2,3,5,7,11,15,22,...
Delta(q) = q (q;q)_infty^24
j(q) = q^-1 + 744 + 196884 q + ...
```

The proof import is not the coefficient list.  It is the sidecar rule:

```text
full modular function packet
  = SL2Z transformation law
  + eta/multiplier balance when needed
  + finite cusp principal part
  + named polar exits
  + certified nonpolar tail
```

Thus an LRC packet carrying q-Pochhammer, eta/Delta, partition, divisor, or
modular-function language must record:

```text
q_cusp_ledger_id
modular_transform_status
eta_multiplier_balance_status
principal_part_order
finite_negative_power_budget
principal_part_coeff_vector
polar_exit_word
q_pochhammer_tail_signature
partition_tail_certificate
log_derivative_divisor_channel
hurwitz_zero_persistence_status
illegal_infinite_polar_tail_flag
```

The LRC analogue of "finite principal part" is "finite named obstruction
debt."  AP/GW boundary equality, q-witness exits, C27 owner-strip exits,
covering-moment exits, and K33/THM-572 state-lift debt are legal polar exits.
An infinite negative q-tail corresponds to an uncontrolled infinite family of
undischarged obstructions and is not a proof object.

## Hurwitz Arithmetic Merge

The same rule governs the arithmetic side.  Hurwitz/Markov/Pell data do not
directly produce an anti-Bohr time.  They name exceptional walls where a scalar
is sharp because a finite address is present: continued-fraction period,
Markov depth, Pell unit, endpoint shell, carry residue, or Vieta mutation word.
HYP-3075's Markov-Pell/cannonball scout is the model:

```text
70^2 = 1^2 + ... + 24^2
70 = Pell P6
29 * 169 - 70^2 = 1
```

The visible square is not the proof coordinate.  The proof coordinate is the
quadratic-unit address and carry between neighboring Markov-Pell wall numbers.
Likewise, HYP-3075's Hurwitz surface orbit is legal because it keeps the seed
`(2,2,2,2)` and Vieta mutation rule.

For LRC14 covering rows, arithmetic sidecars should record:

```text
hurwitz_markov_approximant_class
lagrange_markov_depth
continued_fraction_period_word
pell_wall_unit
pell_cassini_gap
endpoint_shell_address
quadratic_carry_residue
hurwitz_vieta_seed
hurwitz_mutation_depth
destroyed_owner_or_scale_coordinate
required_sidecar_or_exit
```

The transfer is:

```text
best-approximant wall
  -> finite address
  -> endpoint/cusp/owner packet sidecar
  -> legal branch graph entry or named residual debt.
```

## What Remains

### O1. Global finite-address normalizer

Prove that every primitive low-apex, top-balanced covering row either reaches
the HYP-2963/HYP-3082 packet bank or reaches a declared outside-bank normal
form with the same finite-address fields.  This is the missing global
reduction.  It must not quotient by raw residue, raw tournament class,
automatic word, scalar `M` bucket, or q-series tail unless the forgotten
coordinate is reconstructible, dual-annihilated, fiber-constant, or routed to
named residual debt.

### O2. Covering-moment / positive-open discharge

Prove the strict safe-component exit for covering-family packets uniformly. In
the older language this is OPEN-Q-108 / GAP G2 / bounded-core fattening.  In
the Lean formalization status it is the analytic side of the cap floor
`hmeasGP` plus finite-ruler Part A `hpartA` (or the p0/cap positive-margin
route after its certificate/extremality bridge is supplied).  In S59 language,
this is the exponential covering bound, not the additive census.

The q-cusp analogy is useful here because it says the infinite covering tail is
legal only after its finite polar debt is named.  A covering proof should
therefore expose a finite certificate row:

```text
principal debt: covered small-denominator witnesses
tail certificate: p0/moment, gamma/apex-periodic, Haar/Ramanujan, or Node-3
exit: positive witness density, descent, or named residual
```

### O3. K33 / THM-572 state-lift construction

Turn each K33/nonunit residual packet into the complete tournament-state lift
required by THM-572, or give a different terminal discharge.  THM-572 proves
the endpoint contradiction from an `H=7` lift; it does not build the lift.
HYP-3082's `K33-STATE-LIFT` rows are therefore honest debt, not solved rows.

### O4. Branch-closure / no-new-naked-bridge theorem

Promote the HYP-3082 bounded audit into a theorem over the normalizer image:
controlled forgetting is legal exactly when every load-bearing bridge is
protected by a named sidecar, reconstruction, dual-annihilation, descent, or
residual exit.

For parameterized tails, the family theorem should prove:

```text
family extension of a protected branch
  either keeps a reverse verification path
  or descends
  or discharges
  or names new residual debt
```

This is exactly where the modular finite-principal-part rule and Hurwitz seed
rule become proof infrastructure: parameterized infinite tails are legal only
if they are generated by a finite address.

### O5. Integer-vs-real and formal sidecar closure

Formalize the passage from limiting shape/cap arguments back to integer
runners: `rho_K -> rho*`, arc-count/Vmax bounds, shape extraction, and witness
handoff.  The Lean track has a verified conditional assembly reducing LRC14 to
named nodes, and THM-572 is formal once the state lift exists.  What is not
yet formal is the packet normalizer, the covering-moment certificate, the K33
lift, and the sidecar legality theorem saying finite addresses survive every
proof quotient.

The next Lean-facing object should be a record combining:

```text
packet identity
finite principal-part / polar-exit proof
branch graph entry and bridge certificate
terminal discharge certificate
event bridge to goodSet/safeSet or state lift
Part-A or THM-572 readout
```

This is the natural common target for `LRCModularCuspLedger`,
`LRCMedianCenterControl`, `LRCMomentDual`, and `LRCTournamentStateLift`.

S254 instantiates this common target as
`TournamentH7.LRCFiniteAddressBranchClosure`.  The Lean theorem
`lrc14_from_cutting_edge_branch_coverage` proves the exact conditional frontier:
early q-witness / level-7 lift sieve / one-large-speed gates plus low-apex,
top-balanced finite-address packet coverage imply `LRC14Statement`.  After
incoming THM-573, the packet retains both the original multiple-of-14 covering
datum and the sharpened residual count `1..6` multiples of `7`.  The packet now
has formal fields for q-cusp finite principal part, q-Pochhammer tail, Hurwitz
arithmetic address, protected branch certificate, covering-moment dual ledger,
median-center packet, and terminal floor `1/14 <= witnessFloor <= Mreach`.
Incoming HYP-3085-gK8 suggests the covering-moment producer should be a
low-order pairwise `S2` / reflection-`3x3` Perron certificate, and incoming
HYP-3085-covering/K33 supplies concrete nested-refinement versus cross-handoff
rows for the O2/O3 shuttle.  Incoming HYP-3087 supplies the adjacent
hyperoperation/operation-grid address scheduler; it should feed this packet
only after retaining `(p,q)`, danger, level-7 status, destroyed coordinate, and
terminal exit.  No analytic normalizer, covering discharge, or K33 lift is
asserted; those are the remaining producer theorems.

### O6. Boundary census only if the proof routes through equality

The AP/GW census and shell-collapse problem remain mathematically valuable,
but they are not the shortest proof path unless the argument chooses a
boundary-forcing route.  A direct covering bound can bypass full tight-locus
classification.  If a future proof uses "all tight rows are AP/GW" as its
bridge, then it must still prove the HYP-2909 shell-collapse / covering
strictness step; otherwise the census is a side theorem.

## Candidate Branch-Closure Lemma

The theorem shape suggested by this synthesis is:

> **Finite-address branch-closure lemma.** Let `S` be a primitive 13-speed row.
> Suppose the q-witness, apex-majority, and one-large-speed gates reduce `S` to
> a low-apex, top-balanced covering row.  Suppose this row admits a packet
> `P(S)` whose finite-address ledger retains exact scale, endpoint owner,
> route, residual section, bridge-protection sidecar, and any q-cusp/Hurwitz
> arithmetic address used by the proof.  If `P(S)` has no illegal infinite
> polar tail, no naked protected-graph bridge, and every terminal branch is
> discharged by q-witness, boundary, owner-strip, covering-family, or THM-572
> state-lift exit, then `M(S) >= 1/14`.

This lemma is intentionally conditional.  It is a better target than another
scalar separator because each missing hypothesis is now a named proof
obligation.

## Tournament Analysis

Candidate vertices considered: runners, residues, multiples-of-14 blocks,
gaps, fixed circle sections, section boundaries, wall-crossing events, cover
arcs, Fourier modes, q-series coefficients, modular functions, q-cusp
principal parts, Hurwitz/Markov/Pell addresses, Vieta seeds, endpoint owners,
K33 packets, bridge ears, median centers, matroid wall data, and proof
obligations.

Chosen vertices are proof obligations and finite-address branch/cusp/address
packets.  This quotient preserves the LRC predicate `M(S)>=1/14`, exact
q-threshold, multiple-of-14 status, scale ratio, boundary/open status,
endpoint owner, finite polar debt, Hurwitz arithmetic address, bridge
protection, and terminal exit.  It destroys raw runner identity, raw
q-coefficients, raw scalar constants, and naked route labels only after those
losses are irrelevant, reconstructed, dual-annihilated, or named as residual
debt.

Pairwise observable:

```text
retained finite address
covering-branch relevance
low-apex/top-balanced status
finite principal part / no infinite polar debt
endpoint owner and exact scale
protected bridge status
terminal exit specificity
formalization readiness
```

Switch/gauge: orient toward the carrier that retains more of the active
covering-bound proof predicate; tie by fewer destroyed coordinates and earlier
terminal discharge.

Synthesis Hamiltonian path:

```text
finite_address_packet_normalizer
> protected_branch_graph
> covering_moment_exit
> K33_state_lift_exit
> modular_cusp_principal_part_guard
> hurwitz_arithmetic_seed_guard
> route_state_median_scheduler
> apex_majority_gamma_gate
> one_large_interval_peeler
> q_witness_gate
> Lean_sidecar_closure
> raw_scalar_or_census_shadow
```

Known and intended fingerprints:

- HYP-3082 raw scalar-star: `6` nodes, `5` bridges, `5` naked bridges.
- HYP-3082 protected graph: `80` nodes, `83` edges, `69` bridges, `0` naked
  bridges, strongly orientable contracted core.
- HYP-3080 sixth-power rank audit: `0` nontrivial two-lane collisions through
  `250`, `5` three-lane collision certificates through `80`, and `5`
  rank-sensitive edge flips.
- The intended proof-carrier tournament is transitive: singleton score
  histogram, no directed 3-cycles, singleton SCCs, and one Hamiltonian path.
  A directed cycle in an actual packet audit would be useful evidence that two
  carriers each destroy a coordinate the other needs.

## Practical Next Pull

Build the `finite_address_branch_closure` ledger promised by LTT-130.  Each
row should have:

```text
source_row_id
source_row_or_family
contains_multiple_of_14
low_apex_top_balanced_status
normalizer_exit_attempted
finite_address_kind
finite_address_word
exact_M_and_q_threshold
endpoint_owner
destroyed_coordinate
required_sidecar_or_debt
q_cusp_ledger_id
principal_part_order
polar_exit_word
hurwitz_or_pell_address
protected_branch_node
bridge_status_raw
bridge_status_protected
bridge_protection_mode
median_center_kind
terminal_exit
lean_formalization_status
residual_debt_name
```

Concrete first targets:

1. Instantiate one low-apex covering-moment row as a q-cusp ledger with finite
   principal part, polar exit word, product/partition tail certificate, branch
   id, bridge status, and covering discharge route.
2. Turn the p0/moment covering-family certificate into exact rational data:
   feasible dual over miss counts, `L_y <= cap - delta`, and the event bridge
   into `goodSet/safeSet`.
3. Build the K33 state-lift packet record with endpoint owners, incidence
   branch kernel, and the `H(T)=7` payload required by THM-572.
4. Re-run the S250 branch audit on a parameterized low-apex covering family,
   not just the bounded bank, and report whether any bridge becomes naked.
5. Add a Lean `FiniteAddressBranchPacket` record that packages q-cusp,
   branch-bridge, terminal-discharge, and witness-readout fields in one place.

The first useful target is the low-apex covering-moment family, not AP/GW
census.  The second is the K33/THM-572 lift construction.  Those two
obligations are the live proof debt after integrating the incoming redirect
with the older gamma, peeler, packet-bank, branch-closure, q-cusp/Hurwitz, and
Lean tracks.

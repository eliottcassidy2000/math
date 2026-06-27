---
id: HYP-3083
title: LRC14 Hurwitz finite-address branch-closure spine
status: SYNTHESIS / remaining-obligation map; not a proof
source: codex-2026-06-27-S252
tangent: T1167
technique: LTI-232
tournament_technique: LTT-130
related:
  - HYP-2900
  - HYP-2906
  - HYP-2909
  - HYP-3082
  - HYP-3081
  - HYP-3080
  - HYP-3079
  - HYP-3078
  - HYP-3077
  - HYP-3075
  - HYP-2963
  - THM-523
  - THM-525
  - THM-566
  - THM-570
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-3083: LRC14 Hurwitz Finite-Address Branch-Closure Spine

This file completes the S252 synthesis lane for integrating incoming work,
past LRC14 proof work, and the Hurwitz finite-address lesson into a sharper
account of what remains in the LRC14 proof.

The core rule is not that Hurwitz's approximation theorem directly proves
LRC14.  The transferable Hurwitz/Markov/Pell lesson is:

```text
an extremal scalar can be used only after its finite address is retained
```

That same rule now applies to q-cusp tails, sixth-power equalities, median
scheduler centers, proof-graph kernels, endpoint owners, apex-periodic
covering rows, and exact packet addresses.  Any quotient that forgets the
coordinate making a row boundary/open distinguishable must reconstruct it,
dual-annihilate it, prove it fiber-constant, or route it to named residual
debt.

## Critical Reprioritization

The incoming S59 redirect is correct in its main prioritization: the proof
should not put the AP/Goddyn-Wong census on the critical path.  For the
non-strict LRC14 theorem, if no runner is a multiple of `14`, then `t=1/14`
already gives every runner distance at least `1/14`.  The hard case is the
covering case where at least one runner is a multiple of `14`.

Past work refines the covering case:

- THM-523 gives the small-divisor / q-witness gate.
- THM-570 and THM-571 close the apex-majority branch with at least seven
  multiples of `14`.
- HYP-2906 peels a row when `v_max > 13 v_second`.
- HYP-2900 supplies the larger Node-3 equidistribution / induction picture.
- THM-566 rules out a naive fixed bounded-denominator witness atlas.
- HYP-3082 shows that raw scalar quotienting creates naked proof bridges, but
  protected sidecars remove them on the current HYP-2963 bank.

Thus the live residual is:

```text
primitive covering 13-row
  with 1..6 multiples of 14,
  no small-q witness,
  no apex-majority gamma descent,
  no one-large-speed interval peel,
  and no legal quotient that forgets the finite address of its obstruction.
```

That is much smaller than the full census program, but still open.

## Integrated Evidence

1. **Finite address is the common currency.**  HYP-3075 turns
   Hurwitz/Markov/Pell/cannonball coincidences and modular-cusp data into
   address ledgers: continued-fraction period word, Markov/Pell wall, cusp
   principal part, mutation seed, endpoint shell, and required exit.  HYP-3078
   and HYP-3079 make the q-series version precise: positive tails can be
   compressed, but polar debt must be finite and named.  HYP-3080 gives the
   sixth-power version: a two-lane equality is a rigidity gate, while a
   three-lane equality is a native support-six resonance only after lane tuple,
   primitive gcd, shared-term filter, and residue words are retained.

2. **The packet bank already knows many safe exits.**  HYP-2963 reports a
   bounded labelled-packet bank with `21913` packets, no below-threshold rows,
   AP/GW as the zero-open tight rows, and `7235` hard non-AP/GW packets.
   HYP-3024/HYP-3030/HYP-3036 explain why coarse residue/unit gates alone are
   insufficient and why topology plus primitive period decks split the mixed
   residual fibers.  THM-526 warns that a universal single-removal covering
   criterion is false, so the eventual cover must be adaptive and addressed.

3. **Branch closure is the quotient test.**  HYP-3081 imports Robbins'
   no-bridge criterion as a proof-network guardrail: a forgotten coordinate is
   illegal if it becomes the only load-bearing bridge between proof regions.
   HYP-3082 makes this executable on the HYP-2963 bank.  The raw scalar-star
   graph has `6` nodes, `5` bridges, and `5` naked bridges.  After route
   sections, Haar/grid exits, no-lift guards, q-cusp obligations, finalizer
   gates, and named residual exits are protected, the branch graph has `80`
   nodes, `83` edges, `69` bridges, `0` naked bridges, and a strongly
   orientable contracted core.  This is strong evidence for the interface, not
   yet a proof of global coverage.

4. **The formal endpoint is conditional, not missing by syntax.**  THM-572
   blocks a tournament-state packet with `H(T)=7`, but only after the LRC
   residual has been lifted to the required state.  The Lean formalization
   status file keeps the same shape: the sorry-free scaffolding needs analytic
   `hmeasGP`/cap-floor input and `hpartA` finite-ruler/witness handoff, plus
   concrete instantiation of shape, exit, and packet predicates.

## Improved Proof Spine

The improved argument should be stated as a finite-address exponential
covering spine:

```text
S primitive, |S|=13
  -> q-witness gate
  -> covering row with a multiple of 14
  -> apex-majority gamma descent or low-apex residue
  -> one-large-speed interval peeler or top-balanced residue
  -> finite-address packet P(S)
  -> protected branch graph Gamma(P)
  -> strict witness / AP-GW boundary / C27 petal /
     covering moment / K33-THM-572 lift / named residual debt
  -> formal sidecar closure
```

The word "exponential" keeps the useful part of the S59 hyperoperation
reframe: a multiple of `14` is periodic on the `1/14` lift, and the gamma
trick exploits that periodicity.  The Hurwitz/Markov/q-cusp branch supplies
the guardrail: periodicity, q-series tails, route medians, power equalities,
and scalar extrema are legal only after the finite address that made them
special is retained.

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
   THM-523 gives the constructive `1/q` witness.  AP/GW equality atoms live on
   this easy side; they explain the tight floor but are not the hard covering
   proof.

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

Prove the strict safe-component exit for covering-family packets uniformly.
In older language this is OPEN-Q-108 / GAP G2 / bounded-core fattening.  In
the Lean formalization status it is the analytic side of `hmeasGP` plus the
finite-ruler Part A `hpartA` (or the p0/cap positive-margin route after its
certificate/extremality bridge is supplied).  In S59 language, this is the
exponential covering bound, not the additive census.

### O3. K33/state-lift construction

Turn each K33/nonunit residual packet into the complete tournament-state lift
required by THM-572, or give a different terminal discharge.  THM-572 proves
the endpoint contradiction from an `H=7` lift; it does not build the lift.
HYP-3082's `K33-STATE-LIFT` rows are therefore honest debt, not solved rows.

### O4. Branch-closure theorem

Promote the HYP-3082 bounded audit into a theorem over the normalizer image:
controlled forgetting is legal exactly when every load-bearing bridge is
protected by a named sidecar, reconstruction, dual-annihilation, or residual
exit.

### O5. Integer-vs-real and formal sidecar closure

Formalize the passage from limiting shape/cap arguments back to integer
runners: `rho_K -> rho*`, arc-count/Vmax bounds, shape extraction, and witness
handoff.  The Lean track has a verified conditional assembly reducing LRC14 to
named nodes, and THM-572 is formal once the state lift exists.  What is not
yet formal is the packet normalizer, the covering-moment certificate, the K33
lift, and the sidecar legality theorem saying finite addresses survive every
proof quotient.

### O6. Boundary census only if the proof routes through equality

The AP/GW census and shell-collapse problem remain mathematically valuable,
but they are not the shortest proof path unless the argument chooses a
boundary-forcing route.  A direct covering bound can bypass full tight-locus
classification.  If a future proof uses "all tight rows are AP/GW" as its
bridge, then it must still prove the HYP-2909 shell-collapse / covering
strictness step; otherwise the census is a side theorem.

## Tournament Analysis and Assumption Challenge

Do not take runners or arcs as the default vertices.  Candidate vertex sets
considered here were runners, gaps, residues, multiples-of-14 blocks, fixed
circle sections, section boundaries, wall-crossing events, cover arcs, Fourier
modes, Haar rectangles, matroid circuits, endpoint owners, exact-period
packets, q-cusp ledgers, sixth-power collision certificates, K33 packets, and
proof obligations.

The chosen tournament vertices are proof obligations / finite-address branch
carriers.  This preserves the LRC predicate `M(S) >= 1/14`.  It destroys raw
row identity and raw scalar shadows only after the lost coordinate is recorded
as irrelevant, reconstructed, dual-annihilated, or named as residual debt.

Pairwise observable:

```text
retained q-threshold, multiple-of-14 profile, exact scale, boundary/open
status, route, owner topology, primitive period deck, analytic certificate,
branch protection, and formal exit status
```

Switch/gauge:

```text
A -> B when A preserves the target predicate with fewer unnamed destroyed
coordinates, or when A exposes a protected terminal exit that B leaves as a
raw scalar or census shadow.
```

Tie Hamiltonian path for the proof-obligation tournament:

```text
finite_address_packet_normalizer
  > protected_branch_graph
  > covering_moment_exit
  > K33_state_lift_exit
  > apex_majority_gamma_gate
  > one_large_interval_peeler
  > q_witness_gate
  > Lean_sidecar_closure
  > boundary_census_shadow
  > raw_scalar_shadow
```

Known fingerprints from incoming audits:

- HYP-3082 raw scalar-star: `6` nodes, `5` bridges, `5` naked bridges.
- HYP-3082 protected graph: `80` nodes, `83` edges, `69` bridges, `0` naked
  bridges, strongly orientable contracted core.
- HYP-3080 sixth-power rank audit: `0` nontrivial two-lane collisions through
  `250`, `5` three-lane collision certificates through `80`, and `5`
  rank-sensitive edge flips.

Under the displayed switch the intended proof-obligation tournament is
transitive.  A directed cycle in an actual packet audit would be useful
evidence that two carriers each destroy a coordinate the other needs.

## Next Artifact

Build the `finite_address_branch_closure` ledger promised by LTT-130.  Each
row should have:

```text
source_row_or_family
low_apex_top_balanced_status
normalizer_exit_attempted
finite_address_word
preserved_lrc_predicate
destroyed_coordinate
required_sidecar_or_debt
protected_branch_node
bridge_status_raw
bridge_status_protected
terminal_exit
lean_formalization_status
residual_debt_name
```

The first useful target is the low-apex covering-moment family, not AP/GW
census.  The second is the K33/THM-572 lift construction.  Those two
obligations are the live proof debt after integrating the incoming redirect
with the older gamma, peeler, packet-bank, branch-closure, and Lean tracks.

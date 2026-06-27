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

This file completes the S252 synthesis lane for the user's request to integrate
incoming work, past LRC14 proof work, and the Hurwitz finite-address lesson into
a sharper account of what remains in the LRC14 proof.

The intended claim is not that Hurwitz's approximation theorem directly proves
LRC14.  The intended claim is that the sharp Hurwitz/Markov lesson should be
used as a proof-interface rule:

```text
an extremal scalar can be used only after its finite address is retained
```

Incoming work has already produced the adjacent pieces:

- HYP-3075: Hurwitz/Markov/Pell/cannonball and modular-cusp finite-address
  sidecars.
- HYP-3077: legal route-state median hulls can exist while remaining
  nonterminal scheduler centers.
- HYP-3078/HYP-3079: q-series tails are legal only after finite principal
  part, transformation, and zero/pole persistence data are named.
- HYP-3081/HYP-3082: proof quotients must not create naked Robbins bridges;
  the current HYP-2963 packet-bank protected graph has no naked bridges after
  sidecars are attached.

The synthesis target is to turn these into a sharper remaining-obligation map:

```text
global primitive row
  -> finite-address packet bank
  -> protected branch graph
  -> terminal discharge or named residual debt
```

The incoming S59 redirect is correct in its main prioritization: the proof
should not spend its critical path on the AP/Goddyn-Wong census.  For the
non-strict LRC14 theorem, `t=1/14` already handles every row with no multiple
of `14`, and THM-523 strengthens this to the full small-divisor covering
condition for any strict counterexample.  The proof is therefore a covering
bound.

The past work sharpens that slogan.  After THM-571, the covering bound is no
longer the whole "some multiple of 14" world: rows with at least seven multiples
of `14` are discharged by the apex-majority gamma descent.  After HYP-2906,
a single top speed satisfying `v_max > 13 v_second` is discharged by a
connected-safe-arc tooth comparison, and HYP-2900 supplies the larger Node-3
equidistribution/induction picture.  THM-566 rules out a naive uniform
bounded-denominator witness.  Thus the live residual is:

```text
primitive covering 13-row
  with 1..6 multiples of 14,
  no one-large-speed ratio peel,
  no small-q witness,
  no apex-majority gamma descent,
  and no legal quotient that forgets the finite address of its obstruction.
```

That is much smaller, but still open.

## Integrated proof spine

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
```

The phrase "exponential" here keeps the useful part of the S59 hyperoperation
reframe: multiples of `14` are periodic on the `1/14` lift, and the gamma trick
uses that periodicity.  The Hurwitz/Markov/q-cusp branch adds the guardrail:
periodicity, q-series tails, route medians, and scalar extrema are legal only
after the finite address that made them special is retained.

The theorem-facing pieces now line up as follows.

1. **Witness gate.**  If a row omits a multiple of some `q in {2,...,14}`,
   THM-523 gives the constructive `1/q` witness.  The AP/Goddyn-Wong equality
   atoms live on this easy side; they explain the tight floor but are not the
   hard covering proof.

2. **Apex-majority gate.**  If at least seven speeds are multiples of `14`,
   THM-571 gives `M(S)>1/14` using the `14`-phase sieve and, when necessary,
   the `7`-phase gamma descent.  This removes the old "many 14-multiples"
   residual.

3. **Scale-separation gate.**  If the top speed is larger than `13` times the
   second speed, HYP-2906 peels it by a connected seed-safe interval.  HYP-2900
   is the wider analytic Node-3 version: large committed speeds remove the
   expected `1/7` of a positive seed-safe set, pending effective
   Erdos-Turan/Weyl bounds.

4. **Finite-address packet gate.**  The remaining low-apex, top-balanced
   covering residue must emit a packet retaining exact `M`, q-threshold,
   binding denominators, safe topology, endpoint owners, route labels, C27/K33
   status, Ramanujan/q-cusp sidecars, and destroyed-coordinate exits.  This is
   HYP-2963 plus the HYP-3075/HYP-3078/HYP-3079 finite-address discipline.

5. **Protected branch graph.**  HYP-3082 shows that on the current HYP-2963
   packet bank, the raw scalar star has naked bridges but the protected
   branch graph has none after route/section/grid/no-lift/residual sidecars
   are attached.  This changes the next theorem obligation from "find one
   scalar" to "prove every primitive residual reaches the protected graph and
   then discharge each protected exit."

6. **Terminal exits.**  The exits are direct q-witness, strict safe component
   or covering-moment discharge, C27/petal discharge, AP/GW boundary equality,
   K33/THM-572 state-lift contradiction, or named formal residual debt.  THM-572
   is already Lean-formalized as a closure theorem once the missing
   LRC-to-tournament state lift is constructed.

## What remains

The remaining proof obligations are now precise.

### O1. Global finite-address normalizer

Prove that every primitive low-apex, top-balanced covering row either reaches
the HYP-2963/HYP-3082 packet bank or reaches a declared outside-bank normal
form with the same finite-address fields.  This is the missing global
reduction.  It must not quotient by raw residue, raw tournament class,
automatic word, scalar `M` bucket, or q-series tail unless the forgotten
coordinate is reconstructible, dual-annihilated, fiber-constant, or routed to
named residual debt.

### O2. Covering-moment family discharge

Prove the strict safe-component exit for covering-family packets uniformly.
In the older language this is OPEN-Q-108 / GAP G2 / bounded-core fattening.  In
the Lean formalization status it is the analytic side of the cap floor
`hmeasGP` plus finite-ruler Part A `hpartA` (or the p0/cap positive-margin
route after its certificate/extremality bridge is supplied).  In S59 language,
this is the exponential covering bound, not the additive census.

### O3. K33/state-lift construction

Turn each K33/nonunit residual packet into the complete tournament-state lift
required by THM-572, or give a different terminal discharge.  THM-572 proves
the endpoint contradiction from an `H=7` lift; it does not build the lift.
HYP-3082's `K33-STATE-LIFT` rows are therefore honest debt, not solved rows.

### O4. Shell-collapse only if the proof routes through equality

The AP/GW census and shell-collapse problem remain mathematically valuable,
but they are not the shortest proof path unless the argument chooses a
boundary-forcing route.  A direct covering bound can bypass full tight-locus
classification.  If a future proof uses "all tight rows are AP/GW" as its
bridge, then it must still prove the HYP-2909 shell-collapse / covering
strictness step; otherwise the census is a side theorem.

### O5. Formal sidecar closure

The Lean track has a verified conditional assembly reducing LRC14 to named
nodes, and THM-572 is formal once the state lift exists.  What is not yet
formal is the packet normalizer, the covering-moment certificate, the K33 lift,
and the sidecar legality theorem saying finite addresses survive every proof
quotient.

## Tournament Analysis

The tournament vertices for this synthesis are proof obligations, not runners:

```text
q_witness_gate
apex_majority_gamma_gate
one_large_interval_peeler
finite_address_packet_normalizer
protected_branch_graph
covering_moment_exit
K33_state_lift_exit
Lean_sidecar_closure
raw_scalar_or_census_shadow
```

Pairwise observable: which carrier preserves the target predicate
`M(S)>=1/14` while retaining exact address data: q-threshold, multiple-of-14
status, scale ratio, safe topology, endpoint owners, route labels, q-cusp
principal part, bridge protection, and residual exit.

Switch/gauge: orient `A -> B` when `A` preserves all data preserved by `B` and
strictly more address or terminal-discharge information; ties follow the
Hamiltonian path

```text
finite_address_packet_normalizer
> protected_branch_graph
> covering_moment_exit
> K33_state_lift_exit
> apex_majority_gamma_gate
> one_large_interval_peeler
> q_witness_gate
> Lean_sidecar_closure
> raw_scalar_or_census_shadow
```

Fingerprint: the intended proof-carrier tournament is transitive; expected
score histogram `{0:1, ..., 8:1}`, no directed 3-cycles, singleton SCCs, and
one Hamiltonian path.  A directed cycle in an actual packet audit would be a
useful warning that two carriers each destroy a coordinate the other needs.

## Assumption challenge

The synthesis explicitly rejects the default assumption that the useful
vertices are runners or arcs.  Candidate vertex sets considered here were
runners, residues, multiples-of-14 blocks, gaps, fixed circle sections,
section boundaries, wall-crossing events, cover arcs, Fourier modes, q-cusp
principal parts, endpoint owners, K33 packets, matroid wall data, and proof
obligations.  The chosen quotient uses proof obligations and finite-address
packets because that quotient preserves the LRC predicate and exposes exactly
which coordinate is destroyed by each attempted simplification.  It destroys
raw runner identity and most raw scalar shadows, but only after those losses
are either irrelevant to `M(S)>=1/14`, reconstructed by sidecars, or named as
residual debt.

## Practical next pull

Build the `finite_address_branch_closure` ledger promised by LTT-130.  Each
row should have:

```text
source row or family
low-apex/top-balanced status
finite address retained
destroyed coordinate if quotienting
protected branch node
bridge protection mode
terminal exit
Lean/formalization status
residual debt name
```

The first useful target is the low-apex covering-moment family, not AP/GW
census.  The second is the K33/THM-572 lift construction.  Those two obligations
are the live proof debt after integrating the incoming redirect with the older
gamma, peeler, packet-bank, and Lean tracks.

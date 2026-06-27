---
id: HYP-3095
title: LRC14 route sophistication and the hidden observability sheaf
status: SYNTHESIS / proof-route abstraction; not a proof
source: codex-2026-06-27-S256
tangent: T1172
technique: LTI-236
tournament_technique: LTT-134
related:
  - HYP-2976
  - HYP-2990
  - HYP-3083
  - HYP-3085
  - HYP-3087
  - HYP-3088
  - HYP-3089
  - HYP-3090
  - HYP-3091
  - HYP-3092
  - HYP-3093
  - HYP-3094
  - THM-573
  - THM-574
  - THM-575
  - THM-576
  - OPEN-Q-108
---

# HYP-3095: LRC14 Route Sophistication and the Hidden Observability Sheaf

## Claim

The repo's LRC14 route has not merely accumulated more cases.  It has become
more precise about what is legal to forget.  The hidden structure we keep
nearly naming is an **observability sheaf over quotient maps**:

```text
row S
  -> arithmetic observers       (mod 14, mod 7, c-lifts, null polynomial)
  -> arc observers              (safe components, endpoint owners, normalized floor)
  -> cap observers              (pairwise avoidance, triangular ratios)
  -> moment observers           (p0, gK8, low-order S2/Perron modes)
  -> branch observers           (finite address, protected bridge, K33 handoff)
  -> formal witness readout     (M(S) >= 1/14 or named residual debt)
```

Each observer sees a different projection of the same packet.  LRC14 remains
hard because no single projection is the object: the proof must show that the
observers glue without any quotient silently destroying the witness coordinate.

Equivalently, the live conjectural theorem is not "find the right scalar" but:

```text
Every primitive covering 13-row has a finite-address observer family whose
arithmetic, arc, moment, and branch restrictions either glue to a witness at
height 1/14 or expose a named residual debt that is already routed.
```

## Route history as increasing sophistication

### 1. Additive/tight-locus era

The early AP/GW and three-gap work treated the frontier as an additive
geometry problem.  It learned the equality atoms and built census intuition.
That was real information, but later S59/S60 work corrected its role: AP/GW is
the extremal face of the easy non-covering half, not the critical covering
proof.

The first hidden-structure lesson was therefore negative: an elegant
tight-locus classification can be mathematically true while still watching the
wrong half of the theorem.

### 2. Tournament and relational era

The route then moved from scalar tuple invariants to pairwise relations:
tournament spectra, magnitude-aware flips, route labels, owner strips, and
coarse-vs-exact fibers.  This made failures more informative.  A mixed fiber
was no longer "noise"; it identified a missing sidecar.

The main advance was role separation.  The best splitter, the cheapest proof
tooth, and the final route certificate are not required to be the same object.
This is why later work can use a moment functional, an endpoint owner, and a
K33 lift in one proof packet without pretending that one of them subsumes the
others.

### 3. Packet/sidecar era

HYP-2963 and its descendants changed the unit of proof from a runner set to a
finite-address packet.  Exact scale, endpoint owner, boundary/open status,
safe-component stalks, primitive period decks, and residual capacitor exits
became first-class data.

This is the era where the repo learned the governing meta-rule:

```text
a quotient is proof-safe only if the forgotten coordinate is reconstructible,
annihilated by the dual certificate, constant on fibers, or routed to named debt.
```

HYP-2990 later names this as the zipper/no-free-slider principle.

### 4. Finite-address branch era

The q-cusp, Hurwitz, Pell, sixth-power, median, and branch-kernel threads look
eclectic only at scalar level.  At proof-object level they say the same thing:
rare tail behavior is usable only after its finite address and terminal exit
are retained.

HYP-3083 condensed this into the current spine:

```text
primitive covering row
  -> finite-address packet
  -> protected branch graph
  -> terminal discharge or named residual debt
  -> formal witness readout
```

This is the sheaf language in primitive form.  Local data can be compressed
only after the overlap maps are known.

### 5. Covering/exponential redirect

S59 changed the target: LRC14 is the covering bound, not the census.  The
important half is the row containing a multiple of 14, and the exponential /
apex-periodic/gamma route is the proof core.

This turned the older four Farey variations into a real coordinate chart:

```text
additive lane       AP/GW census, three-gap, equality shadow
product lane        14 = 2*7, lift/descent, valuation status
exponential lane    apex periodicity, gamma trick, Node-3 discharge
2-adic lane         Clebsch/cut-space, low-order moment carrier
```

The hidden object is the packet that can be read in all four lanes.

### 6. Broken-field algebra era

The incoming polynomial-method work sharpened why 14 is special.  The paper's
Proposition 4.1 wants `Z/(k+1)Z` to be a field.  At `k=13`, `k+1=14=2*7`;
the field proof breaks into prime-factor lifts.  THM-574 shows the c-lift
family, and THM-573 supplies the level-7 sieve.  HYP-3089/S61 makes the
apex bridge exact: `I(13,7,1)` is covering mod `7`; THM-573 is the `c=7`
lift; and the dyadic `c=2` lift lands on the project's covering condition
mod `14`.  The polynomial-method reflection adds the algebraic shadow: the
units mod 14 have size 6, leaving a degree-7 null direction, while the
low-order moments survive through degree 6.

That is the same number appearing in different languages:

```text
7 lifts at level 7
7 forbidden/copime boundary after THM-573
degree-7 null polynomial after mod 14 stops being a field
S2 / low-order moment carrier before the missing direction
```

The live proof must control the direction the field proof cannot see.

### 7. Raw Conjecture 7.1 failure and normalized repair

THM-575 is not a setback; it is a correction of the observer.  Raw denominator
time is the wrong clock because divisor-loaded apex rows can kill every small
denominator while remaining non-tight.  This failure explains why the witness
route must be normalized in slow/ruler coordinates after the level-7 sieve.

In sheaf terms:

```text
raw time components are not stable under apex loading;
normalized core components are the object that can glue to I(13,7,1),
c=2 lifts, finite-ruler sampling, and moment discharge.
```

The Conjecture 7.1 idea survives only after changing the observed coordinate.

## What is hidden underneath LRC14

The hidden structure is a **controlled-forgetting geometry**.  Its local
coordinates are familiar:

- residue/lift coordinates: mod 14, mod 7, c=2, c=7, `I(13,7,1)`;
- arc coordinates: safe components, endpoint owners, direct and normalized
  witness floors, the `V*` crossover;
- moment coordinates: `p0`, gK8, low-order `S0..S4`, pairwise `S2`, the
  reflection-symmetric Perron mode;
- branch coordinates: finite address, protected bridge, residual exit,
  nested refinement, cross-handoff, K33 lift;
- formal coordinates: integer-vs-real finite-ruler glue and Lean-side records.

The reason the route feels overdetermined is that these are not competing
metaphors.  They are local charts on the same proof object.  We struggle to
describe it perfectly because each chart loses one coordinate the next chart
was invented to restore.

## Current proof obligations through this lens

S258 adds the first concrete starter artifact:
`04-computation/lrc14_observer_gluing_ledger_codex_s258.py` and
`05-knowledge/results/lrc14_observer_gluing_ledger_codex_s258.out`.  It glues
the HYP-3096 direct-witness ledger and the HYP-3097 pair/Pascal scissors fields
on representative rows.  The sample has q-witness and THM-573 terminal exits,
plus seven live residual rows whose direct arcs, CRT status, mod-7 scissors
signatures, and Farey lanes must be kept as separate observer fields.  The
large-apex and divisor-loaded rows show why raw direct largest arcs cannot be
the global scalar chart.

1. **Normalized apex residual.**  Build the THM-573 / THM-575 residual ledger
   for rows with `<=6` multiples of 7.  It should retain `P,E,V`, normalized
   good set `G(P,E)`, component count, largest normalized component,
   `I(13,7,1)=covering mod 7`, `c=7` lift status, dyadic `c=2` lift status
   to covering mod `14`, finite-ruler threshold, and terminal exit.

2. **Degree-7 null-direction control.**  Treat gK8/HYP-3085 as the low-order
   moment chart for the same residual, not as a separate scalar bound.  The
   target is a reflection-Perron / `S2` control statement strong enough to
   bound the concentration that the mod-14 polynomial method cannot eliminate.

3. **Cap chart integration.**  Treat incoming HYP-3090 as the cap /
   pairwise-avoidance chart: the triangular ratios tell us which covering
   cap constants are geometric pair-avoidance shadows and where the remaining
   deviations live.  This chart should constrain, not replace, the normalized
   arc and moment charts.

4. **Shuttle without scalar collapse.**  Use HYP-3094 to keep the
   covering-moment and K33 obligations from merging illegally.  Nested
   refinement and cross-handoff may both have positive safe mass, so the
   branch packet must retain grid class, active binders, and endpoint-owner
   transitions.

5. **Gluing theorem.**  Prove that the arithmetic chart, normalized arc chart,
   cap chart, moment chart, and branch chart have compatible overlaps on the
   residual.  This is the abstract theorem the repo keeps approaching under
   different names: no-free-slider, finite-address branch closure,
   no-naked-bridge, normalized Conjecture 7.1, and c-lift descent.

## Tournament Analysis

This pass deliberately does not use runners as tournament vertices.  The
vertices are proof observers:

```text
finite_address_observability_sheaf
normalized_level7_apex_peel
gK8_moment_perron_chart
covering_K33_shuttle
hyperoperation_grid_address
AP_GW_census_shadow
raw_denominator_pruning
```

Pairwise observable: observer `A` beats observer `B` if `A` preserves every
LRC predicate needed by `B` while also retaining at least one coordinate that
`B` forgets.  Ties are broken by whether the observer names a terminal exit.

Synthetic fingerprint:

```text
score histogram: {6:1, 5:1, 4:1, 3:1, 2:1, 1:1, 0:1}
directed cycles: 0
SCCs: 7 singletons
Hamiltonian paths: 1 under the declared gauge
tie path:
  finite_address_observability_sheaf
  > normalized_level7_apex_peel
  > gK8_moment_perron_chart
  > covering_K33_shuttle
  > hyperoperation_grid_address
  > AP_GW_census_shadow
  > raw_denominator_pruning
```

Interpretation: the current route is not a tournament among runners, arcs, or
residues.  It is a dominance order among allowed observers.  The top object is
not numerically strongest; it is logically least forgetful.

## Assumption challenge

Candidate vertex sets considered: runners, gaps, fixed circle sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
matroid circuits, proof obligations, and observer charts.

Rejected for this synthesis:

- runners: too coarse for c-lift and endpoint-owner debt;
- gaps/sections: good for AP/GW, weak for loaded apex rows;
- residues: good for THM-573/THM-574, blind to normalized arc mass;
- Fourier modes: good for gK8, blind to branch/K33 handoff unless sidecars are
  retained;
- proof obligations alone: too abstract unless each obligation declares which
  LRC predicate it preserves and which coordinate it destroys.

Chosen vertex set: observer charts.  This preserves the LRC predicate "there
exists a witness at height `1/14` or a legal terminal exit" and records the
destroyed coordinate explicitly.  The challenged assumption is that LRC14 has
one natural carrier.  The repo history suggests the carrier is multi-charted;
the proof should be a gluing argument over legally controlled quotients.

## Reframe

The best current slogan is:

```text
LRC14 is a broken-field proof repaired by observability.
```

The field method breaks at 14 and leaves a degree-7 invisible direction.  The
project has spent its history inventing observers for that invisible direction:
safe components, endpoint owners, level-7 lifts, low-order moments, finite
addresses, and branch exits.  The remaining proof is to show that these
observers are not ad hoc repairs but compatible charts on one object.

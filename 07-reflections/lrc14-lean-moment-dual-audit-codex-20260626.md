# LRC14 Lean Moment-Dual Audit - Codex 2026-06-26

status: formalization audit
target: `04-computation/lean/TournamentH7/TournamentH7/LRCMomentDual.lean`
checked: `lake build TournamentH7.LRCMomentDual`

This pass moved the THM-534 moment-dual wrapper from "written but not
build-verified" to checked Lean.  The theorem
`TournamentH7.MomentDual.p0_le_Ly` now proves, from pointwise feasibility

```text
g(n) >= 1[n=0] for n <= 6,
```

that the concrete cover atom obeys

```text
(slowμ (coverSet E)).toReal <= L_y(E,g)
```

where `L_y(E,g) = integral g(missCount E x) dslowμ`.  The axiom printout is
only `[propext, Classical.choice, Quot.sound]`; there is no `sorryAx` in this
local reduction.

## What This Closes

The checked statement is the honest integral step behind the p0 route:
`coverSet E` is exactly the zero-miss event, dual feasibility bounds the
indicator pointwise, and monotonicity of the integral gives the p0 upper bound.
That means Lean is no longer asking us to trust the cover-to-moment relaxation.

This does not prove LRC14.  It says the remaining moment side is now sharply
located: provide a feasible finite dual certificate and a verified inequality
`L_y(E,g) <= cap - delta` for the relevant shape family, and the cover atom
feeds the existing p0/witness floor machinery.

## How Close We Really Are

The proof DAG is close; the unconditional Lean theorem is not one local lemma
away.  Current checked pieces include the concrete `p0`, dense-cover inclusion,
the p0-to-witness floor algebra, finite-ruler Part A arithmetic glue, and now
the moment-dual relaxation.  The hard remaining Lean-facing bridges are:

1. A certificate/extremality bridge: exact rational data for feasible duals
   `g : {0,...,6} -> Q`, plus a proof that the induced `L_y` is at most the
   appropriate cap with positive margin in the binding family.
2. A concrete event bridge from the skeleton's opaque `shapeOf` and
   `witnessG2` to the existing `coverSet`, `safeSet`, `goodSet`, and finite
   phase/ruler objects.
3. The Part-A analytic readout: positive finite witness density must be
   converted to `Mreach >= 1/14` for the actual 13-speed family, not merely for
   an abstract shape.

So the next Lean improvement should not be another abstract sidecar.  It should
be a small certificate record, probably with fields for finite support,
feasibility over `0..6`, exact moment evaluation or an upper bound on `L_y`, and
the margin comparison against `cap`.  `p0_le_Ly` can consume that record
directly.

Incoming-mainline connection: HYP-3069/HYP-3070/HYP-3071/HYP-3073 arrived while
this Lean pass was open.  Their route-center, center-control, cycle-class, and
Dirichlet/current language all point to the same formal interface requirement:
every quotient must name the retained predicate, the destroyed coordinate, and
the certificate column that repairs the loss.  For the Lean moment route, that
means the dual certificate should carry not just a scalar bound, but a named
certificate column tied to `missCount`, cap margin, and the p0-to-witness
handoff.  Otherwise the proof repeats the raw-route projection failure noted in
HYP-3070: a scalar looks centered only because the legal sidecar page was
forgotten.

## Tournament Analysis

For this formalization pass, runners are the wrong tournament vertices.  The
useful vertices are proof obligations and measurable events:

```text
coverSet_zero_miss
dual_feasibility
p0_le_Ly
Ly_cap_certificate
p0_margin_to_witnessG2
finite_ruler_error
PartA_reach
```

Pairwise observable: between two vertices `A` and `B`, record which one removes
an opaque hypothesis needed by the other, and secondarily which one exposes an
exact rational certificate instead of an analytic estimate.

Switch/gauge: orient `A -> B` when formalizing `A` strictly reduces the
hypotheses of `B`.  Ties are broken by preferring vertices whose proof term is
finite/rational over vertices that still require measure or equidistribution.

Tie Hamiltonian path:

```text
coverSet_zero_miss
  -> dual_feasibility
  -> p0_le_Ly
  -> Ly_cap_certificate
  -> p0_margin_to_witnessG2
  -> finite_ruler_error
  -> PartA_reach
```

Fingerprint after this pass: `p0_le_Ly` has left the unchecked SCC and is now a
verified downstream edge from `dual_feasibility` to the p0 route.  The largest
remaining unchecked block is the certificate/extremality vertex together with
the concrete `shapeOf`/event bridge.  The Hamiltonian path above is acyclic at
the proof-obligation level; cycles only appear if we collapse events back to
raw runner sets and lose the miss-count coordinate.

Assumption challenge: I considered runners, gaps, fixed circle sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
matroid circuits, and proof obligations as vertices.  The chosen quotient is
miss-count events plus proof obligations.  It preserves the LRC predicate only
through the chain

```text
p0 margin -> positive witnessG2 -> finite witness density -> Mreach >= 1/14.
```

It destroys endpoint owners and exact lonely-time location.  That destruction
is acceptable for the moment-dual certificate layer, but it would be fatal for
Part A unless the finite-ruler/event bridge reintroduces the lost phase
geometry.

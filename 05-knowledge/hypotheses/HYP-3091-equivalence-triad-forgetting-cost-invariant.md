---
id: HYP-3091
title: Equivalence triad forgetting-cost invariant
status: SYNTHESIS / invariant program; no new theorem claimed
source: codex-2026-06-27-S257
tangent: T1171
tags: [invariants, equinumerosity, equidecomposability, equidistribution, controlled-forgetting, lrc14, tournament-analysis]
related:
  - HYP-2187
  - HYP-2186
  - HYP-2244
  - HYP-2232
  - HYP-2872
  - HYP-2883
  - HYP-2895
  - HYP-2899
  - HYP-2949
  - HYP-3053
  - HYP-3054
  - HYP-3056
  - HYP-3072
  - HYP-3085
  - HYP-3088
  - HYP-3089
  - HYP-3090
  - THM-576
  - OPEN-Q-108
---

# HYP-3091: Equivalence Triad Forgetting-Cost Invariant

This pass searches the older workspace for a more abstract invariant than
raw equidistribution.  The recurring pattern is a triad:

```text
equinumerosity       = same count / same orbit number / same cardinal shadow
equidecomposability = same retained piece-and-sidecar fiber
equidistribution    = legal limiting forgetting after resonance debt is paid
```

The proposed invariant is not any one of these.  It is the cost of moving
between them.

Given a quotient or analogy

```text
q : X -> Y
```

and a target predicate `P`, attach the tuple

```text
F_q(x; P) =
  (cardinal_shadow,
   scissors_fiber,
   observer_cut_orbit,
   distribution_law,
   interaction_order_defect,
   named_residual_debt).
```

This is the "forgetting-cost" invariant.  It measures what has to be retained,
reconstructed, annihilated, descended, boundary-stopped, or named as residual
debt before a quotient becomes proof-safe.

## Why This Exists

The past work has many count coincidences that become useful only after a
side channel is restored.

1. **Royle/even-graph equinumerosity.**  HYP-2187 shows that equal counts are
   a cardinal shadow.  The usable bridge would need to preserve `H`, `beta1`,
   and odd-cycle packet data.  HYP-2872/HYP-2883 then sharpen the warning:
   even graphs are cycle-space addresses, not an `H`-complete obstruction
   closure.

2. **Tournament scissors data.**  HYP-2186 reads `H` as volume and the
   strong-component multiset as the scissors class.  Same `H` can split into
   different piece multisets, so `H` is equinumerosity-like while the
   component multiset is equidecomposability-like.  The `beta1` and packet
   refinements are Dehn-style side channels.

3. **CH/cardinal-shadow method.**  HYP-2232 imports the set-theory lesson:
   cardinal data is incomplete without the model/generic side channel.  In
   finite repo language, a residue count, orbit count, or Hamiltonian-path
   value is not absolute unless the lift/fiber in which it lives is named.

4. **Fixed-path half-tilings.**  HYP-3053 makes A000568 into a cover by
   fixed Hamiltonian path presentations.  The fiber over a class is
   `H(T)/|Aut(T)|`, and rectangle/hourglass cycle laws say when duplicated
   line observations can be compressed to potentials.  The count is the
   shadow; presentation multiplicity and cycle defects are the retained cost.

5. **Baire/Haar events.**  HYP-2949 separates regular-open bulk from finite
   endpoint debt.  Positive Haar mass, endpoint-only residuals, and exact
   LRC loneliness are different predicates until boundary owners are retained.

6. **Covering caps and pairwise avoidance.**  Incoming HYP-3090 and THM-576
   make the cap side sharper: for small `j`, the minimum lonely mass is a
   pairwise avoidance probability, while the binding `k=8,9` rows are exactly
   where order-3 and higher corrections appear.  This is a clean
   interaction-order defect: pairwise equidistribution is almost the right
   invariant, and the remaining constants measure where it stops being exact.

7. **Carrier portfolios.**  HYP-3072 says no small universal scalar portfolio
   covers all active LRC14 obligations.  The useful object is local:
   preserve the LRC predicate, name the destroyed coordinate, and glue local
   portfolios through legal exits.

## Candidate Invariant Fields

### 1. Cardinal Shadow

The visible count:

```text
orbit count, class count, cardinality, H value, cap value, fiber size,
raw safe mass, denominator count, support size.
```

It is useful for navigation and false-positive detection.  It is not a proof
object unless the remaining fields are either irrelevant to `P` or discharged.

### 2. Scissors Fiber

The piece data retained by an equidecomposition:

```text
strong-component multiset,
(H,beta1,packet) fiber,
fixed-path presentation orbit,
support-six packet graph,
owner-strip current,
endpoint/deletion fiber,
rectangle/hourglass cycle class,
route/certificate sidecar set.
```

This is the part that says whether two same-volume or same-count objects can
really be cut and reassembled inside the allowed symmetry group.

### 3. Observer-Cut Orbit

For a quotient `q` and next operation `o`, HYP-3056 suggests:

```text
C_q(x,o) = orbit_Aut_q(x)(boundary slice, incidence word, extended shadow).
```

This is the next-operation version of a Dehn invariant.  If two points have
the same `q(x)` but different `C_q(x,o)` and `P` changes after extension, then
the quotient has forgotten load-bearing information.

### 4. Distribution Law

The limiting law that makes forgetting legal:

```text
Haar regular-open measure,
Weyl or Erdos-Turan rate,
resonance lattice rank and covolume,
primitive-period deck,
squarefree blindness report,
pairwise co-emptiness matrix,
Fourier tail / Abel summation law.
```

Equidistribution is not "ignore the fiber."  It is a certificate that the
forgotten fiber becomes harmless after deleting or naming resonances.

### 5. Interaction-Order Defect

The first order at which a lower-order shadow stops predicting the predicate:

```text
S2 exact but S3/S4 correction active,
pairwise avoidance exact but order-3 cap deviation nonzero,
cycle-space address correct but H nonmonotone,
automatic word route-mixed until exact magnitude or stalk is attached,
line-potential quotient valid until a rectangle/hourglass residue survives.
```

This field is where hidden perspectives tend to live.  It is often more
stable than the scalar value itself.

### 6. Named Residual Debt

The terminal honest label:

```text
AP/GW boundary atom,
K33/THM-572 state lift,
F7/harmonic residual,
Node-3 equidistribution debt,
finite cap deviation constant,
unproved largest-arc normalized-floor producer,
missing formal packet normalizer.
```

The invariant is complete only when every non-distributed defect has a name.

## LRC14 Readout

For LRC14 the triad suggests the following proof posture.

```text
do not ask:
  Which scalar proves the row?

ask:
  Which count shadow is being used?
  Which scissors/fiber sidecars make it predicate-preserving?
  Which resonances prevent equidistribution?
  At what interaction order does the pairwise shadow fail?
  What residual debt remains after legal forgetting?
```

This folds together the current mainline:

- HYP-3088 refutes raw absolute-denominator Conjecture 7.1 and repairs it by
  normalized arc carriers.
- HYP-3089 maps the paper's `I(13,7,1)` bridge to covering mod `7/14` and the
  `V*` crossover.
- HYP-3085/THM-576 say the cap side is controlled by pairwise co-emptiness
  until explicit higher-order deviations appear.
- HYP-3054/HYP-3056 say every quotient must name its next observer-cut orbit.

So the current fundamental invariant is not "uniform distribution" alone.  It
is:

```text
uniform distribution modulo the retained scissors fiber and its first
undistributed observer-cut defect.
```

## Tournament Analysis

Tournament vertices should be invariant carriers, not runners, arcs, raw
classes, or raw denominator charts.

Candidate vertices:

```text
forgetting_cost_tuple
observer_cut_orbit
scissors_fiber
interaction_order_defect
distribution_law
resonance_lattice
boundary_bulk_split
presentation_multiplicity
cardinal_shadow
raw_scalar
```

Pairwise observable:

```text
preserves target predicate,
separates mixed fibers,
survives next operation,
has exact/coboundary/dual/descent discharge,
names residual debt,
detects distribution failure,
low proof cost.
```

Switch/gauge:

```text
A -> B
```

when `A` preserves every predicate-changing distinction preserved by `B` and
adds at least one of: next-operation stability, interaction-order detection,
or a named legal exit.  Ties follow:

```text
forgetting_cost_tuple
> observer_cut_orbit
> scissors_fiber
> interaction_order_defect
> distribution_law
> resonance_lattice
> boundary_bulk_split
> presentation_multiplicity
> cardinal_shadow
> raw_scalar.
```

## Assumption Challenge

This pass explicitly rejects the assumption that the right invariant must be
one scalar, one orbit count, one runner tournament, or one equidistribution
statement.

Alternative vertex sets considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing
events, residues, cover arcs, Fourier modes, Haar tiles, exact-period decks,
Burnside cycle types, fixed-path presentations, observer-cut orbits,
sidecar columns, matroid/cycle-space addresses, proof obligations, and
residual debt names.
```

The chosen quotient preserves only this predicate:

```text
the ability to decide whether a proposed proof quotient is legal for the
target mathematical statement.
```

It destroys direct proof content unless the tuple is instantiated on actual
packet families.  That is why HYP-3091 is a synthesis and protocol, not a
proof of LRC14 or of any cross-domain theorem.

## Next Pull

Build a small `equivalence_triad_probe` over existing HYP-2963 or tournament
rows.  For each candidate quotient, emit:

```text
target_predicate
quotient_name
cardinal_shadow
fiber_id
scissors_fiber_key
observer_kind
observer_cut_orbit_id
distribution_law_id
resonance_lattice_key
interaction_order_first_failure
mixed_fiber_count
separating_sidecar
discharge_mode
residual_debt_name
```

The first useful finite test is not a broad new enumeration.  It is a
side-by-side ledger for three already-known collisions:

1. Royle/even-graph count versus `(H,beta1,packet)` fibers.
2. AP/GW endpoint-only boundary rows versus positive regular-open rows.
3. HYP-3090/THM-576 pairwise cap rows versus the `k=8,9` higher-order
   deviation constants.

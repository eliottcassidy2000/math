---
source: codex-2026-06-01-S539
status: computation plus outside-box mapping atlas
tags:
  - lonely-runner
  - tournament-analysis
  - event-media
  - holes
  - gates
  - kinetic-tournaments
  - circular-arc-covers
  - observer-anchoring
---

# Event-Media Tournaments for LRC

The sector-dual mapping HYP-2022 already answers the user's direct suggestion:
take the `n` evenly spaced sectors as the tournament nodes and let edges change
when runners cross sector boundaries.  This session asks for the next step:
what if the tournament nodes are not runners, not sectors, and not even
positions?  What if the nodes are the *media through which the LRC obstruction
moves*?

That shifts the design question from:

```text
Which objects should be compared?
```

to:

```text
Which finite medium carries the obstruction, and what pairwise relation on that
medium becomes a tournament?
```

This is a helpful discomfort.  A tournament is an absurdly simple object:
between any two nodes, exactly one arrow.  LRC is not simple: it is a one
parameter orbit, a covering problem, a discrepancy problem, a sector occupancy
problem, a Fourier resonance problem, and an endpoint-debt problem.  The art is
choosing a pairwise relation that destroys almost everything except the
obstruction.

## Computation

Artifact:

```text
04-computation/lrc_event_media_tournaments_s539.py
05-knowledge/results/lrc_event_media_tournaments_s539.out
```

The probe uses exact open sector-crossing cells.  Good means:

```text
c_0(t) = c_{n-1}(t) = 0
```

so the two sectors touching the observer are empty.

It compares six functors:

```text
hole_only_bare
  nodes = currently empty sectors only;
  edge = longer survival as a hole beats shorter survival;
  no observer anchoring.

hole_only_anchored
  same nodes, but a hole remembers distance from the observer clasp.

sector_survival_bare
  nodes = all sectors;
  edge = empty/long-lived/low-occupancy sector beats the other;
  no observer anchoring.

sector_survival_anchored
  same tournament, but observer-adjacent sectors and empty status are colored.

gate_priority_bare
  nodes = boundary gates k/n;
  edge = gate with earlier next crossing has higher priority.

gate_priority_anchored
  same gate tournament, but observer gates and adjacent empty statuses are
  colored.
```

Tournament Analysis declaration:

```text
pairwise observables:
  hole survival time, sector occupancy, gate next-crossing priority

switch/gauge:
  bare isomorphism vs observer-anchored coloring

tie Hamiltonian path:
  cyclic sector/gate order

fingerprints:
  class counts and good-only/bad-only/mixed fiber counts
```

The result:

```text
n=4
  hole_only_bare             classes=3   mixed=2
  hole_only_anchored         classes=14  mixed=0
  sector_survival_bare       classes=2   mixed=1
  sector_survival_anchored   classes=24  mixed=0
  gate_priority_bare         classes=3   mixed=1
  gate_priority_anchored     classes=78  mixed=0

n=5
  hole_only_bare             classes=4   mixed=3
  hole_only_anchored         classes=31  mixed=0
  sector_survival_bare       classes=3   mixed=1
  sector_survival_anchored   classes=61  mixed=0
  gate_priority_bare         classes=6   mixed=2
  gate_priority_anchored     classes=373 mixed=0

n=6
  hole_only_bare             classes=6   mixed=3
  hole_only_anchored         classes=68  mixed=0
  sector_survival_bare       classes=10  mixed=2
  sector_survival_anchored   classes=118 mixed=0
  gate_priority_bare         classes=13  mixed=6
  gate_priority_anchored     classes=626 mixed=0
```

The law is stark:

```text
bare media tournaments are tiny and mixed;
observer-anchored media tournaments are pure and larger.
```

This is the event-media version of the metric-retention principle from
HYP-2018/HYP-2020: a quotient only becomes a proof object when it remembers the
part of the original problem that defines the target.

## What the Failed Bare Quotients Teach

The smallest new maps are the most seductive:

```text
hole_only_bare:
  n=4..6 class counts 3,4,6

sector_survival_bare:
  n=4..6 class counts 2,3,10

gate_priority_bare:
  n=4..6 class counts 3,6,13
```

These are much smaller than raw tournament class spaces.  But they are mixed.
They can see "a hole exists" or "which gate fires next" and still not know
whether the hole is the observer's hole.  In hindsight this should be obvious,
but it is exactly the sort of obvious thing a pretty quotient hides.

So:

```text
restriction without target anchoring is compression, not progress.
```

The observer is not decoration.  The observer is the problem.

## The New Node Types

### 1. Holes As Nodes

This is the strongest conceptual inversion.  Since `n-1` runners occupy `n`
sectors, at least one hole always exists.  Instead of tracking runners, track
vacancies.

LRC becomes:

```text
the hole process exhibits two adjacent holes at the observer clasp.
```

Edges compare survival:

```text
hole A -> hole B  iff  A survives longer before a runner enters it.
```

The bare hole tournament is tiny but mixed because it forgets where the holes
are.  The anchored hole tournament is pure because it remembers distance to the
observer.  This suggests an exclusion-process version of LRC:

```text
n-1 indistinguishable particles on n cyclic cells;
speed arithmetic controls which adjacent swap fires next;
LRC asks whether the vacancy pair visits the observer clasp.
```

This is not the same as HYP-2022.  HYP-2022 puts sectors at nodes.  Here the
nodes are *empty sectors*, so the node set itself changes with time.

### 2. Boundary Gates As Nodes

A runner crossing a boundary gate is the atomic event in the sector picture.
For a fixed open cell, every gate has a next failure time:

```text
tau_g(t) = next time some runner crosses boundary g/n.
```

Orient:

```text
g -> h  iff  tau_g < tau_h
```

This is a kinetic tournament.  Edges change when the certificate "gate g fires
before gate h" fails.  That language is standard in kinetic data structures:
objects move continuously, certificates have failure times, and an event queue
updates the structure.  The LRC translation is:

```text
the two observer gates must jointly lose incoming pressure long enough that
their adjacent cells are both empty.
```

The bare gate-priority tournament is again tiny but mixed.  The anchored gate
tournament is pure but pays the largest label tax of the probe.  That is not a
bug: gate priority is future-facing, while loneliness is a present occupancy
condition.  To make a future-priority tournament certify present loneliness,
one must color in the present boundary state.

### 3. Crossing Events As Nodes

The next abstraction is not implemented here but is now well-defined:

```text
node = a crossing event (runner i crosses boundary g);
edge e -> f iff e occurs before f in the cyclic event word, or iff e repairs
more observer debt than f.
```

This turns the LRC clock into an allowable-sequence problem.  The event word is
not an arbitrary word: it is generated by straight-line rotations
`v_i t mod 1`, so it has strong stretchability/arithmetic constraints.  A
counterexample would be a closed event word avoiding the observer-empty factor.

This is the "wiring diagram" version of HYP-2019, but with the tournament on
events rather than runners.

### 4. Danger-Cover Certificates As Nodes

Let `B_i` be the forbidden time set for runner `i`:

```text
B_i = {t : ||v_i t|| < 1/n}.
```

LRC is:

```text
the circular-arc cover {B_i} does not cover the time circle.
```

A tournament can be placed on arcs or endpoints:

```text
B_i -> B_j iff B_i's right endpoint blocks B_j's next uncovered gap before
             B_j blocks B_i's.
```

Or:

```text
endpoint e -> f iff every witness interval destroyed by f was already shadowed
by e.
```

This connects directly to circular-arc graph recognition and forbidden
structure certificates.  The LRC analogue would not ask whether a graph is a
circular-arc graph; it asks whether this highly structured circular-arc family
has a non-cover certificate.  The tournament is a dominance order on certificate
failures.

### 5. Carries And Odometer Digits As Nodes

Sector occupancy can be written as an odometer/carry process.  For each speed
`v`, the map `t -> floor(n frac(vt))` jumps through sector residues.  The
boundary crossings are carries.

Possible tournament:

```text
node = carry digit or residue channel;
carry a -> carry b iff a must fire before b in every lift compatible with the
current sector state.
```

LRC becomes a carry-synchronization problem:

```text
two adjacent observer digits are zero simultaneously.
```

This is where HYP-2016/HYP-2017 prime-power degradation should re-enter.  For
`n=18`, raw parity dies; a carry tournament over the `3`-adic digit stack might
retain the missing memory.

### 6. Rotor-Router / Chip-Firing States As Nodes

The vacancy in the sector model behaves like a chip moving through a directed
cycle with arithmetic scheduling.  Rotor-router theory offers a suggestive
finite-state analogy: recurrent one-chip rotor states are unicycles, and the
operation permutes them in strongly controlled ways.

Possible tournament:

```text
node = local rotor/carry state at a boundary gate;
edge a -> b iff firing a sends the vacancy closer to the observer than firing b.
```

Then LRC asks for recurrence of a particular observer-clasp state, while a
counterexample is a recurrent orbit avoiding it.  This is not the original
rotor-router model, but it is a disciplined metaphor: vacancy motion, cyclic
local orders, and recurrence are the right primitives.

## External Hooks

Three outside languages look especially relevant.

1. Dynamic circular-arc recognition supports insert/delete operations and
   outputs minimal forbidden certificates when an update fails.  That mirrors
   the danger-cover certificate tournament.
   Source: https://arxiv.org/abs/1509.05828

2. Circular-arc graph obstruction theory uses mutually avoiding walks as
   certificates.  That is close in spirit to "a cover can only persist if every
   attempted witness interval is shadowed by another endpoint."
   Source: https://arxiv.org/abs/1408.2639

3. Kinetic data-structure work uses certificate failure times and event queues.
   Gate-priority tournaments are exactly this sort of object: edges compare
   next boundary-crossing failures.
   Source: https://dspace.library.uvic.ca/bitstreams/fc836d25-2016-4d68-b92e-003b95ef608d/download

4. Rotor-router/chip-firing theory gives a clean finite recurrence language for
   chip motion, unicycles, and cyclic local rules.  The sector-vacancy process
   is not a rotor-router system, but it may have a useful rotor-like quotient.
   Source: https://pi.math.cornell.edu/~levine/sand.pdf

## The Assumption Audit

This session breaks six assumptions.

```text
Assumption 1: tournament vertices are runners.
Break: vertices can be holes, gates, events, endpoints, carries, or proof states.

Assumption 2: a tournament represents a geometric state.
Break: it can represent future priority, survival, event ordering, or certificate
       dominance.

Assumption 3: smaller class set means better.
Break: bare hole/gate maps are tiny and mixed.  A small quotient that forgets the
       observer is not a proof quotient.

Assumption 4: the LRC target must be a source.
Break: it can be a pair of empty cells, a non-cover certificate, a vacancy visit,
       a failed gate-pressure SCC, or a recurrent orbit hitting a clasp state.

Assumption 5: boundary events are nuisance walls.
Break: boundary events are the atomic letters of a kinetic tournament.

Assumption 6: labels are cosmetic.
Break: observer anchoring is the gauge that turns compression into purity.
```

The computational moral is the strongest one:

```text
purity = compression + the right anchor.
```

For LRC, the right anchor is usually the observer clasp, sometimes plus
left/right side, adjacent empty status, endpoint ownership, or residue-channel
debt.

## Hypotheses

### H1. Observer Anchoring Is Necessary For Pure Event-Media Quotients

Any event-media tournament whose vertices are invariant under rotation of the
sector circle and which does not mark the observer clasp has mixed fibers for
all sufficiently large `n`.  Reason: rotating a good empty-pair state away from
the observer preserves the bare class but destroys loneliness.

The S539 data already shows the phenomenon for `n=4..6`.

### H2. Hole-Only Anchored Classes Are An Exclusion-Process Target

The anchored hole-only tournament should be the smallest natural event-media
quotient with pure good fibers.  Its class counts in S539 are:

```text
n=4,5,6: 14,31,68.
```

Prediction: these counts grow far slower than sector occupancy classes, because
they only see the vacancy process.  LRC becomes:

```text
every arithmetic exclusion orbit exhibits the observer-adjacent two-hole class
or its compactified boundary.
```

### H3. Gate-Priority Needs A Present-State Color

Gate priority alone is a future event queue and cannot certify present
loneliness.  The anchored gate tournament becomes pure only after coloring the
adjacent empty status.  This is the first clean separation between:

```text
future pressure memory
present target visibility.
```

### H4. A Counterexample Is A Closed Kinetic Tournament Orbit Avoiding The Clasp

In event-media language, an LRC counterexample is not primarily a bad static
tournament.  It is a closed orbit in a kinetic tournament whose edge flips are
gate-certificate failures, avoiding all observer-empty classes.

This suggests using kinetic invariants:

```text
event word forbidden factors,
gate-pressure SCCs,
monodromy of the hole process,
rotor-like recurrence classes.
```

### H5. Danger-Cover Tournaments Should Produce Non-Cover Certificates

The circular-arc cover `{B_i}` has an uncovered interval exactly when LRC holds.
An endpoint/certificate tournament should be designed so that acyclicity or a
source endpoint peels the cover.  A counterexample would require a strongly
connected endpoint-shadow tournament.

This is the cover-dual of HYP-2021's endpoint-owner protection SCC idea.

### H6. Prime-Power Carries Are The Missing Memory At `n=18`

HYP-2017 says single-level parity becomes vacuous at `n=18`.  Event-media
thinking points to the carry stack, not the residue set, as the next object:
nodes should be carry digits/gates in the `3`-adic sector odometer, with edges
ordered by forced firing precedence.

The proof target is not "some residue class wins"; it is:

```text
the carry tournament cannot keep the observer two-cell vacancy out of phase
through a full period.
```

## Conclusion

The sector mapping made the cells into nodes.  This session goes one layer
deeper: make the *changes* into the object.

The most promising outside-box stack is:

```text
real-space layer:
  anchored hole-only tournament

kinetic layer:
  gate-priority tournament with adjacent empty colors

cover layer:
  endpoint-shadow tournament on danger-cover certificates

arithmetic layer:
  carry/odometer tournament over prime-power sector digits

recurrence layer:
  rotor-like vacancy orbit class
```

Then LRC becomes a forced-hit theorem:

```text
Every primitive arithmetic event orbit in this stack hits the observer-clasp
two-hole class, or else creates an impossible closed endpoint/carry SCC.
```

This is stranger than "runners are vertices," but it is not arbitrary.  It
respects the lesson of the computation: the tournament is allowed to be simple,
but the choice of node-set must carry the exact place where the obstruction
lives.


---
source: codex-2026-06-01-S539
status: computational probe plus assumption-challenge manifesto
tags:
  - lonely-runner
  - tournament-analysis
  - section-functors
  - boundary-events
  - quotient-stack
  - assumption-challenge
---

# Section-Boundary Functors and the Vertex-Set Assumption

The user's prompt pushed on exactly the hidden habit that had been shaping
many LRC tournament sessions:

```text
Maybe the tournament vertices are evenly spaced sections of the circle.
Maybe edges change as runners enter or leave those boundaries.
Now go more outside the box than that.
```

The important move is not only "sections instead of runners."  It is the
meta-move:

```text
The vertices of an LRC tournament are a design decision.
```

Once that is made explicit, the A000568 question becomes more flexible and
more proof-shaped.  We still ask which isomorphism classes the arithmetic
clock can exhibit, but we stop assuming the class universe has to be raw
runner tournaments.  The class universe can be fixed sections, fixed
boundaries, empty sectors, crossing events, cover arcs, residue channels, or
proof obligations.

## Computation

Artifacts:

```text
04-computation/lrc_section_boundary_functors_s539.py
05-knowledge/results/lrc_section_boundary_functors_s539.out
```

Tournament Analysis declaration:

```text
vertices:
  fixed circle sections or fixed section boundaries

pairwise observable:
  occupancy count, runner-speed flux, empty/vacuum status, boundary-crossing
  flux

switch/gauge:
  lexicographic pressure comparison plus danger/occupancy colors

tie Hamiltonian path:
  increasing cyclic section index

fingerprints:
  class counts, pure-good/mixed fibers, certification rate, score histograms,
  directed triangles, SCCs, Hamiltonian paths, and a meta-tournament
```

The most direct map is `section_empty_colored_qn`: split the circle into `n`
sections, color observer-adjacent danger sections and occupancy, then orient
edges by section pressure.  In open cells, LRC is exactly exhibition of a
class where those two danger sections are empty.

The computation tested `n=4,5,6` bounded primitive speed sets.  The cleanest
signal:

```text
section_empty_colored_qn has zero mixed good/bad fibers in all three scans.
boundary_flux_q2n also has zero mixed fibers in all three scans.
void_pressure_q2n also has zero mixed fibers in all three scans.
```

The open-cell certification misses are the familiar wall-only AP behavior:
regular-polygon witnesses sit exactly on threshold walls.  That is not a
reason to abandon section vertices.  It is a reason to add the THM-383
compactified boundary layer.

## What Changed Conceptually

The old runner-vertex default says:

```text
Vertices = moving objects.
Edges = pairwise circular order.
Target = observer source or gap threshold.
```

The section functor says:

```text
Vertices = fixed cells of the measuring apparatus.
Edges = which cell has more pressure or more emptiness.
Target = danger cells empty.
```

The boundary functor says:

```text
Vertices = gates between cells.
Edges = which gate currently receives/releases more flux.
Target = the gates adjacent to the observer lose blockers.
```

The void functor says:

```text
Vertices = possible absences.
Edges = which absence is more protected.
Target = a protected absence at the observer.
```

That last line is the useful weirdness.  LRC is usually phrased as finding a
time when all runners are far from the observer.  In void language it is the
existence of an empty interval of prescribed radius around the observer.  The
empty interval is the protagonist.

## Whackier Vertex Sets

Here are the vertex choices future agents should be willing to try before
defaulting back to runner tournaments.

**Event vertices.** A vertex is a wall crossing: runner `i` enters or exits a
danger section.  Edge `e -> f` if event `e` relieves more observer pressure or
creates a larger downstream safe interval than `f`.  LRC becomes a question
about exhibiting a source chain of relief events before a blocker event closes
the corridor.

**Cover-arc vertices.** A vertex is a bad interval
`B_i={t: ||v_i t|| < 1/n}` in the time circle.  Edge `B_i -> B_j` if `B_i`
owns the earlier endpoint, has larger uncovered-overlap debt, or dominates in
the nerve.  LRC becomes non-coverage: the cover tournament must fail to close
all gaps.

**Gap vertices.** A vertex is a gap between adjacent moving points.  Edge
`g -> h` if gap `g` is longer, expanding, closer to the observer, or has
better endpoint protection.  THM-384 becomes a two-gap target inside a
gap-tournament walk.

**Residue-channel vertices.** A vertex is a CRT residue channel, p-adic depth,
or unit orbit.  Edge `a -> b` if channel `a` exports less endpoint debt than
channel `b`.  This ties HYP-2016/HYP-2017's unit-balance laws to tournament
SCC language.

**Fourier-mode vertices.** A vertex is a resonance character in the
inside-debt expansion.  Edge `m -> m'` if `m` has larger cancellation leverage
or smaller denominator penalty.  LRC becomes a forced cancellation-class
problem rather than a direct geometric covering problem.

**Nerve-facet vertices.** A vertex is an intersection pattern in the danger
cover.  Edge `sigma -> tau` if `sigma` implies, contains, or controls the
endpoint debt of `tau`.  The target is a missing top-dimensional closure.

**Proof-obligation vertices.** A vertex is a lemma-shaped obligation: sieve
row, endpoint owner, CRT channel, BLEX source class, boundary event, or cover
hole.  Edge `A -> B` if proving `A` discharges or bounds `B`.  A counterexample
would be an SCC of mutually undischarged obligations.  This is abstract, but
it may be the right tournament space for integrating product-sieve,
permutohedral, BLEX, and endpoint-pressure evidence.

## New Metrics to Measure

The S539 script used standard fingerprints, but this vertex-set reframing
suggests more abstract metrics:

```text
vertex-estrangement:
  how far the chosen vertices are from the original runners

target-locality:
  whether the LRC target is a single vertex mark, a small class family, or a
  global property

label tax:
  bits of coloring/decorations needed to make fibers pure

motion derivative:
  whether edges record current state only, first wall-crossing derivative, or
  higher transition memory

void salience:
  how much of the certificate is carried by empty regions rather than occupied
  regions

compactification debt:
  number and type of wall states needed to recover AP/regular-polygon
  extremals

proof-obligation SCC rank:
  minimal SCC size in the meta-tournament of undischarged proof clauses

functor intersection index:
  size of the class family surviving simultaneous section-empty, BLEX,
  endpoint-owner, CRT-debt, and cover-nerve constraints
```

These metrics are deliberately stranger than score sequence or `H`.  The
point is to keep asking whether the tournament is measuring the proof problem
or only measuring our habit.

## Synthesis

HYP-2018 says metric retention restricts realizable classes.  HYP-2020 says
LRC lives in a restricted quotient stack.  HYP-2021 identifies BLEX as a small
pure oriented threshold quotient.

The upstream sector-duality result HYP-2022 already says fixed sectors are a
highly restricted real-space dual of the resonance picture.  The upstream
event-media result HYP-2023 says holes, gates, certificates, and carries can
also be vertices once the observer anchor is retained.  HYP-2024 adds:

```text
Even the tournament vertex set belongs inside the quotient stack.
```

For the next agents, the standing question should be:

```text
What would the tournament vertices be if runners were forbidden?
```

That one question keeps opening new proof languages: sections, boundaries,
voids, events, covers, residues, Fourier modes, and proof obligations.  Some
will be too lossy, some too large, but the useful ones should make the LRC
target simpler than it was in runner space.

# Projection Defect as Common Residue

**Source:** opus-2026-05-29-S12  
**Computation:** `04-computation/projection_defect_bridge_s12.py`

The "projection defect" angle makes three unrelated threads rhyme:

1. H=63 unlocking through complete Omega.
2. THM-025 breaking real-rootedness.
3. Ghost cycles / HYP-408 in path homology.

Each is asking the same question in a different language:

```text
When we forget coordinates, what residue survives?
```

## Exact Kill vs Near Kill

The two n=8 H=63 classes are exact old-projection kills. Every odd directed
cycle contains one core vertex. Deleting that vertex sends the tournament to the
transitive tournament:

```text
H=63 class: 31 cycles -> delete core -> 0 cycles.
```

This is why the OCF shape is so simple:

```text
Omega(T) = K31
H(T) = 1 + 2*31 = 63.
```

The THM-025 real-rootedness counterexample is different but adjacent. It has no
cycle-family core, yet vertex 3 participates in 92 of 94 odd cycles:

```text
THM-025: 94 cycles -> delete vertex 3 -> 2 cycles.
```

Those two surviving cycles are not noise. They have `alpha=[1,2,1]`, so they are
disjoint. That tiny residue is enough to keep the unique independent triple
alive in the full tournament and to produce the non-real-root polynomial:

```text
I(Omega,x) = 1 + 94x + 10x^2 + x^3.
```

So H=63 is the exact kill case; THM-025 is the near-kill case. The failure does
not come from generic complexity. It comes from an almost-collapsed projection
with a two-cycle residue in precisely the wrong shape.

## The Same Pattern in Homology

HYP-408 and the Ghost Cycle Theorem ask the path-homology version of the same
question. Project away the through-v coordinates. If a 3-cycle generator lives
only through v, does it become a boundary?

The known equivalence says:

```text
codim(pi_old(im d4), pi_old(ker d3)) = 1
iff every through-v-only cycle is a boundary.
```

That is projection defect again. The quotient sees either a genuine old
component or a ghost that disappears into boundaries.

## The Even-Graph Version

For odd n, the tournament-to-even-graph map is another projection:

```text
T_cycle = (I + L(K_n))*T mod 2.
```

For even n this projection is ambiguous because cut and cycle spaces intersect.
This matters: the H=63 classes live at n=8, exactly where the clean odd-n
cycle-space projection is unavailable. The "ambiguity at even n" may not be
incidental; it may be part of why these exact projection-kill objects first
appear there.

## Same Shadow, Different Fibers

The Paley and interval tournaments on 7 vertices give a second outside-facing
version of the same phenomenon. They have the same regular score sequence and
the same number of odd-cycle supports, but the fibers over those supports are
different:

```text
Paley T7:     80 cycles, 36 supports, support_excess 44, max_mult 24
Interval T7:  59 cycles, 36 supports, support_excess 23, max_mult 17
```

So the visible hypergraph shadow is identical, while the multiplicity upstairs
is not. This is not only a tournament fact. In hypergraph language it is a
shadow/fiber defect; in statistical mechanics it changes the hard-core
partition weight without changing the visible support set.

The even-graph projection separates them just as cleanly:

```text
transitive T7 -> 21 edges, degree sequence 6^7
Paley T7      -> 14 edges, degree sequence 4^7
Interval T7   ->  7 edges, degree sequence 2^7
```

That looks like a codeword-weight ladder in the cycle space. The Paley/interval
contrast may therefore belong as much to coding theory and shadow multiplicity
as to tournament extremal theory.

## New Working Principle

The project has many invariants that are really shadows:

- H is the fugacity-2 projection of the full independence polynomial.
- Scores are the cut-space projection of orientation data.
- Even graphs are the cycle-space projection.
- Deletion profiles are old-coordinate projections.
- Path-homology relative groups measure what survives quotienting by a
  subcomplex.

The useful invariant may be the defect, not the projection itself.

Concrete features to track:

- `max_v rho_v`, the maximum fraction of odd cycles killed by deleting v.
- kill vertices with `Omega(T-v)` empty.
- support multiplicity defect: directed cycles minus support sets.
- max support multiplicity and number of repeated supports.
- odd-n even-graph projection degree sequence.
- whether a near-kill residue has `alpha=[1,2,1]` or higher disjoint structure.

Engineering translation: these are cheap Tournament TDA features. They can
distinguish localized inconsistency ("one pivot explains almost all cycles")
from dispersed inconsistency ("many independent cycle packets") before running
expensive homology or full Omega computations.

## Next Experiment

Scan known and newly sampled non-real-root examples at n=9 and n=10. If they
cluster at high `max_v rho_v`, then real-root failure is a projection-residue
phenomenon rather than a generic dense-Omega phenomenon. Then compare the same
statistic against beta_3 / beta_4 path-homology anomalies.

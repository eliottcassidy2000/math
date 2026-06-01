---
source: codex-2026-06-01-S510
status: exploratory synthesis with exact small-clock computation
tags:
  - lonely-runner
  - tournament-analysis
  - A000568
  - tournament-clock
  - isomorphism-classes
  - marked-fibers
---

# LRC and A000568: The Marked Chamber Walk

The user's hypothesis is right in a precise but slightly twisted way.  A000568
is the chamber count for unmarked tournament shapes.  A runner system, under
the half-turn phase comparator, does not roam freely through those chambers; it
draws a closed arithmetic walk through a very small image.  LRC then asks a
marked question on top of that walk: does the observer vertex ever sit in a
`1/N`-safe placement?

So the analogy is:

```text
A000568(N)              ambient unmarked tournament chamber count
tournament clock        arithmetic closed walk through those chambers
observer-fixed quotient marked fiber over the chamber walk
LRC witness             hit of a marked safe target
counterexample          closed arithmetic walk avoiding all marked targets
```

## What S510 measured

`lrc_a000568_marked_chamber_s510.py` first recomputes A000568 by the odd-cycle
Burnside formula.  The large values matter because they set the ambient chamber
scale:

```text
A000568(14) = 28401423719122304
A000568(18) = 1783398846284777975419600287232
```

Then the script canonicalizes small runner clocks exactly.  Initial consecutive
rows are tiny chamber walks:

```text
N5 initial: 12 cells, 4 of 12 unmarked classes
N6 initial: 20 cells, 4 of 56 unmarked classes
N7 initial: 24 cells, 7 of 456 unmarked classes
```

Prime-like or lopsided rows add cells and marked placements, but still only
touch a small circular subfamily:

```text
N7 prime-like: 108 cells, 10 of 456 unmarked classes, 55 marked classes
```

That is the first strong point: LRC does not need all of `G_N`; it needs the
circular arithmetic image of `G_N`.

## Why the marker matters

The unmarked class is not enough.  In `N=5 initial`, the safe witness samples
all land in the regular `H=15` chamber, but that same unmarked chamber appears
at unsafe cell midpoints too.  In `N=6 initial`, the safe samples land in the
near-regular `H=41` score family; the same unmarked score shape also has unsafe
observer placements.

That is the exact reason recent LRC work kept adding observer scores,
two-neighbor brackets, endpoint debt, and pressure labels.  The unmarked
tournament shape says "this configuration is spread."  The marked fiber says
"the stationary vertex is actually safe."

## Large-row shadow

Full isomorphism reduction at `N=14` or `N=18` is not currently the right
computation, but the coarse fingerprint shadow already shows the scale
separation:

```text
n14 initial:     116 cells, 26 coarse fingerprints
n14 row-parent: 1772 cells, 140 fingerprints
n14 gate:       3856 cells, 156 fingerprints
n18 initial:     192 cells, 44 fingerprints
n18 row-parent: 4008 cells, 304 fingerprints
n18 gate:       8512 cells, 344 fingerprints
```

The fingerprints are not classes.  They are deliberately coarse: score
histogram, directed triangles, observer score, largest SCC, source count, and
sink count.  But the number of clock cells itself is an upper bound on visited
phase classes, so even the hard rows occupy a microscopic part of the A000568
ambient space.

## The proof reframe

The deepest version is not "prove LRC by counting tournaments."  It is:

1. characterize the circular arithmetic image inside the A000568 chamber graph;
2. lift it to the observer-marked quotient;
3. mark the safe target set using endpoint-clock data;
4. prove every relevant arithmetic closed walk hits that target or yields a
   pressure certificate.

This unifies the tournament clock, H as a chamber coordinate, and the gauge
bundle.  `H` remains a good coordinate inside the circular image, but endpoint
debt and observer marking decide the actual LRC target.

## Convoluted paths worth keeping

**Burnside resonance.** A000568's Burnside formula only sums odd cycle types:
even permutation cycles cannot fix an oriented complete graph.  Runner-clock
walls collapse when speed differences share denominators.  Both stories are
quotients controlled by parity and gcd data.

**q-deformed pressure.** The repo's `A(n,q)` thread suggests a pressure-weighted
chamber count.  `q=2` counts unmarked tournaments; extra weights could remember
endpoint owners, protectors, or pressure layers.

**Royle-even shadow.** If tournament classes and certain even-graph classes are
two faces of the same quotient, then endpoint-safe masks may be easier to read
as an even-graph shadow of the marked tournament clock.

**Sparse-image proof.** A proof for fixed `N=14` or `N=18` should not classify
A000568(N).  It should classify the much smaller arithmetic image, then prove
target hitting in its marked fiber.

**Disproof signature.** A plausible counterexample should now be described as a
closed arithmetic loop in the marked chamber fiber with no safe target hit and
with a nonempty endpoint-pressure core.  Anything less is only a scalar near
miss.

## Artifacts

```text
04-computation/lrc_a000568_marked_chamber_s510.py
05-knowledge/results/lrc_a000568_marked_chamber_s510.out
05-knowledge/hypotheses/HYP-1979-lrc-a000568-sparse-marked-chamber-walk.md
```

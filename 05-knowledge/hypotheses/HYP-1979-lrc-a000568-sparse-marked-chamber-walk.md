---
id: HYP-1979
status: OPEN
source: codex-2026-06-01-S510
related:
  - HYP-1951
  - HYP-1967
  - HYP-1970
  - HYP-1971
  - HYP-1972
  - HYP-1974
  - HYP-1975
  - HYP-1976
  - HYP-1977
  - THM-373
  - THM-374
---

# HYP-1979: LRC clock images are sparse marked walks inside A000568

## Statement

Refining HYP-1977, the Lonely Runner Conjecture should track how arithmetic
runner clocks occupy a sparse subset of the tournament isomorphism-class space
counted by A000568.  The stationary runner is a marked vertex, and the
endpoint/gap fiber over that marked walk is proof-relevant data.

For total runner count `N`, the half-turn tournament clock maps a speed set to a
closed walk through unlabelled tournament chambers:

```text
speed row V  ->  phase-clock movie T_V(t)  ->  walk in Tour_N / S_N
```

A000568(N) counts the ambient unmarked chamber set.  LRC does not ask only
which unmarked chamber is visited; it asks whether the observer-marked lift of a
visited chamber hits a `1/N`-safe target.  Thus the natural object is a marked
fiber over A000568:

```text
observer-marked clock:  T_V(t) with vertex 0 fixed
target set:             marked chambers containing a safe observer placement
counterexample shape:   arithmetic closed walk avoiding the target set
```

The analogy to A000568 is therefore not "LRC equals an unmarked tournament
count."  It is that LRC lives on a sparse arithmetic image inside the
isomorphism-class metagraph whose chamber count is A000568, while HYP-1977
explains why the endpoint/gap fiber cannot be projected away.

## Evidence

`04-computation/lrc_a000568_marked_chamber_s510.py` computes A000568 by the
odd-cycle Burnside formula, then canonicalizes exact small runner clocks under
all relabelings and under observer-fixing relabelings.

Small exact clocks are extremely sparse in the ambient chamber space:

```text
N5 initial 0..4:      12 cells, 4/12 unmarked classes, 9 marked classes
N6 initial 0..5:      20 cells, 4/56 unmarked classes, 12 marked classes
N7 initial 0..6:      24 cells, 7/456 unmarked classes, 17 marked classes
N7 prime-like:       108 cells, 10/456 unmarked classes, 55 marked classes
```

The marked fiber is not cosmetic.  For `N=5 initial 0..4`, the safe witness
samples all land in the regular `H=15` unmarked chamber, but that chamber also
appears at unsafe cell midpoints.  For `N=6 initial 0..5`, the safe witness
samples land in one of the near-regular `H=41` marked placements while the same
unmarked score class has unsafe placements.  This is exactly the LRC issue:
unmarked tournament shape sees spread, but the stationary vertex label decides
loneliness.

For large LRC rows the exact isomorphism reduction is too large, but coarse
fingerprints show the same scale separation:

```text
n14 initial:     116 clock cells, 26 fingerprints, A000568(14)=28401423719122304
n14 row-parent: 1772 clock cells, 140 fingerprints
n14 gate:       3856 clock cells, 156 fingerprints
n18 initial:     192 clock cells, 44 fingerprints, A000568(18)=1783398846284777975419600287232
n18 row-parent: 4008 clock cells, 304 fingerprints
n18 gate:       8512 clock cells, 344 fingerprints
```

Even allowing that fingerprints are coarse, the cell count itself upper-bounds
the number of visited phase-clock chambers.  The LRC arithmetic image is tiny
relative to A000568(N).

## Interpretation

This hypothesis sharpens HYP-1951 and supplies a sparse-image count for
HYP-1977.  The tournament clock is a walk in `G_N`, but LRC is not solved by
knowing the unmarked walk alone.  The proof object should be:

1. the circular arithmetic image inside `Tour_N / S_N`;
2. the observer-marked fiber over that image;
3. endpoint-clock target labels on the marked fiber;
4. pressure/gauge data from HYP-1975 decorating the same walk.

A counterexample would be a closed arithmetic walk that avoids every marked
target while still satisfying the endpoint-protection constraints.  That makes
LRC look less like a scalar Diophantine inequality and more like a sparse
target-avoidance theorem in an A000568 chamber complex.

## Creative paths

1. **Burnside resonance:** A000568 ignores even permutation cycles because those
   symmetries cannot fix a tournament.  Runner clocks collapse wall events when
   speed differences share denominators.  Both quotient stories are governed by
   cycle parity and gcd data.
2. **q-deformed pressure:** The repo's `A(n,q)` viewpoint suggests weighting
   tournament chambers by endpoint pressure: `q=2` is unmarked counting, while
   larger or multi-parameter weights remember labelled endpoint debt.
3. **Royle-even shadow:** Since even-graph models are equinumerous with
   tournament classes in this repo's metagraph thread, safe endpoint masks may
   be the even-graph shadow of the marked tournament clock.
4. **Sparse-image proof:** A proof need not classify all A000568(N) chambers.
   It may suffice to prove that the circular arithmetic image always intersects
   the marked safe target.
5. **Disproof signature:** A serious counterexample candidate should provide a
   long arithmetic loop in the marked fiber with no safe samples, not merely a
   small scalar gap in one time cell.

## Predictions

1. Speed rows with the same half-turn clock walk can differ in LRC behavior only
   through marked-fiber or endpoint-clock labels.
2. Initial segments should continue to be minimal or near-minimal walks in the
   A000568 chamber space, matching their extremal LRC role.
3. The next useful computation is a labelled endpoint-runner switchboard on top
   of exact small marked fibers, then a fingerprint surrogate for `N=14` and
   `N=18`.
4. Any future proof of fixed `N=14` or `N=18` cases should be expressible as
   "every arithmetic clock walk in the relevant sparse image hits the marked
   target or has an endpoint-pressure certificate."

## Sources

- `04-computation/lrc_a000568_marked_chamber_s510.py`
- `05-knowledge/results/lrc_a000568_marked_chamber_s510.out`
- `07-reflections/lrc-a000568-marked-chamber-walk-s510.md`

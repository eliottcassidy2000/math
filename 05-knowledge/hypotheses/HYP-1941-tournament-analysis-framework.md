---
id: HYP-1941
status: OPEN
source: codex-2026-06-01-S471
related:
  - THM-372
  - THM-370
  - HYP-1900
  - HYP-1903
  - HYP-1904
  - HYP-1921
  - HYP-1930
  - HYP-1931
  - HYP-1940
---

# HYP-1941: Tournament Analysis is the repo's pairwise-data functor

## Statement

For many problems in the repo, the productive abstraction is not the raw
metric, scalar invariant, or weighted graph.  It is a tournament-valued
functor:

```text
objects + pairwise observable + switch functional + tie Hamiltonian path
  -> tournament trajectory.
```

The pairwise observable may be discrete, such as basketball pass counts, or
continuous, such as runner distances on a circle.  The switch functional may
use majority, threshold crossing, metric derivative, isolation, simplex
stress, or two-neighbor deletion relief.  Ties are resolved by a fixed
Hamiltonian path, exactly as in the repo's base-path/staircase model.

The hypothesis is that useful problem analogies preserve the switch-and-residue
structure of this functor, not merely the visual metric.

S500 canonized the finite construction layer as THM-372: a switch value for
each unordered pair plus a fixed tie Hamiltonian path always produces a unique
tournament, and finite-wall switchboards produce piecewise-constant tournament
movies.

## Evidence

`tournament_analysis_framework_s471.py` implements the schema on two families.

Basketball pass matrices:

```text
point-hub:     H=1, transitive
motion-ties:   H=3, one directed 3-cycle
inverted-big:  H=1, transitive with different top role
```

The role order `1->2->3->4->5` acts as a tie Hamiltonian path.  This is the
same abstraction as the fixed base path in the tiling model.

Runner configurations:

```text
initial n=14, t=1/14:
  semicircle      H=24104937, 3cyc=112
  chord-opening   H=24104937, 3cyc=112
  fourier-opening H=24577317, 3cyc=112
  centrality-type rules collapse to H=1

n=14 seven-ladder boundary:
  semicircle      H=19622601, 3cyc=108
  chord-opening   H=21934885, 3cyc=110
  fourier-opening H=20727933, 3cyc=108
  pressure        H=5,        3cyc=2
```

The initial five-runner chord-threshold cuboid trajectory sampled at `t=q/60`
visits:

```text
12 distinct tournaments,
H spectrum {1,5,9,11,13,15}.
```

This verifies that a symmetric continuous metric can be made into a genuine
edge-bit tournament trajectory by thresholding against a base Hamiltonian path.

## Predictions

1. Every major repo thread should be expressible as a choice of pairwise
   observable and switch functional, followed by tournament fingerprints.
2. Scalar methods fail when their switch functional collapses to a transitive
   ranking and discards incidence residue.
3. LRC proof searches should record multiple tournamentizations of the same
   runner movie: semicircle, chord-threshold cuboid, opening/closing,
   pressure, endpoint handoff, and simplex stress.
4. The most useful analogies will be those where the same switch functional
   has a meaningful interpretation in both domains.
5. A future meta-index should catalogue switch functionals by whether they
   preserve cycles, private leaves, endpoint labels, Zeckendorf/path masks, or
   cuboid edge flips.

## Sources

- `04-computation/tournament_analysis_framework_s471.py`
- `05-knowledge/results/tournament_analysis_framework_s471.out`
- `07-reflections/tournament-analysis-framework-s471.md`
- `01-canon/definitions.md`
- THM-372
- THM-370
- HYP-1904
- HYP-1921
- HYP-1930

---
id: THM-826
title: The B3 face cover is a weighted tournament deletion-deck calculus
status: RESERVED / PROOF AND EXACT REPLAY IN PROGRESS
source: codex-2026-07-14-S14
depends_on: [THM-442, THM-550, THM-553, THM-801]
related: [THM-796, THM-809, THM-814, THM-818, HYP-2685, HYP-3234]
---

# THM-826 - the B3 face cover is a weighted tournament deletion-deck calculus

Namespace reservation after rebasing onto live `origin/main` at `0bb0e98a2`.

The proved algebraic core will identify the three B3 face roles with the three
relative positions of a deleted fixed-path vertex with respect to a tournament
chord.  For a tiling bit `t_(a,b)`, the three multiplicities are provisionally

```text
n-a,  a-b-1,  b-1,
```

and sum to `n-2`, the number of vertex deletions under which the chord
survives.  The resulting role masses satisfy a mirror-even internal term and a
mirror-odd boundary current.  This retains the geometric distinction that is
lost when the odd seven-term word
`A+B-C+D-E-F+G` is scalarized.

Still missing at reservation time: a convention-checked coordinate proof, the
exact all-size B3 stratum enumerator, the blue/black one-flip classification,
the Tournament Analysis fingerprints, exhaustive replay, and the final
preservation/destroyed-information boundary.

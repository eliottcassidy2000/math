---
id: HYP-1986
status: OPEN
source: codex-2026-06-01-S516
related:
  - HYP-1990
  - HYP-1989
  - HYP-1985
  - HYP-1982
  - HYP-1981
  - HYP-1979
  - HYP-1977
  - THM-384
  - THM-383
  - THM-382
  - THM-381
  - THM-380
  - THM-369
---

# HYP-1986: LRC proof by compactified source-gap forcing

## Statement

After THM-384, the Lonely Runner Conjecture can be reduced to a forced-visit
problem in a compactified threshold-gap walk.

For a fixed primitive speed system with total denominator `n`, the runner clock
traces an exact cyclic walk through wall-compactified circular states.  Mark the
stationary observer and color the two observer-adjacent gaps by the closed
threshold `>= 1/n`.  The LRC target is the single source-gap fiber

```text
(long, long).
```

A counterexample is therefore an admissible compactified arithmetic walk that
avoids `(long,long)` at every open cell and every equality wall.  The
hypothesis is that no such walk can also satisfy the necessary arithmetic and
pressure constraints from THM-369 and THM-380.

More explicitly: every source-gap-avoiding compactified walk either fails the
denominator sieve, has an empty endpoint core, or has a peelable
owner-compatible pressure graph.  Hence it cannot be a true LRC counterexample.

## Evidence

S516 adds the proof-level local criterion THM-384 and checks it over exact
compactified bounded samples:

```text
total n  max speed  systems  states  mismatches  source-hit
      3         16       79    9475           0          79/79
      4         10      109   12825           0         109/109
      5          9      125   16921           0         125/125
      6          8       56    8568           0          56/56
```

In all sampled states the only good adjacent threshold fiber is `(1,1)`.
The bounded audit also records wall-only witnesses, confirming THM-383's
warning that equality walls must remain in the proof object.

On the selected hard rows from S514, the source-gap lens separates initial
segments from dyadic ladder rows:

```text
n14 initial: open corridors 0, wall-only
n16 initial: open corridors 0, wall-only
n18 initial: open corridors 0, wall-only

n14 row-parent/gate/double-gate: 24/48/96 open corridors
n16 row-parent/gate:             42/84 open corridors
n18 row-parent/gate/double-gate: 64/128/256 open corridors
```

Thus the hard rows are not absent from the source-gap fiber; their witnesses
are short positive corridors whose widths shrink predictably under dyadic
translation.

## Proof Route

1. Use THM-384 to replace the local LRC predicate by the source-gap fiber
   `(long,long)`.
2. Use THM-383's wall compactification so boundary-only witnesses and
   equality-tie states cannot be lost.
3. If the speed system is not denominator-sieve complete, invoke THM-369 to
   obtain a rational lonely time.
4. If a source-gap-avoiding walk remains, construct the terminal endpoint core
   of the closed forbidden arcs.
5. Refine the S514 pressure coordinate from coarse row-time SCCs to the
   owner-compatible endpoint-core pressure graph of THM-380.
6. Prove that the compactified source-gap avoidance relation makes this
   labelled pressure graph peelable.  This contradicts the THM-380 pressure
   cycle required by a pressure-realized nonempty endpoint core.

## Predictions

1. Raw A000568 or half-turn phase classes will continue to be too coarse; the
   useful quotient is the compactified source-gap walk with observer mark.
2. Initial consecutive speed systems are boundary sections, not open-cell
   witnesses.
3. First-even hard ladders should keep producing short positive corridors, with
   widths scaled by the dyadic gate factor.
4. A serious counterexample candidate must be sieve-complete, source-gap
   avoiding, endpoint-core nonempty, and owner-pressure cyclic at once.
5. The pressure obstruction will disappear only after labels are forgotten;
   the labelled owner-compatible core should be peelable in all true arithmetic
   walks.

## Verification

Stored script and output:

```text
04-computation/lrc_gap_threshold_proof_route_s516.py
05-knowledge/results/lrc_gap_threshold_proof_route_s516.out
```

The proof-route tournament in the script ranks `threshold_gap_source`,
`compactified_gap_walk`, and `endpoint_pressure_core` as the top three route
languages.  The tournament is transitive (`H=1`, `c3=0`), which reflects a
strict narrowing rather than a tradeoff cycle.

S518 adds the THM-386 two-gap race refinement: in the source-gap walk, LL can
only be entered from LS and exited to SL.  Thus the compactified forcing route
can be sharpened to proving that at least one wrap-around reset wins the
right-gap race before the left gap decays below the threshold.

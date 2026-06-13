---
id: HYP-1982
status: OPEN
source: codex-2026-06-01-S512
related:
  - HYP-1978
  - HYP-1977
  - HYP-1976
  - HYP-1975
  - HYP-1974
  - HYP-1972
  - HYP-1971
  - HYP-1951
  - HYP-1824
  - THM-373
  - THM-374
  - THM-382
  - THM-383
  - THM-384
  - HYP-1986
---

# HYP-1982: LRC is a forced visit in a threshold-decorated tournament-class fiber

## Statement

The Lonely Runner Conjecture should be read as a constrained walk problem over
the A000568 tournament-isomorphism base, but not over the raw base alone.

The likely object is a fiber:

```text
ordinary tournament iso class
+ circular/tie-wall stratum
+ observer mark
+ 1/n threshold colors on gaps or pair-cells
```

For a fixed speed set, the runner clock traces a finite closed walk in this
decorated class space.  A lonely time is a visit to a good-only fiber.  A proof
of LRC at denominator `n` would show that every admissible runner-clock walk is
forced to hit one of those good-only fibers.

This preserves the user's A000568 hypothesis while correcting the minimal
object: raw unmarked classes are too coarse, but threshold-decorated fibers over
those classes begin to behave like a proof language.

## Evidence

`lrc_iso_class_constraint_s512.py` enumerates A000568 by Burnside and direct
canonicalization through `n<=5`, then tests several runner-to-tournament lifts
on exact wall/cell samples for total `n=3,4`.

S513 canonized the bounded parts of this evidence as THM-382 and THM-383.  The
global forced-visit statement remains open.

S516 sharpened the local gap side as THM-384: once the observer is marked,
LRC-good status is exactly the closed `1/n` threshold color of the two
observer-adjacent gaps.  Thus the remaining open problem is the HYP-1986
compactified source-gap forced-visit theorem.

The open half-turn clock is genuinely class-restricted:

```text
total n=3: open phase classes 2 of A000568(3)=2
total n=4: open phase classes 2 of A000568(4)=4
total n=5: open phase classes 4 of A000568(5)=12
```

But equality walls matter for LRC.  With the fixed tie Hamiltonian path, the
wall-compactified image expands:

```text
total n=4: wall-compactified phase classes 4 of 4
total n=5: wall-compactified phase classes 11 of 12
```

Thus tight LRC witnesses live in a boundary compactification of the circular
menu, not just in the open circular menu.

For total `n=3`, primitive speed sets through `max_speed=16` all have sampled
exact witnesses (`79/79`).  Raw class lifts do not certify them:

```text
phase_half:                     2 classes, 2 mixed, certifies 0/79
phase_marked_observer:          4 classes, 3 mixed, certifies 0/79
gap_rank_marked_adjacent:       3 classes, 2 mixed, certifies 0/79
pair_deficit_origin_marked:     3 classes, 2 mixed, certifies 0/79
```

Threshold-decorated fibers do certify every sampled speed set:

```text
gap_threshold_fiber:            8 classes, 2 good-only, certifies 79/79
pair_deficit_threshold_fiber:   9 classes, 2 good-only, certifies 79/79
```

For total `n=4`, primitive speed sets through `max_speed=10` all have sampled
exact witnesses (`109/109`).  The same separation occurs:

```text
phase_half:                     4 classes, 4 mixed, certifies 0/109
phase_marked_observer:         11 classes, 10 mixed, certifies 0/109
gap_rank_marked_adjacent:       6 classes, 4 mixed, certifies 0/109
pair_deficit_origin_marked:    17 classes, 4 mixed, certifies 0/109
gap_threshold_fiber:           20 classes, 5 good-only, certifies 109/109
pair_deficit_threshold_fiber:  50 classes, 4 good-only, certifies 109/109
```

This is strong negative and positive evidence at once: the A000568 analogy is
real, but it needs threshold color in the fiber.

## Predictions

1. Raw unmarked tournament class will remain too coarse for LRC at larger `n`;
   good and bad clock states will share ordinary A000568 classes.
2. The right `n=14` target is a compressed decorated fiber, not enumeration of
   all A000568(14) classes.  Candidate fibers are gap-threshold classes,
   pair-cell danger-deficit threshold classes, and endpoint-pressure labelled
   classes.
3. The hard `n=14` quotient-ladder rows should be forced into a small region of
   the decorated fiber; a proof would show that every closed walk in that
   region hits a good-only threshold fiber.
4. Pressure DAG failures should correspond to bad-only or mixed decorated
   fibers.  A labelled pressure SCC would be the obstruction to the forced-visit
   theorem.

## Sources

- `04-computation/lrc_iso_class_constraint_s512.py`
- `05-knowledge/results/lrc_iso_class_constraint_s512.out`
- `07-reflections/lrc-as-threshold-decorated-tournament-class-fiber-s512.md`
- HYP-1951
- HYP-1976
- HYP-1975
- HYP-1972

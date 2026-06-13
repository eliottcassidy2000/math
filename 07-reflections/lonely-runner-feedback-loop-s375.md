---
source: codex-2026-05-31-S375
status: research sprint
tags:
  - lonely-runner
  - fourteen-runners
  - fifteen-runners
  - disproof-search
  - scalar-puncture
---

# Lonely Runner Feedback Loop S375

This session enforced the loop:

```text
14-runner dead end -> 15-runner new idea -> disproof construction pressure -> repeat
```

The best new object is not another near-blocker.  It is a repair deficit.

S372 showed that the best `n=14` non-scalar near-blocker is a scalar ramp with
one coordinate punctured.  S375 asked what happens if we try to repair the cells
opened by that puncture.

For the zero-ramp half-turn

```text
(0,0,0,0,0,7,0,0,0,0,0,0,0)
```

there are `56` exposed cells.  Reverting coordinate `6` closes them, of course.
But the best non-reverting move is striking:

```text
coord 3: 0 -> 7
gain on old misses: 56
new total misses: 308
```

So a move can cover the old exposed cells and still make the obstruction much
worse.  The `n=15` analogue has the same taste.  Starting from the best
`5`-layer puncture at coordinate `6`, the best non-reverting repairs cover only
`60` of the `120` old misses and leave at least `340` total misses.

The torsion-shell scans make the same point at larger support:

```text
n=14 half-turn shell support best misses: 56, 112, 126, 182, 168
n=15 {5,10}-shell support best misses:   120, 220, 280, 290
```

Nothing in these local shells points toward a hidden full blocker.  The good
news is that they look structured enough to prove.

The disproof lane stayed leaky too.  Local gated speed-set search found no open
cover.  The best `n=14` examples kept positive gap ratio `0.037879`; the best
`n=15` example kept positive gap ratio `0.030303`; every listed endpoint peel
ended with empty core.

The next proof target should be an exposed-cell repair theorem:

```text
Any non-reverting repair of a scalar puncture creates more exposed cells than it closes.
```

If that can be proved first for the eight `n=14` alpha stencils, it becomes a
branch-and-bound lower bound for the normalized quotient search.  Together with
THM-363/THM-364 and endpoint descent, it gives a clean three-part route:

```text
scalar spine -> puncture repair deficit -> far-field endpoint descent
```

## Artifacts

- `04-computation/lonely_runner_feedback_loop_s375.py`
- `05-knowledge/results/lonely_runner_feedback_loop_s375.out`
- `05-knowledge/hypotheses/HYP-1830-lrc-exposed-cell-repair-deficit.md`

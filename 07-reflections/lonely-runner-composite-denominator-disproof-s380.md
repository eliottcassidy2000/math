# Lonely Runner 14 Composite-Denominator Disproof Attempt

Session: `codex-2026-05-31-S380`

The disproof attempt started from the observation in S376/S377 that `n=14` is
not just another composite denominator.  Its scalar-puncture moat has
curvature concentrated at the `2*7` and `2*7*13` leak layers, while the
speed-side seven-ladder has a very small gap but many exposed endpoints.

The new script is
`04-computation/lonely_runner_14_composite_denominator_disproof_s380.py`, with
stored output in
`05-knowledge/results/lonely_runner_14_composite_denominator_disproof_s380.out`.

## Disproof Target

For `k=13`, threshold `1/14`, a counterexample must be an exact open forbidden
cover of the circle.  By THM-360 it must protect the unit endpoints `a/14`,
so the search was biased toward speed sets with `14`-gates.

The new idea was to go beyond the unit layer and explicitly target the first
composite leak orbit:

```text
9/98, 29/182, 15/182,
```

together with the largest exact gaps of the seven-ladder and single-14-gate
seeds.  This produced a target layer of `187` rational points, with
denominators concentrated at `14,49,56,91,98,182,392,1078,1176,1274` and a few
larger descendants.

## Best New Near-Disproof

The most interesting construction was the target-grid greedy set forced to
contain `(1,14,98)`:

```text
(1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182)
```

This is essentially the `14`-multiple ladder, with the `84` gate skipped.  It
covered all `187` target-layer points on the designed composite grid.  The
exact audit was:

```text
classification = positive_gap
forbidden length = 142/143
max_gap = 5/25872
gap / threshold = 0.002706
unprotected endpoints = 168
first unprotected = 15/196
```

This improves the seven-ladder's gap ratio:

```text
seven-ladder gap/th = 0.005411
14-ladder gap/th   = 0.002706
```

but it does not move toward an open cover in the endpoint sense.  It pays for
smaller open gaps by increasing endpoint leakage.

## Dead End

No open-cover candidate appeared:

```text
open_cover_candidates = 0
boundary_only_candidates = 1
positive_gap_candidates = 89
```

The initial segment remains the only boundary-only object among the audited
best candidates.  Every composite-denominator construction that protected the
unit and `98/182` target layers still exposed a positive gap.

The best positive-gap candidates kept their maximal gap endpoints on the same
composite-denominator descendants:

```text
2156, 2352, 2548
1078, 1176, 1274
490, 686
```

So the anomaly is real, but it behaves like an endpoint-transfer obstruction
rather than a disproof architecture.

## Synthesis

The attempt suggests a sharper negative rule for future disproof searches:

```text
protecting more of the composite target layer lowers open-gap size,
but it increases endpoint-transfer debt.
```

This became HYP-1837.  The next disproof route should stop choosing speeds
first.  It should prescribe an endpoint-protection graph on the unit and
`98/182` leak orbit, solve the finite integer protection constraints, and only
then ask whether any 13-speed realization exists.

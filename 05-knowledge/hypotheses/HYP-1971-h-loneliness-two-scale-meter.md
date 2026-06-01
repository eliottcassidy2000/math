---
id: HYP-1971
status: OPEN
source: codex-2026-06-01-S506
related:
  - HYP-1951
  - HYP-1941
  - HYP-1950
  - HYP-1960
  - THM-370
  - THM-374
---

# HYP-1971: H is a two-scale loneliness meter for runner-clock tournaments

## Statement

For runner-clock tournaments, `H(T)` should be read as a two-scale loneliness
meter, not as a single monotone LRC invariant.

1. Low `H`, and especially `H=1`, measures unanchored bunching.  By THM-374,
   the half-turn circular tournament is transitive exactly when all points lie
   in an open semicircle, equivalently when the maximum circular gap is greater
   than `1/2`.
2. High `H` often measures LRC-style local separation.  Near-regular circular
   configurations have many vertices whose two adjacent circular gaps are both
   at least the LRC threshold `1/n`; the initial segment at `t=1/n` has every
   vertex lonely in this local sense.
3. The Lonely Runner Conjecture needs the marked structure behind the scalar:
   circular order, gap vector, stationary vertex, safe-gap mask, and pressure
   or deletion data.  `H` is a useful scalar shadow of that marked tournament
   movie, but it is not the invariant by itself.

## Evidence

`h_loneliness_meter_s506.py` audits exact clock cells and selected LRC witness
times.  In small circular clocks, `H` is strongly negatively correlated with
maximum circular gap, matching the low-H bunching reading:

- initial `n=5`: Pearson `(H,max_gap) = -0.946`
- primes `n=5`: Pearson `(H,max_gap) = -0.902`
- initial `n=6`: Pearson `(H,max_gap) = -0.912`
- spread `n=7`: Pearson `(H,max_gap) = -0.716`

The same output also shows why the LRC reading is not the low-H endpoint.  For
initial total `n=5..9`, the best marked stationary witness is `t=1/n`; every
vertex is locally lonely, but `H` is high:

```text
n=5: H=15
n=6: H=41
n=7: H=175
n=8: H=629
n=9: H=3267
```

For selected larger rows, the contrast is sharper.  In the initial `n=14` and
`n=18` rows, `t=1/n` has all vertices locally lonely and very large `H`, while
`t=1/(2n)` has `H=1`, a huge empty semicircle, and no marked stationary
loneliness.

The sampled hard-ladder times remain high-H but lose the marked stationary
witness.  For example, the `n=14` selected rows at `3053/25872` and
`4339/51744` both have `H=15649345`, yet the stationary runner is not lonely;
the lonely information is stored in which vertices have two safe adjacent gaps.

## Predictions

1. LRC scripts should report `H`, maximum circular gap, marked stationary
   distance, safe-gap mask, locally lonely vertex set, and pressure/deletion
   DAG data together.
2. High-H clock cells should correlate better with the number of local lonely
   vertices than with the marked stationary witness alone.
3. Low-H clock cells should remain exactly the right detector for unanchored
   semicircle bunching, but should not be used as a proxy for LRC success.
4. Counterexample-shaped LRC searches should be judged by marked tournament
   structure first: stationary vertex, adjacent safe gaps, endpoint protection,
   and pressure DAG/SCC status.

## Sources

- `04-computation/h_loneliness_meter_s506.py`
- `05-knowledge/results/h_loneliness_meter_s506.out`
- `07-reflections/h-as-loneliness-meter-s506.md`
- HYP-1951
- THM-370
- THM-374

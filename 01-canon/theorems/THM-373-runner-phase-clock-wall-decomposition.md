---
id: THM-373
name: runner-phase-clock-wall-decomposition
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S500
depends_on:
  - THM-372
---

# THM-373: The runner phase comparator gives a finite closed tournament walk

## Statement

Let `s_0,...,s_{n-1}` be distinct integer speeds.  Put runner `i` at

```text
x_i(t) = s_i t mod 1
```

on the unit circle.  Define the half-turn phase comparator by

```text
i -> j  iff  frac((s_j - s_i)t) in (0,1/2),
```

with collision and antipodal ties completed by a fixed Hamiltonian tie path.

Then:

1. the only times at which the arc between `i` and `j` can change are

```text
t = m / (2 |s_i - s_j|)
```

for integers `m`;

2. on `[0,1)` the comparator has finitely many wall times;

3. between consecutive wall times the resulting tournament is constant;

4. the tournament path is closed, i.e. the limiting clock repeats with period
   `1`.

## Proof

The relative phase of runners `i` and `j` is

```text
frac((s_j - s_i)t).
```

The strict half-turn decision changes only when this value crosses one of the
boundary values `0` or `1/2`.  Thus

```text
(s_j - s_i)t in Z        or        (s_j - s_i)t in Z + 1/2.
```

Equivalently,

```text
2 |s_j - s_i| t in Z,
```

so the wall times for the pair are exactly `m/(2|s_i-s_j|)`.

There are finitely many pairs and finitely many such residues modulo `1` for
each pair, so there are finitely many wall times in `[0,1)`.  On the complement
of this finite wall set, every pairwise strict comparison is constant.  By
THM-372, the resulting tournament is constant on each open cell between
consecutive wall times.

Finally, each `s_i` is an integer, so `x_i(t+1)=x_i(t)` for every `i`.  The full
configuration, and therefore the completed tournament, has period `1`.

## Significance

This formalizes the "tournament clock" from the recent Tournament Analysis
sessions.  A runner system is not only a scalar lonely-runner problem; under
the half-turn phase gauge it is a closed walk through the tournament cube, with
exact rational wall times controlled by speed differences.

## Related

- THM-372
- HYP-1931
- HYP-1951
- `04-computation/tournament_clock_s24.py`
- `07-reflections/tournament-clock-runner-walks-in-Gn-s24.md`

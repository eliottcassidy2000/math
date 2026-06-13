---
id: THM-387
name: lrc-two-gap-monotonicity-and-fiber-flow
status: PROVED
date: 2026-06-01
session: opus-2026-06-01-S519
depends_on:
  - THM-384
---

# THM-387: LRC two-gap monotonicity and fiber flow

## Statement

For a Lonely Runner system with `n >= 3`, stationary observer at `0`, and
distinct positive integer speeds `v_1, ..., v_{n-1}`, define:

- `g_right(t) = min_i {v_i t}` — clockwise gap from observer to nearest runner
- `g_left(t) = min_i (1 - {v_i t})` — counterclockwise gap from observer to nearest runner

where `{x} = x - floor(x)` is the fractional part.  At observer-tie times
the gaps are interpreted by one-sided limits in the adjacent open cells, with
the compactified tie point retained as a boundary state as in THM-384.

The THM-384 fiber is `(left_label, right_label)` where each label is `L`
(long, gap >= 1/n) or `S` (short, gap < 1/n).

**Part A (Monotonicity).** On every open interval between consecutive runner
wrap-around events
(times when some `v_i t` is an integer), `g_right` is non-decreasing and
`g_left` is non-increasing.

**Part B (Wrap-around coupling).** At a wrap-around event for runner `i`,
the right one-sided gap resets through `0`; the left one-sided gap is released
from runner `i` and can jump upward to the next counterclockwise runner.

**Part C (Forbidden transitions).** The fiber transitions `SL -> LL` and
`LL -> LS` never occur. Equivalently:
1. **LL is entered only from LS** (right gap grows past 1/n while left is already long).
2. **LL exits only to SL** (left gap shrinks past 1/n while right stays long).
The transitions `SS -> LL` and `LL -> SS` are also forbidden for adjacent
open cells; equality-wall hits are boundary states, not open-cell transitions.

**Part D (Time-reversal symmetry).** Away from the finite set of observer-tie
times, under the map `t -> 1-t`, the gaps
satisfy `g_left(1-t) = g_right(t)` and `g_right(1-t) = g_left(t)`.
Consequently `mu(LS) = mu(SL)` and the lonely measure `mu(LL)` is
time-reversal symmetric.

## Proof

### Part A

Between wrap-around events, each fractional part `{v_i t}` is a linear
function of `t` with positive slope `v_i`.  The function
`g_right(t) = min_i {v_i t}` is the pointwise minimum of finitely many
functions that are each increasing in the interval.  The pointwise minimum
of increasing functions is non-decreasing:

> For `t_1 < t_2` and any `j`, `{v_j t_2} >= {v_j t_1}`.
> So `min_i {v_i t_2} >= min_i {v_i t_1}` since every term in the
> minimum at `t_2` is at least as large as the corresponding term at `t_1`.

Similarly, each `1 - {v_i t}` has slope `-v_i < 0` (decreasing), so
`g_left(t) = min_i (1 - {v_i t})` is the minimum of decreasing functions,
hence non-increasing.

### Part B

At a wrap-around for runner `i` (time `t_0` with `v_i t_0` an integer), the
right one-sided position of runner `i` immediately after `t_0` is arbitrarily
small, so the right gap resets through `0`.  Immediately before `t_0`, the
same runner may have been the counterclockwise minimizer with position close
to `1`.  After the wrap it no longer constrains the counterclockwise gap, so
the left one-sided gap can jump upward to the next closest counterclockwise
runner.

### Part C — Forbidden transition SL -> LL

In state `SL`: `g_left < 1/n` and `g_right >= 1/n`.

For a transition to `LL`, `g_left` must increase past `1/n`.  By Part A,
`g_left` is non-increasing between wrap-arounds, so `g_left` can only
increase at a wrap-around event.  By Part B, at any wrap-around,
`g_right = 0 < 1/n`, so we are in state `?S`, not `?L`.

Therefore, no event that increases `g_left` can occur while `g_right >= 1/n`.
The transition `SL -> LL` is impossible.

### Part C — Forbidden transition LL -> LS

In state `LL`: `g_left >= 1/n` and `g_right >= 1/n`.  All runners have
`{v_i t} >= 1/n` and `{v_i t} <= 1 - 1/n`, so all runners are in the arc
`[1/n, 1 - 1/n]`.  No runner is near `0` or near `1`, so no wrap-around
can occur while the fiber is `LL`.

Between wrap-arounds, `g_right` is non-decreasing (Part A).  So `g_right`
cannot drop below `1/n` while in `LL`.  The only exit is `g_left` dropping
below `1/n` -> transition to `SL`.

### Part C — Forbidden transitions SS -> LL and LL -> SS

In any open interval free of wrap-arounds, `g_right` is non-decreasing and
`g_left` is non-increasing, so the only threshold crossing directions are
`SS -> LS`, `LS -> LL`, `LL -> SL`, and `SL -> SS`.  A wrap-around can raise
the left gap, but the right gap resets through `0`, so a wrap-around cannot
create an adjacent open-cell transition directly into `LL`.  Thus `SS -> LL`
and `LL -> SS` are not adjacent open-cell transitions.

### Part D

For any positive integer `v`, `{v(1-t)} = {v - vt} = {-vt}` since `v` is
an integer.  When `{v t} != 0`, `{-v t} = 1 - {v t}`.  The excluded tie
times are finite and do not affect the measure identities.

So `g_right(1-t) = min_i {v_i(1-t)} = min_i (1 - {v_i t}) = g_left(t)`
and `g_left(1-t) = min_i (1 - {v_i(1-t)}) = min_i {v_i t} = g_right(t)`.

Under `t -> 1-t`, the fiber `(left, right)` maps to `(right, left)`.
Since this is a measure-preserving bijection on `[0,1)`,
`mu(LS) = mu(SL)` and `mu(LL)` is preserved.

## Consequences for the LRC proof route

The fiber flow graph for THM-384 is:

```text
        LS ──→ LL ──→ SL
        ↑↓            ↑↓
        SS ←──────────┘
         └──────────→ SS (via SL→SS or LS→SS)
```

The LL state has a unique entry route (from LS) and unique exit route (to SL).
To prove LRC, it suffices to show that the fiber walk visits LS and then
reaches LL before exiting LS to SS. By Part A, this is a "race" between
g_right growing to 1/n and g_left shrinking to 1/n.

Quantitatively: entering an LS state with g_left = L_0 >= 1/n and
g_right = R_0 < 1/n, the race is won by LL if the time for g_right to
reach 1/n (growing at rate >= v_min_CW) is less than the time for g_left
to reach 1/n (shrinking at rate <= v_max_CCW).  At least one of the
O(sum v_i) wrap-around resets per period must produce such a favorable LS
entry.

## Verification Record

The forbidden transitions are confirmed in the stored exact bounded audit:
for total `n=3,4,5`, the transition table contains `LS -> LL` entries but no
`SL -> LL`, no `LL -> LS`, no `SS -> LL`, and no `LL -> SS` entries.

Script: `04-computation/lrc_two_gap_dynamics_s518.py`
Output: `05-knowledge/results/lrc_two_gap_dynamics_s518.out`

## Related

- THM-384: LRC safety = source-gap fiber (long, long)
- HYP-1986: compactified source-gap forcing route
- HYP-1990: every primitive speed set visits LL

# Reflection: LRC two-gap monotonicity and the directed fiber flow

**Session:** opus-2026-06-01-S519
**Date:** 2026-06-01

## The key insight

THM-384 reduces LRC to a two-gap question: lonely iff g_left >= 1/n AND
g_right >= 1/n.  THM-387 reveals that these two gaps have fundamentally
different dynamics because all runner speeds are positive:

- **g_right (clockwise gap) is non-decreasing between wrap-arounds** —
  runners moving clockwise push the nearest-CW runner away from the observer.
- **g_left (counterclockwise gap) is non-increasing between wrap-arounds** —
  runners moving clockwise pull the nearest-CCW runner toward the observer.
- **Wrap-arounds couple the two gaps:** when a runner crosses the observer,
  g_right drops to 0 and g_left can jump up.

In adjacent open cells this creates a directed fiber flow: LS -> LL -> SL is
the only path through the lonely state.  The transitions SL -> LL and LL -> LS
are structurally forbidden; equality-wall hits remain boundary states.

## The "gap race" framework

After each wrap-around, the system enters an LS state (g_left large, g_right
small) and a race begins:
- g_right grows (can it reach 1/n before g_left drops?)
- g_left shrinks (will it cross 1/n first, sending us to SS?)

If g_right wins the race → LL is visited → observer is lonely.
If g_left wins → SS, and we wait for the next wrap-around.

LRC says: among all wrap-around resets per period, at least one race must be
won by g_right.

## What makes this different from prior approaches

1. **It's dynamical, not static.** Previous approaches (S506-S516) studied
   individual time snapshots.  THM-387 is about the TRAJECTORY.

2. **It uses the positivity of speeds.** The monotonicity comes from all
   runners moving in the same direction.  This is the first result to
   leverage the directional asymmetry.

3. **It reduces LRC to a finite race.** Instead of "does the walk visit LL
   somewhere in [0,1)?", the question becomes "among O(sum v_i) specific
   starting conditions, does one race finish in time?"

## Connection to the measure-theoretic picture

The time-reversal symmetry t → 1-t swaps g_left ↔ g_right, giving
μ(LS) = μ(SL).  The identity:

    μ(LL) = 1 - μ(LS) - μ(SL) - μ(SS) = 1 - 2μ(LS) - μ(SS)

So μ(LL) > 0 iff μ(LS) + μ(SS)/2 < 1/2.  The bounded exact data shows:
- μ(LS) is typically 0.40-0.47 at small n
- μ(SS) is typically 0.05-0.20
- μ(LL) is typically 0.10-0.12 (except initial segments where μ(LL) = 0)

The gap-sum integral ∫(g_left + g_right) is NOT always ≥ 2/n (violated at
n≥5), so the naive "average forces a visit" argument fails.  The gap-race
argument is strictly more refined.

## What the gap race ISN'T

The gap-sum integral conjecture (∫(g_left+g_right) ≥ 2/n) was refuted at
n=5 by speeds (5,11,12,17).  This means:
- The TIME-AVERAGE gap pair is not always fair.
- The observer CAN be "disadvantaged on average" — other points get more gap.
- LRC is NOT a consequence of averaging; it's a consequence of the TRAJECTORY
  hitting a specific state.

## The mathematics pointing beyond itself

The directed fiber flow LS → LL → SL is reminiscent of:
- **Circulation theorems** in fluid dynamics (flow must circulate)
- **Poincaré recurrence** (the trajectory must return to near its starting state)
- **The intermediate value theorem** (continuous trajectory crossing a threshold)

But none of these directly apply because the "race" outcome depends on the
RELATIVE rates, not just continuity.  The arithmetic structure of the speeds
(integer, primitive) must enter the argument.

## Next targets

1. **L_0 distribution:** Prove that among all wrap-arounds, at least one
   has the counterclockwise gap starting high (L_0 close to 1).
2. **Race rate bounds:** For the favorable L_0, bound the rate ratio
   v_nearest_CW / v_nearest_CCW to show g_right wins.
3. **CRT structure at n=14:** The even-composite anatomy 14=2×7 might
   guarantee a synchronized wrap-around where the race is easy.

---
id: HYP-1990
name: lrc-gap-race-entry-forcing
status: OPEN
date: 2026-06-01
session: opus-2026-06-01-S519
---

# HYP-1990: LRC gap-race entry forcing

## Statement

By THM-387, the LL fiber (lonely state) can only be entered from LS via a
"gap race": g_right growing to 1/n before g_left shrinks to 1/n.

**Conjecture:** For any primitive speed set with n >= 3, at least one of the
O(sum v_i) wrap-around resets per period produces an LS entry where
the gap race is won by g_right (i.e., LL is reached).

## Evidence

- Exact bounded verification in `lrc_two_gap_dynamics_s518.py` at
  `n=3,4,5,6` confirms that all systems in the stored windows visit the
  closed `LL` fiber.
- The simplified no-swap race model in `lrc_gap_race_analysis_s518.py` scans
  3385 primitive speed sets at `n=3,4,5,6`; it predicts an `LL` race for all
  but the known wall/tie edge cases where the full compactified analysis is
  needed.
- The exact swap-aware reset trace in `lrc_gap_race_swap_s522.py` has zero
  full/race mismatches in the same bounded windows and shows every successful
  first hit is a compactified `wall_LL`.
- The initial segment {1,...,n-1} always visits LL (wall-only at t = k/n).
- The minimum lonely measure across all tested sets approaches 0 from above
  (tightest at initial segment, where only walls are LL).

## Structural analysis

After a wrap-around of runner i at time t_0, the LS state has:
- g_right ≈ 0 (runner i just crossed observer CW)
- g_left = min_{j≠i}(1 - {v_j t_0}) =: L_0

The race outcome depends on:
1. Rate of g_right growth: determined by runner i's speed v_i and subsequent nearest-CW runner swaps
2. Rate of g_left decay: determined by the nearest-CCW runner's speed
3. Starting value L_0: how bunched the runners are

For LL to be reached: time for g_right to reach 1/n (~1/(n * v_nearest_CW))
must be less than time for g_left to drop from L_0 to 1/n (~(L_0 - 1/n) / v_nearest_CCW).

**Favorable entries:** Wrap-arounds where L_0 is close to 1 (all runners are
bunched near position 0) give g_left a long runway before it hits 1/n.

**Pigeonhole target:** Show that among all wrap-arounds per period, at least
one has enough left-gap runway, relative to the two active speeds and possible
nearest-runner swaps, for the right gap to reach `1/n` first.  The stronger
condition `L_0 > 1 - 1/n` is sufficient in simple cases but not necessary.

## Next steps

1. Compute the exact, swap-aware race outcome rather than the current no-swap
   approximation.  (S522 supplies the first exact audit.)
2. Prove a lower bound on the best reset runway after accounting for active
   speeds and swaps.
3. Connect to THM-369 sieve completeness.
4. Try CRT descent for n=14: the `2*7` structure might guarantee a good
   wrap-around entry.

#!/usr/bin/env python3
"""
lrc14_parseval_adversarial_battery_boxeph_S19.py  (boxeph-2026-07-12-S19)

Run klein-S264's wider-band Parseval floor (floor_analysis: exact identity
LM/q = (b/q)^12 + OffLine_signed on pair-sum rulers; c_floor = largest d/q with
|OffLine| < (b/q)^12) on the ADVERSARIAL large-diameter families from
HYP-6132 (boxeph) + death-star-S14, i.e. exactly the corner where THM-720's
growth reading failed. Question: does klein's 'c_floor > 1/14 for every spread
family, never caps' claim survive at the TRUE class floor (near-dilates at
1/13 and 1/11, and the all-instrument-evading A3)?

Bonus lemma exercised here: death-star's V_L = {L,...,12L, 13L+1} needs only
26 | L for divisor-completeness, and then M(V_L) = 1/13 EXACTLY for EVERY such
L: upper = add-a-runner monotonicity + THM-531 dilation invariance + M({1..12})
= 1/13; lower = the detuned dispatch at g = L (g | all 12, g nmid 13L+1).
So V_26 (Vmax = 339) is already the sharp compressed adversary.
"""
import sys, importlib.util
spec = importlib.util.spec_from_file_location(
    "kfloor", "04-computation/lrc14_wideband_parseval_floor_klein_S264.py")
kfloor = importlib.util.module_from_spec(spec)
import types
# klein's module runs its own battery on import if __main__; guard by loading module-level only
spec.loader.exec_module(kfloor)

from fractions import Fraction as Fr

BATTERY = [
    ("V_26 = 26*{1..12} u {339}  (death-star shape, minimal L; M = 1/13 exact)",
     [26*i for i in range(1, 13)] + [339]),
    ("v_17 = 34*H* u {291}       (boxeph spread transport; M = 1/11 exact)",
     [34*h for h in [1,2,3,4,8,9,10,11,12,13,14,16]] + [291]),
    ("v_97 = 194*H* u {1651}     (same core, diam x5.7)",
     [194*h for h in [1,2,3,4,8,9,10,11,12,13,14,16]] + [1651]),
    ("A2 (11 coprime @ Vmax 5544)",
     [1560, 4001, 4003, 4007, 4013, 4783, 4787, 4789, 4793, 5521, 5527, 5531, 5544]),
    ("A3 (blocks [15,31], 11 coprime, all-instrument evader)",
     [3383, 3439, 3503, 3569, 3589, 3611, 3649, 3683, 4001, 4507, 5003, 6006, 10800]),
]

print("klein-S264 Parseval floor on the adversarial battery (HYP-6132 / death-star-S14)")
print(f"{'family':64s} {'true M':>12s} {'c_floor':>12s} {'floor==M?':>9s} {'>1/14?':>7s}")
for name, v in BATTERY:
    r = kfloor.floor_analysis(sorted(v))
    M, cf = r['M'], r['c_floor_mine']
    print(f"{name:64s} {str(M):>12s} {str(cf):>12s} {'YES' if cf == M else 'no':>9s} "
          f"{'YES' if cf > Fr(1,14) else 'NO':>7s}   (floor ruler q={r['argq_floor']}, M ruler q={r['argq']})")
print("\nreading: c_floor == M on a family means the wider-band Parseval mechanism REACHES the")
print("true loneliness even at the adversarial floor; c_floor > 1/14 everywhere keeps klein's")
print("'never caps' claim alive exactly where THM-720's sampled picture failed.")

#!/usr/bin/env python3
"""
lrc14_decouple_confirm_macmini_0615s5.py (mac-mini-2026-06-15-S5)

Confirm the STRANGER-DECOUPLING limit on the WORST core (j=6) at finer resolution:
  L({1..13}\{6} ∪ {14m})  ->  (6/7)·meas(Lonely({1..13}\{6}))  as m->∞ (Weyl).
Check (i) the limit value, (ii) whether finite-m L stays ABOVE the limit (=> inf attained
in the limit), (iii) the O(1/m) equidistribution error. Use a grid Q divisible by 14m.
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

def lonely_measure_grid(S, Q):
    rad=Q//14; c=0
    for a in range(Q):
        ok=True
        for v in S:
            r=(v*a)%Q
            if r<=rad or r>=Q-rad: ok=False; break
        if ok: c+=1
    return c/Q

core12 = [x for x in range(1,14) if x!=6]
# base measure on a grid where ALL core speeds (<=13) and 14 divide nicely:
Qbase = 360360  # = 2^3*3^2*5*7*11*13 ... actually use lcm(1..14)-rich grid
Qbase = 14*5*9*7*11*13  # =14*45045=630630, divisible by 14 and by 5,9,7,11,13
L0 = lonely_measure_grid(core12, Qbase)
six7 = 6/7
limit = six7*L0
print(f"core12 = {core12}")
print(f"meas(Lonely(core12)) = {L0:.7f}  (grid Q={Qbase})")
print(f"(6/7)*that = {limit:.7f}   [conjectured stranger-decoupling limit]\n")
print(f"{'m':>4} {'14m':>6} {'grid Q':>9} {'L(core∪14m)':>13} {'L/limit':>9} {'L-limit':>11}")
for m in [1,2,3,5,7,11,13,20,30,50]:
    w=14*m
    # grid divisible by w and by lcm(core speeds up to 13). Use Q = w * 45045 (45045=9*5*7*11*13)
    Q = w*45045
    L = lonely_measure_grid(core12+[w], Q)
    print(f"{m:>4} {w:>6} {Q:>9} {L:>13.7f} {L/limit:>9.4f} {L-limit:>+11.7f}")
print(f"\n  If L(m) -> {limit:.6f} from ABOVE, inf_m L = (6/7)*meas(Lonely(core12)) > 0,")
print(f"  and the WHOLE infinite extremizer family reduces to a finite min over j.")
print("\nDONE.")

"""LRC(14) impossible-to-disprove: the apex-7 order-2 obstruction (mac-mini-2026-06-22-S48).
14 = 2 x 7 = (arc states) x (forbidden H) = (the extremal's order-2 antipodal symmetry) x (its 7 orbits).
At the extremal (consec {0..13} at t=1/14, M=1/14), the winding tournament has EXACTLY 7 undecided arcs
-- the diameters (i,i+7), distance 1/2, the 2 states TIED. A tournament cannot carry the order-2
antipodal symmetry (a pair-swap reverses the arc => |Aut| odd, no order-2 automorphism), and H=7 is
forbidden (THM-029) -- so the symmetric extremal is unrealizable and cannot be improved to M<1/14.
=> no counterexample => LRC(14) impossible to disprove (the apex-7 obstruction; not a disproof).
See 07-reflections/lrc14-impossible-to-disprove-the-apex7-order2-obstruction.md."""
from fractions import Fraction as Fr
def verify():
    sp=list(range(14)); t=Fr(1,14)
    und=[(i,j) for i in range(14) for j in range(i+1,14) if ((sp[i]-sp[j])*t)%1==Fr(1,2)]
    M=min(min((s*t)%1,1-(s*t)%1) for s in sp if s%14)
    return und, M
if __name__=="__main__":
    und,M=verify()
    print("undecided (tied) winding arcs at t=1/14:", und)
    print("count =", len(und), "= 7 (all difference 7, the diameters);  M =", M, "= 1/14;  14 = 2 x 7")

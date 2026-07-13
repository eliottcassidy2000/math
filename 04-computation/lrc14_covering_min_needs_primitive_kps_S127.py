# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont57: the covering-min 14/183 is a PRIMITIVE-DC statement -- load-bearing witnesses.
#
# DILATION invariance M(cv)=M(v) (THM-531) but DC is NOT dilation-invariant. So an IMPRIMITIVE DC family
# c*v can have M(c*v)=M(v) < 14/183 while its primitive reduction v is NON-DC (=> THM-366 gives only >= 1/14).
# WITNESSES (all DC, all M < 14/183):
#   {2,4,..,26} = 2*{1..13}  (2*AP)  : DC, M = 1/14   (prim reduction {1..13} = tight AP, non-DC)
#   2*GW {2,..,22,26,48}             : DC, M = 1/14   (prim reduction = GW, non-DC)
#   {2,4,..,24,182} = 2*{1..12,91}   : DC, M = 7/92   (prim reduction {1..12,91}, non-DC, killer 91=7*13 only)
# => the covering-min crux is "PRIMITIVE DC => M >= 14/183"; without 'primitive' it is FALSE (imprimitive DC
#    reaches the global floor 1/14). This is WHY the fleet's "primitive q-covering" qualifier is essential, and
#    it CONTROLS opus-S253's open "large-s trade": a core with M_core=1/13 and binding speed s>=2 is a dilation
#    s*{1..12}, whose covering killer must be even (resonant + mult of lcm(13,14)=182) => imprimitive => barred.
from math import gcd
from fractions import Fraction as F
from functools import reduce

def norm(x): r = x - int(x); r = r + 1 if r < 0 else r; return min(r, 1 - r)
def Mexact(v, qcap=None):
    if qcap is None: qcap = 4 * max(v) + 4
    best = F(0); bestt = None
    for q in range(2, qcap):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(F(vi * p, q)) for vi in v)
                if m > best: best, bestt = m, F(p, q)
    return best, bestt
def is_cov(v, N=14): return all(any(x % d == 0 for x in v) for d in range(2, N + 1))
def g(v): return reduce(gcd, v)

def main():
    tgt = F(14, 183)
    print(f"PRIMITIVE covering-min = 14/183 = {float(tgt):.6f}; global LRC floor 1/14 = {float(F(1,14)):.6f}\n")
    print(f"{'family':<30} | gcd | DC? | prim-reduction | prim DC? | {'M':>7} | base q | bucket")
    fams = {
        "2*{1..13} = {2,4,..,26}": [2*i for i in range(1,14)],
        "2*GW {2,..,22,26,48}   ": [2*x for x in [1,2,3,4,5,6,7,8,9,10,11,13,24]],
        "2*{1..12,91}={2..24,182}": [2*i for i in range(1,13)]+[182],
        "deep well {1..12,182}  ": list(range(1,13))+[182],
        "AP {1..13} (primitive) ": list(range(1,14)),
    }
    for name, v in fams.items():
        gg = g(v); dc = is_cov(v); pv = [x//gg for x in v]; pdc = is_cov(pv)
        M, t = Mexact(v)
        bucket = ("PRIM-DC (crux, odd base)" if (gg==1 and dc) else
                  "PRIM non-DC (THM-366)" if gg==1 else
                  "IMPRIM DC -> reduces to prim non-DC")
        pvs = str(pv) if max(pv) < 30 else str(pv[:4])+"..+["+str(pv[-1])+"]"
        print(f"{name:<30} | {gg:>3} | {str(dc):>4} | {pvs:<14} | {str(pdc):>7} | {float(M):.5f} | {t.denominator:>5}{'(even)' if t.denominator%2==0 else '(odd)'} | {bucket}")
    print()
    print("READING: imprimitive DC (gcd>1) dilates to a NON-DC primitive => even base => down to 1/14. Only")
    print("PRIMITIVE DC lives at the odd doorway (base Phi6=183) with M>=14/183. The crux is primitive-only, and")
    print("primitivity is exactly what bars the low-M large-binding-speed cores (opus-S253's isolated large-s trade).")

if __name__ == "__main__":
    main()

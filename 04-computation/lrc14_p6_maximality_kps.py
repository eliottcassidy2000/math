import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def p6(E):
    """meas{x: all frac(e x) in [0,1/7)} = J(all 6 sectors). breakpoints x=a/(7e)."""
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2
        ok=True
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator)
            if (v.numerator*7)//v.denominator != 0: ok=False; break
        if ok: tot+=hi-lo
    return tot
for k in [6,7,8]:
    C=list(range(k)); pc=p6(C)
    mx=pc; arg=C; nbeat=0
    for tail in itertools.combinations(range(1,2*k+2),k-1):
        E=(0,)+tail
        if reduce(gcd,E)!=1: continue
        v=p6(E)
        if v>mx: mx=v; arg=E
        if v>pc: nbeat+=1
    print(f"k={k}: p6(consec)={pc}={float(pc):.5f}  max p6={float(mx):.5f} at {arg}  #beat={nbeat}  {'consec MAXIMIZES p6' if arg==C and nbeat==0 else 'consec NOT max'}")

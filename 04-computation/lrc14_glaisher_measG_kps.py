import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def measG_direct(E):
    """meas{x: frac(e x) in [0,1/7) for all e in E}, exact."""
    E=sorted(set(e for e in E if e>0)); bps={Fraction(0),Fraction(1)}
    for e in E:
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; ok=True
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator)
            if (v.numerator*7)//v.denominator!=0: ok=False; break
        if ok: tot+=hi-lo
    return tot
def glaisher_data(E):
    """odd-part b -> max 2-adic level a_b (only b's that are themselves in E)."""
    S=set(e for e in E if e>0); d={}
    for e in S:
        b=e; a=0
        while b%2==0: b//=2; a+=1
        if b in S:  # anchor present
            d[b]=max(d.get(b,0),a)
    return d
def measG_glaisher(E):
    """meas{x: frac(b x) in [0, 1/(7*2^a_b)) for all odd b in glaisher_data}, exact."""
    d=glaisher_data(E)
    if not d: return Fraction(1)
    bs=sorted(d); bps={Fraction(0),Fraction(1)}
    for b in bs:
        for a in range(0,b+1): bps.add(Fraction(a,b))  # frac(bx) breakpoints
        thr=Fraction(1,7*(2**d[b]))
        for a in range(0,b+1): bps.add(Fraction(a,b)+thr/b if (Fraction(a,b)+thr/b)<=1 else Fraction(1))
    # exact: integrate over fine grid of breakpoints
    allb=set(bps)
    for b in bs:
        thr=Fraction(1,7*(2**d[b]))
        for m in range(0,b+1):
            allb.add(Fraction(m,b)); 
            x=Fraction(m,b)+thr/b
            if 0<=x<=1: allb.add(x)
    allb=sorted(allb); tot=Fraction(0)
    for i in range(len(allb)-1):
        lo,hi=allb[i],allb[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; ok=True
        for b in bs:
            v=b*mid; v=v-(v.numerator//v.denominator)
            if v>=Fraction(1,7*(2**d[b])): ok=False; break
        if ok: tot+=hi-lo
    return tot
print("=== VERIFY: meas(G_E) [direct] == meas via Glaisher (odd b, level a_b)?  (odd-part-complete sets) ===")
tests=[ list(range(1,14)), [1,2,3,4,5,6,7], [1,3,5,7,9,11,13], [1,2,4,8,3,6,5], list(range(1,8)) ]
for E in tests:
    if not all((b in set(E)) for b in glaisher_data(E)): pass
    dG=measG_direct(E); gG=measG_glaisher(E)
    print(f"  E={E}: direct={dG} glaisher={gG}  [{'MATCH' if dG==gG else 'DIFFER'}]  data={glaisher_data(E)}")
print("\n=== meas(G) for consec {1..13} (the n=14 AP) ===")
E=list(range(1,14)); print(f"  meas(G_consec13) = {measG_direct(E)} = {float(measG_direct(E)):.6f}; glaisher data={glaisher_data(E)}")

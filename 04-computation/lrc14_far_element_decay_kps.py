import sys
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def p0(E):
    """meas(S7) = meas{x: all 7 sectors [j/7,(j+1)/7) hit by {frac(e x)}}."""
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if len(hit)==7: tot+=hi-lo
    return tot
capf={8:0.38153,9:0.49426,10:0.6044}
for k in [8,9,10]:
    core=list(range(k-1))  # {0,1,...,k-2}
    consec=list(range(k)); pc=p0(consec)
    print(f"k={k}: consec p0={float(pc):.5f} cap={capf[k]}  | family core{{0..{k-2}}} U {{w}}, w grows:")
    vals=[]
    for w in range(k-1, 4*k+10):
        E=core+[w]
        if len(set(E))!=k: continue
        if reduce(gcd,E)!=1: continue
        v=p0(E); vals.append((w,v))
    # print key w and the trend
    mx=max(vals,key=lambda t:t[1])
    for w,v in vals:
        if w<=k+3 or w%5==0 or (w,v)==mx:
            tag=" <-- MAX" if (w,v)==mx else (" =consec" if w==k-1 else "")
            print(f"    w={w:3d}: p0={float(v):.5f}{tag}")
    print(f"    => family max p0={float(mx[1]):.5f} at w={mx[0]}  (consec={float(pc):.5f}); decreasing tail? last p0={float(vals[-1][1]):.5f} at w={vals[-1][0]}")
    print()

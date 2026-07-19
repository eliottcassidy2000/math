# Robustness: does the SHARP horn ever FALSE-certify? It exhibits a real witness, so it can't,
# but let's confirm it certifies the tightest covering families as lonely (they ARE lonely).
from fractions import Fraction as F
def danger(v):
    w=F(1,14*v); return [(F(j,v)-w,F(j,v)+w) for j in range(v)]
def sub(safe,arcs):
    cuts=[]
    for lo,hi in arcs:
        lo%=1; hi%=1
        if lo<hi: cuts.append((lo,hi))
        else: cuts.append((lo,F(1))); cuts.append((F(0),hi))
    for clo,chi in sorted(cuts):
        new=[]
        for lo,hi in safe:
            if chi<=lo or clo>=hi: new.append((lo,hi)); continue
            if clo>lo: new.append((lo,clo))
            if chi<hi: new.append((chi,hi))
        safe=new
    return safe
def horn_certifies(speeds):
    # Remove all but the largest.  The residual is closed and the danger teeth
    # are open, so equality at 1/(7*kmax) certifies as well.
    speeds=sorted(speeds); kmax=speeds[-1]; rest=speeds[:-1]
    safe=[(F(0),F(1))]
    for v in rest: safe=sub(safe,danger(v))
    L=max((hi-lo for lo,hi in safe), default=F(0))
    return L >= F(1,7*kmax), L, kmax
# deep well and its tower (all lonely, M=14/183)
for fam in [list(range(1,13))+[182], list(range(1,13))+[364],
            [1,2,3,4,5,6,7,8,9,10,11,13,24],      # GW second tight (non-covering)
            list(range(1,14)),                      # tight AP (non-covering)
            [1,2,3,4,5,6,7,8,9,10,11,13,84]]:       # worst |core|=1 body
    ok,L,km=horn_certifies(fam)
    print(f"{fam}: sharp-horn certifies={ok}  L={float(L):.6g} kmax={km} 1/(7kmax)={1/(7*km):.6g}")

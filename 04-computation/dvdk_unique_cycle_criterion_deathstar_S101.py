"""S101: the sharp DvdK-free criterion for GMC(2) via one-dimensional coprime intervals.
CT(f_F^m*) is a SINGLE nonzero term (=> DvdK-free, coefficient-independent) exactly when the
lowest face has a UNIQUE minimal balanced channel = a unique tournament-zeta PRIMITIVE CYCLE."""
from functools import reduce
from math import gcd
import itertools

def min_channels(charges, cap=30):
    ch=sorted(set(charges)); n=len(ch)
    for m in range(1,cap+1):
        out=[]
        def rec(i,rem,acc,c):
            if i==n-1:
                if c+rem*ch[i]==0: out.append(tuple(acc+[rem]))
                return
            for k in range(rem+1): rec(i+1,rem-k,acc+[k],c+k*ch[i])
        rec(0,m,[],0)
        if out: return m,out
    return None,[]

print("[A] Unique minimal cycle => single nonzero CT(f^m*) => DvdK-FREE (support-only, any coeffs):")
for c in [[-1,1],[-1,2],[-1,1,2],[-1,2,3],[-1,3,5],[-1,0,1],[-6,10,15],[-2,3,7],[-1,1,3,5]]:
    m,ch=min_channels(c); print(f"   charges {c}: m*={m}, #cycles={len(ch)} -> {'FREE' if len(ch)==1 else 'HARD'}  cycle {ch[0] if ch else None}")

print("\n[B] Scan straddling supports (sizes 3-4, range +-4): fraction DvdK-free:")
pool=[-4,-3,-2,-1,1,2,3,4]; free=0; hard=[]; tot=0
for r in (3,4):
    for c in itertools.combinations(pool,r):
        if min(c)<0<max(c):
            tot+=1; m,ch=min_channels(list(c))
            if len(ch)==1: free+=1
            else: hard.append(c)
print(f"   {free}/{tot} = {100*free//tot}% DvdK-FREE; {len(hard)} HARD (>=2 coincident minimal cycles)")
print(f"   hard stratum examples: {hard[:6]} -- thin arithmetic-coincidence set (e.g. two antipodal pairs)")

print("\n[C] Connection: minimal balanced channel = fundamental PRIMITIVE CYCLE of the tournament/walk")
print("    zeta (THM-1926, zeta starts at u^ell = the shortest cycle). Unique shortest cycle => single")
print("    nonzero leading term. >=2 coincident shortest cycles = degenerate zeta = the ONLY DvdK-hard case.")

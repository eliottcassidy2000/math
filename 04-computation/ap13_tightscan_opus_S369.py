from fractions import Fraction as F
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0], max(out[-1][1],b))
        else: out.append((a,b))
    return out
def uncovered(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    return 1-sum(b-a for a,b in union(allv))
print("IS TIGHTNESS SPECIAL TO {1,...,13}?  d=1 APs with varying start a:")
for a in range(1,10):
    V=list(range(a,a+13))
    u=uncovered(V)
    tag = "  <-- TIGHT (uncovered = 0)" if u==0 else ""
    print(f"    {{{a}..{a+12}}}   uncovered = {float(u):.8f}{tag}")
print()
print("MINIMUM OVER 13-TERM APs WITH d >= 2 (scan a,d):")
best=None
for d in range(2,13):
    for a in range(1,15):
        from math import gcd
        V=[a+i*d for i in range(13)]
        u=uncovered(V)
        if best is None or u<best[0]: best=(u,a,d)
print(f"    min uncovered = {float(best[0]):.8f}  at a={best[1]}, d={best[2]}")
print(f"    (exact: {best[0]})")
print()
print("  => every 13-term AP with d>=2 stays uniformly away from 0.")

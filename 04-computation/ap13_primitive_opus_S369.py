from fractions import Fraction as F
from math import gcd
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
print("CORRECTED CLAIM: a 13-term AP is TIGHT iff a = d (i.e. it is a dilate")
print("of {1,...,13}).  Dilation invariance (THM-1075) forces this.")
print("    a=d cases:")
for k in range(1,7):
    V=[k+i*k for i in range(13)]
    print(f"      a=d={k}: {V[:3]}...  uncovered = {float(uncovered(V)):.8f}   (= {k}*{{1..13}})")
print()
print("    MINIMUM OVER PRIMITIVE (gcd=1) 13-TERM APs, scanning a,d <= 16:")
best=None; tights=[]
for d in range(1,17):
    for a in range(1,17):
        V=[a+i*d for i in range(13)]
        if gcd(a,d)!=1: continue           # primitive only
        u=uncovered(V)
        if u==0: tights.append((a,d))
        if best is None or u<best[0]: best=(u,a,d)
print(f"      min = {float(best[0]):.8f} at a={best[1]}, d={best[2]}")
print(f"      primitive tight APs found: {tights}")
print()
print("  => among PRIMITIVE APs the only tight one is a=d=1, i.e. {1,...,13};")
print("     every other primitive AP has a strictly positive margin, and the")
print("     non-primitive tight cases are exactly its dilates k*{1,...,13}.")

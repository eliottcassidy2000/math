# opus-2026-07-17-S376 -- HYP-7590: THE SIMULTANEOUS BLOCKING TRADE-OFF.
#
# THM-1110 showed no single modulus binds -- every q > 14 is individually
# blockable within the 13-speed budget -- so the constraint is SIMULTANEOUS
# blocking, whose only cheap route is divisibility (q | v kills every p).
#
# THE TENSION TO TEST.  Blocking many moduli forces speeds divisible by many
# things, i.e. LARGE and highly structured.  But a large speed has a very fine
# comb, which is close to equidistributed and therefore close to INDEPENDENT of
# the others -- and independence is the BEST case for loneliness (it is what
# makes the uncovered measure ~ (1-2*lam)^13 rather than 0).  So the adversary
# may be forced to buy blocking with measure.  If uncovered measure is bounded
# BELOW as a function of how much blocking is done, that is a real theorem.
from fractions import Fraction as F
from functools import reduce
from math import gcd
import random
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
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def uncovered(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    return 1-sum(b-a for a,b in union(allv))
def q0(V):
    q=1
    while any(v%q==0 for v in V): q+=1
    return q

print("(1) UNCOVERED MEASURE vs q0 -- does blocking more moduli cost measure?")
random.seed(376)
buckets={}
for _ in range(400):
    V=sorted(random.sample(range(1,260),13))
    if reduce(gcd,V)!=1: continue
    z=q0(V); u=float(uncovered(V))
    buckets.setdefault(min(z,20), []).append(u)
print("    q0   n    min uncovered   median   max")
for z in sorted(buckets):
    b=sorted(buckets[z])
    if len(b)<3: continue
    print(f"    {z:3d} {len(b):4d}   {b[0]:.6f}      {b[len(b)//2]:.6f}  {b[-1]:.6f}")

print()
print("(2) THE lcm FAMILIES -- maximal blocking, what does it cost them?")
def lcm2(a,b): return a*b//gcd(a,b)
for Q in [6,8,10,12]:
    L=reduce(lcm2,range(1,Q+1))
    V=sorted([L]+[3,5,11,13,17,19,23,29,31,37,41,43])
    print(f"    Q={Q:3d}  lcm={L:7d}  q0={q0(V):4d}  uncovered = {float(uncovered(V)):.6f}")

print()
print("(3) ADVERSARIAL MINIMISATION of uncovered measure over PRIMITIVE families")
print("    (the tight family {1..13} gives 0; is it isolated?)")
random.seed(3762)
best=(1.0,None)
for trial in range(25):
    V=sorted(random.sample(range(1,120),13))
    if reduce(gcd,V)!=1: continue
    cur=float(uncovered(V)); stall=0
    for step in range(220):
        W=list(V); i=random.randrange(13)
        W[i]=max(1,W[i]+random.choice([-4,-2,-1,1,2,4]))
        W=sorted(set(W))
        if len(W)!=13 or reduce(gcd,W)!=1: continue
        u=float(uncovered(W))
        if u<=cur:
            if u<cur: stall=0
            V,cur=W,u
        else: stall+=1
        if stall>140: break
    if cur<best[0]: best=(cur,list(V))
print(f"    minimum uncovered found: {best[0]:.8f}")
print(f"      V = {best[1]}")
print(f"      q0 = {q0(best[1]) if best[1] else '-'}")

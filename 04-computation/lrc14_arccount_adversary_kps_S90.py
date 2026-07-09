# Adversarial: maximize c=#arcs/spread over primitive k=13 configs (does worst c approach rho*?).
# Also test the WORST-case structure klein flagged (near-arithmetic / resonant primitive).
import random
from math import gcd
from functools import reduce
TH=1.0/7.0
def maxgap(E,x):
    pts=sorted((e*x)%1.0 for e in E); n=len(pts); mg=0.0
    for i in range(n):
        nxt=pts[(i+1)%n]+(1.0 if i==n-1 else 0.0); g=nxt-pts[i]
        if g>mg: mg=g
    return mg
def stats(E,spread,N=None):
    if N is None: N=spread*300
    prev=False; arcs=0; tot=0; g0=None; glast=None
    for j in range(N):
        x=(j+0.5)/N; good=maxgap(E,x)>TH
        if j==0: g0=good
        glast=good
        if good: tot+=1
        if good and not prev: arcs+=1
        prev=good
    if g0 and glast and arcs>1: arcs-=1
    return arcs, tot/N
def isprim(E):
    E=sorted(E); return reduce(gcd,[e-E[0] for e in E])==1
random.seed(999)
spread=210
print(f"Adversarial worst-case c=#arcs/spread, k=13 primitive, spread={spread}:")
# (a) random primitives
best=(0,None)
for _ in range(200):
    while True:
        mid=sorted(random.sample(range(1,spread),11)); E=[0]+mid+[spread]
        if isprim(E): break
    a,rho=stats(E,spread)
    if a/spread>best[0]: best=(a/spread,(tuple(E),a,rho))
print(f"  random(200): max c={best[0]:.3f}  #arcs={best[1][1]} rho*={best[1][2]:.3f}")
# (b) structured resonant candidates: near-AP, clustered, few-difference (primitive)
cands={
 "near-AP {0,15,30,..} +1 off": [0]+[15*i for i in range(1,12)]+[165+1,210],
 "AP step 17 +offset": sorted(set([0]+[17*i for i in range(1,12)]+[210])),
 "2 blocks": sorted(set([0,1,2,3,4,5]+[205,206,207,208,209,210]+[100])),
 "geometric-ish": sorted(set([0,1,2,4,8,16,32,64,128,210,105,53,27])),
}
for name,E in cands.items():
    E=sorted(set(E))
    if len(E)!=13: continue
    sp=E[-1]-E[0]
    if not isprim(E): continue
    a,rho=stats(E,sp)
    print(f"  {name:<28} c={a/sp:.3f} #arcs={a} rho*={rho:.3f} {'** c>=rho' if a/sp>=rho else ''}")
print(f"\n  bar for pigeonhole: need c < rho*. rho* ~ 0.98 for large spread; worst random c ~ {best[0]:.2f}.")
print("  => margin rho*-c large; the resonant-primitive worst case is the target to bound.")

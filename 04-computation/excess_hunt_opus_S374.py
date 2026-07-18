# adversarial hunt on the EXCESS  min-den - q0  (the salvageable conjecture)
from math import gcd
from functools import reduce
import random
def lonely_at(V,p,q):
    for v in V:
        r=(v*p)%q
        if min(r,q-r)*14 < q: return False
    return True
def works(V,q): return any(lonely_at(V,p,q) for p in range(1,q) if gcd(p,q)==1)
def min_den(V,Qmax=600):
    for q in range(1,Qmax+1):
        if works(V,q): return q
    return None
def q0(V):
    q=1
    while any(v%q==0 for v in V): q+=1
    return q
def exc(V):
    m=min_den(V)
    return (m-q0(V)) if m else 999
print("(7) ADVERSARIAL HUNT ON THE EXCESS  E(V) = min-den - q0")
random.seed(37421)
best=(0,None)
for trial in range(40):
    V=sorted(random.sample(range(2,400),13))
    if reduce(gcd,V)!=1: continue
    cur=exc(V); stall=0
    for step in range(400):
        W=list(V); i=random.randrange(13)
        W[i]=max(1,W[i]+random.choice([-6,-4,-2,-1,1,2,4,6]))
        W=sorted(set(W))
        if len(W)!=13 or reduce(gcd,W)!=1: continue
        e=exc(W)
        if e>=cur:
            if e>cur: stall=0
            V,cur=W,e
        else: stall+=1
        if stall>200: break
    if cur>best[0]: best=(cur,list(V))
print(f"    max excess found: {best[0]}")
print(f"      V = {best[1]}")
if best[1]:
    print(f"      q0 = {q0(best[1])}, min-den = {min_den(best[1])}")
print()
print("(8) THE lcm CONSTRUCTION'S EXCESS (my S374 refutation family)")
def lcm(a,b): return a*b//gcd(a,b)
for Q in [10,20,30,40]:
    L=reduce(lcm,range(1,Q+1))
    V=sorted([L]+[3,5,11,13,17,19,23,29,31,37,41,43])
    print(f"    Q={Q:3d}: q0={q0(V):4d}  min-den={min_den(V,300)}  excess={exc(V)}")

from fractions import Fraction as F
from functools import reduce
from math import gcd
from itertools import combinations
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def maximizer(V):
    Dn=set()
    for a,b in combinations(V,2):
        Dn.add(a+b)
        if a!=b: Dn.add(abs(a-b))
    for v in V: Dn.add(2*v)
    best=F(0); arg=None
    for q in sorted(Dn):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,arg=g,F(p,q)
    return best,arg
def active(V,t,g):
    up=[];dn=[]
    for v in V:
        x=v*t; a=int(x); r=x-a
        if r<0: a-=1; r+=1
        if r==g: up.append((v,a))
        if 1-r==g: dn.append((v,a+1))
    return up,dn
def q0(V):
    q=1
    while any(v%q==0 for v in V): q+=1
    return q
def info(V):
    g,t=maximizer(V); up,dn=active(V,t,g)
    if not up or not dn: return g,None,None,None
    vi,ai=up[0]; vj,aj=dn[0]; D=vi*aj-vj*ai
    return g,D,vi+vj,(vi,vj)
TARGET=F(1,14)
print("(4) ADVERSARIAL MINIMISATION RESTRICTED TO THE HARD STRATUM (q0 > 14)")
print("    -- the branch where the sieve is too weak and a D>=2 pair must rescue")
random.seed(3921)
def make_hard():
    # force every q in 1..14 to divide some speed
    base=[8,9,10,11,12,13,14]
    rest=[random.randint(2,90) for _ in range(6)]
    V=sorted(set(base+rest))
    while len(V)<13: V=sorted(set(V+[random.randint(2,90)]))
    return V[:13]
best=(F(1),None)
n=0
for trial in range(30):
    V=make_hard()
    if len(V)!=13 or reduce(gcd,V)!=1 or q0(V)<=14: continue
    cur,_,_,_=info(V); stall=0; n+=1
    for step in range(140):
        W=list(V); i=random.randrange(13)
        W[i]=max(2,W[i]+random.choice([-4,-2,-1,1,2,4]))
        W=sorted(set(W))
        if len(W)!=13 or reduce(gcd,W)!=1 or q0(W)<=14: continue
        g2,_,_,_=info(W)
        if g2<=cur:
            if g2<cur: stall=0
            V,cur=W,g2
        else: stall+=1
        if stall>90: break
    if cur<best[0]: best=(cur,list(V))
g,D,s,pr=info(best[1])
print(f"    {n} hard-stratum starts")
print(f"    minimum gap over the hard stratum: {best[0]} = {float(best[0]):.8f}")
print(f"      vs 1/14 = {float(TARGET):.8f}  -> margin factor {float(best[0]/TARGET):.3f}x")
print(f"      V = {best[1]}")
print(f"      q0={q0(best[1])}  D={D}  s={s}  pair={pr}  s/D={F(s,D) if D else '-'}  (need s<=14D)")
print()
print("    => the sieve at q0 would give only 1/q0 = "
      f"{float(F(1,q0(best[1]))):.6f} < 1/14; the D>=2 pair delivers {float(g):.6f} instead")

# opus-2026-07-19-S392 -- HYP-7780: THE D >= 2 BRANCH.
#
# THM-1210 split LRC(14) into two mechanisms: D = 1 (the classical sieve, q0 <= 14,
# containing the extremal families at s = 14) and D >= 2 (the hard stratum).
# In the D >= 2 branch LRC(14) needs  s = v_i + v_j <= 14*D  at the active pair,
# and measured ratios s/D sit at 5-6.5 -- a factor >2 of slack.  So the branch
# breaks only if some family drives g = D/s down toward 1/14.
#
# DIRECT ATTACK: adversarially MINIMISE g(V) over primitive 13-families and read
# the boundary in (D, s) coordinates.  If the minimum is 1/14, attained only at
# the known tight families, that both supports LRC(14) and characterises the
# extremal set in the new coordinates.
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
def report(V):
    g,t=maximizer(V); up,dn=active(V,t,g)
    if not up or not dn: return g,t,None,None,None
    vi,ai=up[0]; vj,aj=dn[0]; D=vi*aj-vj*ai
    return g,t,D,vi+vj,(vi,vj)

TARGET=F(1,14)
print(f"    LRC(14) threshold 1/14 = {float(TARGET):.8f}")
print()
print("(1) ADVERSARIAL MINIMISATION of g over PRIMITIVE 13-families")
random.seed(392)
best=(F(1),None)
for trial in range(22):
    V=sorted(random.sample(range(1,45),13))
    if reduce(gcd,V)!=1: continue
    cur,_ ,_,_,_ = report(V); stall=0
    for step in range(160):
        W=list(V); i=random.randrange(13)
        W[i]=max(1,W[i]+random.choice([-3,-2,-1,1,2,3]))
        W=sorted(set(W))
        if len(W)!=13 or reduce(gcd,W)!=1: continue
        g2,_,_,_,_ = report(W)
        if g2<=cur:
            if g2<cur: stall=0
            V,cur=W,g2
        else: stall+=1
        if cur==TARGET: break
        if stall>110: break
    if cur<best[0]: best=(cur,list(V))
g,t,D,s,pr = report(best[1])
print(f"    minimum gap found: {best[0]} = {float(best[0]):.8f}   ({'== 1/14' if best[0]==TARGET else 'BELOW 1/14!' if best[0]<TARGET else 'above 1/14'})")
print(f"      V = {best[1]}")
print(f"      t*={t}  D={D}  s={s}  pair={pr}  ratio s/D = {F(s,D) if D else None}")

print()
print("(2) THE KNOWN EXTREMALS in (D,s) coordinates, and near-misses")
for nm,V in [("{1,...,13}",list(range(1,14))),
             ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("{1..12,14}",list(range(1,13))+[14]),
             ("{1..12,15}",list(range(1,13))+[15]),
             ("{2,...,14}",list(range(2,15))),
             ("{1..11,13,25}",[1,2,3,4,5,6,7,8,9,10,11,13,25])]:
    g,t,D,s,pr=report(V)
    print(f"    {nm:16s} g={str(g):9s} ({float(g):.6f})  D={D}  s={s}  pair={pr}"
          f"  s/D={F(s,D) if D else '-'}   {'<== EXTREMAL' if g==TARGET else ''}")

print()
print("(3) WHAT THE MINIMISERS LOOK LIKE -- do they all contain a long run 1..k?")
V=best[1]
runs=1; mx=1
for a,b in zip(V,V[1:]):
    runs = runs+1 if b==a+1 else 1
    mx=max(mx,runs)
print(f"    longest consecutive run in the minimiser: {mx}")
print(f"    tight family {list(range(1,14))[:4]}... has run 13; {[1,2,3,4,5,6,7,8,9,10,11,13,24][:4]}... has run 11")

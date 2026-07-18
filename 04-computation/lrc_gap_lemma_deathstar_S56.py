from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
random.seed(1)
def M_exact(fam):
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best
def covers(fam,S): return all(any(v%q==0 for v in fam) for q in S)
def is_AP(W): W=sorted(W); d=W[0]; return W==[d*i for i in range(1,13)]
C=2366; base=list(range(1,13)); pool=[v for v in range(1,31) if v%13 and v%14]
worst=None; tested=0; viol=[]
def check(W):
    global worst,tested,viol
    W=sorted(set(W))
    if len(W)!=12 or is_AP(W) or any(v%13==0 or v%14==0 for v in W): return
    if not covers(W,range(2,13)): return
    M=M_exact(W); delta=M-F(1,13)
    if delta<=0: return
    tested+=1; margin=delta-F(max(W),C)
    if worst is None or margin<worst[0]: worst=(margin,tuple(W),float(delta),max(W))
    if margin<=0: viol.append((tuple(W),float(delta),max(W)))
# k=1 exhaustive
for p in range(12):
    for x in pool:
        W=base[:]; W[p]=x; check(W)
# k=2 sampled
for _ in range(8000):
    p=random.sample(range(12),2); x=random.sample(pool,2)
    W=base[:]; W[p[0]]=x[0]; W[p[1]]=x[1]; check(W)
print(f"valid non-AP cores tested: {tested}, violations (delta<=max/2366): {len(viol)}")
if worst:
    m,W,d,mx=worst
    print(f"tightest: W={list(W)} delta={d:.5f} max={mx} max/2366={mx/C:.5f} margin={float(m):.5f} => {'HOLDS' if float(m)>0 else 'VIOLATED'}")
for v in viol[:5]: print("  VIOL:",v)

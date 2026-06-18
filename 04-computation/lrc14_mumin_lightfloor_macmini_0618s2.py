# Light bounded-spread minima for k=8..13 (small spread, fast exact) to report concrete floors.
from fractions import Fraction as F
import itertools
from math import gcd
TWO7=F(2,7)
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def meas(arcs): return sum((b-a for a,b in arcs), F(0))
def good_set_exact(E):
    E=sorted(set(E)); k=len(E)
    diffs=set()
    for a in range(k):
        for b in range(a+1,k): diffs.add(E[b]-E[a])
    bps={F(0),F(1)}
    for d in diffs:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1)
    good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=sorted(((E[t]*xm)%1, E[t]) for t in range(k))
        order=[e for _,e in pts]; floors=[int((e*xm)//1) for e in order]
        for idx in range(k):
            e_cur=order[idx]; f_cur=floors[idx]
            if idx<k-1: e_nx=order[idx+1]; f_nx=floors[idx+1]; wrap=F(0)
            else: e_nx=order[0]; f_nx=floors[0]; wrap=F(1)
            A=F(e_nx-e_cur); Cc=F(f_cur-f_nx)+wrap
            if A==0:
                if Cc>TWO7: good.append((x0,x1))
                continue
            xb=(TWO7-Cc)/A
            if A>0: lo=max(x0,xb); hi=x1
            else: lo=x0; hi=min(x1,xb)
            if lo<hi: good.append((lo,hi))
    return merge(good)
def mu(E): return meas(good_set_exact(E))
def prim(E):
    g=0
    for e in E: g=gcd(g,e)
    return [e//g for e in E] if g>1 else list(E)
# exhaustive small spread for k=8..13
for k in range(8,14):
    bb={8:14,9:14,10:15,11:15,12:15,13:14}[k]
    best=(F(2),None)
    for rest in itertools.combinations(range(1,bb+1),k-1):
        E=[0]+list(rest)
        if prim(E)!=E: continue
        m=mu(E)
        if m<best[0]: best=(m,E)
    print(f"k={k}: bounded-spread<= {bb} exhaustive min mu = {best[0]} = {float(best[0]):.6f} at E={best[1]}")

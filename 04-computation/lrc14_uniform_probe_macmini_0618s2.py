from fractions import Fraction as F
import random
random.seed(99)
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
from math import gcd
def prim(E):
    g=0
    for e in E: g=gcd(g,e)
    return [e//g for e in E] if g>1 else list(E)

# Candidate uniform bounds. For sorted 0=e1<...<ek=s, define:
#   d2 = e2 (smallest positive offset)
#   dtop = s - e_{k-1}  (top gap)
#   dmin = min consecutive gap
# Conjecture to test: mu(E) >= 5/(7 * d2) ?  (near x=a/d2 the offset e2 sweeps; near small the cluster forms)
# Test which of these is a valid uniform LOWER bound.
worst={'5/(7 e2)':99,'5/(7 dtop)':99}
for trial in range(40000):
    k=random.randint(4,13)
    sp=random.choice([k,k+1,2*k,3*k,5*k,30,60,120,250])
    rest=random.sample(range(1,sp+1),min(k-1,sp))
    E=prim([0]+rest)
    if len(set(E))<3: continue
    E=sorted(set(E)); s=max(E)
    e2=E[1]; dtop=s-E[-2]
    m=mu(E)
    r1=float(m*F(7*e2,5))   # m / (5/(7 e2)) = m*7e2/5
    r2=float(m*F(7*dtop,5))
    if r1<worst['5/(7 e2)']: worst['5/(7 e2)']=r1; w1=E
    if r2<worst['5/(7 dtop)']: worst['5/(7 dtop)']=r2; w2=E
print('min ratio mu/(5/(7 e2)) =',round(worst['5/(7 e2)'],3),'(>=1 means valid bound) at',w1)
print('min ratio mu/(5/(7 dtop))=',round(worst['5/(7 dtop)'],3),'at',w2)

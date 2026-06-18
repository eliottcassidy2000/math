# Bulletproof L2: for many random E, verify (0, 5/(7s)) subset Good and that it is TIGHT
# (the near-0 arc ends exactly when the top point e_max*x first creates a <=2/7 max-gap).
from fractions import Fraction as F
import random
from math import gcd
random.seed(5)
TWO7=F(2,7)
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
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
fails=0; checked=0
for _ in range(3000):
    k=random.randint(2,13)
    sp=random.choice([k,2*k,5*k,50,200])
    rest=random.sample(range(1,sp+1),min(k-1,sp))
    E=sorted(set([0]+rest)); s=max(E)
    g=good_set_exact(E)
    w=F(5,7*s)
    a0,b0=g[0]
    if not (a0==0 and b0>=w):
        fails+=1
        if fails<=3: print("FAIL", E, "first arc",(a0,b0),"need>=",w)
    checked+=1
print(f"L2 check: {checked} random sets, {fails} failures of '(0,5/(7s)) subset Good'.")
print("L2 is", "CONFIRMED EXACT" if fails==0 else "VIOLATED")

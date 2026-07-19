from fractions import Fraction as F
from itertools import combinations
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def Mfull(V):
    D=set()
    for a,b in combinations(V,2):
        D.add(a+b)
        if a!=b: D.add(abs(a-b))
    for v in V: D.add(2*v)
    best=F(0); arg=None
    for q in sorted(D):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,arg=g,F(p,q)
    return best,arg
def act(V,t,g):
    up=[];dn=[]
    for v in V:
        x=v*t; a=int(x); r=x-a
        if r<0: a-=1; r+=1
        if r==g: up.append((v,a))
        if 1-r==g: dn.append((v,a+1))
    return up,dn
LO=F(1,14); HI=F(2,27)
print("THE LADDER  {1,...,11,13,12m}")
print("   m    x=12m   M          D    s    14D   slack=14D-s   in gap (1/14,2/27)?")
for m in range(2,11):
    x=12*m
    V=sorted(set(list(range(1,12))+[13,x]))
    if len(V)!=13: continue
    M,t=Mfull(V); up,dn=act(V,t,M)
    if up and dn:
        vi,ai=up[0]; vj,aj=dn[0]; D=vi*aj-vj*ai; s=vi+vj
    else: D=s=None
    ing = LO<M<HI
    print(f"   {m:2d}  {x:6d}   {str(M):9s} {str(D):4s} {str(s):4s} {str(14*D) if D else '-':5s}"
          f" {str(14*D-s) if D else '-':11s}  {'YES' if ing else ''}")
print()
print("  (M = 1/14 is the extremal; anything strictly between 1/14 and 2/27 populates")
print("   the stability gap, and slack = 14D - s is the distance from violating LRC(14))")

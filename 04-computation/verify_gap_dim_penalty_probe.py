from fractions import Fraction as F
import itertools
from math import gcd
from functools import reduce
def N_at(E,x):
    h=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator); h.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in h)
def dist_p(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[F(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        p[N_at(E,(lo+hi)/2)]+=(hi-lo)
    return p
def Ly(E,k):
    p=dist_p(E)
    g=[(F((t-1)*(t-2)*(t-4)*(t-5),40) if k==8 else F(-(t-2)*(t-3)*(t-6),36)) for t in range(7)]
    return sum(p[t]*g[t] for t in range(7))
def prim(E): return reduce(gcd,sorted(set(E)))==1
def isAP(E):
    E=sorted(set(E))
    if len(E)<=2: return True
    d=E[1]-E[0]; return all(E[i+1]-E[i]==d for i in range(len(E)-1))
# k=8 exhaustive non-AP max<=14
known={8:F(297,980),10:F(3307,5880)}
for k,W in [(8,14),(10,15)]:
    best=None;cnt=0
    for combo in itertools.combinations(range(1,W+1),k-1):
        E=(0,)+combo
        if not prim(E) or isAP(E): continue
        cnt+=1; L=Ly(E,k)
        if best is None or L>best[0]: best=(L,E)
    over = best[0]>known[k]
    print(f"k={k} W<={W}: {cnt} sets, max={float(best[0]):.6f} at {list(best[1])}; exceeds claimed? {over}  (claim {float(known[k]):.6f})")

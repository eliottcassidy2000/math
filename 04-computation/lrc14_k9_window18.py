import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
def N_at(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)
def dist_p(E):
    E=sorted(set(E)); bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        p[N_at(E,(lo+hi)/2)]+=(hi-lo)
    return p
g=[Fraction(-(t-2)*(t-3)*(t-6),36) for t in range(7)]
def L_y(E): 
    p=dist_p(E); return sum(p[t]*g[t] for t in range(7))
def is_ap(E):
    E=sorted(E); d=E[1]-E[0]
    return all(E[i+1]-E[i]==d for i in range(len(E)-1))
cap=0.49426; W=18; best=(Fraction(-1),None); cnt=0
for combo in itertools.combinations(range(1,W+1),8):
    E=(0,)+combo
    if reduce(gcd,E[1:])!=1: continue
    if is_ap(E): continue
    cnt+=1; Lv=L_y(E)
    if Lv>best[0]: best=(Lv,E)
Lv,E=best
print(f"k=9 W={W}: scanned {cnt} non-AP sets; max L_y={float(Lv):.6f} at {list(E)} margin={cap-float(Lv):.6f}")

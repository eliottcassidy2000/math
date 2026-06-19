import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def N_at(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)
def dist_p(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        p[N_at(E,(lo+hi)/2)]+=hi-lo
    return p
def g_poly(k,t):
    if k==8: return Fraction((t-1)*(t-2)*(t-4)*(t-5),40)
    if k in(9,10): return Fraction(-(t-2)*(t-3)*(t-6),36)
    return Fraction((t-3)*(t-4),12)
def L_y(E,k): p=dist_p(E); return sum(p[t]*g_poly(k,t) for t in range(7))
def excess(E):
    E=set(E); return len({a+b for a in E for b in E})-(2*len(E)-1)
def primitive(E): return reduce(gcd,E)==1
caps={8:Fraction(2402,6300),9:0.49426,10:0.6044}  # use float caps
capf={8:0.38153,9:0.49426,10:0.6044}
# For k=9,10: max L_y over excess>=1 (non-full-AP) primitive sets, bounded spread.
for k,box in [(9,15),(10,14)]:
    C=list(range(k)); Lc=L_y(C,k)
    best=(Fraction(0),None); cnt=0
    for tail in itertools.combinations(range(1,box+1),k-1):
        E=(0,)+tail
        if not primitive(E): continue
        if excess(E)==0: continue
        cnt+=1
        L=L_y(E,k)
        if L>best[0]: best=(L,E)
    print(f"k={k}: consec={float(Lc):.5f} cap={capf[k]} | max non-AP={float(best[0]):.5f} at {best[1]} | "
          f"gap(consec-nonAP)={float(Lc-best[0]):.5f} | margin(cap-nonAP)={capf[k]-float(best[0]):.5f} | scanned {cnt}")

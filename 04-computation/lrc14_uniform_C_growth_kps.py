import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def G0(y):
    f=y-(y.numerator//y.denominator)
    return Fraction(6,7)*f if f<=Fraction(1,7) else Fraction(6,49)-(f-Fraction(1,7))/7
def cells_of(Ep):
    Ep=sorted(set(Ep)); bps={Fraction(0),Fraction(1)}
    for e in Ep:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); out=[]
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in Ep:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        m=[j for j in range(1,7) if j not in hit]
        if len(m)==1:
            if out and out[-1][2]==m[0] and out[-1][1]==lo: out[-1]=(out[-1][0],hi,m[0])
            else: out.append((lo,hi,m[0]))
    return out
def wD(cells,w):
    return sum(G0(w*b-Fraction(s,7))-G0(w*a-Fraction(s,7)) for (a,b,s) in cells)
def supw(Ep,wmax):
    cells=cells_of(Ep); best=(Fraction(0),None)
    for w in range(2,wmax):
        if reduce(gcd,tuple(Ep)+(w,))!=1: continue
        v=abs(wD(cells,w))
        if v>best[0]: best=(v,w)
    return best
print("=== growth of sup_w |w*Delta_w| with core span (does the uniform bound hold?) ===")
print("  core (span) : sup_w|wD|  at w   | #N1cells")
fams = [
  ("consec8", list(range(8))),
  ("odd-struct s11", [0,1,3,5,7,9,10,11]),
  ("odds5+ends", [0,1,3,5,7,9,11]),
  ("odds7", [0,1,3,5,7,9,11,13]),
  ("odds9", [0,1,3,5,7,9,11,13,15,17]),
  ("dilate2 consec", [2*i for i in range(8)]),
  ("dilate3 consec", [3*i for i in range(8)]),
  ("wide {0,1,2,20,21,22,40}", [0,1,2,20,21,22,40]),
  ("multiscale", [0,1,2,30,31,32,60,61]),
  ("AP-step2 long", [0,2,4,6,8,10,12,14,16,18]),
  ("sidon-ish", [0,1,3,7,12,20,30,44]),
]
for name,Ep in fams:
    if reduce(gcd,tuple(x for x in Ep if x))==0: pass
    cells=cells_of(Ep)
    wm = 6*max(Ep)+30
    best=supw(Ep, wm)
    print(f"  {name:24s} span={max(Ep):3d} : sup|wD|={float(best[0]):.4f}={best[0]}  w={best[1]}  #cells={len(cells)}")

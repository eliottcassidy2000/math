import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd, comb
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
def Ly(p,k): return sum(p[t]*g_poly(k,t) for t in range(7))
def Lyinf(k):
    y={8:[1,-1,1,Fraction(-9,10),Fraction(3,5)],9:[1,Fraction(-13,18),Fraction(4,9),Fraction(-1,6)],10:[1,Fraction(-13,18),Fraction(4,9),Fraction(-1,6)]}[k]
    s=Fraction(0)
    for r,yr in enumerate(y):
        s+=yr*comb(6,r)*Fraction((7-r),7)**k
    return s
caps={8:Fraction(2402,6300),9:None,10:None}  # approx; use floats
capf={8:0.38153,9:0.49426,10:0.60440}
for k in [8,9,10]:
    C=list(range(k)); pc=dist_p(C); p0c=pc[0]; Lyc=Ly(pc,k)
    print(f"k={k}: consec p0={float(p0c):.5f} L_y={float(Lyc):.5f}  L_y^inf={float(Lyinf(k)):.5f}  cap={capf[k]}")
    # search: maximize p0 and L_y over bounded E (0 in E, entries<=2k+4), find if any beats consec p0
    maxp0=p0c; argp0=C; maxLy=Lyc; argLy=C; nbeat_p0=0; beat_ex=None
    rng=list(range(1,2*k+5))
    cnt=0
    for tail in itertools.combinations(rng,k-1):
        E=(0,)+tail
        if reduce(gcd,E)!=1: continue
        cnt+=1
        p=dist_p(E); p0=p[0]; L=Ly(p,k)
        if p0>maxp0: maxp0=p0; argp0=E
        if p0>p0c: 
            nbeat_p0+=1
            if beat_ex is None: beat_ex=(E,float(p0))
        if L>maxLy: maxLy=L; argLy=E
    print(f"   scanned {cnt} sets (entries<= {2*k+4}): max p0={float(maxp0):.5f} at {argp0}; #sets with p0>consec={nbeat_p0} ex={beat_ex}")
    print(f"   max L_y={float(maxLy):.5f} at {argLy}  (consec is {'MAX' if argLy==C else 'NOT max'} for L_y)")

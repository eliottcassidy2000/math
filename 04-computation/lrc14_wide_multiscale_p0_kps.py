import sys, itertools, random
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def p0(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if len(hit)==7: tot+=hi-lo
    return tot
def primitive(E): return reduce(gcd,tuple(E))==1
cap9=Fraction(1979,4004)
print(f"cap_9 = 1979/4004 = {float(cap9):.5f}")
print("=== p_0(E) for WIDE / MULTISCALE k=9 sets (does the CONCLUSION p_0<cap hold despite C growing?) ===")
tests = [
  [0,1,2,30,31,32,60,61,90],          # 3-cluster multiscale + far
  [0,1,2,30,31,32,60,61,62],          # 3 full clusters
  [0,1,3,5,7,9,10,11,90],             # the C=2.54 core + far
  [0,1,2,20,21,22,40,41,80],
  [0,1,2,100,101,102,200,201,202],    # very wide 3-cluster
  [0,1,2,3,300,301,302,303,600],
  [0,5,10,15,20,25,30,35,40],         # wide AP (dilated) -> = consec by scale-inv
  [0,1,2,3,4,5,6,7,8],                # consec (reference, the max)
]
mx=(Fraction(0),None)
for E in tests:
    if len(set(E))!=9: 
        print(f"  {E}: size!=9"); continue
    pr=primitive(E); v=p0(E)
    if v>mx[0]: mx=(v,E)
    print(f"  {E}: p_0={float(v):.5f}={v}  prim={pr}  {'<cap OK' if v<cap9 else '*** >= cap ***'}  margin={float(cap9-v):.4f}")
# broader random wide multiscale sample
print("\n=== random wide/multiscale primitive k=9 sample: max p_0 ===")
rng=random.Random(20260619); worst=(Fraction(0),None)
for _ in range(3000):
    nclust=rng.choice([2,3,4]); base=0; E=[]
    scale=rng.choice([20,30,50,100])
    for c in range(nclust):
        for j in range(rng.choice([2,3])):
            E.append(c*scale+j)
    E=sorted(set(E))
    while len(E)<9: 
        x=rng.randint(0,nclust*scale+5)
        if x not in E: E.append(x)
    E=sorted(E)[:9]
    if len(E)!=9 or not primitive(E): continue
    v=p0(E)
    if v>worst[0]: worst=(v,tuple(E))
print(f"  max p_0 over random wide multiscale = {float(worst[0]):.5f} at {worst[1]}  (cap_9={float(cap9):.5f}, margin {float(cap9-worst[0]):.4f})")

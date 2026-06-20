import sys, itertools, random
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def dist_p(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        p[sum(1 for j in range(1,7) if j not in hit)]+=hi-lo
    return p
def Phi(E):
    p=dist_p(E); return p[0]+Fraction(1,7)*p[1]
def prim(E): return reduce(gcd,tuple(E))==1
# Q(8) = max over bounded 8-sets of Phi (from earlier: 0.36210 at consec_8)
Q8=Fraction(141899,635376)  # placeholder; compute consec_8 Phi exactly
consec8=list(range(8)); Q8=Phi(consec8)
print(f"Phi(consec_8) = {Q8} = {float(Q8):.5f}  (= Q(8), the conjectured plateau max)")
print("=== is Phi(F) <= Phi(consec_8) for ALL 8-sets F, including WIDE? ===")
worst=(Fraction(0),None); cnt=0
# bounded scan
for tail in itertools.combinations(range(1,14),7):
    F=(0,)+tail
    if not prim(F): continue
    cnt+=1; v=Phi(F)
    if v>worst[0]: worst=(v,F)
print(f"  bounded (span<=13, {cnt} sets): max Phi = {float(worst[0]):.5f} at {worst[1]}  {'= consec OK' if worst[1]==tuple(consec8) else 'NOT consec'}")
# wide random scan
wworst=(Fraction(0),None)
rng=random.Random(7)
for _ in range(4000):
    F=sorted(rng.sample(range(1,60),7)); F=[0]+F
    if len(set(F))!=8 or not prim(F): continue
    v=Phi(tuple(F))
    if v>wworst[0]: wworst=(v,tuple(F))
print(f"  wide random (span<=60): max Phi = {float(wworst[0]):.5f} at {wworst[1]}  {'<= Q(8) OK' if wworst[0]<=Q8 else '*** EXCEEDS Q(8) ***'}")
# multiscale clusters
mworst=(Fraction(0),None)
for _ in range(2000):
    sc=rng.choice([20,30,50]); ncl=rng.choice([2,3]); F=[]
    for c in range(ncl):
        for j in range(rng.choice([2,3,4])): F.append(c*sc+j)
    F=sorted(set(F))
    while len(F)<8:
        x=rng.randint(0,ncl*sc+5)
        if x not in F: F.append(x)
    F=tuple(sorted(F)[:8])
    if len(F)!=8 or not prim(F): continue
    v=Phi(F)
    if v>mworst[0]: mworst=(v,F)
print(f"  multiscale clusters: max Phi = {float(mworst[0]):.5f} at {mworst[1]}  {'<= Q(8) OK' if mworst[0]<=Q8 else '*** EXCEEDS Q(8) ***'}")
print(f"\n  => if all <= Q(8)={float(Q8):.5f}: the plateau is maximized at bounded consec (clean inductive ingredient).")

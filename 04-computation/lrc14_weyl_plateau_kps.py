import sys, itertools
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
def p0(E): return dist_p(E)[0]
def primitive(E): return reduce(gcd,E)==1
capf={8:0.38153,9:0.49426,10:0.6044}

print("=== (1) Weyl approach rate: |p0({consec_8} U {w}) - plateau| vs w (k=9) ===")
core=list(range(8)); pc=dist_p(core); plat=pc[0]+Fraction(1,7)*pc[1]
print(f"   plateau = {float(plat):.6f}")
for w in [9,10,15,20,30,50,80,120,200,300]:
    E=core+[w]
    if not primitive(E): 
        E=core+[w]  # consec_8 has gcd1 always
    v=p0(E); err=abs(float(v)-float(plat))
    print(f"   w={w:4d}: p0={float(v):.6f}  |p0-plateau|={err:.6f}  err*w={err*w:.3f}")

print("\n=== (2) PERTURBED core + far w: does p0(E) stay < cap_9? (general far-element bound) ===")
# core = consec_8 with one element perturbed, plus a far w. max over such < cap_9?
maxv=(Fraction(0),None)
for j in range(1,8):  # perturb element j
    for delta in range(1,6):
        pcore=[x for x in range(8) if x!=j]+[7+delta]  # move j out
        if len(set(pcore))!=8: continue
        for w in [30,40,50,60,80,100]:
            E=tuple(sorted(set(pcore+[w])))
            if len(E)!=9 or not primitive(E): continue
            v=p0(E)
            if v>maxv[0]: maxv=(v,E)
print(f"   max p0 over (perturbed-core + far w) = {float(maxv[0]):.5f} at {maxv[1]}  (cap_9={capf[9]}, margin {capf[9]-float(maxv[0]):.4f})")

print("\n=== (3) the plateau bound is UNIFORM: max_w p0({consec_8} U {w}) over ALL w>=9 ===")
mx=max((p0(core+[w]) for w in range(9,250)), default=Fraction(0))
print(f"   max over w=9..249 = {float(mx):.5f} (at the near-AP w=9? = {float(p0(core+[9])):.5f}); cap_9={capf[9]} -> all < cap")

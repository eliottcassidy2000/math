#!/usr/bin/env python3
"""cont.43: (1) k=8 base compact+tail structure (parity with k=9 cont.42): deg-3 majorant
UPPER bound on Phi, must stay <= cap_9=0.4943; hardest at consec {0..7}; check no large-d
family exceeds. (2) decorrelation rate: |J9 - J9_iid|*diam for block+far tail -> constant C."""
from fractions import Fraction as F
def pvec(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return p
def bound8(E):  # deg-3 majorant value (UPPER bound on Phi)
    p=pvec(E)
    def q3(N): return 1 - F(2,3)*N + F(47,252)*N*(N-1) - F(5,252)*N*(N-1)*(N-2)
    return sum(p[N]*q3(N) for N in range(8))
def J9(E):
    p=pvec(E); return sum(p[n]*n*(7-n) for n in range(8))
cap9=F(1979,4004)
from itertools import combinations
print(f"=== (1) k=8 base: deg-3 bound must stay <= cap_9 = {float(cap9):.4f} (hardest at consec {{0..7}}) ===")
print("compact exhaustive (8-sets {0}uS, S subset {1..d}, |S|=7, diam d), MAX bound per d:")
gmax=(F(-1),None)
for d in range(7,18):
    best=(F(-1),None)
    for S in combinations(range(1,d),6):
        E=[0]+list(S)+[d]
        b=bound8(E)
        if b>best[0]: best=(b,E)
    if best[0]>gmax[0]: gmax=best
    print(f"  d={d:2d}: max bound = {float(best[0]):.4f}  margin {float(cap9-best[0]):+.4f}  {'*** EXCEEDS' if best[0]>cap9 else ''}")
# structured extremal candidates at large diameter (MISTAKE-138)
print("  structured large-d candidates (max bound; exceed => break k=8):")
cand=[("block+far d=50",[0,1,2,3,4,5,6,50]),("block+far d=200",[0,1,2,3,4,5,6,200]),
      ("2AP+far d=80",[0,2,4,6,8,10,12,80]),("mod7 {1,8..50}",[1,8,15,22,29,36,43,50]),
      ("block+2far",[0,1,2,3,4,5,100,101])]
cmax=gmax[0]
for nm,E in cand:
    b=bound8(E)
    if b>cmax: cmax=b
    print(f"    {nm:20s} bound={float(b):.4f}  {'*** EXCEEDS cap_9' if b>cap9 else 'ok'}")
print(f"  => k=8 GLOBAL max bound (compact+structured) = {float(cmax):.4f} <= cap_9={float(cap9):.4f}: {cmax<=cap9}, margin {float(cap9-cmax):+.4f}")
print()
print("=== (2) decorrelation rate |J9 - J9_iid|*diam for block+far {0..7,d} ===")
Jiid=F(0)  # use float iid from cont.42
Jiid_f=8.4560
print(f"  J9_iid = {Jiid_f}")
for d in [20,40,80,160,320,640]:
    E=[0,1,2,3,4,5,6,7,d]; j=float(J9(E))
    print(f"  d={d:4d}: J9={j:.4f}  |J9-Jiid|*d = {abs(j-Jiid_f)*d:.2f}")
print("  (|J-Jiid|*d bounded => rate C/d; D0 = C/(Jiid-thr) = C/3.71 for the tail bound)")

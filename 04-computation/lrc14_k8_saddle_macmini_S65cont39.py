#!/usr/bin/env python3
"""cont.39: does the k=8 cubic base (THM-714) have the SAME isolated-saddle structure
as k=9 (THM-716)? The deg-3 requirement is a cubic in the moments. Test: (a) is the
optimal-deg3-majorant value J3(E) minimized at consec {1..8}? (b) is it an isolated
saddle (low-'coverage' and high-'sync' both fail to beat it)?"""
from fractions import Fraction as F
from itertools import combinations
import random
# optimal deg-3 majorant q3 = 1 - 2N/3 + 47 N(N-1)/252 - 5 N(N-1)(N-2)/252 (THM-714)
def bound3(E):
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
    # bound = sum p[N] q3(N)  (an UPPER bound on Phi; k=8 closes iff bound <= cap_9)
    def q3(N): return 1 - F(2,3)*N + F(47,252)*N*(N-1) - F(5,252)*N*(N-1)*(N-2)
    return sum(p[N]*q3(N) for N in range(8))
cap9=F(1979,4004)
print(f"k=8: optimal deg-3 majorant bound vs cap_9={float(cap9):.4f} (closes iff bound<=cap9)")
print("consec-shift curve (want unique MAX at shift 1 = {1..8}, the hardest):")
worst=(F(-1),None)
for sh in range(0,10):
    E=list(range(sh,sh+8)); b=bound3(E)
    if b>worst[0]: worst=(b,sh)
    print(f"  cshift{sh} {{{sh}..{sh+7}}}: bound={float(b):.4f}  margin {float(cap9-b):+.4f}")
print(f"  => hardest (max bound) at shift {worst[1]}, bound={float(worst[0]):.4f}")
# adversarial: does anything push bound ABOVE cap9 (break k=8)?
random.seed(8)
gmax=(F(-1),None)
def climb():
    E=sorted(random.sample(range(0,26),8)); best=bound3(E)
    for _ in range(500):
        i=random.randrange(8); nd=random.choice([-2,-1,1,2])
        E2=E[:]; E2[i]+=nd
        if len(set(E2))!=8 or min(E2)<0 or max(E2)>55: continue
        v=bound3(sorted(E2))
        if v>best: best=v; E=sorted(E2)
    return E,best
print("\nadversarial MAX-bound hill-climbs (break k=8 if bound > cap_9):")
for t in range(8):
    E,b=climb()
    if b>gmax[0]: gmax=(b,E)
    print(f"  run {t}: bound={float(b):.4f}  {'*** BREAKS' if b>cap9 else 'ok'} at {E}")
print(f"\nGLOBAL max bound: {float(gmax[0]):.4f} at {gmax[1]}  (cap_9={float(cap9):.4f})")
print(f"k=8 SADDLE CONFIRMED (consec-max, nothing exceeds cap_9): {gmax[0]<=cap9}, margin {float(cap9-gmax[0]):+.4f}")

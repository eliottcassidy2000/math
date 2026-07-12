from math import gcd
from fractions import Fraction as F
from itertools import combinations
import random
def lcm(a,b): return a//gcd(a,b)*b
def T56(E):
    nz=[abs(e) for e in E if e]; has0=len(nz)<len(E)
    L=1
    for e in nz: L=lcm(L,e)
    D=7*L; pts=set([0,D])
    for e in nz:
        st=L//e; pts.update(range(0,D+1,st))
    pts=sorted(pts); pn=[0]*8; b0=1 if has0 else 0
    for t1,t2 in zip(pts,pts[1:]):
        s=t1+t2; hit=b0
        for e in nz: hit|=1<<((7*e*s//(2*D))%7)
        pn[7-bin(hit).count("1")]+=t2-t1
    p=[F(x,D) for x in pn]; T=[sum(p[j:]) for j in range(8)]
    return T[5],T[6]
def norm(E):
    g=0
    for e in E: g=gcd(g,e)
    return tuple(sorted(e//g for e in E)) if g>1 else tuple(sorted(E))
k=9; M=7*k-6  # = 57
bT6=F(2,M); bT5=F(5,M)
print(f"k={k}: conjecture T6 <= 2/{M} = {float(bT6):.6f}, T5 <= 5/{M} = {float(bT5):.6f}")
# search: exhaustive mod-7 (all coherence classes) + adversarial + all AP steps
maxT6=None;maxT6s=None;maxT5=None;maxT5s=None;viol=0
cands=[]
for step in range(1,20):
    for off in range(1,step+1 if step>1 else 2):
        cands.append(tuple(off+step*i for i in range(k)))
# all mod-q coherent for q=2..14
for q in range(2,15):
    for r in range(q):
        cands.append(tuple(r+q*i for i in range(k)) if r>0 else tuple(q*(i+1) for i in range(k)))
random.seed(3)
for _ in range(6000): cands.append(tuple(sorted(random.sample(range(1,120),k))))
seen=set()
for E in cands:
    E=norm(E)
    if len(set(E))<k or 0 in E: continue
    if E in seen: continue
    seen.add(E)
    t5,t6=T56(E)
    if t6>bT6 or t5>bT5: viol+=1
    if maxT6 is None or t6>maxT6: maxT6=t6;maxT6s=E
    if maxT5 is None or t5>maxT5: maxT5=t5;maxT5s=E
print(f"searched {len(seen)} families:")
print(f"  MAX T6 = {maxT6} ~ {float(maxT6):.6f} at {maxT6s}  {'<= bound' if maxT6<=bT6 else '*** EXCEEDS 2/57 ***'}")
print(f"  MAX T5 = {maxT5} ~ {float(maxT5):.6f} at {maxT5s}  {'<= bound' if maxT5<=bT5 else '*** EXCEEDS 5/57 ***'}")
print(f"  violations of (T6<=2/57 AND T5<=5/57): {viol}")
print(f"  => BUNCH = 2T5+4T6 <= 2*5/57+4*2/57 = 18/57 = 6/19 = {float(F(18,57)):.5f}: "
      f"{'holds' if maxT6<=bT6 and maxT5<=bT5 else 'CHECK'}")

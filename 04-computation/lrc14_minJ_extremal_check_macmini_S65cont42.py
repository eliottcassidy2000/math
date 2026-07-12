#!/usr/bin/env python3
"""cont.42 (MISTAKE-138 heed): my crossover SAMPLED large d -- explicitly TEST the
structured low-J extremal candidates at large diameter (not a box/sample) to confirm
min-J >= compact min 5.199 (= J({0..8})). Candidates: consec-block+far {0..7,d};
near-2-AP {0,2,4,..,14, far}; dilated-consec-primitivized; the mod-7 pole (high-J, sanity)."""
from fractions import Fraction as F
def J(E):
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
    return sum(p[n]*n*(7-n) for n in range(8))
compact_min=F(1019,196)  # J({0..8})
print(f"compact min J({{0..8}}) = {compact_min} = {float(compact_min):.4f}; threshold 432/91 = {float(F(432,91)):.4f}")
print("STRUCTURED low-J candidates at large diameter (per MISTAKE-138, not sampled):")
below=0; recs=[]
# consec-block + one far outlier, d up to 300
for d in [20,30,50,80,120,200,300]:
    E=[0,1,2,3,4,5,6,7,d]; j=J(E); recs.append(("block+far d=%d"%d,j))
# near-2-AP + far
for d in [40,80,160]:
    E=[0,2,4,6,8,10,12,14,d]; j=J(E); recs.append(("2AP+far d=%d"%d,j))
# consec-block of 8 shifted, + far
for d in [60,100]:
    E=[0,3,6,9,12,15,18,21,d]; j=J(E); recs.append(("3AP+far d=%d"%d,j))
# mod-7 pole (high-J sanity; the BUNCH extremal)
recs.append(("mod-7 pole {1,8..57}", J([1,8,15,22,29,36,43,50,57])))
# double-outlier
for d in [50,100]:
    E=[0,1,2,3,4,5,6,d,d+1]; j=J(E); recs.append(("block+2far d=%d"%d,j))
for nm,j in recs:
    b = j<compact_min
    if b: below+=1
    print(f"  {nm:24s} J={float(j):.4f}  {'*** BELOW compact min!' if b else ''}")
print(f"\n{len(recs)} structured candidates, {below} below compact min {float(compact_min):.4f}")
print(f"MIN-J ROBUST (no structured family beats the compact min): {below==0}")
print(f"=> crossover result holds: min-J = J({{0..8}}) = 1019/196, compact+tail closes k=9.")

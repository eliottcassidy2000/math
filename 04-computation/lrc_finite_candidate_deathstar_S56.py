from fractions import Fraction as F
from math import gcd
def M_exact(fam):
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best
def candidates_vmax(W):
    # stability: v_max <= max(W)/(13 delta); covering(miss 13,14): v_max = mult of 182 in [max(vmax_min), bound]
    M=M_exact(W); delta=M-F(1,13)
    if delta<=0: return M,[]
    ub=max(W)/(13*delta)   # v_max <= this
    lo=max(182, max(W)+1)  # v_max > max(W) (largest) and >=182
    out=[k for k in range(182, int(ub)+1, 182) if k>=lo]
    return M,out

print("Finite check: for each violating non-AP core W={1..11,X}, check ALL candidate v_max families:")
allclosed=True
for X in [16,17,18,19,20,22]:
    W=sorted([1,2,3,4,5,6,7,8,9,10,11,X])
    if len(W)!=12: continue
    MW,cands=candidates_vmax(W)
    results=[]
    for vm in cands:
        V=sorted(W+[vm])
        if len(V)!=13: continue
        MV=M_exact(V)
        results.append((vm,MV,MV<F(1,13)))
    allge = all(not r[2] for r in results)
    print(f"  W={{1..11,{X}}}: M(W)={float(MW):.4f} delta={float(MW-F(1,13)):.5f}  candidate v_max={cands}")
    for vm,MV,lt in results:
        print(f"      V=W+{{{vm}}}: M(V)={MV}={float(MV):.4f}  M<1/13? {lt}")
    if not allge: allclosed=False
print(f"\n=> all violating cores eliminated by finite check (every candidate gives M>=1/13)? {allclosed}")
print("   Mechanism: stability(v_max<=max/(13d)) + covering(182|v_max) => FINITE candidates per core, all fail.")

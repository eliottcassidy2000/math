from fractions import Fraction as F
from math import gcd
def M_arg(fam):
    Q=2*max(fam)+2; best=F(0); arg=None
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q); arg=(a,q)
    return best,arg
W=[1,2,3,4,5,6,7,8,9,10,11,24]
MW,_=M_arg(W); delta=MW-F(1,13)
ub=int(max(W)/(13*delta))
cands=[k for k in range(182,ub+1,182) if k>max(W)]
print(f"core W={W}: M(W)={MW}={float(MW):.5f} delta={float(delta):.5f}  candidate window [182,{ub}] => v_max in {cands}")
allge=True
for vm in cands:
    V=sorted(W+[vm]); MV,(a,q)=M_arg(V)
    lt = MV<F(1,13)
    if lt: allge=False
    print(f"  V=W+{{{vm}}}: M(V)={MV}={float(MV):.5f}  @t={a}/{q}  M<1/13? {lt}  {'*** COUNTEREXAMPLE ***' if lt else ''}")
print(f"\ncore {'ELIMINATED (all M>=1/13)' if allge else 'IS A COUNTEREXAMPLE to the inverse theorem!'}")
# also the OTHER violating cores
for W in [[1,2,3,4,5,6,7,8,9,10,12,22],[1,2,3,4,5,6,7,8,9,11,20,24]]:
    MW,_=M_arg(W); delta=MW-F(1,13)
    if delta<=0: continue
    ub=int(max(W)/(13*delta)); cands=[k for k in range(182,ub+1,182) if k>max(W)]
    res=[(vm,M_arg(sorted(W+[vm]))[0]) for vm in cands]
    ok=all(mv>=F(1,13) for _,mv in res)
    print(f"core {W}: cands={cands} -> {[(vm,float(mv)) for vm,mv in res]} all>=1/13:{ok}")

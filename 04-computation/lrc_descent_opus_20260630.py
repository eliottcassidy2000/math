"""
The 2-adic descent product. Chain of odd cores O_j = {s/2^j : v_2(s)=j}. Per-level apex gap
g(O_j mod 7) = min_{k!=0} |sum_{x in O_j} omega^{kx}|^2 (THM-590). Product, per-level min, cusp (Z_7) hits.
"""
import math, cmath
w=cmath.exp(2j*cmath.pi/7)
def g(Oset):  # apex gap of a subset of Z_7 (as a set of residues, with multiplicity)
    res=Oset  # list of residues mod 7
    best=None
    for k in range(1,7):
        val=abs(sum(w**((k*x)%7) for x in res))**2
        best=val if best is None else min(best,val)
    return best
def descent(S):
    chain=[]; cur=list(S)
    while cur:
        O=[x for x in cur if x%2==1]
        E=[x//2 for x in cur if x%2==0]
        if O: chain.append(O)
        cur=E
        if len(chain)>20: break
    return chain
def analyze(S,name):
    ch=descent(S)
    rows=[]
    prod=1.0; pmin=1.0; cusp_levels=[]
    for j,O in enumerate(ch):
        res=[x%7 for x in O]
        gv=g(res)
        prod*=gv if gv>1e-12 else 0.0
        if gv>1e-9: pmin=min(pmin,gv)
        isZ7 = set(res)==set(range(7))
        if isZ7: cusp_levels.append(j)
        rows.append((j,O,sorted(set(res)),round(gv,4),'CUSP Z_7' if isZ7 else ''))
    print(f"\n{name}: S={S}")
    for r in rows:
        print(f"   level {r[0]}: O={r[1]}  mod7={r[2]}  g={r[3]}  {r[4]}")
    print(f"   product prod g = {prod:.5f} | per-level min g (nonzero) = {pmin:.5f} | cusp at levels {cusp_levels}")
    return prod,pmin,cusp_levels
# the AP (extremal, M=1/14, NON-covering)
analyze(list(range(1,14)),"AP {1..13} (extremal)")
# the covering-min construction {1..12, 182}
analyze(list(range(1,13))+[182],"covering {1..12,182} (M=14/183)")
# the other covering construction {1..11,13,84}
analyze(list(range(1,12))+[13,84],"covering {1..11,13,84} (M=7/89)")
# a generic covering set: must contain mult of each q in 2..14
analyze([2,3,4,5,6,7,8,9,10,11,12,13,14],"{2..14} (covering, all q present)")
analyze([1,2,3,4,5,6,7,8,9,10,11,12,14],"{1..12,14} (covering via 14)")

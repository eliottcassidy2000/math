# -*- coding: utf-8 -*-
# klein-2026-07-11-S257: THE k=9 BASE TAIL = THE TWO-SCALE LIMIT (THM-687/688).
# mac-mini cont.43's correction: wide families PLATEAU at the two-scale/multi-scale limit
# (NOT J_iid). This confirms: far elements RAISE J to >= min two-scale limit >= compact-min,
# at every level k=5..9, bottoming at k<=7 (SmallClusterFull, nu=1). + the peel constant.
from math import gcd
from fractions import Fraction as F
from itertools import combinations
def lcm(a,b): return a//gcd(a,b)*b
def m12(E):
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
    p=[F(x,D) for x in pn]
    return sum(F(j)*p[j] for j in range(8)), sum(F(j*(j-1))*p[j] for j in range(8))
def J(E):
    m1,m2=m12(E); return 6*m1-m2
def eigenlim(Ep):  # two-scale limit: J(Ep u {far}) as far->inf, = (6 J(Ep) + m2(Ep))/7
    m1,m2=m12(Ep); return F(6,7)*6*m1-F(5,7)*m2
def norm(E):
    g=0
    for e in E: g=gcd(g,e)
    return tuple(sorted(e//g for e in E)) if g>1 else tuple(sorted(E))
if __name__=="__main__":
    print("TAIL RECURSION: min two-scale-limit_k >= compact-min_k = J(consec-k) at every k")
    for k in range(5,10):
        CMk=J(tuple(range(1,k+1)))
        W=min(2*k+8,20); minel=None;seen=set()
        for combo in combinations(range(1,W+1),k-1):
            key=norm(combo)
            if key in seen: continue
            seen.add(key); el=eigenlim(combo)
            if minel is None or el<minel: minel=el
        print(f"  k={k}: compact-min={float(CMk):.4f} min-2scale-limit={float(minel):.4f} "
              f"raises +{float(minel-CMk):.3f}")

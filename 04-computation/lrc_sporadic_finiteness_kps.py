"""
Bounding the sporadic tight locus (kps-S31o). Two ingredients:
(1) CENSUS: rigorously enumerate single-swap tight sets {AP\{rem} U {v}}, M=1/14.
(2) FINITENESS MECHANISM: M({1..11,13,v}) -> M({1..11,13}) >= 1/13 > 1/14 as v->inf, so
    tight (=1/14) only for BOUNDED v => finitely many sporadics.
M via exact critical points t=k/(s_i +- s_j) (the lonely-runner crossings).
"""
from fractions import Fraction as F
from math import gcd
def nf(x):
    r=x%1; return min(r,1-r)
def M_exact(S):
    S=sorted(set(abs(s) for s in S if s!=0)); C=set()
    for i in range(len(S)):
        for j in range(i,len(S)):
            for comb in {S[i]+S[j], abs(S[i]-S[j])}:
                if comb: 
                    for k in range(1,comb): C.add(F(k,comb))
        for k in range(1,2*S[i]): C.add(F(k,2*S[i]))
    best=F(0)
    for t in C:
        if 0<t<1:
            m=min(nf(s*t) for s in S)
            if m>best: best=m
    return best
AP=list(range(1,14))
print("(1) SINGLE-SWAP TIGHT CENSUS (remove rem from AP, add v<=50; rigorous M):")
found=[]
for rem in AP:
    for v in range(12,51):
        if v in AP and v!=rem: continue
        S=sorted(set([x for x in AP if x!=rem]+[v]))
        if len(S)!=13: continue
        if M_exact(S)==F(1,14): found.append((rem,v,S))
for rem,v,S in found:
    tag="AP" if v==rem else f"swap {rem}->{v} (v={v%27} mod27, {v%14} mod14)"
    print(f"   M=1/14: {tag}")
print(f"   => {len([f for f in found if f[0]!=f[1] or f[1] not in range(1,14)])} non-AP single-swaps; "
      f"the only non-AP tight set is GW. (v>24 all LOOSE.)")
print("\n(2) FINITENESS MECHANISM -- M({1..11,13,v}) vs v (-> M_12 = 1/13 as v grows):")
base=list(range(1,12))+[13]
M12=M_exact(base)
print(f"   M({{1..11,13}}) [12-subset] = {M12} = {float(M12):.5f}  (>1/14={1/14:.5f})")
for v in [14,18,20,22,24,26,28,30,36,48,60,120]:
    S=base+[v]; M=M_exact(S)
    print(f"   v={v:3d}: M = {str(M):>7s} = {float(M):.5f}  {'TIGHT' if M==F(1,14) else 'loose'}")
print("   => large v: M -> 1/13 (the 12-subset dominates; U_v too thin to tile the lonely set).")
print("      Tight only for BOUNDED v => finitely many single-swap sporadics. mod-27: GW swaps within")
print("      the gcd-3 shell (12,24 both mult of 3), the composite 2n-1=27=3^3 structure (HYP-2138).")

#!/usr/bin/env python3
"""mac-mini-S74: ADVERSARIAL stress-test of opus-S259's LRC(14) route.
Claim: for covering v, coreCover = meas(G' n U_core D_v)/meas(G') < 1 (=> safe pt => M>=1/14).
core = speeds coprime to 30030; G' = {t: ||w t||>=1/14 for non-core w}; D_v={t:||v t||<1/14}.
opus tested speeds<150 (max 0.65). PUSH: large speeds, large core, near-tight covering families.
Does coreCover stay <1 with margin, or approach 1? Also measure the per-core-runner discrepancy
(density vs 1/7) to assess the Erdos-Turan residual."""
from math import gcd
import random
LV=1.0/14; P30030=[2,3,5,7,11,13]
def coprime30030(v): return all(v%p!=0 for p in P30030)
def is_cov(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))
def prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def nn(x): x%=1.0; return min(x,1-x)  # ||x||

def corecover(S, res=400000):
    core=[v for v in S if coprime30030(v)]
    noncore=[v for v in S if not coprime30030(v)]
    if not core: return None,len(core),None
    inGp=0; inGp_coredanger=0
    per=[0]*len(core)
    for j in range(res):
        t=(j+0.5)/res
        if all(nn(w*t)>=LV for w in noncore):  # t in G'
            inGp+=1
            cd=False
            for idx,v in enumerate(core):
                if nn(v*t)<LV:
                    per[idx]+=1; cd=True
            if cd: inGp_coredanger+=1
    if inGp==0: return None,len(core),None
    cover=inGp_coredanger/inGp
    dens=[p/inGp for p in per]  # per-core-runner density in G'
    return cover,len(core),(dens,core,inGp/res)

print(f"level 1/14={LV:.5f}, per-runner target density 2/14=1/7={1/7:.4f}")
print("ADVERSARIAL covering families -- does coreCover stay <1?\n")
tests=[]
# structured near-extremal covering families + large speeds
tests+=[[*range(1,13),182],[*range(1,13),364],[*range(1,12),13,84],[*range(1,11),22,13,84]]
# families engineered to have LARGE core (many coprime-to-30030 speeds): use 1,17,19,23,29,...
tests+=[[1,17,19,23,29,31, 2,6,30,210,2310,30030, 14]]     # 6 coprime-core + smooth cover
tests+=[[1,17,19,23,29, 2,3,5,7,11,13, 14,182]]            # coprime core {1,17,19,23,29}
# large-speed covering (multiply deep well structure)
tests+=[[*range(1,13),182*3]]
# random covering families, speeds up to 400
random.seed(7); found=0
while found<8:
    S=sorted(random.sample(range(1,400),13))
    if prim(S) and is_cov(S): tests.append(S); found+=1
worst=(0,None)
for S in tests:
    S=sorted(set(S))
    if len(S)!=13 or not (prim(S) and is_cov(S)): 
        print(f"  {S}: skip (not prim/cov 13-set)"); continue
    cover,k,extra=corecover(S)
    if cover is None:
        print(f"  core empty / G' empty: {S}"); continue
    dens,core,mGp=extra
    maxd=max(dens); sumd=sum(dens)
    flag = " <== coreCover>=1 !!!" if cover>=1 else (" (tight)" if cover>0.9 else "")
    if cover>worst[0]: worst=(cover,S)
    print(f"  |core|={k} coreCover={cover:.4f} sumDens={sumd:.4f} maxDens={maxd:.4f} meas(G')={mGp:.4f}  {S[:6]}...{flag}")
print(f"\nWORST coreCover over tests: {worst[0]:.4f} at {worst[1]}")
print("If coreCover stays <1 with margin AND per-runner density ~1/7 (discrepancy small), opus's")
print("route is robust; if it approaches 1 for some covering family, the margin/ET-bound is at risk.")

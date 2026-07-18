#!/usr/bin/env python3
"""r4_feasibility_kps_S128c61.py -- kind-pasteur S128 cont.61.
SIZE the r=4 finite horn honestly: how many quadruples actually pass the necessary
condition sum_i frac_i >= 1?  If that tail is small the run is feasible; if not, say so.
Sample of cores, exact counts.  PRINT DATA ONLY."""
import sys, itertools
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
KB=400
CORES=[sorted(c) for c in itertools.combinations(range(1,13),9)]
print("### frac distribution and the size of the sum>=1 tail ###")
print("  core                         |bits| #killers  mean frac  #heavy(>=.25)  est. quads passing")
import random
random.seed(6111)
tot_est=0
for P in random.sample(CORES,10):
    M=max(P); lo=13*M+1
    bits=[i for i,(q,a) in enumerate(QS) if all(la(p*a,q)>=-(-q//14) for p in P)]
    n=len(bits)
    fr=[]
    for k in range(lo,KB):
        b=sum(1 for i in bits if la(k*QS[i][1],QS[i][0])<-(-QS[i][0]//14))
        fr.append((b/n,k))
    fr.sort(reverse=True)
    fs=[f for f,_ in fr]
    heavy=sum(1 for f in fs if f>=0.25)
    mean=sum(fs)/len(fs)
    # exact count of 4-subsets with sum >= 1 is expensive; estimate by sampling
    hit=0; TR=20000
    for _ in range(TR):
        s=sum(random.choice(fs) for _ in range(4))
        if s>=1.0: hit+=1
    import math
    tot4=math.comb(len(fs),4)
    est=tot4*hit/TR
    tot_est+=est
    print("  %-28s %-6d %-9d %-10.4f %-14d %.3e"%(str(P),n,len(fs),mean,heavy,est))
print()
print("  extrapolated over all 220 cores: ~%.2e quadruples pass the necessary condition"%(tot_est/10*220))
print("  (each needs a 4-way bitmask OR + covering test)")
print()
print("### verdict on feasibility ###")
est_total=tot_est/10*220
print("  passing quadruples ~ %.2e ; at ~3e5 tests/sec in Python that is ~%.1f hours"%(
    est_total, est_total/3e5/3600))
print("  brute force (no pruning) would be ~2e9 per the raw count, so the prune helps by ~%.0fx"%(2e9/max(est_total,1)))
print("DONE")

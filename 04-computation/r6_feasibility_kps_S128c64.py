#!/usr/bin/env python3
"""r6_feasibility_kps_S128c64.py -- kind-pasteur S128 cont.64.
Size the r=6 finite horn at KB = 333 (= max T 308.4 + 25) before running it.
A sextuple can only be uncertified if its six kill-sets COVER the core's safe (q,a) set,
which requires sum frac >= 1 -- but with SIX killers that condition is far easier to meet
than with four or five, so the prune may no longer save us.  Measure.  PRINT DATA ONLY."""
import sys, itertools, random, math
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
KB=333
C7=[sorted(c) for c in itertools.combinations(range(1,13),7)]
random.seed(64)
print("### r=6 finite horn sizing at KB = %d ###"%KB)
print("  core                        |bits| #killers  mean frac  P(sum6>=1)  est sextuples passing")
tot=0.0
sample=random.sample(C7,8)
for P in sample:
    M=max(P); lo=13*M+1
    bits=[i for i,(q,a) in enumerate(QS) if all(la(p*a,q)>=-(-q//14) for p in P)]
    n=len(bits)
    fs=[]
    for k in range(lo,KB):
        b=sum(1 for i in bits if la(k*QS[i][1],QS[i][0])<-(-QS[i][0]//14))
        fs.append(b/n)
    mean=sum(fs)/len(fs)
    hit=0; TR=20000
    for _ in range(TR):
        if sum(random.choice(fs) for _ in range(6))>=1.0: hit+=1
    p=hit/TR
    est=math.comb(len(fs),6)*p
    tot+=est
    print("  %-27s %-6d %-9d %-10.4f %-11.4f %.3e"%(str(P),n,len(fs),mean,p,est))
allest=tot/len(sample)*len(C7)
print()
print("  extrapolated over all %d cores: ~%.3e sextuples pass the necessary condition"%(len(C7),allest))
print("  at ~3e5 tests/sec that is ~%.1f hours"%(allest/3e5/3600))
print()
print("### comparison across r ###")
print("  r=4: KB=400, 1.43e8 passing, ran in ~25 min")
print("  r=5: KB=235, 2.64e8 passing, ran in ~9 min")
print("  r=6: KB=333, %.3e passing"%allest)
print("  ratio r=6 / r=5: %.1fx"%(allest/2.64e8))
print("DONE")

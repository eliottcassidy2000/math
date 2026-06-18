# LRC(14) — LEAN, FAST statistics for the parked-runner-binding question on genuine S3.
# kps-S4-wf. Smaller speeds, single Mval per set, deterministic-ish generation.
import sys
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
def pr(*a): print(*a, flush=True)
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def Mval(S):
    b=F(0); bt=None
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v; bt=t
    return b,bt
def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def is_prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def binding_set(S,t,M): return [v for v in S if nrm(v*t)==M]
def owners(S): return {q:[v for v in S if v%q==0] for q in range(2,15)}

random.seed(7)
S3=[]; tries=0
# smaller maxbig=120 to keep Mval cheap; k in {2,3}
while len(S3)<800 and tries<400000:
    tries+=1
    kbig=random.choice([2,2,3])
    ksmall=13-kbig
    P=sorted(random.sample(range(1,14),ksmall))
    L=sorted(random.sample(range(14,121),kbig))
    S=sorted(set(P+L))
    if len(S)!=13: continue
    if not is_prim(S): continue
    if not is_cov(S): continue
    if sum(1 for v in S if v>13)<2: continue
    if max(S)<13*min(S): continue
    S3.append(S)
pr(f"genuine-S3 lean sample: {len(S3)} sets (k in 2,3; big<=120)")
n_break=0; minM=F(1); minMset=None
n_uniq14=0; n_parked_binding=0
n_smallpair=0; n_anybig_binding=0
worst=[]
for S in S3:
    M,t=Mval(S)
    if M<F(1,14): n_break+=1
    if M<minM: minM=M; minMset=S
    B=binding_set(S,t,M)
    ow=owners(S)
    if len(ow[14])==1:
        n_uniq14+=1
        w=ow[14][0]
        if w in B: n_parked_binding+=1
    # is the binding set ALL small (<=13)?
    if all(v<=13 for v in B): n_smallpair+=1
    if any(v>13 for v in B): n_anybig_binding+=1
    worst.append((M,S,B,ow[14]))
worst.sort(key=lambda x:x[0])
pr(f"LRC breaks: {n_break}")
pr(f"min M = {minM} ({float(minM):.5f}; 1/14={float(F(1,14)):.5f}) on {minMset}")
pr(f"unique-mult14 sets: {n_uniq14}; parked runner BINDING at tau*: {n_parked_binding} ({100*n_parked_binding/max(1,n_uniq14):.1f}%)")
pr(f"sets whose binding runners are ALL small (<=13): {n_smallpair}/{len(S3)} ({100*n_smallpair/len(S3):.1f}%)")
pr(f"sets with SOME big (>13) binding runner: {n_anybig_binding}/{len(S3)} ({100*n_anybig_binding/len(S3):.1f}%)")
pr("\n12 hardest:")
for M,S,B,m14 in worst[:12]:
    pr(f"  M={M} ({float(M):.5f}) binding={B} mult14={m14} S={S}")

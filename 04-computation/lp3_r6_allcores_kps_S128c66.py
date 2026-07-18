#!/usr/bin/env python3
"""lp3_r6_allcores_kps_S128c66.py -- kind-pasteur S128 cont.66 (background).
Heuristic triple-moment audit on every r=6 core.  For each of the 792 seven-speed cores,
test one top-cardinality, three greedy, and fourteen seeded-random sextuple candidates.
LP3 = max sum n_d s.t. sum n_d d = S1, sum n_d C(d,2) = S2, sum n_d C(d,3) = S3, n_d >= 0.
LP3 < n certifies the particular sextuple tested, but the candidate generator is not exhaustive:
this script does not close a core or the r=6 branch.  PRINT DATA ONLY."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
KB=333; R=6
def Cc(d,k):
    if d<k: return 0
    num=1
    for t in range(k): num*=(d-t)
    for t in range(1,k+1): num//=t
    return num
DS=list(range(1,R+1))
ROWS=[[F(d) for d in DS],[F(Cc(d,2)) for d in DS],[F(Cc(d,3)) for d in DS]]
COMBOS=list(itertools.combinations(range(len(DS)),3))
def lp3(S1,S2,S3):
    rhs=[F(S1),F(S2),F(S3)]; best=None
    for combo in COMBOS:
        M=[[ROWS[i][j] for j in combo]+[rhs[i]] for i in range(3)]
        ok=True
        for c in range(3):
            piv=None
            for rr in range(c,3):
                if M[rr][c]!=0: piv=rr; break
            if piv is None: ok=False; break
            M[c],M[piv]=M[piv],M[c]
            pv=M[c][c]; M[c]=[x/pv for x in M[c]]
            for rr in range(3):
                if rr!=c and M[rr][c]!=0:
                    f=M[rr][c]; M[rr]=[M[rr][t]-f*M[c][t] for t in range(4)]
        if not ok: continue
        sol=[M[i][3] for i in range(3)]
        if any(x<0 for x in sol): continue
        t=sum(sol)
        if best is None or t>best: best=t
    return best
def mom(ms):
    S1=sum(bin(m).count("1") for m in ms)
    S2=sum(bin(ms[i]&ms[j]).count("1") for i in range(R) for j in range(i+1,R))
    S3=sum(bin(ms[i]&ms[j]&ms[k]).count("1") for i in range(R) for j in range(i+1,R) for k in range(j+1,R))
    return S1,S2,S3
random.seed(661)
C7=[sorted(c) for c in itertools.combinations(range(1,13),7)]
worst=None; sampled_nonclosing_cores=0; sampled_candidates=0
for ci,P in enumerate(C7):
    M=max(P); lo=13*M+1
    bits=[i for i,(q,a) in enumerate(QS) if all(la(p*a,q)>=-(-q//14) for p in P)]
    n=len(bits)
    if n==0 or lo>=KB: continue
    masks={}
    for k in range(lo,KB):
        km=0
        for jj,i in enumerate(bits):
            q,a=QS[i]
            if la(k*a,q)<-(-q//14): km|=(1<<jj)
        masks[k]=km
    ks=sorted(masks,key=lambda k:-bin(masks[k]).count("1"))
    cands=[ks[:R]]
    for st in ks[:3]:
        S=[st]; cur=masks[st]
        while len(S)<R:
            nx=max((k for k in ks if k not in S),key=lambda k:bin(masks[k]&~cur).count("1"))
            S.append(nx); cur|=masks[nx]
        cands.append(S)
    for _ in range(14): cands.append(random.sample(ks[:60],R))
    bestb=None
    for S in cands:
        sampled_candidates+=1
        ms=[masks[k] for k in S]
        b=lp3(*mom(ms))
        if b is None: continue
        if bestb is None or b>bestb: bestb=b
    if bestb is None: continue
    marg=float(bestb)-n
    if marg>=0: sampled_nonclosing_cores+=1
    if worst is None or marg>worst[0]: worst=(marg,tuple(P),n,float(bestb))
    if ci%150==0: print("  ... core %d/792  worst margin so far %+.1f"%(ci,worst[0]))
m,P,n,b=worst
print()
print("### r=6, all 792 cores, sampled sextuple candidates ###")
print("  sampled candidate instances: %d"%sampled_candidates)
print("  WORST sampled margin (LP3 - n) = %+.1f  at core %s (n=%d, LP3=%.1f)"%(m,list(P),n,b))
print("  cores with a sampled candidate having LP3 >= n: %d"%sampled_nonclosing_cores)
print("  cores with LP3 < n on every sampled candidate: %d"%(792-sampled_nonclosing_cores))
print("  NOTE: candidates are heuristic, not all sextuples; these counts do not close any core.")
print("DONE")

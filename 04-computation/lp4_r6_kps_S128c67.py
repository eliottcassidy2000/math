#!/usr/bin/env python3
"""lp4_r6_kps_S128c67.py -- kind-pasteur S128 cont.67 (background).
STEP 1: add S4 (quadruple overlaps) to the moment LP and see how many r=6 cores survive.
   max sum n_d  s.t.  sum n_d*C(d,j) = S_j for j = 1..4,  n_d >= 0,  d = 1..6.
Four equality constraints => basic solutions have <= 4 nonzero n_d => exact by enumerating
4-subsets of {1..6} (15 of them) and a 4x4 rational solve.  Coverage requires optimum >= n.
Report cores containing a positive sampled candidate.  Negative samples do not
certify a core because only 18 heuristic sextuples are tested.  PRINT DATA ONLY."""
import sys, itertools, random, json
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
def lpk(S,K):
    """S = [S1..SK]; max sum n_d s.t. sum n_d C(d,j) = S[j-1], j=1..K"""
    rows=[[F(Cc(d,j)) for d in DS] for j in range(1,K+1)]
    rhs=[F(x) for x in S[:K]]
    best=None
    for combo in itertools.combinations(range(len(DS)),K):
        M=[[rows[i][j] for j in combo]+[rhs[i]] for i in range(K)]
        ok=True
        for c in range(K):
            piv=None
            for rr in range(c,K):
                if M[rr][c]!=0: piv=rr; break
            if piv is None: ok=False; break
            M[c],M[piv]=M[piv],M[c]
            pv=M[c][c]; M[c]=[x/pv for x in M[c]]
            for rr in range(K):
                if rr!=c and M[rr][c]!=0:
                    f=M[rr][c]; M[rr]=[M[rr][t]-f*M[c][t] for t in range(K+1)]
        if not ok: continue
        sol=[M[i][K] for i in range(K)]
        if any(x<0 for x in sol): continue
        t=sum(sol)
        if best is None or t>best: best=t
    return best
def moms(ms,K=4):
    out=[]
    for j in range(1,K+1):
        s=0
        for combo in itertools.combinations(range(len(ms)),j):
            x=ms[combo[0]]
            for t in combo[1:]: x&=ms[t]
            s+=bin(x).count("1")
        out.append(s)
    return out
C7=[sorted(c) for c in itertools.combinations(range(1,13),7)]
print("### S4 fixed-sextuple certificate: full-core / sampled-sextuple audit ###")
pos3=[]; pos4=[]; worst3=None; worst4=None
for ci,P in enumerate(C7):
    rng=random.Random((67<<16)+ci)
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
    for _ in range(14): cands.append(rng.sample(ks[:60],R))
    b3=None; b4=None
    for S in cands:
        ms=[masks[k] for k in S]
        mm=moms(ms,4)
        v3=lpk(mm,3); v4=lpk(mm,4)
        if v3 is not None and (b3 is None or v3>b3): b3=v3
        if v4 is not None and (b4 is None or v4>b4): b4=v4
    if b3 is not None:
        m3=b3-n
        if m3>=0: pos3.append(tuple(P))
        if worst3 is None or m3>worst3[0]: worst3=(m3,tuple(P))
    if b4 is not None:
        m4=b4-n
        if m4>=0: pos4.append(tuple(P))
        if worst4 is None or m4>worst4[0]: worst4=(m4,tuple(P))
    if ci%200==0: print("  ... core %d/792  positive-sample(LP3)=%d positive-sample(LP4)=%d"%(ci,len(pos3),len(pos4)))
print()
print("  LP3: worst sampled margin %+.1f ; positive-sample cores: %d"%(float(worst3[0]),len(pos3)))
print("  LP4: worst sampled margin %+.1f ; positive-sample cores: %d"%(float(worst4[0]),len(pos4)))
print("  sampled-positive reduction LP3->LP4: %d -> %d"%(len(pos3),len(pos4)))
print("  negative sampled cores certified: 0")
json.dump([list(p) for p in pos4],open('05-knowledge/results/r6_lp4_positive_sample_cores_kps_S128c67.json','w'))
print("  positive-sample list written to 05-knowledge/results/r6_lp4_positive_sample_cores_kps_S128c67.json")
if pos4[:6]:
    print("  positive-sample cores:",[list(p) for p in pos4[:6]])
    print("  sampled structure: contain 1: %d/%d ; max=12: %d/%d"%(
        sum(1 for p in pos4 if 1 in p),len(pos4),
        sum(1 for p in pos4 if max(p)==12),len(pos4)))
print("DONE")

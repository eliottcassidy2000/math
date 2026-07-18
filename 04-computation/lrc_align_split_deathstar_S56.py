from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
random.seed(4)
def M_arg(fam):
    Q=2*max(fam)+2; best=F(0); args=[]
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            fr=F(m,q)
            if fr>best: best=fr; args=[(a,q)]
            elif fr==best: args.append((a,q))
    return best,args
def covers(fam,S): return all(any(v%q==0 for v in fam) for q in S)
def is_AP(W): W=sorted(W); d=W[0]; return W==[d*i for i in range(1,13)]
def candidates(W,M):
    delta=M-F(1,13)
    if delta<=0: return []
    ub=int(max(W)/(13*delta))
    return [k for k in range(182,ub+1,182) if k>max(W)]
# Decomposition test on valid non-AP cores: aligned (13|some D_W) vs non-aligned.
# Claim: non-aligned => all candidates give M(V)=M(W) (misalignment).  aligned => empty window (stability).
base=list(range(1,13)); pool=[v for v in range(1,36) if v%13 and v%14]
n_align_empty=0; n_align_nonempty=0; n_nonalign_covered=0; n_nonalign_MVeqMW=0; fails=[]
tested=0
def test(W):
    global n_align_empty,n_align_nonempty,n_nonalign_covered,n_nonalign_MVeqMW,tested
    W=sorted(set(W))
    if len(W)!=12 or is_AP(W) or not covers(W,range(2,13)) or any(v%13==0 or v%14==0 for v in W): return
    tested+=1
    M,args=M_arg(W); denoms=[q for a,q in args]; aligned=any(q%13==0 for q in denoms)
    cands=candidates(W,M)
    if aligned:
        if cands: n_align_nonempty+=1
        else: n_align_empty+=1
    else:
        if not cands: return
        alleq=all(M_arg(sorted(W+[vm]))[0]==M for vm in cands)
        if alleq: n_nonalign_MVeqMW+=1
        else:
            n_nonalign_covered+=1; fails.append((W,cands))
for k in (1,2):
    for pos in combinations(range(12),k):
        combos=[tuple(random.choice(pool) for _ in range(k)) for _ in range(250)] if k==2 else [(x,) for x in pool]
        for nv in combos:
            W=base[:]
            for p,x in zip(pos,nv): W[p]=x
            test(W)
print(f"valid non-AP cores tested: {tested}")
print(f"  ALIGNED (13|D_W): empty window(stability kills)={n_align_empty}, nonempty window={n_align_nonempty}")
print(f"  NON-ALIGNED: M(V)=M(W) for all candidates (misalignment kills)={n_nonalign_MVeqMW}, far COVERED (M(V)<M(W))={n_nonalign_covered}")
print(f"  cores where a candidate actually lowered M below M(W): {len(fails)}")
for W,c in fails[:5]: print("   ",W,c)

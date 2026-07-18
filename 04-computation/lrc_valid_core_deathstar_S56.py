from fractions import Fraction as F
from math import gcd
def M_arg(fam):
    Q=2*max(fam)+2; best=F(0); arg=None
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q); arg=(a,q)
    return best,arg
def covers(fam,S): return all(any(v%q==0 for v in fam) for q in S)
def is_AP(W): W=sorted(W); d=W[0]; return W==[d*i for i in range(1,13)]
# VALID non-AP compact cores: cover 2..12, miss 13,14, non-AP. Candidate far v_max=182k in stability window.
# Test: M(V)>=1/13? and record the WITNESS time (mechanism).
def candidates(W):
    (M,_)=M_arg(W); delta=M-F(1,13)
    if delta<=0: return M,delta,[]
    ub=int(max(W)/(13*delta)); lo=max(182,max(W)+1)
    return M,delta,[k for k in range(182,ub+1,182) if k>=lo]
print("VALID non-AP compact cores (cover 2..12, miss 13,14): finite check + witness:")
cores=[]
# {2..12, X}
for X in range(15,40):
    if X%13==0 or X%14==0: continue
    W=sorted(list(range(2,13))+[X])
    if len(W)==12 and covers(W,range(2,13)) and not is_AP(W): cores.append(W)
# {1..12} with one middle element bumped, keeping coverage
for rep,newv in [(6,18),(6,20),(4,16),(10,22),(8,20)]:
    W=list(range(1,13)); W[rep-1]=newv; W=sorted(set(W))
    if len(W)==12 and covers(W,range(2,13)) and not is_AP(W) and not any(v%13==0 or v%14==0 for v in W): cores.append(W)
allok=True
for W in cores[:14]:
    M,delta,cands=candidates(W)
    if not cands: 
        print(f"  W={W}: M(W)={float(M):.4f} delta={float(delta):.5f} NO candidates (window empty)"); continue
    res=[]
    for vm in cands:
        V=sorted(W+[vm]); MV,(a,q)=M_arg(V); res.append((vm,MV,q))
    ok=all(mv>=F(1,13) for _,mv,_ in res)
    if not ok: allok=False
    ws=", ".join(f"vmax={vm}:M={float(mv):.4f}@den{q}" for vm,mv,q in res)
    print(f"  W={W} (M(W)={float(M):.4f}): {ws}  all>=1/13:{ok}")
print(f"\n=> valid non-AP compact cores all eliminated (M>=1/13)? {allok}")
print("   witness denominators reveal the mechanism.")

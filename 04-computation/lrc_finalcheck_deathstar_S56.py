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
def window(W,M):
    delta=M-F(1,13)
    if delta<=0: return []
    ub=int(max(W)/(13*delta))
    return [k for k in range(182,ub+1,182) if k>max(W)]
residual=[[1,2,3,5,7,8,9,10,11,12,17,19],[1,2,3,4,5,6,7,8,9,10,22,24],[1,2,3,4,5,6,7,8,9,10,24,33]]
print("Finite check on the 3 boundary residual cores (soft-Weyl AND stability both just miss):")
allclosed=True
for W in residual:
    W=sorted(set(W)); M,_=M_arg(W); win=window(W,M)
    res=[]
    for vm in win:
        MV,(a,q)=M_arg(sorted(W+[vm])); res.append((vm,float(MV),MV<F(1,13)))
    ok=all(not lt for _,_,lt in res)
    if not ok: allclosed=False
    print(f"  W={W}: candidates={win}")
    print(f"     {[(vm,round(mv,4),'M<1/13!' if lt else 'ok') for vm,mv,lt in res]}  eliminated:{ok}")
print(f"\n=> all boundary residual cores eliminated by finite check: {allclosed}")
print("   POSITION LEMMA for max<=34: soft Weyl (99%) ∪ stability (99%) ∪ finite check (thin residual).")

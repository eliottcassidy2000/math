# death-star-S57: EXHAUSTIVE finite check of cover-gap uniqueness on the very-near-tight fragmented residual, max<=34.
# For every near-tight non-AP core (covers 2..12, misses 13,14) where BOTH lenses fail
# (soft Weyl C>464mu AND stability delta<=max/2366), verify coverGap(W,182k)>=1/13 for all far elements 182k in window.
# Also records min coverGap and the maximizer denominator D_W (alignment: 13|D_W only for the AP).
from fractions import Fraction as F
from math import gcd, floor, ceil
from itertools import combinations
import multiprocessing as mp
TH=F(1,13); NEAR=F(1,13)+F(34,2366)
PAIRS=[(1,12),(2,11),(3,10),(4,9),(5,8),(6,7)]
POOL=[v for v in range(1,35) if v%13 and v%14]
def d13ok(W):
    rs=set(v%13 for v in W)
    if 0 in rs: return True
    return all((r in rs or s in rs) for r,s in PAIRS)
def M_D(fam):
    Q=2*max(fam)+2; best=F(0); D=None
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(r,q)>best: best=F(r,q); D=q
            if best>NEAR: return None,None
    return best,D
def good_components(W):
    ivs=[]
    for w in W:
        for k in range(0,w+1):
            lo=(F(k)-TH)/w; hi=(F(k)+TH)/w
            a=lo if lo>0 else F(0); b=hi if hi<1 else F(1)
            if b>a: ivs.append((a,b))
    ivs.sort(); mg=[]
    for a,b in ivs:
        if mg and a<=mg[-1][1]:
            if b>mg[-1][1]: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    comps=[]; prev=F(0)
    for a,b in mg:
        if a>prev: comps.append((prev,a))
        if b>prev: prev=b
    if prev<1: comps.append((prev,F(1)))
    return comps
def ndist(x):
    f=x-floor(x); return f if f<=F(1,2) else 1-f
def max_norm_on(a,b,vmax):
    va=vmax*a; vb=vmax*b
    if vb-va>=1: return F(1,2)
    m=max(ndist(va),ndist(vb))
    if ceil(va-F(1,2))+F(1,2)<=vb: return F(1,2)
    return m
def coverGap(comps,vmax):
    if not comps: return F(0)
    return max(max_norm_on(a,b,vmax) for a,b in comps)
def isAP(W):
    W=sorted(W); d=W[0]; return d>0 and all(W[i]==d*(i+1) for i in range(12))
def work(first):
    rest=[v for v in POOL if v>first]; res=[]; fails=[]; near=0
    for tail in combinations(rest,11):
        W=(first,)+tail
        if not d13ok(W): continue
        if not all(any(v%q==0 for v in W) for q in range(2,13)): continue
        M,D=M_D(W)
        if M is None: continue
        near+=1
        if isAP(W): continue
        delta=M-TH; comps=good_components(W); C=len(comps); mu=sum(b-a for a,b in comps)
        if C<=F(4643,10)*mu or delta>F(max(W),2366): continue   # covered by a lens
        # RESIDUAL: finite check over window
        ub=int(max(W)/(13*delta)) if delta>0 else 0
        mincg=F(1)
        for vm in range(((max(W)//182)+1)*182, ub+1, 182):
            g=coverGap(comps,vm); mincg=min(mincg,g)
            if g<TH: fails.append((tuple(W),vm))
        res.append((tuple(W),str(M),C,float(mu),float(delta),float(mincg),D))
    return res,fails,near
if __name__=='__main__':
    firsts=[f for f in POOL if f<=23]
    with mp.Pool(4) as p:
        out=p.map(work,firsts)
    RES=[]; FAILS=[]; near=0
    for r,f,n in out: RES+=r; FAILS+=f; near+=n
    print("near-tight non-AP cores (max<=34): %d ; RESIDUAL (both lenses fail): %d"%(near,len(RES)),flush=True)
    print("coverGap<1/13 FAILURES in residual: %d"%len(FAILS),flush=True)
    RES.sort(key=lambda r:r[5])
    for r in RES: print("  %-36s M=%-6s C=%2d mu=%.5f delta=%.5f minCoverGap=%.4f D_W=%d (13|D_W?%s)"%(str(list(r[0])),r[1],r[2],r[3],r[4],r[5],r[6],r[6]%13==0),flush=True)
    if not FAILS:
        m=min((r[5] for r in RES),default=1.0)
        print("=> ALL residual cores have coverGap>=1/13 (min over all = %.4f, 1/13=%.4f). COMPACT residual CLOSED by finite check."%(m,float(TH)),flush=True)
    else:
        for w in FAILS[:20]: print("  *** FAIL:",list(w[0]),"vmax=",w[1],flush=True)

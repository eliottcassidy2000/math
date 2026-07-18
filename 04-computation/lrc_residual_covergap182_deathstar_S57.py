# death-star-S57: characterize the very-near-tight FRAGMENTED residual and its coverGap(182).
# Residual = non-AP core (covers 2..12, misses 13,14) with C>464mu AND delta<=max/2366 (both lenses fail).
# Question: is coverGap(W,182) barely >=1/13 (hard kernel) or =1/2 (far element too fine -> tractable)?
from fractions import Fraction as F
from math import gcd, floor, ceil
from itertools import combinations
TH=F(1,13); NEAR=F(1,13)+F(34,2366)
def M_val(fam):
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
def maxlen(comps): return max((b-a for a,b in comps),default=F(0))
def isAP(W):
    W=sorted(W); d=W[0]; return d>0 and all(W[i]==d*(i+1) for i in range(12))
pool=[v for v in range(1,25) if v%13 and v%14]   # max<=24 for speed
residual=[]; nearcnt=0
for W in combinations(pool,12):
    if not all(any(v%q==0 for v in W) for q in range(2,13)): continue
    M,D=M_val(W)
    if M is None: continue
    nearcnt+=1
    if isAP(W): continue
    delta=M-TH; comps=good_components(W); C=len(comps); mu=sum(b-a for a,b in comps)
    softweyl = C<=F(4643,10)*mu       # C<=464.3 mu
    stability = delta>F(max(W),2366)
    if softweyl or stability: continue   # covered by a rigorous lens
    # RESIDUAL: both lenses fail. Compute coverGap(182) and window min.
    cg182=coverGap(comps,182)
    # min coverGap over window 182k
    ub=int(max(W)/(13*delta)) if delta>0 else 0
    mincg=cg182;
    for k in range(1, ub//182+2):
        vm=182*k
        if vm>max(W) and vm<=ub:
            g=coverGap(comps,vm); mincg=min(mincg,g)
    residual.append((tuple(W),str(M),C,float(mu),float(delta),float(cg182),float(mincg),float(maxlen(comps)),D))
print("near-tight non-AP cores (max<=24):",nearcnt-1 if nearcnt else 0,"  RESIDUAL (both lenses fail):",len(residual))
print("%-34s %-6s %2s %8s %8s %9s %9s %9s %3s"%("W","M","C","mu","delta","cGap182","minGapWin","maxlen","D"))
for r in residual:
    print("%-34s %-6s %2d %8.5f %8.5f %9.4f %9.4f %9.5f %3d"%(str(list(r[0])),r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8]))
if residual:
    m=min(r[6] for r in residual)
    print("\nMIN coverGap over window across all residual cores: %.5f (1/13=%.5f)"%(m,float(TH)))
    print("=> if all >=1/13 with margin, residual is tractable; 1/2 means far element too fine to cover.")
else:
    print("=> RESIDUAL EMPTY at max<=24: soft Weyl + stability cover ALL non-AP cores.")

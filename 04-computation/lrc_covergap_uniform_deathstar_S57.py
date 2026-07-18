# death-star-S57: verify the structure of the uniform cover-gap bound before proving it.
# For near-tight non-AP cores W: compute mu, C, vmax_min, min coverGap over the far-element window.
# Check (i) claim: coverGap>=1/13 for ALL valid far elements; (ii) Prop B soft-Weyl split C <= 2.551*vmax*mu.
from fractions import Fraction as F
from math import gcd, floor, ceil
from itertools import combinations
TH=F(1,13); NEAR=F(1,13)+F(34,2366)
def M_and_args(fam):
    Q=2*max(fam)+2; best=F(0); args=[]
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam); fr=F(r,q)
            if fr>NEAR: return None,None
            if fr>best: best=fr; args=[(a,q)]
            elif fr==best: args.append((a,q))
    return best,args
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
    m=ndist(va); d=ndist(vb)
    if d>m: m=d
    if ceil(va-F(1,2))+F(1,2)<=vb: return F(1,2)
    return m
def coverGap(comps,vmax):
    best=F(0)
    for (a,b) in comps:
        v=max_norm_on(a,b,vmax)
        if v>best: best=v
    return best
def is_dilated_AP(W):
    W=sorted(W); d=W[0]; return all(W[i]==d*(i+1) for i in range(12))
# enumerate near-tight non-AP cores (perturbations of the AP and dilated APs), check the two claims
bad=[]; frag=[]; nonfrag=0; total=0; minmargin=(F(1),None)
import random
seeds=[[i for i in range(1,13)]]+[[2*i for i in range(1,13)]]+[[3*i for i in range(1,12)]+[x] for x in range(34,35)]
# systematic: replace 1-2 elements of {1..12} with values up to 34
base=list(range(1,13))
cands=[]
for r in range(1,3):
    for pos in combinations(range(12),r):
        for repl in combinations([x for x in range(1,35) if x%13 and x%14 and x not in base],r):
            W=base[:]
            for idx,p in enumerate(pos): W[p]=repl[idx]
            if len(set(W))<12: continue
            cands.append(tuple(sorted(W)))
cands=list(set(cands))
for W in cands:
    if not all(any(v%q==0 for v in W) for q in range(2,13)): continue   # covers 2..12
    M,args=M_and_args(W)
    if M is None: continue
    if is_dilated_AP(W): continue
    delta=M-TH
    if delta<=0: continue
    total+=1
    comps=good_components(W); C=len(comps); mu=sum(b-a for a,b in comps)
    miss=[q for q in range(2,14) if not any(v%q==0 for v in W)]; L=1
    for q in miss: L=L*q//gcd(L,q)
    ub=int(max(W)/(13*delta)); start=((max(W)//L)+1)*L
    if start>ub: continue
    # min coverGap over ALL valid far elements
    mincg=F(1); argcg=None
    for vmax in range(start,ub+1,L):
        cg=coverGap(comps,vmax)
        if cg<mincg: mincg=cg; argcg=vmax
    if mincg<minmargin[0]: minmargin=(mincg,W,argcg)
    if mincg<TH: bad.append((W,argcg,float(mincg)))     # claim violation = counterexample!
    # Prop B split at the binding (smallest) far element
    vm=start
    if C<=F(2551,1000)*vm*mu: nonfrag+=1
    else: frag.append((W,C,float(mu),vm,float(coverGap(comps,vm))))
print("near-tight non-AP cores tested: %d"%total,flush=True)
print("CLAIM violations (min coverGap < 1/13 = counterexamples): %d"%len(bad),flush=True)
if bad:
    for w in bad[:10]: print("  *** VIOLATION",w)
print("min coverGap over all cores/far-elements: %.5f (1/13=%.5f) at %s"%(float(minmargin[0]),float(TH),str(minmargin[1]) if minmargin[1] else None),flush=True)
print("Prop B (C<=2.551*vmax*mu, non-fragmented, soft Weyl closes): %d/%d"%(nonfrag,total),flush=True)
print("FRAGMENTED (C>2.551*vmax*mu, reduce to displacement kernel): %d"%len(frag),flush=True)
for w in frag[:12]: print("  frag: W=%s C=%d mu=%.4f vmax=%d coverGap=%.4f"%(list(w[0]),w[1],w[2],w[3],w[4]),flush=True)

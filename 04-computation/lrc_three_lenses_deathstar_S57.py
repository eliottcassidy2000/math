# death-star-S57: do THREE rigorous lenses close cover-gap uniqueness for ALL near-tight non-AP cores (max<=34)?
#  L1 soft Weyl:  C <= 464.3 mu           => coverGap(182) >= avg >= 1/13
#  L2 stability:  delta > max/2366        => component half-width > far-arc => coverGap >= 1/13
#  L3 maximizer:  some maximizer t0=p/D_W has ||182 t0|| >= 1/13  => coverGap >= ||182 t0|| >= 1/13  (t0 in G_W)
# Report cores where ALL THREE fail. If none => uniform proof of the compact far-element inverse theorem.
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
def maximizers(fam):
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
    C=0; prev=F(0); mu=F(0)
    for a,b in mg:
        if a>prev: C+=1; mu+=a-prev
        if b>prev: prev=b
    if prev<1: C+=1; mu+=1-prev
    return C,mu
def nd(x):
    f=x-floor(x); return f if f<=F(1,2) else 1-f
def isAP(W):
    W=sorted(W); d=W[0]; return d>0 and all(W[i]==d*(i+1) for i in range(12))
def work(first):
    rest=[v for v in POOL if v>first]; near=0; allfail=[]
    for tail in combinations(rest,11):
        W=(first,)+tail
        if not d13ok(W): continue
        if not all(any(v%q==0 for v in W) for q in range(2,13)): continue
        M,args=maximizers(W)
        if M is None: continue
        near+=1
        if isAP(W): continue
        delta=M-TH
        # L3 maximizer (cheapest, try first)
        if max(nd(182*F(a,q)) for a,q in args)>=TH: continue
        # L2 stability
        if delta>F(max(W),2366): continue
        # L1 soft Weyl
        C,mu=good_components(W)
        if C<=F(4643,10)*mu: continue
        allfail.append((tuple(W),str(M),C,float(mu),float(delta)))
    return near,allfail
if __name__=='__main__':
    firsts=[f for f in POOL if f<=23]
    with mp.Pool(4) as p:
        out=p.map(work,firsts)
    near=0; AF=[]
    for n,af in out: near+=n; AF+=af
    print("near-tight non-AP cores (max<=34): %d"%near,flush=True)
    print("cores where ALL THREE lenses fail: %d"%len(AF),flush=True)
    for r in AF[:30]: print("  ALL-FAIL:",list(r[0]),"M=%s C=%d mu=%.5f delta=%.5f"%(r[1],r[2],r[3],r[4]),flush=True)
    if not AF:
        print("=> UNIFORM: soft Weyl OR stability OR maximizer closes EVERY near-tight non-AP core (max<=34).",flush=True)
        print("=> The far-element inverse theorem (compact) is PROVED by three rigorous lenses -- no finite check.",flush=True)

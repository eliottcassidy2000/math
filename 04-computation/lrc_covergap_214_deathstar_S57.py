# death-star-S57: CORRECT cover-gap enumeration at COVERING 2..14 (far element 182|v_max).
# Class: primitive covering-2..14 13-families with M<1/13, AP-core-adjacent (core misses 13,14).
# Core W = 12-subset of {1..MAXW}\{mults of 13,14}, covers 2..12, near-tight (M<=NEAR), non-AP.
# Far element vmax = 182k (covers 13 AND 14, so V covers 2..14), max(W)<vmax<=max(W)/(13 delta) [stability].
# coverGap(W,vmax)=max_{G_W}||vmax t|| : M(V)<1/13 iff coverGap<1/13 (far element covers G_W). Exact.
# Confirms the CORRECT compact-case closure (retracts the wrong 'mults of 13' cont22 run, MISTAKE-161).
from fractions import Fraction as F
from math import gcd, floor, ceil
from itertools import combinations
import multiprocessing as mp, sys
MAXW=int(sys.argv[1]) if len(sys.argv)>1 else 34
TH=F(1,13); NEAR=F(1,13)+F(MAXW,2366)
POOL=[v for v in range(1,MAXW+1) if v%13 and v%14]
def M_and_args(fam):
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam); fr=F(r,q)
            if fr>NEAR: return None
            if fr>best: best=fr
    return best
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
def far_covers(comps,vmax):
    for (a,b) in comps:
        if max_norm_on(a,b,vmax)>=TH: return False
    return True
def is_dilated_AP(W):
    W=sorted(W); d=W[0]; return all(W[i]==d*(i+1) for i in range(12))
def work(first):
    counter=[]; near=0; checks=0
    rest=[v for v in POOL if v>first]
    for tail in combinations(rest,11):
        W=(first,)+tail
        if not all(any(v%q==0 for v in W) for q in range(2,13)): continue   # covers 2..12
        M=M_and_args(W)
        if M is None: continue
        near+=1
        if is_dilated_AP(W): continue
        delta=M-TH
        if delta<=0: continue
        ub=int(max(W)/(13*delta)); start=((max(W)//182)+1)*182     # far element = 182k > max(W)
        if start>ub: continue
        comps=good_components(W)
        for vmax in range(start,ub+1,182):
            checks+=1
            if far_covers(comps,vmax):
                counter.append((tuple(W),vmax,str(M)))
    return counter,near,checks
if __name__=='__main__':
    firsts=[f for f in POOL if f<=MAXW-11]
    with mp.Pool(4) as p:
        res=p.map(work,firsts)
    C=[]; near=0; checks=0
    for c,n,ch in res: C+=c; near+=n; checks+=ch
    print("MAXW=%d: near-tight non-AP cores (cover 2..12, miss 13,14)=%d; far-element(182k) cover-gap checks=%d"%(MAXW,near,checks))
    print("COUNTEREXAMPLES (coverGap<1/13 => M(V)<1/13, non-AP core, covering 2..14): %d"%len(C))
    if not C:
        print("=> NO counterexample. COVERING-2..14 compact case (max<=%d) CLOSED: every near-tight non-AP core has coverGap>=1/13 for all 182k."%MAXW)
    else:
        for w in C[:20]: print("  *** COUNTEREXAMPLE:",list(w[0]),"vmax=",w[1],"M=",w[2])

# death-star-S56: COMPLETE closure of the compact case max<=34, PARALLEL (4 workers).
# Same WLOG reductions as lrc_residual_complete; adds a cheap necessary pre-filter D13 and splits by min element.
#  D13 gate: M(W) >= (max_{a=1..12} min_w dist(w*a mod 13))/13. If W has NO multiple of 13, near-tight (M<=NEAR<2/13)
#   FORCES min_w dist(w*a,13)<=1 for every a, i.e. residues mod 13 hit all 6 antipodal pairs {1,12}..{6,7}.
#   (If W has a multiple of 13, D13=0 and the gate passes trivially.)  Necessary => safe pre-filter.
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import multiprocessing as mp
NEAR=F(1,13)+F(34,2366)
def cmask(v):
    m=0
    for i,q in enumerate(range(2,11)):
        if v%q==0: m|=(1<<i)
    return m
MASK=[cmask(v) for v in range(35)]; FULL=(1<<9)-1
PAIRS=[(1,12),(2,11),(3,10),(4,9),(5,8),(6,7)]
def M_and_args(fam):
    Q=2*max(fam)+2; best=F(0); args=[]
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            fr=F(r,q)
            if fr>NEAR: return None,None
            if fr>best: best=fr; args=[(a,q)]
            elif fr==best: args.append((a,q))
    return best,args
def M_exact(fam):
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(r,q)>best: best=F(r,q)
    return best
def is_dilated_AP(W):
    W=sorted(W); d=W[0]
    return all(W[i]==d*(i+1) for i in range(12))
def safe_at(vmax,args):
    for (a,q) in args:
        r=(vmax*a)%q
        if 13*min(r,q-r)>=q: return True
    return False
def work(k):
    counter=[]; near=0; checks=0; scanned=0
    for tail in combinations(range(k+1,35),11):
        W=(k,)+tail; scanned+=1
        res=[v%13 for v in W]                     # D13 gate
        if 0 not in res:
            rs=set(res)
            if not all((r in rs or (13-r) in rs) for (r,_) in PAIRS): continue
        om=0                                       # covering 2..10
        for v in W: om|=MASK[v]
        if om!=FULL: continue
        M,args=M_and_args(W)                       # near-tight?
        if M is None: continue
        near+=1
        if is_dilated_AP(W): continue
        delta=M-F(1,13)
        if delta<=0: continue
        miss=[q for q in range(2,14) if not any(v%q==0 for v in W)]
        L=1
        for q in miss: L=L*q//gcd(L,q)
        ub=int(max(W)/(13*delta)); start=((max(W)//L)+1)*L
        for vmax in range(start,ub+1,L):
            checks+=1
            if not safe_at(vmax,args) and M_exact(list(W)+[vmax])<F(1,13):
                counter.append((tuple(W),vmax))
    return counter,near,checks,scanned
if __name__=='__main__':
    with mp.Pool(4) as p:
        results=p.map(work, range(1,24))
    C=[]; near=0; checks=0; scanned=0
    for c,n,ch,sc in results:
        C+=c; near+=n; checks+=ch; scanned+=sc
    print("DONE: %d subsets scanned across 4 workers."%scanned,flush=True)
    print("near-tight non-AP covering cores found; total candidate (W,vmax) far-element checks = %d"%checks,flush=True)
    print("COUNTEREXAMPLES (M(V)<1/13, non-AP core): %d"%len(C),flush=True)
    if not C:
        print("=> NO counterexample. Every near-tight non-AP core with max<=34 has M(V)>=1/13 for ALL far elements.",flush=True)
        print("=> COMPACT CASE max<=34 CLOSED: boxeph inverse theorem holds for all cores with max<=34.",flush=True)
    else:
        for w in C[:30]: print("  *** COUNTEREXAMPLE:",list(w[0]),"vmax=",w[1],flush=True)

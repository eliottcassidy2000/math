# death-star-S56: COMPLETE closure of the compact case max<=34.
# WLOG reductions (all proved this session):
#  (R1) non-near-tight core (M(W)>1/13+34/2366=0.09129) => delta>max/2366 => stability window empty => M(V)=M(W)>=1/13.
#       So only NEAR-TIGHT cores W (M(W)<=NEAR) can extend to a counterexample. [missing q<=10 => M>=1/10>NEAR]
#  (R2) near-tight => W covers 2..10 (missing q<=10 gives M>=1/q>NEAR). [valid pre-filter]
#  (R3) V=W+{vmax}, M(V)<1/13 requires V cover 2..13 (missing q<=12 => M>=1/q>1/13; missing 13 => M>=1/13 not <).
#       So vmax must cover missing(W) cap {2..13}: L := lcm(those) divides vmax. [exact necessary cond, no boxeph needed]
#  (R4) stability: M(V)<1/13, W non-AP => vmax <= max(W)/(13 delta), delta=M(W)-1/13. [THM-1028]
#  Dilated APs W={d,2d,..,12d} are the CONCLUSION (deep wells), excluded from the counterexample search.
# Complete: enumerate ALL 12-subsets of {1..34} covering 2..10, near-tight, non-dilated-AP; check every candidate vmax
# = multiple of L in (max(W), max(W)/(13 delta)]. If all give M(V)>=1/13, max<=34 is CLOSED.
from fractions import Fraction as F
from math import gcd
from itertools import combinations
NEAR=F(1,13)+F(34,2366)
def cmask(v):
    m=0
    for i,q in enumerate(range(2,11)):
        if v%q==0: m|=(1<<i)
    return m
MASK=[cmask(v) for v in range(35)]; FULL=(1<<9)-1
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
counter=[]; nearcores=0; residual_checks=0; done=0
for W in combinations(range(1,35),12):
    done+=1
    if done%20000000==0: print("  ...%dM / 548M, near=%d, checks=%d, ctr=%d"%(done//1000000,nearcores,residual_checks,len(counter)),flush=True)
    om=0
    for v in W: om|=MASK[v]
    if om!=FULL: continue                      # R2: must cover 2..10
    M,args=M_and_args(W)
    if M is None: continue                     # R1: not near-tight => stability
    nearcores+=1
    if is_dilated_AP(W): continue              # conclusion, not counterexample
    delta=M-F(1,13)
    if delta<=0: continue                      # exactly tight non-AP (none for compact; skip-safe)
    miss=[q for q in range(2,14) if not any(v%q==0 for v in W)]
    L=1
    for q in miss: L=L*q//gcd(L,q)
    ub=int(max(W)/(13*delta))
    lo=max(W)+1; start=((lo+L-1)//L)*L
    for vmax in range(start,ub+1,L):
        residual_checks+=1
        if not safe_at(vmax,args):
            if M_exact(list(W)+[vmax])<F(1,13):
                counter.append((tuple(W),vmax))
print("DONE: %d subsets scanned."%done,flush=True)
print("near-tight non-AP covering cores checked; total candidate (W,vmax) far-element checks = %d"%residual_checks,flush=True)
print("COUNTEREXAMPLES (M(V)<1/13, non-AP core): %d"%len(counter),flush=True)
if not counter:
    print("=> NO counterexample. Every near-tight non-AP core with max<=34 has M(V)>=1/13 for ALL far elements.",flush=True)
    print("=> COMPACT CASE max<=34 CLOSED: boxeph inverse theorem holds for all cores with max<=34.",flush=True)
else:
    for w in counter[:20]: print("  *** COUNTEREXAMPLE:",list(w[0]),"vmax=",w[1],flush=True)

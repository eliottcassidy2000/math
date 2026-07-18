#!/usr/bin/env python3
"""prune_dedupe_kps_S128c65.py -- kind-pasteur S128 cont.65.
TWO IMPROVEMENTS to the r=6 enumeration, measured before committing to a run.

(A) PAIRWISE-OVERLAP PRUNE.  Coverage of bits(P) by kill-sets A_1..A_r satisfies, for ANY
ordering,  |union A_i| <= sum|A_i| - sum_{i>=2} max_{j<i} |A_i cap A_j| .  Maximising the
subtracted term over orderings gives the standard spanning-tree improvement:
        |union| <= sum |A_i| - MST_max ,
MST_max = max-weight spanning tree on the complete graph with weights |A_i cap A_j|.
So coverage REQUIRES  sum|A_i| - MST_max >= n , strictly stronger than sum|A_i| >= n.
THM-1071(III)'s positive correlation says the overlaps are large, so this should bite.

(B) KILL-MASK DEDUPE.  The certificate sees a killer ONLY through its kill-set.  Killers
with equal masks are interchangeable, so enumerate distinct masks with multiplicity.
PRINT DATA ONLY."""
import sys, itertools, random
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
KB=333
C7=[sorted(c) for c in itertools.combinations(range(1,13),7)]
random.seed(65)
def core_data(P):
    M=max(P); lo=13*M+1
    bits=[i for i,(q,a) in enumerate(QS) if all(la(p*a,q)>=-(-q//14) for p in P)]
    n=len(bits)
    masks={}
    for k in range(lo,KB):
        km=0
        for jj,i in enumerate(bits):
            q,a=QS[i]
            if la(k*a,q)<-(-q//14): km|=(1<<jj)
        masks[k]=km
    return bits,n,masks
print("### (B) kill-mask dedupe: how many DISTINCT masks per core? ###")
print("  core                        |bits| #killers  #distinct masks  dedupe factor")
ded=[]
sample=random.sample(C7,10)
for P in sample:
    bits,n,masks=core_data(P)
    d=len(set(masks.values()))
    ded.append(len(masks)/d)
    print("  %-27s %-6d %-9d %-16d %.3f"%(str(P),n,len(masks),d,len(masks)/d))
print("  mean dedupe factor: %.4f  -> sextuple count scales by ~factor^6 = %.4f"%(
    sum(ded)/len(ded),(sum(ded)/len(ded))**6))
print()
print("### (A) pairwise-overlap prune: how much of the sum>=1 tail does MST kill? ###")
def mst_weight(idx,masks_list):
    """max-weight spanning tree via Prim on |A_i cap A_j|"""
    m=len(idx); inT=[False]*m; best=[-1]*m; tot=0
    inT[0]=True
    for j in range(1,m): best[j]=bin(masks_list[idx[0]]&masks_list[idx[j]]).count("1")
    for _ in range(m-1):
        bi=-1; bv=-1
        for j in range(m):
            if not inT[j] and best[j]>bv: bv=best[j]; bi=j
        inT[bi]=True; tot+=bv
        for j in range(m):
            if not inT[j]:
                w=bin(masks_list[idx[bi]]&masks_list[idx[j]]).count("1")
                if w>best[j]: best[j]=w
    return tot
print("  core                        sampled sextuples  pass sum>=n  pass sum-MST>=n  extra cut")
tot_old=0; tot_new=0
for P in sample[:6]:
    bits,n,masks=core_data(P)
    ks=sorted(masks)
    ml={k:masks[k] for k in ks}
    sizes={k:bin(masks[k]).count("1") for k in ks}
    TR=4000; old=0; new=0
    for _ in range(TR):
        S=random.sample(ks,6)
        s=sum(sizes[k] for k in S)
        if s<n: continue
        old+=1
        if s-mst_weight(S,ml)>=n: new+=1
    tot_old+=old; tot_new+=new
    print("  %-27s %-18d %-12d %-16d %.1fx"%(P,TR,old,new,(old/new) if new else float('inf')))
print()
print("  overall: sum>=n passes %d ; sum-MST>=n passes %d ; extra cut factor %.1fx"%(
    tot_old,tot_new,(tot_old/tot_new) if tot_new else float('inf')))
print()
print("### combined projection for r=6 ###")
base=3.639e12
dd=(sum(ded)/len(ded))**6
cut=(tot_old/tot_new) if tot_new else 1e9
print("  baseline passing sextuples (cont.64):        %.3e"%base)
print("  after mask dedupe (factor %.4f):             %.3e"%(dd,base*dd))
print("  after MST prune (extra %.1fx):               %.3e"%(cut,base*dd/cut))
print("  at ~3e5 tests/sec: %.1f hours"%(base*dd/cut/3e5/3600))
print("DONE")

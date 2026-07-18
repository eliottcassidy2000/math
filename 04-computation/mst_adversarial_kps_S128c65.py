#!/usr/bin/env python3
"""mst_adversarial_kps_S128c65.py -- kind-pasteur S128 cont.65.
THE DECISIVE TEST.  Random sextuples all failed sum - MST >= n.  But the cases that matter
are ADVERSARIAL: the killers with the largest kill-sets.  If
        max over sextuples of ( sum|A_i| - MST_max )  <  n
for every core, then NO sextuple can cover bits(P), so no r=6 family is uncertified and the
whole enumeration is unnecessary.
Search that max by (i) the top-6 by size, (ii) greedy max-marginal, (iii) local search.
The bound  |union A_i| <= sum|A_i| - MST_max  is valid for any sets (each A_i beyond the
first loses at least its best overlap with a predecessor; maximising over orderings is
exactly the max spanning tree).  PRINT DATA ONLY."""
import sys, itertools, random
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
KB=333
C7=[sorted(c) for c in itertools.combinations(range(1,13),7)]
def mstw(ms):
    m=len(ms); inT=[False]*m; best=[-1]*m; tot=0; inT[0]=True
    for j in range(1,m): best[j]=bin(ms[0]&ms[j]).count("1")
    for _ in range(m-1):
        bi=-1; bv=-1
        for j in range(m):
            if not inT[j] and best[j]>bv: bv=best[j]; bi=j
        inT[bi]=True; tot+=bv
        for j in range(m):
            if not inT[j]:
                w=bin(ms[bi]&ms[j]).count("1")
                if w>best[j]: best[j]=w
    return tot
def score(ks,masks):
    ms=[masks[k] for k in ks]
    return sum(bin(m).count("1") for m in ms)-mstw(ms)
random.seed(651)
print("### max over sextuples of (sum|A_i| - MST) vs n = |bits(P)| ###")
print("  core                        n     best score found   union of best   margin  covers?")
worst=None; anycover=False
for P in C7:
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
    ks=sorted(masks,key=lambda k:-bin(masks[k]).count("1"))
    best=None
    cand=[tuple(ks[:6])]
    # greedy: maximise marginal new coverage
    for start in ks[:6]:
        S=[start]; cur=masks[start]
        while len(S)<6:
            nx=max((k for k in ks if k not in S), key=lambda k: bin(masks[k]&~cur).count("1"))
            S.append(nx); cur|=masks[nx]
        cand.append(tuple(S))
    # random restarts + local search on the score
    for _ in range(60):
        S=list(random.sample(ks[:60],6))
        improved=True
        while improved:
            improved=False
            base=score(S,masks)
            for pos in range(6):
                for repl in ks[:80]:
                    if repl in S: continue
                    T=list(S); T[pos]=repl
                    sc=score(T,masks)
                    if sc>base: S=T; base=sc; improved=True; break
                if improved: break
        cand.append(tuple(S))
    for S in cand:
        sc=score(list(S),masks)
        if best is None or sc>best[0]: best=(sc,S)
    sc,S=best
    un=0
    for k in S: un|=masks[k]
    uc=bin(un).count("1")
    if uc==n: anycover=True
    if worst is None or (sc-n)>worst[0]: worst=(sc-n,tuple(P),n,sc,S)
print("  worst margin (score - n) over all %d cores: %d"%(len(C7),worst[0]))
print("    at core %s : n = %d, best score = %d, sextuple %s"%(list(worst[1]),worst[2],worst[3],worst[4]))
print("  any sextuple actually covering bits(P): %s"%anycover)
print()
if worst[0]<0:
    print("  => sum|A_i| - MST < n for EVERY core and every sextuple searched.")
    print("     Coverage is impossible, so NO r=6 family is uncertified -- no enumeration needed.")
else:
    print("  => some sextuple passes the MST condition; enumeration still required there.")
print("DONE")

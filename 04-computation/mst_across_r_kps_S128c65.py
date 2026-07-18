#!/usr/bin/env python3
"""mst_across_r_kps_S128c65.py -- kind-pasteur S128 cont.65.
Is coverage IMPOSSIBLE for structural reasons?  Test max( sum|A_i| - MST ) vs n across
r = 4, 5, 6.  If the bound falls short of n at every r, it explains why the r=4 and r=5
enumerations found zero uncertified families, and it would replace enumeration entirely.
Leaner search: top-by-size, greedy max-marginal, and a short local search.  PRINT DATA ONLY."""
import sys, itertools, random
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
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
random.seed(652)
print("### max( sum|A_i| - MST ) vs n, across r ###")
print("   r  cores  KB    worst margin (score-n)   max union / n      any cover?")
for r,size,KB in [(4,9,400),(5,8,235),(6,7,333)]:
    CS=[sorted(c) for c in itertools.combinations(range(1,13),size)]
    CS=random.sample(CS,min(40,len(CS)))
    wm=None; bestuc=0; bestn=1; cover=False
    for P in CS:
        M=max(P); lo=13*M+1
        bits=[i for i,(q,a) in enumerate(QS) if all(la(p*a,q)>=-(-q//14) for p in P)]
        n=len(bits)
        if n==0 or lo>=KB: continue
        masks={}
        for k in range(lo,KB):
            km=0
            for jj,i in enumerate(bits):
                q,a=QS[i]
                if la(k*a,q)<-(-q//14): km|=(1<<jj)
            masks[k]=km
        ks=sorted(masks,key=lambda k:-bin(masks[k]).count("1"))
        cands=[ks[:r]]
        for st in ks[:4]:
            S=[st]; cur=masks[st]
            while len(S)<r:
                nx=max((k for k in ks if k not in S),key=lambda k:bin(masks[k]&~cur).count("1"))
                S.append(nx); cur|=masks[nx]
            cands.append(S)
        for _ in range(12):
            cands.append(random.sample(ks[:50],r))
        for S in cands:
            ms=[masks[k] for k in S]
            sc=sum(bin(m).count("1") for m in ms)-mstw(ms)
            un=0
            for m in ms: un|=m
            uc=bin(un).count("1")
            if uc==n: cover=True
            if uc/n>bestuc/max(bestn,1): bestuc,bestn=uc,n
            if wm is None or sc-n>wm[0]: wm=(sc-n,tuple(P),n,sc)
    print("   %d  %-6d %-5d %-24d %-18s %s"%(r,len(CS),KB,wm[0],"%d/%d = %.3f"%(bestuc,bestn,bestuc/bestn),cover))
print()
print("  interpretation: score - n < 0 everywhere means sum|A_i| - MST never reaches n,")
print("  so the kill-sets provably cannot cover bits(P) -- coverage impossible, no enumeration needed.")
print("  max union / n shows how far the BEST sextuple actually gets toward covering.")
print("DONE")

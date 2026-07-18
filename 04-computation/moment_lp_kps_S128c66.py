#!/usr/bin/env python3
"""moment_lp_kps_S128c66.py -- kind-pasteur S128 cont.66.
THE TRIPLE-OVERLAP CORRECTION, as a MOMENT LP.
Let n_d = #points of bits(P) lying in exactly d of the kill-sets.  Then
   |union A_i| = sum_d n_d ,   S1 = sum_d n_d * d  = sum|A_i| ,
   S2 = sum_d n_d * C(d,2) = sum_{i<j} |A_i cap A_j| ,
   S3 = sum_d n_d * C(d,3) = sum_{i<j<k} |A_i cap A_j cap A_k| .
So an UPPER bound on the union is the LP
        maximise sum_d n_d   s.t.  the three moment equations,  n_d >= 0, d = 1..r.
Coverage requires that optimum >= n.  With 3 equality constraints the optimum sits at a
basic solution with <= 3 nonzero n_d, so it can be solved EXACTLY by enumerating triples of
d-values -- no LP solver needed.
Compare against the MST bound on the same heuristic candidate families.  The cores and
sextuples/quadruples below are seeded samples, so a negative displayed margin is a certificate
only for that displayed candidate, not a uniform r-branch theorem.  PRINT DATA ONLY."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
def C(d,k):
    if d<k: return 0
    num=1
    for t in range(k): num*=(d-t)
    for t in range(1,k+1): num//=t
    return num
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
def moments(ms):
    r=len(ms)
    S1=sum(bin(m).count("1") for m in ms)
    S2=sum(bin(ms[i]&ms[j]).count("1") for i in range(r) for j in range(i+1,r))
    S3=sum(bin(ms[i]&ms[j]&ms[k]).count("1") for i in range(r) for j in range(i+1,r) for k in range(j+1,r))
    return S1,S2,S3
def lp_bound(S1,S2,S3,r,use3=True):
    """max sum n_d s.t. moment equations, n_d>=0 ; exact via basic solutions"""
    ds=list(range(1,r+1))
    best=None
    rows=[[F(d) for d in ds],[F(C(d,2)) for d in ds]]
    rhs=[F(S1),F(S2)]
    if use3:
        rows.append([F(C(d,3)) for d in ds]); rhs.append(F(S3))
    k=len(rows)
    for combo in itertools.combinations(range(len(ds)),k):
        # solve k x k system
        A=[[rows[i][j] for j in combo] for i in range(k)]
        b=list(rhs)
        M=[row[:]+[b[i]] for i,row in enumerate(A)]
        ok=True
        for c in range(k):
            piv=None
            for rr in range(c,k):
                if M[rr][c]!=0: piv=rr; break
            if piv is None: ok=False; break
            M[c],M[piv]=M[piv],M[c]
            pv=M[c][c]
            M[c]=[x/pv for x in M[c]]
            for rr in range(k):
                if rr!=c and M[rr][c]!=0:
                    f=M[rr][c]
                    M[rr]=[M[rr][t]-f*M[c][t] for t in range(k+1)]
        if not ok: continue
        sol=[M[i][k] for i in range(k)]
        if any(x<0 for x in sol): continue
        tot=sum(sol)
        if best is None or tot>best: best=tot
    return best
random.seed(66)
print("### moment LP vs MST, on sampled sextuple/quadruple candidates ###")
print("   r  core                     n    |union|  MST bound  LP(S1,S2)  LP(S1,S2,S3)  rejects candidate?")
for r,size,KB in [(4,9,400),(5,8,235),(6,7,333)]:
    CS=[sorted(c) for c in itertools.combinations(range(1,13),size)]
    CS=random.sample(CS,14)
    worstMST=None; worstLP=None
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
        for st in ks[:3]:
            S=[st]; cur=masks[st]
            while len(S)<r:
                nx=max((k for k in ks if k not in S),key=lambda k:bin(masks[k]&~cur).count("1"))
                S.append(nx); cur|=masks[nx]
            cands.append(S)
        for _ in range(10): cands.append(random.sample(ks[:50],r))
        for S in cands:
            ms=[masks[k] for k in S]
            S1,S2,S3=moments(ms)
            mst=S1-mstw(ms)
            lp2=lp_bound(S1,S2,0,r,use3=False)
            lp3=lp_bound(S1,S2,S3,r,use3=True)
            un=0
            for m in ms: un|=m
            uc=bin(un).count("1")
            if worstMST is None or mst-n>worstMST[0]: worstMST=(mst-n,)
            if worstLP is None or (lp3 is not None and float(lp3)-n>worstLP[0]):
                worstLP=(float(lp3)-n,tuple(P),n,uc,mst,float(lp2) if lp2 else -1,float(lp3))
    m,P,n,uc,mst,lp2,lp3=worstLP
    print("   %d  %-24s %-4d %-8d %-11d %-10.1f %-13.1f %s"%(
        r,str(list(P)),n,uc,mst,lp2,lp3,"YES" if lp3<n else "no"))
    print("      worst margins:  MST %+d   LP2 %+.1f   LP3 %+.1f"%(worstMST[0],lp2-n,lp3-n))
print()
print("  margin < 0 certifies noncoverage for that sampled candidate only.")
print("  NOTE: the candidate generator is heuristic; this output proves no uniform r-branch closure.")
print("DONE")

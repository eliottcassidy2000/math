"""Leaner covering-min IP + M_exact verification of the found set."""
import math
from fractions import Fraction
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix
def M_exact(S,Qmax):
    best=Fraction(0)
    for q in range(2,Qmax+1):
        bb=0
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            m=q
            for s in S:
                r=(s*a)%q; d=r if r<=q-r else q-r
                if d<m:m=d
                if m<=bb:break
            if m>bb:bb=m
        v=Fraction(bb,q)
        if v>best:best=v
    return best
def cov_min_IP(n,Smax,Qmax):
    cands=list(range(1,Smax+1)); idx={s:i for i,s in enumerate(cands)}; N=len(cands)
    witnesses=[(q,a) for q in range(2,Qmax+1) for a in range(1,q) if math.gcd(a,q)==1]
    vs=sorted({Fraction(m,q) for q in range(n,Qmax+1) for m in range(1,q//(n-1)+2) if Fraction(1,n)<=Fraction(m,q)<Fraction(1,n-1)})
    covrows=[[idx[s] for s in cands if s%q==0] for q in range(2,n+1)]
    c=np.ones(N)
    for v in vs:
        nr=[]
        ok=True
        for (q,a) in witnesses:
            thr=float(v)*q+1e-9
            near=[idx[s] for s in cands if min((s*a)%q,q-(s*a)%q)<=thr]
            if not near: ok=False; break
            nr.append(near)
        if not ok: continue
        rows=covrows+nr; data=[];ri=[];ci=[]
        for r,row in enumerate(rows):
            for j in row: data.append(1);ri.append(r);ci.append(j)
        A=csr_matrix((data,(ri,ci)),shape=(len(rows),N))
        res=milp(c=c,constraints=[LinearConstraint(A,np.ones(len(rows)),np.full(len(rows),np.inf))],integrality=np.ones(N),bounds=Bounds(0,1))
        if res.success and round(res.fun)<=n-1:
            S=sorted(cands[i] for i in range(N) if res.x[i]>0.5)
            return v,S
    return None,None
for n,Sm,Qm in [(7,18,22),(9,40,40)]:
    v,S=cov_min_IP(n,Sm,Qm)
    if v is None: print(f"n={n}: no feasible v found in range (raise bounds)"); continue
    Mv=M_exact(S,3*(n*n-n+1))
    print(f"n={n}: IP covmin candidate v={v}={float(v):.5f}, set={S}; VERIFY M_exact={Mv}={float(Mv):.5f}; match={v==Mv}")

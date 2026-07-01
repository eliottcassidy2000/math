#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S72 -- LAZY-CUT (cutting-plane) prover for the covering-min at the CONSTRUCTION scale.
Extends the exact set-cover ILP (HYP-3731/3733/3778) to V=n(n-1), where the full ILP times out, by adding
danger-arc-covering cuts lazily until infeasibility.

Question: is there a PRIMITIVE COVERING (n-1)-set S with all speeds <= V=n(n-1) and M(S) < r := n/Phi6 ?
(a STRICT beater of the construction {1..n-2,n(n-1)} which has M=n/Phi6 exactly.)
Loop: solve feasibility [cardinality + covering + primitivity + accumulated cuts]; if a candidate S appears,
compute M(S) EXACTLY; if M(S)<r it's a BEATER (overturn); else take its M-witness t* (min_v||v t*||=M(S)>=r)
and add the VALID cut "a strict beater must strictly-danger t*": sum_{v: ||v t*||<r} x_v >= 1.  Any strict
beater satisfies every cut (its open danger arcs cover all t), so INFEASIBLE => no strict beater =>
covering-min(speeds<=V) = n/Phi6, RIGOROUS.
"""
import sys, time, numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix

def Mexact(S):
    Sg=sorted(set(S)); cand=set()
    for i in range(len(Sg)):
        for j in range(len(Sg)):
            for d in (Sg[i]-Sg[j],Sg[i]+Sg[j]):
                if d>0:
                    for k in range(1,d): cand.add(F(k,d))
    best=F(0); at=None
    for t in cand:
        g=min(min((v*t)%1,1-((v*t)%1)) for v in Sg)
        if g>best: best=g; at=t
    return best, at
def primes_upto(V): return [p for p in range(2,V+1) if all(p%d for d in range(2,int(p**.5)+1))]

def prove(n, maxit=9000, tbudget=1500, verbose=True):
    V=n*(n-1); r=F(n, n*n-n+1); k=n-1
    cols=list(range(1,V+1)); idx={v:i for i,v in enumerate(cols)}
    rowidx=[]; lbs=[]; ubs=[]     # each row = list of column-indices (sparse)
    def addrow(vs, lb, ub):
        rowidx.append([idx[v] for v in vs]); lbs.append(lb); ubs.append(ub)
    addrow(cols, k, k)
    for q in range(2,n+1): addrow([v for v in cols if v%q==0], 1, np.inf)
    for p in primes_upto(V): addrow([v for v in cols if v%p!=0], 1, np.inf)
    integ=np.ones(V); bnds=Bounds(np.zeros(V), np.ones(V)); c=np.zeros(V)
    cuts=0; t0=time.time()
    for it in range(maxit):
        # build csr from indptr (O(nnz), no dense)
        indptr=[0]; indices=[]
        for rr in rowidx:
            indices.extend(rr); indptr.append(len(indices))
        data=np.ones(len(indices))
        A=csr_matrix((data, np.array(indices), np.array(indptr)), shape=(len(rowidx), V))
        lc=LinearConstraint(A, np.array(lbs,dtype=float), np.array(ubs,dtype=float))
        res=milp(c=c, constraints=[lc], integrality=integ, bounds=bnds)
        if not res.success:
            return ("INFEASIBLE", cuts, r)
        S=[cols[i] for i in range(V) if res.x[i]>0.5]
        M,tstar=Mexact(S)
        if M<r:
            return ("BEATER", S, M, r)
        a,d=tstar.numerator, tstar.denominator
        dang=[v for v in cols if min((v*a)%d, d-((v*a)%d)) < r*d]
        addrow(dang, 1, np.inf); cuts+=1
        if verbose and cuts%40==0: print(f"    n={n}: {cuts} cuts, last M={float(M):.5f} (target {float(r):.5f}), {time.time()-t0:.0f}s", flush=True)
        if time.time()-t0>tbudget: return ("TIMEOUT", cuts, r)
    return ("MAXIT", cuts, r)

if __name__=="__main__":
    for n in ([int(sys.argv[1])] if len(sys.argv)>1 else [12,13,14]):
        r=F(n,n*n-n+1)
        print(f"\n=== n={n}: target covering-min r = n/Phi6 = {r} = {float(r):.6f}, V=n(n-1)={n*(n-1)} ===", flush=True)
        out=prove(n)
        if out[0]=="INFEASIBLE":
            print(f"  INFEASIBLE after {out[1]} cuts => NO covering (n-1)-set with speeds<={n*(n-1)} has M<{r}")
            print(f"  => covering-min(n={n}, speeds<=n(n-1)) = {r} RIGOROUS (construction optimal at construction scale).")
        elif out[0]=="BEATER":
            print(f"  !!!! BEATER FOUND: S={sorted(out[1])} M={out[2]}={float(out[2]):.6f} < {out[3]} -- OVERTURN")
        else:
            print(f"  {out[0]} ({out[1]} cuts) -- inconclusive within budget")

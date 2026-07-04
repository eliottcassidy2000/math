#!/usr/bin/env python3
"""
The covering-min as a RATIONAL min-max with EQUIOSCILLATION (mac-mini-2026-07-04-S40).
Inspired by arXiv:1612.00337 (AAA rational approximation: greedy/adaptive/barycentric min-max).
The LRC view M(S)=max_t min_i ||v_i t|| is a max-min; its optimizer t* is a RATIONAL (CF/Ostrowski point),
and at t* the binding runners EQUIOSCILLATE (>=2 runners at exactly M). Tests:
 (1) for covering families: t*, its continued fraction, the binding (equioscillating) runner set.
 (2) GREEDY Stern-Brocot descent for t* (AAA-greedy analog): start coarse, refine by mediant toward max-min.
 (3) equioscillation => a de la Vallee-Poussin-style LOWER bound: >=2 binding runners at t*=a/q* pin M=|.|/q*.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
def gcd_all(xs): return reduce(gcd,xs)
def nd(x):
    x=x%1; return x if x<=1-x else 1-x
def cf(fr):
    a=[]; x=fr
    for _ in range(10):
        if x==0: break
        r=1/x; ai=r.numerator//r.denominator; a.append(ai); x=r-ai
    return a
def M_and_t(sp):
    vs=sorted(set(sp)); Q=set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1,len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best=F(0); bt=None
    for q in Q:
        if q<2: continue
        for a in range(1,q):
            t=F(a,q); m=min(nd(v*t) for v in sp)
            if m>best: best,bt=m,t
    return best,bt
def binding(sp,t,M):
    return [v for v in sp if nd(v*t)==M]
def is_cov(sp,n): return all(any(v%q==0 for v in sp) for q in range(2,n+1))
def greedy_sternbrocot(sp, depth=40):
    """AAA-greedy analog: descend the Stern-Brocot tree, at each node pick the child (mediant side) with
       larger min_i ||v_i t||; return the best t found (a rational hiding spot)."""
    lo=(F(0),F(1)); hi=(F(1),F(1))  # not used; do interval [0,1] mediant search over Farey
    a,b,c,d=0,1,1,1  # interval (a/b, c/d) starts (0/1,1/1)
    best=F(0); bestt=None
    for _ in range(depth):
        med=F(a+c,b+d)
        for t in [med]:
            m=min(nd(v*t) for v in sp)
            if m>best: best,bestt=m,t
        # pick side with larger min at its mediant-ish probe
        lm=min(nd(v*F(a+(a+c),b+(b+d))) for v in sp) if b+d>0 else 0
        rm=min(nd(v*F((a+c)+c,(b+d)+d)) for v in sp)
        if lm>=rm: c,d=a+c,b+d
        else: a,b=a+c,b+d
    return best,bestt

if __name__=="__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    n=14
    fams=[("deep well",list(range(1,13))+[182]),("2..14",list(range(2,15))),
          ("cov A",[2,3,4,5,6,7,8,9,10,11,12,13,14]),("cov B",[1,2,4,6,7,8,9,10,11,12,13,14,182])]
    out("(1)+(3) covering families: t*, CF, EQUIOSCILLATION (binding runners), greedy vs exact:")
    out(f"{'family':>12} {'M':>9} {'t*':>10} {'CF(t*)':>14} {'#binding':>8} {'binding runners':>22}")
    for name,S in fams:
        if not is_cov(S,n):
            out(f"  {name}: not covering, skip"); continue
        M,t=M_and_t(S)
        bind=binding(S,t,M)
        out(f"{name:>12} {str(M):>9} {str(t):>10} {str(cf(t)):>14} {len(bind):>8} {str(sorted(bind)):>22}")
    out("\n(2) GREEDY Stern-Brocot descent (AAA-greedy analog) vs exact M:")
    rng=random.Random(40); ok=0; tot=0
    for name,S in fams:
        if not is_cov(S,n): continue
        Mex,_=M_and_t(S); Mg,tg=greedy_sternbrocot(S)
        tot+=1
        hit = (Mg==Mex)
        if hit: ok+=1
        out(f"   {name:>12}: exact M={float(Mex):.5f}  greedy M={float(Mg):.5f} at t={tg}  {'HIT' if hit else 'miss (local max)'}")
    # random covering families: does greedy reach >= 14/183?
    below=0; cnt=0
    for _ in range(3000):
        hi=rng.choice([20,40,80,200]); S=sorted(set(rng.sample(range(1,hi),13)))
        if len(S)!=13 or gcd_all(S)!=1 or not is_cov(S,14): continue
        cnt+=1; Mg,_=greedy_sternbrocot(S)
        if Mg<F(14,183): below+=1
    out(f"\n   greedy on {cnt} covering families: reached >=14/183 in {cnt-below}, below in {below}")
    out("=> equioscillation (>=2 binding runners at t*=a/q*) is the de la Vallee-Poussin lower-bound structure;")
    out("   greedy Stern-Brocot (AAA-style) finds t* when it doesn't stick in a local max => the certificate is")
    out("   a RATIONAL t* whose CF convergents are the binding denominators.")

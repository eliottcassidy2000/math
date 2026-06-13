#!/usr/bin/env python3
"""The MASTER OBJECT: covering-depth distribution p_k = meas{t: depth(t)=k}, depth =
#runners within delta of origin. Lonely = {depth=0}, measure p_0. (A) singleton-wall
exponent: p_0(delta') as delta'->collapse (predict exponent 1). (B) moments. (C) the
p_0=0 collapse family = AP + additive chains (1,3,4,7),(1,3,4,5,9). opus-2026-06-03-S599.
Convention: n runners, gap delta=1/(n+1)."""
from fractions import Fraction as F
from math import log
def dist(x): x%=1; return min(x,1-x)
def depth_dist(V,delta):
    # partition [0,1) by breakpoints ||v t||=delta i.e. t=(k +- delta)/v; depth const on each cell
    eps=set([F(0)])
    for v in V:
        for k in range(v+1):
            for s in(1,-1): eps.add((F(k)+s*delta)/v % 1)
    pts=sorted(eps); from_collections=None
    p={}; L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        d=sum(1 for v in V if dist(v*mid)<delta)
        p[d]=p.get(d,F(0))+ln
    return p
def p0(V,delta): return depth_dist(V,delta).get(0,F(0))
def main():
    print("(B) MASTER OBJECT: depth distribution p_k, moments (sum=1, E[depth]=2n*delta)")
    for V,n in [((1,2,3,4),4),((1,3,4,7),4),((1,2,3,4,5),5),((1,3,4,5,9),5)]:
        d=F(1,n+1); p=depth_dist(V,d)
        tot=sum(p.values()); Ed=sum(k*v for k,v in p.items())
        print(f"  V={V} (n={n}, delta=1/{n+1}): p_k={ {k:str(v) for k,v in sorted(p.items())} }; sum={tot}; E[depth]={Ed}={float(Ed):.3f} (2n*delta={float(2*n*d):.3f}); p_0={p.get(0,0)} COLLAPSE={p.get(0,F(0))==0}")
    print()
    print("(C) the p_0=0 COLLAPSE family (additive chains, top=sum of two below):")
    for V,n,why in [((1,2,3,4,5),5,"AP"),((1,3,4,7),4,"4=1+3,7=3+4"),((1,3,4,5,9),5,"4=1+3,5=1+4,9=4+5"),((1,2,4,7),4,"non-chain control")]:
        d=F(1,n+1); pp=p0(V,d)
        print(f"  V={V} ({why}): p_0={pp} collapse={pp==0}")
    print()
    print("(A) SINGLETON-WALL EXPONENT: p_0(delta') as delta'->1/(n+1) from below; fit alpha (p_0~eps^alpha)")
    for V,n in [((1,2,3,4),4),((1,2,3,4,5),5),((1,3,4,7),4)]:
        dc=F(1,n+1); 
        vals=[]
        for ek in [F(1,50),F(1,100),F(1,200),F(1,400)]:
            dp=dc-ek; pp=p0(V,dp); vals.append((float(ek),float(pp)))
        # exponent between consecutive
        exps=[]
        for i in range(len(vals)-1):
            (e1,p1),(e2,p2)=vals[i],vals[i+1]
            if p1>0 and p2>0: exps.append(log(p1/p2)/log(e1/e2))
        print(f"  V={V} (n={n}): p_0(eps)={[('%.4f'%p) for e,p in vals]} at eps={[e for e,p in vals]}; exponents alpha~{[round(x,2) for x in exps]}")
if __name__=='__main__': main()

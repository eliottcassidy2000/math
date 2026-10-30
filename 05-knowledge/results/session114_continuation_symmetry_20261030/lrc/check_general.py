#!/usr/bin/env python3
"""Exact audits of GENERAL_CORE.md; stdlib only, explicit exception gates."""
from fractions import Fraction as F
from math import gcd
import json
from check_resonant_tail import survivor, maximum, distance

checks=0

def check(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(str(label))

def least_unit(q):return next(u for u in range(2,q) if gcd(u,q)==1)

def core(q):return tuple(range(1,q))+(q+1,)

def windows(q,r):
    eps=F(1,q)-r; out=[]
    for a in range(1,q):
        if gcd(a,q)!=1:continue
        p=pow(a,-1,q);n=q-p
        if p==1:p=q+1
        if n==1:n=q+1
        out.append((F(a,q)-eps/p,F(a,q)+eps/n))
    return out

def terminal_B(V,M,S):
    B=F(0)
    for a,b in S:
        P=[v for v in V if (v*a)%1==M]
        N=[v for v in V if (v*b)%1==1-M]
        check(P and N,('terminal owners',V,M,a,b))
        B+=F(1,max(P))+F(1,max(N))
    return B

def explicit_witness(q,m):
    s=least_unit(q)
    if 2*m<=s:return F(1,q+2)
    a=pow(s,-1,q);h=(a*s-1)//q
    return F(a*m+h,q*m+s)

def main():
    visits=0; rows=0; general=[]
    for q in range(3,41):
        W=core(q);s=least_unit(q);r0=F(1,q+2);rho=F(1,q)
        for lam in (F(1,100),F(1,3),F(3,4),F(99,100),F(1)):
            r=r0+lam*(rho-r0)
            check(survivor(W,r)==windows(q,r),('complete top chamber',q,r))
        # Prime character of least missing unit, checked directly.
        check(all(s%d for d in range(2,s)),('least unit prime',q,s))
        if q<=24:
            for m in range(1,9):
                V=W+(q*m,);want=max(r0,F(m,q*m+s))
                got,points,n=maximum(V); visits+=n;rows+=1
                check(got==want,('global spectrum',q,m,s,got,want))
                if 2*m>s:
                    t=explicit_witness(q,m)
                    check(points==sorted([t,1-t]),('all strict maximizers',q,m,points,t))
                    S=survivor(V,want)
                    check(S==[(t,t) for t in points],('literal maximizing survivor',q,m))
                    check(terminal_B(V,want,S)==2*(F(1,s)+F(1,q*m)),('strict derivative',q,m))
                else:
                    t=F(1,q+2)
                    check(min(distance(v,t) for v in V)==want,('floor witness',q,m))
            for v in sorted(set((q+2,q+3,2*q-1,2*q+1))):
                if v%q:
                    got,points,n=maximum(W+(v,));visits+=n;rows+=1
                    check(got==rho,('nonmultiple spectrum',q,v))
                    check(points==[F(a,q) for a in range(1,q) if gcd(a,q)==1],('nonmultiple maximizers',q,v))
    # Large structural controls, using explicit witnesses and exact integer phases.
    for q,m in ((30,3),(30,4),(210,5),(210,6),(2310,6),(2310,7),(30030,8),(30030,9),(12,10**30+7)):
        s=least_unit(q);rho=max(F(1,q+2),F(m,q*m+s));t=explicit_witness(q,m)
        a,b=t.numerator,t.denominator
        mn=min(min(v*a%b,b-v*a%b) for v in core(q)+(q*m,))
        check(F(mn,b)==rho,('large witness',q,m,s))
        if 2*m>s:check(q*m*(F(1,q)-rho)/s==rho,('large upper equality',q,m))
        general.append({'q':q,'m':m,'s':s,'M':str(rho),'witness':str(t),'tight':rho==F(1,q+2)})
    print(json.dumps({'status':'passed','checks':checks,'core_window_q':'3..40 at five exact radii','global_maxima':'q3..24,m1..8 plus four nonmultiple tails','global_rows':rows,'candidate_clock_visits':visits,'large_structural_controls':general,'proof':'GENERAL_CORE.md','unrestricted_LRC':'OPEN'},indent=2,sort_keys=True))

if __name__=='__main__':main()

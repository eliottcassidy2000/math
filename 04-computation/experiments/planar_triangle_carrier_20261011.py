#!/usr/bin/env python3
"""Exact algebra/ramification controls, not a search proving planar JC."""
from fractions import Fraction as F
from math import gcd
import json
import sympy as s

checks=0


def check(ok,message):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(message)


def signature(a,b,c):
    d=gcd(gcd(a,b),c)
    a,b,c=a//d,b//d,c//d
    w,k=2*a+3*b+6*c,a+b+3*c
    g0,g1,gi=gcd(w,k),gcd(w,c),gcd(w,k-c)
    numerator=w+2-g0-g1-gi
    check(numerator%2 == 0,"integral genus")
    genus=numerator//2
    orders=(w//g0,w//g1,w//gi)
    chi=sum((F(1,e) for e in orders),F(-1))
    check(2-2*genus == w*chi,"orbifold Riemann-Hurwitz")
    check(gcd(gcd(w,k),c)==1,"connected primitive Kummer cover")
    check(genus>=2 and chi<0,"hyperbolic invariant")
    return {"primitive":[a,b,c],"power":d,"degree":w,"orders":orders,"genus":genus}


def main():
    signatures={}
    check([F(1,2)+F(1,3)+F(1,m)-1 for m in (5,6,7)]==[F(1,30),F(0),F(-1,42)],
          "spherical Euclidean hyperbolic triangle boundary")
    for a in range(1,25):
        for b in range(25):
            for c in range(1,13):
                info=signature(a,b,c)
                if b==0 and c==1:
                    check(info['genus']==(a+1)//2+1,"inherited hyperelliptic genus")
    p,y,r,z,tau,x,t=s.symbols('p y r z tau x t')
    D=p**3-y**2
    def bracket(f,g):
        return s.expand(-D/p*(s.diff(f,p)*s.diff(g,y)-s.diff(f,y)*s.diff(g,p)))
    small=0
    for a in range(2,6):
        for b in range(4):
            for c in range(1,4):
                I=p**a*y**b*D**c
                w,k=2*a+3*b+6*c,a+b+3*c
                check(s.cancel(I.subs({p:r**2/z,y:r**3/z})-
                               r**w*(1-z)**c/z**k)==0,"birational Kummer identity")
                dp=bracket(p,I)
                dy=bracket(y,I)
                check(s.expand(s.diff(I,p)*dp+s.diff(I,y)*dy)==0,"invariant Hamiltonian")
                ds=s.expand((dy*p-y*dp)/p**2)
                dtau=s.expand(dp-2*(y/p)*ds)
                # Independent chart differentiation: I=tau^c r^b(r^2+tau)^(a+b+2c).
                J=tau**c*r**b*(r*r+tau)**(a+b+2*c)
                check(s.cancel(ds.subs({p:r*r+tau,y:r*(r*r+tau)})-
                               tau*s.diff(J,tau))==0,"chart s coefficient")
                check(s.cancel(dtau.subs({p:r*r+tau,y:r*(r*r+tau)})+
                               tau*s.diff(J,r))==0,"chart tau coefficient")
                leading=s.Poly(s.expand(dp.subs({p:r*r+tau,y:r*(r*r+tau)})),tau).coeff_monomial(tau**c)
                check(s.expand(leading-2*c*r**(w-2*c+1))==0,"nonzero leading displacement")
                small+=1
    for abc in [(2,0,1),(3,0,1),(2,1,1),(2,0,2),(2,0,3),(4,2,2),(1,0,1)]:
        signatures[str(abc)]=signature(*abc)
    # Literal monodromy permutations give a second finite genus calculation.
    for a in range(1,9):
        for b in range(9):
            for c in range(1,6):
                if gcd(gcd(a,b),c) != 1:
                    continue
                w,k=2*a+3*b+6*c,a+b+3*c
                cycle_counts=[]
                for shift in (k,-c,c-k):
                    unseen=set(range(w))
                    cycles=0
                    while unseen:
                        v=next(iter(unseen))
                        while v in unseen:
                            unseen.remove(v)
                            v=(v+shift)%w
                        cycles+=1
                    cycle_counts.append(cycles)
                genus=signature(a,b,c)['genus']
                check(sum(cycle_counts)-w==2-2*genus,"monodromy Euler characteristic")
    level=s.symbols('level')
    univariate=[]
    for H in [p**2, p**3, p**2*(p-1),p**2*(p-1)**2,
              p**2*(p-1)*(p+2),p**4*(p*p+1)**3]:
        h=s.degree(H,p)
        factors=s.factor_list(H,p)[1]
        odd=sum(s.degree(factor,p) for factor,m in factors if m%2)
        numerator=s.Poly(p**3*H-level,p,domain=s.QQ.frac_field(level))
        check(s.gcd(numerator,numerator.diff()).degree()==0,"generic squarefree numerator")
        check(s.gcd(numerator,s.Poly(H,p,domain=s.QQ.frac_field(level))).degree()==0,
              "no cancellation at poles")
        genus=(int(h)+int(odd)+2)//2
        check((h+odd)%2==0 and genus>=2,"univariate genus")
        displacement=bracket(p,H*D).subs({p:r*r+tau,y:r*(r*r+tau)})
        leading=s.Poly(s.expand(displacement),tau).coeff_monomial(tau)
        check(s.expand(leading-2*r**5*H.subs(p,r*r))==0,"univariate leading displacement")
        univariate.append({"H":str(H),"degree":int(h),"odd_poles":int(odd),"genus":genus})
    # Same primitive curve, misleading larger raw Newton triangle.
    check(signature(2,0,2)['genus']==signature(1,0,1)['genus'],"strip invariant power")
    # Full weighted-homogeneous factor patterns with 1..3 distinct factors.
    weighted=0
    for a in range(2,9):
        for b in range(6):
            for mults in [(1,), (2,), (1,1),(1,2),(2,2),(1,1,1),(1,2,3)]:
                d=gcd(a,b)
                for c in mults:
                    d=gcd(d,c)
                aa,bb=a//d,b//d
                cs=[c//d for c in mults]
                total=sum(cs)
                w,k=2*aa+3*bb+6*total,aa+bb+3*total
                n=len(cs)*w-gcd(w,k)-gcd(w,k-total)-sum(gcd(w,c) for c in cs)
                check(n>0 and n%2==0,"weighted homogeneous genus bound")
                if len(cs)==1:
                    check(n//2+1==signature(a,b,mults[0])['genus'],"triangle specialization")
                weighted+=1
    commutator=bracket(p*p*D,p*p*y*D)
    check(s.expand(commutator+p*p*D*D*(5*p**3-2*y*y))==0,"mixed leading commutator")
    check(2*15-2 == 2*22-11-2-2-1,"commutator genus fifteen")
    # Genus-zero non-LND flow is a required hostile.
    lam=s.symbols('lam')
    X=x/(1-lam*x)
    T=t*(1-lam*x)**2
    check(s.cancel(s.diff(X,x)*s.diff(T,t)-s.diff(X,t)*s.diff(T,x)-1)==0,"rational hostile volume")
    check(s.cancel(X**2*T-x*x*t)==0,"rational hostile invariant")
    print(json.dumps({"status":"FINITE-EXACT; universal genus/nonrationality in companion proof",
                      "checks":checks,"ramification_universe":[[1,24],[0,24],[1,12]],
                      "symbolic_carriers":small,"examples":signatures,
                      "univariate_controls":univariate,
                      "weighted_factor_patterns":weighted,"leading_commutator_genus":15},indent=2,sort_keys=True))


if __name__=='__main__':
    main()

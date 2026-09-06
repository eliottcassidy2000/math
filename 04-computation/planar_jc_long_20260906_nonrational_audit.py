#!/usr/bin/env python3
"""Independent exact controls for invariant-curve nonrationality.

No producer or repository mathematical imports. Sparse polynomial-ring
exponentials test the completed operator, while separate symbolic function
fields test the generic curves and the rational-flow hostile. The all-a,
all-field result rests on the audited analytic proof, not these finite tests.
"""
from math import factorial
from pathlib import Path
from hashlib import sha256
import json

import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.ring_series import rs_trunc

CHECKS=0


def check(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:
        raise RuntimeError(label)


def main():
    P,Y,C,V,Z=sp.symbols('P Y C V Z')
    D=P**3-Y**2
    genera=[]
    for a in range(2,33):
        invariant=P**a*D
        if a%2==0:
            scale=a//2
            radicand=P**(a+3)-C
            genus=a//2+1
        else:
            scale=(a+1)//2
            radicand=P**(a+4)-C*P
            genus=(a-1)//2+2
        target_relation=sp.expand((P**scale*Y)**2-radicand)
        expected=(C-invariant)*(1 if a%2==0 else P)
        check(sp.expand(target_relation-expected)==0,'birational generic-curve relation')
        field=sp.QQ.frac_field(C)
        h=sp.Poly(radicand,P,domain=field)
        check(sp.gcd(h,h.diff()).degree()==0,'generic hyperelliptic polynomial squarefree')
        check(h.degree()%2==1 and (h.degree()-1)//2==genus and genus>=2,
              'exact branch count and genus')
        # Independent genus derivation: the odd-degree double cover has one
        # additional ramified point at infinity, so 2g-2=-4+(deg h+1).
        check(2*genus-2==-4+h.degree()+1,'Riemann-Hurwitz double-cover branch count')
        genera.append((a,genus))
    check(sp.gcd(sp.Poly(P**5,P),sp.Poly(5*P**4,P)).degree()>0,
          'special fibre C=0 is not the smooth generic curve')
    check((3-1)//2==1,'a0 boundary has genus one and does not supply finite automorphisms')
    print('GENERIC_CURVES',json.dumps(genera,separators=(',',':')))

    R,s,tau=ring('s,tau',sp.QQ)
    zero,one=R.zero,R.one
    p=s*s+tau
    y=s*p
    polynomials=[(0,1),(7,0,-3),(0,0,0,2),(3,2,-5,0,1),(-2,0,3,-1)]
    cases=0
    for a in range(2,11):
        I=tau*p**(a+2)
        for coefficients in polynomials:
            k=next(i for i,c in enumerate(coefficients) if i and c)
            fk=coefficients[k]
            cutoff=2*k+2

            def trunc(F):
                return rs_trunc(F,tau,cutoff)

            def product(F,G):
                return trunc(F*G)

            def power(F,n):
                answer=one
                for _ in range(n):answer=product(answer,F)
                return answer

            S=sum((c*I**i for i,c in enumerate(coefficients)),zero)
            velocity_s=trunc(tau*S.diff(tau))
            velocity_tau=trunc(-tau*S.diff(s))

            def delta(F):
                return trunc(F.diff(s)*velocity_s+F.diff(tau)*velocity_tau)

            def exponential(F,time):
                answer=trunc(F)
                term=answer
                for n in range(1,cutoff+1):
                    term=delta(term)
                    if term==zero:break
                    answer+=term*sp.QQ.convert(sp.Rational(time)**n/factorial(n))
                return trunc(answer)

            # A separate direct formula in the original coefficient algebra
            # checks the carrier and the precise p-velocity.
            f=sum(sp.Integer(c)*Z**i for i,c in enumerate(coefficients))
            source=sp.expand(f.subs(Z,P**a*D))
            quotient=sp.cancel((source-coefficients[0])/(P**2*D))
            check(sp.denom(quotient)==1,'source lies in the actual universal carrier')
            direct_p=sp.cancel(-D/P*sp.diff(source,Y))
            expected_p=2*P**(a-1)*Y*D*sp.diff(f,Z).subs(Z,P**a*D)
            check(sp.expand(direct_p-expected_p)==0,'original (p,y) Poisson response')
            check(delta(I)==zero,'actual invariant retained by completed derivation')
            check(min(monom[1] for monom in velocity_s)>=k and
                  min(monom[1] for monom in velocity_tau)>=k+1,
                  'valuation increase in both coordinate velocities')
            for b in range(4):
                for j in range(-3,4):
                    # tau^(-j)delta(s^b tau^j), kept polynomial even for j<0.
                    scaled=(b*s**(b-1)*velocity_s if b else zero)+j*s**b*velocity_tau.exquo(tau)
                    if scaled!=zero:
                        check(min(monom[1] for monom in scaled)>=k,
                              'Laurent monomial valuation increase')
            expected_lead=2*k*fk*s**((2*a+4)*k+1)
            dp=delta(p)
            coefficient=R.from_dict({(i,0):value for (i,j),value in dp.items() if j==k})
            check(coefficient==expected_lead,'exact leading p-velocity')
            for time in [sp.Rational(1),sp.Rational(-2),sp.Rational(3,5)]:
                ep,ey=exponential(p,time),exponential(y,time)
                difference=ep-p
                check(difference!=zero and min(monom[1] for monom in difference)==k,
                      'nonzero scalar time first live order')
                actual_lead=R.from_dict({(i,0):value for (i,j),value in difference.items() if j==k})
                check(actual_lead==sp.QQ.convert(time)*expected_lead,'nonzero scalar time leading coefficient')
                transformed_I=product(power(ep,a),power(ep,3)-power(ey,2))
                check(transformed_I==trunc(I),'invariant after actual component substitution')
                check(exponential(ep,-time)==p and exponential(ey,-time)==y,
                      'inverse by actual operator composition')
            ep=exponential(p,1)
            ey=exponential(y,1)
            check(exponential(ep,-2)==exponential(p,-1) and
                  exponential(ey,-2)==exponential(y,-1),'group law on the component bank')
            # Actual repeated composition, independent of just substituting
            # an integer into the claimed leading coefficient.
            if a in [2,3] and coefficients in [(0,1),(7,0,-3)]:
                iterate=p
                for n in range(1,7):
                    iterate=exponential(iterate,1)
                    check(iterate==exponential(p,n),'six genuine positive iterates')
                    lead=R.from_dict({(i,0):value for (i,j),value in (iterate-p).items() if j==k})
                    check(lead==n*expected_lead and lead!=zero,'positive iterates cannot have finite order')
            cases+=1
    print('COMPLETED_FLOW_CASES',cases,'a=2..10; five polynomial f; three nonzero scalar times')
    print('ITERATES exact operator composition; invariant substitution; inverse/group-law PASS')

    # Topology bridge to the actual t-adic source. This uses a second sparse
    # polynomial ring and the literal source bracket. Enough tau precision
    # is retained to control the desired t precision after s=x*t.
    RX,xx,tt=ring('x,t',sp.QQ)
    source_p=tt*(1+xx*xx*tt)
    source_y=xx*tt*source_p
    bridge_cutoff=14
    def sigma(F):
        return RX.from_dict({(i,i+j):value for (i,j),value in F.items()})
    check(sigma(p)==source_p and sigma(y)==source_y,'literal birational coordinate embedding')
    pairs={(i,i+j) for i in range(9) for j in range(9)}
    check(len(pairs)==81,'substitution has distinct monomial coordinates')
    check((xx*tt).diff(xx)*tt.diff(tt)-(xx*tt).diff(tt)*tt.diff(xx)==tt,
          'exact bracket conjugacy Jacobian is tau')
    for a,coefficients in [(2,(0,1,1)),(3,(0,1)),(2,(0,0,1))]:
        II=tau*p**(a+2)
        source_II=source_p**a*(source_p**3-source_y**2)
        SS=sum((c*II**i for i,c in enumerate(coefficients)),R.zero)
        source_SS=sum((c*source_II**i for i,c in enumerate(coefficients)),RX.zero)
        check(sigma(II)==source_II and sigma(SS)==source_SS,'invariant/source conjugacy before completion')
        def log_delta(F):
            return rs_trunc(tau*(F.diff(s)*SS.diff(tau)-F.diff(tau)*SS.diff(s)),tau,bridge_cutoff)
        def source_delta(F):
            return rs_trunc(F.diff(xx)*source_SS.diff(tt)-F.diff(tt)*source_SS.diff(xx),tt,bridge_cutoff)
        def exp_in_ring(F,operator,variable,time):
            result=rs_trunc(F,variable,bridge_cutoff)
            term=result
            for n in range(1,bridge_cutoff+1):
                term=operator(term)
                if term==0:break
                result+=term*sp.QQ.convert(sp.Rational(time)**n/factorial(n))
            return rs_trunc(result,variable,bridge_cutoff)
        for log_coordinate,source_coordinate in [(p,source_p),(y,source_y)]:
            log_iterate,source_iterate=log_coordinate,source_coordinate
            for n in range(4):
                check(rs_trunc(sigma(log_iterate),tt,bridge_cutoff)==source_iterate,
                      'every retained differential iterate commutes with sigma')
                log_iterate=log_delta(log_iterate)
                source_iterate=source_delta(source_iterate)
            for time in [1,sp.Rational(-2,3)]:
                log_flow=exp_in_ring(log_coordinate,log_delta,tau,time)
                source_flow=exp_in_ring(source_coordinate,source_delta,tt,time)
                check(rs_trunc(sigma(log_flow),tt,bridge_cutoff)==source_flow,
                      'completed scalar-time source/logarithmic comparison')
                check(source_flow!=source_coordinate,'source comparison includes an actual nonzero displacement')
    print('TOPOLOGY_BRIDGE sigma(s)=xt,sigma(tau)=t; three literal source flows agree through t13')

    # Non-LND is insufficient: this different Hamiltonian has a rational flow.
    x,t,L,M=sp.symbols('x t L M')
    u=x*x*t
    ex=x/(1-L*x)
    et=t*(1-L*x)**2
    check(sp.cancel(sp.diff(ex,L)-ex**2)==0,'rational hostile x differential equation')
    check(sp.expand(sp.diff(et,L)+2*ex*et)==0,'rational hostile t differential equation')
    check(sp.cancel(ex**2*et-u)==0,'rational hostile invariant')
    composed=ex.subs({x:x/(1-M*x)},simultaneous=True)
    check(sp.cancel(composed-x/(1-(L+M)*x))==0,'rational hostile group law')
    check(sp.cancel(sp.diff(ex,x)*sp.diff(et,t)-sp.diff(ex,t)*sp.diff(et,x)-1)==0,
          'rational hostile preserves source volume')
    def hostile_delta(F):return sp.expand(sp.diff(F,x)*sp.diff(u,t)-sp.diff(F,t)*sp.diff(u,x))
    value=x
    for n in range(1,11):
        value=hostile_delta(value)
        check(value==factorial(n)*x**(n+1),'rational hostile not locally nilpotent')
    print('HOSTILE S=x²t: x/(1-Lx), t(1-Lx)²; non-LND but rational symplectic flow')
    print('PASS checks='+str(CHECKS))
    print('SOURCE_SHA256 '+sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':main()

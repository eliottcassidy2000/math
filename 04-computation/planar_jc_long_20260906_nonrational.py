#!/usr/bin/env python3
"""Exact controls for the high-genus non-rational scalar-time theorem.

No repository mathematical imports. The universal conclusion is analytic;
this finite universe checks the identities, models, and hostile boundary.
"""
from hashlib import sha256
from pathlib import Path
from math import factorial
import sympy as S

CHECKS = 0
def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def main():
    p,y,c,s,tau,z,x,t,lam = S.symbols('p y c s tau z x t lam')
    D=p**3-y**2
    pp=s*s+tau
    def trunc(h, n):
        return S.Add(*[
            coef*tau**j[0] for j,coef in S.Poly(S.expand(h),tau).terms() if j[0]<n])
    def logbracket(h,q):
        return S.expand(tau*(S.diff(h,s)*S.diff(q,tau)-S.diff(h,tau)*S.diff(q,s)))
    def sourcebracket(h,q):
        return S.expand(S.diff(h,x)*S.diff(q,t)-S.diff(h,t)*S.diff(q,x))
    def flow(h,q,time,n):
        current=trunc(h,n); answer=current
        for j in range(1,n):
            current=trunc(logbracket(current,q),n)
            if current==0:break
            answer+=time**j*current/factorial(j)
        return trunc(answer,n)
    models=[]
    for a in range(2,21):
        I=p**a*D
        r=a//2
        power=r if a%2==0 else r+1
        model=p**(a+3)-c if a%2==0 else p**(a+4)-c*p
        identity=S.expand(p**(2*power)*(p**3-c/p**a)-model)
        check(identity==0,'birational invariant model '+str(a))
        check(S.Poly(model,p,domain=S.QQ.frac_field(c)).gcd(S.Poly(S.diff(model,p),p,domain=S.QQ.frac_field(c))).degree()==0,
              'geometrically squarefree model '+str(a))
        degree=S.degree(model,p); genus=(degree-1)//2
        check(degree%2==1 and genus==(a+1)//2+1 and genus>=2,'branch count genus '+str(a))
        dp=S.cancel(-D*S.diff(I,y)/p)
        dy=S.cancel(D*S.diff(I,p)/p)
        check(S.expand(dp-2*p**(a-1)*y*D)==0,'literal p bracket '+str(a))
        check(S.expand(S.diff(I,p)*dp+S.diff(I,y)*dy)==0,'exact invariant '+str(a))
        models.append((a,int(degree),int(genus)))
    forms=[z, z*z, 7+3*z*z-z**3, 2-2*z**3]
    for a in range(2,8):
        I=tau*pp**(a+2)
        for f in forms:
            k=min(j[0] for j,coef in S.Poly(f-f.subs(z,0),z).terms() if coef)
            fk=S.expand(f).coeff(z,k)
            q=trunc(f.subs(z,I),k+2)
            dp=logbracket(pp,q)
            leading=2*k*fk*s**((2*a+4)*k+1)
            check(trunc(dp,k)==0 and S.expand(dp).coeff(tau,k)==leading,'universal leading p formula')
            ds=logbracket(s,q);dt=logbracket(tau,q)
            check(trunc(ds,k)==0 and trunc(dt,k+1)==0,'valuation-raising vector field')
            ep=flow(pp,q,lam,k+2)
            check(S.expand(ep-pp).coeff(tau,k)==lam*leading,'nonzero scalar-time coefficient')
            check(trunc(logbracket(I,q),k+2)==0,'invariant under truncated carrier')
    # Independent literal source versus logarithmic bracket controls.
    px=t*(1+x*x*t); yx=x*t*px; Dx=S.expand(px**3-yx**2)
    check(S.expand(Dx-t**3*(1+x*x*t)**2)==0,'source cusp identity')
    for a in range(2,5):
        Ix=px**a*Dx
        check(S.expand(sourcebracket(px,Ix)-2*px**(a-1)*yx*Dx)==0,'literal source bracket')
    # The source comparison is injective on all formal series by its
    # monomial rule; this finite bank checks addresses and intertwining.
    addresses=[(i,i+j) for i in range(9) for j in range(9-i)]
    check(len(addresses)==len(set(addresses))==45,'injective source comparison monomial addresses')
    for a in (2,3):
        Ilog=tau*pp**(a+2)
        for hlog,hsource in ((pp,px),(s*pp,yx)):
            transformed=logbracket(hlog,Ilog).subs({s:x*t,tau:t},simultaneous=True)
            check(S.expand(transformed-sourcebracket(hsource,px**a*Dx))==0,'source and logarithmic derivations intertwine')
    # Finite group law and invariant controls by applying the derivation
    # directly to series (no substitution of growing rational expressions).
    for a in (2,3):
        I=tau*pp**(a+2);q=trunc(I,5)
        for h in (pp,s*pp,I):
            left=flow(flow(h,q,1,5),q,2,5)
            right=flow(h,q,3,5)
            check(S.expand(left-right)==0,'finite scalar additive group law')
            check(S.expand(flow(flow(h,q,1,5),q,-1,5)-trunc(h,5))==0,'finite inverse control')
        check(S.expand(flow(I,q,1,5)-trunc(I,5))==0,'finite invariant preservation')
    # Positive rational hostile to the invalid non-LND inference.
    u=x*x*t; xp=x/(1-lam*x);tp=t*(1-lam*x)**2
    check(S.cancel(xp*xp*tp-u)==0,'rational hostile preserves invariant')
    check(S.cancel(S.diff(xp,x)*S.diff(tp,t)-S.diff(xp,t)*S.diff(tp,x)-1)==0,'rational hostile preserves volume')
    check(S.cancel(xp/(1+lam*xp)-x)==0,'rational hostile has inverse')
    check(S.cancel(tp*(1+lam*xp)**2-t)==0,'rational hostile inverse t')
    current=x
    for n in range(1,9):
        current=sourcebracket(current,u)
        check(S.expand(current-factorial(n)*x**(n+1))==0,'non-LND hostile iterate')
    print('UNIVERSE model_a=2..20; carrier_a=2..7; f=z,z^2,7+3z^2-z^3,2-2z^3')
    print('MODELS (a,odd_degree,genus)',models)
    print('FLOW_CONTROLS a=2,3 through tau^4; scalar times1,2,3,-1; invariant and group law')
    print('HOSTILE Hamiltonian_x^2t non_LND_rational_symplectic_flow; eight_nonzero_iterates')
    print('CHECKS',CHECKS,'PASS')
    print('SCOPE finite identity controls; all-a all-f nonrationality uses analytic genus proof; JC2 OPEN')
    print('SOURCE_SHA256',sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':main()

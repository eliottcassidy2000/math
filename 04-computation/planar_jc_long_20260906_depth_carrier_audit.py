#!/usr/bin/env python3
"""Independent generic quotient-rule audit of the maximal depth carrier."""
from pathlib import Path
from hashlib import sha256
import sympy as S

def main():
    p,y,r,rp,ry=S.symbols('p y r rp ry');D=p**3-y*y
    x=y*p/D;u=y*y/D
    dp=-p*D*ry;dy=D*(2*r+p*rp)
    derived_x=S.diff(x,p)*dp+S.diff(x,y)*dy
    derived_u=S.diff(u,p)*dp+S.diff(u,y)*dy
    expected_x=p*(1+2*u)*(2*r+p*rp)+p*y*(2+3*u)*ry
    expected_u=(1+u)*y*(4*r+2*p*rp+3*y*ry)
    checks=0
    def check(ok,label):
        nonlocal checks
        checks+=1
        if not ok:raise RuntimeError(label)
    for variable in (r,rp,ry):
        check(S.cancel(S.diff(derived_x-expected_x,variable))==0,'generic quotient x coefficient')
        check(S.cancel(S.diff(derived_u-expected_u,variable))==0,'generic quotient u coefficient')
    check(S.cancel(derived_u-p**3*y*(4*r+2*p*rp+3*y*ry)/D)==0,'affine source residue exact')
    check(S.gcd(p,D)==1 and S.gcd(p**3*y,D)==1,'necessary denominator coprimality')
    check(S.expand(2*p*S.diff(D,p)+3*y*S.diff(D,y)-6*D)==0,'Euler descends to cusp')
    for a in range(9):
        for b in range(9-a):
            R=p**a*y**b;z=S.Symbol('z')
            ER=4*R+2*p*S.diff(R,p)+3*y*S.diff(R,y)
            check(S.expand(ER.subs({p:z*z,y:z**3})-(4+2*a+3*b)*z**(2*a+3*b))==0,
                  'cusp weight diagonal action')
    ss,tau=S.symbols('s tau');pl=ss*ss+tau;yl=ss*pl
    hostile_u=expected_u.subs({r:1,rp:0,ry:0}).subs({p:pl,y:yl},simultaneous=True)
    check(S.limit(tau*hostile_u,tau,0)==4*ss**5,'depth preserving does not lower depth')
    check(S.rem(D,p,p)==-y*y,'fixed-source D outside maximal depth carrier')
    print('UNIVERSE generic_R_Rp_Ry coefficient identities;45 cusp monomials; two carrier hostiles')
    print('METHOD quotient differentiation of x=yp/D,u=y²/D; no source-derivative producer import')
    print('CHECKS',checks,'PASS')
    print('SOURCE_SHA256',sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=='__main__':main()

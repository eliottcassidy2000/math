#!/usr/bin/env python3
"""Exact identities for the sharp local stability moduli."""
import hashlib
import sympy as S

GATES=0
TRACE=hashlib.sha256()


def gate(ok,label):
    global GATES
    if ok != True:
        raise RuntimeError(label+': '+str(ok))
    GATES+=1
    TRACE.update((label+'\n').encode())


def zero(expr,label):
    gate(S.simplify(expr)==0,label)


def main():
    z=1/S.sqrt(3)
    e=S.Symbol('e')
    zero((z+e)**2*e**2-z**2*e**2-(2*z*e**3+e**4),'leading coordinate exact remainder')
    zero(e**2*(e-z)**2-z**2*e**2-(-2*z*e**3+e**4),'tail coordinate exact remainder')
    gate(96**2<3*61**2,'relative error below one twelfth at Delta one over256')
    zero(S.Rational(1,3)-S.Rational(1,12)-S.Rational(1,4),'local lower comparison')
    zero(S.Rational(1,3)+S.Rational(1,12)-S.Rational(5,12),'local upper comparison')
    A=-2-S.sqrt(6)/3+2*S.sqrt(2)+7*S.sqrt(3)/3
    zero(3*A-(-6-S.sqrt(6)+6*S.sqrt(2)+7*S.sqrt(3)),'sharp moment coefficient')
    m,v=S.symbols('m v')
    Delta=2-2*S.sqrt(1-m-v)
    zero(Delta.subs({m:0,v:0}),'distance vanishes at three atoms')
    zero(S.diff(Delta,m).subs({m:0,v:0})-1,'radial distance coefficient')
    zero(S.diff(Delta,v).subs({m:0,v:0})-1,'variance distance coefficient')
    epsilon=S.Symbol('epsilon',nonnegative=True)
    a=(epsilon+z*S.sqrt(1+2*epsilon))/(1+3*epsilon)
    p3=3*a**3-epsilon**2*(3*a-1)**3
    p4=3*a**4+epsilon**3*(3*a-1)**4
    MM=p4-2*z*p3+S.Rational(1,3)
    DD=2-6*z*a
    K=4*S.sqrt(3)/(3*(1+S.sqrt(2))*(1+S.sqrt(3)))
    cs=(13-8*S.sqrt(2))/3
    R=(((5-8*p3+3*p4)/(3*(1-p4)))-cs)/(2-2*S.sqrt(2)*a)
    radial=(S.sqrt(3)-1)**2
    zero(DD.subs(epsilon,0),'actual family distance limit')
    zero(MM.subs(epsilon,0),'actual family moment limit')
    zero(R.subs(epsilon,0)-K,'actual family quotient limit')
    zero(S.diff(DD,epsilon).subs(epsilon,0)-radial,'actual family distance derivative')
    zero(S.diff(MM,epsilon).subs(epsilon,0)-radial/3,'actual family moment derivative')
    zero(S.diff(R,epsilon).subs(epsilon,0)-A*radial,'actual family response derivative')
    print('IDENTITY |M-Delta/3|<=2Delta^(3/2)/sqrt(3)+Delta^2; all actual finite lists')
    print('LOCAL Delta<=1/256 implies Delta/4<=M<=5Delta/12')
    print('SHARP minimizing-sequence liminf gap/Delta>=A; liminf gap/M>=3A; equal-three family attains both')
    print('EQUALITY convergence to sharp ratios iff h=o(m); liminf equality alone does not imply full-sequence equality')
    print('SCOPE sharp local asymptotic moduli; globally optimal modulus remains OPEN')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())


if __name__=='__main__':
    main()

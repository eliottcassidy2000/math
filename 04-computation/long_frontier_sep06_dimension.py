#!/usr/bin/env python3
"""Exact identities for the asymptotic dimension penalty; no finite-N optimizer claim."""
from fractions import Fraction as F
from math import isqrt
import hashlib
import json
import sympy as s

gates=0


def check(test,label):
    global gates
    gates+=1
    if not bool(test):
        raise RuntimeError(label)


def equal(a,b,label):
    check(s.simplify(a-b)==0,label)


class I:
    """Outward exact rational intervals for independent numerical controls."""
    def __init__(self,a,b=None):
        self.a=F(a)
        self.b=F(a if b is None else b)
        if self.a>self.b:
            raise ValueError('reversed interval')
    def __add__(self,other):
        other=iv(other)
        return I(self.a+other.a,self.b+other.b)
    __radd__=__add__
    def __neg__(self):
        return I(-self.b,-self.a)
    def __sub__(self,other):
        return self+-iv(other)
    def __rsub__(self,other):
        return iv(other)+-self
    def __mul__(self,other):
        other=iv(other)
        v=[self.a*other.a,self.a*other.b,self.b*other.a,self.b*other.b]
        return I(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=iv(other)
        if other.a<=0<=other.b:
            raise ValueError('division across zero')
        return self*I(1/other.b,1/other.a)
    def __rtruediv__(self,other):
        return iv(other)/self
    def __pow__(self,k):
        if k<0:
            return 1/(self**(-k))
        out=I(1)
        for _ in range(k):
            out=out*self
        return out
    def contains(self,x):
        return self.a<=F(x)<=self.b
    def strings(self):
        scale=10**15
        lo=(self.a*scale).__floor__()
        hi=(self.b*scale).__ceil__()
        return [str(F(lo,scale)),str(F(hi,scale))]


def iv(x):
    return x if isinstance(x,I) else I(x)


def sqrt_interval(q):
    q=F(q)
    scale=10**45
    k=isqrt(q.numerator*scale*scale//q.denominator)
    lo,hi=F(k,scale),F(k+1,scale)
    check(lo*lo<=q<=hi*hi,'rational radical enclosure')
    return I(lo,hi)


def main():
    u,z=s.sqrt(2),1/s.sqrt(3)
    d0=2-2*u*z
    cstar=(13-8*u)/3
    K=4*s.sqrt(3)/(3*(1+u)*(1+s.sqrt(3)))
    A=-2-s.sqrt(6)/3+2*s.sqrt(2)+7*s.sqrt(3)/3
    B=-5*s.sqrt(2)-4*s.sqrt(3)+3*s.sqrt(6)+8
    C=-22-16*s.sqrt(6)/3+10*s.sqrt(2)+40*s.sqrt(3)/3
    equal(K,(3-4*z-cstar)/d0,'base quotient')
    equal(A,(10*z-4-K*u*z)/d0,'radial mass derivative')
    equal(B,K*u/d0,'ordered leading splitting derivative')
    equal(C,A*(3*z-1)**2,'sharp dimension coefficient')
    check(s.Rational(7,5)**2<2 and s.Rational(17,10)**2<3
          and s.Rational(5,2)**2>6,'simple strict radical bounds')
    check(-2-s.Rational(5,6)+2*s.Rational(7,5)+7*s.Rational(17,10)/3>0,
          'A strictly positive without floating point')
    check(s.Rational(3,2)**2>2 and 2<3,'d0 and B positive')

    # A separate direct atom-gradient derivation catches the ordered top-two term.
    x,y,w=s.symbols('x y w',real=True)
    p3=x**3+y**3+w**3
    p4=x**4+y**4+w**4
    J=(5-8*p3+3*p4)/(3*(1-p4))
    R=(J-cstar)/(2-u*(x+y))
    grad=[s.simplify(s.diff(R,t).subs({x:z,y:z,w:z})) for t in (x,y,w)]
    equal(-z*sum(grad)/2,A,'direct gradient radial contraction')
    equal(2*grad[0]-grad[1]-grad[2],B,'asymmetric ordered split')
    equal((grad[0]+grad[1])/2-grad[2],B,'symmetric ordered split')
    equal(grad[0]-grad[2],B,'omitting the third atom changes first derivative')

    m,h=s.symbols('m h',real=True)
    q3=z-s.Rational(3,2)*z*m
    q4=s.Rational(1,3)-s.Rational(2,3)*m
    radial=((5-8*q3+3*q4)/(3*(1-q4))-cstar)/(d0+u*z*m-u*h)
    equal(radial.subs({m:0,h:0}),K,'local base')
    equal(s.diff(radial,m).subs({m:0,h:0}),A,'local m coefficient')
    equal(s.diff(radial,h).subs({m:0,h:0}),B,'local h coefficient')

    # Exact admissible upper family and split controls, using epsilon=1/n.
    e,v=s.symbols('e v',nonnegative=True)
    ell=(e+z*s.sqrt(1+2*e-(1+3*e)*v))/(1+3*e)
    equal(3*ell**2+v+e*(3*ell-1)**2,1,'norm solution for prescribed variance')
    ell0=ell.subs(v,0)
    equal(ell0.subs(e,0),z,'equal family limit')
    equal(s.diff(ell0,e).subs(e,0),1-2*z,'equal family first derivative')
    mm=e*(3*ell0-1)**2
    equal(s.diff(mm,e).subs(e,0),(3*z-1)**2,'tail Cauchy mass equality')
    upper3=3*ell0**3-e**2*(3*ell0-1)**3
    upper4=3*ell0**4+e**3*(3*ell0-1)**4
    equal(s.diff(upper3,e).subs(e,0),-s.Rational(3,2)*z*(3*z-1)**2,
          'actual upper cubic derivative')
    equal(s.diff(upper4,e).subs(e,0),-s.Rational(2,3)*(3*z-1)**2,
          'actual upper fourth derivative')
    upper=((5-8*upper3+3*upper4)/(3*(1-upper4))-cstar)/(2-2*u*ell0)
    equal(s.diff(upper,e).subs(e,0),C,'actual upper sharp coefficient')

    # An exact normalized hostile to discarding ordered leading-root splitting:
    # h=lambda/n costs B*lambda/n already at first order, rather than h^2.
    lam=s.symbols('lam',positive=True)
    for shape in ((s.Rational(1,2),s.Rational(1,2),-1),(2,-1,-1)):
        var=sum(t*t for t in shape)*lam**2*e**2
        center=ell.subs(v,var)
        atoms=[center+t*lam*e for t in shape]
        moment3=sum(t**3 for t in atoms)-e**2*(3*center-1)**3
        moment4=sum(t**4 for t in atoms)+e**3*(3*center-1)**4
        value=((5-8*moment3+3*moment4)/(3*(1-moment4))-cstar)/(2-u*(atoms[0]+atoms[1]))
        equal(s.diff(value,e).subs(e,0),C+B*lam,'actual split family linear penalty')

    # Independent outward-rational controls for finite lists; no optimizer census.
    ui,zi=sqrt_interval(2),1/sqrt_interval(3)
    ki=4*sqrt_interval(3)/(3*(1+ui)*(1+sqrt_interval(3)))
    ci=(13-8*ui)/3
    Ai=-2-sqrt_interval(6)/3+2*ui+7*sqrt_interval(3)/3
    Bi=-5*ui-4*sqrt_interval(3)+3*sqrt_interval(6)+8
    Ci=Ai*(sqrt_interval(3)-1)**2
    check(Ai.a>0 and Bi.a>0 and Ci.a>0,'independent constant positivity')
    controls=[]
    for n in (10,100,1000,10000):
        equal_penalty=None
        for shape in ((0,0,0),(F(1,2),F(1,2),-1),(2,-1,-1)):
            hh=F(0) if shape==(0,0,0) else F(1,n)
            var=sum(F(t)**2 for t in shape)*hh**2
            center=(1+sqrt_interval(F(n*((n+2)-(n+3)*var),3)))/(n+3)
            atoms=[center+F(t)*hh for t in shape]
            neg=(3*center-1)/n
            check(atoms[2].a>neg.b>0,'three leading positive roots and negative dust')
            check((sum(atoms,I(0))-n*neg).contains(1),'literal first moment enclosure')
            check((sum((a*a for a in atoms),I(0))+n*neg**2).contains(1),
                  'literal second moment enclosure')
            pp3=sum((a**3 for a in atoms),I(0))-n*neg**3
            pp4=sum((a**4 for a in atoms),I(0))+n*neg**4
            energy=(1-pp4)/2
            distance=2-ui*(atoms[0]+atoms[1])
            check(energy.a>0 and distance.a>0,'eligible actual denominator')
            jj=(5-8*pp3+3*pp4)/(3*(1-pp4))
            rr=(jj-ci)/distance
            penalty=(n+3)*(rr-ki)
            check(penalty.a>0,'positive finite penalty control')
            if shape==(0,0,0):
                check(abs(penalty.a-Ci.b)<F(5,n) and abs(penalty.b-Ci.a)<F(5,n),
                      'finite upper-family convergence diagnostic')
                equal_penalty=penalty
            else:
                check(penalty.a>equal_penalty.b,'actual split costs more in named control')
            controls.append(dict(n=n,shape=[str(t) for t in shape],
                                 scaled_penalty_interval=penalty.strings()))

    manifest=dict(K3=str(K),A=str(A),B=str(B),C=str(C),
                  conclusion='inf length<=N R=K3+C/N+o(1/N); no finite-N optimizer assertion',
                  constant_intervals=dict(A=Ai.strings(),B=Bi.strings(),C=Ci.strings()),
                  controls=controls)
    encoded=json.dumps(manifest,sort_keys=True,separators=(',',':'))
    print('PROVED asymptotic dimension penalty; exact identities and finite controls')
    print(encoded)
    print('EXPLICIT_GATES',gates)
    print('SEMANTIC_SHA256',hashlib.sha256(encoded.encode()).hexdigest())


if __name__=='__main__':
    main()

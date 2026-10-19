#!/usr/bin/env python3
"""Finite-exact controls for uniform neighbor-count coefficient bounds.

The arbitrary-r theorem is proved in the companion note. This checker uses
actual dense univariate integer polynomials for the division recurrences,
then checks odd multiplication independently against ordinary group addition.
It also checks first jets without factoring or canceling polynomial pairs.
"""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd
import json
import sys
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ

sys.set_int_max_str_digits(0)
COUNTS = Counter()


def check(ok, label, kind='exact_integer'):
    COUNTS[kind] += 1
    if not ok:
        raise RuntimeError(label)


def normalized_division(x, aa, bb):
    g = x**3 + aa*x + bb
    base = {1:x**0, 2:x**0,
            3:3*x**4 + 6*aa*x*x + 12*bb*x - aa*aa,
            4:2*(x**6 + 5*aa*x**4 + 20*bb*x**3 - 5*aa*aa*x*x
                 - 4*aa*bb*x - aa**3 - 8*bb*bb)}
    @lru_cache(None)
    def fn(n):
        if n in base:
            return base[n]
        m=n//2
        if n % 2:
            first=fn(m+2)*fn(m)**3
            second=fn(m-1)*fn(m+1)**3
            return 16*g*g*first-second if m % 2==0 else first-16*g*g*second
        return fn(m)*(fn(m+2)*fn(m-1)**2-fn(m-2)*fn(m+1)**2)
    return fn


def norm(p):
    return sum(abs(int(c)) for c in p.values())


def peval(p, v):
    out=Q(0)
    for j in range(p.degree(), -1, -1):
        out=out*v+int(p.get((j,), 0))
    return out


class Jet:
    """Exact first jet, with division only by units at h=0."""
    def __init__(self, a=0, b=0):
        self.a,self.b=Q(a),Q(b)
    @staticmethod
    def coerce(v):
        return v if isinstance(v,Jet) else Jet(v)
    def __add__(self,v):
        v=self.coerce(v);return Jet(self.a+v.a,self.b+v.b)
    __radd__=__add__
    def __neg__(self):
        return Jet(-self.a,-self.b)
    def __sub__(self,v):
        return self+-self.coerce(v)
    def __rsub__(self,v):
        return self.coerce(v)+-self
    def __mul__(self,v):
        v=self.coerce(v);return Jet(self.a*v.a,self.a*v.b+self.b*v.a)
    __rmul__=__mul__
    def __truediv__(self,v):
        v=self.coerce(v)
        check(v.a!=0,'nonzero jet denominator','exact_symbolic_jet')
        return Jet(self.a/v.a,(self.b*v.a-self.a*v.b)/(v.a*v.a))
    def __rtruediv__(self,v):
        return self.coerce(v)/self
    def __pow__(self,n):
        if n<0:return (1/self)**(-n)
        out=Jet(1)
        for _ in range(n):out=out*self
        return out


def pjet(p):
    return Jet(int(p.get((0,),0)),int(p.get((1,),0)))


def add(P,R,a2,a4):
    if P is None:return R
    if R is None:return P
    x,y=P;xx,yy=R
    if x==xx:
        if y==-yy:return None
        slope=(3*x*x+2*a2*x+a4)/(2*y)
    else:slope=(yy-y)/(xx-x)
    outx=slope*slope-a2-x-xx
    return outx,-y+slope*(x-outx)


def odd_multiple(x,y,fn,d):
    return (x-4*y*y*fn(d-1)*fn(d+1)/fn(d)**2,
            y*fn(2*d)/fn(d)**4)


def inverse(X,Y,k):
    K=k*k
    u=(k*Y+X-K)/(X*(X-K))
    z=(2*u*u*X-K-2*u-K*u*u)/k
    return u,z


def spread(t):
    return 4*t*(t*t-1)/(t*t+1)**2


def check_port(eps):
    R,h=ring('h',ZZ)
    t=1+h;D=t*t+1
    A=t*t+2*t-1 if eps==1 else -t*t+2*t+1
    aa=-81*D**4-27*A**4
    bb=-54*A**6+486*A*A*D**4
    x=9*D*(t*t-4*t+1)-3*A*A
    y=-216*(t*t if eps==1 else t)*(t-1)*D
    check(y*y==x**3+aa*x+bb,'polynomial short model','exact_symbolic')
    fn=normalized_division(x,aa,bb)
    H=max(norm(x),norm(aa),norm(bb),1)
    polyrows=[]
    for n in range(1,15):
        g=fn(n);gn=norm(g)
        check(g.degree()<=2*n*n,'degree bound')
        check(gn<=H**(n*n)*(n+1)**(4*n*n),'substituted norm bound')
        if n%2:
            check(g[(0,)]!=0,'nonzero odd cusp denominator')
        polyrows.append({'n':n,'degree':g.degree(),'norm_bits':gn.bit_length()})
    check(fn(3)[(0,)] == -5184**2,'cusp F3 value')
    x0,y0=pjet(x),pjet(y)
    D0,A0=pjet(D),pjet(A)
    k0=A0/D0
    jetrows=[]
    for d in (3,5,7):
        xx,yy=odd_multiple(x0,y0,lambda n:pjet(fn(n)),d)
        check(xx.a==-48 and yy.a==0,'odd cusp specialization','exact_symbolic_jet')
        if eps==1:
            px=-3*A0*A0;py=27*A0*D0*D0
            slope=(-py-yy)/(px-xx)
            newx=slope*slope-xx-px
            xx,yy=newx,-yy+slope*(xx-newx)
        X=(xx+3*A0*A0)/(9*D0*D0)
        Y=yy/(27*D0**3)
        u,z=inverse(X,Y,k0);v=eps*u
        check((v.a,v.b,z.a)==(1,d,2*eps),'first jet v and z','exact_symbolic_jet')
        jetrows.append({'d':d,'v_constant':str(v.a),'v_linear':str(v.b),'z_constant':str(z.a)})
    rationalrows=[]
    for hv in (Q(1,64),Q(1,1024),Q(1,8)):
        tv=1+hv;dv=tv*tv+1;av=tv*tv+2*tv-1 if eps==1 else -tv*tv+2*tv+1
        kv=av/dv;K=kv*kv
        B=((tv*tv-4*tv+1)/dv,-8*(tv*tv if eps==1 else tv)*(tv-1)/(dv*dv))
        xx,yy=peval(x,hv),peval(y,hv)
        evaluated=lambda n:peval(fn(n),hv)
        multiple=None
        for d in range(1,8):
            multiple=add(multiple,B,-K,-1)
            if d not in (3,5,7):continue
            sx,sy=odd_multiple(xx,yy,evaluated,d)
            converted=((sx+3*av*av)/(9*dv*dv),sy/(27*dv**3))
            check(converted==multiple,'division multiple vs ordinary additions','exact_rational')
            point=add(multiple,(Q(0),-kv),-K,-1) if eps==1 else multiple
            if point is None or point[0] in (0,K):
                rationalrows.append({'h':str(hv),'d':d,'chart':'excluded'})
                continue
            u,z=inverse(*point,kv);v=eps*u
            root=abs(z)/(v*v+1)
            check(root*root==1+eps*spread(tv)-eps*spread(v),'exact contact square','exact_rational')
            admissible=v>1 and 0<spread(v)<Q(1,2) and spread(v)!=spread(tv)
            if hv==Q(1,1024):
                check(tv<v<Q(5,4),'small cusp ordering','exact_rational')
                check(admissible,'small cusp success','exact_rational')
                n,m=v.numerator,v.denominator;div=2 if n%2 and m%2 else 1
                a,b,c=(n*n-m*m)//div,2*n*m//div,(n*n+m*m)//div
                check(gcd(a,b)==1 and a*a+b*b==c*c and c*c>4*a*b,
                      'primitive admissible output','exact_integer')
            rationalrows.append({'h':str(hv),'d':d,'admissible':admissible,
                                 'parameter_bits':max(abs(v.numerator).bit_length(),v.denominator.bit_length())})
    return {'port':eps,'polynomials':polyrows,'jets':jetrows,'rational_probes':rationalrows}


def cutoff_checks():
    # Explicit unit-constant polynomials exercise the cutoff itself.  These
    # controls do not stand in for constructing the full arbitrary-r atlas.
    examples=[(1,-1000,7),(7,-1000000,-500,12),(3,0,0,-999999999),(1,)]
    out=[]
    for coeff in examples:
        M=sum(abs(a) for a in coeff)
        b=max(4,1+M.bit_length());eta=Q(1,2**b)
        for hv in (Q(0),eta/2,eta):
            value=sum(Q(a)*hv**i for i,a in enumerate(coeff))
            check(value>=Q(coeff[0],2),'dyadic cutoff positivity','exact_rational')
        delta=spread(1+eta)
        check(delta>=eta and 2**(b+2)*delta>Q(5,2),'cutoff period budget','exact_rational')
        for r in (1,2,3,8):
            W=Q(5,7)*(6+b);N=(20*r*W.numerator)//W.denominator;Mindex=2*r*N
            check(Mindex>0 and Mindex%2==0,'even recurrence budget')
        out.append({'coefficient_bits':M.bit_length(),'b':b,'J':b+2})
    return out


if __name__=='__main__':
    ports=[check_port(eps) for eps in (1,-1)]
    out={'status':'FINITE-EXACT','ports':ports,'cutoff_controls':cutoff_checks(),'checks':dict(COUNTS)}
    print(json.dumps(out,indent=2,sort_keys=True))

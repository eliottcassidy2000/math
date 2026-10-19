#!/usr/bin/env python3
"""Finite-exact controls for the fixed-count contact atlas.

This deliberately does not duplicate the predecessor's complete census.
Checks remain active under python -O. The theorem for arbitrary fixed r
uses the formal-coordinate argument in the companion note.
"""
from collections import Counter
from fractions import Fraction as F
from math import gcd
import json
import sys
from sympy import QQ
from sympy.polys.fields import field

sys.set_int_max_str_digits(0)
CHECKS = Counter()

def check(ok, label, kind='exact_rational'):
    CHECKS[kind] += 1
    if not ok:
        raise RuntimeError(label)

def spread(t):
    return 4*t*(t*t-1)/(t*t+1)**2

def add(A,B,k):
    if A is None: return B
    if B is None: return A
    x,y=A;xx,yy=B
    if x==xx and y==-yy: return None
    m=(3*x*x-2*k*k*x-1)/(2*y) if A==B else (yy-y)/(xx-x)
    X=m*m+k*k-x-xx
    return X,-y+m*(x-X)

def neg(A):
    return None if A is None else (A[0],-A[1])

def inverse(A,k):
    if A is None: return None
    X,Y=A;K=k*k
    if X==0 or X==K: return None
    u=(k*Y+X-K)/(X*(X-K))
    z=(2*u*u*X-K-2*u-K*u*u)/k
    return u,z

def anchor(t,eps,sheet=1):
    k=(t*t+2*t-1)/(t*t+1) if eps==1 else (-t*t+2*t+1)/(t*t+1)
    K=k*k;u=eps*t;z=sheet*(t*t+1)
    X=(k*z+K+2*u+K*u*u)/(2*u*u)
    return k,(X,(X-K)*(X*u-1)/k)

def normalize(t):
    n,m=t.numerator,t.denominator
    divisor=2 if n%2 and m%2 else 1
    a,b,c=(n*n-m*m)//divisor,2*n*m//divisor,(n*n+m*m)//divisor
    if a%2==0: a,b=b,a
    check(a>0 and b>0 and a*a+b*b==c*c,'positive PPT')
    check(gcd(a,b)==gcd(a,c)==gcd(b,c)==1,'primitive PPT')
    check(F(2*a*b,c*c)==spread(t),'preserved spread')
    return (a,b,c)

def polynomial_eval(p,x):
    return p.evaluate(0,x)

def evaluate(g,x):
    return polynomial_eval(g.numer,x)/polynomial_eval(g.denom,x)

def symbolic_controls():
    FF,t=field('t',QQ)
    rows=[]
    for eps in (1,-1):
        k,Q=anchor(t,eps,1 if eps==1 else -1)
        P=(FF.zero,k)
        B=add(Q,P,k) if eps==1 else Q
        check(evaluate(B[0],QQ(1))==-1,'B tends to smooth torsion','exact_symbolic')
        check(evaluate(B[1].diff(t),QQ(1))==-2,'B Y slope -2','exact_symbolic')
        S=add(B,B,k)
        w=-S[0]/S[1]
        check(evaluate(w,QQ(1))==0,'S tends to identity','exact_symbolic')
        check(evaluate(w.diff(t),QQ(1))==-1,'formal w slope -1','exact_symbolic')
        R=Q
        for j in range(1,4):
            R=add(R,S,k)
            u,z=inverse(R,k);v=eps*u;d=2*j+1
            check(evaluate(v,QQ(1))==1,'v_d cusp value','exact_symbolic')
            check(evaluate(v.diff(t),QQ(1))==d,'v_d derivative d','exact_symbolic')
            check(evaluate(z,QQ(1))==2*eps,'z_d cusp value','exact_symbolic')
            check(z*z==k*k*(u*u+1)**2-4*u*(u*u-1),'section quartic identity','exact_symbolic')
            root=z/(v*v+1)
            check(root*root==1+eps*spread(t)-eps*spread(v),'section square identity','exact_symbolic')
            coeff=[evaluate(v,QQ(1))]
            deriv=v
            fac=1
            for order in range(1,4):
                deriv=deriv.diff(t);fac*=order
                coeff.append(evaluate(deriv,QQ(1))/fac)
            rows.append({'port':eps,'d':d,'parameter_degrees':[v.numer.degree(),v.denom.degree()],
                         'taylor_through_h3':[str(x) for x in coeff]})
    return rows

def rational_cusp_controls():
    rows=[]
    for t in [F(1025,1024),F(4097,4096),F(1000001,1000000)]:
        for eps in (1,-1):
            k,Q=anchor(t,eps,1 if eps==1 else -1)
            B=add(Q,(F(0),k),k) if eps==1 else Q
            S=add(B,B,k);R=Q;last=t;petals=[]
            for d in (3,5,7):
                R=add(R,S,k);u,z=inverse(R,k);v=eps*u
                check(last<v<F(5,4),'ordered cusp parameters')
                check(0<spread(t)<spread(v)<F(1,2),'admissible changed spread')
                root=abs(z)/(v*v+1)
                check(root*root==1+eps*spread(t)-eps*spread(v),'contact certificate')
                T=normalize(v);check(T not in petals,'distinct primitive cusp petals')
                petals.append(T);last=v
                rows.append({'t':str(t),'port':eps,'d':d,'hypotenuse_bits':T[2].bit_length()})
    return rows

def fallback_controls():
    # Interior probes deliberately do not purport to be a complete census.
    rows=[]
    for t in [F(9,8),F(6,5),F(5,4),F(13,10)]:
        f0=spread(t);check(0<f0<F(1,2),'admitted compact anchor')
        for eps in (1,-1):
            k,Q=anchor(t,eps,1);A=add((F(0),k),(F(0),k),k)
            J=0;D=2+f0 if eps==1 else F(2)
            while 2**J*f0<D: J+=1
            r=3;W=F(5,7)*(4+J);budget=2*r*(20*r*W.numerator//W.denominator)
            R=None;found=[];indices=[]
            for j in range(1,budget//2+1):
                R=add(R,A,k)
                for sign in (1,-1):
                    pair=inverse(add(Q,R if sign==1 else neg(R),k),k)
                    if pair is None: continue
                    u,z=pair;v=eps*u
                    if v<=1 or z<=0: continue
                    f1=spread(v)
                    if not(0<f1<F(1,2)) or f1==f0 or f1 in found: continue
                    root=z/(v*v+1)
                    check(root*root==1+eps*f0-eps*f1,'fallback square')
                    normalize(v);found.append(f1);indices.append(sign*2*j)
                    if len(found)==r:break
                if len(found)==r:break
            check(len(found)==r,'three distinct fallback contacts')
            check(all(abs(q)<=budget and q%2==0 for q in indices),'fixed-count budget')
            rows.append({'t':str(t),'port':eps,'r':r,'J':J,'budget':budget,'successful_q':indices})
    return rows

if __name__=='__main__':
    out={'status':'FINITE-EXACT', 'symbolic':symbolic_controls(),
         'cusp_probes':rational_cusp_controls(),'fallback_probes':fallback_controls()}
    out['checks']=dict(CHECKS)
    print(json.dumps(out,indent=2,sort_keys=True))

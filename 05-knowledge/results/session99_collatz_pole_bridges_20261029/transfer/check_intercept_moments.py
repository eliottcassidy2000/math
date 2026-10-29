#!/usr/bin/env python3
import json
from fractions import Fraction as F
import sympy as s
import mpmath as mp

def need(v,msg):
    if not v: raise RuntimeError(msg)
m,h,y,t=s.symbols('m h y t'); b=s.symbols('b')
uv=s.symbols('u1 v1 u2 v2 u3 v3'); u1,v1,u2,v2,u3,v3=uv
U=3*b*u1*u2*u3; V=b*u1*u2*u3+v1*u2*u3+v2*u3+v3
mom={(0,0):1,(1,0):1,(2,0):m,(0,1):y,(1,1):h,(0,2):t}
def integrate(p,q):
    ans=s.Rational(3**p,2**(p+q+1)-1)-s.Rational(3**p,2**(3*(p+q+1))-1)
    for exps,coeff in s.Poly(s.expand(U**p*V**q),b,*uv).terms():
        term=coeff*s.Rational(1,2**(3*(exps[0]+1))-1)
        for i in range(3): term*=mom[exps[1+2*i],exps[2+2*i]]
        ans+=term
    return s.expand(ans)
S=1+m+m*m
rootm=s.Rational(648,511)+s.Rational(9,511)*m**3
need(s.expand(integrate(1,1)-(rootm/3+h*S/21))==0,'UV expansion')
need(s.expand(integrate(0,2)-(rootm/9+t*S/7+s.Rational(2,63)*h*S+s.Rational(2,7)*y*h*(m+2)))==0,'V2 expansion')
thresholds=[]
for r in range(1,5):
    B=F(3**r,2**(3*(r+1))-1); L=F(3**r,2**(r+1)-1)-B
    value=L+8*B
    need((value<2)==(r<=3),'moment threshold')
    thresholds.append({'r':r,'L_plus_8B':str(value),'below_two':value<2})
mp.mp.dps=50
M=mp.findroot(lambda x:mp.mpf(648)/511+mp.mpf(9)/511*x**3-x,1.3)
SS=1+M+M*M; H=(M/3)/(1-SS/21)
T=(M/9+mp.mpf(2)/63*H*SS+mp.mpf(1)/6*H*(M+2))/(1-SS/7)
print(json.dumps({'status':'PASS','symbolic_mixed_equations':2,'thresholds':thresholds,'approximations':{'EU2':str(M),'EUV':str(H),'EV2':str(T)}},indent=2))

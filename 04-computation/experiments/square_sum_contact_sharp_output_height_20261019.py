#!/usr/bin/env python3
"""Small exact controls; no large orbit search or canonical-height numerics."""
from fractions import Fraction as F
from math import gcd
import json
import sys
import sympy as sp

checks = []
def check(condition, label):
    if not condition:
        raise RuntimeError(label)
    checks.append(label)

t, x, K, r, s, a, b = sp.symbols('t x K r s a b')
k = (t*t+2*t-1)/(t*t+1)
spread = lambda v: 4*v*(v*v-1)/(v*v+1)**2
check(sp.cancel(k*k-1-spread(t)) == 0, 'K=1+f(t)')
XQ=(t*t+4*t+1)/(t*t+1)
YQ=8*t*(t+1)/(t*t+1)**2
check(sp.cancel(YQ*YQ-(XQ-k*k)*(XQ*XQ-1)) == 0, 'Q is a section')
forward=(k*(t*t+1)+k*k+2*t+k*k*t*t)/(2*t*t)
check(sp.cancel(forward-XQ) == 0, 'Q is the prescribed anchor image')
check(sp.cancel(spread(-1/t)-spread(t)) == 0, 'alternate chart preserves spread')

cubic=(x-K)*(x*x-1)
Delta=16*sp.discriminant(cubic,x)
check(sp.factor(Delta-64*(K*K-1)**2) == 0, 'discriminant')
b2=-4*K; b4=-2
c4=b2*b2-24*b4
check(sp.expand(c4-16*(K*K+3)) == 0, 'c4')
j=sp.cancel(c4**3/Delta)
check(sp.cancel(j-64*(K*K+3)**3/(K*K-1)**2) == 0, 'j invariant')
jt=sp.cancel(j.subs(K,k*k)); jn,jd=sp.fraction(jt)
check(sp.degree(jn,t)==24 and sp.degree(jd,t)==22, 'j(t) has degree 24')
check(jn.subs(t,1)!=0 and jd.subs(t,1)==0, 'j(t) has a pole at t=1')
check(k.subs(t,F(5,4))==F(49,41), 'smooth nontorsion-test specialization')
check(spread(F(5,4))==F(720,1681), 'admissible nontorsion-test specialization')

fn=x**4+2*x*x-8*K*x+4*K*K+1
gn=4*(x-K)*(x*x-1)
slope_squared=(3*x*x-2*K*x-1)**2/(4*cubic)
check(sp.cancel(slope_squared+K-2*x-fn/gn)==0, 'duplication formula')
check(sp.cancel((fn/gn).subs(x,0)-K-1/(4*K))==0, 'x(2P)')
check(sp.expand(fn.subs(x,1)-4*(K-1)**2)==0, 'root x=1')
check(sp.expand(fn.subs(x,-1)-4*(K+1)**2)==0, 'root x=-1')
check(sp.expand(fn.subs(x,K)-(K*K-1)**2)==0, 'root x=K')
check(sp.factor(sp.resultant(fn,gn,x)-4096*(K*K-1)**4)==0,
      'duplication resultant')
Fn=b*b*r**4+2*b*b*r*r*s*s-8*a*b*r*s**3+(4*a*a+b*b)*s**4
Gn=4*b*s*(b*r-a*s)*(r*r-s*s)
check(sp.cancel(Fn/Gn-(fn/gn).subs({K:a/b,x:r/s}))==0,
      'homogeneous duplication')
for kval in [F(7,5),F(9,7),F(101,99),F(1,3),F(-2,3)]:
    for xv in [F(0),F(1),F(-1),kval,F(17,19),F(-23,5)]:
        aa,bb=kval.numerator,kval.denominator
        rr,ss=xv.numerator,xv.denominator
        ff=int(Fn.subs({a:aa,b:bb,r:rr,s:ss}))
        gg=int(Gn.subs({a:aa,b:bb,r:rr,s:ss}))
        bound=16*max(abs(aa),bb)**2*max(abs(rr),ss)**4
        check(max(abs(ff),abs(gg))<=bound and (ff,gg)!=(0,0),
              'homogeneous size bound K=%s x=%s'%(kval,xv))

for v in [F(3,2),F(7,3),F(8),F(101,100),F(10001,9999)]:
    n,m=v.numerator,v.denominator
    d=2 if n%2 and m%2 else 1
    triple=((n*n-m*m)//d,2*n*m//d,(n*n+m*m)//d)
    A,B,C=triple
    H=max(n,m)
    check(A*A+B*B==C*C and gcd(gcd(A,B),C)==1,
          'primitive triple v=%s'%v)
    check(H*H<=2*C and C<=2*H*H, 'primitive height bounds v=%s'%v)
    u=-1/v
    check(max(abs(u.numerator),u.denominator)==H,
          'alternate chart exact rational height v=%s'%v)

for jj in [3,4,7,20,100]:
    left,right=F(jj+2,jj+1),F(2*jj+2)
    check(left.numerator==jj+2 and left.denominator==jj+1,
          'left reduced parameter j=%s'%jj)
    check(max(right.numerator,right.denominator)==2*jj+2,
          'right rational height j=%s'%jj)
    check((jj+2)**2+(jj+1)**2==2*jj*jj+6*jj+5,
          'left input hypotenuse j=%s'%jj)
    check((2*jj+2)**2+1==4*jj*jj+8*jj+5,
          'right input hypotenuse j=%s'%jj)

out={'status':'PASS','exact_checks':len(checks),'checks':checks,
     'scope':'Finite algebraic controls only; no numerical asymptotic evidence.'}
body=json.dumps(out,indent=2)+'\n'
if len(sys.argv)>1:
    with open(sys.argv[1],'w') as f:
        f.write(body)
else:
    print(body,end='')

#!/usr/bin/env python3
"""Exact bounded controls for the maximal actual depth-preserving carrier."""
from pathlib import Path
import hashlib
import json
import sympy as sp

x,t,p,y,u,s,tau=sp.symbols('x t p y u s tau')
D=p**3-y**2
source={p:t*(1+x*x*t),y:x*t*t*(1+x*x*t),u:x*x*t}
logchart={p:s*s+tau,y:s*(s*s+tau),u:s*s/tau,x:s/tau}
CHECKS=0


def check(ok,label):
    global CHECKS
    if not ok:raise RuntimeError(label)
    CHECKS+=1


def emit(tag,*values):
    print(tag,json.dumps(values,separators=(',',':')),flush=True)


def src(f):
    return sp.expand(sp.sympify(f).subs(source,simultaneous=True))


def bracket(f,g):
    return sp.expand(sp.diff(f,x)*sp.diff(g,t)-sp.diff(f,t)*sp.diff(g,x))


def val(f,variable):
    f=sp.expand(f)
    if f==0:return None
    return min(int(term.as_powers_dict().get(variable,0)) for term in sp.Add.make_args(f))


for R in [sp.Integer(1),p,y,p*p,p*y,y*y,D,1+p+y+p*y]:
    rp,ry=sp.diff(R,p),sp.diff(R,y)
    images={p:-p*D*ry,y:D*(2*R+p*rp),
            u:(1+u)*y*(4*R+2*p*rp+3*y*ry),
            x:p*(1+2*u)*(2*R+p*rp)+p*y*(2+3*u)*ry}
    S=src(p*p*R)
    for f in (x,u,p,y):
        check(bracket(src(f),S)==src(images[f]),'literal actual generator image')
    for a in range(4):
        for b in range(4-a):
            for c in range(3):
                for e in range(3-c):
                    f=x**a*u**b*p**c*y**e
                    symbolic=sum(sp.diff(f,g)*images[g] for g in (x,u,p,y))
                    check(src(symbolic)==bracket(src(f),S),'whole monomial Leibniz response')
                    check(all(monom[0]+monom[1]<=a+b
                              for monom,_ in sp.Poly(sp.expand(symbolic),x,u,p,y).terms()),
                          'actual depth representative')
                    vv=val(sp.expand(symbolic.subs(logchart,simultaneous=True)),tau)
                    check(vv is None or vv>=-(a+b),'independent Laurent depth')
    emit('CARRIER',str(R),'60 actual depth controls PASS')

# Complete bounded coefficient necessity: monomials p^a*y^b of total degree
# at most six. Both polynomial generator images hold iff a>=2 or constant.
for a in range(7):
    for b in range(7-a):
        S=p**a*y**b
        vp=sp.cancel(-D*sp.diff(S,y)/p)
        vy=sp.cancel(D*sp.diff(S,p)/p)
        admitted=sp.denom(vp)==1 and sp.denom(vy)==1
        check(admitted==(a>=2 or a+b==0),'complete bounded necessity monomial')

# Cusp quotient comparison with p=z^2, y=z^3; actual Euler factors are n+4.
z=sp.Symbol('z')
for a in range(7):
    for b in range(7-a):
        R=p**a*y**b
        ER=4*R+2*p*sp.diff(R,p)+3*y*sp.diff(R,y)
        check(sp.expand(ER.subs({p:z*z,y:z**3})-(4+2*a+3*b)*z**(2*a+3*b))==0,
              'cusp eigenvalue exact')
        check(4+2*a+3*b>0,'no zero Euler eigenvalue')
        DR=sp.expand(D*R)
        response=sp.expand(4*DR+2*p*sp.diff(DR,p)+3*y*sp.diff(DR,y))
        check(sp.rem(sp.Poly(response,y,p),sp.Poly(D,y,p)).is_zero,
              'source intersection sufficiency')

R=sp.Integer(1)
du=4*(1+u)*y
logdu=sp.expand(du.subs(logchart,simultaneous=True))
check(val(logdu,tau)==-1,'depth-preserving is not depth-lowering')
check(sp.expand(tau*logdu).coeff(tau,0)==4*s**5,'sharp hostile depth symbol')
dx=src(2*p*(1+2*u))
check(val(dx,t)==1,'no D-adic gain in maximal carrier')
check(sp.rem(sp.Poly(4*p**3*y,y,p),sp.Poly(D,y,p)).as_expr()!=0,
      'p2 fails every affine source')
check(sp.denom(sp.cancel(2*y*D/p))==p,'D fixed-source carrier fails depth zero')
emit('HOSTILES','p² preserves all depths but fails every affine source',
     'D preserves source H=0 but fails depth zero',
     'non-LND does not decide an individual scalar specialization')
emit('PASS',CHECKS)
emit('SOURCE_SHA256',hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

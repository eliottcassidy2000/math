#!/usr/bin/env python3
import json
import sympy as S
v,t,r,c,a=S.symbols('v t r c a',nonzero=True)
l=2-3*v-t
checks=[]
def check(name,expr):
    if S.cancel(expr)!=0:raise RuntimeError(name)
    checks.append(name)
for w in (1,2,3,5,8):
    x=c/l; y=v/x; z=t/x**2; q=r/x**w
    check('nonzero fiber w%d'%w,x*(2-3*x*y-x*x*z)-c)
    check('recover v w%d'%w,x*y-v)
    check('recover t w%d'%w,x*x*z-t)
    check('recover r w%d'%w,x**w*q-r)
    x0=a;y0=v/a;z0=(2-3*v)/a**2;q0=r/a**w
    check('zero off-hyperplane w%d'%w,x0*(2-3*x0*y0-x0*x0*z0))
    check('zero extra invariant w%d'%w,x0**w*q0-r)
print(json.dumps({'status':'PASS','checks':len(checks),'check_names':checks,'extra_weights':[1,2,3,5,8],'scope':'Exact reconstruction controls for padded seed. General varieties and non-slice conclusion are proved in the companion corollary.'},indent=2,sort_keys=True))

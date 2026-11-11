#!/usr/bin/env python3
"""Finite controls for degree inflation; proofs cover all N in PROOF.md."""
import json
import sympy as s
x,y=s.symbols('x y')
checks=0
def require(v):
    global checks
    checks+=1
    if not v: raise AssertionError('degree control failed')
def subst(f,a,b): return s.expand(f.subs({x:a,y:b}, simultaneous=True))
def jac(p,q): return s.expand(s.diff(p,x)*s.diff(q,y)-s.diff(p,y)*s.diff(q,x))
rows=[]
for a in range(2,6):
    p=x+y**a
    for b in range(2,5):
        q=y+p**b
        require(jac(p,q)==1)
        inv_y=y-x**b
        inv_x=x-inv_y**a
        require(subst(inv_x,p,q)==x)
        require(subst(inv_y,p,q)==y)
        for n in range(2,6):
            r=s.expand(p+q**n)
            require(jac(r,q)==1)
            require(s.Poly(r,x,y).total_degree()==n*a*b)
            require(s.expand(r-q**n-p)==0)
            rows.append({'a':a,'b':b,'N':n,'degree':n*a*b})
print(json.dumps({'status':'PASS','checks':checks,'families':len(rows),
                  'max_degree':max(z['degree'] for z in rows)},indent=2))

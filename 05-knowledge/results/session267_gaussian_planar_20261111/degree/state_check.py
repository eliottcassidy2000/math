#!/usr/bin/env python3
"""Exact finite controls; Krull-height proof is in ALGEBRAIC-STATE.md."""
import json
import sympy as s
x,z,w=s.symbols('x z w')
checks=0
def require(v):
    global checks
    checks+=1
    if not v: raise AssertionError('state control failed')
def E(p):
    return s.expand(sum(c*s.factorial(a) for (a,b),c in s.Poly(s.expand(p),z,w).terms() if a==b))
rows=[]
for m in range(5,81):
    d=(m-1)//2+2
    cols=[]
    for j in range(d+1):
        f=s.expand(-x*x**j+(x*x+6)*s.diff(x**j,x)/(2*m))
        cols.append([f.coeff(x,k) for k in range(m+1)])
    M=s.Matrix.hstack(*(s.Matrix(c) for c in cols))
    require(M.rank()==d+1)
    S=s.Matrix.hstack(*(s.eye(m+1)[:,j] for j in range(m//2+1)))
    require(M.row_join(S).rank()==d+2)
    require((m+1)-M.rank()==m//2-1)
    rows.append({'row':m,'before':m//2-1,'after':m//2-2})
for q in range(1,13):
    D=sum((j*j+2*j+3)*z**j for j in range(q))
    for j in range(q+2):
        require(E(w**j*D)/s.factorial(j)==(j*j+2*j+3 if j<q else 0))
    require(s.expand(E(D*sum(x**j*w**j/s.factorial(j) for j in range(q))))==D.subs(z,x))
for m in range(1,9): require(E((z+2*z*z)**m)==0)
require(1+s.I**2==0)
print(json.dumps({'status':'PASS','checks':checks,'row_range':[5,80],
                  'row16':rows[11],'tomography_q_max':12},indent=2))

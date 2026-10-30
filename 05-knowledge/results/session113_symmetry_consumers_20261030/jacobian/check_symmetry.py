#!/usr/bin/env python3
"""Exact finite tests plus symbolic identities; no floating roots of unity."""
import json
import sympy as s
def need(ok):
    if not ok:
        raise RuntimeError("exact symmetry check failed")
x,y,t=s.symbols('x y t')
def bracket(p,q): return s.expand(s.diff(p,x)*s.diff(q,y)-s.diff(p,y)*s.diff(q,x))
def project(p,r,m,w=1):
    return s.Add(*(c*x**i*y**j for (i,j),c in s.Poly(p,x,y).terms() if (i-w*j-r)%m==0))
def aliases(m,w,D):
    return [(r,i,j) for r in (1,-w) for i in range(D+1) for j in range(D+1-i)
            if (i-w*j-r)%m==0 and i-w*j!=r]
checks=0
for m in range(3,81):
    need(not aliases(m,1,m-2))
    need(aliases(m,1,m-1))
    p,q=x+y**(m-1),y
    need(bracket(p,q)==1)
    need(project(p,1,m)==p and project(q,-1,m)==q)
    checks+=1
weighted=0
for w in range(1,13):
    for D in range(1,18):
        m=w*D+2
        need(not aliases(m,w,D))
        weighted+=1
# Sharp weighted subsequence m=wD+1, where m>w+1.
for w in range(1,13):
    for D in range(2,18):
        m=w*D+1
        p,q=x+y**D,y
        need(bracket(p,q)==1 and project(p,1,m,w)==p and project(q,-w,m,w)==q)
# Twisted one-variable identity and its planar direct check.
A=1+2*t-3*t**2; B=4-t+5*t**3
for w in range(1,8):
    z=x**w*y
    need(s.expand(bracket(x*A.subs(t,z),y*B.subs(t,z))-(A*B+w*t*s.diff(A,t)*B+t*A*s.diff(B,t)).subs(t,z))==0)
# Balanced identity is (tAB)' exactly.
need(s.expand(bracket(x*A.subs(t,x*y),y*B.subs(t,x*y))-s.diff(t*A*B,t).subs(t,x*y))==0)
P=x+y**2; Q=y+(x+y**2)**2
need(bracket(P,Q)==1)
for m in range(6,31):
    U,V=project(P,1,m),project(Q,-1,m)
    need(U==x and s.expand(V-y-2*x*y**2)==0)
    need(bracket(U,V)==1+4*x*y)
    channels={r:project(P,r,m) for r in range(m)}
    channelsq={r:project(Q,r,m) for r in range(m)}
    for k in range(m):
        total=sum(bracket(channels[r],channelsq[(k-r)%m]) for r in range(m))
        need(s.expand(total-(1 if k==0 else 0))==0)
    omitted=sum(bracket(channels[r],channelsq[-r%m]) for r in range(m) if r!=1)
    need(s.expand(omitted)==-4*x*y)
# Hostile controls: origin-Jacobian remains 1, projection still loses Keller away from origin.
need(bracket(x,y+2*x*y**2).subs({x:0,y:0})==1)
need(bracket(x,y+2*x*y**2).subs({x:1,y:1})==5)
result={'status':'PASS','balanced_moduli':checks,'weighted_no_alias_tests':weighted,
        'averaging_moduli':25,'exact_defect':'4*x*y','scope':'finite controls support the accompanying proofs; not a JC proof'}
print(json.dumps(result,indent=2))

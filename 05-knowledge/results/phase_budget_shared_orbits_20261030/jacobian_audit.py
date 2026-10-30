#!/usr/bin/env python3
"""Exact symbolic check of weighted two-shear projection and resonance."""
from pathlib import Path
import json
import sympy as sp
x,y=sp.symbols('x y')
def need(b,msg):
    if not b:raise RuntimeError(msg)
def project(poly,weight,w,m=None):
    return sp.Add(*(c*x**i*y**j for (i,j),c in sp.Poly(poly,x,y).terms()
                   if ((i-w*j)==weight if m is None else (i-w*j-weight)%m==0)))
def jac(P,Q):return sp.expand(sp.diff(P,x)*sp.diff(Q,y)-sp.diff(P,y)*sp.diff(Q,x))
count=cyclic=0
examples=[]
for w in range(1,5):
    for r in range(2,7):
        for s in range(2,9):
            P=x+y**r;Q=y+P**s
            need(jac(P,Q)==1,'original polynomial automorphism')
            pp=project(P,1,w);qq=project(Q,-w,w)
            need(pp==x,'first projected component')
            resonance=(s+w)%(w*r+1)==0
            need((jac(pp,qq)!=1)==resonance,'torus arithmetic resonance iff')
            if resonance:
                j=(s+w)//(w*r+1)
                expected=1+sp.binomial(s,j)*r*j*x**(s-j)*y**(r*j-1)
                need(sp.expand(jac(pp,qq)-expected)==0,'explicit missing determinant channel')
            count+=1
            # Moduli with no alias in first component: exact finite-cyclic criterion.
            for m in range(2,17):
                if (w*r+1)%m==0:continue
                cp=project(P,1,w,m);cq=project(Q,-w,w,m)
                js=[j for j in range(1,s+1) if (s+w-(w*r+1)*j)%m==0]
                need((jac(cp,cq)!=1)==bool(js),'cyclic resonance iff under no P alias')
                cyclic+=1
            if (w,r,s) in ((1,2,2),(1,2,5),(1,2,3),(2,2,3)):
                examples.append(dict(w=w,r=r,s=s,projection=str(qq),determinant=str(jac(pp,qq))))
out=dict(status='exact symbolic controls; elementary all-parameter proof in GEOMETRY_JC.md',
         torus_cases=count,cyclic_cases=cyclic,examples=examples)
Path(__file__).with_name('jacobian_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

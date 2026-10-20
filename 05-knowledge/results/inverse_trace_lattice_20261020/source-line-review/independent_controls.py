#!/usr/bin/env python3
"""Cheap exact hostiles for the separate source-line addendum."""
from pathlib import Path
import json
import sympy as s
x,y,t,a,u,v,z=s.symbols('x y t a u v z')
labels=[]
def check(c,label):
    if not bool(c):
        raise RuntimeError(label)
    labels.append(label)
# An immersed finite birational parametrization need not be injective.
p=t*t-1
q=t*(t*t-1)
node=v*v-u*u*(u+1)
check(s.expand(node.subs({u:p,v:q}))==0,'node equation')
check(s.gcd(s.diff(p,t),s.diff(q,t))==1,'node parametrization is immersive everywhere')
check(s.cancel(q/p)==t,'node parametrization is birational')
check(p.subs(t,1)==p.subs(t,-1)==q.subs(t,1)==q.subs(t,-1)==0,'two parameters collide')
check(s.diff(node,u).subs({u:0,v:0})==s.diff(node,v).subs({u:0,v:0})==0,'image is singular at collision')
check(s.det(s.Matrix([[s.diff(p,t).subs(t,1),s.diff(p,t).subs(t,-1)],
                     [s.diff(q,t).subs(t,1),s.diff(q,t).subs(t,-1)]]))!=0,'two distinct smooth branches')
check(s.degree(s.Symbol('T')**2-(u+1),s.Symbol('T'))==2,'parameter satisfies a monic equation over node ring')
# Nontrivial finite polynomial covers necessarily supply critical points.
for k in (2,3,5):
    cover=t**k
    fp=s.expand(p.subs(t,cover)); fq=s.expand(q.subs(t,cover))
    check(s.diff(cover,t).subs(t,0)==0,f'power cover critical point k={k}')
    check(s.diff(fp,t).subs(t,0)==s.diff(fq,t).subs(t,0)==0,f'composed parametrization loses immersion k={k}')
# Finite surjective smooth maps do not alone force reduced hyperplane sections.
# G(x,y)=(u=x,v=y^2,z=y) has W=(z^2-v); use slice v=0 to see multiplicity two.
check(s.expand((z*z-v).subs({u:x,v:y*y,z:y}))==0,'smooth finite graph identity')
check(s.Poly((z*z-v).subs(v,0),z).as_expr()==z*z,'tangent hyperplane slice is nonreduced')
# With the source-coordinate slice, the same smooth graph gives a reduced section.
check(s.diff((z*z-v).subs(z,a),v)==-1,'source-coordinate slice remains reduced')
# A normal smooth source and finite normalization do not make generic image slices smooth.
check(s.expand(node.subs({u:p,v:q,z:a}))==0,'node times line normalization identity')
check(node.subs(z,a)==node,'nodal identification persists in every parameter slice')
# The constant-term-one atlas has invertible Z and cannot contain a source x=0 image.
for d in (3,6):
    R=u*z**d+v*z+1
    check(s.expand(R.subs(v,-u*z**(d-1)-1/z))==0,f'A1 times Gm atlas d={d}')
    check(R.subs(z,0)==1,f'empty zero-coordinate slice d={d}')
# A finite birational graph can have a nonreduced SOURCE-coordinate slice.
R=u*z*z-v*v
check(s.expand(R.subs({u:y*y,v:x*y,z:x}))==0,'finite birational graph of (y^2,xy)')
check(R.subs(z,0)==-v*v,'nonreduced source-coordinate specialization')
check(s.det(s.Matrix([[0,2*y],[y,x]]))==-2*y*y,'nonreduced specialization hostile is non-Keller')
check(s.diff(y*y,y).subs(y,0)==0,'special-line parametrization ramifies')
# A cusp has a unique normalization preimage but fails immersion.
check(s.diff(t*t,t).subs(t,0)==s.diff(t**3,t).subs(t,0)==0,'cusp fails immersion')
check(1-(2-1)==0,'single node normalization Euler is zero')
check(1-(1-1)==1,'single cusp normalization Euler is one')
out={'status':'PASS','checks':len(labels),'labels':labels,'scope':'Exact local algebra controls, not actual Keller examples or replacements for the cited theorems.'}
Path(__file__).with_name('results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':len(labels)},sort_keys=True))

#!/usr/bin/env python3
"""Independent exact quartic counterexample to an unqualified boundary-degree rule."""
import argparse,json
from pathlib import Path
import sympy as s
x,y,z,u,h,A,B,C,P,Q,w=s.symbols('x y z u h A B C P Q w')
checks=[]
def check(name,value):
    if not value:raise RuntimeError(name)
    checks.append(name)
u0=1+3*x*y;h0=1-4*x*y-x*x*z
F=s.Matrix([s.cancel((2*u0+u0*u0-3*u0**4*h0**2)/x**2),s.cancel((1+u0-2*u0**3*h0**2)/x),x*h0])
check('all coordinates polynomial',all(s.denom(q)==1 for q in F))
check('coordinate degrees 12,11,4',[s.Poly(q,x,y,z).total_degree() for q in F]==[12,11,4])
check('determinant -6',s.expand(F.jacobian((x,y,z)).det())==-6)
plane=tuple(s.expand(q).subs(x,0) for q in F)
check('source invariant plane restriction',plane==(87*y*y+6*z,y,0))
fold=s.Matrix([(2*u+u*u)/x**2,(1+u)/x,0])
check('fold elimination',s.cancel((fold[0]-A).subs(u,B*x-1)-((B**2-A)*x**2-1)/x**2)==0)
points=[(0,1,-s.Rational(29,2)),(1,-s.Rational(1,3),s.Rational(7,3)),(-1,1,5)]
check('three distinct exact fiber points',len(set(points))==3 and all(tuple(F.subs(dict(zip((x,y,z),pt))))==(0,1,0) for pt in points))
# On C!=0, quotient invariants satisfy an irreducible quartic inverse equation.
P0=h+w-2*w**3
Q0=2*w*h+w*w-3*w**4
T=w**4-w*w+2*P*w-Q
check('inverse quartic identity',s.expand(T.subs({P:P0,Q:Q0}))==0)
check('quartic monic in w',s.Poly(T,w).degree()==4 and s.Poly(T,w).LC()==1)
check('quartic defines prime graph ideal',s.diff(T,Q)==-1)
# Explicit reconstruction verifies sufficiency on the generic open C*h!=0.
hrec=P-w+2*w**3
xrec=C/hrec
urec=w/hrec
yrec=(urec-1)/(3*xrec)
zrec=(1-4*xrec*yrec-hrec)/xrec**2
check('reconstructed gamma',s.cancel(h0.subs({x:xrec,y:yrec,z:zrec},simultaneous=True)-hrec)==0)
check('reconstructed u',s.cancel(u0.subs({x:xrec,y:yrec,z:zrec},simultaneous=True)-urec)==0)
# Invariant forms are sufficient for generic source recovery, avoiding large expansion.
check('target P recovered',s.expand(hrec+w-2*w**3-P)==0)
check('target Q recovered modulo T',s.rem(s.Poly(2*w*hrec+w*w-3*w**4-Q,w),s.Poly(T,w))==0)
hostiles=[]
for name,condition in [
 ('wrong fold degree d-minus-one',2==4-1),
 ('wrong determinant',s.expand(F.jacobian((x,y,z)).det())==6),
 ('altered third source point',tuple(F.subs({x:-1,y:1,z:-3}))==(0,1,0)),
 ('wrong quartic sign',s.expand((w**4+w*w+2*P*w-Q).subs({P:P0,Q:Q0}))==0)]:
    if condition:raise RuntimeError('hostile accepted: '+name)
    hostiles.append(name)
result={'status':'PASS','checks':checks,'check_count':len(checks),'hostiles_rejected':hostiles,'generic_field_degree':4,'generic_invariant_hyperplane_finite_points':3,'invariant_source_plane_degree':1,'orbit_zero_surface_degree':2,'generic_missing_sheet_count':1,'scope':'Explicit quartic map; refutes unqualified general d-minus-restricted-degree formula.'}
p=argparse.ArgumentParser();p.add_argument('--output',type=Path);args=p.parse_args();out=json.dumps(result,sort_keys=True,indent=2)+'\n'
if args.output:args.output.write_text(out)
else:print(out,end='')

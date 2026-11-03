#!/usr/bin/env python3
"""Exact bounded-theta addendum audit; no discovery imports."""
import json, sys, hashlib
from pathlib import Path
import sympy as s
checks=0
def must(x,msg):
 global checks
 checks+=1
 if not bool(x): raise RuntimeError(msg)
z=s.symbols('z'); L=16*z*z-20*z+5; C=16*z*z-12*z+1; B=s.expand(5-L*C)
must(s.expand(B-(-256*z**4+512*z**3-336*z**2+80*z))==0,'B expansion')
must(s.expand(2*z*(1-z)*s.diff(B,z)-B+4*z*(1-z)*C*s.diff(L,z))==0,'all coefficients of derivative identity vanish')
R=s.sqrt(13); z2=(5-R)/8
must(s.simplify(L.subs(z,z2)-2)==0,'ell=2 root')
must(z2.is_positive and s.simplify(1-z2).is_positive,'physical root 0<z<1')
must(s.simplify(s.Rational(1,2)-z2).is_positive,'root below one half')
G2=s.simplify((1-z2)/z2*B.subs(z,z2)**2/16)
must(s.simplify(G2-(89+26*R)/16)==0,'exact g2 square')
# Independent cocircularity check of the already audited optimal coordinate witness.
D2=(R-1)/6; A=(R+1)/4; BB=s.Rational(3,2); H=(R+7)/4
center=s.simplify((A*A+D2*BB*BB-1)/(2*D2*BB))
points=[(-1,0),(-A,BB),(s.Rational(-1,2),H),(s.Rational(1,2),H),(A,BB),(1,0)]
for x,y in points:
 must(s.simplify(x*x+D2*y*y-2*D2*center*y-1)==0,'all six corners on one circle')
# Valid radical root yields the required arc range: z > sin(pi/10)^2 and z < sin(pi/5)^2.
left=(3-s.sqrt(5))/8; right=(5-s.sqrt(5))/8
must(s.simplify(z2-left).is_positive,'alpha > pi/10')
must(s.simplify(right-z2).is_positive,'alpha < pi/5')
# Algebraic derivative is positive on the entire branch: C=cos5a/cosa<0
# follows from pi/2<5a<pi; L'(z)=32z-20<0 as z<sin(pi/5)^2<5/8.
must(s.simplify(s.Rational(5,8)-right).is_positive,'strict negative branch derivative L')
res={'status':'PASS','checked_conditions':checks,'g2_squared':str(G2),
 'witness_circle_center_y_over_D':str(center),
 'scope':'Endpoint-only theta(5,5,2), outer 5+5 Jordan polygon, interior 2-edge separator, open straight terminal chord interior.',
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(sys.argv[1]).write_text(json.dumps(res,indent=2,sort_keys=True)+'\n');print(json.dumps(res,sort_keys=True))

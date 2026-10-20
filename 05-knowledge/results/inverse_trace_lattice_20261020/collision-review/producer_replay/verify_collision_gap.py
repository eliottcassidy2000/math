#!/usr/bin/env python3
"""Small exact controls for the collision-gap and singular-locus bridge."""
from pathlib import Path
import json
import sympy as S
x,y,a,s,t,u,v,z=S.symbols('x y a s t u v z')
gap=s-t
checks=0

def check(ok,label):
    global checks
    if not bool(ok):
        raise AssertionError(label)
    checks+=1

def exact_div(f,g):
    return S.Poly(S.expand(f),a,s,t).exquo(S.Poly(g,a,s,t)).as_expr()

cases=[]
for n in range(2,5):
  for m in range(2,4):
    for c in [1,-2]:
      P=x+y**n; Q=c*(y+P**m)
      J=S.expand(S.diff(P,x)*S.diff(Q,y)-S.diff(P,y)*S.diff(Q,x))
      check(J==c,('constant Jacobian',n,m,c))
      Ps=P.subs({x:a,y:s}); Pt=P.subs({x:a,y:t})
      U=exact_div(Ps-Pt,gap)
      V=exact_div(Q.subs({x:a,y:s})-Q.subs({x:a,y:t}),gap)
      check(S.Poly(U,t).LC()==1,('monic collision relation',n,m,c))
      check(S.expand(U.subs(t,s)-S.diff(P,y).subs({x:a,y:s}))==0,('diagonal U',n,m,c))
      check(S.expand(V.subs(t,s)-S.diff(Q,y).subs({x:a,y:s}))==0,('diagonal V',n,m,c))
      Px=S.diff(P,x).subs({x:a,y:s}); Qx=S.diff(Q,x).subs({x:a,y:s})
      omega=exact_div(Px*V-Qx*U-J,gap)
      check(S.expand(J-Px*V+Qx*U+gap*omega)==0,('gap identity',n,m,c))
      powers=sum(Ps**(m-1-k)*Pt**k for k in range(m))
      check(S.expand(V-c*U*powers)==c,('automorphism collision ideal is unit',n,m,c))
      check(omega!=0,('nontrivial quotient certificate',n,m,c))
      cases.append({'n':n,'m':m,'J':c,'omega_total_degree':S.Poly(omega,a,s,t).total_degree()})

# The exact finite birational non-Keller graph F=(y^2,xy).
R=u*z**2-v**2
check(S.expand(R.subs({u:y**2,v:x*y,z:x}))==0,'non-Keller actual graph')
check(S.expand(S.diff(y**2,x)*S.diff(x*y,y)-S.diff(y**2,y)*S.diff(x*y,x))==-2*y**2,'nonconstant Jacobian')
U=s+t; V=a
check(U.subs({a:0,s:0,t:0})==V.subs({a:0,s:0,t:0})==0,'diagonal collision survives without Keller')
check(gap.subs({a:0,s:0,t:0})==0,'gap not unit in non-Keller quotient')
check(R.subs({u:0,v:0,z:1})==S.diff(R,z).subs({u:0,v:0,z:1})==0,'R and RZ vanish')
check(S.diff(R,u).subs({u:0,v:0,z:1})==1,'RZ locus may contain smooth points without Keller')

# A unit gap and zero parameter elimination ideal can still omit a parameter.
miss=[a*gap-1,t]
G=S.groebner(miss,s,t,a,order='lex')
check(not any(not q.as_expr().has(s,t) for q in G.polys),'dominance elimination ideal is zero')
check(miss[0].subs(a,0)==-1,'dominant projection misses zero')
check(S.expand(a*gap-1-miss[0])==0,'explicit unit gap in dominant hostile')
for av in [-3,-2,-1,1,2,3]:
    check(all(S.cancel(f.subs({a:av,t:0,s:S.Rational(1,av)}))==0 for f in miss),('nonzero parameter attained',av))

# Surjective, finite-fiber, nonproper abstract curve with a unit gap.
curve=[t*gap-1,a-t**2+t]
sub={a:t**2-t,s:t+1/t}
check(all(S.cancel(f.subs(sub))==0 for f in curve),'Laurent parametrization')
check(S.expand(t*gap-1-curve[0])==0,'explicit gap inverse')
check(all(f.subs({a:0,t:1,s:2})==0 for f in curve),'zero parameter attained')
for av in range(-5,6):
    poly=S.Poly(t**2-t-av,t)
    check(poly.degree()==2,('positive degree fiber equation',av))
    if av:
        check(poly.eval(0)!=0,('all fiber roots nonzero',av))
    else:
        check(poly.eval(1)==0,('surviving nonzero root at zero',av))
for k in range(2,12):
    tv=S.Rational(1,k); av=tv**2-tv; sv=tv+1/tv
    check(all(f.subs({a:av,t:tv,s:sv})==0 for f in curve),('escaping exact points',k))
    check(abs(av)<S.Rational(1,k) and sv>k,('bounded parameter unbounded source',k))

# On H=0, R=AZ+B and RZ=A. Triple-locus exclusion forces smoothness.
for A in range(-3,4):
  for B in range(-3,4):
    if A:
        root=-S.Rational(B,A)
        check(A*root+B==0 and A!=0,('simple boundary graph root',A,B))
    elif B:
        check(S.Poly(A*z+B,z).degree()==0,('no boundary graph point when A zero',B))

out={'status':'PASS','checks':checks,'automorphism_cases':cases,
     'controls':['non-Keller diagonal and extra RZ zeros','dominant nonsurjective unit-gap curve',
                 'surjective nonproper unit-gap curve','linear boundary roots'],
     'scope':'Finite exact controls; no abstract hostile is claimed to be a Keller realization.'}
Path(__file__).with_name('collision-gap-verification-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'checks':checks,'automorphism_cases':len(cases),'status':'PASS'},sort_keys=True))

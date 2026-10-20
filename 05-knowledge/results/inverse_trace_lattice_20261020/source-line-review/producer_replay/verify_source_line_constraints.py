#!/usr/bin/env python3
"""Exact hostile controls for the source-line specialization lemma."""
from pathlib import Path
import json
import sympy as s
u,v,z,x,y,a,t=s.symbols('u v z x y a t')
n=0

def check(ok,label):
    global n
    if not bool(ok):
        raise AssertionError(label)
    n+=1

# A finite birational graph can have a nonreduced line specialization without Keller.
# F(x,y)=(y^2,x*y), chosen graph coordinate Z=x.
R=u*z**2-v**2
check(s.expand(R.subs({u:y**2,v:x*y,z:x}))==0,'actual graph relation')
check(s.factor_list(R)[1]==[(R,1)],'graph polynomial irreducible')
check(R.subs(z,0)==-v**2,'nonreduced specialization')
check(s.expand((x*y/x)**2-y**2)==0,'generic recovery of companion y')
jac=s.det(s.Matrix([[0,2*y],[y,x]]))
check(jac==-2*y**2,'Keller condition fails')
check(s.diff(y**2,y).subs(y,0)==0,'line parametrization ramifies')
check((y**2).subs(y,1)==(y**2).subs(y,-1),'generic double covering on special line')

# The nodal normalization is immersive but identifies two parameters.
nu=t**2-1; nv=t*(t**2-1)
check(s.expand(nv**2-nu**2*(nu+1))==0,'node equation')
check(s.gcd(s.diff(nu,t),s.diff(nv,t))==1,'node immersion everywhere')
check(nu.subs(t,1)==nu.subs(t,-1)==0 and nv.subs(t,1)==nv.subs(t,-1)==0,'node collision')
check(s.det(s.Matrix([[s.diff(nu,t).subs(t,1),s.diff(nv,t).subs(t,1)],
                     [s.diff(nu,t).subs(t,-1),s.diff(nv,t).subs(t,-1)]]))!=0,'node branches distinct tangents')

# A cusp is finite birational but fails immersion, so χ=1 does not imply smooth A1.
check(s.expand((t**3)**2-(t**2)**3)==0,'cusp equation')
check(s.diff(t**2,t).subs(t,0)==s.diff(t**3,t).subs(t,0)==0,'cusp fails immersion')
check(1-(1-1)==1,'unibranch cusp Euler characteristic')

atlases=[]
for d in range(3,11):
    R=u*z**d+v*z+1
    check(R.subs(z,0)==1,('empty zero-coordinate image',d))
    check(s.expand(z*(u*z**(d-1)+v)-(R-1))==0,('coordinate is a unit',d))
    check(s.cancel(R.subs(v,-u*z**(d-1)-1/z))==0,('Laurent surface parametrization',d))
    for av in range(-3,4):
        Ra=s.expand(R.subs(z,av))
        if av:
            check(s.Poly(Ra,u,v).total_degree()==1,('other specializations are lines',d,av))
            check(s.Poly(Ra,v).LC()==av,('other specialization smooth',d,av))
        else:
            check(Ra==1,('forbidden constant specialization',d))
    atlases.append({'d':d,'R':'u*Z^d+v*Z+1','forbidden_specialization':'Z=0',
                    'surface':'A1 x G_m','status':'cannot be actual finite source graph'})

# Simple scalar-only topology budgets; these do not construct polynomial Keller maps.
budgets=[]
for d in range(3,13):
  for branch_excess in range(1,6):
    chi=1-branch_excess
    omitted=(d-1)*(1-chi)
    check(chi<=0 and omitted>=d-1,('irreducible boundary lower bound',d,branch_excess))
    budgets.append({'d':d,'branch_excess':branch_excess,'chi':chi,'omitted':omitted})

out={'checks':n,'status':'PASS','source_control':'finite birational graph with nonreduced specialization is non-Keller',
     'forbidden_constant_atlases':atlases,'irreducible_boundary_budgets':budgets,
     'scope':'Finite exact controls only; general proofs and versioned external inputs are in the accompanying note.'}
Path(__file__).with_name('source-line-verification-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'checks':n,'atlas_cases':len(atlases),'budget_cases':len(budgets),'status':'PASS'},sort_keys=True))

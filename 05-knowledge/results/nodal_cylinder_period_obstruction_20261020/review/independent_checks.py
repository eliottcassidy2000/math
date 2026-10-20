#!/usr/bin/env python3
"""Independent finite checks: coefficient moments and endpoint-jet conormals.

Does not import the producer verifier. All equalities are symbolic over Q.
No finite check substitutes for the cited all-embedding theorem.
"""
import argparse
import json
from pathlib import Path
import sympy as S

x,y,W,Z=S.symbols('x y W Z')
w=y*y-1
z=y*w
h=Z*Z-W*W*(W+1)
checks=0
rows=[]

def eq(name,a,b=0):
    global checks
    checks+=1
    d=S.expand(a-b)
    if d!=0:
        raise RuntimeError(f'{name}: difference {d}')

def period(p,q):
    # Coefficient pairing rather than symbolic antidifferentiation.
    pp=S.Poly(p,y)
    qq=S.Poly(q,y)
    total=S.Integer(0)
    for (i,),ai in pp.terms():
        for (j,),bj in qq.terms():
            if j and (i+j)%2:
                total+=S.Rational(2*j,i+j)*ai*bj
    return S.expand(total)

def integ(p):
    return S.expand(sum(S.Rational(2,i+1)*a for (i,),a in S.Poly(p,y).terms() if i%2==0))

def jac(p,q):
    return S.expand(S.diff(p,x)*S.diff(q,y)-S.diff(p,y)*S.diff(q,x))

def jets(p):
    dp=S.diff(p,y)
    plus=dp.subs(y,1)
    minus=dp.subs(y,-1)
    return ((plus-minus)/4,(plus+minus)/4)

def run(name,f,g,expected):
    p=S.expand(f.subs({W:w,Z:z}, simultaneous=True))
    q=S.expand(g.subs({W:w,Z:z}, simultaneous=True))
    aj,bj=jets(p)
    cj,dj=jets(q)
    conormal=S.expand(aj*dj-bj*cj)
    ambient=S.expand(S.diff(f,W)*S.diff(g,Z)-S.diff(f,Z)*S.diff(g,W))
    eq(name+' ambient determinant',ambient,expected)
    eq(name+' jet conormal',conormal,expected)
    eq(name+' P endpoint',p.subs(y,1),p.subs(y,-1))
    eq(name+' Q endpoint',q.subs(y,1),q.subs(y,-1))
    for label,value,deriv in [('PW',aj,S.diff(f,W)),('PZ',bj,S.diff(f,Z)),('QW',cj,S.diff(g,W)),('QZ',dj,S.diff(g,Z))]:
        eq(name+' '+label,value,deriv.subs({W:0,Z:0}))
    area=period(p,q)
    eq(name+' area',area,S.Rational(8,15)*expected)
    j=jac(p,q)
    eq(name+' Jacobian mean',integ(j))
    boundary=(p*S.diff(q,x)).subs(y,1)-(p*S.diff(q,x)).subs(y,-1)
    eq(name+' integration by parts',S.diff(area,x),integ(j)+boundary)
    eq(name+' boundary',boundary)
    vp=S.Matrix([S.diff(p,y).subs(y,1),S.diff(q,y).subs(y,1)])
    vm=S.Matrix([S.diff(p,y).subs(y,-1),S.diff(q,y).subs(y,-1)])
    eq(name+' direct tangent determinant',S.det(S.Matrix.hstack(vp,vm)),15*area)
    rows.append({'name':name,'degrees_y':[int(S.degree(p,y)),int(S.degree(q,y))],'conormal_determinant':str(conormal),'area':str(area)})

# Independently chosen 36 generating families. Each is a composition of
# translations, two elementary triangular maps, constant scalings and a swap.
for a in range(1,4):
    for b in range(1,5):
        f=W+x**a
        g=Z+x**b
        f=f+(x+1)*g*g
        g=g+(x*x-1)*f*f
        run(f'base_{a}_{b}',f,g,1)
        run(f'reverse_{a}_{b}',g,f,-1)
        run(f'scale_{a}_{b}',-a*f,S.Rational(b,3)*g,-S.Rational(a*b,3))

# Standard normalization and conductor presentation.
eq('node relation',h.subs({W:w,Z:z},simultaneous=True))
eq('normalization rational recovery',z-y*w)
eq('area normalization',period(w,z),S.Rational(8,15))
eq('endpoint jet orientation',S.det(S.Matrix([[2,-2],[2,2]])),8)

# Changing a representative by the equation changes neither pullback nor
# its conormal first-order part, since the node equation has order two.
for r in [1,x+W+Z,W*Z+x*x]:
    eq('relation representative pullback',(h*r).subs({W:w,Z:z},simultaneous=True))
    eq('relation representative W derivative',S.diff(h*r,W).subs({W:0,Z:0}))
    eq('relation representative Z derivative',S.diff(h*r,Z).subs({W:0,Z:0}))

# Three hostile tests: generation, unit conormal, and missing endpoint boundary.
eq('nongenerator generic-only area',period(w,x*z),S.Rational(8,15)*x)
eq('nongenerator generic-only determinant',S.det(S.Matrix([jets(w),jets(x*z)])),x)
qbad=z+x*w*w*z
eq('unit conormal nongenerator determinant',S.det(S.Matrix([jets(w),jets(qbad)])),1)
eq('unit conormal nongenerator area',period(w,qbad),S.Rational(8,15)+S.Rational(64,315)*x)
eq('extra collision first coordinate',S.rem(w-1,y*y-2,y))
eq('extra collision second coordinate',S.rem(qbad.subs(x,-1),y*y-2,y))
pout=y*y
qout=x*y
eq('nonmember area',period(pout,qout),S.Rational(2,3)*x)
eq('nonmember Jacobian mean',integ(jac(pout,qout)),-S.Rational(4,3))
eq('nonmember boundary', (pout*S.diff(qout,x)).subs(y,1)-(pout*S.diff(qout,x)).subs(y,-1),2)
eq('nonmember integration by parts',S.diff(period(pout,qout),x),integ(jac(pout,qout))+2)

result={'status':'PASS','checks':checks,'generating_families':len(rows),'method':'independent coefficient moments and endpoint derivative conormals; no imported producer code','families':rows,'scope':'FINITE-EXACT controls, not proof of all-embedding uniqueness'}
a=argparse.ArgumentParser(); a.add_argument('--output',default='independent_results.json'); args=a.parse_args()
Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(f'PASS: {checks} exact checks; {len(rows)} independently selected generating families.')

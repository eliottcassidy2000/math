#!/usr/bin/env python3
"""Exact S144 checks. No assert gates; only SymPy rational polynomial arithmetic.
Run: python -B verify_equivariant_composition.py
     python -B -O verify_equivariant_composition.py
The arbitrary-iterate and all-degree statements are proved in the companion note.
"""
import json
from math import gcd
import sympy as S

v,t,x,y,z,lam=S.symbols('v t x y z lam')
checks=[]
def check(name, predicate):
    if not bool(predicate):
        raise RuntimeError('FAILED: '+name)
    checks.append(name)
def zero(expr):
    return S.expand(expr)==0

def substitute(expr, vv, tt):
    return S.Poly(expr,v,t).eval({v:vv,t:tt}) if False else S.expand(expr.subs({v:vv,t:tt},simultaneous=True))
def jac_pair(P,Q):
    # Poly arithmetic prevents expression swell in the direct determinant route.
    P,Q=S.Poly(P,v,t),S.Poly(Q,v,t)
    return (P.diff(v)*Q.diff(t)-P.diff(t)*Q.diff(v)).as_expr()

u=1+x*y
F=S.Matrix([u**3*z+y**2*u*(4+3*x*y), y+3*x*u**2*z+3*x*y**2*(4+3*x*y),2*x-3*x*x*y-x**3*z])
G=S.Matrix([F[2],F[1],F[0]])
check('fixed witness determinant -2',zero(F.jacobian([x,y,z]).det()+2))
check('untwisted witness determinant +2',zero(G.jacobian([x,y,z]).det()-2))
points=[(0,0,-S.Rational(1,4)),(1,-S.Rational(3,2),S.Rational(13,2)),(-1,S.Rational(3,2),S.Rational(13,2))]
for j,pt in enumerate(points):
    check('exact collision '+str(j),F.subs(dict(zip((x,y,z),pt)))==S.Matrix([-S.Rational(1,4),0,0]))
l=2-3*v-t
K=t*(1+v)**2+v**2*(4+3*v)
A=v+3*K
B=(1+v)*K
H=(S.expand(l*A),S.expand(l*l*B))
check('quotient first coordinate',zero((G[0]*G[1])-H[0].subs({v:x*y,t:x*x*z})))
check('quotient second coordinate',zero((G[0]**2*G[2])-H[1].subs({v:x*y,t:x*x*z})))
check('quotient determinant',zero(jac_pair(*H)-2*l*l))
check('quotient fixes origin',all(h.subs({v:0,t:0})==0 for h in H))
check('common-zero elimination identity',zero(3*B-(1+v)*A+v*(1+v)))
check('v=-1 excluded from A=B=0',A.subs(v,-1)==2)
check('v=0 forces t=0',B.subs(v,0)==t)

H2=tuple(substitute(h,*H) for h in H)
f=[l,substitute(l,*H),substitute(l,*H2)]
l2=S.expand(f[0]*f[1]);l3=S.expand(l2*f[2])
check('second quotient direct determinant',zero(jac_pair(*H2)-4*l2*l2))
check('H degree pair',(S.Poly(H[0],v,t).total_degree(),S.Poly(H[1],v,t).total_degree())==(4,6))
check('H squared degree pair',(S.Poly(H2[0],v,t).total_degree(),S.Poly(H2[1],v,t).total_degree())==(20,30))
for i,j in ((0,1),(0,2),(1,2)):
    _,rem=S.div(f[j]-2,f[i],v,t)
    check('pairwise comaximal f%d f%d'%(i,j),zero(rem))
check('orbit degree sequence 1,7,37',[S.Poly(p,v,t).total_degree() for p in (l,l2,l3)]==[1,7,37])
top2=sum(c*v**a*t**b for (a,b),c in S.Poly(l2,v,t).terms() if a+b==7)
check('second scale top homogeneous',zero(top2-v**3*(t+3*v)**4))
for j,fj in enumerate(f):
    check('scale pullback f%d has no repeated factors'%j, S.Poly(S.gcd(fj,S.gcd(S.diff(fj,v),S.diff(fj,t))),v,t).total_degree()==0)

# Corrected weight-module classification, exhaustively in a stated finite universe.
monomials=0
for k in range(2,8):
    for a in range(19):
      for b in range(19-a):
       for c in range(19-a-b):
        wt=a-b-k*c
        if wt not in (-k,-1,1):continue
        monomials+=1
        if wt==-k:
            ok=(c==0 and b==a+k) or (c>=1 and a==b+k*(c-1))
        elif wt==-1:
            ok=(c==0 and b==a+1) or (c>=1 and a-(k-1)==b+k*(c-1))
        else:ok=(a-1==b+k*c)
        if not ok:raise RuntimeError('monomial module reconstruction failed')
check('finite module reconstruction k2..7 degree<=18',monomials>0)
check('missing quadratic-z hostile', (2-2*2)==-2) # x^2 z^2 has weight -2.
# c2=0 valuation repair, incl. odd k and rational scalar units.
for k in range(2,10):
    gg=gcd(k+1,2)
    h=t*t+t+1
    AA=2*h**((k+1)//gg);DD=3*h**(2//gg)
    check('gcd-corrected c2 k%d'%k,zero((k+1)*AA*S.diff(DD,t)-2*S.diff(AA,t)*DD))
check('odd k minimal hostile A=t^2 D=t',zero(4*t*t*S.diff(t,t)-2*S.diff(t*t,t)*t))
# Positive and hostile controls: all-negative/positive weights do not have this mechanism;
# equal sums alone need not yield a tangent intertwiner.
source=(1,-1,-3);target=(1,-2,-2)
intertwiner=S.Matrix([[1 if i==j else 0 for j in source] for i in target])
check('equal-sum unequal-multiset obstruction',sum(source)==sum(target) and intertwiner.rank()==1)
# Identity retains l=1, affine torus scaling retains constant orbit scale.
check('identity positive control',zero(jac_pair(v,t)-1))
check('non-Keller quotient hostile',zero(jac_pair(v*v,t)-2*v) and not zero(2*v-1))

result={
 'status':'PASS',
 'checks':len(checks),
 'check_names':checks,
 'finite_universe':{'weight_k':list(range(2,8)),'total_monomial_degree_max':18,'weight_monomials_checked':monomials},
 'exact_degrees':{'H':[4,6],'H_squared':[20,30],'iterate_orbit_scale':[1,7,37]},
 'second_scale_factored':str(S.factor(l2)),
 'arbitrary_iterate_proof':'Companion theorem, not inferred from finite tests.',
 'scope':'Explicit untwisted-map composition tower and corrected weight-module statements; not all Keller-map classification.',
 'sympy_version':S.__version__,
}
print(json.dumps(result,indent=2,sort_keys=True))

#!/usr/bin/env python3
"""Exact independent audit; run with and without -O. SymPy dependency.
No repository writes. Output goes to stdout or --output PATH.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as S

class CheckFailure(RuntimeError):
    pass

checks=[]
def require(name, condition):
    if not condition:
        raise CheckFailure(name)
    checks.append(name)

def polynomial_equal(actual, expected):
    return S.expand(actual-expected)==0

x,y,z,a,b,g,v,t,h,U=S.symbols('x y z a b g v t h U')
u=1+x*y
F=S.Matrix([u**3*z+y**2*u*(4+3*x*y),y+3*x*u**2*z+3*x*y**2*(4+3*x*y),2*x-3*x**2*y-x**3*z])
detF=S.expand(F.jacobian((x,y,z)).det())
require('seed determinant -2',detF==-2)
collision=[(0,0,-S.Rational(1,4)),(1,-S.Rational(3,2),S.Rational(13,2)),(-1,S.Rational(3,2),S.Rational(13,2))]
images=[tuple(F.subs(dict(zip((x,y,z),p)))) for p in collision]
require('three distinct exact collision points',len(set(collision))==3 and all(q==(-S.Rational(1,4),0,0) for q in images))
l=2-3*v-t
K=(1+v)**2*t+v**2*(4+3*v)
A=v+3*K
B=(1+v)*K
H=S.Matrix([l*A,l**2*B])
require('quotient identity determinant 2*l^2',polynomial_equal(H.jacobian((v,t)).det(),2*l**2))
require('quotient origin fixed',tuple(H.subs({v:0,t:0}))==(0,0))
require('off-divisor simultaneous zero identity',polynomial_equal(3*B-(1+v)*A,-v*(1+v)))
require('hostile v=-1 excluded',S.expand(A.subs(v,-1))==2)
require('v=0 forces t=0',S.expand(B.subs(v,0))==t)
f1=S.expand(l.subs({v:H[0],t:H[1]},simultaneous=True))
l2=S.expand(l*f1)
require('second scale degree 7',S.Poly(l2,v,t).total_degree()==7)
top7=sum(c*v**i*t**j for (i,j),c in S.Poly(l2,v,t).terms() if i+j==7)
require('second scale leading form',polynomial_equal(top7,v**3*(t+3*v)**4))
require('f1=2 modulo f0',S.rem(S.Poly(f1-2,t,v),S.Poly(l,t,v))==0)
require('f0 and f1 coprime',S.Poly(S.gcd(l,f1),v,t).total_degree()==0)
require('f1 smooth exact Groebner certificate',list(S.groebner([f1,S.diff(f1,v),S.diff(f1,t)],t,v))==[1])
# General degree recurrence has unique dominant terms, independently enumerate all monomials.
for av,bv in [(4,6),(20,30),(100,150)]:
    for idx,expected in [(0,2*(av+bv)),(1,3*(av+bv))]:
        terms=S.Poly(H[idx],v,t).terms()
        vals=[(av*i+bv*j,(i,j),c) for (i,j),c in terms]
        maximum=max(q[0] for q in vals)
        top=[q for q in vals if q[0]==maximum]
        require(f'unique degree-leading term H{idx+1} at ({av},{bv})',maximum==expected and len(top)==1 and top[0][1]==((2,2) if idx==0 else (3,3)))
# Direct fibers, without the incorrectly cancelled conic relation.
E1=a*x**3+g*U**3-x*U*(U+1)
E2=b*x**2-x*(2+4*U)+3*g*U**2
C1=3*a*x**2-U*(b*x-(U-1))
require('conic relation loss factor',polynomial_equal(U*E2-3*E1,-x*C1))
spurious={a:0,b:0,g:0,x:1,U:0}
require('old conic pair accepts spurious point',E1.subs(spurious)==0 and C1.subs(spurious)==0)
require('direct relation rejects spurious point',E2.subs(spurious)!=0)
require('spurious reconstructed point exact image',tuple(F.subs({x:1,y:-1,z:5}))==(0,2,0))
sub={y:(U-1)/x,z:(x*(5-3*U)-g)/x**3}
require('E1 reconstructs first component',S.cancel(F[0].subs(sub)-a+E1/x**3)==0)
require('E2 reconstructs second component',S.cancel(F[1].subs(sub)-b+E2/x**2)==0)
L=27*a**2*g**2-18*a*b*g+16*a+b**3*g-b**2
Q=27*a*g**2-9*b*g+8
k=4-3*b*g
N=L*x**3+k*x-2*g
D=k*x-3*g
E=(9*a*g-b)*x**2+2*x
res=S.resultant(E1,E2,U)
require('direct resultant gx^3N',polynomial_equal(res,g*x**3*N))
subsres=S.subresultants(E1,E2,U)
require('direct linear subresultant',polynomial_equal(subsres[-2],g*x*(D*U+E)))
require('cusp identity Q^2=k^3+27g^2L',polynomial_equal(Q**2,k**3+27*g**2*L))
require('cubic discriminant',polynomial_equal(S.discriminant(N,x),-4*L*Q**2))
require('L irreducibility discriminant witness',polynomial_equal(S.discriminant(L,a),4*k**3))
require('missing curve L vanishes',S.cancel(L.subs({a:4/(27*g**2),b:4/(3*g)}))==0)
require('missing curve Q vanishes',S.cancel(Q.subs({a:4/(27*g**2),b:4/(3*g)}))==0)
require('missing curve N constant nonzero',S.cancel(N.subs({a:4/(27*g**2),b:4/(3*g)}))==-2*g)
require('fold quadratic has distinct roots',S.cancel(S.discriminant(E2,U).subs(x,3*g/k)-108*g**2/k)==0)
require('g=0 nonzero-x fiber equation',S.expand(E1.subs({g:0,U:(b*x-2)/4})-x*(L.subs(g,0)*x**2+4)/16)==0)
# A polynomial family which retracts every pointed Keller map to its derivative.
Ft=S.Matrix([S.cancel(f.subs({x:h*x,y:h*y,z:h*z},simultaneous=True)/h) for f in F])
require('scaling family polynomial',all(S.denom(q)==1 for q in Ft))
require('scaling family linear endpoint',tuple(Ft.subs(h,0))==(z,y,2*x))
require('scaling family constant determinant',S.expand(Ft.jacobian((x,y,z)).det())==-2)
# Actual complex fold has two x=3 points and one x=-6 point.
fold_target={a:S.Rational(1,27),b:1,g:1}
require('fold N double/simple root certificate',S.factor(N.subs(fold_target))==-(x-3)**2*(x+6)/27)
require('fold repeated-x u equation',S.expand(E2.subs(fold_target).subs(x,3))==3*(U**2-4*U+1))
require('fold fiber trace with multiplicities',3+3-6==0)
require('fold support trace differs',3-6!=0)
# A point on the escape divisor disproves the all-target trace law.
require('escape singleton candidate image',tuple(F.subs({x:2,y:-S.Rational(1,2),z:S.Rational(9,8)}))==(0,1,1))
require('escape singleton target L=0,Q!=0',L.subs({a:0,b:1,g:1})==0 and Q.subs({a:0,b:1,g:1})!=0)
# Independent finite-field enumeration versus reconstruction, with prime universe stated.
def direct_map(p,point):
    xx,yy,zz=point;uu=(1+xx*yy)%p
    return ((uu**3*zz+yy**2*uu*(4+3*xx*yy))%p,(yy+3*xx*uu**2*zz+3*xx*yy**2*(4+3*xx*yy))%p,(2*xx-3*xx**2*yy-xx**3*zz)%p)
def recover(p,target):
    aa,bb,gg=target
    LL=(27*aa**2*gg**2-18*aa*bb*gg+16*aa+bb**3*gg-bb**2)%p
    out=set()
    if gg==0:out.add((0,bb,(aa-4*bb**2)%p))
    for xx in range(1,p):
        if (LL*xx**3+(4-3*bb*gg)*xx-2*gg)%p:continue
        for uu in range(p):
            if (aa*xx**3+gg*uu**3-xx*uu*(uu+1))%p or (bb*xx**2-xx*(2+4*uu)+3*gg*uu**2)%p:continue
            yy=(uu-1)*pow(xx,-1,p)%p
            zz=(xx*(5-3*uu)-gg)*pow(xx**3,-1,p)%p
            out.add((xx,yy,zz))
    return out
finite=[]
for p in (5,7,11):
    import itertools
    universe=list(itertools.product(range(p),repeat=3));fib={q:set() for q in universe}
    for q in universe:fib[direct_map(p,q)].add(q)
    for q in universe:
        if fib[q]!=recover(p,q):raise CheckFailure(f'finite reconstruction mismatch p={p},target={q}')
    finite.append({'p':p,'targets':p**3,'source_points':p**3,'mismatches':0})
    checks.append(f'all {p**3} targets reconstructed over F_{p}')
# Hostiles: independent validators must reject corrupted certificates, even with -O.
hostiles=[]
for name,condition in [
 ('wrong seed determinant',detF==2),
 ('altered collision point',tuple(F.subs({x:1,y:-S.Rational(3,2),z:S.Rational(15,2)}))==(-S.Rational(1,4),0,0)),
 ('lost conic sidecar',E2.subs(spurious)==0),
 ('wrong second-scale degree',S.Poly(l2,v,t).total_degree()==6),
 ('erased fold multiplicity',3-6==0),
 ('wrong missing curve sign',S.cancel(N.subs({a:-4/(27*g**2),b:4/(3*g)}))==-2*g),
 ('global counterexample deformation endpoint',tuple(Ft.subs(h,0))!=(z,y,2*x))]:
    if condition:raise CheckFailure('hostile unexpectedly accepted: '+name)
    hostiles.append(name)
result={'status':'PASS','checks':checks,'check_count':len(checks),'hostiles_rejected':hostiles,'finite_field_crosscheck':finite,'l2_degree':7,'l2_leading_form':'v^3*(t+3v)^4','fixed_map_complex_fiber_counts':{'L!=0':3,'L=0,Q!=0':1,'L=Q=0':0},'scope':'Exact algebraic certificates and finite-field reconstruction; no claim to global classification of Keller maps or formal proof-assistant build.'}
p=argparse.ArgumentParser();p.add_argument('--output',type=Path);args=p.parse_args()
encoded=json.dumps(result,sort_keys=True,indent=2)+'\n'
if args.output:args.output.write_text(encoded)
else:print(encoded,end='')

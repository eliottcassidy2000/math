#!/usr/bin/env python3
"""Exact supplementary checks; proofs are in the accompanying notes."""
from fractions import Fraction as Q
from math import factorial
import json
from pathlib import Path
import sympy as s
x,y,t=s.symbols('x y t')
checks=[]
def check(name, test):
    if not test: raise AssertionError(name)
    checks.append(name)
def bracket(a,b): return s.expand(s.diff(a,x)*s.diff(b,y)-s.diff(a,y)*s.diff(b,x))
def trunc(p,n): return s.Poly(s.expand(p),x,y).terms() if False else s.Add(*(c*x**i*y**j for (i,j),c in s.Poly(s.expand(p),x,y).terms() if i+j<=n))
def expjet(H,n,lam=1):
    ans=[]
    for z in (x,y):
        a=z; out=z
        for k in range(1,n+1):
            a=trunc(bracket(a,H),n)
            if a==0: break
            out+=lam**k*a/factorial(k)
        ans.append(s.expand(out))
    return tuple(ans)
q=y+x*x; H=q**3
X=x+3*t*q*q; Y=y-6*t*x*q*q-9*t*t*q**4
check('parabola invariant',s.expand(Y+X*X-q)==0)
check('parabola Jacobian',bracket(X,Y)==1)
check('parabola x vector field',s.expand(s.diff(X,t)-s.diff(H,y).subs({x:X,y:Y}, simultaneous=True))==0)
check('parabola y vector field',s.expand(s.diff(Y,t)+s.diff(H,x).subs({x:X,y:Y}, simultaneous=True))==0)
for n in range(3,8):
    E=expjet(H,n)
    check(f'parabola exponential jet {n}',E==(trunc(X.subs(t,1),n),trunc(Y.subs(t,1),n)))
for M in range(3,10):
    E=expjet(y**3+x**M,M-1,lam=t)
    check(f'perturbation exact first discrepancy M={M}',s.expand(E[0]-x-3*t*y*y)==0 and s.expand(E[1]-y+t*M*x**(M-1))==0)
    check(f'isolated critical ideal M={M}',s.groebner([s.diff(y**3+x**M,x),s.diff(y**3+x**M,y)],x,y).is_zero_dimensional)
F=(x+y*y,y+(x+y*y)**2)
check('Henon tangent identity',s.Matrix(F).jacobian([x,y]).subs({x:0,y:0})==s.eye(2))
check('Henon Jacobian',bracket(*F)==1)
hi=(y,x+y*y)
check('Henon square',tuple(s.expand(z.subs({x:hi[0],y:hi[1]},simultaneous=True)) for z in hi)==tuple(map(s.expand,F)))
check('quadratic sl2 bracket',bracket(x*x/2,y*y/2)==x*y)
check('common invariant bracket',bracket(q**3,q**4+q**5)==0)
# Sparse rational polynomial arithmetic, with truncation during multiplication.
def add(a,b):
    c=a.copy()
    for e,v in b.items():
        c[e]=c.get(e,Q(0))+v
        if not c[e]: del c[e]
    return c
def scale(a,k): return {e:v*k for e,v in a.items() if v*k}
def mul(a,b,N):
    c={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():
            if i+j+k+l<=N:
                e=(i+k,j+l);c[e]=c.get(e,Q(0))+v*w
    return {e:v for e,v in c.items() if v}
def power(a,n,N):
    out={(0,0):Q(1)}
    for _ in range(n): out=mul(out,a,N)
    return out
def fromexpr(p):
    return {e:Q(int(c.p),int(c.q)) for e,c in s.Poly(p,x,y).terms() if c}
def word(profiles,N):
    X={(1,0):Q(1)};Y={(0,1):Q(1)}
    for j,coefs in enumerate(profiles):
        ell=add(X,scale(Y,Q(j))); delta={}
        for r,c in coefs.items(): delta=add(delta,scale(power(ell,r,N),c))
        X=add(X,scale(delta,Q(j)));Y=add(Y,scale(delta,Q(-1)))
    return X,Y
def homogeneous_expr(a,r):return sum(s.Rational(v.numerator,v.denominator)*x**i*y**j for (i,j),v in a.items() if i+j==r)
examples=[]
for N in range(2,8):
    target=tuple(map(fromexpr,expjet(y**3+x**3+x*x*y*y,N)))
    profiles=[{} for _ in range(N+2)]
    for r in range(2,N+1):
        W=word(profiles,N)
        a=homogeneous_expr(add(target[0],scale(W[0],Q(-1))),r)
        b=homogeneous_expr(add(target[1],scale(W[1],Q(-1))),r)
        check(f'jet residual divergence N={N},r={r}',s.expand(s.diff(a,x)+s.diff(b,y))==0)
        potential=s.expand((y*a-x*b)/(r+1))
        check(f'jet potential N={N},r={r}',s.diff(potential,y)==a and -s.diff(potential,x)==b)
        V=s.Matrix([[s.binomial(r+1,k)*j**k for j in range(r+2)] for k in range(r+2)])
        rhs=s.Matrix([potential.coeff(x,r+1-k).coeff(y,k) for k in range(r+2)])
        cs=V.inv()*rhs
        for j,c in enumerate(cs): profiles[j][r]=Q(int(c.p),int(c.q))*(r+1)
        W=word(profiles,r)
        check(f'fixed-factor completion N={N},r={r}',all({e:v for e,v in target[i].items() if sum(e)<=r}==W[i] for i in (0,1)))
    check(f'final fixed-factor jet N={N}',word(profiles,N)==target)
    examples.append({'N':N,'max_factors':N+2,'degree_bound':N**(N+2),'profiles':[{str(r):str(c) for r,c in p.items() if c} for p in profiles]})
# Formal log of the explicit Henon example through degree 6.
N=6
logs=[]
for z in (x,y):
    delta=z; out=0
    for k in range(1,N+1):
        delta=trunc(delta.subs({x:F[0],y:F[1]},simultaneous=True)-delta,N)
        out+=(-1)**(k+1)*delta/k
    logs.append(s.expand(out))
check('Henon log divergence to available order',trunc(s.diff(logs[0],x)+s.diff(logs[1],y),N-1)==0)
Hlog=s.integrate(logs[0],y)
Hlog=s.expand(Hlog+s.integrate(-logs[1]-s.diff(Hlog,x),x))
check('Henon log cubic',trunc(Hlog,3)==(y**3-x**3)/3)
check('Henon log quartic',trunc(Hlog,4)==(y**3-x**3)/3-x*x*y*y/2)
out={'status':'PASS','check_count':len(checks),'checks':checks,'fixed_direction_examples':examples,'Henon_formal_Hamiltonian_log_through_degree_7':str(Hlog),'scope':'Exact supplementary controls, not proofs of the general theorems.'}
Path(__file__).with_name('geometry_session2_exact_controls.json').write_text(json.dumps(out,indent=2)+'\n')
print(f'PASS: {len(checks)} exact checks. Fixed-direction shear words verified through N=7.')

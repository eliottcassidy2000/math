#!/usr/bin/env python3
"""Finite exact controls of the boundary-tangent Hamiltonian family."""
import json
import sympy as sp
s,z,x,t,p,y,L=sp.symbols('s z x t p y L')
h=s*s+z
checks=0

def check(ok,label):
    global checks
    if not bool(ok): raise RuntimeError(label)
    checks+=1

def zero(f,label): check(sp.cancel(sp.expand(f))==0,label)
def clean(f): return sp.cancel(sp.expand(f))
def ev(f): return sp.cancel(f.subs(z,-s*s))
def sigma(f): return sp.cancel(f.subs({s:x*t,z:t},simultaneous=True))
def delta(H,f): return clean(z*(sp.diff(f,s)*sp.diff(H,z)-sp.diff(f,z)*sp.diff(H,s)))
def boundary(v,f): return clean(v*sp.diff(f,s))
def partial(op,f,time=1,n=3):
    out=f;term=f
    for j in range(1,n+1):
        term=op(term);out+=sp.Rational(time)**j*term/sp.factorial(j)
    return clean(out)

def val(f,var):
    f=sp.expand(f)
    if f==0: return 10**9
    return min(m[0] for m,c in sp.Poly(f,var).terms() if c!=0)

As=[sp.Integer(0),sp.Integer(1),s,z,h,h*h,s+z,1+s+s*z,h*(1+s+z)]
records=[]
for i,A in enumerate(As):
    H=sp.expand(z*h*A)
    ds=delta(H,s);dz=delta(H,z)
    zero(ds-z*(h*A+z*A+z*h*sp.diff(A,z)),f'ds formula {i}')
    zero(dz+z*z*(2*s*A+h*sp.diff(A,s)),f'dtau formula {i}')
    check(val(ds,z)>=1 and val(dz,z)>=2,f'common convergence {i}')
    dh=delta(H,h)
    zero(dh-z*h*(2*s*A+2*s*z*sp.diff(A,z)-z*sp.diff(A,s)),f'boundary ideal formula {i}')
    v=s**4*ev(A)
    zero(ev(ds)-v,f'boundary ds {i}')
    zero(ev(dz)+2*s*v,f'boundary dtau {i}')
    check(v==0 or val(v,s)>=4,f'boundary convergence {i}')
    dxsource=sigma(delta(H,s/z));dtsource=sigma(dz)
    check(sp.denom(dxsource)==1 and sp.denom(dtsource)==1,f'source polynomial {i}')
    check(val(dxsource,t)>=1 and val(dtsource,t)>=2,f'source convergence {i}')
    for a in range(4):
        for b in range(-2,4):
            f=s**a*z**b
            zero(ev(delta(H,f))-boundary(v,ev(f)),f'intertwining {i}:{a},{b}')
    remainder=sp.rem(A,h,z)
    zero(remainder-ev(A),f'pointwise divisibility {i}')
    check((v==0)==(remainder==0),f'pointwise kernel {i}')
    records.append({'A':str(A),'boundary_velocity':str(v),'pointwise_fixed':v==0})

# Finite partial operators check exponential intertwining and word preservation.
for wi,word in enumerate([[(1,1),(2,1)],[(2,1),(1,1)],[(4,2),(7,-1)]]):
    for initial in [s,z,h,s/z]:
        common=initial;target=ev(initial)
        for index,time in word:
            A=As[index];H=z*h*A;v=s**4*ev(A)
            common=partial(lambda f:delta(H,f),common,time,n=2)
            target=partial(lambda f:boundary(v,f),target,time,n=2)
        zero(ev(common)-target,f'word intertwining {wi}:{initial}')
        if initial==h: zero(ev(common),f'word boundary set {wi}')

# Distinct noncommuting invariants occur in the enlarged family.
H1=z*h;H2=z*h*s
comm=delta(H1,delta(H2,s))-delta(H2,delta(H1,s))
check(sp.expand(comm)!=0,'noncommuting pair')
zero(ev(comm)-(boundary(s**4,boundary(s**5,s))-boundary(s**5,boundary(s**4,s))),
     'boundary commutator transport')

# Exact embedding of the older universal carrier family: A=h^3 R(h,sh).
D=p**3-y**2
for R in [sp.Integer(1),p,y,1+p+y,D]:
    old=(p*p*D*R).subs({p:h,y:s*h},simultaneous=True)
    A=h**3*R.subs({p:h,y:s*h},simultaneous=True)
    zero(old-z*h*A,f'universal inclusion {R}')
    zero(ev(A),f'universal boundary fixedness {R}')
check(val(sigma(z*h),t)==2,'strict enlargement source order')

# A=1 gives nontrivial but setwise boundary motion.
v=s**4
series=s;term=s
for j in range(1,5):
    term=boundary(v,term);series+=L**j*term/sp.factorial(j)
closed=sp.series(s*(1-3*L*s**3)**sp.Rational(-1,3),s,0,14).removeO()
zero(series-closed,'explicit boundary algebraic flow')
check(ev(delta(z*h,s))!=0,'setwise is not pointwise')

# Dropping h permits an exact polynomial shear H=t^2.
H=z*z
zero(delta(H,s/z)-2*z,'shear first derivative')
zero(delta(H,delta(H,s/z)),'shear termination')
zero(delta(H,z),'shear fixes t')
q=1+x*x*t
shearq=q.subs(x,x+2*L*t)-q
check(sp.expand(shearq)!=0,'shear does not preserve q identically')
check(sp.expand(shearq.subs({x:-1/s,t:-s*s},simultaneous=True))!=0,
      'shear fails actual boundary')
for f in [t,t*t,1+t,t+t**3]:
    g=f.subs(t,-s*s)
    check(g!=0 and sp.expand(g*(2-s*g))!=0,f'shear restriction gate {f}')

print(json.dumps({'status':'PASS','checks':checks,'A_controls':records,
 'finite_words':3,'partial_exponential_order':2,
 'boundary_input_monomials_per_A':24,
 'claims':'Finite controls; universal theorem uses the written proof and cited fixed-line theorem.'},indent=2,sort_keys=True))

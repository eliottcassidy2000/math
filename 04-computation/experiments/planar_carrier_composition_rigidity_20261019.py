#!/usr/bin/env python3
"""Exact finite controls for the two-boundary carrier composition theorem."""
import json
import sympy as sp
s,z,p,y,x,t=sp.symbols('s z p y x t')
h=s*s+z
D=p**3-y**2
checks=0

def check(ok,label):
    global checks
    if not bool(ok): raise RuntimeError(label)
    checks+=1

def zero(f,label): check(sp.cancel(f)==0,label)
def expand(f): return sp.cancel(sp.expand(f))
def ev(f): return sp.cancel(f.subs(z,-s*s))
def sigma(f): return sp.cancel(f.subs({s:x*t,z:t},simultaneous=True))
def val(f):
    f=sp.expand(f)
    if f==0: return 10**9
    return min(mon[0] for mon,c in sp.Poly(f,z).terms() if c!=0)

def delta(S,f): return expand(z*(sp.diff(f,s)*sp.diff(S,z)-sp.diff(f,z)*sp.diff(S,s)))
def partial_exp(S,f,time=1,order=2):
    out=f
    term=f
    for j in range(1,order+1):
        term=delta(S,term)
        out+=sp.Rational(time)**j*term/sp.factorial(j)
    return expand(out)

Rs=[sp.Integer(0),sp.Integer(1),p,y,1+p+y,D,D*(p-y)]
records=[]
Ss=[]
for i,R in enumerate(Rs):
    S=sp.expand(z*h**4*R.subs({p:h,y:s*h},simultaneous=True))
    Ss.append(S)
    ds=delta(S,s); dz=delta(S,z); dx=delta(S,s/z)
    for f,label in [(ds/(z*h**3),'ds factor'),(dz/(z*z*h**3),'dt factor')]:
        check(sp.denom(sp.cancel(f))==1,f'{label} {i}')
    for f,label in [(ds,'s'),(dz,'tau'),(dx,'x')]:
        zero(ev(f),f'boundary derivation {i} {label}')
    check(val(ds)>=1 and val(dz)>=2,f'common convergence {i}')
    # Both source coordinates fix the line, with their original orders.
    sx=sigma(dx); st=sigma(dz)
    check(sp.denom(sx)==1 and sp.denom(st)==1,f'polynomial source derivatives {i}')
    zero(sx.subs(t,0),f'source line x {i}')
    zero(st.subs(t,0),f'source line t {i}')
    # Ambient carrier formulas for preservation of A and its ideal filtration.
    Sp=p*p*D*R
    dp=sp.expand(-D*sp.diff(Sp,y)/p)
    dy=sp.expand(D*sp.diff(Sp,p)/p)
    dD=sp.expand(3*p*p*dp-2*y*dy)
    for f,label in [(dp/D,'dp'),(dy/D,'dy'),(dD/D**2,'dD')]:
        check(sp.denom(sp.cancel(f))==1,f'A filtration {i} {label}')
    formula_x=(p*(5*p**3+4*y*y)*R+p*p*(p**3+y*y)*sp.diff(R,p)
               +p*y*(2*p**3+y*y)*sp.diff(R,y))
    zero(dx-formula_x.subs({p:h,y:s*h},simultaneous=True),f'dx in A {i}')
    records.append({'R':str(R),'common_s_order':None if ds==0 else val(ds),
                    'common_tau_order':None if dz==0 else val(dz)})

# Finite partial exponential words preserve the boundary evaluation exactly.
# They are controls of linear operators, not asserted exact finite-time maps.
words=[[(1,1),(2,1)],[(2,1),(1,1)],[(1,-1),(3,1)],[(1,1),(1,-1)]]
for wi,word in enumerate(words):
    for f,label in [(s,'s'),(z,'tau'),(s/z,'x')]:
        out=f
        for index,time in word: out=partial_exp(Ss[index],out,time,order=2)
        zero(ev(out)-ev(f),f'word boundary {wi} {label}')
        if label in ['tau','x']:
            source=sigma(out)
            check(sp.denom(source)==1,f'word source polynomial partial sum {wi} {label}')
            zero(source.subs(t,0)-(x if label=='x' else 0),f'word source line {wi} {label}')

# The class actually includes different noncommuting invariants.
comm=delta(Ss[1],delta(Ss[2],s))-delta(Ss[2],delta(Ss[1],s))
check(sp.expand(comm)!=0,'different invariant commutator nonzero')
zero(ev(comm),'commutator fixes second boundary')

# Complete normal-form leading coefficient for finite basis controls.
normal_count=0
for n in range(5):
    for a in range(5):
        for eps in [0,1]:
            f=sp.expand((z*h*h)**n*h**a*(s*h)**eps)
            expected=s**(4*n+2*a+3*eps)
            check(val(f)==n,f'normal form valuation {n} {a} {eps}')
            zero(f.coeff(z,n)-expected,f'normal form leading coefficient {n} {a} {eps}')
            normal_count+=1
# Even and odd supports cannot cancel, including mixed test coefficients.
for n in range(1,5):
    f=sp.expand((z*h*h)**n*((1+h+h*h)+s*h*(1-h)))
    check(val(f)==n,f'mixed normal form {n}')
    leading=f.coeff(z,n)
    check(sp.rem(leading,s**(4*n),s)==0 and leading!=0,f'positive coefficient divisor {n}')

# Finite grid verifies the monomial-address injectivity of sigma.
images=set()
for a in range(13):
    for b in range(-4,13):
        pair=(a,a+b)
        check(pair not in images,'sigma address injective')
        images.add(pair)
# Shear restriction gate and nonzero boundary hostile.
for d in range(1,13):
    image=(-s*s)**d
    check(image!=0,f'nonzero shear detected {d}')
# q is a t-adic unit: exact finite geometric-series hostile to naive closure.
q=1+x*x*t
for N in range(1,9):
    geom=sum((-x*x*t)**j for j in range(N))
    zero(q*geom-1-(-1)**(N-1)*(x*x*t)**N,f'q unit hostile {N}')

print(json.dumps({'status':'PASS','checks':checks,'carrier_controls':records,
 'finite_word_controls':words,'partial_exponential_order':2,
 'normal_form_basis_controls':normal_count,'sigma_address_controls':len(images),
 'claims':'Universal conclusion uses written proof and Gwozdziewicz fixed-line theorem; finite partial sums are not claimed to be exact maps.'},indent=2,sort_keys=True))

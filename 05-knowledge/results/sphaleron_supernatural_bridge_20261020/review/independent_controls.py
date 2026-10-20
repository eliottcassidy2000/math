#!/usr/bin/env python3
"""Independent finite referee; no imports from the source verifier."""
import json
import math
from collections import Counter
from itertools import product
import sympy as sp

checks=0

def need(ok,label):
    global checks
    checks+=1
    if not ok: raise RuntimeError(label)

x=sp.Symbol('x')

def polynomial(bits):
    return sp.Poly(sum(((bits>>i)&1)*x**i for i in range(bits.bit_length())),x,modulus=2)

def bits(poly):
    return sum((int(c)%2)<<i for i,c in enumerate(reversed(poly.all_coeffs())))

def multiply_by_x(a,f):
    n=f.bit_length()-1
    answer=a<<1
    if answer&(1<<n): answer^=f
    return answer

def mul(a,b,f):
    # Independent quotient-ring linear-operator algorithm.
    out=0
    while b:
        if b&1: out^=a
        a=multiply_by_x(a,f); b>>=1
    return out

def eval_poly(p,a,f):
    out=0
    for i in range(p.bit_length()-1,-1,-1):
        out=mul(out,a,f)^((p>>i)&1)
    return out

def element_x_order(f):
    n=f.bit_length()-1
    current=1
    for t in range(1,1<<n):
        current=multiply_by_x(current,f)
        if current==1: return t
    raise RuntimeError('Unit order exceeds quotient ring size')

# Every constant-one monic polynomial, degrees 1..9, including inseparable
# repeated factors. Compare factor formula with literal full basis evolution.
registers=[]
for n in range(1,10):
    count=0; maxorder=0; repeated=0
    for f in range((1<<n)|1,1<<(n+1),2):
        factors=sp.factor_list(polynomial(f))[1]
        prediction=1
        for p,e in factors:
            r=element_x_order(bits(p))
            need(r%2==1,'irreducible factor root order odd')
            thick=1
            while thick<e: thick*=2
            prediction=math.lcm(prediction,r*thick)
        observed=element_x_order(f)
        need(prediction==observed,'repeated-factor operator formula')
        image=[1<<i for i in range(n)]
        for t in range(observed): image=[multiply_by_x(a,f) for a in image]
        need(image==[1<<i for i in range(n)],'whole basis identity')
        repeated+=int(any(e>1 for _,e in factors)); count+=1; maxorder=max(maxorder,observed)
    registers.append({'degree':n,'polynomials':count,'repeated_factor_polynomials':repeated,'maximum_order':maxorder})

# Normal-basis shift has x^n-1 as minimal polynomial: 1,e_1,... form a
# literal cyclic orbit. Fixed vectors also independently count necklaces.
shift_rows=[]
for n in range(1,13):
    mask=(1<<n)-1
    shift=lambda a: ((a<<1)&mask)|(a>>(n-1))
    need(len({1<<i for i in range(n)})==n,'distinct basis orbit')
    curr=1
    for k in range(1,n+1):
        curr=shift(curr)
        need((curr==1)==(k==n),'cyclic shift exact order')
        fixed=0
        for v in range(1<<n):
            a=v
            for _ in range(k): a=shift(a)
            fixed+=a==v
        need(fixed==2**math.gcd(n,k),'cyclic shift fixed count')
    shift_rows.append({'n':n,'operator_order':n,'fixed_states':2})

# Scheme/cover arrow directions: a map of n-covers to m-covers corresponds
# to an embedding F_(2^m) -> F_(2^n). Root counts supply independent check.
field_moduli={}
for n in range(1,7):
    field_moduli[n]=next(f for f in range((1<<n)|1,1<<(n+1),2) if polynomial(f).is_irreducible)
cover_rows=[]
for n,m in product(range(1,7),repeat=2):
    equivariant=0
    for image0 in range(m):
        f=[(image0+i)%m for i in range(n)]
        equivariant+=all(f[(i+1)%n]==(f[i]+1)%m for i in range(n))
    roots=sum(eval_poly(field_moduli[m],a,field_moduli[n])==0 for a in range(1<<n))
    prediction=m if n%m==0 else 0
    need(equivariant==roots==prediction,'cover/finite-field hom counts')
    cover_rows.append({'source_degree':n,'target_degree':m,'morphisms':prediction})

# Independent degree route: the regular target (0,-1,0,0) has exactly one
# preimage: mu=theta=pi/2, phi=0 modulo 2pi, since all three factors of
# sin(mu)*sin(theta)*cos(phi)=1 must equal one. The local degree is +1.
mu,th,ph=sp.symbols('mu th ph',real=True)
y=sp.Matrix([-sp.sin(mu)*sp.sin(th)*sp.sin(ph),-sp.sin(mu)*sp.sin(th)*sp.cos(ph),sp.sin(mu)*sp.cos(mu)*(sp.cos(th)-1),sp.cos(mu)**2+sp.sin(mu)**2*sp.cos(th)])
point={mu:sp.pi/2,th:sp.pi/2,ph:0}
need(y.subs(point)==sp.Matrix([0,-1,0,0]),'regular target')
local=sp.Matrix.hstack(y,y.diff(mu),y.diff(th),y.diff(ph)).subs(point).det()
need(local==1,'positive local degree at unique preimage')

midpoints=[]
for n in range(1,65):
    roots=[a for a in range(n) if (2*a-1)%n==0]
    need(len(roots)==(1 if n%2 else 0),'midpoint parity obstruction')
    midpoints.append({'modulus':n,'solutions':len(roots)})

print(json.dumps({'status':'INDEPENDENT FINITE-EXACT PASS','checks':checks,'all_constant_one_registers_degrees_1_to_9':registers,'cyclic_shifts':shift_rows,'cover_and_field_morphism_pairs':cover_rows,'degree_check':{'method':'unique regular preimage and positive local determinant','local_determinant':int(local),'degree':1},'midpoint_moduli_tested':64,'scope':'No physical dynamics, continuous path identification, or infinite theorem inferred from finite data.'},indent=2,sort_keys=True))

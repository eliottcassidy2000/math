#!/usr/bin/env python3
"""Bounded exact controls for BRIDGE.md; universal claims have written proofs."""
import json
from math import gcd, comb
from pathlib import Path
import sympy as sp

COUNTS = {}
def gate(condition, label):
    if not condition:
        raise RuntimeError(label)
    COUNTS[label.split(':')[0]] = COUNTS.get(label.split(':')[0], 0) + 1

def norm(v,p):
    for x in v:
        if x % p:
            inv=pow(x,-1,p)
            return tuple((a*inv)%p for a in v)
    raise RuntimeError('zero projective vector')

def act(v,r,p):
    return norm(tuple(pow(x,r,p) for x in v),p)

primes=(2,3,5,7,11,13,17,19,23,29,31)
for p in primes:
    line={(1,b) for b in range(p)} | {(0,1)}
    conic={(1,v,w) for v in range(p) for w in range(p) if (v-w*w)%p==0}
    conic |= {(0,1,w) for w in range(p) if w*w%p==0}
    def ver(v):
        a,b=v
        return norm((a*a,b*b,a*b),p)
    gate({ver(v) for v in line}==conic,'conic_chart:all_points')
    for r in range(1,13):
        gate(all(ver(act(v,r,p))==act(ver(v),r,p) for v in line),'conic_equivariance:all_points')
        fixed=sum(act(v,r,p)==v for v in conic)
        gate(fixed==2+gcd(p-1,r-1),'fixed_count:prime_field')
        for t in range(4):
            for n in range(1,6):
                raw=sum(act(v,r**(t+n),p)==act(v,r**t,p) for v in conic)
                gate(raw==2+gcd(p-1,r**t*(r**n-1)),'tail_period:prime_field')
    for m in range(2,13):
        for n in range(1,m):
            if gcd(m,n)!=1:
                continue
            a=norm((m,n),p)
            for r in range(1,6):
                reduced=norm((pow(m,r),pow(n,r)),p)
                gate(act(a,r,p)==reduced,'integer_reduction:commutes')
                d=m-n
                dr=d*sum(m**(r-1-j)*n**j for j in range(r))
                gate(dr==m**r-n**r,'exceptional_difference:retained')

# Independent explicit orbit inventory for the small example in the note.
def orbit_type(x,p=7,r=2):
    seen={}; seq=[]
    while x not in seen:
        seen[x]=len(seq);seq.append(x);x=act(x,r,p)
    return seen[x],len(seq)-seen[x]
inv={}
for a in {(1,b) for b in range(7)}|{(0,1)}:
    key=orbit_type(a);inv[key]=inv.get(key,0)+1
gate(inv=={(0,1):3,(0,2):2,(1,1):1,(1,2):2},'example:Q7_square_map')

# Polynomial identities over Z; no generic-fiber-only numerical test.
t,u,v=sp.symbols('t u v')
def mul(r,z):
    return sum(comb(r,j)*t**(j-1)*z**j for j in range(1,r+1))
def zero(expr):
    return sp.Poly(sp.expand(expr),t,u,v).is_zero
for r in range(1,7):
    f=mul(r,u)
    gate(zero(1+t*f-(1+t*u)**r),'power_polynomial:identity')
    gate(zero(f.subs(t,0)-r*u),'additive_specialization:integer_scalar')
    gate(zero(sp.diff(f,u).subs(u,0)-r),'tangent:integer_scalar')
    fu,fv=mul(r,u),mul(r,v)
    gate(zero(mul(r,u+v+t*u*v)-(fu+fv+t*fu*fv)),'group_law:preserved')
    for s in range(1,7):
        gate(zero(mul(r,mul(s,u))-mul(r*s,u)),'composition:semigroup')
for p in (2,3,5,7,11,13):
    poly=sp.Poly(sp.expand(mul(p,u)-t**(p-1)*u**p),t,u,modulus=p)
    gate(poly.is_zero,'prime_power:mod_p_formula')
    gate(sp.Poly(mul(p,u).subs(t,0),u,modulus=p).is_zero,'zero_special_fiber:not_frobenius')
    gate(not sp.Poly(u**p,u,modulus=p).is_zero,'frobenius_special_fiber:nonzero')
    relative=(u+v+t*u*v)**p-(u**p+v**p+t**p*u**p*v**p)
    gate(sp.Poly(sp.expand(relative),t,u,v,modulus=p).is_zero,'relative_frobenius:base_twist')

receipt={
    'status':'PASS',
    'scope':{
        'conic_and_dynamics':'p in [2,3,5,7,11,13,17,19,23,29,31]; r=1..12; t=0..3; period-divisor n=1..5',
        'integer_reduction':'coprime 1<=n<m<=12; powers=1..5; same primes',
        'polynomials':'exact Z[t,u,v], r,s=1..6; characteristic-p checks p=2,3,5,7,11,13',
        'limits':'No extension-field enumeration, root-surface recertification, parent-automaton construction, or external novelty test.'
    },
    'gates':dict(sorted(COUNTS.items())),
    'total_gates':sum(COUNTS.values()),
    'Q7_r2_tail_period_points':{str(k):v for k,v in sorted(inv.items())},
    'sympy_version':sp.__version__,
}
print(json.dumps(receipt,sort_keys=True,indent=2))

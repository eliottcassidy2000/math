#!/usr/bin/env python3
"""Exact finite controls for the proved metabelian word/jet theorem.

The tests supplement the proofs; they do not establish general minimality.
Run with SymPy installed: python3 check_metabelian_jets.py
"""
import hashlib
import json
from pathlib import Path
import sympy as s

x,y,t,c,u,z,w=s.symbols('x y t c u z w')
checks=[]

def check(label, assertion):
    assert assertion, label
    checks.append(label)

def eq(a,b):
    return s.cancel(a-b)==0

def M(a):
    return (x/(1-a*x),y*(1-a*x)**2)

def T(b):
    return (x,y+b*x**2)

def compose(after,before):
    return tuple(s.cancel(h.subs({x:before[0],y:before[1]}, simultaneous=True)) for h in after)

def same(a,b):
    return all(eq(p,q) for p,q in zip(a,b))

def r(a):
    return x**2/(1-a*x)**4

check('M additive subgroup',same(compose(M(u),M(t)),M(t+u)))
check('T additive subgroup',same(compose(T(u),T(c)),T(c+u)))
check('conjugated vertical shear',same(compose(M(-t),compose(T(c),M(t))),(x,y+c*r(t))))
check('canonical coordinate bracket',eq(s.diff(1/x,x)*s.diff(-y*x**2,y)-s.diff(1/x,y)*s.diff(-y*x**2,x),1))
check('M Jacobian one',eq(s.det(s.Matrix(M(t)).jacobian([x,y])),1))
check('G=x^2*y vector field',same(tuple(s.diff(v,t).subs(t,0) for v in M(t)),(x**2,-2*x*y)))

rows=[]
for N in range(3,15):
    p=N-2
    if p==1:
        nodes=[s.Integer(1)]
    else:
        harmonic=sum(s.Rational(1,j) for j in range(1,p))
        nodes=list(map(s.Integer,range(1,p)))+[-1/harmonic]
    Q=s.prod(z-a for a in nodes).expand()
    Q0=Q.subs(z,0)
    c0=s.diff(Q,z).subs(z,0)/(4*Q0)
    coeff=[-Q0/(4*a*a*s.diff(Q,z).subs(z,a)) for a in nodes]
    check(f'N={N}: distinct nonzero sites',len(set(nodes))==p and all(a!=0 for a in nodes))
    check(f'N={N}: all residues nonzero',all(b!=0 for b in coeff))
    check(f'N={N}: zero-site coefficient',c0==(-s.Rational(1,4) if N==3 else 0))
    moments=[s.simplify(sum(b*a**k for a,b in zip(nodes,coeff))+(c0 if k==0 else 0)) for k in range(p+2)]
    target=[s.Integer(0),s.Rational(1,4)]+[s.Integer(0)]*(p-1)
    check(f'N={N}: jet moments',moments[:p+1]==target)
    check(f'N={N}: exact next moment',eq(moments[p+1],-Q0/4) and moments[p+1]!=0)
    error=s.binomial(p+4,3)*moments[p+1]
    word=[]
    if c0:
        word.append(('T',c0))
    prev=s.Integer(0)
    for a,b in zip(nodes,coeff):
        word.extend([('M',a-prev),('T',b)])
        prev=a
    word.append(('M',-prev))
    check(f'N={N}: no zero-time factors',all(v!=0 for _,v in word))
    check(f'N={N}: sharp factor count',len(word)==(4 if N==3 else 2*N-3))
    if N<=7:
        f=s.cancel(c0*x*x+sum(b*r(a) for a,b in zip(nodes,coeff)))
        denominator=s.fraction(f)[1]
        check(f'N={N}: reduced denominator degree',s.degree(denominator,x)==4*p)
        check(f'N={N}: rational jet and first error',eq(s.series(f,x,0,N+2).removeO(),x**3+error*x**(N+1)))
        point=(x,y)
        for typ,value in word:
            point=compose(M(value) if typ=='M' else T(value),point)
        check(f'N={N}: chronological word realizes residues',same(point,(x,y+f)))
    rows.append({'N':N,'nonzero_sites':p,'nodes':[str(a) for a in nodes],
        'coefficients':[str(b) for b in coeff],'zero_site_coefficient':str(c0),
        'word':[[typ,str(value)] for typ,value in word], 'word_length':len(word),
        'reduced_denominator_degree':4*p,'first_error_degree':N+1,'first_error_coefficient':str(error)})

# Arbitrary distinct nonzero sites attain the pole minimum, but usually need c0.
for p in range(1,9):
    nodes=list(map(s.Integer,range(1,p+1)))
    Q=s.prod(z-a for a in nodes).expand(); Q0=Q.subs(z,0)
    c0=s.diff(Q,z).subs(z,0)/(4*Q0)
    coeff=[-Q0/(4*a*a*s.diff(Q,z).subs(z,a)) for a in nodes]
    for k in range(p+1):
        moment=s.simplify(sum(b*a**k for a,b in zip(nodes,coeff))+(c0 if k==0 else 0))
        check(f'arbitrary sites p={p}, moment {k}',moment==(s.Rational(1,4) if k==1 else 0))
    check(f'arbitrary sites p={p}: c0 reciprocal formula',eq(c0,-sum(1/a for a in nodes)/4))

# Pure point-adic density: solve independent finite Vandermonde systems.
for N in range(2,11):
    nodes=list(map(s.Integer,range(N-1)))
    coeffs=[s.Rational((-1)**k*(k+1),k+3) for k in range(N-1)]
    V=s.Matrix([[s.binomial(k+3,3)*a**k for a in nodes] for k in range(N-1)])
    weights=V.inv()*s.Matrix(coeffs)
    check(f'density N={N}: arbitrary rational target jet',V*weights==s.Matrix(coeffs))

out={'status':'PASS','check_count':len(checks),'field':'exact rationals (finite controls only)',
     'scope':'N=3..14 sharp cubic jets; exact rational words N<=7; arbitrary-site and density controls',
     'checks':checks,'sharp_cubic_approximants':rows,
     'proof_warning':'All-order minimality and closure are proved in the companion note, not inferred from these finite checks.'}
path=Path(__file__).with_name('metabelian_jets_exact_controls.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(f'PASS: {len(checks)} exact controls. Output: {path.name}')
print('Output SHA256:',hashlib.sha256(path.read_bytes()).hexdigest())

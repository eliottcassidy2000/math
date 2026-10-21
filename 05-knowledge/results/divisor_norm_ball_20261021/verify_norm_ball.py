#!/usr/bin/env python3
"""Exact controls for DIVISOR_NORM_BALL.md; stdlib only; checks survive -O."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd, isqrt
import json

checks = 0
def need(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)

def source(k):
    return [v for v in product(range(k + 1), range(2), range(2))
            if v not in [(0, 0, 0), (k, 1, 1)]]

def encode(k, v):
    a, b, c = v
    return (a - k*c, k*(b-c))  # actual basis (1,omega), not (1,k omega)

def decode(k, z):
    x, t = z
    if t == k:
        return (x, 1, 0)
    if t == -k:
        return (x+k, 0, 1)
    if t == 0 and x > 0:
        return (x, 0, 0)
    if t == 0 and x < 0:
        return (x+k, 1, 1)
    raise ValueError('outside the nonzero norm ball')

def norm(z):
    x, t = z
    return x*x-x*t+t*t

def conjugate(z):
    x, t = z
    return (x-t, -t)

def neg(z):
    return (-z[0], -z[1])

def rotor(z):
    x, t = z
    return (-t, x-t)

def complement(k, v):
    a, b, c = v
    return (k-a, 1-b, 1-c)

def exchange(v):
    a, b, c = v
    return (a, c, b)

def factors(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out

def chi(p):
    return 0 if p == 3 else (1 if p % 3 == 1 else -1)

def class_number(f):
    if f == 1:
        return 1
    h = Fraction(f, 3)
    for p in factors(f):
        h *= Fraction(p-chi(p), p)
    need(h.denominator == 1, 'class number integrality')
    return h.numerator

def forms(D):
    out = []
    for A in range(1, isqrt((-D)//3)+1):
        for B in range(-A, A+1):
            if (B*B-D) % (4*A):
                continue
            C = (B*B-D)//(4*A)
            if A > C or gcd(gcd(A, B), C) != 1:
                continue
            if (abs(B) == A or A == C) and B < 0:
                continue
            out.append((A, B, C))
    return out

def ambiguous(t):
    A, B, C = t
    return B == 0 or abs(B) == A or A == C

def label(v):
    out = []
    for letter, e in zip('pqr', v):
        if e:
            out.append(letter if e == 1 else f'{letter}^{e}')
    return ''.join(out) or '1'

ball_rows = []
for k in range(1, 65):
    src = source(k)
    image = {encode(k, v) for v in src}
    # The completed square bounds |y|<=1 and |x|<=k; deliberately enlarge both.
    ball = {(x,k*y) for x in range(-2*k,2*k+1) for y in range(-3,4)
            if 0 < norm((x,k*y)) <= k*k}
    need(len(src) == len(image) == 4*k+2, f'injective {k}')
    need(image == ball, f'exhaustive norm ball {k}')
    need(encode(k,(0,0,0)) == encode(k,(k,1,1)) == (0,0), 'endpoints')
    for v in src:
        z = encode(k,v)
        a,b,c = v
        need(decode(k,z) == v, 'source decode')
        need(encode(k,complement(k,v)) == neg(z), 'C transport')
        need(encode(k,exchange(v)) == conjugate(z), 'Q transport')
        need(encode(k,complement(k,exchange(v))) == neg(conjugate(z)), 'CQ transport')
        need(2*z[0]-z[1] == 2*a-k*(b+c), 'real chart')
        need(norm(z) == norm(rotor(z)) == norm(conjugate(z)), 'norm invariance')
        need((rotor(z) in ball) == (z[0] % k == 0), 'rotor membership')
    fixed_C = [v for v in src if complement(k,v)==v]
    fixed_Q = [v for v in src if exchange(v)==v]
    fixed_CQ = [v for v in src if complement(k,exchange(v))==v]
    need(not fixed_C and len(fixed_Q)==2*k, 'fixed C,Q')
    need(len(fixed_CQ)==(2 if k%2==0 else 0), 'CQ parity')
    core = [v for v in src if encode(k,v)[0] % k == 0]
    expected_core = {(k,0,0),(0,1,0),(0,0,1),(0,1,1),(k,1,0),(k,0,1)}
    need(set(core)==expected_core, 'literal rotor core')
    need(all(norm(encode(k,v)) == k*k for v in core), 'core outer shell')
    need(all(rotor(encode(k,v)) in {encode(k,u) for u in core} for v in core), 'core closed')
    real_sign = Counter((2*a-k*(b+c)>0)-(2*a-k*(b+c)<0) for a,b,c in src)
    imag_sign = Counter(b-c for a,b,c in src)
    need([imag_sign[i] for i in (-1,0,1)] == [k+1,2*k,k+1], 'imag split')
    need([real_sign[i] for i in (-1,0,1)] == ([2*k,2,2*k] if k%2==0 else [2*k+1,0,2*k+1]), 'real split')
    if k>=2:
        need(sorted(z for z in ball if norm(z)==1)==[(-1,0),(1,0)], 'minimal vectors')
        F,S,U = len(src),sum(max(v)<=1 for v in src),3
        need(F-S-U==4*(k-2), 'balance defect')
    if k <= 8:
        ball_rows.append({'k':k,'size':len(ball),'Q_fixed':len(fixed_Q),'CQ_fixed':len(fixed_CQ),'rotor_core':len(core)})

# Direct integer divisors do not use the exponent-box enumerator.
for p,q,r in [(2,3,5),(3,2,7),(5,2,3)]:
    for k in range(1,6):
        N=p**k*q*r
        ds=set()
        for d in range(1,isqrt(N)+1):
            if N%d==0:
                ds.update([d,N//d])
        ds.discard(1); ds.discard(N)
        encoded={p**a*q**b*r**c for a,b,c in source(k)}
        need(ds==encoded, 'direct integer divisor set')

# Entire-lattice rotor/conductor identities, not restricted to the norm ball.
for k in range(1,17):
    for x,y in product(range(-3*k,3*k+1),range(-3,4)):
        z=(x,k*y)
        need((rotor(z)[1]%k==0)==(x%k==0), 'whole-lattice rotation test')
        need((x%k==0 and z[1]%k==0)==(z[1]%k==0 and rotor(z)[1]%k==0), 'conductor identity')

literal=[]
for v in sorted(source(2),key=lambda v:(sum(v),v),reverse=False):
    z=encode(2,v)
    literal.append({'divisor':label(v),'exponents':list(v),'eisenstein_basis':list(z),
                    'rectangular_basis':[z[0]-z[1]//2,z[1]//2], 'norm':norm(z),
                    'C':label(complement(2,v)),'Q':label(exchange(v)),
                    'CQ':label(complement(2,exchange(v)))})
need(Counter(row['norm'] for row in literal)=={1:2,3:2,4:6}, 'ten-point norm census')
native={(1,0),(2,1),(1,1),(1,2),(0,1),(-1,1),(-1,0),(-2,-1),(-1,-1),(-1,-2),(0,-1),(1,-1)}
need(Counter(map(norm,native))=={1:6,3:6}, 'inherited native alphabet')
need(all(rotor(z) in native for z in native), 'native rotor closed')
need(rotor((1,0)) not in {encode(2,v) for v in source(2)}, 'smallest rotor hostile')
need(any(norm(encode(2,v)) not in set(map(norm,native)) for v in source(2)), 'norm-preserving target hostile')

small_expected={
 2:[(1,0,3)],
 4:[(1,0,12),(3,0,4)],
 8:[(1,0,48),(3,0,16),(4,4,13),(7,2,7)],
 16:[(1,0,192),(3,0,64),(4,4,49),(7,-4,28),(7,4,28),(12,12,19),(13,-8,16),(13,8,16)]}
class_rows=[]
for f in range(1,129):
    fs=forms(-3*f*f)
    need(len(fs)==class_number(f), f'form count versus class formula {f}')
    if f%2==0:
        need(all(ambiguous(t) for t in fs)==(f in (2,4,8)), f'even classification {f}')
        fac=factors(f)
        value=2**(fac.get(2,0)-1)*3**fac.get(3,0)
        for p,e in fac.items():
            if p>=5:
                value*=p**(e-1)*(p-chi(p))
                need((p-chi(p))%3==0, 'odd-prime obstruction')
        need(value==len(fs), 'even class number product')
    if f in small_expected:
        need(fs==small_expected[f], 'complete literal reduced form list')
        class_rows.append({'conductor':f,'discriminant':-3*f*f,'idoneal_integer':3*(f//2)**2,
                           'class_number':len(fs),'ambiguous_classes':sum(map(ambiguous,fs)),
                           'forms':[list(t) for t in fs]})
need((2,0,6) not in forms(-48), 'reject imprimitive false class')
need(class_number(16)==8 and not all(map(ambiguous,forms(-768))), '2-power order not exponent-two')
for a in range(4,129):
    f=2**a
    A,B,C=(13,8,16) if a==4 else (16,8,1+3*4**(a-3))
    need(B*B-4*A*C==-3*f*f, 'uniform dyadic witness discriminant')
    need(gcd(gcd(A,B),C)==1, 'uniform dyadic witness primitive')
    need(0<B<A<C, 'strictly reduced nonambiguous dyadic witness')

result={'status':'PASS','live_checks':checks,'scope':{'norm_balls_k':[1,64],
 'whole_lattice_k':[1,16],'direct_divisor_assignments':3,'direct_divisor_k':[1,5],
 'reduced_forms_conductors':[1,128],'dyadic_witness_exponents':[4,128]},
 'ball_rows':ball_rows,'literal_k2_table':literal,'even_idoneal_rows':class_rows,
 'native_alphabet_norm_histogram':{'1':6,'3':6},
 'k2_norm_histogram':{'1':2,'3':2,'4':6},
 'uniform_theorem_basis':'written proofs, not extrapolation from these finite checks'}
print(json.dumps(result,indent=2,sort_keys=True))

#!/usr/bin/env python3
"""Independent exact finite audit; checks stay active under python -O.

No writes outside this audit directory; frozen producer implementation is imported
only for cross-comparison. Reference is direct proper-divisor accumulation;
independent implementation uses trial factorization and canonical squarefull pairs.
"""
import hashlib
import importlib.util
import json
from math import gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('producer', ROOT/'source-snapshot/aliquot_fibers.py')
producer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(producer)
checks = 0

def require(condition, context):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(context)

def factor(n):
    out = {}
    for p in range(2, isqrt(n)+1):
        if p*p > n:
            break
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out

def aliquot(n):
    value = 1
    for p,e in factor(n).items():
        value *= sum(p**k for k in range(e+1))
    return value-n

def prime(n):
    return n >= 2 and factor(n) == {n:1}

def cube_root(n):
    r = 0
    while (r+1)**3 <= n:
        r += 1
    return r

def squarefull_pairs(B):
    out = {}
    for v in range(1, cube_root(B)+1):
        if any(e != 1 for e in factor(v).values()):
            continue
        for u in range(1, isqrt(B//v**3)+1):
            n = u*u*v**3
            require(n not in out, ('nonunique squarefull form',n))
            out[n] = (u,v)
    return out

M = 1000
B = (M-1)**2
reference_s = [0]*(B+1)
for d in range(1,B//2+1):
    for n in range(2*d,B+1,d):
        reference_s[n] += d
reference = {m:[] for m in range(2,M+1)}
bound_sources = 0
for n in range(1,B+1):
    m = reference_s[n]
    if 2 <= m <= M:
        reference[m].append(n)
    if m <= 1:
        continue
    bound_sources += 1
    require(n <= (m-1)**2, ('global bound',n,m))
    r = isqrt(n)
    is_prime_square = r*r == n and reference_s[r] == 1
    require((n == (m-1)**2) == is_prime_square, ('sharpness iff',n,m))
    if not is_prime_square:
        require(4*n <= (m-1)**2, ('non-prime-square bound',n,m))
    if n%2 == 0 and n >= 6:
        require(n <= 2*m-6, ('even bound',n,m))
    twice_square = n%2 == 0 and isqrt(n//2)**2 == n//2
    require(((m+n)%2 == 1) == (r*r == n or twice_square), ('sigma parity',n,m))
    if n%2 and m%2 == 0:
        require(r*r == n, ('odd source even target',n,m))
        if not is_prime_square:
            require(m*m > r**3, ('strict composite root bound',r,m))
for n in range(1,10001):
    require(aliquot(n) == reference_s[n], ('factor oracle mismatch',n))

pairs = squarefull_pairs(B)
require(len(pairs) <= 3*(M-1), 'squarefull count bound')
squarefull_values = {n:aliquot(n) for n in pairs if n > 1}
computed = {}
for m in range(2,M+1):
    found = set()
    for a in range(2,m):
        S = reference_s[a]
        T = a+S
        numerator = m-T
        if numerator <= 0 or numerator%S:
            continue
        p = numerator//S
        if p >= 2 and reference_s[p] == 1 and gcd(a,p) == 1:
            found.add(a*p)
    found.update(n for n,t in squarefull_values.items() if t == m)
    computed[m] = sorted(found)
    require(computed[m] == reference[m], ('independent inversion mismatch',m))
    require(producer.fiber(m) == reference[m], ('producer inversion mismatch',m))

def census(X):
    fibers = {m:[] for m in range(2,X+1,2)}
    for n in sorted({4} | set(range(6,max(5,2*X-6)+1,2))):
        m = aliquot(n)
        if m in fibers:
            fibers[m].append(n)
    for a in range(3,X,2):
        if prime(a):
            fibers[a+1].append(a*a)
        elif a**3 <= X*X:
            m = aliquot(a*a)
            if m in fibers:
                fibers[m].append(a*a)
    return {m:sorted(set(v)) for m,v in fibers.items()}

for X in list(range(0,41))+[99,100,121,500,999,1000]:
    expected = {m:reference[m] for m in range(2,X+1,2)}
    require(census(X) == expected, ('independent prefix',X))
    require(producer.even_prefix(X) == expected, ('producer prefix',X))

# Every possible m is checked over the finite named source range; prime/gcd
# exclusions are applied independently, without using the hyperbola identity.
hyperbola_cases = 0
for a in range(2,51):
    expected_by_m = {m:[] for m in range(2,M+1)}
    for m,sources in reference.items():
        for n in sources:
            if n%a:
                continue
            pq = n//a
            fpq = factor(pq)
            if len(fpq) == 2 and all(e == 1 for e in fpq.values()) and gcd(a,pq) == 1:
                expected_by_m[m].append(n)
    for m in range(2,M+1):
        require(producer.two_prime_core(m,a) == expected_by_m[m], ('hyperbola family',a,m))
        hyperbola_cases += 1

# Hostile controls force survival of each exceptional coordinate.
require(aliquot(4)==3 and 4>2*3-6, 'even linear boundary must exclude4')
require(reference[5]==[], 'five untouchable')
require(121 in reference[12] and 121>2*12, 'prime square beyond linear cutoff')
require(32 in reference[31] and gcd(*pairs[32])>1, 'squarefull u,v noncoprime')
require(27 in reference[13] and isqrt(27)**2!=27, 'squarefull is not square')
require(aliquot(9)!=1+3+3, 'Goldbach repeated prime fails')
require(aliquot(4)!=aliquot(2)*2+(aliquot(2)+2), 'linear core needs gcd')
require(aliquot(2)==aliquot(3)==1 and aliquot(4)!=aliquot(6), 'fiber quotient not multiplicative')
require(producer.fiber(0)==[1] and producer.fiber(-1)==[], 'zero/negative targets')
try:
    producer.fiber(1)
except ValueError:
    pass
else:
    raise RuntimeError('infinite prime fiber must not be silently finite')

serialized = json.dumps(computed,sort_keys=True,separators=(',',':')).encode()
report = {
    'status':'PASS; finite-exact audit, no asymptotic conjecture closure',
    'total_explicit_checks':checks,
    'full_inverse_targets':[2,M], 'full_source_oracle':[1,B],
    'composite_bound_sources':bound_sources,
    'independent_factor_oracle_check':[1,10000],
    'squarefull_pairs_including1':len(pairs),
    'squarefull_bound':3*(M-1),
    'prefix_endpoints':list(range(0,41))+[99,100,121,500,999,1000],
    'hyperbola_a':[2,50], 'hyperbola_target':[2,M], 'hyperbola_cases':hyperbola_cases,
    'fiber_sha256':hashlib.sha256(serialized).hexdigest(),
    'untouchables_through1000':[m for m in range(2,M+1) if not reference[m]],
    'frozen_producer_sha256':hashlib.sha256((ROOT/'source-snapshot/aliquot_fibers.py').read_bytes()).hexdigest(),
    'blockers':[]
}
print(json.dumps(report,indent=2,sort_keys=True))

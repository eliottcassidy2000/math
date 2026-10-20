#!/usr/bin/env python3
"""Exact tests for Carmichael power retractions. Standard library only.
No Python asserts: identical checks remain active under python -O.
"""
from collections import Counter
from itertools import product
from math import comb, gcd, prod
import json

checks = Counter()

def check(condition, lane, detail=None):
    checks[lane] += 1
    if not condition:
        raise RuntimeError(f'{lane}: {detail}')

def factor(n):
    out = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0)+1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out

def vp(x,p):
    if x == 0:
        return 10**9
    v = 0
    while x % p == 0:
        v += 1
        x //= p
    return v

def floor_log(p,x):
    h = 0
    while p**(h+1) <= x:
        h += 1
    return h

def ceil_log(p,x):
    h = 0
    while p**h < x:
        h += 1
    return h

def hensel(a, exponent, p, k):
    """Digit correction using derivative exponent*t^(exponent-1)-1=-1 mod p."""
    t = a % p
    for j in range(1,k):
        residual = (pow(t, exponent, p**(j+1))-t)//(p**j)
        t += (residual % p)*p**j
    return t

def crt(values, moduli):
    M = prod(moduli)
    return sum(v*(M//m)*pow(M//m,-1,m) for v,m in zip(values,moduli)) % M

def section(a,n,k):
    return pow(a,n**(k-1),n**k)

def carry(a,b,n):
    m = n*n
    return ((pow(a+b,n,m)-pow(a,n,m)-pow(b,n,m)) % m)//n

def poly_mul(a,b,m):
    c = [0]*len(a)
    for i,x in enumerate(a):
        for j in range(len(a)-i):
            c[i+j] = (c[i+j]+x*b[j]) % m
    return tuple(c)

def poly_pow(a,e,m):
    ans = (1,)+(0,)*(len(a)-1)
    while e:
        if e&1:
            ans = poly_mul(ans,a,m)
        a = poly_mul(a,a,m)
        e >>= 1
    return ans

def multiplicative_order(p,l):
    z,k = p % l,1
    while z != 1:
        z = z*p % l
        k += 1
    return k

# Exact finite universe: every exponent 2..1200, all its residue classes.
# Early exit on the first failing residue is logically sufficient for a negative.
identity_n = []
for n in range(2,1201):
    fac = factor(n)
    korselt = all(e == 1 and (n-1) % (p-1) == 0 for p,e in fac.items())
    actual = all(pow(a,n,n) == a for a in range(n))
    check(actual == korselt,'pointwise_korselt',n)
    if actual:
        identity_n.append(n)

carmichaels = [561,1105,1729,2465,2821,6601]
for n in carmichaels:
    ps = list(factor(n))
    for a in range(n):
        check(pow(a,n,n) == a,'carmichael_residues',(n,a))
        # All residue classes, precisions 1..5, independent Hensel/CRT route.
        for k in range(1,6):
            local = [hensel(a,n,p,k) for p in ps]
            expected = crt(local,[p**k for p in ps])
            got = section(a,n,k)
            check(got == expected,'section_vs_hensel_crt',(n,a,k))
            check(pow(got,n,n**k) == got,'fixed_section',(n,a,k))
            if k > 1:
                check(got % n**(k-1) == section(a,n,k-1),'precision_compatibility')
            for p,t in zip(ps,local):
                check(t == hensel(a,p,p,k),'classical_teichmuller_identification')
    for k in range(2,10):
        m = n**k
        r = k-2
        before = pow(1+n,n**r,m)
        check(before != 1 and pow(before,n,m) == 1,'arithmetic_sharpness',(n,k))

# Exhaustive global n=561, k=2: all 314721 inputs, including zero-divisors.
n,k = 561,2
m = n*n
roots = []
for x in range(m):
    y = pow(x,n,m)
    check(y == section(x % n,n,2),'global_complete_retraction',x)
    check(pow(y,n,m) == y,'global_complete_idempotence',x)
    if y == x:
        roots.append(x)
check(len(roots) == n,'global_root_count',len(roots))

# Exact unit versus zero-residue hitting-time formula, complete local rings.
local_cases = [(3,561,5),(11,561,3),(17,561,3),(5,1105,4),(13,1105,3)]
for p,e,k in local_cases:
    m = p**k
    for x in range(m):
        t = hensel(x,e,p,k)
        s = min(k,vp(x-t,p))
        if x % p:
            expected = k-s
        else:
            expected = 0
            while e**expected*s < k:
                expected += 1
        y,hits = x,0
        while y != t:
            y = pow(y,e,m)
            hits += 1
            if hits > k:
                raise RuntimeError('unexpected nonconvergence')
        check(hits == expected,'exact_local_hitting_time',(p,e,k,x,hits,expected))

# Complete truncated polynomial rings for prime and nonprime exponents,
# including p=2 and higher v_p(e). No Carmichael assumption in this theorem.
jet_cases = [(3,561,1,5),(3,561,2,4),(3,561,3,2),
             (3,21,2,3),(5,1105,1,5),(5,45,2,2),
             (2,6,2,4),(2,6,3,3),(3,9,2,3),(2,4,2,4)]
jet_rows = []
for p,e,a,b in jet_cases:
    check((e-1) % (p-1) == 0 and e % p == 0,'jet_hypotheses')
    s = vp(e,p)
    h = floor_log(p,b-1)
    bound = (a+h+s-1)//s
    m = p**a
    worst = 0
    for x in product(range(m),repeat=b):
        target = (hensel(x[0],e,p,a),)+(0,)*(b-1)
        y,hits = x,0
        while y != target:
            y = poly_pow(y,e,m)
            hits += 1
            if hits > bound:
                raise RuntimeError(f'jet bound failed: {p,e,a,b,x}')
        worst = max(worst,hits)
        check(poly_pow(x,e**bound,m) == target,'complete_mixed_jet_ring')
    check(worst == bound,'mixed_jet_sharpness',(p,e,a,b,worst,bound))
    witness = (1,1)+(0,)*(b-2)
    check(poly_pow(witness,e**(bound-1),m) != (1,)+(0,)*(b-1),'one_plus_t_witness')
    jet_rows.append({'p':p,'exponent':e,'arithmetic_precision':a,'jet_length':b,
                     'ring_size':m**b,'worst_iterations':worst})

# Coefficient-path independent of convolution. Exact valuation of binomials
# establishes the sharp leading witness even at large jet lengths.
for p,e in [(3,561),(11,561),(17,561),(5,1105),(2,6),(3,9),(2,4)]:
    s = vp(e,p)
    for r in range(0,5):
        E = e**r
        for j in range(1,min(p**(s*r),90)+1):
            check(vp(comb(E,j),p) == s*r-vp(j,p),'binomial_valuation',(p,e,r,j))
    for a in range(1,5):
        for b in range(2,45):
            h = floor_log(p,b-1)
            r = (a+h+s-1)//s
            check(all(comb(e**r,j) % p**a == 0 for j in range(1,b)),
                  'all_jet_coefficients_vanish',(p,e,a,b,r))
            check(any(comb(e**(r-1),j) % p**a != 0 for j in range(1,b)),
                  'prior_jet_coefficient_survives',(p,e,a,b,r))

# Composite carry cocycle and local normalization; all 561^2 input pairs.
n = 561
ps = list(factor(n))
nonzero_carries = 0
for a in range(n):
    for b in range(n):
        c = carry(a,b,n)
        nonzero_carries += bool(c)
        for p in ps:
            expected = pow(n//p,-1,p)*carry(a % p,b % p,p) % p
            check(c % p == expected,'carry_crt_normalization',(a,b,p))
        if a < 25 and b < 25:
            for d in range(5):
                check((c+carry(a+b,d,n)-carry(b,d,n)-carry(a,b+d,n)) % n == 0,
                      'carry_cocycle')
check(nonzero_carries > 0,'carry_nonadditivity')
check(comb(n,3) % n != 0,'formal_freshmans_dream_hostile')

# First Hasse jet and first equal-characteristic nonadditivity.
jet_hostiles = []
for n in carmichaels:
    for p in factor(n):
        check(all(comb(n,j) % p == 0 for j in range(1,p)),'lower_hasse_jets_zero')
        check(comb(n,p) % p == (n//p) % p != 0,'first_hasse_jet')
        b = p+1
        one = (1,)+(0,)*(b-1)
        t = (0,1)+(0,)*(b-2)
        one_t = (1,1)+(0,)*(b-2)
        lhs = poly_pow(one_t,n,p)
        rhs = tuple((x+y) % p for x,y in zip(poly_pow(one,n,p),poly_pow(t,n,p)))
        check(lhs != rhs and lhs[p] == (n//p) % p,'nilpotent_nonadditivity')
        jet_hostiles.append({'n':n,'p':p,'first_hasse_order':p,'coefficient':lhs[p]})

# Finite-field counts are verified independently on actual field elements
# in three tiny extensions, using polynomial representatives and exact products.
def field_mul(x,y,p,f,relation):
    z = [0]*(2*f-1)
    for i,a in enumerate(x):
        for j,b in enumerate(y):
            z[i+j] = (z[i+j]+a*b)%p
    for d in range(2*f-2,f-1,-1):
        c = z[d]
        for j in range(f):
            z[d-f+j] = (z[d-f+j]+c*relation[j])%p
    return tuple(z[:f])

def field_pow(x,e,p,f,relation):
    z = (1,)+(0,)*(f-1)
    while e:
        if e&1:
            z = field_mul(z,x,p,f,relation)
        x = field_mul(x,x,p,f,relation)
        e >>= 1
    return z

field_rows = []
# T^2+1 over F3; T^3-T-1 over F3; T^2+1 over F11 are irreducible.
for p,f,rel in [(3,2,(2,0)),(3,3,(1,1,0)),(11,2,(10,0))]:
    q,n = p**f,561
    elems = list(product(range(p),repeat=f))
    image = [field_pow(x,n,p,f,rel) for x in elems]
    fixed = sum(x == y for x,y in zip(elems,image))
    check(fixed == 1+gcd(n-1,q-1),'extension_fixed_count',(p,f,fixed))
    check((len(set(image)) == q) == (gcd(n,q-1) == 1),'extension_permutation')
    check((fixed == q) == ((n-1) % (q-1) == 0),'extension_identity')
    # This additionally confirms the supplied quotient really is a field.
    for x in elems[1:]:
        check(field_pow(x,q-1,p,f,rel) == (1,)+(0,)*(f-1),'extension_field_control')
    field_rows.append({'p':p,'f':f,'field_size':q,'fixed_points':fixed,'image_size':len(set(image))})

first_hidden = {}
for p in factor(561):
    first_hidden[str(p)] = min(multiplicative_order(p,l) for l in factor(561//p))

print(json.dumps({'status':'PASS','total_checks':sum(checks.values()),
    'checks':dict(sorted(checks.items())),
    'universes':{'pointwise_korselt_exponents':[2,1200], 'carmichaels':carmichaels,
                 'global_complete_ring':'Z/561^2', 'global_inputs':561**2,
                 'section_precisions':[1,5],'local_cases':local_cases,
                 'composite_carry_pairs':561**2},
    'mixed_jet_cases':jet_rows,'finite_field_cases':field_rows,
    'first_noninjective_extension_degree_for_561':first_hidden,
    'nonzero_composite_carry_pairs':nonzero_carries,
    'first_hasse_jet_examples':jet_hostiles},indent=2,sort_keys=True))

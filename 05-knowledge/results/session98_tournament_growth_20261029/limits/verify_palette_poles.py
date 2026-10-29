#!/usr/bin/env python3
"""Self-contained exact palette theorem controls; explicit non-removable gates."""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product
from math import comb, factorial, gcd
from pathlib import Path
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def partitions(n, least=1):
    if n == 0:
        yield ()
    for d in range(least, n + 1):
        for tail in partitions(n-d, d):
            yield (d,) + tail


def centralizer(mu):
    result = 1
    for d, c in Counter(mu).items():
        result *= d**c * factorial(c)
    return result


def fixed_formula(mu, q, t):
    u = sum((d-1)//2 for d in mu)
    u += sum(gcd(x,y) for i,x in enumerate(mu) for y in mu[i+1:])
    v = sum(d % 2 == 0 for d in mu)
    return q**u * t**v


def cycle_type(p):
    seen = set()
    cycles = []
    for i in range(len(p)):
        if i in seen:
            continue
        length = 0
        while i not in seen:
            seen.add(i)
            length += 1
            i = p[i]
        cycles.append(length)
    return tuple(sorted(cycles))


def exhaustive(n, q, t):
    tau = list(range(t))
    for i in range(t, q, 2):
        tau += [i+1, i]
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    index = {pair:k for k,pair in enumerate(pairs)}
    ps = list(permutations(range(n)))
    maps = [[(index[tuple(sorted((p[i],p[j])))], p[i]>p[j]) for i,j in pairs] for p in ps]
    fixed = [0] * len(ps)
    canonical = set()
    for colors in product(range(q), repeat=len(pairs)):
        images = []
        for h, mapping in enumerate(maps):
            image = [None] * len(pairs)
            for r,(dest,reverse) in enumerate(mapping):
                image[dest] = tau[colors[r]] if reverse else colors[r]
            image = tuple(image)
            images.append(image)
            if image == colors:
                fixed[h] += 1
        canonical.add(min(images))
    for p,value in zip(ps,fixed):
        require(value == fixed_formula(cycle_type(p),q,t), f'fixed count mismatch {n,q,t,p}')
    require(sum(fixed) == factorial(n)*len(canonical), 'direct orbit/Burnside mismatch')
    formula = sum((F(fixed_formula(mu,q,t),centralizer(mu)) for mu in partitions(n)), F(0))
    require(formula == len(canonical), 'partition count/canonical orbit mismatch')
    return len(canonical)


def core(mu,q,t):
    m = sum(mu)
    a = m-len(mu)
    fixed = fixed_formula(mu,q,t)
    # Fourth powers make the quadratic inequality exact integer arithmetic.
    require(fixed**4 * q**(m*(m-2)) <= q**(4*comb(m,2)), 'quadratic core bound')
    require(2*a >= m, 'core defect bound')
    C = F(factorial(m)*fixed,centralizer(mu)*q**comb(m,2))
    return m,a,C


palettes = [(q,t) for q in range(2,5) for t in range(q+1) if (q-t)%2 == 0]
receipts = []
for q,t in palettes:
    values = [exhaustive(n,q,t) for n in range(5)]
    for n in range(13):
        full = F(factorial(n), q**comb(n,2))*sum((F(fixed_formula(mu,q,t),centralizer(mu)) for mu in partitions(n)),F(0))
        extracted = F(0)
        for m in range(n+1):
            for mu in partitions(m,2):
                _,a,C = core(mu,q,t)
                extracted += C*comb(n,m)*F(1,q**(a*(n-m)))
        require(full == extracted, f'core extraction {n,q,t}')
    for k in range(7):
        top = (2,)*k if t else (3,)*k
        m,a,C = core(top,q,t)
        expected = F(factorial(2*k),2**k*factorial(k))*t**k*q**(k*k) if t else F(factorial(3*k),3**k*factorial(k))*q**(3*k*k+k)
        require(C*q**(a*m) == expected, 'leading amplitude')
        target = k if t else 2*k
        bank = [mu for size in range(m+1) for mu in partitions(size,2) if size-len(mu)==target and fixed_formula(mu,q,t)]
        require([mu for mu in bank if sum(mu)==m] == [top], 'unique maximal core')
        require(len(bank) == len(list(partitions(k))), 'bank partition count')
    receipts.append({'q':q,'t':t,'unlabeled_n_0_to_4':values,'fixed_count_exhaustive_n_max':4,'core_coefficient_n_max':12,'pole_k_max':6})

report = {'status':'PASS','palettes':receipts,'arithmetic':'exact integers and Fractions','limits':'finite checks do not prove global convergence or infinite pole sets'}
Path(__file__).with_name('palette_poles_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

"""Independent direct subset DP and marked-endpoint monomial compiler."""
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from math import factorial, prod


def hp(n, rev):
    @lru_cache(None)
    def f(m, v):
        if m == 1 << v:
            return 1
        return sum(f(m ^ (1 << v), i) for i in range(n)
                   if i != v and (m >> i) & 1
                   and ((i < v) ^ ((min(i, v), max(i, v)) in rev)))
    return sum(f((1 << n) - 1, i) for i in range(n))


def inventory(k, rev, gaps):
    out = Counter()
    for p in permutations(range(k)):
        seq = (-1,) + p + (k,)
        slots, forced = [], []
        for u, v in zip(seq, seq[1:]):
            if u > v and (v, u) not in rev:
                break
            if u < v:
                row = tuple(j for j, g in enumerate(gaps) if u < g <= v)
                slots.append(row)
                if (u, v) in rev:
                    forced.append(row)
        else:
            c = [sum(j in row for row in slots) for j in range(len(gaps))]
            for mask in range(1 << len(forced)):
                d = [x - sum((mask >> i) & 1 for i, row in enumerate(forced)
                             if j in row) for j, x in enumerate(c)]
                out[tuple(d)] += (-1)**mask.bit_count()
    return {k: v for k, v in out.items() if v}


def evaluate(coeffs, exponents):
    return sum(c * prod(b**e for b, e in zip(key, exponents))
               for key, c in coeffs.items())


def closed(a, b, c):
    x = (2**(a+c+1)*3**b + (2**(a+1)+2**(c+1))*2**b + 1
         - (a == 0)*2**c - (c == 0)*2**a)
    y = (2**(a+c+1)*3**b + (2**(a+c+2)+1)*2**b + 1
         - 2**(a+c+1) - (b == 0))
    return x, y


X, N = {(0, 2), (1, 3)}, {(0, 3), (1, 2)}
models = [inventory(4, rev, (1, 2, 3)) for rev in (X, N)]
expected = [{(1,1,1):1,(1,2,2):2,(0,1,2):-1,(2,2,1):2,(2,1,0):-1,(2,3,2):2},
            {(1,1,1):1,(1,0,1):-1,(1,2,1):1,(2,2,2):4,(2,1,2):-2,(2,3,2):2}]
if models != expected:
    raise ValueError('monomial table')
for a, b, c in product(range(4), repeat=3):
    endpoints = (0, a+1, a+b+2, a+b+c+3)
    counts = tuple(hp(endpoints[-1]+1, {(endpoints[i],endpoints[j]) for i,j in rev})
                   for rev in (X,N))
    if counts != closed(a,b,c) or counts != tuple(evaluate(m,(a,b,c)) for m in models):
        raise ValueError((a,b,c,counts))
for r in range(1,5):
    for p in permutations(range(r)):
        m = inventory(2*r, {(l,r+v) for l,v in enumerate(p)}, (r,))
        if m[(r+1,)] != factorial(r) or max(k[0] for k in m) != r+1:
            raise ValueError((r,p,m))
for b in range(1,100):
    if closed(0,b,0)[1] - closed(0,b,0)[0] != 2**b:
        raise ValueError('binary carry')
    for u in range(2,10):
        if closed(u-1,b+1,u+2)[0] != 3*closed(u,b,u+1)[0] - 2:
            raise ValueError('width exchange')
print('128 subset DP counts, 33 perfect matchings, and exact carry/exchange controls passed')

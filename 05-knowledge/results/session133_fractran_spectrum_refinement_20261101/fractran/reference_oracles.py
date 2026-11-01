"""Frozen independent direct-permutation/orbit and partition oracles from Session132.
Trimmed to the imported standard-library-only checks; not inputs to any program.
"""

from collections import Counter

from fractions import Fraction

from functools import lru_cache

from itertools import combinations, permutations, product

from math import factorial, gcd

from pathlib import Path

import json

def require(test, message):
    if not test:
        raise RuntimeError(message)

def odd_partitions(n, least=1):
    if n == 0:
        yield ()
    for d in range(least, n + 1, 2):
        for tail in odd_partitions(n-d, d):
            yield (d,) + tail

def cycle_data(parts):
    multiplicities = Counter(parts)
    z = 1
    for d, m in multiplicities.items():
        z *= d ** m * factorial(m)
    exponent = sum((d-1)//2 for d in parts)
    exponent += sum(gcd(d, e) for d, e in combinations(parts, 2))
    return z, exponent, len(parts)

def partition_count(n, q=2, x=1):
    """F_n(q,x), with integer summands in n! F_n(q,x)."""
    total = 0
    fact = factorial(n)
    num_partitions = 0
    for parts in odd_partitions(n):
        z, exponent, r = cycle_data(parts)
        quotient, rem = divmod(fact, z)
        require(rem == 0, 'cycle conjugacy class size must be integral')
        total += quotient * q**exponent * x**r
        num_partitions += 1
    value, rem = divmod(total, fact)
    require(rem == 0, 'Burnside average must be integral')
    return value, num_partitions

def pair_action(n, p):
    pairs = list(combinations(range(n), 2))
    positions = {edge:i for i,edge in enumerate(pairs)}
    out = []
    for u,v in pairs:
        a,b = p[u],p[v]
        out.append((positions[tuple(sorted((a,b)))], int(a>b)))
    return out

def fixed_count_from_signed_pair_cycles(action, converse=False):
    unseen = set(range(len(action)))
    free = 0
    while unseen:
        start = next(iter(unseen))
        i = start
        parity = 0
        while i in unseen:
            unseen.remove(i)
            i, swap = action[i]
            parity ^= swap ^ converse
        require(i == start, 'pair action is a permutation')
        if parity:
            return 0
        free += 1
    return 2**free

def permutation_burnside(n):
    plain = twisted = 0
    for p in permutations(range(n)):
        action = pair_action(n,p)
        plain += fixed_count_from_signed_pair_cycles(action)
        twisted += fixed_count_from_signed_pair_cycles(action,True)
    a, ra = divmod(plain, factorial(n))
    s, rs = divmod(twisted, factorial(n))
    require(ra == rs == 0, 'direct signed permutation Burnside integral')
    return a,s

def transform(mask, action):
    out = 0
    for i,(j,swap) in enumerate(action):
        out |= (((mask>>i)&1)^swap)<<j
    return out

def permute_subset(mask, p):
    out = 0
    for i,j in enumerate(p):
        out |= ((mask>>i)&1)<<j
    return out

def orbit_number(n, group):
    return len({min(p[v] for p in group) for v in range(n)})

def strong(n, mask):
    if n == 0:
        return False
    reach = [1<<i for i in range(n)]
    for k,(u,v) in enumerate(combinations(range(n),2)):
        if (mask>>k)&1:
            reach[u] |= 1<<v
        else:
            reach[v] |= 1<<u
    for k in range(n):
        for i in range(n):
            if (reach[i]>>k)&1:
                reach[i] |= reach[k]
    return all(x == (1<<n)-1 for x in reach)

class Census:
    def __init__(self,n):
        self.n = n
        self.ps = list(permutations(range(n)))
        self.actions = [pair_action(n,p) for p in self.ps]
        self.full = (1<<(n*(n-1)//2))-1
        unseen = set(range(self.full+1))
        self.reps = []
        self.canonical = {}
        while unseen:
            t = min(unseen)
            orb = {transform(t,a) for a in self.actions}
            for u in orb:
                self.canonical[u] = t
            unseen -= orb
            self.reps.append(t)
        self.sc = [t for t in self.reps if self.canonical[t^self.full] == t]
        self.merged = [t for t in self.reps if t <= self.canonical[t^self.full]]
        self.auts = {}
        self.antis = {}
        for t in self.reps:
            self.auts[t] = [p for p,a in zip(self.ps,self.actions) if transform(t,a) == t]
            self.antis[t] = [p for p,a in zip(self.ps,self.actions) if transform(t,a) == (t^self.full)]
            require(len(self.auts[t])%2 == 1, 'tournament automorphism order odd')
            require(len(self.antis[t]) in (0,len(self.auts[t])), 'anti-isomorphism torsor')

    def merged_canonical(self,t):
        return min(self.canonical[t],self.canonical[t^self.full])

    def summary(self):
        return dict(n=self.n, labelled=self.full+1, A=len(self.reps), S=len(self.sc),
                    M=len(self.merged), C=sum(strong(self.n,t) for t in self.reps),
                    D=sum(strong(self.n,t) for t in self.sc),
                    automorphism_order_histogram=dict(sorted(Counter(len(self.auts[t]) for t in self.reps).items())),
                    rooted=sum(orbit_number(self.n,self.auts[t]) for t in self.reps),
                    merged_rooted=sum(orbit_number(self.n,self.auts[t]+self.antis[t]) for t in self.merged))

def pair_tournament_census(n, colors=1):
    """Independent orbit enumeration of two tournament layers, optional 2-coloring."""
    require(colors in (1,2), 'supported direct census coloring')
    ps = list(permutations(range(n)))
    acts = [pair_action(n,p) for p in ps]
    tmax = 1<<(n*(n-1)//2)
    cmax = 1 if colors == 1 else 1<<n
    ttrans = [[transform(t,a) for a in acts] for t in range(tmax)]
    ctrans = [[permute_subset(c,p) for p in ps] for c in range(cmax)]
    unseen = set(product(range(tmax),range(tmax),range(cmax)))
    count = 0
    while unseen:
        a,b,c = next(iter(unseen))
        unseen -= {(ttrans[a][k],ttrans[b][k],ctrans[c][k]) for k in range(len(ps))}
        count += 1
    return count

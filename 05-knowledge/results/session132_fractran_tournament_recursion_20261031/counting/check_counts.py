#!/usr/bin/env python3
"""Exact tournament counters; independent cycle, permutation, and orbit routes.

Python 3 standard library. No assert is used: -O retains every check.
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


def cycle_insertion(n, q=2, x=1):
    """Integer recurrence, tracking all cross-cycle gcd interaction exponents."""
    @lru_cache(None)
    def B(k, h):
        if k == 0:
            return 1
        total = 0
        for d in range(1, k+1, 2):
            choices = factorial(k-1)//factorial(k-d)
            exponent = (d-1)//2 + h[(d-1)//2]
            child_h = tuple(h[(j-1)//2]+gcd(d,j)
                            for j in range(1, k-d+1, 2))
            total += choices * x * q**exponent * B(k-d, child_h)
        return total
    numerator = B(n, (0,)*((n+1)//2))
    value, rem = divmod(numerator, factorial(n))
    require(rem == 0, 'integer recursive average')
    return value, B.cache_info().currsize


def twisted_partition_count(n):
    """Independent doubled-part cycle formula in S_n, not F_{floor(n/2)}."""
    singleton = n % 2
    total = 0
    for odd_parts in odd_partitions(n//2):
        parts = tuple(2*d for d in odd_parts) + (1,)*singleton
        counts = Counter(parts)
        z = 1
        for d, m in counts.items():
            z *= d**m * factorial(m)
        e = sum(d//2 for d in parts)
        e += sum(gcd(d,f) for d,f in combinations(parts,2))
        q, r = divmod(factorial(n),z)
        require(r == 0, 'twisted conjugacy class size')
        total += q * 2**e
    q, r = divmod(total, factorial(n))
    require(r == 0, 'twisted average integral')
    return q


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


def extension_mask(n, t, subset):
    """New vertex n dominates exactly subset of old vertices."""
    old = {e:i for i,e in enumerate(combinations(range(n),2))}
    u = 0
    for i,(a,b) in enumerate(combinations(range(n+1),2)):
        bit = 1-((subset>>a)&1) if b == n else ((t>>old[(a,b)])&1)
        u |= bit<<i
    return u


def verify_extension(parent,child,merged=False):
    n = parent.n
    total = 0
    unpointed = Fraction(0)
    targets = Counter()
    for t in (parent.merged if merged else parent.reps):
        full = (1<<n)-1
        unseen = set(range(1<<n))
        while unseen:
            subset = min(unseen)
            orb = {permute_subset(subset,p) for p in parent.auts[t]}
            if merged:
                orb |= {permute_subset(subset^full,p) for p in parent.antis[t]}
            unseen -= orb
            total += 1
            u = extension_mask(n,t,subset)
            rep = child.merged_canonical(u) if merged else child.canonical[u]
            group = child.auts[rep] + (child.antis[rep] if merged else [])
            vertex_orbits = orbit_number(n+1,group)
            unpointed += Fraction(1,vertex_orbits)
            targets[rep] += 1
    expected = child.summary()['merged_rooted' if merged else 'rooted']
    require(total == expected, 'extension roots vs vertex orbits')
    require(unpointed == len(child.merged if merged else child.reps), 'weighted unpointing')
    for rep,mult in targets.items():
        group = child.auts[rep] + (child.antis[rep] if merged else [])
        require(mult == orbit_number(n+1,group), 'every child once for each vertex orbit')
    return dict(parent_n=n, merged=merged, extension_orbits=total,
                weighted_unpointed=int(unpointed),
                naive_divide_by_nplus1=str(Fraction(total,n+1)))


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


def main():
    A = []; S = []; M = []; rows = []
    for n in range(26):
        a,partitions = partition_count(n)
        ar,states = cycle_insertion(n)
        require(a == ar, 'partition vs recursive A '+str(n))
        s,_ = partition_count(n//2,4,2**(n%2))
        sr,sstates = cycle_insertion(n//2,4,2**(n%2))
        require(s == sr == twisted_partition_count(n), 'three independent S routes '+str(n))
        require((a+s)%2 == 0, 'merged parity')
        A.append(a); S.append(s); M.append((a+s)//2)
        rows.append(dict(n=n,A=a,S=s,M=(a+s)//2,
                         odd_partitions=partitions,recursive_states=states,
                         self_converse_recursive_states=sstates))
    # Published OEIS data: independent external controls, not computation inputs.
    pinned = json.loads((Path(__file__).parent/'sources'/'oeis_values.json').read_text())
    for name,computed in [('A000568',A),('A002785',S),('A059735',M)]:
        for n,value in pinned[name].items():
            require(computed[int(n)] == value, 'OEIS '+name+' n='+n)
    C = [0]*26; D = [0]*26
    for n in range(1,26):
        C[n] = A[n] - sum(C[k]*A[n-k] for k in range(1,n))
        D[n] = S[n] - (A[n//2] if n%2 == 0 else 0)
        D[n] -= sum(D[j]*A[(n-j)//2] for j in range(1,n) if (n-j)%2 == 0)
        require(C[n] >= D[n] >= 0, 'strong-count positivity')
        rows[n].update(C=C[n],D=D[n],strong_merged=(C[n]+D[n])//2)
        require((C[n]+D[n])%2 == 0, 'strong merged parity')
    direct_burnside = []
    for n in range(9):
        a,s = permutation_burnside(n)
        require((a,s) == (A[n],S[n]), 'signed pair-cycle direct Burnside n='+str(n))
        direct_burnside.append(dict(n=n,A=a,S=s,permutations=factorial(n)))
    censuses = [Census(n) for n in range(6)]
    summaries = [c.summary() for c in censuses]
    for row in summaries:
        n = row['n']
        require([row[k] for k in ('A','S','M','C','D')] == [A[n],S[n],M[n],C[n],D[n]],
                'full independent isomorphism census '+str(n))
    extensions = [verify_extension(censuses[n],censuses[n+1],merged)
                  for n in range(5) for merged in (False,True)]
    pair_checks = []
    for n in range(5):
        for colors in (1,2):
            actual = pair_tournament_census(n,colors)
            expected = partition_count(n,4,colors)[0]
            require(actual == expected, 'ordered pair orbit census')
            pair_checks.append(dict(n=n,vertex_colors=colors,count=actual))
    print(json.dumps(dict(status='PASS',range='0..25 inclusive',counts=rows,
                         direct_signed_permutation_burnside=direct_burnside,
                         independent_tournament_census=summaries,
                         pointed_extension_checks=extensions,
                         independent_pair_tournament_census=pair_checks,
                         oeis_compared_terms={k:len(v) for k,v in pinned.items()}),indent=2,sort_keys=True))


if __name__ == '__main__':
    main()

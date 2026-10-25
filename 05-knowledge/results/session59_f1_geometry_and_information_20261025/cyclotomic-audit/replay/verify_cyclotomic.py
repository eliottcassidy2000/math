#!/usr/bin/env python3
"""Exact independent coordinate enumeration versus derived cycle formulas.

No third-party dependencies; no assertions (normal and -O runs are identical).
The distinguished Python value None denotes monoid zero; integer 0 is a unit.
"""
from collections import Counter
from itertools import combinations, permutations, product
from math import comb, gcd, isqrt, lcm
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def N(k, m):
    require(k >= 0 and m >= 1, 'N domain')
    return ((m + 1) ** k - 1) // m


def normalize(x, m):
    anchor = next((a for a in x if a is not None), None)
    require(anchor is not None, 'zero vector cannot be projectivized')
    return tuple(None if a is None else (a - anchor) % m for a in x)


def points(k, m):
    # Independent construction: normalize each raw nonzero homogeneous vector.
    return sorted({normalize(x, m) for x in product([None] + list(range(m)), repeat=k)
                   if any(a is not None for a in x)},
                  key=lambda x: tuple(-1 if a is None else a for a in x))


def act(x, m, b, sigma, r):
    y = [None] * len(x)
    for j, i in enumerate(sigma):
        if x[j] is not None:
            y[i] = (b[i] + r * x[j]) % m
    return normalize(y, m)


def cycles(sigma):
    unseen = set(range(len(sigma)))
    result = []
    while unseen:
        start = min(unseen)
        c, i = [], start
        while i in unseen:
            unseen.remove(i)
            c.append(i)
            i = sigma[i]
        result.append(c)
    return result


def cycle_fixed(k, m, b, sigma, r):
    # For C=(i0,...,i_{ell-1}), recurse i0->i1->...->i0.
    data = []
    for c in cycles(sigma):
        order = c[1:] + c[:1]
        B, A = 0, 0
        for i in order:
            B = r * B + b[i]
            A = r * A + 1
        g = gcd(m, r ** len(c) - 1)
        data.append((g, A, B))
    total = 0
    for t in range(m):
        v = 1
        for g, A, B in data:
            v *= 1 + (g if (A * t - B) % g == 0 else 0)
        total += v - 1
    require(total % m == 0, 'fixed quotient not integral')
    return total // m


def compose(left, right, m):
    # left after right: (b,sigma,r)(c,tau,s)=(b+r sigma(c),sigma tau,rs)
    b, sigma, r = left
    c, tau, s = right
    pc = [0] * len(b)
    for j, i in enumerate(sigma):
        pc[i] = c[j]
    return (tuple((b[i] + r * pc[i]) % m for i in range(len(b))),
            tuple(sigma[tau[i]] for i in range(len(b))), (r * s) % m)


def phi(n):
    return sum(gcd(a, n) == 1 for a in range(1, n + 1))


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def mobius(n):
    result, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            result = -result
            if n % p == 0:
                return 0
            while n % p == 0:
                n //= p
        p += 1
    return -result if n > 1 else result


def jordan(d, n):
    return sum(mobius(d // e) * e ** n for e in divisors(d))


def unit_orbit_formula(k, m):
    a = sum(N(k, gcd(m, r - 1)) for r in range(m) if gcd(r, m) == 1)
    require(a % phi(m) == 0, 'unit Burnside not integral')
    burnside = a // phi(m)
    total = 0
    for s in range(1, k + 1):
        for d in divisors(m):
            j = jordan(d, s - 1)
            require(j % phi(d) == 0, 'Jordan orbit count not integral')
            total += comb(k, s) * j // phi(d)
    require(burnside == total, 'Jordan and Burnside mismatch')
    return total


def orbit_count(pts, maps):
    unseen = set(pts)
    total = 0
    while unseen:
        total += 1
        todo = [unseen.pop()]
        while todo:
            x = todo.pop()
            for f in maps:
                y = f(x)
                if y in unseen:
                    unseen.remove(y)
                    todo.append(y)
    return total


def tail_period(x, f):
    seen, y = {}, x
    while y not in seen:
        seen[y] = len(seen)
        y = f(y)
    return seen[y], len(seen) - seen[y]


def exact_tail_period(k, m, r, t, n):
    ans = 0
    for d in divisors(n):
        v = N(k, gcd(m, r ** t * (r ** d - 1)))
        if t:
            v -= N(k, gcd(m, r ** (t - 1) * (r ** d - 1)))
        ans += mobius(n // d) * v
    return ans


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def power(z, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, z)
    return out


def gap(js, primes, exponents):
    z = (1, 0)
    for j, pi, a in zip(js, primes, exponents):
        z = mul(z, mul(power(pi, a + j), power((pi[0], -pi[1]), a - j)))
    return abs(2 * z[0] * z[1])


def main():
    stats = Counter()
    samples = {}
    # Counts, support fibers, power fixed points, canonical unit action.
    for m in range(1, 9):
        for k in range(1, 6):
            pts = points(k, m)
            require(len(pts) == N(k, m), f'count {k,m}')
            support = Counter(tuple(i for i, a in enumerate(x) if a is not None) for x in pts)
            for S, count in support.items():
                require(count == m ** (len(S) - 1), f'support {k,m,S}')
            stats['count_cases'] += 1
            for r in range(0, 2 * m + 2):
                f = lambda x: act(x, m, (0,) * k, tuple(range(k)), r)
                actual = sum(f(x) == x for x in pts)
                expected = N(k, gcd(m, r - 1))
                require(actual == expected, f'power fixed {k,m,r}')
                stats['power_fixed_cases'] += 1
            if k <= 4:
                maps = [lambda x, r=r: act(x, m, (0,) * k, tuple(range(k)), r)
                        for r in range(m) if gcd(r, m) == 1]
                actual = orbit_count(pts, maps)
                require(actual == unit_orbit_formula(k, m), f'unit orbit {k,m}')
                stats['unit_orbit_cases'] += 1
    # Exhaust ALL diagonal twists, permutations and powers in this bounded grid.
    for m in range(1, 5):
        for k in range(1, 4):
            pts = points(k, m)
            for sigma in permutations(range(k)):
                for b in product(range(m), repeat=k):
                    for r in range(m + 1):
                        actual = sum(act(x, m, b, sigma, r) == x for x in pts)
                        require(actual == cycle_fixed(k, m, b, sigma, r), f'twist fixed {k,m,b,sigma,r}')
                        stats['exhaustive_twist_cases'] += 1
    # Larger composite moduli, 4-coordinate permutations, deterministic twists.
    for m in range(1, 9):
        k, pts = 4, points(4, m)
        for sigma in permutations(range(k)):
            for b in [(0,0,0,0), tuple(i % m for i in range(k)),
                      tuple((i*i+1) % m for i in range(k))]:
                for r in range(m + 1):
                    actual = sum(act(x, m, b, sigma, r) == x for x in pts)
                    require(actual == cycle_fixed(k, m, b, sigma, r), f'large twist {m,b,sigma,r}')
                    stats['larger_twist_cases'] += 1
    # Semilinear composition directly on points, and cyclic orbit Burnside.
    for m in range(2, 9):
        k, pts = 3, points(3, m)
        for r in range(m):
            left = ((0,1,2), (1,2,0), r)
            right = ((1,0,1), (1,0,2), 2)
            joined = compose(left, right, m)
            for x in pts:
                a = act(act(x, m, *right), m, *left)
                require(a == act(x, m, *joined), f'composition {m,r,x}')
                stats['composition_point_cases'] += 1
            if gcd(r,m) != 1:
                continue
            f = lambda x: act(x, m, *left)
            periods = Counter(tail_period(x, f)[1] for x in pts)
            order = lcm(*periods)
            identity = ((0,) * k, tuple(range(k)), 1 % m)
            current, total = identity, 0
            for _ in range(order):
                total += cycle_fixed(k, m, *current)
                current = compose(left, current, m)
            require(total % order == 0, 'cyclic Burnside integral')
            oc = orbit_count(pts, [f])
            require(total // order == oc, f'cyclic Burnside {m,r}')
            stats['cyclic_burnside_cases'] += 1
            if m in (4,6,8) and r == m - 1:
                samples[f'twisted_m{m}_k3_r{r}'] = {'b':[0,1,2], 'sigma':[1,2,0], 'order':order,
                    'fixed':cycle_fixed(k,m,*left), 'orbit_count':oc,
                    'cycle_counts':{str(n):count//n for n,count in sorted(periods.items())}}
    # Full functional graphs test exact tails/periods, not just bounded iterates.
    for m in range(1, 13):
        for k in range(1, 4):
            pts = points(k, m)
            for r in range(1, 2*m+1):
                f = lambda x: act(x, m, (0,) * k, tuple(range(k)), r)
                profile = Counter(tail_period(x, f) for x in pts)
                maxt = max(t for t,n in profile)
                maxn = max(n for t,n in profile)
                for t in range(maxt + 2):
                    for n in range(1, maxn + 2):
                        require(profile[t,n] == exact_tail_period(k,m,r,t,n),
                                f'tail period {k,m,r,t,n}')
                        stats['tail_period_cells'] += 1
                stats['functional_graph_cases'] += 1
                if (m,k,r) in [(12,3,2),(8,3,3),(6,3,2)]:
                    samples[f'dynamics_m{m}_k{k}_r{r}'] = {f'{t},{n}':v for (t,n),v in sorted(profile.items())}
    # Twisted finite-field descent: equal totals, unequal support fibers.
    for q in range(2, 6):
        for k in range(1, 4):
            for sigma in permutations(range(k)):
                d = lcm(*(len(c) for c in cycles(sigma)))
                m = q ** d - 1
                require(cycle_fixed(k,m,(0,)*k,sigma,q) == (q**k-1)//(q-1),
                        f'twisted descent {q,k,sigma}')
                stats['twisted_descent_formula_cases'] += 1
                if m <= 8:
                    pts = points(k,m)
                    actual = sum(act(x,m,(0,)*k,sigma,q) == x for x in pts)
                    require(actual == (q**k-1)//(q-1), 'twisted descent raw')
                    stats['twisted_descent_raw_cases'] += 1
    descent = {}
    for sigma in [(0,1),(1,0)]:
        pts = points(2,8)
        profile = Counter(sum(a is not None for a in x) for x in pts
                          if act(x,8,(0,0),sigma,3) == x)
        descent[str(sigma)] = dict(sorted(profile.items()))
    samples['q3_quadratic_descent_support_sizes'] = descent
    # Coefficient map multiplicity: image of a generator must have order dividing n.
    for n in range(1,13):
        for m in range(1,13):
            require(sum(n*a % m == 0 for a in range(m)) == gcd(n,m), 'coefficient homs')
            stats['coefficient_hom_cases'] += 1
    # Arithmetic bridge controls inherited Gaussian formulas, not new search.
    pis = [(2,1),(3,2)]
    g65 = {str(j):gap(j,pis,[1,1]) for j in [(1,0),(0,1),(1,1),(1,-1)]}
    require(list(g65.values()) == [4056,3000,3696,2016], 'e65 values')
    require(gap([1],[(2,1)],[2]) == 600, 'e25 j1')
    require(gap([2],[(2,1)],[2]) == 336, 'e25 j2')
    require(all((e*e-d)>=0 and isqrt(e*e-d)**2 == e*e-d and
                isqrt(e*e+d)**2 == e*e+d
                for e,d in [(65,d) for d in g65.values()] + [(25,600),(25,336)]), 'endpoint control')
    samples['arithmetic_e65'] = g65
    samples['higher_exponent_e25'] = {'j1':600, 'j2':336, 'gap_count':2, 'C4_projective_zero_space_count':1}
    samples['unit_orbit_table_k3'] = {str(m):unit_orbit_formula(3,m) for m in range(1,17)}
    report = {'status':'PASS', 'method':'exact enumeration with explicit RuntimeError gates; no assert statements',
              'statistics':dict(sorted(stats.items())), 'samples':samples}
    payload = json.dumps(report, indent=2, sort_keys=True) + '\n'
    (ROOT / 'verification.json').write_text(payload)
    print(payload, end='')

if __name__ == '__main__':
    main()

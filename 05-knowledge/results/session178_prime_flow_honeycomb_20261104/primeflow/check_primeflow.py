"""Standard-library exact checks; run from any directory. No network required."""

def require(condition, message='exact check failed'):
    if not condition:
        raise AssertionError(message)
from itertools import combinations
from math import lcm, prod
from fractions import Fraction
import json
from pathlib import Path

def primes(limit):
    return [p for p in range(3, limit + 1, 2) if all((p % d for d in range(2, int(p ** 0.5) + 1)))]

def valuation(n, p):
    require(n > 0)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def order(p, r):
    for d in range(1, r):
        if pow(p, d, r) == 1:
            return d
    raise AssertionError()

def ilog(n, p):
    require(n >= 1)
    t = 0
    while n >= p:
        n //= p
        t += 1
    return t

def sigma(p, a):
    return (p ** (a + 1) - 1) // (p - 1)

def graph(S, q):
    E = {}
    for p in S:
        for r in S:
            if p == r:
                continue
            d = order(p, r)
            if p == q and d % 4 == 0 or (p != q and d % 2 == 0):
                continue
            b = 2 if p == q else 1
            E[p, r] = {'d': d, 'k': r if d == 1 else d, 'ell': lcm(b, r if d == 1 else d), 'c': 0 if d == 1 else valuation(p ** d - 1, r), 'b': b}
    return E

def connected(S, E):
    for p in S:
        reach = {p}
        while True:
            new = reach | {r for t, r in E if t in reach}
            if new == reach:
                break
            reach = new
        if reach != set(S):
            return False
    return True

def capacity(z, U, r):
    if z['ell'] > U + 1:
        return 0
    if z['d'] == 1:
        return ilog((U + 1) // z['b'], r)
    return z['c'] + ilog((U + 1) // z['ell'], r)

def certificate(S, q):
    E = graph(S, q)
    require(connected(S, E))
    C = max((sum((z['c'] for (p, t), z in E.items() if t == r)) for r in S))
    B = 2
    while C + (len(S) - 1) * (ilog(B + 1, 3) + 1) >= B:
        B *= 2
    U = dict.fromkeys(S, B - 1)
    rounds = []
    while True:
        E = {e: z for e, z in E.items() if z['ell'] <= U[e[0]] + 1}
        if not connected(S, E):
            return {'q': q, 'C': C, 'B': B, 'rounds': rounds, 'reason': 'not strongly connected', 'remaining_edges': sorted(E)}
        caps = {r: sum((capacity(z, U[p], r) for (p, t), z in E.items() if t == r)) for r in S}
        V = {p: min(U[p], 1 + 4 * ((caps[p] - 1) // 4) if p == q else 2 * (caps[p] // 2)) for p in S}
        rounds.append({'before': U, 'raw_column_caps': caps, 'after_role_rounding': V})
        if any((V[p] < (1 if p == q else 2) for p in S)):
            return {'q': q, 'C': C, 'B': B, 'rounds': rounds, 'reason': 'positive role exponent impossible'}
        require(V != U, 'This particular example should be rejected')
        U = V

def run():
    P = primes(71)
    counts = {'lte': 0, 'capacities': 0, 'lcm_patterns': 0}
    for p in P:
        for r in P:
            if p == r:
                continue
            d = order(p, r)
            c = 0 if d == 1 else valuation(p ** d - 1, r)
            k = r if d == 1 else d
            for a in range(1, 101):
                expected = 0 if (a + 1) % k else (0 if d == 1 else c) + valuation(a + 1, r)
                require(valuation(sigma(p, a), r) == expected)
                counts['lte'] += 1
            for euler in (False, True):
                q = p if euler else -1
                E = graph((p, r), q)
                for U in (2, 5, 10, 15, 31, 64, 100):
                    cap = capacity(E[p, r], U, r) if (p, r) in E else 0
                    for a in range(1, U + 1):
                        if euler and a % 4 != 1 or (not euler and a % 2):
                            continue
                        require(valuation(sigma(p, a), r) <= cap)
                        counts['capacities'] += 1
    for p in P[:8]:
        for euler in (False, True):
            q = p if euler else -1
            E = graph(P[:8], q)
            out = {r: z for (t, r), z in E.items() if t == p}
            b = 2 if euler else 1
            for a in range(1, 101):
                if euler and a % 4 != 1 or (not euler and a % 2):
                    continue
                targets = {r for r in out if sigma(p, a) % r == 0}
                L = lcm(b, *(out[r]['k'] for r in targets))
                require(L <= a + 1 and (a + 1) % L == 0)
                require(targets == {r for r, z in out.items() if L % z['k'] == 0})
                counts['lcm_patterns'] += 1
    S = (3, 5, 7, 11, 13)
    certs = [certificate(S, q) for q in (5, 13)]
    require(certs[0]['C'] == 3 and certs[0]['B'] == 16)
    require(certs[0]['rounds'][0]['raw_column_caps'][13] == 1)
    require(certs[1]['C'] == 4 and certs[1]['B'] == 32)
    require(certs[1]['rounds'][0]['after_role_rounding'] == {3: 4, 5: 2, 7: 4, 11: 4, 13: 1})
    require(not any((p == 5 for p, r in certs[1]['remaining_edges'])))
    require(prod((Fraction(p, p - 1) for p in S)) == Fraction(1001, 384) > 2)
    require(sigma(3, 4) % 11 == 0 and sigma(3, 2) % 13 == 0)
    require(lcm(order(3, 11), order(3, 13)) == 15 > 11)
    require(order(3, 4561) == 15 and all((4561 % d for d in range(2, 68))))
    counts['cut_bookkeeping'] = 0
    for k in (2, 3, 5, 7):
        p = 2 ** k - 1
        a = {2: k - 1, p: 1}
        support = (2, p)
        w = {(u, r): valuation(sigma(u, a[u]), r) for u in support for r in support}
        w[p, 2] -= 1
        require(all((sum((w[u, r] for u in support)) == a[r] for r in support)))
        for t in support:
            T = {t}
            outside = set(support) - T
            O = prod((r ** sum((w[u, r] for u in T)) for r in outside))
            I = prod((r ** sum((w[u, r] for u in outside)) for r in T))
            ratio = Fraction(sigma(t, a[t]), t ** a[t])
            require(ratio == 2 ** int(p in T) * Fraction(O, I))
            require(1 < ratio < 2)
            require(O < I < 2 * O if p in T else I < O < 2 * I)
            counts['cut_bookkeeping'] += 1
    out = {'checks': counts, 'support': S, 'abundancy_ceiling': '1001/384', 'certificates': certs, 'status': 'all exact checks passed'}
    Path(__file__).with_name('RESULTS.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    run()

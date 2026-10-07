#!/usr/bin/env python3
"""A golden analogue of THM-4581's pair chain.

Map on O_2 = Z_2[phi] (2 is inert; residue field F_4 = {0, 1, phi, phi^2}):
    G(x) = x/2                 if x = 0 mod 2
    G(x) = (phi x + rho(x))/2  otherwise, rho(x) = the {0,1}-coefficient representative of phi*x mod 2.
G is 4-to-1 and Haar-preserving; residues of G^n(v) are i.i.d. uniform on F_4 for Haar v.
Relation u = phi^k v + e (phi a unit, so e stays in O).  The chain (k, e) is driven by r = v mod 2.
Checks: (1) the chain table against direct orbits of integral u, v (exact, many steps);
        (2) BFS: every small state (0, e) reaches (0, 0);
        (3) Monte Carlo: merge probability and tail of the merge time, from (0, 1).
"""
import random, sys

def fib(k):
    a, b = 0, 1
    for _ in range(k): a, b = b, a + b
    return a

def phipow(k):
    if k >= 0:
        return (1, 0) if k == 0 else (fib(k - 1), fib(k))
    m = -k; s = -1 if m % 2 else 1
    return (s * fib(m + 1), -s * fib(m))

def mul(x, y):
    a, b = x; c, d = y
    return (a * c + b * d, a * d + b * c + b * d)
def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def half(x):
    assert x[0] % 2 == 0 and x[1] % 2 == 0, x
    return (x[0] // 2, x[1] // 2)
def res(x): return (x[0] % 2, x[1] % 2)
PHI = (0, 1)

def G(x):
    if res(x) == (0, 0):
        return half(x)
    rho = res(mul(PHI, x))
    return half(add(mul(PHI, x), rho))

def chain_step(k, e, r):
    """r = v mod 2 in F_4 (as a 0/1 pair); returns new (k, e)."""
    pk = phipow(k)
    ures = res(add(mul(pk, r), e))
    if r == (0, 0) and ures == (0, 0):
        return k, half(e)
    rho_v = res(mul(PHI, r))
    if r != (0, 0) and ures != (0, 0):
        rho_u = res(mul(PHI, ures))
        return k, half(sub(add(mul(PHI, e), rho_u), mul(pk, rho_v)))
    if r == (0, 0):            # v even, u odd: k+1
        rho_u = res(mul(PHI, ures))
        return k + 1, half(add(mul(PHI, e), rho_u))
    # v odd, u even: k-1
    return k - 1, half(sub(e, mul(phipow(k - 1), rho_v)))

def check_table(trials=300, steps=200):
    random.seed(5)
    for _ in range(trials):
        v = (random.randrange(-10 ** 30, 10 ** 30), random.randrange(-10 ** 30, 10 ** 30))
        k = random.randrange(-4, 5)
        e = (random.randrange(-50, 50), random.randrange(-50, 50))
        u = add(mul(phipow(k), v), e)
        for _ in range(steps):
            r = res(v)
            k, e = chain_step(k, e, r)
            u, v = G(u), G(v)
            assert u == add(mul(phipow(k), v), e)
    print(f'(1) chain table matches direct orbits: {trials} pairs x {steps} steps, 0 mismatches')

def bfs_absorb(E=12, depth=60, K=20):
    from collections import deque
    R4 = [(0, 0), (1, 0), (0, 1), (1, 1)]
    bad = []
    for a in range(-E, E + 1):
        for b in range(-E, E + 1):
            if (a, b) == (0, 0): continue
            start = (0, (a, b)); seen = {start}; dq = deque([(start, 0)]); ok = False
            while dq and not ok:
                (k, e), d = dq.popleft()
                if d >= depth: continue
                for r in R4:
                    s = chain_step(k, e, r)
                    if s == (0, (0, 0)):
                        ok = True; break
                    if abs(s[0]) <= K and max(abs(s[1][0]), abs(s[1][1])) <= 50 * E and s not in seen:
                        seen.add(s); dq.append((s, d + 1))
            if not ok: bad.append((a, b))
    print(f'(2) BFS: states (0,e), |a|,|b| <= {E}: {(2*E+1)**2-1 - len(bad)} reach (0,0); failures: {bad[:10]}')

def monte_carlo(paths=4000, T=40000, start=(0, (1, 0))):
    random.seed(7)
    R4 = [(0, 0), (1, 0), (0, 1), (1, 1)]
    times = []
    for _ in range(paths):
        k, e = start; t = 0
        while t < T:
            k, e = chain_step(k, e, R4[random.randrange(4)]); t += 1
            if k == 0 and e == (0, 0): break
            if abs(k) > 400: t = T; break
        times.append(t if (k, e) == (0, (0, 0)) else None)
    merged = [t for t in times if t is not None]
    print(f'(3) Monte Carlo from (k,e)={start}: merged {len(merged)}/{paths} by T={T}; median merge time '
          f'{sorted(merged)[len(merged)//2] if merged else None}')
    for TT in (100, 400, 1600, 6400, 25600):
        q = sum(1 for t in times if t is None or t > TT) / paths
        print(f'     T={TT:6d}  P(no merge by T)={q:.4f}  sqrt(T)*q={q * TT ** 0.5:.2f}')

if __name__ == '__main__':
    check_table()
    bfs_absorb()
    monte_carlo(int(sys.argv[1]) if len(sys.argv) > 1 else 3000)

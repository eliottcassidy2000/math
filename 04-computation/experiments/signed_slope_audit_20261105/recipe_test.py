#!/usr/bin/env python3
"""Independent test of THM-4882 (3)(ii): both-odd slopes, the cyclic-corrected clock.

Own implementation of Session 198's O/Z calculus, the clock word, the constants R, C_lin, C_cyc, the periodic-orbit
computation of K_cyc (cyclic run decomposition), the transfer formula with the boundary term -eps_0 E, and the recipe
seed -> O -> W^h for chosen targets.  The final point is cross-checked against the orbit (D must equal g(U)), and
small final words are written out for an exact subset-DP count (hamcount.exe).
"""
import sys, random
from fractions import Fraction
from math import gcd
sys.setrecursionlimit(100000)

def g(n):
    # iterative to handle huge n: unwind the recurrence
    acc = 0
    while n > 1:
        if n % 2 == 0: n //= 2
        else:
            m = (n + 1) // 2
            acc += m - (m % 2)
            n = m
    return acc + 1

def O(p): D, E, U, d = p; return (D, 0, 2*U, d)
def Z(p): D, E, U, d = p; return (D - E + U, d, 2*U - d, d)
def word_apply(p, w):
    for c in w: p = O(p) if c == 'O' else Z(p)
    return p

def clock(a, b):
    k = 1
    while pow(2, k, b) != 1: k += 1
    B = 1 << k; A = a * (B - 1) // b; assert A * b == a * (B - 1)
    eps = [(A >> i) & 1 for i in range(k)]
    W = ''.join('Z' if e else 'O' for e in eps)
    R = sum(eps[i] << (k - 1 - i) for i in range(k))
    C_lin = sum(eps[j] * eps[i] << (i - j - 1) for i in range(k) for j in range(i)) + sum(eps[i-1] * eps[i] for i in range(1, k))
    C_cyc = C_lin + eps[k-1] * eps[0]
    return dict(k=k, B=B, A=A, eps=eps, W=W, R=R, C_lin=C_lin, C_cyc=C_cyc, K_lin=a*R - b*C_lin, K_cyc=a*R - b*C_cyc)

def periodic_K(a, b, c):
    """K_cyc / b via the periodic orbit u_0 = R/(B-1) with cyclic ending counts, and the run decomposition."""
    eps = c['eps']; k = c['k']
    u = Fraction(c['R'], c['B'] - 1); us = []; terms = []
    for i in range(k):
        assert 0 < u < 1 and eps[i] == (1 if 2*u >= 1 else 0)
        us.append(u); terms.append(eps[i] * (u - eps[i-1]))   # eps[-1] = eps[k-1]: cyclic
        u = 2*u - eps[i]
    assert u == us[0]
    delta = sum(terms)
    # run decomposition over the cyclic word
    runs = []
    for i in range(k):
        if eps[i] and not eps[i-1]:
            j = 1
            while eps[(i+j) % k]: j += 1
            v = us[(i+j) % k]
            part = v + (1 - v) / 2**j
            assert part == sum(terms[(i+h) % k] for h in range(j)), "run contribution mismatch"
            runs.append(part)
    assert sum(runs) == delta
    return delta, runs

def main():
    rng = random.Random(7)
    slopes = [(3,5), (5,7), (7,11), (5,9), (7,9), (9,11), (7,13), (9,13), (3,5), (11,13), (9,17), (13,15), (11,19), (15,17)]
    seen = set()
    for (a, b) in slopes:
        if (a, b) in seen: continue
        seen.add((a, b))
        assert a % 2 == 1 and b % 2 == 1 and gcd(a, b) == 1 and b < 2*a < 2*b
        c = clock(a, b)
        assert c['A'] % 2 == 1 and c['W'][0] == 'Z' and c['W'][-1] == 'Z' and 'O' in c['W'] and c['B']//2 <= c['A'] < c['B'] - 1
        delta, runs = periodic_K(a, b, c)
        assert delta == Fraction(c['K_cyc'], b), (delta, c)
        assert c['K_cyc'] > 0 and c['K_cyc'] % 2 == 0 and c['K_lin'] == c['K_cyc'] + b
        # transfer with boundary term on random inputs (d = 1, E in {0,1}) and on symbolic basis vectors
        for _ in range(100):
            D = rng.randrange(1, 10**6); U = rng.randrange(1, 10**6); E = rng.randrange(0, 2)
            p = (D, E, U, 1); q = word_apply(p, c['W'])
            assert q == (D + c['A']*U - c['C_lin'] - c['eps'][0]*E, 1, c['B']*U - c['R'], 1), "transfer formula"
            assert a*q[2] - b*q[0] == a*U - b*D - c['K_lin'] + b*c['eps'][0]*E, "response drop"
        # recipe on targets
        targets = [0, 2, -2, 4, -4, 6, -6, 100, -100, -1000, -10000, 20000, -30000] + [rng.randrange(-10**4, 10**4) * 2 for _ in range(10)]
        results = []
        for r in targets:
            # residue condition: 2 a U0 - b D0 = r + b (mod K_cyc) on a strong orbit point (odd n0, D0 = g(n0)), with S - b - r >= K_cyc
            K = c['K_cyc']; found = None
            for n0 in range(3, 1 << 22, 2):
                D0 = g(n0); S = 2*a*n0 - b*D0
                if (2*a*n0 - b*D0 - r - b) % K == 0 and S - b - r >= K:
                    found = n0; break
            assert found is not None, (a, b, r)
            n0 = found; D0 = g(n0); S = 2*a*n0 - b*D0; h = (S - b - r) // K
            assert h >= 1
            p = O((D0, 1, n0, 1))
            p = word_apply(p, c['W'])            # first period: drop K_lin
            assert a*p[2] - b*p[0] == S - c['K_lin']
            for _ in range(h - 1): p = word_apply(p, c['W'])
            D, E, U, d = p
            assert E == 1 and d == 1 and U % 2 == 1 and D % 2 == 1
            assert a*U - b*D == r, (a, b, r, a*U - b*D)
            assert g(U) == D, "final point is not on the orbit"   # independent cross-check of the whole transfer
            results.append((r, n0, h, U.bit_length()))
        print(f"slope {a}/{b}: k={c['k']} A={c['A']} W={c['W']} R={c['R']} C_lin={c['C_lin']} C_cyc={c['C_cyc']} K_lin={c['K_lin']} K_cyc={c['K_cyc']} runs={[str(x) for x in runs]} | {len(results)} targets hit exactly; (target, seed n0, h, bits of final U): {results[:7]} ...")
    # Write small final words for exact DP counting: small seeds n0 and 1..3 periods; expected r = S - b - h K_cyc
    with open('small_words.txt', 'w') as f:
        cnt = 0
        for (a, b) in [(3,5), (5,7), (7,11), (5,9), (7,9)]:
            c = clock(a, b); K = c['K_cyc']
            for n0 in [3, 5, 7, 9, 11, 13, 15]:
                S = 2*a*n0 - b*g(n0)
                for h in [1, 2, 3]:
                    p = word_apply(O((g(n0), 1, n0, 1)), c['W'] * h)
                    r = S - b - h*K
                    assert a*p[2] - b*p[0] == r
                    if (p[2] - 1).bit_length() <= 18:
                        f.write(" ".join(str(x) for x in (a, b, r, n0, h, p[0], p[1], p[2], p[3], bin(p[2]-1)[2:])) + chr(10)); cnt += 1
    print(f"{cnt} small words written to small_words.txt (a b r n0 h D E U d word)")

main()

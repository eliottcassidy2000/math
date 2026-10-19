"""Goldbach / twin-prime singular series as a sum over Farey fractions, truncated by denominator Q versus by Stern-Brocot depth d.
S(N) = sum_{q>=1} mu(q)^2/phi(q)^2 * c_q(N),  c_q(N) = Ramanujan sum = sum over a/q in lowest terms of e(aN/q)  (Hardy-Littlewood).
Twin primes: N = 2 (singular series 2 C_2 = 1.3203...). Goldbach: even N. We compare the two truncation orders."""
from math import gcd, cos, pi, log
from sympy import totient, mobius, factorint
import sys
def ramanujan(q, N):
    # c_q(N) = sum_{d | gcd(q,N)} d mu(q/d)
    g = gcd(q, N); s = 0
    for d in range(1, g+1):
        if g % d == 0: s += d * int(mobius(q // d))
    return s
def term(q, N):
    m = int(mobius(q))
    if m == 0: return 0.0
    return ramanujan(q, N) / float(totient(q))**2
def exact_series(N):
    # Euler product: 2 C2 prod_{p | N, p>2} (p-1)/(p-2) for even N, 0 for odd N
    if N % 2: return 0.0
    C2 = 1.0; p = 3
    from sympy import primerange
    for p in primerange(3, 200000): C2 *= 1 - 1/(p-1)**2
    S = 2*C2
    for p in factorint(N):
        if p > 2: S *= (p-1)/(p-2)
    return S
def sb_depth_fractions(dmax):
    """Farey fractions a/q in (0,1) by Stern-Brocot depth; depth of a/q = (sum of CF partial quotients) - 1 (root 1/2 at depth 1)"""
    out = {}
    frontier = [((0,1),(1,1))]
    for d in range(1, dmax+1):
        nf = []
        for (a,b),(c,e) in frontier:
            m = (a+c, b+e); out.setdefault(d, []).append(m)
            nf.append(((a,b), m)); nf.append((m, (c,e)))
        frontier = nf
    return out
for N in (2, 4, 30, 210, 1000, 2310):
    S = exact_series(N)
    # denominator truncation
    partial_q = {}; s = 0.0
    for q in range(1, 2001):
        s += term(q, N)
        if q in (10, 30, 100, 300, 1000, 2000): partial_q[q] = s
    # Stern-Brocot depth truncation: include a/q when its depth <= d (the term for q is split over its phi(q) fractions a/q: each contributes mu(q)^2/phi(q)^2 * e(aN/q))
    fr = sb_depth_fractions(16)
    partial_d = {}; s = term(1, N)   # q = 1 (the fraction 0/1 ~ 1/1) first
    for d in range(1, 17):
        for (a, q) in fr[d]:
            m = int(mobius(q))
            if m == 0: continue
            s += cos(2*pi*a*N/q) / float(totient(q))**2
        if d in (4, 8, 12, 16): partial_d[d] = s
    print(f"N={N}: exact S(N) = {S:.6f}; by denominator Q: " + ", ".join(f"Q={q}: {v:.5f}" for q, v in partial_q.items()))
    print(f"        by Stern-Brocot depth d (max denominator F_(d+2)): " + ", ".join(f"d={d}: {v:.5f}" for d, v in partial_d.items()))

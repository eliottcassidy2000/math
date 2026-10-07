#!/usr/bin/env python3
"""Audit A, item 6(d)/(f): S19's lag-D Mersenne switch vs the pair chain on actual integers.
p = 2*3^(a-1) - 1 (post-run of M_a), q_D = 2*3^(a-1-D) - 1 (post-run of M_(a-D)); p = 3^D q_D + 3^D - 1.
Chain from (D, 3^D - 1) driven by q_D's Terras parities. Compare: absorption time t, S19's switch (first i with
U^i(p) == U^(i+D)(q_D)) and its template total (T-time of p at U^i(p)), sigma(M_a) == sigma(M_(a-D))."""
import sys
from a3_drift import step, v2
def U(x):
    y = 3 * x + 1
    return y >> v2(y)
def sigma(n):
    s = 0
    while n != 1:
        n = U(n); s += 1
    return s
def T(x): return x >> 1 if x % 2 == 0 else (3 * x + 1) >> 1

def chain_absorb(p, q, D, tmax):
    k, N = D, 3 ** D - 1
    assert p == 3 ** D * q + N
    u, v = p, q
    for t in range(tmax):
        if (k, N) == (0, 0):
            assert u == v
            return t, u
        if u == 1 or v == 1:   # integer orbits reached 1 before absorption: stop (cycle region)
            return None, None
        beta = v & 1
        k, N = step(k, N, beta)
        u, v = T(u), T(v)
        assert u == (3 ** k * v + N if k >= 0 else None) or k < 0
    return None, None

def s19_switch(p, q, D):
    """first i with U^i(p) == U^(i+D)(q), stopped at 1; returns (i, template total of p at index i, value)"""
    P = [p]; Ptot = [0]
    x = p; tot = 0
    while x != 1:
        y = 3 * x + 1; e = v2(y); x = y >> e; tot += e; P.append(x); Ptot.append(tot)
    Q = [q]; x = q
    while x != 1:
        x = U(x); Q.append(x)
    for i in range(len(P)):
        j = i + D
        if j >= len(Q): return None
        if P[i] == Q[j] and P[i] != 1:
            return i, Ptot[i], P[i]
    return None

stats = {'absorbed': 0, 'switch': 0, 'abs_and_switch': 0, 'abs_no_switch': 0, 'switch_no_abs': 0,
         'tt_eq_t+v2': 0, 'tt_ne': 0, 'sigma_eq_given_abs': 0, 'sigma_ne_given_abs': 0, 'q_mod16': set(), 'odd_counts_ok': 0}
D = 1
A = int(sys.argv[1]) if len(sys.argv) > 1 else 1201
for a in range(5, A + 1, 2):
    p = 2 * 3 ** (a - 1) - 1; q = 2 * 3 ** (a - 1 - D) - 1
    stats['q_mod16'].add(q % 16)
    t, z = chain_absorb(p, q, D, 10 ** 6)
    sw = s19_switch(p, q, D)
    if t is not None: stats['absorbed'] += 1
    if sw is not None: stats['switch'] += 1
    if t is not None and sw is not None:
        stats['abs_and_switch'] += 1
        i, tt, val = sw
        if tt == t + v2(z): stats['tt_eq_t+v2'] += 1
        else: stats['tt_ne'] += 1
    elif t is not None: stats['abs_no_switch'] += 1
    elif sw is not None: stats['switch_no_abs'] += 1
    if t is not None:
        if sigma(2 ** a - 1) == sigma(2 ** (a - 1) - 1): stats['sigma_eq_given_abs'] += 1
        else: stats['sigma_ne_given_abs'] += 1
print(f"[mersenne lag-1, odd a in [5,{A}]] {stats}")

#!/usr/bin/env python3
"""Imaginary quadratic class numbers as ranking excesses of Kronecker tournaments (S40).

For a negative fundamental discriminant D with Kronecker character chi_D, let Gamma_D be the oriented graph on Z/|D|
with an arc i -> j whenever chi_D(j - i) = +1 (chi_D(-1) = -1, so each pair with gcd(j - i, D) = 1 carries exactly one
arc). For D = -p, p = 3 (mod 4) prime, Gamma_D is the Paley tournament P_p. Ordering the vertices 0 < 1 < ... < |D|-1,

    forward - backward = sum_(d=1)^(|D|-1) (|D| - d) chi_D(d) = -sum_d d chi_D(d) = (2/w_D) |D| h(D)

(Dirichlet's class number formula h(D) = -(w_D / (2|D|)) sum_(a=1)^(|D|) a chi_D(a); w_D = 6, 4, 2 for D = -3, -4,
D < -4). This script checks the identity against class numbers counted from reduced forms, for every prime
p = 3 (mod 4) below PMAX and every negative fundamental discriminant down to -DMAX; checks Dirichlet's half-window form
(out-degree minus in-degree of vertex 0 on {1, ..., (p-1)/2} equals (2 - (2/p)) h(-p) for p > 3); and computes the exact
maximum acyclic subgraph MAS(P_p) for small p by subset dynamic programming (vertex 0 first, by vertex-transitivity),
comparing it with the natural order.
Usage: python paley_class_number_ranking_20261017.py [PMAX] [DMAX] [MASMAX]
Session opus-2026-10-17-S40.
"""
import sys, json, os
import numpy as np
from numba import njit


@njit(cache=True)
def class_number(D):
    """Number of reduced primitive positive definite forms of discriminant D < 0."""
    h = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a) != 0:
                continue
            c = (b * b - D) // (4 * a)
            if c < a:
                continue
            if b < 0 and (a == c):
                continue
            g = a
            x, y = abs(b), c
            while x:
                g, x = x, g % x
            while y:
                g, y = y, g % y
            if g != 1:
                continue
            h += 1
        a += 1
    return h


@njit(cache=True)
def kronecker_table(D):
    """chi_D(a) for 0 <= a < |D| (fundamental D < 0), via multiplicativity over the prime factors of D."""
    n = -D
    chi = np.zeros(n, dtype=np.int8)
    for a in range(n):
        chi[a] = kron(D, a)
    return chi


@njit(cache=True)
def kron(D, a):
    """Kronecker symbol (D / a) for a >= 0 (D a discriminant)."""
    if a == 0:
        return 1 if (D == 1 or D == -1) else 0
    res = 1
    # factor out 2
    while a % 2 == 0:
        a //= 2
        if D % 2 == 0:
            return 0
        r = D % 8
        if r == 3 or r == 5:
            res = -res
    # now a odd: Jacobi symbol (D / a)
    if a == 1:
        return res
    x = D % a
    y = a
    t = 1
    while x != 0:
        while x % 2 == 0:
            x //= 2
            r = y % 8
            if r == 3 or r == 5:
                t = -t
        x, y = y, x
        if x % 4 == 3 and y % 4 == 3:
            t = -t
        x = x % y
    if y == 1:
        return res * t
    return 0


@njit(cache=True)
def excess(chi):
    n = chi.shape[0]
    s = 0
    for d in range(1, n):
        s += (n - d) * chi[d]
    return s


@njit(cache=True)
def mas_dp(adj, n):
    """Maximum number of forward arcs over all linear orders of a vertex-transitive tournament with vertex 0 first.
    adj[u] = bitmask of out-neighbours of u. DP over subsets of {1, ..., n-1} placed after 0."""
    m = n - 1
    size = 1 << m
    dp = np.full(size, -1, dtype=np.int16)
    dp[0] = 0
    # inbits[v]: in-neighbours of v among 1..n-1 as a mask over the m remaining vertices (bit v-1)
    inmask = np.zeros(n, dtype=np.int64)
    for v in range(1, n):
        msk = 0
        for u in range(1, n):
            if (adj[u] >> v) & 1:
                msk |= 1 << (u - 1)
        inmask[v] = msk
    from0 = np.zeros(n, dtype=np.int64)
    for v in range(1, n):
        from0[v] = (adj[0] >> v) & 1
    for S in range(size):
        base = dp[S]
        if base < 0:
            continue
        for v in range(1, n):
            bit = 1 << (v - 1)
            if S & bit:
                continue
            x = S & inmask[v]
            c = 0
            while x:
                x &= x - 1
                c += 1
            val = base + c + from0[v]
            T = S | bit
            if val > dp[T]:
                dp[T] = val
    return dp[size - 1]


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def is_fundamental(D):
    n = -D
    if D % 4 == 1 or (D % 4 + 4) % 4 == 1:
        r = n
        k = 2
        while k * k <= r:
            if r % (k * k) == 0:
                return False
            k += 1
        return True
    if D % 4 == 0:
        m = D // 4
        if m % 4 in (2, 3) or (m % 4 + 4) % 4 in (2, 3):
            r = abs(m)
            k = 2
            while k * k <= r:
                if r % (k * k) == 0:
                    return False
                k += 1
            return True
    return False


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    DMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    MASMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 23
    out = {}
    # 1. Paley tournaments
    bad = 0; count = 0; half_bad = 0; one = []
    for p in range(3, PMAX + 1, 4):
        if not is_prime(p):
            continue
        chi = kronecker_table(-p)
        h = class_number(-p)
        w = 6 if p == 3 else 2
        F_minus_B = excess(chi)
        if F_minus_B * w != 2 * p * h:
            bad += 1
        if p > 3:
            half = int(np.sum(chi[1:(p - 1) // 2 + 1]))
            two = 1 if p % 8 == 7 else -1
            if half != (2 - two) * h:
                half_bad += 1
        if h == 1:
            one.append(p)
        count += 1
    print('Paley tournaments P_p, p = 3 (mod 4) prime <= %d: %d primes; forward - backward = p h(-p) (2/w) fails at %d; '
          'half-window imbalance = (2 - (2/p)) h(-p) fails at %d' % (PMAX, count, bad, half_bad))
    print('primes with h(-p) = 1 (natural order has the least possible excess p, w-corrected at 3):', one)
    out['paley'] = dict(PMAX=PMAX, primes=count, identity_failures=bad, half_window_failures=half_bad, class_number_one=one)
    # 2. all negative fundamental discriminants
    bad = 0; count = 0
    for n in range(3, DMAX + 1):
        D = -n
        if not is_fundamental(D):
            continue
        chi = kronecker_table(D)
        h = class_number(D)
        w = 6 if D == -3 else (4 if D == -4 else 2)
        if excess(chi) * w != 2 * n * h:
            bad += 1
        count += 1
    print('Kronecker oriented graphs, fundamental D in [-%d, -3]: %d discriminants; identity fails at %d' % (DMAX, count, bad))
    out['kronecker'] = dict(DMAX=DMAX, discriminants=count, identity_failures=bad)
    # 3. exact MAS for small Paley tournaments
    rows = []
    for p in range(3, MASMAX + 1, 4):
        if not is_prime(p):
            continue
        chi = kronecker_table(-p)
        adj = np.zeros(p, dtype=np.int64)
        for u in range(p):
            m = 0
            for v in range(p):
                if v != u and chi[(v - u) % p] == 1:
                    m |= 1 << v
            adj[u] = m
        arcs = p * (p - 1) // 2
        h = class_number(-p)
        w = 6 if p == 3 else 2
        natural = (arcs + 2 * p * h // w) // 2
        mas = int(mas_dp(adj, p))
        rows.append(dict(p=p, arcs=arcs, h=h, natural_forward=natural, MAS=mas, natural_excess=2 * natural - arcs,
                         optimal_excess=2 * mas - arcs))
        print('P_%d: arcs %d, h(-p) = %d, natural order forward %d (excess %d), MAS %d (excess %d)'
              % (p, arcs, h, natural, 2 * natural - arcs, mas, 2 * mas - arcs), flush=True)
    print('P_43 (cited, Chindelevitch-Harutyunyan arXiv:2607.26690 Thm 7.1): MAS 543 of 903; natural order forward %d'
          % ((903 + 43 * class_number(-43)) // 2))
    out['mas'] = rows
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'paley_class_number_ranking_20261017.json'), 'w'),
              indent=1)


if __name__ == '__main__':
    main()

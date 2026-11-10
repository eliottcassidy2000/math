#!/usr/bin/env python3
"""AUDIT: identities that tie Theorems A and B to canon the packet does not cite (auditor's own code).
 (i)   THM-4757 (remarks after its section 3; audited 2026-10-22):  sum_L (n-L)! N_L t^L = (1-t)^(n-1) F(T,(1+t)/(1-t))
       and [r^(n-1-L)] W(r) = 2^(-L) (n-L)! N_L, with N_L the signed count of directed linear forests with L arcs,
       and THM-4761 (1): kappa_L = N_L/(n)_L is the Walsh chaos of degree L of kappa = H 2^(n-1)/n!, for every n.
       Together:  2^(n-1) (N_rho H)(T) = sum_L (n-L)! N_L rho^L  -- Theorem A (1)-(2) for all n.
       Checked here by brute force: N_L from linear forests (via ordered path covers), Walsh levels from the cube.
 (ii)  THM-059: E_p[H_n] = F_(n-1)(r)/(2r)^(n-1) with F the master polynomial, and for n = 2k+1
       sum_j gamma_(n,j) w^j = sum_j b_(k,j) (2j+1)! w^j (1-4w)^(k-j)   (central factorial triangle b).
 (iii) THM-4756: F(T,x) = H(T) A_n(x) (mod 3)  ==>  gamma_j(T) = H(T) gamma_(n,j) (mod 3) for every j.
"""
from fractions import Fraction as Fr
from itertools import permutations
from math import comb, factorial
import random
def pairs(n): return [(i, j) for i in range(n) for j in range(i + 1, n)]
def adj_from(n, code):
    a = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs(n)):
        if code >> k & 1: a[i][j] = 1
        else: a[j][i] = 1
    return a
def profile(n, adj):
    a = [0] * n
    for P in permutations(range(n)): a[sum(adj[P[t]][P[t + 1]] for t in range(n - 1))] += 1
    return a
def gamma(a):
    d = len(a) - 1; c = list(a); g = []
    for j in range(d // 2 + 1):
        g.append(c[j])
        for i in range(d - 2 * j + 1): c[j + i] -= g[j] * comb(d - 2 * j, i)
    assert not any(c); return g
def eulerian(n): return [sum((-1) ** i * comb(n + 1, i) * (k + 1 - i) ** n for i in range(k + 2)) for k in range(n)]
def padd(A, B):
    m = max(len(A), len(B)); return [(A[i] if i < len(A) else 0) + (B[i] if i < len(B) else 0) for i in range(m)]
def pmul(A, B):
    C = [0] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        for j, y in enumerate(B): C[i + j] += x * y
    return C
def ppow(A, e):
    R = [1]
    for _ in range(e): R = pmul(R, A)
    return R
ok = True
# (i) N_L by brute force: N_L = sum over directed linear forests with L arcs of prod S_e, S = A - A^T.
#     A directed linear forest with L arcs on n vertices = set of n-L directed paths covering V; count via permutations:
#     each forest corresponds to (n-L)! orderings of its paths, i.e. to pairs (permutation, set of L "kept" steps).
print("== (i) THM-4757 / THM-4761 form of Theorem A (all labelled tournaments, n = 3, 4, 5) ==")
for n in (3, 4, 5):
    m = n * (n - 1) // 2; N = 1 << m; d = n - 1; bad = 0
    Hall = []
    for code in range(N):
        adj = adj_from(n, code); Hall.append(sum(1 for P in permutations(range(n)) if all(adj[P[t]][P[t + 1]] for t in range(d))))
    hat = Hall[:]; h = 1
    while h < N:
        for i in range(0, N, 2 * h):
            for j in range(i, i + h): x, y = hat[j], hat[j + h]; hat[j], hat[j + h] = x + y, x - y
        h *= 2
    for code in range(N):
        adj = adj_from(n, code); S = [[adj[i][j] - adj[j][i] for j in range(n)] for i in range(n)]
        NL = [Fr(0)] * n                      # sum over (perm, kept step set K) of prod_{t in K} S / (n-L)!
        for P in permutations(range(n)):
            e = [1]                           # elementary symmetric polynomials of the step signs
            for t in range(d):
                s = S[P[t]][P[t + 1]]; e = [(e[k] if k < len(e) else 0) + (s * e[k - 1] if k >= 1 else 0) for k in range(len(e) + 1)]
            for L in range(n): NL[L] += e[L]
        NL = [NL[L] / factorial(n - L) for L in range(n)]
        lev = [Fr(0)] * (m + 1)
        for Sx in range(N):
            if hat[Sx]: lev[bin(Sx).count("1")] += Fr(hat[Sx] * (-1) ** bin(Sx & code).count("1"), N)
        a = profile(n, adj)
        for L in range(n):
            if lev[L] * 2 ** d != factorial(n - L) * NL[L]: bad += 1          # THM-4761: Walsh level L of H = n! kappa_L / 2^(n-1)
        for rho in (Fr(1, 3), Fr(-2, 5)):
            p = (1 + rho) / 2; q = (1 - rho) / 2
            if 2 ** d * sum(a[k] * p ** k * q ** (d - k) for k in range(n)) != sum(factorial(n - L) * NL[L] * rho ** L for L in range(n)): bad += 1
        W = [0]
        for k in range(n): W = padd(W, [a[k] * x for x in pmul(ppow([Fr(1, 2), 1], k), ppow([Fr(-1, 2), 1], d - k))])
        for L in range(n):
            if W[d - L] != Fr(factorial(n - L), 2 ** L) * NL[L]: bad += 1     # THM-4757: [r^(n-1-L)] W = 2^-L (n-L)! N_L
    ok &= (bad == 0)
    print(f"   n={n}: {N} tournaments, failures {bad}")
# (ii) THM-059 master polynomial in the w coordinate
print("== (ii) THM-059: sum_j gamma_(n,j) w^j = sum_j b_(k,j) (2j+1)! w^j (1-4w)^(k-j), n = 2k+1 ==")
b = {(0, 0): 1}
for k in range(1, 8):
    for j in range(k + 1): b[(k, j)] = b.get((k - 1, j - 1), 0) + (j + 1) ** 2 * b.get((k - 1, j), 0)
for k in range(1, 8):
    n = 2 * k + 1; g = gamma(eulerian(n)); rhs = [0]
    for j in range(k + 1): rhs = padd(rhs, [b[(k, j)] * factorial(2 * j + 1) * x for x in pmul(ppow([0, 1], j), ppow([1, -4], k - j))])
    good = (rhs == g); ok &= good
    print(f"   n={n}: gamma = {g}; central-factorial row {[b[(k, j)] for j in range(k + 1)]}; identity holds: {good}")
# (iii) THM-4756 consequence
print("== (iii) THM-4756 ==> gamma_j(T) = H(T) gamma_(n,j) mod 3 ==")
rng = random.Random(7)
for n in (3, 4, 5, 6, 7):
    m = n * (n - 1) // 2; gn = gamma(eulerian(n))
    codes = range(1 << m) if n <= 5 else [rng.getrandbits(m) for _ in range(40 if n == 6 else 8)]
    bad = 0; cnt = 0
    for code in codes:
        a = profile(n, adj_from(n, code)); g = gamma(a); cnt += 1
        if any((g[j] - g[0] * gn[j]) % 3 for j in range(len(g))): bad += 1
    ok &= (bad == 0)
    print(f"   n={n}: {cnt} tournaments ({'all' if n <= 5 else 'seeded sample'}), failures {bad}")
print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")

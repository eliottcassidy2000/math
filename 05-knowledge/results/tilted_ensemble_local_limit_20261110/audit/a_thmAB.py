#!/usr/bin/env python3
"""AUDIT check (a): Theorems A and B of the tilted-ensemble packet, exact arithmetic.
Written from scratch by the auditor; shares no code with the packet.

Conventions: vertices 0..n-1, arc code bit k (pairs in lexicographic order (i<j)) = 1 iff i -> j.
"""
from fractions import Fraction as Fr
from itertools import permutations
from math import factorial, comb
import random, sys

def pair_list(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def adj_from_code(n, code):
    adj = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pair_list(n)):
        if code >> k & 1: adj[i][j] = 1
        else: adj[j][i] = 1
    return adj

def forward_profile(n, adj):
    """a[k] = number of permutations P with exactly k steps P(t)->P(t+1) that are arcs (brute force over n! perms)."""
    a = [0] * n
    for P in permutations(range(n)):
        a[sum(adj[P[t]][P[t + 1]] for t in range(n - 1))] += 1
    return a

def ham_dp(n, adj):
    """Hamiltonian path count by a start-vertex subset DP (cnt[S][v] = paths visiting exactly S, starting at v)."""
    out = [sum(1 << w for w in range(n) if adj[v][w]) for v in range(n)]
    cnt = [[0] * n for _ in range(1 << n)]
    for v in range(n): cnt[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if S >> v & 1 and S != (1 << v):
                rest = S ^ (1 << v); tot = 0; ww = out[v] & rest
                while ww:
                    w = (ww & -ww).bit_length() - 1; ww &= ww - 1
                    tot += cnt[rest][w]
                cnt[S][v] = tot
    return sum(cnt[(1 << n) - 1])

def gamma_coords(a):
    """Solve sum_k a_k x^k y^(d-k) = sum_j g_j (xy)^j (x+y)^(d-2j) by triangular elimination. Returns (g, remainder)."""
    d = len(a) - 1; c = list(a); g = []
    for j in range(d // 2 + 1):
        gj = c[j]; g.append(gj)
        for i in range(d - 2 * j + 1):
            c[j + i] -= gj * comb(d - 2 * j, i)
    return g, c

def poly_mul(A, B):
    C = [0] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        for j, y in enumerate(B): C[i + j] += x * y
    return C
def poly_pow(A, e):
    R = [1]
    for _ in range(e): R = poly_mul(R, A)
    return R
def poly_add(A, B):
    n = max(len(A), len(B)); return [(A[i] if i < len(A) else 0) + (B[i] if i < len(B) else 0) for i in range(n)]
def strip(A):
    A = list(A)
    while len(A) > 1 and A[-1] == 0: A.pop()
    return A
P_ = [0, 1]; Q_ = [1, -1]            # p and q = 1 - p as polynomials in p
def poly_in_w(A):
    """write a polynomial in p that is symmetric under p -> 1-p as a polynomial in w = p - p^2 (exact); None if impossible"""
    A = strip(A); out = {}
    while A != [0]:
        deg = len(A) - 1
        if deg % 2: return None
        k = deg // 2; c = A[-1] * (-1) ** k; out[k] = c
        A = strip(poly_add(A, [-c * x for x in poly_pow([0, 1, -1], k)]))
    return [out.get(k, 0) for k in range(max(out) + 1)] if out else [0]

ok_all = True
def check(cond, msg):
    global ok_all
    if not cond:
        ok_all = False; print("   FAIL:", msg)

# ---------------------------------------------------------------- Theorem A, exhaustive n <= 5, sample at n = 6
print("== Theorem A: all labelled tournaments, n = 2..5; 12 seeded random tournaments at n = 6 ==")
rng = random.Random(987654321)
for n in range(2, 7):
    m = n * (n - 1) // 2; d = n - 1
    codes = range(1 << m) if n <= 5 else [rng.getrandbits(m) for _ in range(12)]
    # H for every labelled tournament (needed for the brute-force noise average)
    Hall = [ham_dp(n, adj_from_code(n, c)) for c in range(1 << m)]
    cnt = 0
    for code in codes:
        adj = adj_from_code(n, code); a = forward_profile(n, adj)
        check(sum(a) == factorial(n), "sum a_k != n!")
        check(a == a[::-1], "a_k not palindromic")
        check(a[d] == Hall[code], "a_{n-1} != H (perm count vs DP)")
        # (1): brute-force noise average at rational p (two values, one with rho < 0)
        bydist = [0] * (m + 1)
        for c2 in range(1 << m): bydist[bin(code ^ c2).count("1")] += Hall[c2]
        for p in (Fr(2, 3), Fr(1, 5), Fr(1, 2), Fr(1)):
            q = 1 - p
            brute = sum(bydist[t] * q ** t * p ** (m - t) for t in range(m + 1))
            formula = sum(a[k] * p ** k * q ** (d - k) for k in range(d + 1))
            check(brute == formula, f"(1) fails n={n} code={code} p={p}")
            # (3) w-expansion
            g, rem = gamma_coords(a)
            check(all(x == 0 for x in rem), "gamma expansion has a remainder")
            check(brute == sum(gj * (p * q) ** j for j, gj in enumerate(g)), "w-expansion fails")
        g, rem = gamma_coords(a)
        check(g[0] == Hall[code], "gamma_0 != H")
        check(sum(Fr(gj, 4 ** j) for j, gj in enumerate(g)) == Fr(factorial(n), 2 ** (n - 1)), "sum gamma_j 4^-j != n!/2^(n-1)")
        # (2): W_T(r) by the repository definition sum_P prod_steps (r + s), s = +1/2 on arcs, -1/2 against arcs
        for r in (Fr(3, 2), Fr(-2), Fr(1, 2), Fr(7, 5)):
            W = Fr(0)
            for k in range(d + 1): W += a[k] * (r + Fr(1, 2)) ** k * (r - Fr(1, 2)) ** (d - k)
            rho = 1 / (2 * r); p = (1 + rho) / 2; q = (1 - rho) / 2
            check(0 <= p <= 1, "p out of range")
            noise = sum(bydist[t] * q ** t * p ** (m - t) for t in range(m + 1))
            check(W == (2 * r) ** d * noise, f"(2) fails n={n} code={code} r={r}")
        cnt += 1
    print(f"   n={n}: {cnt} tournaments checked: (1) at p=2/3,1/5,1/2,1; (2) at r=3/2,-2,1/2,7/5; (3) existence, gamma_0=H, sum gamma_j/4^j = n!/2^(n-1)")

# W by direct product over permutations (independent of a_k), n = 4, all tournaments, one r
n = 4; m = 6
for code in range(1 << m):
    adj = adj_from_code(n, code); r = Fr(5, 3)
    Wdir = sum(( (r + Fr(1,2)) if adj[P[0]][P[1]] else (r - Fr(1,2)) ) * ( (r + Fr(1,2)) if adj[P[1]][P[2]] else (r - Fr(1,2)) ) * ( (r + Fr(1,2)) if adj[P[2]][P[3]] else (r - Fr(1,2)) ) for P in permutations(range(4)))
    a = forward_profile(n, adj)
    check(Wdir == sum(a[k] * (r + Fr(1, 2)) ** k * (r - Fr(1, 2)) ** (3 - k) for k in range(4)), "W direct product != a_k form")

# ---------------------------------------------------------------- Walsh remark (n <= 5)
print("== Remark after Theorem A: Walsh levels (n = 3, 4, 5, all labelled tournaments) ==")
for n in (3, 4, 5):
    m = n * (n - 1) // 2; N = 1 << m; d = n - 1
    Hall = [ham_dp(n, adj_from_code(n, c)) for c in range(N)]
    # Walsh-Hadamard transform: hat[S] = sum_c H(c) (-1)^{|S & c|}  (exact integers; true coefficient = hat/N)
    hat = Hall[:]
    h = 1
    while h < N:
        for i in range(0, N, 2 * h):
            for j in range(i, i + h):
                x, y = hat[j], hat[j + h]; hat[j], hat[j + h] = x + y, x - y
        h *= 2
    odd_levels_vanish = all(hat[S] == 0 for S in range(N) if bin(S).count("1") % 2)
    check(odd_levels_vanish, "odd Walsh levels do not vanish")
    maxlev = max(bin(S).count("1") for S in range(N) if hat[S])
    check(maxlev <= d, "Walsh degree exceeds n-1")
    bad = 0
    for code in range(N):
        lev = [Fr(0)] * (m + 1)
        for S in range(N):
            if hat[S]:
                lev[bin(S).count("1")] += Fr(hat[S] * (-1) ** bin(S & code).count("1"), N)
        a = forward_profile(n, adj_from_code(n, code)); g, _ = gamma_coords(a)
        # gamma_j = (-4)^j sum_{i>=j} C(i,j) H^{=2i}
        for j in range(len(g)):
            if g[j] != (-4) ** j * sum(comb(i, j) * lev[2 * i] for i in range(j, m // 2 + 1)): bad += 1
        # coefficient of r^(d-2j) in W_T(r) equals 2^(d-2j) H^{=2j}
        W = [0]
        for k in range(d + 1):
            W = poly_add(W, [a[k] * x for x in poly_mul(poly_pow([Fr(1, 2), 1], k), poly_pow([Fr(-1, 2), 1], d - k))])
        for j in range(d // 2 + 1):
            if W[d - 2 * j] != 2 ** (d - 2 * j) * lev[2 * j]: bad += 1
        for e in range(d + 1):
            if (d - e) % 2 == 1 and W[e] != 0: bad += 1
    check(bad == 0, f"Walsh remark fails at n={n}")
    print(f"   n={n}: odd levels vanish: {odd_levels_vanish}; max level {maxlev}; gamma/Walsh and W-coefficient identities failures: {bad}")

# ---------------------------------------------------------------- Theorem B, exhaustive n = 3..6
print("== Theorem B: exhaustive over all labelled tournaments, n = 3..6 (polynomial identities in p) ==")
claimed = {3: [1, 2], 4: [1, 8], 5: [1, 22, 16], 6: [1, 52, 136], 7: [1, 114, 720, 272]}
def eulerian_bruteforce(n):
    A = [0] * n
    for P in permutations(range(n)):
        A[sum(1 for t in range(n - 1) if P[t] > P[t + 1])] += 1
    return A
second_moment = {}
for n in range(3, 7):
    m = n * (n - 1) // 2; d = n - 1
    S1 = [0] * (m + 1); S2 = [0] * (m + 1); S0 = [0] * (m + 1)
    for code in range(1 << m):
        H = ham_dp(n, adj_from_code(n, code)); b = m - bin(code).count("1")      # b = number of back arcs (j -> i, i<j)
        S0[b] += 1; S1[b] += H; S2[b] += H * H
    mean = [0]; sq = [0]; tot = [0]
    for b in range(m + 1):
        basis = poly_mul(poly_pow(P_, m - b), poly_pow(Q_, b))
        mean = poly_add(mean, [S1[b] * x for x in basis]); sq = poly_add(sq, [S2[b] * x for x in basis]); tot = poly_add(tot, [S0[b] * x for x in basis])
    check(strip(tot) == [1], "probabilities do not sum to 1")
    A = eulerian_bruteforce(n)
    eul = [0]
    for k in range(n): eul = poly_add(eul, [A[k] * x for x in poly_mul(poly_pow(P_, d - k), poly_pow(Q_, k))])
    g, rem = gamma_coords(A)
    gam = [0]
    for j, gj in enumerate(g): gam = poly_add(gam, [gj * x for x in poly_pow([0, 1, -1], j)])
    e1 = strip(mean) == strip(eul); e2 = strip(mean) == strip(gam)
    check(e1 and e2, f"Theorem B identity fails at n={n}")
    check(g == claimed[n], f"gamma vector differs from README at n={n}: {g}")
    check(g[1] == 2 ** n - 2 * n, "gamma_1 != 2^n - 2n")
    check(all(x == 0 for x in rem), "Eulerian gamma remainder")
    w2 = poly_in_w(sq); second_moment[n] = w2
    print(f"   n={n}: Eulerian {A}; gamma {g}; E_p[H]=Eulerian form: {e1}; =gamma form: {e2}; gamma_1=2^n-2n: {g[1] == 2**n - 2*n}; E_p[H^2] in w: {w2}")
for n in (7, 8, 9):
    g, rem = gamma_coords(eulerian_bruteforce(n))
    print(f"   n={n}: gamma vector of A_n by brute-force Eulerian numbers: {g}; gamma_1 = 2^n - 2n: {g[1] == 2**n - 2*n}; all >= 0: {all(x >= 0 for x in g)}")
    if n in claimed: check(g == claimed[n], "n=7 gamma differs from README")
readme_sq = {3: [1, 8], 4: [1, 40, 16], 5: [1, 152, 544, 56, 48], 6: [1, 552, 6288, 6944, 3072]}
for n in readme_sq: check(second_moment[n] == readme_sq[n], f"second moment polynomial differs at n={n}")
print("   second-moment polynomials n=3..6 equal README:", all(second_moment[n] == readme_sq[n] for n in readme_sq))
# THM-589 value: E[H^2] at w=1/4 equals n! W(n) / 4^(n-1), W = 8, 32, 158, 928 for n = 3..6
Wn = {3: 8, 4: 32, 5: 158, 6: 928}
for n in Wn:
    val = sum(Fr(c, 4 ** j) for j, c in enumerate(second_moment[n]))
    check(val == Fr(factorial(n) * Wn[n], 4 ** (n - 1)), f"uniform second moment != n! W(n)/4^(n-1) at n={n}")
print("   uniform second moments equal n! W(n)/4^(n-1) with W = 8, 32, 158, 928:", True if ok_all else "see failures")

# strict monotonicity on a grid + endpoints (Theorem B consequence), n = 3..40
import decimal
mono = True
for n in range(3, 41):
    # gamma by recurrence-free route: Eulerian numbers via explicit formula A(n,k) = sum_i (-1)^i C(n+1,i) (k+1-i)^n
    A = [sum((-1) ** i * comb(n + 1, i) * (k + 1 - i) ** n for i in range(k + 2)) for k in range(n)]
    g, rem = gamma_coords(A)
    if any(rem) or min(g) < 0 or g[0] != 1 or g[1] != 2 ** n - 2 * n: mono = False
    vals = [sum(gj * (Fr(t, 400) * (1 - Fr(t, 400))) ** j for j, gj in enumerate(g)) for t in range(200, 401, 10)]
    if not all(vals[i] > vals[i + 1] for i in range(len(vals) - 1)): mono = False
    if vals[0] != Fr(factorial(n), 2 ** (n - 1)) or vals[-1] != 1: mono = False
check(mono, "monotonicity / endpoints")
print("== Theorem B consequences n=3..40: gamma >= 0, gamma_0 = 1, gamma_1 = 2^n - 2n, endpoints n!/2^(n-1) and 1, strictly decreasing on a grid:", mono)
print("ALL CHECKS PASSED" if ok_all else "SOME CHECKS FAILED")

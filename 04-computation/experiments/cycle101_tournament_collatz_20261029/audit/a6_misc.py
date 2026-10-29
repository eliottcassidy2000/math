#!/usr/bin/env python3
"""Audit of the remaining checkable items:
  THM-4849 (F): Gauss-period identity mod 2^q - 1; T-trace = 2 U-trace + k on random rational cycles;
  HYP-9446 (5): 5 | D_q iff q = 3 mod 8, Paley-pair 5-integrality at q = 11, 19, 43;
  results note: H(R7) = 175, H(P7) = 189; THM-4848 (v) P7 example: the clock cubic of the 79-part of
  coker M_3 for the QR and the NQR labelling (computed from the left null space mod 79).
"""
import os, sys, random
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a3_runobs import circulant, observer
from sympy import primerange, factorint

def say(*a):
    print(" ".join(str(x) for x in a), flush=True)

def C_of(w, a=3):
    k = len(w); s = 0; t = 0
    for j in range(k):
        t += a ** (k - 1 - j) * 2 ** s
        s += w[j]
    return t

def T(x):
    return x / 2 if x.numerator % 2 == 0 else (3 * x + 1) / 2

# (F) Gauss periods
ok_all = True
qs = [q for q in primerange(7, 300) if q % 4 == 3]
for q in qs:
    QR = {(x * x) % q for x in range(1, q)}
    NQR = set(range(1, q)) - QR
    N = 2 ** q - 1
    cS = sum(2 ** j for j in QR) % N
    cN = sum(2 ** j for j in NQR) % N
    ok = (cS + cN) % N == N - 1 and (cS * cN) % N == (q + 1) // 4 % N
    ok &= (cS * cS + cS + (q + 1) // 4) % N == 0
    ok_all &= ok
say(f"(F) Gauss-period identity c_QR + c_NQR = -1, c_QR c_NQR = (q+1)/4, c^2 + c + (q+1)/4 = 0 mod 2^q - 1 for all "
    f"{len(qs)} primes q = 3 mod 4, 7 <= q < 300: {ok_all}")
say(f"    q = 7: c_QR = {sum(2 ** j for j in (1, 2, 4))}, c_NQR = {sum(2 ** j for j in (3, 5, 6))}, product mod 127 = {(22 * 104) % 127}")

# (F) T-trace = 2 U-trace + k on random words
rnd = random.Random(3)
cnt = 0
for _ in range(500):
    k = rnd.randint(1, 6)
    w = tuple(rnd.randint(1, 6) for _ in range(k))
    H = sum(w)
    D = 2 ** H - 3 ** k
    if D == 0:
        continue
    rots = [w[i:] + w[:i] for i in range(k)]
    U = sum(Fr(C_of(r), D) for r in rots)
    x = Fr(C_of(w), D)
    y = x
    tt = Fr(0)
    for _ in range(H):
        tt += y
        y = T(y)
    assert y == x
    assert tt == 2 * U + k, (w,)
    cnt += 1
say(f"(F) T-trace = 2 U-trace + k on {cnt} random rational cycles (signed denominators included): OK")

# HYP-9446 (5)
bad = [q for q in primerange(3, 2000) if ((2 ** q - 3 ** ((q - 1) // 2)) % 5 == 0) != (q % 8 == 3)]
say(f"HYP-9446 (5): 5 | D_q iff q = 3 mod 8 for primes q < 2000: {not bad}")
def paley_red(q, S):
    k = (q - 1) // 2
    s = sorted(S)
    w = tuple(((s[(i + 1) % k] - s[i]) % q) for i in range(k))
    D = 2 ** q - 3 ** k
    x = Fr(C_of(w), D)
    return x.denominator, D
for q in (11, 19, 43, 59, 67, 83):
    QR = {(x * x) % q for x in range(1, q)}
    NQR = set(range(1, q)) - QR
    dq, D = paley_red(q, QR)
    dn, _ = paley_red(q, NQR)
    say(f"   q={q}: D_q = {D}; QR cycle reduced denominator D/{D // dq}; NQR cycle D/{D // dn}")

# H(R7), H(P7)
def ham(adj, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for m in range(1, 1 << n):
        for v in range(n):
            c = dp[m][v]
            if c:
                for u in range(n):
                    if not (m >> u) & 1 and adj[v][u]:
                        dp[m | (1 << u)][u] += c
    return sum(dp[-1])
say(f"H(R7) = {ham(circulant({1, 2, 3}, 7), 7)}, H(P7) = {ham(circulant({1, 2, 4}, 7), 7)}, 2^8 - 3^4 = {2 ** 8 - 3 ** 4}")

# THM-4848 (v): clock cubic on the 79-part of coker M_3, QR and NQR labellings
def left_null_mod(M, p):
    """basis of {y : y M = 0 mod p}"""
    n = len(M)
    A = [[M[i][j] % p for i in range(n)] for j in range(n)]  # A = M^T; y M = 0 <=> M^T y^T = 0
    # row reduce A, find null space of A
    rows = [r[:] for r in A]
    piv_cols = []
    r = 0
    for c in range(n):
        pr = next((i for i in range(r, n) if rows[i][c]), None)
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        inv = pow(rows[r][c], -1, p)
        rows[r] = [x * inv % p for x in rows[r]]
        for i in range(n):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [(x - f * y) % p for x, y in zip(rows[i], rows[r])]
        piv_cols.append(c)
        r += 1
    free = [c for c in range(n) if c not in piv_cols]
    basis = []
    for fc in free:
        v = [0] * n
        v[fc] = 1
        for i, pc in enumerate(piv_cols):
            v[pc] = (-rows[i][fc]) % p
        basis.append(v)
    return basis

def clock_charpoly(S, p=79):
    n, k = 7, 3
    M, tup, idx = observer(n, circulant(S, n), k)
    Bn = left_null_mod(M, p)
    d = len(Bn)
    # translation permutation matrix R_tau with R_tau e_t = e_(tau t); then (B R_tau)[i][t] = B[i][tau t]
    def tau(t, s=1):
        return tuple((x + s) % n for x in t)
    BR = [[Bi[idx[tau(t, 1)]] for t in tup] for Bi in Bn]
    # solve BR = C B: choose d pivot columns where B is invertible
    import itertools
    from sympy import Matrix
    Bm = Matrix(Bn)
    # find independent columns
    cols = []
    rnk = 0
    for c in range(len(tup)):
        trial = cols + [c]
        if Matrix([[Bn[i][j] for j in trial] for i in range(d)]).rank(iszerofunc=lambda x: x % p == 0) > rnk:
            cols = trial
            rnk += 1
            if rnk == d:
                break
    Bsub = Matrix([[Bn[i][j] for j in cols] for i in range(d)])
    BRsub = Matrix([[BR[i][j] for j in cols] for i in range(d)])
    C = (BRsub * Bsub.inv_mod(p)) % p
    # check BR = C B on all columns
    CB = (C * Bm) % p
    okc = all((CB[i, j] - BR[i][j]) % p == 0 for i in range(d) for j in range(len(tup)))
    from sympy import symbols, Poly
    X = symbols('X')
    cp = Poly(C.charpoly(X).as_expr(), X)
    coeffs = [int(c) % p for c in cp.all_coeffs()]
    return d, okc, coeffs

for lab, S in (("QR {1,2,4}", {1, 2, 4}), ("NQR {3,5,6}", {3, 5, 6})):
    d, okc, cf = clock_charpoly(S)
    say(f"THM-4848 (v): P7 {lab}, k = 3, p = 79: dim of coker M_3 tensor F_79 = {d}; B R_tau = C B verified: {okc}; "
        f"charpoly of the translation clock C = X^3 + {cf[1]} X^2 + {cf[2]} X + {cf[3]} (mod 79)")

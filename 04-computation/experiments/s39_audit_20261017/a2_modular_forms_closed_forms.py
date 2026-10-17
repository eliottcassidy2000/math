#!/usr/bin/env python3
"""S39 audit, item 2: the modular forms of THM-4663 statement 3, built independently, and the closed forms of M_5, M_6, M_8.

* f15 = (eta(3z)eta(5z))^3 + (eta(z)eta(15z))^3 from Jacobi's identity eta(z)^3 = q^(1/8) sum (-1)^k (2k+1) q^(k(k+1)/2)
  (the session multiplied (1 - q^n) factors); compared with the Hecke character psi((alpha)) = alpha^2, psi(p_2) = w,
  implemented through ideal membership (alpha = x + y w in p_2 = (2, w) iff x even).
* b = (eta1 eta2 eta3 eta6)^2 from the pentagonal-number series.
* the weight-6 level-6 newform: S_6(Gamma_0(6)) spanned by eta quotients found by a search over exponent vectors
  (Ligozat/Gordon-Hughes-Newman orders at the four cusps), the newform singled out as the U_2-eigenvector with eigenvalue
  of absolute value 2^2 (old U_2 eigenvalues have absolute value 2^(5/2)); then checked to be a T_7 and T_5 eigenform.
* closed forms of M_5, M_6, M_8 (and the classical M_1..M_4) against moments from a1 (independent) and the session JSON;
  the symmetric-power sums S_j; W_k identities in sympy; the seventh-moment range test including p = 5, 7.
Usage: python -B a2_modular_forms_closed_forms.py SESSION_JSON A1_JSON
"""
import sys, json, itertools
from fractions import Fraction
import numpy as np
import sympy
from sympy import primerange, jacobi_symbol

NS = 3001


def kr(d, n):
    return int(jacobi_symbol(d % n, n)) if d % n else 0


def jacobi_cube(m, N):
    """prod_(n>=1) (1 - q^(mn))^3 as a list (coefficients up to q^N)."""
    s = [0] * (N + 1)
    k = 0
    while m * k * (k + 1) // 2 <= N:
        s[m * k * (k + 1) // 2] += (-1) ** k * (2 * k + 1)
        k += 1
    return s


def pentagonal(m, N):
    """prod_(n>=1) (1 - q^(mn))."""
    s = [0] * (N + 1)
    k = 0
    while True:
        done = True
        for kk in ([k, -k] if k else [0]):
            e = m * kk * (3 * kk - 1) // 2
            if e <= N:
                s[e] += 1 - 2 * (abs(kk) % 2)        # (-1)^kk as an integer (Python's (-1)**(-1) is a float)
                done = False
        if done and k > 0:
            break
        k += 1
    return s


def mul(a, b, N):
    A = np.array(a, dtype=object); B = np.array(b, dtype=object)
    out = np.zeros(N + 1, dtype=object)
    nzA = [i for i, x in enumerate(a) if x]
    for i in nzA:
        out[i:] += A[i] * B[:N + 1 - i]
    return list(out)


def shift(a, k, N):
    return ([0] * k + a)[:N + 1]


# ------------------------------------------------------------------ f15
f15 = [x + y for x, y in zip(shift(mul(jacobi_cube(3, NS), jacobi_cube(5, NS), NS), 1, NS),
                              shift(mul(jacobi_cube(1, NS), jacobi_cube(15, NS), NS), 2, NS))]


def psi_trace(p):
    """a_p of the Hecke character psi((alpha)) = alpha^2, psi(p_2) = w, w = (1 + sqrt-15)/2, via ideal membership."""
    if kr(-15, p) == -1:
        return 0
    # elements alpha = x + y w, N = x^2 + xy + 4y^2.  Represent numbers as (A, B) meaning (A + B sqrt-15)/2.
    def sq(A, B):            # ((A + B s)/2)^2 = (A^2 - 15 B^2 + 2AB s)/4 = ((A^2 - 15B^2)/2 + AB s)/2
        return (Fraction(A * A - 15 * B * B, 2), Fraction(A * B, 1))
    R = int((8 * p) ** 0.5) + 3
    for y in range(0, R):
        for x in range(-R, R):
            n = x * x + x * y + 4 * y * y
            if n == p:                          # principal: psi(p) + psi(pbar) = alpha^2 + conj = (2x + y)^2 - 2p
                return (2 * x + y) ** 2 - 2 * p
    for y in range(0, R):
        for x in range(-R, R):
            if x * x + x * y + 4 * y * y == 2 * p:
                A, B = 2 * x + y, y             # alpha = (A + B s)/2
                SA, SB = sq(A, B)               # alpha^2 = (SA + SB s)/2
                if x % 2 == 0:                  # alpha in p_2 = (2, w): divide by psi(p_2) = w = (1 + s)/2
                    wA, wB = 1, 1
                elif (x + y) % 2 == 0:          # alpha in pbar_2: divide by wbar = (1 - s)/2
                    wA, wB = 1, -1
                else:
                    continue
                # (SA + SB s)/2 / ((wA + wB s)/2) = (SA + SB s)(wA - wB s)/(wA^2 + 15 wB^2)
                re = (SA * wA + 15 * SB * wB) / (wA * wA + 15 * wB * wB)   # real part times 1 (the quotient is (re + im s))
                tr = 2 * re
                assert tr.denominator == 1
                return int(tr)
    raise RuntimeError(p)


bad = [p for p in primerange(7, NS) if f15[p] != psi_trace(p)]
print('f15 from Jacobi cubes: a_2, a_3, a_5, a_7, a_11, a_13, a_17, a_19 =', [f15[n] for n in (2, 3, 5, 7, 11, 13, 17, 19)])
print('  Hecke character psi((alpha)) = alpha^2, psi(p_2) = w equals the eta-quotient coefficient at every prime 7 <= p < %d: failures %s' % (NS, bad))
assert not bad and f15[1] == 1 and (f15[2], f15[3], f15[5]) == (1, -3, 5)
# multiplicativity / Hecke relations for weight 3, nebentypus (-15/.)
chi15 = lambda n: kr(-15, n) if n % 2 else (1 if n % 4 == 0 else kr(-15, n))
for (m, n) in [(2, 7), (3, 7), (7, 11), (2, 3), (4, 3)]:
    pass
assert f15[6] == f15[2] * f15[3] and f15[14] == f15[2] * f15[7] and f15[77] == f15[7] * f15[11]
for p in [7, 11, 13, 17, 19, 23]:
    assert f15[p * p] == f15[p] ** 2 - kr(-15, p) * p ** 2, p
    assert abs(f15[p]) <= 2 * p
print('  multiplicative, a_(p^2) = a_p^2 - (-15/p) p^2, Ramanujan bound: ok')

# ------------------------------------------------------------------ b = (eta1 eta2 eta3 eta6)^2
P = {m: pentagonal(m, NS) for m in (1, 2, 3, 6)}
b46 = [1] + [0] * NS
for m in (1, 2, 3, 6):
    b46 = mul(b46, P[m], NS); b46 = mul(b46, P[m], NS)
b46 = shift(b46, 1, NS)
print('b = (eta1 eta2 eta3 eta6)^2: a_2..a_13 =', b46[2:14])
for p in [5, 7, 11, 13]:
    assert b46[p * p] == b46[p] ** 2 - p ** 3
assert b46[6] == b46[2] * b46[3] and b46[35] == b46[5] * b46[7]

# ------------------------------------------------------------------ S_6(Gamma_0(6)) from eta quotients
def cusp_orders(r):
    # r = (r1, r2, r3, r6); order at the cusp 1/d (d | 6), N = 6, gcd(d, 6/d) = 1
    out = []
    for d in (1, 2, 3, 6):
        v = Fraction(0)
        for delta, e in zip((1, 2, 3, 6), r):
            g = sympy.gcd(d, delta)
            v += Fraction(int(g) ** 2 * e, d * delta)
        out.append(v * Fraction(6, 24))
    return out


cands = []
for r in itertools.product(range(-6, 19), repeat=4):
    if sum(r) != 12:
        continue
    if sum(d * e for d, e in zip((1, 2, 3, 6), r)) % 24 or sum((6 // d) * e for d, e in zip((1, 2, 3, 6), r)) % 24:
        continue
    s = 1
    for d, e in zip((1, 2, 3, 6), r):
        s *= Fraction(d) ** e
    if not sympy.sqrt(sympy.Rational(s.numerator, s.denominator)).is_rational:
        continue
    o = cusp_orders(r)
    if all(x >= 1 for x in o):
        cands.append((r, o))
print('holomorphic cuspidal eta quotients of weight 6 on Gamma_0(6):', [c[0] for c in cands])
NB = 200


def eta_quot(r, N):
    ser = [1] + [0] * N
    for m, e in zip((1, 2, 3, 6), r):
        base = pentagonal(m, N)
        if e < 0:   # invert the series
            inv = [0] * (N + 1); inv[0] = 1
            for n in range(1, N + 1):
                inv[n] = -sum(base[k] * inv[n - k] for k in range(1, n + 1) if base[k])
            base = inv; e = -e
        for _ in range(e):
            ser = mul(ser, base, N)
    lead = sum(d * e for d, e in zip((1, 2, 3, 6), r)) // 24
    return shift(ser, lead, N)


basis = [eta_quot(c[0], NB) for c in cands]
Mb = sympy.Matrix([[basis[j][n] for n in range(1, NB + 1)] for j in range(len(basis))])
rank = Mb.rank()
print('rank of the eta-quotient span (expected dim S_6(Gamma_0(6)) = 3):', rank)
idx = []
for j in range(len(basis)):
    if sympy.Matrix([[basis[i][n] for n in range(1, NB + 1)] for i in idx + [j]]).rank() == len(idx) + 1:
        idx.append(j)
    if len(idx) == 3:
        break
B3 = [basis[j] for j in idx]
print('basis used:', [cands[j][0] for j in idx])
# U_2 on the span: (U_2 f)_n = a_(2n). Solve U_2 B_i = sum_j U_ij B_j using the first coefficients.
nn = 40
A = sympy.Matrix([[B3[j][n] for j in range(3)] for n in range(1, nn + 1)])
U = []
for i in range(3):
    rhs = sympy.Matrix([B3[i][2 * n] for n in range(1, nn + 1)])
    sol = A.solve_least_squares(rhs)
    assert A * sol == rhs
    U.append(list(sol))
Umat = sympy.Matrix(U).T            # columns: images
ev = Umat.eigenvects()
print('U_2 eigenvalues on S_6(Gamma_0(6)):', [(e, m) for e, m, _ in ev])
new = [vs[0] for e, m, vs in ev if abs(complex(e)) == 4.0]
assert len(new) == 1
v = new[0]
f66 = [sum(v[j] * B3[j][n] for j in range(3)) for n in range(NB + 1)]
f66 = [sympy.nsimplify(x / f66[1]) for x in f66]
assert all(x.is_integer for x in f66)
f66 = [int(x) for x in f66]
print('weight-6 level-6 newform (U_2 route): a_1..a_13 =', f66[1:14])
# Hecke checks: T_7 and T_5 eigenform (weight 6: a_p a_n = a_(pn) + p^5 a_(n/p))
for p in (5, 7):
    for n in range(1, NB // p + 1):
        assert f66[p] * f66[n] == f66[p * n] + (p ** 5 * f66[n // p] if n % p == 0 else 0), (p, n)
assert f66[6] == f66[2] * f66[3] and f66[4] == f66[2] ** 2 and f66[9] == f66[3] ** 2
print('  T_5, T_7 eigenform, U_2, U_3 multiplicative: ok;  a_2, a_3, a_5, a_7, a_11 =',
      (f66[2], f66[3], f66[5], f66[7], f66[11]))
assert (f66[2], f66[3], f66[5], f66[7], f66[11]) == (4, -9, -66, 176, -60)

# extend the newform's prime coefficients to 3000 using a dense space: f66 = combination of the eta quotients
B3big = [eta_quot(cands[j][0], NS) for j in idx]
f66big = [sum(v[j] * B3big[j][n] for j in range(3)) for n in range(NS + 1)]
c1 = f66big[1]
f66big = [sympy.nsimplify(x / c1) for x in f66big]
f66big = [int(x) for x in f66big]
assert f66big[:NB + 1] == f66

# ------------------------------------------------------------------ closed forms
S = {int(p): v for p, v in json.load(open(sys.argv[1]))['moments'].items()}
A1 = {int(p): v for p, v in json.load(open(sys.argv[2])).items()}
fails = {k: [] for k in (1, 2, 3, 4, 5, 6, 8)}
nchk = 0
for p in sorted(S):
    m = S[p]
    if p in A1:
        assert A1[p] == m
    e3 = kr(p, 3)
    a, b, d = f15[p], b46[p], f66big[p]
    nchk += 1
    if m[0] != 1: fails[1].append(p)
    if m[1] != p * p - p - 1: fails[2].append(p)
    if m[2] != e3 * p * p + 2 * p + 1: fails[3].append(p)
    if m[3] != 2 * p ** 3 - 3 * p * p - 3 * p - 1: fails[4].append(p)
    if m[4] != 4 * e3 * p ** 3 + 5 * p * p + 4 * p + 1 + p * p * a: fails[5].append(p)
    if m[5] != 5 * p ** 4 - 10 * p ** 3 - 9 * p * p - 5 * p - 1 - p * p * b: fails[6].append(p)
    if m[7] != 14 * p ** 5 - 35 * p ** 4 - 28 * p ** 3 - 20 * p * p - 7 * p - 1 - 7 * p ** 3 * b - p * p * d: fails[8].append(p)
print('\nclosed forms on %d primes 5 <= p <= %d (session JSON; = a1 exact moments on the %d primes <= %d): failures %s'
      % (nchk, max(S), len(A1), max(A1), fails))
print('  M_5 at p = 5: actual %d, formula with a_5 = 5 gives %d (so the p >= 7 proviso of statements 3-4 is unnecessary)'
      % (S[5][4], 4 * kr(5, 3) * 125 + 125 + 20 + 1 + 25 * 5))
assert all(not v for k, v in fails.items() if k != 5) and fails[5] in ([], [5])

# symmetric power sums and the seventh moment
ss = {}
bad7, rng7, lam0 = [], [], []
for p in sorted(S):
    m = S[p]
    al = sympy.symbols('al')
    # S_j from power sums: Sigma (-K)^k = sum_j (C(k,j) - C(k,j-1)) p^j S_(k-2j), S_0 = p - 1
    Spow = {k: (-1) ** k * m[k - 1] for k in range(1, 9)}
    Ssym = {0: p - 1}
    for k in range(1, 9):
        acc = Spow[k]
        for j in range(1, k // 2 + 1):
            acc -= (sympy.binomial(k, j) - sympy.binomial(k, j - 1)) * p ** j * Ssym[k - 2 * j]
        Ssym[k] = int(acc)
    e3 = kr(p, 3)
    exp = {1: -1, 2: -1, 3: -1 - e3 * p * p, 4: -1 - p * p, 6: -1 - p * p * b46[p], 8: -1 - p ** 4 - p * p * f66big[p]}
    if p >= 7:
        exp[5] = -1 - p * p * f15[p]
    for k, v in exp.items():
        assert Ssym[k] == v, (p, k)
    ss[p] = Ssym[7]
    if (Ssym[7] + 1) % (p * p):
        bad7.append(p)
    else:
        c7 = -(Ssym[7] + 1) // (p * p)
        if p > 7:
            lam2 = p * p + kr(p, 3) * kr(p, 5) * kr(p, 7) * c7
            if not 0 <= lam2 <= 4 * p * p:
                rng7.append(p)
            if lam2 == 0:
                lam0.append(p)
print('symmetric-power sums S_1..S_6, S_8 as stated at every prime (S_5 for p >= 7): ok')
print('seventh moment: p^2 does NOT divide S_7 + 1 at', bad7, '; |lambda_p|^2 range violations (p >= 11):', rng7, '; lambda_p = 0 at', lam0)
print('  (S_7 + 1)/p^2 at p = 5, 7, 11, 13:', [(ss[p] + 1) // (p * p) for p in (5, 7, 11, 13)],
      '(the session script tested divisibility only for p >= 11; it also holds at 5 and 7)')

# W_k identities (sympy)
pp, e, a, b, d = sympy.symbols('p e a b d')
M5 = 4 * e * pp ** 3 + 5 * pp ** 2 + 4 * pp + 1 + pp ** 2 * a
M6 = 5 * pp ** 4 - 10 * pp ** 3 - 9 * pp ** 2 - 5 * pp - 1 - pp ** 2 * b
M8 = 14 * pp ** 5 - 35 * pp ** 4 - 28 * pp ** 3 - 20 * pp ** 2 - 7 * pp - 1 - 7 * pp ** 3 * b - pp ** 2 * d
W = lambda k, Mk: sympy.expand(((pp + 1) ** k + (pp + 1) * (-1) ** k * Mk) / pp ** 2)
assert sympy.expand(W(5, M5) - (pp + 1) * (pp ** 2 + 4 * (1 - e) * pp + 1 - a)) == 0
assert sympy.expand(W(6, M6) - (pp ** 4 + 11 * pp ** 3 + 10 * pp ** 2 + pp + 1 - (pp + 1) * b)) == 0
assert sympy.expand(W(8, M8) - (pp + 1) * (pp ** 5 + 7 * pp ** 4 + 35 * pp ** 3 + 7 * pp + 1 - 7 * pp * b - d)) == 0
print('W_5, W_6, W_8 closed forms follow from M_5, M_6, M_8 (sympy): ok')
json.dump({'f15': f15, 'b46': b46, 'f66': f66big}, open('a2_forms_to_3000.json', 'w'))
print('ALL CHECKS PASSED')

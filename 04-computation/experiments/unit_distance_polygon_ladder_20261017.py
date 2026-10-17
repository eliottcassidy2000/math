#!/usr/bin/env python3
"""The Kloosterman polygon ladder of the finite unit-distance graphs G_q = Cay(F_(q^2), mu_(q+1)) (S39).

W_k(q) = number of closed k-walks from a vertex of G_q = number of k-tuples of norm-one steps summing to zero
       = q^(-2) [ (q+1)^k + (q+1) sum_(a in F_q^*) (-K_q(a))^k ]   (spectrum = {q+1} u {-K_q(N(b))}, THM-4628).

Checks, all exact unless marked:
  1. the classical moments M_1..M_4 and the closed forms of M_5 (weight-3 CM form of Q(sqrt-15)), M_6 (eta(z)^2 eta(2z)^2
     eta(3z)^2 eta(6z)^2), M_8 (that form and the weight-6 level-6 newform), on every prime 5 <= p <= 2999 (exact moments from
     kloosterman_moments_exact_20261017.py); the seventh moment against Evans' |lambda_p|^2 shape (a range test only);
  2. the resulting closed forms of W_k(p), k <= 8, against traces of powers of the adjacency matrix of G_p, p <= 17;
  3. the pentagon count of G_p (Harary-Manvel for regular graphs) against a direct cycle count, p <= 13;
  4. prime powers q = 49, 121, 169: moments over F_q against the same closed forms with Frobenius powers;
  5. NUMERICAL (25 digits): the five- and six-Bessel moments and the 5-step planar walk density against L-values.
Usage: python unit_distance_polygon_ladder_20261017.py moments.json
Session opus-2026-10-17-S39.
"""
import sys, json
from math import isqrt
import numpy as np
import sympy
from sympy import jacobi_symbol, factorint, primerange, divisor_sigma, Matrix, Rational

MOM = json.load(open(sys.argv[1]))['moments']
M = {int(p): v for p, v in MOM.items()}
PS = sorted(M)


def kr(d, n):
    return int(jacobi_symbol(d % n, n)) if d % n else 0


# ---------------------------------------------------------------- the three modular forms
def _mul(x, y):                    # O_K = Z[w], w = (1 + sqrt(-15))/2, w^2 = w - 4
    a, b = x; c, d = y
    return (a * c - 4 * b * d, a * d + b * c + b * d)


def a15(p):
    """p-th coefficient of the weight-3 CM newform of level 15: psi((alpha)) = alpha^2, psi(p_2) = w (p_2^2 = (w))."""
    if p == 2: return 1
    if p == 3: return -3
    if p == 5: return 5
    if kr(-15, p) == -1: return 0
    r = isqrt(p) + 3
    for y in range(0, r):
        for x in range(-r, r):
            if x * x + x * y + 4 * y * y == p:          # principal class
                return (2 * x + y) ** 2 - 2 * p
    r = isqrt(8 * p) + 3
    for b in range(-r, r):                              # nonprincipal class: (alpha) = p p_2 or p conj(p_2), N(alpha) = 2p
        for a in range(-r, r):
            if a * a + a * b + 4 * b * b == 2 * p:
                al2 = _mul((a, b), (a, b))
                for ws in [(1, -1), (0, 1)]:            # divide by w or by conj(w): x / w = x conj(w) / 4
                    num = _mul(al2, ws)
                    if num[0] % 4 == 0 and num[1] % 4 == 0:
                        return 2 * (num[0] // 4) + num[1] // 4
    raise RuntimeError(p)


NS = 3001


def eta_series(exps, N=NS):
    s = np.zeros(N + 1, dtype=object); s[0] = 1
    for m, e in exps:
        for n in range(1, N // m + 1):
            for _ in range(e):
                s[m * n:] = s[m * n:] - s[:N + 1 - m * n].copy()
    return s


def mulser(a, b, N=NS):
    c = np.zeros(N + 1, dtype=object)
    for i in range(N + 1):
        if a[i]:
            c[i:] += a[i] * b[:N + 1 - i]
    return c


P46 = eta_series([(1, 2), (2, 2), (3, 2), (6, 2)])
f46 = np.zeros(NS + 1, dtype=object); f46[1:] = P46[:NS]          # eta(z)^2 eta(2z)^2 eta(3z)^2 eta(6z)^2 = 6.4.a.a
sig = [0] + [int(divisor_sigma(n, 1)) for n in range(1, NS + 1)]


def E2(d):
    s = np.zeros(NS + 1, dtype=object); s[0] = 1
    for n in range(1, NS // d + 1):
        s[d * n] = -24 * sig[n]
    return s


# S_6(Gamma_0(6)) = f46 * M_2(Gamma_0(6)) (dimension 3: two oldforms from level 3, one newform)
basis6 = [mulser(f46, E2(1) - d * E2(d)) for d in (2, 3, 6)]
Bm = Matrix([[basis6[j][n] for j in range(3)] for n in range(1, 13)])
TBm = Matrix([[basis6[j][5 * n] + (5 ** 5 * basis6[j][n // 5] if n % 5 == 0 else 0) for j in range(3)] for n in range(1, 13)])
T5 = (Bm.T * Bm).LUsolve(Bm.T * TBm)
newvec = [v for val, mult, vecs in T5.eigenvects() if mult == 1 for v in vecs]
assert len(newvec) == 1
f66 = sum(newvec[0][j] * basis6[j] for j in range(3))
f66 = [Rational(x) / Rational(f66[1]) for x in f66]
assert all(x.q == 1 for x in f66)
f66 = [int(x) for x in f66]
print('weight-6 level-6 newform: a_2..a_7 =', f66[2:8], '(T_5 eigenvalue %s, old eigenvalue %s)' %
      (f66[5], [val for val, mult, vecs in T5.eigenvects() if mult == 2][0]))
print('eta(z)^2 eta(2z)^2 eta(3z)^2 eta(6z)^2: a_2..a_7 =', [int(f46[n]) for n in range(2, 8)])

A15 = {p: a15(p) for p in PS}
B6 = {p: int(f46[p]) for p in PS}
D6 = {p: f66[p] for p in PS}

# the Hecke-character coefficients against Broadhurst's eta-quotient (eta(3z)eta(5z))^3 + (eta(z)eta(15z))^3
_A = mulser(eta_series([(3, 3)]), eta_series([(5, 3)]))
_B = mulser(eta_series([(1, 3)]), eta_series([(15, 3)]))
f315 = np.zeros(NS + 1, dtype=object); f315[1:] += _A[:NS]; f315[2:] += _B[:NS - 1]
bad15 = [p for p in primerange(2, NS) if int(f315[p]) != a15(p)]
print('level-15 CM coefficients = coefficients of (eta3 eta5)^3 + (eta1 eta15)^3 at every prime < %d; differences at' % NS, bad15)
assert not bad15

# ---------------------------------------------------------------- 1. moment closed forms
fails = {k: [] for k in range(1, 9)}
for p in PS:
    m = M[p]
    if m[0] != 1: fails[1].append(p)
    if m[1] != p * p - p - 1: fails[2].append(p)
    if m[2] != kr(p, 3) * p * p + 2 * p + 1: fails[3].append(p)
    if m[3] != 2 * p ** 3 - 3 * p * p - 3 * p - 1: fails[4].append(p)
    if p >= 7 and m[4] != 4 * kr(p, 3) * p ** 3 + 5 * p * p + 4 * p + 1 + p * p * A15[p]: fails[5].append(p)
    if m[5] != 5 * p ** 4 - 10 * p ** 3 - 9 * p * p - 5 * p - 1 - p * p * B6[p]: fails[6].append(p)
    if m[7] != 14 * p ** 5 - 35 * p ** 4 - 28 * p ** 3 - 20 * p * p - 7 * p - 1 - 7 * p ** 3 * B6[p] - p * p * D6[p]:
        fails[8].append(p)
print('\n1. closed forms of M_k(p) on all %d primes 5 <= p <= %d (p >= 7 for k = 5): failures' % (len(PS), PS[-1]),
      {k: v for k, v in fails.items() if k != 7})
assert all(not v for v in fails.values())
# seventh moment: S_7 = -1 - p^2 c_7(p) with c_7(p) = (p/105)(|lambda_p|^2 - p^2) (Evans, Yun; Broadhurst's form).
viol, lam0 = [], []
for p in PS:
    if p <= 7: continue
    S5 = -1 - p * p * A15[p]; S3 = -1 - kr(p, 3) * p * p; S1 = -1
    S7 = -M[p][6] - 6 * p * S5 - 14 * p * p * S3 - 14 * p ** 3 * S1
    assert (S7 + 1) % (p * p) == 0
    c7 = -(S7 + 1) // (p * p)
    lam2 = p * p + kr(p, 3) * kr(p, 5) * kr(p, 7) * c7
    if not 0 <= lam2 <= 4 * p * p: viol.append(p)
    if lam2 == 0: lam0.append(p)
print('   seventh moment: S_7 + 1 divisible by p^2 for every p; implied |lambda_p|^2 in [0, 4p^2] fails at', viol,
      '; lambda_p = 0 at', lam0, '(RANGE TEST ONLY: lambda_p is not computed independently)')


# ---------------------------------------------------------------- 2. closed forms of W_k against adjacency traces
def W_closed(p):
    e3 = kr(p, 3)
    W = {2: p + 1, 3: (p + 1) * (1 - e3), 4: 3 * p * (p + 1),
         6: p ** 4 + 11 * p ** 3 + 10 * p * p + p + 1 - (p + 1) * B6[p],
         8: (p + 1) * (p ** 5 + 7 * p ** 4 + 35 * p ** 3 + 7 * p + 1 - 7 * p * B6[p] - D6[p])}
    if p >= 7:
        W[5] = (p + 1) * (p * p + 4 * (1 - e3) * p + 1 - A15[p])
    return W


def W_from_moments(p, k):
    num = (p + 1) ** k + (p + 1) * (-1) ** k * M[p][k - 1]
    assert num % (p * p) == 0
    return num // (p * p)


for p in PS:
    for k, v in W_closed(p).items():
        assert v == W_from_moments(p, k), (p, k)
print('\n2. closed forms of W_k(p), k in {2,3,4,5,6,8}, agree with the moment formula on all primes; W_5:',
      {p: W_closed(p)[5] for p in [7, 11, 13, 17, 19]})


def unit_graph(p):
    n = next(x for x in range(2, p) if kr(x, p) == -1)        # F_(p^2) = F_p(sqrt n)
    steps = [(x, y) for x in range(p) for y in range(p) if (x * x - n * y * y) % p == 1]
    assert len(steps) == p + 1
    A = np.zeros((p * p, p * p), dtype=np.int64)
    for x in range(p):
        for y in range(p):
            for dx, dy in steps:
                A[x * p + y, ((x + dx) % p) * p + (y + dy) % p] = 1
    return A, steps


for p in [5, 7, 11, 13, 17]:
    A, steps = unit_graph(p)
    Ak = np.eye(p * p, dtype=np.int64)
    got = []
    for k in range(1, 9):
        Ak = Ak @ A
        tr = int(np.trace(Ak))
        assert tr % (p * p) == 0
        got.append(tr // (p * p))
        assert got[-1] == W_from_moments(p, k), (p, k)
    print('   p = %2d: adjacency traces / p^2 = %s  (= moment formula)' % (p, got))


# ---------------------------------------------------------------- 3. pentagons
def count_c5_through0(A):
    nbr = [np.nonzero(A[v])[0] for v in range(A.shape[0])]
    cnt = 0
    for a in nbr[0]:
        for b in nbr[a]:
            if b == 0: continue
            for c in nbr[b]:
                if c in (0, a): continue
                for d in nbr[c]:
                    if d in (0, a, b): continue
                    if A[d, 0]:
                        cnt += 1
    return cnt // 2                     # each 5-cycle through 0 is traversed in two directions


for p in [5, 7, 11, 13]:
    A, _ = unit_graph(p)
    c5 = p * p * count_c5_through0(A) // 5
    W = W_closed(p) if p >= 7 else {3: W_from_moments(p, 3), 5: W_from_moments(p, 5)}
    formula = p * p * (W[5] - 5 * p * W[3]) // 10             # Harary-Manvel, (p+1)-regular graph
    unified = p * p * (p + 1) * (p * p + 1 - (1 - kr(p, 3)) * p - A15[p]) // 10 if p >= 7 else None
    print('\n3. p = %2d: 5-cycles in G_p by direct count %d, by the trace formula %d, closed form %s' % (p, c5, formula, unified))
    assert c5 == formula and (unified is None or unified == c5)


# ---------------------------------------------------------------- 3b. theta graphs and the four-bar linkage curve
def radial_walks(p, k, A):
    """R_k(d): number of k-walks from 0 to d, as a vector over the vertices of G_p."""
    e0 = np.zeros(p * p, dtype=object); e0[0] = 1
    v = e0
    for _ in range(k):
        v = A.T.astype(object).dot(v)
    return v


def ap_30a1(p):
    cnt = 1
    for x in range(p):
        for y in range(p):
            if (y * y + x * y + y - (x ** 3 + x + 2)) % p == 0:
                cnt += 1
    return p + 1 - cnt


def fourbar_count(p):
    """R_3(4) = #{(u1,u2,u3) in mu^3 : u1 + u2 + u3 = 2} (three unit steps from 0 to the point 2, N(2) = 4)."""
    n = next(x for x in range(2, p) if kr(x, p) == -1)
    mu = [(x, y) for x in range(p) for y in range(p) if (x * x - n * y * y) % p == 1]
    mus = set(mu)
    return sum(1 for (a, b) in mu for (c, d) in mu if ((2 - a - c) % p, (-b - d) % p) in mus)


for p in [7, 11, 13, 17, 19]:
    A, _ = unit_graph(p)
    R2, R3 = radial_walks(p, 2, A), radial_walks(p, 3, A)
    th223 = int(sum(R2[i] * R2[i] * R3[i] for i in range(p * p)))
    W = {k: W_from_moments(p, k) for k in (3, 5)}
    pred = 2 * W[5] + (p * p - 1) * W[3] - (p + 1) * (p + 1 - ap_30a1(p))
    print('\n3b. p = %2d: rooted hom(theta(2,2,3), G_p) = %d; 2W_5 + (p^2-1)W_3 - (p+1)(p+1-a_p(30a1)) = %d' % (p, th223, pred))
    assert th223 == pred
badfb = [p for p in primerange(7, 400) if fourbar_count(p) != p + 1 - ap_30a1(p)]
print('   four-bar count R_3(4) = p + 1 - a_p(30a1) for every prime 7 <= p < 400; failures:', badfb)
assert not badfb
# j-invariant of the four-bar curve: v = t u through (0,0) gives w^2 = t^4 - 6t^3 - 5t^2 - 6t + 1
a, b, c, d, e = 1, -6, -5, -6, 1
I = 12 * a * e - 3 * b * d + c * c
J = 72 * a * c * e + 9 * b * c * d - 27 * a * d * d - 27 * e * b * b - 2 * c ** 3
jq = Rational(6912 * I ** 3, 4 * I ** 3 - J * J)
c4, c6, Dl = -71, -1837, -2160                     # 30a1 = [1,0,1,1,2]: b2 = 1, b4 = 3, b6 = 9, b8 = 0
assert c4 ** 3 - c6 ** 2 == 1728 * Dl
print('   j(four-bar quartic) = %s = j(30a1) = c4^3/Delta = %s' % (jq, Rational(c4 ** 3, Dl)))
assert jq == Rational(c4 ** 3, Dl)


def kl_float(p):
    g = sympy.primitive_root(p)
    pw = np.array([pow(g, i, p) for i in range(p - 1)])
    f = np.exp(2j * np.pi * pw / p)
    K = np.fft.ifft(np.fft.fft(f) ** 2).real
    Kf = np.zeros(p); Kf[pw] = K
    return Kf


def radial_float(p, k, Kf):
    """R_k as a function of the norm class c (c = 0: the vertex 0), rounded from Kloosterman sums."""
    a = np.arange(1, p)
    base = (-Kf[a]) ** k
    R = np.zeros(p)
    R[0] = ((p + 1) ** k + (p + 1) * base.sum()) / p ** 2
    for c in range(1, p):
        R[c] = ((p + 1) ** k + (base * (-Kf[(a * c) % p])).sum()) / p ** 2
    Ri = np.rint(R)
    assert np.abs(R - Ri).max() < 1e-3
    return [int(x) for x in Ri]


thbad = []
for p in primerange(7, 300):
    Kf = kl_float(p)
    R = {k: radial_float(p, k, Kf) for k in (1, 2, 3)}
    th = lambda x, y, z: R[x][0] * R[y][0] * R[z][0] + (p + 1) * sum(R[x][c] * R[y][c] * R[z][c] for c in range(1, p))
    e3 = kr(p, 3)
    W5, W3 = W_closed(p)[5], W_closed(p)[3]
    checks = [th(2, 2, 2) == p ** 3 + 7 * p * p + 4 * p - 2, th(1, 2, 3) == 3 * (1 - e3) * p * (p + 1),
              th(1, 3, 3) == 9 * p * p * (p + 1),
              th(2, 2, 3) == 2 * W5 + (p * p - 1) * W3 - (p + 1) * (p + 1 - ap_30a1(p)),
              R[3][4 % p] == p + 1 - ap_30a1(p)]
    if not all(checks):
        thbad.append((p, checks))
print('   theta(2,2,2), theta(1,2,3), theta(1,3,3) polynomial and the theta(2,2,3) four-bar identity at every prime '
      '7 <= p < 300 (radial Kloosterman evaluation, rounded); failures:', thbad)
assert not thbad

# ---------------------------------------------------------------- 3c. the four-bar family = the Gamma_1(6) surface
from sympy import symbols, expand, Poly, factor, cancel
tt, cc = symbols('t c')
Dlt = expand((tt ** 2 + (1 + cc) * tt + 1) ** 2 - 4 * cc * tt * (tt + 1) ** 2)
co = [Poly(Dlt, tt).coeff_monomial(tt ** k) for k in (4, 3, 2, 1, 0)]
Ic = expand(12 * co[0] * co[4] - 3 * co[1] * co[3] + co[2] ** 2)
Jc = expand(72 * co[0] * co[2] * co[4] + 9 * co[1] * co[2] * co[3] - 27 * co[0] * co[3] ** 2 - 27 * co[4] * co[1] ** 2 - 2 * co[2] ** 3)
Dc_ = factor(4 * Ic ** 3 - Jc ** 2)
jc = factor(cancel(6912 * Ic ** 3 / (4 * Ic ** 3 - Jc ** 2)))
print('\n3c. four-bar family w^2 = Delta_c(t): 4I^3 - J^2 =', Dc_, '; j(c) =', jc)
assert expand(4 * Ic ** 3 - Jc ** 2 - 6912 * cc ** 2 * (cc - 9) * (cc - 1) ** 3) == 0
assert all(Ic.subs(cc, r) != 0 for r in (0, 1, 9))          # multiplicative fibres I_2, I_3, I_1 at 0, 1, 9; I_6 at infinity
fam_bad, mod6_bad = [], []
for p in primerange(7, 102):
    Kf = kl_float(p)
    R3 = radial_float(p, 3, Kf)
    for cv in range(1, p):
        if cv % p in (1, 9 % p):
            continue
        Ec = 2 + sum(1 + kr(((tv * tv + (1 + cv) * tv + 1) ** 2 - 4 * cv * tv * (tv + 1) ** 2) % p, p) for tv in range(p))
        if R3[cv] != Ec:
            fam_bad.append((p, cv))
        if R3[cv] % 6:
            mod6_bad.append((p, cv))
print('   R_3(c) = #E_c(F_p) for every prime 7 <= p <= 101 and every c not in {0, 1, 9}: failures', fam_bad[:5],
      '; R_3(c) = 0 mod 6: failures', mod6_bad[:5])
assert not fam_bad and not mod6_bad

# ---------------------------------------------------------------- 3d. moments of the four-bar family
mom_bad = []
for p in primerange(11, 400):
    S1 = S2 = 0
    for cv in range(1, p):
        if cv in (1, 9 % p):
            continue
        ac = p + 1 - (2 + sum(1 + kr(((tv * tv - 4 * cv) * (tv - cv - 1) * (tv - cv + 3)) % p, p) for tv in range(p)))
        S1 += ac; S2 += ac * ac
    if S1 != -2 - 2 * kr(p, 3) or S2 != (p + 1) * (p - 4) - B6[p]:
        mom_bad.append(p)
    if p < 102:
        Kf = kl_float(p); R3 = radial_float(p, 3, Kf)
        assert R3[1] == 3 * p and R3[9 % p] == p + 1 - kr(p, 3)
print('   sum_c a_p(E\'_c) = -2 - 2(p/3) and sum_c a_p(E\'_c)^2 = (p+1)(p-4) - b_p (b = eta1^2 eta2^2 eta3^2 eta6^2) '
      'at every prime 11 <= p < 400; failures:', mom_bad, '; R_3(1) = 3p, R_3(9) = p + 1 - (p/3) for p <= 101')
assert not mom_bad


# ---------------------------------------------------------------- 4. prime powers q = p^2
def moments_Fq2(p, kmax=8):
    """Kloosterman moments over F_q, q = p^2, F_q = F_p(sqrt n); trace F_q -> F_p is 2x for x + y sqrt n."""
    n = next(x for x in range(2, p) if kr(x, p) == -1)
    q = p * p
    els = [(x, y) for x in range(p) for y in range(p) if (x, y) != (0, 0)]
    def inv(z):
        x, y = z; d = (x * x - n * y * y) % p; di = pow(d, -1, p)
        return (x * di % p, (-y) * di % p)
    def mul2(z, w):
        return ((z[0] * w[0] + n * z[1] * w[1]) % p, (z[0] * w[1] + z[1] * w[0]) % p)
    ze = np.exp(2j * np.pi / p)
    K = {}
    for a in els:
        s = 0
        for x in els:
            t = mul2(a, inv(x))
            s += ze ** ((2 * (x[0] + t[0])) % p)
        K[a] = s.real
    return [round(sum(v ** k for v in K.values())) for k in range(1, kmax + 1)]


for p in [7, 11, 13]:
    q = p * p
    Mq = moments_Fq2(p)
    chi = kr(-15, p)
    A2 = A15[p] ** 2 - 2 * chi * p * p                     # Frobenius^2 trace, weight 3, nebentypus (-15/.)
    B2 = B6[p] ** 2 - 2 * p ** 3                            # weight 4, trivial character
    D2 = D6[p] ** 2 - 2 * p ** 5                            # weight 6
    e3 = 1                                                  # (q/3) = (p/3)^2 = 1
    pred = {5: 4 * e3 * q ** 3 + 5 * q * q + 4 * q + 1 + q * q * A2,
            6: 5 * q ** 4 - 10 * q ** 3 - 9 * q * q - 5 * q - 1 - q * q * B2,
            8: 14 * q ** 5 - 35 * q ** 4 - 28 * q ** 3 - 20 * q * q - 7 * q - 1 - 7 * q ** 3 * B2 - q * q * D2,
            3: e3 * q * q + 2 * q + 1, 4: 2 * q ** 3 - 3 * q * q - 3 * q - 1, 2: q * q - q - 1, 1: 1}
    ok = {k: Mq[k - 1] == v for k, v in pred.items()}
    print('\n4. q = %d: M_k(q) = %s; closed forms with Frobenius-squared traces hold: %s' % (q, Mq, ok))
    assert all(ok.values())


# ---------------------------------------------------------------- 5. Bessel moments and L-values (NUMERICAL)
import mpmath as mp
mp.mp.dps = 25


def coeffs_from_primes(ap_fn, chi_fn, w, bad, N):
    a = [0] * (N + 1); a[1] = 1
    for n in range(2, N + 1):
        v = 1
        for p, k in factorint(n).items():
            if p in bad:
                v *= ap_fn(p) ** k
            else:
                A = [1, ap_fn(p)]
                for j in range(2, k + 1):
                    A.append(ap_fn(p) * A[-1] - chi_fn(p) * p ** (w - 1) * A[-2])
                v *= A[k]
        a[n] = v
    return a


def L_values(a, level, w, points):
    Aq = mp.sqrt(level) / (2 * mp.pi)
    def Lam(s, eps):
        t = mp.mpf(0)
        for n in range(1, len(a)):
            if a[n]:
                x = 2 * mp.pi * n / mp.sqrt(level)
                t += a[n] * ((Aq / n) ** s * mp.gammainc(s, x) + eps * (Aq / n) ** (w - s) * mp.gammainc(w - s, x))
        return t
    s0 = mp.mpf(w) + mp.mpf('1.5')
    direct = mp.nsum(lambda n: 0, [1, 1])
    euler = mp.mpf(1)
    return {s: Lam(mp.mpf(s), 1) / (Aq ** s * mp.gamma(s)) for s in points}, Lam


a15n = coeffs_from_primes(a15, lambda p: 1 if p == 2 else kr(-15, p), 3, (3, 5), 400)
L15, Lam15 = L_values(a15n, 15, 3, [1, 2])
# root number check: Euler product at s = 3.5 against the approximate functional equation with eps = +1
s = mp.mpf('3.5'); E = mp.mpf(1)
for p in primerange(2, 20000):
    E *= 1 / (1 - a15(p) * mp.mpf(p) ** (-s)) if p in (3, 5) else \
        1 / (1 - a15(p) * mp.mpf(p) ** (-s) + (1 if p == 2 else kr(-15, p)) * mp.mpf(p) ** (2 - 2 * s))
afe = Lam15(s, 1) / ((mp.sqrt(15) / (2 * mp.pi)) ** s * mp.gamma(s))
print('\n5. level 15: Euler product L(3.5) = %s, functional equation (eps = +1) gives %s' % (mp.nstr(E, 12), mp.nstr(afe, 12)))
assert abs(E - afe) < 1e-8
a46n = [0] + [int(P46[n - 1]) for n in range(1, 301)]
L46, _ = L_values(a46n, 6, 4, [1, 2, 3])
J5 = mp.quadosc(lambda x: x * mp.besselj(0, x) ** 5, [0, mp.inf], omega=5)
IK4 = mp.quad(lambda t: t * mp.besseli(0, t) * mp.besselk(0, t) ** 4, [0, 1, 5, 20, mp.inf])
I2K3 = mp.quad(lambda t: t * mp.besseli(0, t) ** 2 * mp.besselk(0, t) ** 3, [0, 1, 5, 20, 60, mp.inf])
IK5 = mp.quad(lambda t: t * mp.besseli(0, t) * mp.besselk(0, t) ** 5, [0, 1, 5, 20, mp.inf])
I2K4 = mp.quad(lambda t: t * mp.besseli(0, t) ** 2 * mp.besselk(0, t) ** 4, [0, 1, 5, 20, mp.inf])
pi = mp.pi
rows = [('int_0^oo x J0(x)^5 dx  (5-step planar walk density p5\'(0))', J5, 6 * L15[1] / pi ** 2, '6 L(f15,1)/pi^2'),
        ('int_0^oo x J0(x)^5 dx  (Wick rotation)', J5, 30 * IK4 / pi ** 4, '30/pi^4 int t I0 K0^4'),
        ('int_0^oo t I0 K0^4 dt', IK4, pi ** 2 * L15[1] / 5, 'pi^2 L(f15,1)/5'),
        ('int_0^oo t I0^2 K0^3 dt', I2K3, 3 * L15[2] / 4, '3 L(f15,2)/4'),
        ('L(f15,1)/L(f15,2)', L15[1] / L15[2], mp.sqrt(15) / (2 * pi), 'sqrt(15)/(2 pi)'),
        ('int_0^oo t I0 K0^5 dt', IK5, pi ** 2 * L46[2] / 2, 'pi^2 L(f6,2)/2'),
        ('int_0^oo t I0^2 K0^4 dt', I2K4, pi ** 2 * L46[1] / 2, 'pi^2 L(f6,1)/2'),
        ('int_0^oo t I0^2 K0^4 dt', I2K4, 3 * L46[3] / 2, '3 L(f6,3)/2')]
for name, lhs, rhs, form in rows:
    print('   %-58s = %s ; %-18s = %s ; |diff| = %s' % (name, mp.nstr(lhs, 20), form, mp.nstr(rhs, 20), mp.nstr(abs(lhs - rhs), 3)))
    assert abs(lhs - rhs) < mp.mpf(10) ** -18
print('\nALL CHECKS PASSED')

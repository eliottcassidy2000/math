"""Audit S52 (independent code) - part C: pdf_referee.py claims, by different methods.

C1  q-series: tau(n) from the Jacobi triple product (prod (1-q^n)^3 = sum (-1)^k (2k+1) q^(k(k+1)/2)) to the 8th power,
    E4, E6 from divisor sums; Delta = (E4^3 - E6^2)/1728 cross-check; Leech theta = E4^3 - 720 Delta vs
    (65520/691)(sigma_11 - tau) for n <= 3000.
C2  extremal rank-48 theta from the basis {E4^6, E4^3 E6^2, E6^4} (not E4^3 Delta, Delta^2); shells; density;
    the PDF's covering-radius volume bound (24!/pi^24)^(1/48); Cohn-Elkies remark (interpolation of published LP data).
C3  L(Delta, s) via the Mellin integral split at TWO different points t0 (functional equation check), against the
    Manin/Zagier period relations r^+(X) ~ (36/691)(X^10 - 1) - X^2 (X^2-1)^3, r^-(X) ~ 4X^9 - 25X^7 + 42X^5 - 25X^3 + 4X.
C4  zeta_Leech(10) by the theta-function analytic continuation (Riemann's method, using ONLY the Leech theta
    coefficients, no L(Delta, s)): Z(2w) = pi^w/Gamma(w) [int_1^oo (theta(it)-1)(t^(w-1) + t^(11-w)) dt - 1/w - 1/(12-w)].
C5  zeta_Leech(26): closed form vs the theta-integral vs a direct partial sum to n = 3000 with an exact divisor-sum tail
    and a Deligne bound for the cusp tail.  zeta_E8(10) likewise.
C6  McLaughlin slice: brute-force binary even forms of det 15 with min >= 4; Prop. 6.1 over all fundamental
    discriminants -4000 < D < 0; the Leech count of norm-4 x with |v - x|^2 = 4 for a norm-6 v (2 * 276).
"""
import math, itertools, time
from fractions import Fraction as Fr
import mpmath as mp

mp.mp.dps = 60
FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


NMAX = 3000


def sigma_table(k, M):
    s = [0] * (M + 1)
    for d in range(1, M + 1):
        dk = d ** k
        for m in range(d, M + 1, d):
            s[m] += dk
    return s


def mul(a, b, M):
    c = [0] * (M + 1)
    nz = [(i, x) for i, x in enumerate(a) if x]
    for i, x in nz:
        for j in range(0, M + 1 - i):
            if b[j]:
                c[i + j] += x * b[j]
    return c


# C1 -------------------------------------------------------------------------------------------------------------
print('C1  q-series')
jac = [0] * (NMAX + 1)
k = 0
while k * (k + 1) // 2 <= NMAX:
    jac[k * (k + 1) // 2] += (-1) ** k * (2 * k + 1)
    k += 1
p8 = [1] + [0] * NMAX
for _ in range(8):
    p8 = mul(p8, jac, NMAX)          # dense x sparse
tau = [0] + p8[:NMAX]                # Delta = q * prod(1 - q^n)^24
check(tau[1:13] == [1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920, 534612, -370944],
      f'tau(1..12) from the Jacobi triple product: {tau[1:7]}...')
check(all(tau[m * n] == tau[m] * tau[n] for m in range(2, 40) for n in range(2, 40) if math.gcd(m, n) == 1 and m * n <= NMAX),
      'tau multiplicative on coprime pairs (m, n < 40)')
check(all(tau[p * p] == tau[p] ** 2 - p ** 11 for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53)),
      'Hecke relation tau(p^2) = tau(p)^2 - p^11 for p <= 53')
s3 = sigma_table(3, 60)
s5 = sigma_table(5, 60)
E4s = [1] + [240 * s3[n] for n in range(1, 61)]
E6s = [1] + [-504 * s5[n] for n in range(1, 61)]
E4c = mul(mul(E4s, E4s, 60), E4s, 60)
E6q = mul(E6s, E6s, 60)
check(all(E4c[n] - E6q[n] == 1728 * tau[n] for n in range(61)), 'Delta = (E4^3 - E6^2)/1728 through q^60')
s11 = sigma_table(11, NMAX)
s3b = sigma_table(3, NMAX)
E4 = [1] + [240 * s3b[n] for n in range(1, NMAX + 1)]
E4sq = mul(E4, E4, NMAX)
E4cube = mul(E4sq, E4, NMAX)
leech = [E4cube[n] - 720 * tau[n] for n in range(NMAX + 1)]
check(leech[:5] == [1, 0, 196560, 16773120, 398034000], f'Leech theta E4^3 - 720 Delta: {leech[:5]}')
check(all(691 * leech[n] == 65520 * (s11[n] - tau[n]) for n in range(1, NMAX + 1)),
      'Leech coefficients a(n) = (65520/691)(sigma_11(n) - tau(n)) for all n <= 3000 (two independent routes)')
print(f'     [{time.time() - T0:.0f}s]')

# C2 -------------------------------------------------------------------------------------------------------------
print('C2  extremal rank-48 theta series')
M = 8
E4m = E4s[:M + 1]
E6m = E6s[:M + 1]
b1 = mul(mul(E4m, E4m, M), mul(E4m, mul(E4m, mul(E4m, E4m, M), M), M), M)        # E4^6
b2 = mul(mul(mul(E4m, E4m, M), E4m, M), mul(E6m, E6m, M), M)                    # E4^3 E6^2
b3 = mul(mul(E6m, E6m, M), mul(E6m, E6m, M), M)                                 # E6^4
# solve x b1 + y b2 + z b3 with coefficients (1, 0, 0) at q^0, q^1, q^2
# exact rational solve
from sympy import Matrix, Rational
S = Matrix([[b1[i], b2[i], b3[i]] for i in range(3)]).LUsolve(Matrix([1, 0, 0]))
th48 = [sum(Rational(c) * b[i] for c, b in zip(S, (b1, b2, b3))) for i in range(M + 1)]
check([int(t) for t in th48[:7]] == [1, 0, 0, 52416000, 39007332000, 6609020221440, int(th48[6])] and all(t.q == 1 for t in th48),
      f'extremal theta (from E4^6, E4^3 E6^2, E6^4): shells {[int(t) for t in th48[:7]]}')
print(f'     coefficients in that basis: {list(S)}; a_6 (norm 12) = {int(th48[6])}')
check(sum(int(t) for t in th48[3:6]) == 6648079969440 and sum(int(t) for t in th48[3:6]) // 2 == 3324039984720,
      'norms 6, 8, 10 total 6,648,079,969,440 vectors = 3,324,039,984,720 antipodal pairs (a LOWER bound for the facet count, as the PDF says)')
dens = (3 * mp.pi / 2) ** 24 / mp.factorial(24)
check(abs(dens - mp.mpf('2.31782953e-8')) < mp.mpf('1e-16'), f'density (3 pi/2)^24/24! = {mp.nstr(dens, 12)}')
Rvol = (mp.factorial(24) / mp.pi ** 24) ** (mp.mpf(1) / 48)
check(abs(Rvol - mp.mpf('1.766471')) < 1e-6 and abs(4 * Rvol ** 2 - mp.mpf('12.481675')) < 1e-6,
      f'PDF Sec. 6 covering bound R_cov >= (24!/pi^24)^(1/48) = {mp.nstr(Rvol, 8)}, 4R^2 >= {mp.nstr(4 * Rvol ** 2, 9)} (correct)')
# Cohn-Elkies: Afkhami-Jeddi-Cohn-Hartman-de Laat-Tajdini (JHEP 2020) Table 1, LP bound in (log2 density)/d
lp = {8: -0.24737, 16: -0.33192, 32: -0.40382, 64: -0.46101, 128: -0.50432}


def quad_interp(d, pts):
    xs_ = [math.log2(p) for p in pts]
    ys_ = [lp[p] for p in pts]
    x = math.log2(d)
    tot = 0.0
    for i in range(3):
        li = 1.0
        for j in range(3):
            if j != i:
                li *= (x - xs_[j]) / (xs_[i] - xs_[j])
        tot += ys_[i] * li
    return tot


leech_dens = mp.pi ** 12 / mp.factorial(12)
f24 = quad_interp(24, (8, 16, 32))
check(abs(f24 - float(mp.log(leech_dens, 2)) / 24) < 5e-4,
      f'interpolation sanity: LP (log2 density)/d at d = 24 interpolates to {f24:.5f}; Leech (LP sharp) = {float(mp.log(leech_dens, 2)) / 24:.5f}')
f48a = quad_interp(48, (16, 32, 64))
f48b = quad_interp(48, (32, 64, 128))
lpd48 = [2 ** (48 * f) for f in (f48a, f48b)]
p48 = float(mp.log(dens, 2)) / 48
check(min(lpd48) > 10 * float(dens),
      f'LP bound at d = 48 interpolates to (log2 density)/48 = {f48a:.5f} / {f48b:.5f}, i.e. density ~ {lpd48[0]:.2e} / {lpd48[1]:.2e}, about {lpd48[0] / float(dens):.0f}x the extremal density ({p48:.5f}); the PDF value 3.01e-8 (and its "23% gap") is inconsistent with the published LP numerics')
print(f'     [{time.time() - T0:.0f}s]')

# C3 -------------------------------------------------------------------------------------------------------------
print('C3  L(Delta, s) with two split points; Manin/Zagier period relations')
NT = 80


def Lambda(s, t0):
    s = mp.mpf(s)
    tot = mp.mpf(0)
    for n in range(1, NT + 1):
        x = 2 * mp.pi * n
        tot += tau[n] * (mp.gammainc(s, x * t0) / x ** s + mp.gammainc(12 - s, x / t0) / x ** (12 - s))
    return tot


def Lval(s, t0=mp.mpf(1)):
    return Lambda(s, t0) * (2 * mp.pi) ** s / mp.gamma(s)


lam = {k: Lambda(k, mp.mpf('1.0')) for k in range(1, 12)}
lam2 = {k: Lambda(k, mp.mpf('1.37')) for k in range(1, 12)}
check(all(abs(lam[k] - lam2[k]) < mp.mpf(10) ** -45 * abs(lam[k]) for k in lam), 'Lambda(k) independent of the split point (t0 = 1 vs 1.37) to 45 digits, k = 1..11')
check(all(abs(lam[k] - lam[12 - k]) < mp.mpf(10) ** -45 * abs(lam[k]) for k in range(1, 12)), 'Lambda(k) = Lambda(12 - k)')
r = {k - 1: lam[k] for k in lam}    # r_j = int_0^oo Delta(it) t^j dt = Lambda(j + 1)
ratio_even = [r[0] / r[4], r[2] / r[4]]
ratio_odd = [r[1] / r[5], r[3] / r[5]]
check(abs(ratio_even[0] - mp.mpf(36) / 691 * 70) < mp.mpf(10) ** -40 and abs(ratio_even[1] - mp.mpf(70) / 45) < mp.mpf(10) ** -40,
      f'r_0 : r_2 : r_4 = 36/691 : 1/45 : 1/70 (Zagier r^+), i.e. r_0/r_4 = {mp.nstr(ratio_even[0], 20)} = 2520/691, r_2/r_4 = 14/9')
check(abs(ratio_odd[0] - mp.mpf(12) / 5) < mp.mpf(10) ** -40 and abs(ratio_odd[1] - mp.mpf(5) / 4) < mp.mpf(10) ** -40,
      f'r_1 : r_3 : r_5 = 2/5 : 5/24 : 1/6 (Zagier r^-), i.e. r_1/r_5 = 12/5, r_3/r_5 = 5/4')
L5 = Lval(5)
L6 = Lval(6)
L13 = Lval(13)
direct13 = mp.fsum(mp.mpf(tau[n]) / mp.mpf(n) ** 13 for n in range(1, NMAX + 1))
check(abs(L13 - direct13) < mp.mpf(10) ** -20, f'L(Delta, 13) = {mp.nstr(L13, 25)} (Mellin) = direct Dirichlet sum to n = 3000 (|diff| = {mp.nstr(abs(L13 - direct13), 3)}, tail < 1e-20)')
print(f'     L(Delta, 5) = {mp.nstr(L5, 30)}   L(Delta, 6) (central) = {mp.nstr(L6, 30)}')
check(abs(L5 - mp.mpf('0.666709188434003643')) < mp.mpf(10) ** -17, 'session value L(Delta, 5) = 0.666709188434003643 confirmed')
print(f'     [{time.time() - T0:.0f}s]')

# C4/C5 ---------------------------------------------------------------------------------------------------------
print('C4/C5  Epstein zeta functions by Riemann\'s theta method (no L-functions)')
NA = 120


def epstein_theta(coef, dimhalf, w):
    """Z(2w) = sum' |v|^(-2w) for an even unimodular lattice of dim 2*dimhalf with theta coefficients coef[n] (norm 2n):
    Z(2w) = pi^w/Gamma(w) [sum_n coef[n] (G(w) + G(dimhalf - w)) - 1/w - 1/(dimhalf - w)],
    G(a) = int_1^oo e^(-2 pi n t) t^(a-1) dt = Gamma(a, 2 pi n)/(2 pi n)^a."""
    w = mp.mpf(w)
    tot = mp.mpf(0)
    for n in range(1, NA + 1):
        if coef[n]:
            x = 2 * mp.pi * n
            tot += coef[n] * (mp.gammainc(w, x) / x ** w + mp.gammainc(dimhalf - w, x) / x ** (dimhalf - w))
    return mp.pi ** w / mp.gamma(w) * (tot - 1 / w - 1 / (dimhalf - w))


z10 = epstein_theta(leech, 12, 5)
z10_closed = -mp.mpf(65520) / 691 * mp.mpf(2) ** -5 * L5
check(abs(z10 - z10_closed) < mp.mpf(10) ** -40,
      f'zeta_Leech(10) (theta continuation, Leech coefficients only) = {mp.nstr(z10, 25)} = -(65520/691) 2^-5 L(Delta, 5) = {mp.nstr(z10_closed, 25)}')
check(abs(z10 - mp.mpf('-1.97552397007037')) < 1e-13, 'session value -1.97552397007037 confirmed')
z0 = epstein_theta(leech, 12, mp.mpf('1e-30'))
print(f'     sanity: Z(s) -> -1 as s -> 0: Z(2e-30) = {mp.nstr(z0, 12)}')
z26_theta = epstein_theta(leech, 12, 13)
z26_closed = mp.mpf(65520) / 691 * mp.mpf(2) ** -13 * (mp.zeta(13) * mp.zeta(2) - L13)
# direct: sum_{n <= N} a(n) (2n)^-13 + tail; tail = (65520/691) 2^-13 [sum_{n>N} sigma_11(n) n^-13 - sum_{n>N} tau(n) n^-13]
part = mp.fsum(mp.mpf(leech[n]) / (2 * mp.mpf(n)) ** 13 for n in range(1, NMAX + 1))
# sum_{n>N} sigma_11(n)/n^13 = sum_d d^-13 * sum_{m > N/d} m^-2 = sum_d d^-13 psi'(floor(N/d) + 1)
tail_sig = mp.fsum(mp.mpf(d) ** -13 * mp.psi(1, NMAX // d + 1) for d in range(1, NMAX + 1)) + mp.zeta(13, NMAX + 1) * mp.pi ** 2 / 6
cusp_bound = mp.fsum(mp.mpf(2) * mp.sqrt(n) * mp.mpf(n) ** mp.mpf('-7.5') for n in range(NMAX + 1, 20 * NMAX)) + mp.mpf(10) ** -25
direct = part + mp.mpf(65520) / 691 * mp.mpf(2) ** -13 * tail_sig
err = mp.mpf(65520) / 691 * mp.mpf(2) ** -13 * cusp_bound
check(abs(z26_theta - z26_closed) < mp.mpf(10) ** -40 and abs(direct - z26_closed) <= err,
      f'zeta_Leech(26): closed form {mp.nstr(z26_closed, 22)} = theta method {mp.nstr(z26_theta, 22)}; direct sum (n <= 3000, exact sigma tail) {mp.nstr(direct, 22)}, |diff| = {mp.nstr(abs(direct - z26_closed), 3)} <= cusp-tail bound {mp.nstr(err, 3)}')
check(abs(z26_closed - mp.mpf('0.00749948328893')) < 1e-14 and abs(L13 - mp.mpf('0.99720986806')) < 1e-11,
      'session values zeta_Leech(26) = 0.00749948328893 and L(Delta, 13) = 0.99720986806 confirmed')
e8 = [1] + [240 * s3[n] for n in range(1, 61)] + [240 * s3b[n] for n in range(61, NA + 1)]
ze8 = epstein_theta(e8, 4, 5)
check(abs(ze8 - mp.mpf(5) / 4 * mp.pi ** 2 * mp.zeta(5)) < mp.mpf(10) ** -40,
      f'zeta_E8(10) = (5/4) pi^2 zeta(5) = {mp.nstr(ze8, 15)} by the theta method')
print(f'     [{time.time() - T0:.0f}s]')

# C6 -------------------------------------------------------------------------------------------------------------
print('C6  McLaughlin slice')
found = set()
for a2 in range(2, 60, 2):
    for c2 in range(a2, 200, 2):
        for b in range(-a2, a2 + 1):
            if a2 * c2 - b * b != 15:
                continue
            mn = min(a2 * x * x + 2 * b * x * y + c2 * y * y for x in range(-6, 7) for y in range(-6, 7) if (x, y) != (0, 0))
            if mn >= 4:
                found.add((a2, b, c2))


def gauss_reduce(A, B, C):
    # Gram [[A, B], [B, C]] -> reduced (|2B| <= A <= C) by Lagrange-Gauss steps
    while True:
        if A > C:
            A, C = C, A
            continue
        if abs(2 * B) > A:
            k = round(B / A)
            B, C = B - k * A, C - 2 * k * B + k * k * A
            continue
        return (A, abs(B), C)


red = {gauss_reduce(*f) for f in found}
check(red == {(4, 1, 4)} and len(found) >= 2, f'Gram matrices [[2a, b], [b, 2c]] with det 15 and min >= 4 (2a <= 2c <= 198, |b| <= 2a): {len(found)} found, all Gauss-reduce to [[4, +-1], [+-1, 4]]')


def is_fund(D):
    if D % 4 == 1:
        m = -D
        return all(m % (p * p) for p in range(2, math.isqrt(m) + 1))
    if D % 4 == 0:
        m = D // 4
        return m % 4 in (2, 3) and all((-m) % (p * p) for p in range(2, math.isqrt(-m) + 1))
    return False


hits = []
for D in range(-4000, 0):
    if not is_fund(D):
        continue
    for b in range(-2, 3):
        if b * b - 16 == D:
            hits.append((D, b))
check(sorted(hits) == [(-15, -1), (-15, 1)], f'fundamental D in (-4000, 0) with a reduced form (2, b, 2): {sorted(hits)} (b = 0, +-2 give D = -16, -12, not fundamental): Prop. 6.1 true, its proof skips even D')
# Leech: norm-4 x with |v - x|^2 = 4 for a fixed norm-6 v (in raw coordinates norm 8*N)
# extended binary Golay code from the generator polynomial of the [23,12] cyclic code, extended by parity
g = [1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1]   # 1 + x^2 + x^4 + x^5 + x^6 + x^10 + x^11 (coefficients low to high)
words = set()
gens = []
for sh in range(12):
    w = [0] * 23
    for i, c in enumerate(g):
        if c:
            w[(i + sh) % 23] ^= 1
    gens.append(w + [sum(w) % 2])
for mask in range(1 << 12):
    w = [0] * 24
    for i in range(12):
        if (mask >> i) & 1:
            w = [a ^ b for a, b in zip(w, gens[i])]
    words.add(tuple(w))
wt = {}
for w in words:
    wt[sum(w)] = wt.get(sum(w), 0) + 1
check(len(words) == 4096 and wt == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}, f'extended Golay code: weight distribution {dict(sorted(wt.items()))}')
octads = [w for w in words if sum(w) == 8]
minv = []
for o in octads:
    S = [i for i in range(24) if o[i]]
    for signs in itertools.product((1, -1), repeat=7):
        sg = list(signs) + [1 if signs.count(-1) % 2 == 0 else -1]
        v = [0] * 24
        for i, s in zip(S, sg):
            v[i] = 2 * s
        minv.append(tuple(v))
# shape (-+3, +-1^23): start from (-3 at i, +1 elsewhere) and flip the signs on a Golay codeword
for w in words:
    for i in range(24):
        minv.append(tuple((-1) ** c * (-3 if j == i else 1) for j, c in enumerate(w)))
for i, j in itertools.combinations(range(24), 2):
    for si in (4, -4):
        for sj in (4, -4):
            v = [0] * 24
            v[i], v[j] = si, sj
            minv.append(tuple(v))
check(len(minv) == 196560 and len(set(minv)) == 196560 and all(sum(t * t for t in v) == 32 for v in minv),
      f'Leech minimal vectors constructed: {len(set(minv))} distinct, raw norm 32 (= 4 after /8)')
# a norm-6 vector: (2^12 0^12) on a dodecad with an even (0 mod 4) pattern: take x = 2 c for a dodecad c, adjust signs
dodec = next(w for w in words if sum(w) == 12)
v = [2 * c for c in dodec]
s = sum(v)
# require sum = 0 mod 8: 24 = 0 mod 8 ok; check membership: v.m must be = 0 mod 8 for all minimal m (unimodular, raw/8 integral)
ok_member = all(sum(a * b for a, b in zip(v, m)) % 8 == 0 for m in minv[:5000])
cntx = sum(1 for m in minv if sum(a * b for a, b in zip(v, m)) == 24)
check(sum(t * t for t in v) == 48 and ok_member and cntx == 552,
      f'v = 2 * (dodecad) has norm 6; minimal x with x.v = 3 (|v - x|^2 = 4): {cntx} = 2 * 276 (Co_3 on 276 pairs {{x, v - x}}, stabiliser McL:2)')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')

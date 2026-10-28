"""opus-2026-10-28-S67: the owner's Borwein ladder and the twelve behind {2,3,4,5,7,13}, {2,3,11}, {2,5,7}.

Sections
  A  odd harmonic crossings O_m = 1 + 1/3 + ... + 1/(2m+1): thresholds 2, 3, 4, 5, 6, 7 and 2*pi (interval arithmetic)
  B  Borwein integrals: exact first-failure values (rational multiples of pi) and the mantissas 2.31, 2.33, 4.21
  C  the crossing formula m(n) = floor(exp(2n - gamma)/4): rigorous window test for 2 <= n <= NMAX (HYP-9436)
  D  von Staudt-Clausen as orbifold Euler characteristics; triangle groups (2,3,p); genus of X(p)
  E  sigma = 144: the two prime triples {2,3,11}, {2,5,7}; the vSC sinc products (Borwein excess = -B_2k)
  F  the owner's {231, 233, 420} in prior repository work (Collatz 231 -> 233; square-sigma aliquot targets)

Every assertion is an `assert`; run with python3 (and python3 -O for the receipts-only pass).
"""
import sys
import time
from fractions import Fraction as Fr
from math import comb, factorial, isqrt, gcd

import mpmath as mp
from sympy import bernoulli, isprime, divisor_sigma, factorint

T0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000


def say(*a):
    print(*a, flush=True)


def O_exact(m):
    return sum(Fr(1, 2 * j + 1) for j in range(m + 1))


# ---------------------------------------------------------------- A
say("=" * 78)
say("A. Odd harmonic crossings: m(T) = least m with O_m = sum_{j<=m} 1/(2j+1) > T")
say("=" * 78)
mp.mp.dps = 40
iv = mp.iv
iv.dps = 40
LO = lambda x: mp.make_mpf(x._mpi_[0])
HI = lambda x: mp.make_mpf(x._mpi_[1])


def crossing(T_lo, T_hi, label):
    """least m with O_m > T, certified with interval arithmetic: T in [T_lo, T_hi]."""
    s_f = 0.0
    m = -1
    # fast float scan to a neighbourhood, then certify with intervals
    while True:
        m += 1
        s_f += 1.0 / (2 * m + 1)
        if s_f > float(T_hi) - 1e-9:
            break
    m0 = max(0, m - 3)
    S = iv.mpf(0)
    for j in range(m0 + 1):
        S += iv.mpf(1) / (2 * j + 1)
    k = m0
    while True:
        if S.a > T_hi:          # certainly O_k > T
            break
        if S.b > T_lo:          # undecided
            raise RuntimeError("interval too wide at " + label)
        k += 1
        S += iv.mpf(1) / (2 * k + 1)
    # certify O_{k-1} <= T
    Sprev = S - iv.mpf(1) / (2 * k + 1)
    assert Sprev.b < T_lo
    return k, S, Sprev


table = {}
for T in [2, 3, 4, 5, 6, 7]:
    k, S, Sp = crossing(mp.mpf(T), mp.mpf(T), str(T))
    table[T] = k
    say("  T = %-6s m(T) = %-7d denominator 2m+1 = %-7d O_m - T in [%s, %s]"
        % (T, k, 2 * k + 1, mp.nstr(LO(S) - T, 8), mp.nstr(HI(S) - T, 8)))
assert [table[T] for T in [2, 3, 4, 5, 6, 7]] == [7, 56, 418, 3091, 22845, 168803]
twopi = iv.pi * 2
k2pi, S2pi, Sp2pi = crossing(twopi.a, twopi.b, "2pi")
say("  T = 2*pi  m(T) = %d denominator = %d ; O_m - 2pi in [%s, %s] ; O_(m-1) - 2pi <= %s"
    % (k2pi, 2 * k2pi + 1, mp.nstr(LO(S2pi) - HI(twopi), 12), mp.nstr(HI(S2pi) - LO(twopi), 12),
       mp.nstr(HI(Sp2pi) - LO(twopi), 6)))
assert k2pi == 40249 and 2 * k2pi + 1 == 80499
say("  => owner's 7, 56, 418 = m(2), m(3), m(4); owner's 15, 113 = 2m+1 at T = 2, 3; owner's 80499 = 2m+1 at T = 2*pi.")
say("     (the T = 4 denominator is 837 = 2*418 + 1; 80499 is NOT on the integer ladder 15, 113, 837, 6183, 45691, 337607)")

# ---------------------------------------------------------------- B
say("")
say("=" * 78)
say("B. Borwein integrals I_m(t) = int_0^oo (sin tx / x) prod_{j=1}^m sinc(x/(2j+1)) dx = pi/2 - pi*P(S'_m > t)")
say("   F_m(xi) = int_R prod_{j=0}^m sinc(x/(2j+1)) cos(xi x) dx = pi*P(xi-1 <= S'_m <= xi+1)")
say("=" * 78)


def tail_prob(a, t):
    """exact P(sum U_j > t), U_j uniform on [-a_j, a_j] (Fractions), inclusion-exclusion over corners."""
    m = len(a)
    A = sum(a)
    c = A - t
    if c <= 0:
        return Fr(0)
    # enumerate subsets with positive remainder only (pruned DFS on sorted widths)
    w = sorted((2 * x for x in a))
    tot = Fr(0)
    stack = [(0, Fr(0), 0)]
    while stack:
        i, used, size = stack.pop()
        if i == len(w):
            r = c - used
            if r > 0:
                tot += (-1) ** size * r ** m
            continue
        stack.append((i + 1, used, size))
        if c - used - w[i] > 0:
            stack.append((i + 1, used + w[i], size + 1))
    den = factorial(m)
    for x in a:
        den *= 2 * x
    return tot / den


# B1. the classical t = 1 case and the published value
a7 = [Fr(1, 2 * j + 1) for j in range(1, 8)]
for m in range(1, 7):
    assert tail_prob(a7[:m], 1) == 0
P7 = tail_prob(a7, 1)
I7 = Fr(1, 2) - P7          # I_7(1) / pi
say("  I_7(1)/pi = %s/%s" % (I7.numerator, I7.denominator))
assert I7 == Fr(467807924713440738696537864469, 935615849440640907310521750000)   # Borwein-Borwein (2001)
say("  matches the published Borwein value 467807924713440738696537864469/935615849440640907310521750000 * pi")


def first_failure_value(m, xi_lo, xi_hi=None):
    """F_m(xi) at the first failure via the closed form pi*eps^m*(2m+1)!!/(2m)!!, eps = O_m - xi (interval)."""
    S = iv.mpf(0)
    for j in range(m + 1):
        S += iv.mpf(1) / (2 * j + 1)
    xi = iv.mpf([xi_lo, xi_hi if xi_hi is not None else xi_lo])
    eps = S - xi
    assert eps.a > 0 and eps.b <= iv.mpf(1) / (2 * m + 1)
    # log10 of (2m+1)!!/(2m)!! = (2m+1)!/(2^m m!)^2
    lg = iv.log10(iv.pi) + m * iv.log10(eps) + iv.log10(iv.factorial(2 * m + 1)) \
        - 2 * (m * iv.log10(2) + iv.log10(iv.factorial(m)))
    e = int(mp.floor(lg.a))
    man = iv.power(10, lg - e)
    return man, e, eps


rows = []
for (m, xi, lab) in [(7, 2, "t=1 (classic, owner 15)"), (56, 3, "t=2 (2cos x variant, owner 113)"),
                     (418, 4, "t=3 (837)"), (40249, None, "xi=2pi (owner 80499)")]:
    if xi is None:
        man, e, eps = first_failure_value(m, twopi.a, twopi.b)
    else:
        man, e, eps = first_failure_value(m, mp.mpf(xi))
        # exact rational cross-check for the integer thresholds
        P = tail_prob([Fr(1, 2 * j + 1) for j in range(1, m + 1)], xi - 1)
        eps_ex = O_exact(m) - xi
        closed = eps_ex ** m
        for j in range(1, m + 1):
            closed *= (2 * j + 1)
        closed /= factorial(m) * 2 ** m
        assert P == closed, "single-corner closed form differs from inclusion-exclusion"
        val = mp.mpf(P.numerator) / P.denominator * mp.pi
        assert LO(man) <= val / mp.mpf(10) ** e <= HI(man)
    rows.append((lab, m, 2 * m + 1, man, e))
    say("  %-34s m = %-6d D = %-6d first-failure value = [%s, %s] x 10^%d"
        % (lab, m, 2 * m + 1, mp.nstr(LO(man), 11), mp.nstr(HI(man), 11), e))
mans = [r[3] for r in rows]
assert mp.nstr(LO(mans[0]), 3) == "2.31" == mp.nstr(HI(mans[0]), 3) and rows[0][4] == -11
assert mp.nstr(LO(mans[1]), 3) == "2.33" == mp.nstr(HI(mans[1]), 3) and rows[1][4] == -138
assert mp.nstr(LO(mans[3]), 3) == "4.21" == mp.nstr(HI(mans[3]), 3) and rows[3][4] == -226577
say("  => owner's 231 and 233 are the 3-digit mantissas at the breaks 15 and 113;")
say("     the 2*pi break (80499) has mantissa 4.2102 (owner's 420 agrees to two digits).")

# B2. how a double-precision computation reads the third mantissa
fw = 0.0
for j in range(40249 + 1):
    fw += 1.0 / (2 * j + 1)
eps_f = mp.mpf(fw) - 2 * mp.pi
lgf = mp.log10(mp.pi) + 40249 * mp.log10(eps_f) + mp.log10(mp.factorial(80499)) \
    - 2 * (40249 * mp.log10(2) + mp.log10(mp.factorial(40249)))
manf = mp.power(10, lgf - mp.floor(lgf))
say("  naive forward double-precision summation gives eps = %s and mantissa %s (truncates to 4.20)"
    % (mp.nstr(eps_f, 12), mp.nstr(manf, 6)))
assert 4.20 <= manf < 4.21

# B3. telescoping identity pi/2 - I_m(t) = sum_k F_m(t + 1 + 2k), checked exactly for small m and rational t
for m in range(1, 9):
    a = [Fr(1, 2 * j + 1) for j in range(1, m + 1)]
    for t in [Fr(0), Fr(1, 3), Fr(1, 2), Fr(1), Fr(3, 2)]:
        lhs = tail_prob(a, t)                       # (pi/2 - I_m(t))/pi
        rhs = Fr(0)
        k = 0
        while t + 2 * k <= sum(a):
            # F_m(xi)/pi = P(S' >= xi-1) - P(S' > xi+1), xi = t+1+2k
            rhs += tail_prob(a, t + 2 * k) - tail_prob(a, t + 2 * k + 2)
            k += 1
        assert lhs == rhs
say("  telescoping identity pi/2 - I_m(t) = sum_{k>=0} F_m(t+1+2k): exact for m <= 8, t in {0,1/3,1/2,1,3/2}")

# B4. light quadrature sanity check of the 2.31e-11 deviation (independent of the probabilistic formula)
mp.mp.dps = 30
f = lambda x: mp.sinc(x) * mp.sinc(x / 3) * mp.sinc(x / 5) * mp.sinc(x / 7) * mp.sinc(x / 9) * \
    mp.sinc(x / 11) * mp.sinc(x / 13) * mp.sinc(x / 15)
L = 700
pts = [mp.pi * k for k in range(0, int(L / mp.pi) + 1)]
quad = mp.quad(f, pts)
tail_bound = mp.mpf(2027025) / (7 * mp.mpf(pts[-1]) ** 7)
dev_q = mp.pi / 2 - quad
say("  quadrature on [0, %s]: pi/2 - I_7(1) = %s (tail <= %s); closed form 2.3100572725e-11"
    % (mp.nstr(pts[-1], 6), mp.nstr(dev_q, 8), mp.nstr(tail_bound, 2)))
assert abs(dev_q - mp.mpf("2.3100572725e-11")) < 1e-15 + tail_bound
mp.mp.dps = 40

# ---------------------------------------------------------------- C
say("")
say("=" * 78)
say("C. Crossing formula: m(xi) = floor(X), X = exp(2 xi - gamma)/4, whenever frac(X) >= X/(24 floor(X)^2)")
say("   (O_m = (1/2)ln(4M) + gamma/2 + r(M), M = m+1, 0 < r(M) < 1/(48 M^2) for M >= 1: GKP (9.89) bounds on H_N)")
say("=" * 78)
# r(M) bounds checked directly for small M (the GKP envelope gives them for all M >= 2)
for M in list(range(1, 400)) + [1000, 5000]:
    Om = mp.fsum([mp.mpf(1) / (2 * j + 1) for j in range(M)])
    r = Om - (mp.log(4 * M) + mp.euler) / 2
    assert 0 < r < 1 / (48 * mp.mpf(M) ** 2), M
say("  0 < r(M) < 1/(48 M^2) verified directly for M < 400, M = 1000, 5000")
# consistency with Section A
for T, k in table.items():
    X = mp.exp(2 * T - mp.euler) / 4
    assert int(mp.floor(X)) == k
X2pi = mp.exp(4 * mp.pi - mp.euler) / 4
assert int(mp.floor(X2pi)) == 40249
say("  formula reproduces m(2..7) and m(2*pi) = 40249 (X = %s)" % mp.nstr(X2pi, 10))

worst = None
exceptions = []
fracs = []
tC = time.time()
DMAX = int(0.8686 * NMAX) + 60          # one working precision; gamma computed once
mp.mp.dps = DMAX
G = +mp.euler
for n in range(2, NMAX + 1):
    X = mp.exp(2 * n - G) / 4
    fl = mp.floor(X)
    fr = X - fl
    window = X / (24 * fl ** 2)
    # precision guard: |error in X| < 10^(digits before the point + 5 - DMAX); require a clear margin
    guard = mp.mpf(10) ** (int(0.8686 * n) + 5 - DMAX)
    ratio = fr / window
    fracs.append(float(fr))
    if worst is None or ratio < worst[1]:
        worst = (n, ratio, fr)
    if not (fr - guard > window):
        exceptions.append(n)
mp.mp.dps = 40
say("  n = 2..%d: exceptions to the window test: %s ; tightest n = %d with frac(X)/window = %s (frac = %s)"
    % (NMAX, exceptions, worst[0], mp.nstr(worst[1], 6), mp.nstr(worst[2], 6)))
say("  => m(n) = floor(exp(2n - gamma)/4) PROVED for every 2 <= n <= %d (%.1fs)" % (NMAX, time.time() - tC))
assert exceptions == []
# equidistribution of frac(X_n) (heuristic only): Kolmogorov-Smirnov statistic
fs = sorted(fracs)
N = len(fs)
D = max(max((i + 1) / N - x, x - i / N) for i, x in enumerate(fs))
say("  heuristic: frac(X_n) KS distance from uniform D = %.4f (N = %d; 1%%-critical ~ %.4f)" % (D, N, 1.63 / N ** 0.5))
expected = mp.fsum([1 / (24 * mp.exp(2 * n - mp.euler) / 4) for n in range(2, 60)])
say("  heuristic expected number of exceptions over all n >= 2 (uniform fractional parts): %s" % mp.nstr(expected, 4))

# C2. the owner's 2*pi generalised: thresholds xi = k*pi, 1 <= k <= KMAX (pre-registered HYP-9436 (P3))
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 500
DK = int(0.8686 * 2 * mp.pi * KMAX) + 80
mp.mp.dps = DK
G = +mp.euler
PI = +mp.pi
exc_pi = []
first = []
for k in range(1, KMAX + 1):
    X = mp.exp(2 * k * PI - G) / 4
    fl = mp.floor(X)
    fr = X - fl
    window = X / (24 * fl ** 2)
    guard = mp.mpf(10) ** (int(0.8686 * 2 * float(PI) * k) + 5 - DK)
    if not (fr - guard > window):
        exc_pi.append(k)
    if k <= 3:
        first.append((k, int(fl), 2 * int(fl) + 1))
mp.mp.dps = 40
say("  thresholds xi = k*pi, k = 1..%d: window exceptions %s ; m(k pi), 2m+1 for k = 1, 2, 3: %s" % (KMAX, exc_pi, first))
assert exc_pi == [] and first[1] == (2, 40249, 80499)

# C3. OEIS A092315 lists two formulas: floor(X) (Kirchner 2020) and floor(X + 1/8) (Hasler 2017, 'conjectured').
#     They differ exactly when frac(X) >= 7/8; there the window test proves floor(X), so the +1/8 form fails.
mp.mp.dps = DMAX
G = +mp.euler
hasler_fail = []
for n in range(2, NMAX + 1):
    X = mp.exp(2 * n - G) / 4
    if X - mp.floor(X) >= mp.mpf(7) / 8:
        hasler_fail.append(n)
mp.mp.dps = 80
def O_dg(m):                       # O_m = (psi(m + 3/2) - psi(1/2))/2, an independent evaluation
    return (mp.digamma(m + mp.mpf(3) / 2) - mp.digamma(mp.mpf(1) / 2)) / 2
for m in [0, 1, 7, 56, 418]:
    assert abs(O_dg(m) - mp.fsum([mp.mpf(1) / (2 * j + 1) for j in range(m + 1)])) < mp.mpf(10) ** -70
checks = []
for n in hasler_fail[:3]:
    X = mp.exp(2 * n - mp.euler) / 4
    F = int(mp.floor(X))
    assert O_dg(F - 1) < n < O_dg(F)           # a(n) = floor(X), not floor(X) + 1
    checks.append((n, F, mp.nstr(O_dg(F - 1) - n, 3), mp.nstr(O_dg(F) - n, 3)))
mp.mp.dps = 40
say("  Hasler's floor(X + 1/8) differs from the proved value at %d of the n <= %d; first: %s" % (len(hasler_fail), NMAX, hasler_fail[:8]))
say("  digamma cross-check (n, a(n), O_(a-1) - n, O_a - n): %s" % checks)
assert hasler_fail[0] == 17 and checks[0][1] == 81897532160124

# ---------------------------------------------------------------- D
say("")
say("=" * 78)
say("D. von Staudt-Clausen as orbifold Euler characteristics")
say("=" * 78)


def vsc(n2k):
    return [p for p in range(2, n2k + 2) if isprime(p) and n2k % (p - 1) == 0]


def B(n):
    b = bernoulli(n)
    return Fr(int(b.p), int(b.q))


exact_k = []
for k in range(1, 61):
    v = vsc(2 * k)
    chi = 2 - len(v) + sum(Fr(1, p) for p in v)          # sphere with cone points v
    Nk = sum(Fr(1, p) for p in v) + B(2 * k)
    assert Nk.denominator == 1                            # von Staudt-Clausen
    assert (chi + B(2 * k)).denominator == 1              # chi = -B_2k mod 1
    if chi == -B(2 * k):
        exact_k.append(2 * k)
    if k <= 7:
        say("  2k = %-3d vSC = %-17s chi_orb = %-12s -B_2k = %-12s N_k = %s"
            % (2 * k, v, chi, -B(2 * k), Nk))
say("  chi_orb(S^2(vSC(2k))) = -B_2k (mod 1) for all 2k <= 120; exact equality only at 2k in %s" % exact_k)
assert exact_k == [4, 6, 8, 10]
# |B_2k| > 1 for 2k >= 14 and |chi| <= |vSC| - 2 + ... => no exact equality beyond (proof in the note); check |B| growth
for k in range(7, 61):
    assert abs(B(2 * k)) > 1

tri_exact, tri_mod = [], []
for p in range(5, 400):
    if not isprime(p):
        continue
    chi = Fr(1, p) - Fr(1, 6)
    if chi == -B(p - 1):
        tri_exact.append(p)
    if (chi + B(p - 1)).denominator == 1:
        tri_mod.append(p)
        assert vsc(p - 1) == [2, 3, p]
    else:
        assert vsc(p - 1) != [2, 3, p]
say("  triangle groups (2,3,p): chi = 1/p - 1/6 = -B_(p-1) exactly for p in %s" % tri_exact)
say("  chi = -B_(p-1) mod 1 iff vSC(p-1) = {2,3,p}: p in %s" % tri_mod[:14])
assert tri_exact == [5, 7, 11]
for p, g in [(5, 0), (7, 3), (11, 26)]:
    G = p * (p * p - 1) // 2
    assert 2 * g - 2 == G * B(p - 1)
    assert g == 1 + (p * p - 1) * (p - 6) // 24
    say("  X(%d): 2g - 2 = |PSL2(F_%d)| * B_%d = %d * %s = %d  (g = %d)" % (p, p, p - 1, G, B(p - 1), 2 * g - 2, g))

# ---------------------------------------------------------------- E
say("")
say("=" * 78)
say("E. The twelve: {2,3,4,5,7,13}, {2,3,11}, {2,5,7}")
say("=" * 78)
d12 = sorted(d + 1 for d in range(1, 13) if 12 % d == 0)
assert d12 == [2, 3, 4, 5, 7, 13] and set(d12) == {4, 2, 5, 3, 7, 13}
assert [p for p in d12 if isprime(p)] == vsc(12) == [2, 3, 5, 7, 13]
assert B(12) == Fr(-691, 2730) and B(10) == Fr(5, 66) and vsc(10) == [2, 3, 11]
say("  {d+1 : d | 12} = %s ; primes = vSC(12) = %s ; B_12 = -691/2730 ; vSC(10) = [2, 3, 11], B_10 = 5/66" % (d12, vsc(12)))
# sigma^{-1}(144) and prime sets with prod (p+1) = 144
pre = [n for n in range(1, 145) if divisor_sigma(n) == 144]
assert pre == [66, 70, 94, 115, 119]
primes = [p for p in range(2, 144) if isprime(p) and 144 % (p + 1) == 0]
sets = []


def rec(i, prod, chosen):
    if prod == 144:
        sets.append(list(chosen))
        return
    for j in range(i, len(primes)):
        q = primes[j] + 1
        if 144 % (prod * q) == 0:
            rec(j + 1, prod * q, chosen + [primes[j]])


rec(0, 1, [])
say("  sigma^{-1}(144) = %s ; prime sets with prod(p+1) = 144: %s" % (pre, sets))
assert sorted(map(tuple, sets)) == [(2, 3, 11), (2, 5, 7), (2, 47), (5, 23), (7, 17)]
assert [s for s in sets if len(s) == 3] in ([[2, 3, 11], [2, 5, 7]], [[2, 5, 7], [2, 3, 11]])
subs = [s for s in sets if set(s) <= set(vsc(12))]
assert subs == [[2, 5, 7]]
say("  the only 3-prime solutions are {2,3,11} = vSC(10) and {2,5,7}, the unique subset of vSC(12) with sigma = 144")
say("  triangle groups: chi(2,3,11) = %s = -B_10 ; chi(2,5,7) = %s" % (Fr(1, 2) + Fr(1, 3) + Fr(1, 11) - 1,
                                                                       Fr(1, 2) + Fr(1, 5) + Fr(1, 7) - 1))
# the vSC sinc products: int_0^oo sinc(x) prod_{p in vSC(2k)} sinc(x/p) dx = pi/2 - pi*P(sum U_p > 1)
for k in range(1, 7):
    v = vsc(2 * k)
    a = [Fr(1, p) for p in v]
    P = tail_prob(a, 1)
    excess = sum(a) - 1
    assert excess == -B(2 * k)                     # N_k = 1 for 2k <= 12
    say("  2k = %-3d vSC = %-17s excess sum(1/p) - 1 = %-10s deficit (pi/2 - I)/pi = %s"
        % (2 * k, v, excess, P))
    if 2 * k in (4, 8):
        assert P == Fr(1, 43200)
    if 2 * k == 12:
        assert P == Fr(691 ** 5 - 271 ** 5, 3840 * 2730 ** 4)
        assert 691 - 271 == 420 == 2 * 2730 // 13
say("  => for 2k <= 12 the vSC sinc product is a Borwein integral with excess -B_2k: exact (pi/2) for 2k = 2, 6, 10;")
say("     deficit pi/43200 for 2k = 4, 8; deficit pi*(691^5 - 271^5)/(3840*2730^4) for 2k = 12 (271 = 691 - 420, 420 = 2*2730/13)")
# J-homomorphism / Eisenstein constants (classical)
for k, want in [(5, 264), (6, 65520)]:
    val = B(2 * k) / (4 * k)
    assert abs(val.denominator) == want
say("  denominators of B_10/20 and B_12/24 (= |im J| in pi_19^s and pi_23^s): 264 = 2^3*3*11 and 65520 = 2^4*3^2*5*7*13")

# E2. the twelve behind the labels: (p-1) | 12, (p+1) | 12 and (p+1) | 24 lists; genus of X_0(p); Monster exponents
P200 = [p for p in range(2, 200) if isprime(p)]
def genus_X0(p):
    if p == 2 or p == 3:
        return 0
    nu2 = 1 + (1 if p % 4 == 1 else -1)
    nu3 = 1 + (1 if p % 3 == 1 else -1)
    g = Fr(p + 1, 12) - Fr(nu2, 4) - Fr(nu3, 3)
    assert g.denominator == 1
    return int(g)
gz = [p for p in P200 if genus_X0(p) == 0]
m12 = [p for p in P200 if 12 % (p - 1) == 0]
p12 = [p for p in P200 if 12 % (p + 1) == 0]
p24 = [p for p in P200 if 24 % (p + 1) == 0]
assert gz == m12 == vsc(12) == [2, 3, 5, 7, 13]
assert p12 == [2, 3, 5, 11] and p24 == [2, 3, 5, 7, 11, 23]
assert all(12 % (q + 1) == 0 for q in (2, 3, 11)) and all(12 % (q - 1) == 0 for q in (2, 5, 7))
say("  genus X_0(p) = 0 exactly for p in %s = {p : (p-1) | 12} = vSC(12)" % gz)
say("  {p : (p+1) | 12} = %s contains {2,3,11}; {p : (p-1) | 12} contains {2,5,7}; {p : (p+1) | 24} = %s (M24 eta products)" % (p12, p24))
MON = {2: 46, 3: 20, 5: 9, 7: 6, 11: 2, 13: 3, 17: 1, 19: 1, 23: 1, 29: 1, 31: 1, 41: 1, 47: 1, 59: 1, 71: 1}
BAB = {2: 41, 3: 13, 5: 6, 7: 2, 11: 1, 13: 1, 17: 1, 19: 1, 23: 1, 31: 1, 47: 1}
cube = sorted(p for p, e in MON.items() if e >= 3)
ratio = {p: MON[p] - BAB.get(p, 0) for p in MON if MON[p] - BAB.get(p, 0) > 0}
sq_ratio = sorted(p for p, e in ratio.items() if e >= 2)
assert cube == [2, 3, 5, 7, 13] == sq_ratio
assert {p: ratio[p] for p in (2, 3, 5, 7, 11, 13)} == {2: 5, 3: 7, 5: 3, 7: 4, 11: 1, 13: 2}
assert all(MON[p] >= 1 + 24 // (p - 1) for p in cube)          # room for the extraspecial p^(1+24/(p-1))
say("  Monster: p^3 | |M| exactly for p in %s; primes with exponent >= 2 in |M|/|B| = %s; exponents 2^5 3^7 5^3 7^4 11 13^2 (Session43 cycle)" % (cube, sq_ratio))

# ---------------------------------------------------------------- F
say("")
say("=" * 78)
say("F. {231, 233, 420} in prior repository work")
say("=" * 78)


def syr(n):
    n = 3 * n + 1
    while n % 2 == 0:
        n //= 2
    return n


w = [231]
for _ in range(17):
    w.append(syr(w[-1]))
assert w == [231, 347, 521, 391, 587, 881, 661, 31, 47, 71, 107, 161, 121, 91, 137, 103, 155, 233]
for t in range(1, 4):
    x = 231 + (1 << 28) * t
    for _ in range(17):
        x = syr(x)
    assert x == 233 + 258280326 * t and 258280326 == 2 * 3 ** 17
say("  Collatz: 17 Syracuse steps 231 -> 233 (crossroads233_20260926_graph.md); m_17(231 + 2^28 t) = 233 + 2*3^17 t checked t = 1..3")


def s(n):
    return int(divisor_sigma(n)) - n


assert s(345) == 231 and isqrt(int(divisor_sigma(345))) ** 2 == 576 == int(divisor_sigma(345))
assert s(364) == 420 and int(divisor_sigma(364)) == 784 == 28 ** 2
# 233 is not s(n) for any n with square sigma (s(n) >= sqrt(n) + 1 for composite n, so n < 233^2)
hits = [n for n in range(2, 233 ** 2) if s(n) == 233 and isqrt(int(divisor_sigma(n))) ** 2 == int(divisor_sigma(n))]
say("  aliquot: 231 = s(345), sigma(345) = 24^2 ; 420 = s(364), sigma(364) = 28^2 ; square-sigma preimages of 233: %s" % hits)
assert hits == []
say("")
say("ALL ASSERTIONS PASSED (%.1fs)" % (time.time() - T0))

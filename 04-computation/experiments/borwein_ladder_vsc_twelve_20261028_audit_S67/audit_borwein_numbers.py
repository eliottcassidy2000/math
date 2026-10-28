"""Independent referee re-computation (audit of opus-2026-10-28-S67) -- Borwein ladder numerics.

Written from scratch; it does not import or copy the session script.

Methods (all independent of the session code):
  * O_m = sum_{j=0}^{m} 1/(2j+1) as an EXACT rational by binary splitting (python ints);
  * pi enclosed by EXACT rational bounds from a hard-coded 100-digit decimal (checked against mpmath);
  * first-failure values from the Borwein-Borwein (2001) Corollary 1 / Wikipedia "general formula"
    (sign-sum over gamma in {+-1}^n), not from the session's probabilistic derivation;
  * crossing-index window test with mpmath INTERVAL arithmetic (iv), not floating point with a guard;
  * the n = 17 OEIS check by the DLMF 5.11(ii) enveloping bound for psi (rigorous) AND mpmath digamma;
  * m(3 pi) by math.fsum over correctly-rounded doubles with an explicit error bound AND by the envelope.

Run:  python audit_borwein_numbers.py [NMAX=3000] [KMAX=500]
"""
import sys
import time
import math
from fractions import Fraction as Fr
from itertools import product as iproduct

import mpmath as mp

T0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 500
iv = mp.iv
FAIL = []
LO = lambda x: mp.make_mpf(x._mpi_[0])       # exact lower endpoint (mpf) of an mpmath interval
HI = lambda x: mp.make_mpf(x._mpi_[1])       # exact upper endpoint


def say(*a):
    print(*a, flush=True)


def check(cond, label):
    say(("  [PASS] " if cond else "  [FAIL] ") + label)
    if not cond:
        FAIL.append(label)


# ------------------------------------------------------------------ exact helpers
def bs(a, b):
    """sum_{j=a}^{b-1} 1/(2j+1) = P/Q exactly, Q = prod (2j+1)."""
    if b - a == 1:
        return 1, 2 * a + 1
    mid = (a + b) // 2
    p1, q1 = bs(a, mid)
    p2, q2 = bs(mid, b)
    return p1 * q2 + p2 * q1, q1 * q2


def O_pq(m):
    return bs(0, m + 1)


PI_STR = ("3.14159265358979323846264338327950288419716939937510"
          "58209749445923078164062862089986280348253421170679")   # 100 decimals
mp.mp.dps = 120
assert mp.nstr(mp.pi, 115).replace(".", "")[:101] == PI_STR.replace(".", "")[:101]   # truncation, not rounding
PI_LO = Fr(int(PI_STR.replace(".", "")), 10 ** 100)
PI_HI = PI_LO + Fr(1, 10 ** 100)


def ratio_rel(p, q, T):
    """(p/q - T) as mpf (T rational)."""
    T = Fr(T)
    num = p * T.denominator - T.numerator * q
    return mp.mpf(num) / (mp.mpf(q) * T.denominator)


def float_scan(T):
    s, m = 0.0, -1
    while True:
        m += 1
        s += 1.0 / (2 * m + 1)
        if s > T:
            return m


# ======================================================================== 1
say("=" * 100)
say("1. Odd harmonic crossings m(T) = least m with O_m > T, exact rationals (binary splitting)")
say("=" * 100)
mp.mp.dps = 40
crossings = {}
for T in [1, 2, 3, 4, 5, 6, 7]:
    m0 = float_scan(T)
    while True:
        p, q = O_pq(m0)
        above = p > T * q                                   # O_m0 > T ?
        prev_below = (p - T * q) * (2 * m0 + 1) <= q         # O_{m0-1} <= T ? (O_0 = 1 = T at T = 1)
        if above and prev_below:
            break
        m0 = m0 + 1 if not above else m0 - 1
    crossings[T] = m0
    say("  T = %d : m = %-7d 2m+1 = %-7d O_m - T = %s ; O_(m-1) - T = %s"
        % (T, m0, 2 * m0 + 1, mp.nstr(ratio_rel(p, q, T), 12),
           mp.nstr(ratio_rel(p, q, T) - mp.mpf(1) / (2 * m0 + 1), 12)))
check([crossings[T] for T in range(2, 8)] == [7, 56, 418, 3091, 22845, 168803],
      "m(2..7) = 7, 56, 418, 3091, 22845, 168803 (exact rationals)")
check(crossings[1] == 1, "m(1) = 1")

# k*pi thresholds with exact rational pi bounds: m(pi), m(2pi)
for k, claim in [(1, 75), (2, 40249)]:
    m = claim
    p, q = O_pq(m)
    lo, hi = k * PI_LO, k * PI_HI
    above = Fr(p, 1) > hi * q                                   # O_m > k*pi_hi
    pm1 = Fr(p * (2 * m + 1) - q, q * (2 * m + 1))              # O_{m-1}
    below = pm1 < lo
    mp.mp.dps = 50
    d_hi = mp.mpf(p) / q - k * mp.mpf(PI_LO.numerator) / PI_LO.denominator
    d_lo = mp.mpf(pm1.numerator) / pm1.denominator - k * mp.mpf(PI_LO.numerator) / PI_LO.denominator
    say("  T = %d*pi : m = %d, 2m+1 = %d ; O_m - T = %s ; O_(m-1) - T = %s"
        % (k, m, 2 * m + 1, mp.nstr(d_hi, 12), mp.nstr(d_lo, 12)))
    check(above and below, "m(%d pi) = %d certified by exact rationals and rational pi bounds (100 digits)" % (k, claim))
    if k == 2:
        P2PI, Q2PI = p, q

# ======================================================================== 2
say("")
say("=" * 100)
say("2. First-failure values  pi * eps^m * (2m+1)!!/(2m)!!  (BB 2001 Cor. 1 / BBB 2008 Thm 2 'First Bite')")
say("=" * 100)


def log10_ratio_dfact(m, exact=True):
    """log10((2m+1)!!/(2m)!!)."""
    if exact:
        num = 1
        for j in range(1, m + 1):
            num *= 2 * j + 1
        den = math.factorial(m) * (1 << m)
        return mp.log10(mp.mpf(num)) - mp.log10(mp.mpf(den))
    return (mp.loggamma(2 * m + 2) - 2 * mp.loggamma(m + 1) - 2 * m * mp.log(2)) / mp.log(10)


def man_exp(L):
    e = int(mp.floor(L))
    return mp.power(10, L - e), e


mp.mp.dps = 80
results = {}
for T in [2, 3, 4]:
    m = crossings[T]
    p, q = O_pq(m)
    epsnum, epsden = p - T * q, q                    # eps = O_m - T exactly
    L = mp.log10(mp.pi) + m * (mp.log10(mp.mpf(epsnum)) - mp.log10(mp.mpf(epsden))) + log10_ratio_dfact(m)
    L2 = mp.log10(mp.pi) + m * (mp.log10(mp.mpf(epsnum)) - mp.log10(mp.mpf(epsden))) + log10_ratio_dfact(m, exact=False)
    man, e = man_exp(L)
    results[T] = (man, e)
    say("  T = %d (denominator %d): first-failure value = %s x 10^%d   (loggamma route agrees to %s)"
        % (T, 2 * m + 1, mp.nstr(man, 14), e, mp.nstr(abs(L - L2), 3)))
# 2 pi: eps = O_m - 2 pi with exact O_m and rational pi bounds -> an enclosure of the value
m = 40249
vals = []
for pib in (PI_HI, PI_LO):                     # eps_lo uses pi_hi
    eps = Fr(P2PI, Q2PI) - 2 * pib
    L = mp.log10(mp.pi) + m * (mp.log10(mp.mpf(eps.numerator)) - mp.log10(mp.mpf(eps.denominator))) \
        + log10_ratio_dfact(m, exact=False)
    vals.append(man_exp(L))
# pi itself inside log10(pi): use the mpmath value (error 1e-80, negligible)
say("  T = 2pi (denominator 80499): first-failure value in [%s, %s] x 10^%d"
    % (mp.nstr(vals[0][0], 14), mp.nstr(vals[1][0], 14), vals[0][1]))
say("     two-sided Poisson discrepancy sum_Z - int_R = 2 x value = %s x 10^%d"
    % (mp.nstr(2 * vals[0][0], 12), vals[0][1]))
check(vals[0][1] == vals[1][1] == -226577 and mp.nstr(vals[0][0], 11) == mp.nstr(vals[1][0], 11) == "4.2101759371",
      "2 pi value = 4.2101759371e-226577 (enclosure)")
check(mp.nstr(results[2][0], 11) == "2.3100572725" and results[2][1] == -11, "T=2 value 2.3100572725e-11")
check(mp.nstr(results[3][0], 11) == "2.3324296128" and results[3][1] == -138, "T=3 value 2.3324296128e-138")
check(mp.nstr(results[4][0], 11) == "4.5195214277" and results[4][1] == -1320, "T=4 value 4.5195214277e-1320")

# Published anchors for 1/15: BB 2001 eq.(10) and BBB 2008 eq.(5)
p7, q7 = O_pq(7)
eps7 = Fr(p7, q7) - 2
P7 = eps7 ** 7
for j in range(1, 8):
    P7 *= (2 * j + 1)
P7 /= math.factorial(7) * 2 ** 7
check(Fr(1, 2) - P7 == Fr(467807924713440738696537864469, 935615849440640907310521750000),
      "I_7(1)/pi = 467807924713440738696537864469/935615849440640907310521750000 (BB 2001 eq. 10)")
check(P7 == Fr(6879714958723010531, 935615849440640907310521750000),
      "deficit/pi = 6879714958723010531/935615849440640907310521750000 (BBB 2008 eq. 5)")

# Independent evaluation by the BB 2001 Theorem 2 / Wikipedia general formula (all 2^7 sign vectors)
a0, a = Fr(1), [Fr(1, 2 * k + 1) for k in range(1, 8)]
n = len(a)
S = Fr(0)
for g in iproduct((1, -1), repeat=n):
    b = a0 + sum(gi * ai for gi, ai in zip(g, a))
    eg = 1
    for gi in g:
        eg *= gi
    sg = (b > 0) - (b < 0)
    S += eg * b ** n * sg
prod_a = Fr(1)
for ai in a:
    prod_a *= ai
Cn = S / (2 ** n * math.factorial(n) * prod_a)
check(Cn / 2 == Fr(467807924713440738696537864469, 935615849440640907310521750000),
      "general sign-sum formula (BB 2001 Thm 2) reproduces the published 1/15 value independently")

# 2cos variant at 113: how many sign vectors gamma give b_gamma < 0 (a0 = 2, a_k = 1/(2k+1), k=1..56)?
a56 = [Fr(1, 2 * k + 1) for k in range(1, 57)]
s56 = sum(a56)
excess = s56 - 2
neg = sum(1 for k in range(56) if 2 * a56[k] < excess)       # a single flip that keeps b < 0
check(excess > 0 and neg == 0, "1/113 case: exactly one corner (all gamma_k = -1) has b_gamma < 0; excess = %s"
      % mp.nstr(mp.mpf(excess.numerator) / excess.denominator, 10))

# Poisson normalisation check against BBB 2008 Example 4: s_7 + 1/2 - i_7 = 2 pi (7 - 2 pi)^6 / (6! 2^7)
import sympy as sp
x = sp.symbols('x')
bbb = (sp.Rational(43141, 15360) - sp.Rational(5887, 23040)) * x - sp.Rational(16807, 3840) * x ** 2 \
    + sp.Rational(2401, 768) * x ** 3 - sp.Rational(343, 288) * x ** 4 + sp.Rational(49, 192) * x ** 5 \
    - sp.Rational(7, 240) * x ** 6 + sp.Rational(1, 720) * x ** 7
mine = sp.expand(2 * x * (7 - 2 * x) ** 6 / (720 * 128))
check(sp.expand(bbb - mine) == 0,
      "one-sided sum-minus-integral = F(2 pi) (single cosine term): BBB Example 4 exact s_7, i_7 reproduced")

# ======================================================================== 3
say("")
say("=" * 100)
say("3. 'Naive double precision' mantissa at 80499 (sensitivity to summation order)")
say("=" * 100)
import numpy as np
m = 40249
terms = [1.0 / (2 * j + 1) for j in range(m + 1)]
fw = 0.0
for t in terms:
    fw += t
bw = 0.0
for t in reversed(terms):
    bw += t
fs = math.fsum(terms)
npsum = float(np.sum(np.array(terms)))
mp.mp.dps = 60
exactO = mp.mpf(P2PI) / Q2PI
lr = log10_ratio_dfact(m, exact=False)
rows = []
for name, sval in [("forward (python loop)", fw), ("reverse (smallest first)", bw),
                   ("math.fsum (correctly rounded)", fs), ("numpy.sum (pairwise)", npsum)]:
    for piname, twopi in [("exact 2pi", 2 * mp.pi), ("double 2*math.pi", mp.mpf(2 * math.pi))]:
        eps = mp.mpf(sval) - twopi
        L = mp.log10(mp.pi) + m * mp.log10(eps) + lr
        man, e = man_exp(L)
        rows.append((name, piname, sval, eps, man, e))
        say("  %-30s %-17s O = %.17g  eps = %s  mantissa = %s x 10^%d"
            % (name, piname, sval, mp.nstr(eps, 12), mp.nstr(man, 6), e))
say("  exact O_40249 - 2pi = %s" % mp.nstr(exactO - 2 * mp.pi, 15))
check(abs(rows[0][4] - mp.mpf("4.208")) < 0.001, "forward double sum + exact 2pi gives 4.208 (session's claim reproduced)")

# ======================================================================== 4
say("")
say("=" * 100)
say("4. Crossing-index window test with INTERVAL arithmetic: m(n) = floor(X), X = exp(2n - gamma)/4,")
say("   certified when frac(X) >= X/(24 floor(X)^2)  [O_m = ln(4M)/2 + gamma/2 + r(M), 0 < r < 1/(48M^2)]")
say("=" * 100)
tC = time.time()
exc, undecided, fracs, hasler = [], [], [], []
worst = None
for nn in range(2, NMAX + 1):
    iv.dps = int(0.8686 * nn) + 30
    mp.mp.dps = iv.dps
    X = iv.exp(iv.mpf(2 * nn) - iv.euler) / 4
    Fl = int(LO(X))                      # exact floor of the lower endpoint (python int)
    if int(HI(X)) != Fl:
        undecided.append(nn)
        continue
    fr = X - iv.mpf(Fl)                  # interval for frac(X)
    win = X / (24 * iv.mpf(Fl) ** 2)     # interval for the window
    if not (LO(fr) > HI(win)):
        exc.append(nn)
    ratio = LO(fr) / HI(win)
    if worst is None or ratio < worst[1]:
        worst = (nn, ratio)
    fracs.append(float((LO(fr) + HI(fr)) / 2))
    if LO(fr) >= mp.mpf(7) / 8:
        hasler.append(nn)
    elif HI(fr) >= mp.mpf(7) / 8:
        undecided.append(("hasler", nn))
iv.dps = 30
say("  n = 2..%d : window exceptions %s ; undecided %s ; tightest n = %d with frac/window >= %s  (%.1fs)"
    % (NMAX, exc, undecided, worst[0], mp.nstr(worst[1], 6), time.time() - tC))
check(exc == [] and undecided == [], "m(n) = floor(exp(2n - gamma)/4) certified for 2 <= n <= %d (interval arithmetic)" % NMAX)
check(abs(worst[1] - 101.84) < 0.01 and worst[0] == 2, "tightest case n = 2, frac/window = 101.8")
fs_ = sorted(fracs)
N = len(fs_)
D = max(max((i + 1) / N - v, v - i / N) for i, v in enumerate(fs_))
say("  KS distance of frac(X_n) from uniform: D = %.4f (N = %d; 1%% critical 1.63/sqrt(N) = %.4f)" % (D, N, 1.63 / math.sqrt(N)))
check(abs(D - 0.0192) < 0.00005 and D < 1.63 / math.sqrt(N), "KS D = 0.0192 < 0.0298")
mp.mp.dps = 30
heur = mp.fsum([1 / (6 * mp.exp(2 * nn - mp.euler)) for nn in range(2, 200)])
heur61 = mp.fsum([1 / (6 * mp.exp(2 * nn - mp.euler)) for nn in range(61, 400)])
say("  heuristic expected exceptions: all n >= 2: %s ; n >= 61: %s" % (mp.nstr(heur, 4), mp.nstr(heur61, 3)))
say("  Hasler floor(X + 1/8) differs from floor(X) at %d of the n <= %d; first %s" % (len(hasler), NMAX, hasler[:10]))
check(len(hasler) == 375 and hasler[:8] == [17, 31, 45, 57, 59, 63, 86, 118], "375 Hasler failures, first 17, 31, 45, 57, 59, 63, 86, 118")

# OEIS DATA (read 2026-10-28) vs floor(X)
OEIS = [1, 7, 56, 418, 3091, 22845, 168803, 1247297, 9216353, 68100150, 503195828, 3718142207,
        27473561357, 203003686105, 1500005624923, 11083625711270, 81897532160124, 605145459495140,
        4471453748222756, 33039822589391675, 244133102611731230, 1803913190804074903]
mp.mp.dps = 60
fl = [int(mp.floor(mp.exp(2 * k - mp.euler) / 4)) for k in range(1, 23)]
fl8 = [int(mp.floor(mp.exp(2 * k - mp.euler) / 4 + mp.mpf(1) / 8)) for k in range(1, 23)]
check(fl == OEIS, "OEIS A092315 DATA a(1..22) = floor(exp(2n - gamma)/4) for every listed term")
diff8 = [k + 1 for k in range(22) if fl8[k] != OEIS[k]]
say("  OEIS DATA terms where Hasler's +1/8 formula disagrees with the DATA itself: n = %s" % diff8)
check(diff8 == [17], "the entry's own DATA contradicts the +1/8 formula exactly at n = 17 (within n <= 22)")

# n = 17 certification by the DLMF 5.11(ii) envelope (rigorous) and by digamma
iv.dps = 50


def H_iv(Nint):
    Nv = iv.mpf(Nint)
    L = iv.log(Nv) + iv.euler + 1 / (2 * Nv) - 1 / (12 * Nv ** 2) + 1 / (120 * Nv ** 4)
    return iv.mpf([LO(L - 1 / (252 * Nv ** 6)), HI(L)])


def O_iv(mint):
    M = mint + 1
    return H_iv(2 * M) - H_iv(M) / 2


a17 = 81897532160124
up = O_iv(a17) - 17
dn = O_iv(a17 - 1) - 17
say("  n = 17, a = %d : O_a - 17 in [%s, %s] ; O_(a-1) - 17 in [%s, %s]"
    % (a17, mp.nstr(LO(up), 8), mp.nstr(HI(up), 8), mp.nstr(LO(dn), 8), mp.nstr(HI(dn), 8)))
mp.mp.dps = 60
Odg = lambda mm: (mp.digamma(mm + mp.mpf(3) / 2) - mp.digamma(mp.mpf(1) / 2)) / 2
say("     digamma: O_a - 17 = %s ; O_(a-1) - 17 = %s" % (mp.nstr(Odg(a17) - 17, 6), mp.nstr(Odg(a17 - 1) - 17, 6)))
X17 = mp.exp(34 - mp.euler) / 4
check(LO(up) > 0 and HI(dn) < 0 and int(mp.floor(X17 + mp.mpf(1) / 8)) == a17 + 1,
      "a(17) = 81897532160124 certified (envelope + intervals); floor(X+1/8) = %d" % (a17 + 1))

# k*pi thresholds, k = 1..KMAX, interval window test
exc_pi, und_pi, firsts = [], [], []
for k in range(1, KMAX + 1):
    iv.dps = int(2.7288 * k) + 30
    mp.mp.dps = iv.dps
    X = iv.exp(2 * k * iv.pi - iv.euler) / 4
    Fl = int(LO(X))
    if int(HI(X)) != Fl:
        und_pi.append(k)
        continue
    fr = X - iv.mpf(Fl)
    win = X / (24 * iv.mpf(Fl) ** 2)
    if not (LO(fr) > HI(win)):
        exc_pi.append(k)
    if k <= 3:
        firsts.append((k, Fl, 2 * Fl + 1, float((LO(fr) + HI(fr)) / 2)))
say("  k*pi, k = 1..%d: window exceptions %s, undecided %s ; (k, m, 2m+1, frac): %s" % (KMAX, exc_pi, und_pi, firsts))
check(exc_pi == [] and und_pi == [] and [f[1] for f in firsts] == [75, 40249, 21553437],
      "m(k pi) = floor(X) for k <= %d; m(pi), m(2pi), m(3pi) = 75, 40249, 21553437" % KMAX)

# m(3 pi) by a second, elementary route: math.fsum of correctly rounded doubles, explicit error bound
m3 = 21553437
arr = 1.0 / (2.0 * np.arange(m3 + 1, dtype=np.float64) + 1.0)
s_hi = math.fsum(arr)                    # O_m3 (rounded)
s_lo = s_hi - float(arr[-1])             # O_{m3-1} (approximately)
err = 2.0 ** -53 * 10.0 + 2.0 ** -52 * 10.0
mp.mp.dps = 30
d_up = mp.mpf(s_hi) - 3 * mp.pi
d_dn = mp.mpf(s_lo) - 3 * mp.pi
say("  fsum route: O_m - 3pi = %s ; O_(m-1) - 3pi = %s ; rigorous rounding error < %.1e" % (mp.nstr(d_up, 6), mp.nstr(d_dn, 6), err))
check(d_up > err and d_dn < -err, "m(3 pi) = 21553437 confirmed by fsum with explicit error bound")
iv.dps = 40
up = O_iv(m3) - 3 * iv.pi
dn = O_iv(m3 - 1) - 3 * iv.pi
say("  envelope route: O_m - 3pi in [%s, %s] ; O_(m-1) - 3pi in [%s, %s]" % (mp.nstr(LO(up), 8), mp.nstr(HI(up), 8), mp.nstr(LO(dn), 8), mp.nstr(HI(dn), 8)))
check(LO(up) > 0 and HI(dn) < 0, "m(3 pi) = 21553437 confirmed by the DLMF envelope in interval arithmetic")

# ======================================================================== 5
say("")
say("=" * 100)
say("5. The remainder r(M) = O_(M-1) - ln(4M)/2 - gamma/2: algebra and bounds")
say("=" * 100)
Msym = sp.symbols('M', positive=True)


def Hasy(Nv):
    return sp.log(Nv) + sp.EulerGamma + 1 / (2 * Nv) - sp.Rational(1, 12) / Nv ** 2 + sp.Rational(1, 120) / Nv ** 4


r_expr = sp.simplify(Hasy(2 * Msym) - Hasy(Msym) / 2 - sp.log(4 * Msym) / 2 - sp.EulerGamma / 2)
check(sp.simplify(r_expr - (sp.Rational(1, 48) / Msym ** 2 - sp.Rational(7, 1920) / Msym ** 4)) == 0,
      "r(M) = 1/(48M^2) - 7/(1920M^4) - eps_2M + eps_M/2 (algebra)")
# exact-rational/high-precision direct check for M = 1..3000
mp.mp.dps = 50
run = mp.mpf(0)
ok = True
for M in range(1, 3001):
    run += mp.mpf(1) / (2 * M - 1)                  # O_(M-1)
    r = run - (mp.log(4 * M) + mp.euler) / 2
    if not (0 < r < 1 / (48 * mp.mpf(M) ** 2)):
        ok = False
        say("   bound fails at M =", M)
    lo = 1 / (48 * mp.mpf(M) ** 2) - mp.mpf(7) / (1920 * mp.mpf(M) ** 4) - 1 / (252 * (2 * mp.mpf(M)) ** 6)
    hi = 1 / (48 * mp.mpf(M) ** 2) - mp.mpf(7) / (1920 * mp.mpf(M) ** 4) + 1 / (504 * mp.mpf(M) ** 6)
    if not (lo < r < hi):
        ok = False
        say("   envelope fails at M =", M)
check(ok, "0 < r(M) < 1/(48M^2) and the two-sided envelope hold for M = 1..3000 (50-digit direct sums)")
# GKP/DLMF envelope for H_N itself at small N (exact H_N)
okH = True
for Nn in range(1, 200):
    Hn = mp.fsum([mp.mpf(1) / k for k in range(1, Nn + 1)])
    L = mp.log(Nn) + mp.euler + mp.mpf(1) / (2 * Nn) - mp.mpf(1) / (12 * Nn ** 2) + mp.mpf(1) / (120 * Nn ** 4)
    e = L - Hn
    if not (0 < e < 1 / (252 * mp.mpf(Nn) ** 6)):
        okH = False
check(okH, "0 < eps_N < 1/(252 N^6) (GKP (9.89); DLMF 5.11.2 + 5.11(ii)) checked directly for N = 1..199")

# ======================================================================== 6
say("")
say("=" * 100)
say("6. Normalisation of (1)-(2) by quadrature (m = 4), independent of all closed forms")
say("=" * 100)


def tail_exact(avals, t):
    """P(sum U_j > t), U_j uniform on [-a_j, a_j]: inclusion-exclusion over all subsets (own code)."""
    n = len(avals)
    A = sum(avals)
    c = A - t
    tot = Fr(0)
    for mask in range(1 << n):
        used = sum(2 * avals[j] for j in range(n) if mask >> j & 1)
        r = c - used
        if r > 0:
            tot += (-1) ** bin(mask).count("1") * r ** n
    den = math.factorial(n)
    for v in avals:
        den *= 2 * v
    return tot / den


mp.mp.dps = 25
a4 = [Fr(1, 3), Fr(1, 5), Fr(1, 7), Fr(1, 9)]
f4 = lambda xx, xi: mp.sinc(xx) * mp.sinc(xx / 3) * mp.sinc(xx / 5) * mp.sinc(xx / 7) * mp.sinc(xx / 9) * mp.cos(xi * xx)
for xi in [Fr(3, 2), Fr(7, 4)]:
    pts = [mp.pi * k for k in range(0, 64)]
    q = 2 * mp.quad(lambda xx: f4(xx, mp.mpf(xi.numerator) / xi.denominator), pts)
    tail = 2 * 945 / (4 * pts[-1] ** 4)
    pr = tail_exact(a4, xi - 1) - tail_exact(a4, xi + 1)          # P(xi-1 < S' <= xi+1), exact rational
    predv = mp.pi * mp.mpf(pr.numerator) / pr.denominator
    say("  F_4(%s): quadrature %s ; pi*P(xi-1 <= S' <= xi+1) = %s ; |tail| <= %s"
        % (xi, mp.nstr(q, 12), mp.nstr(predv, 12), mp.nstr(tail, 2)))
    check(abs(q - predv) < tail + 1e-12, "F_m(xi) = pi*P(xi-1 <= S' <= xi+1) at xi = %s (quadrature)" % xi)

say("")
say("FAILURES: %s" % (FAIL if FAIL else "none"))
say("total %.1fs" % (time.time() - T0))

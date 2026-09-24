#!/usr/bin/env python3
"""collatz_procgen_20260924_fixedpt_chains.py

Lane: fixed points (session collatz-procgen-20260922, 2026-09-24).
The correct fixed-point theorems for Collatz, and the "fixed point chain" and its growth.

Shortcut map T(x) = x/2 (x even), (3x+1)/2 (x odd).  Inverse branches D(x) = 2x, E(x) = (2x-1)/3.
A forward parity word w = (w_0..w_(p-1)) with a ones has T^p(x) = (3^a x + c_w)/2^p on its cylinder,
c_() = 0, c_(z e) = 3^e c_z + e 2^|z|, and its inverse word W = G_(w_0) o ... o G_(w_(p-1)) (G_0 = D,
G_1 = E) has the unique fixed point x_w = c_w/(2^p - 3^a).

Sections (printed):
  F1  Banach in Z_2: D, E contract by 1/2; every word's inverse map has exactly one fixed point, the
      rational periodic point x_w; Banach iteration reaches it; T^p(x_w) = x_w with parity word w
  F2  growth of the fixed-point chain by period p <= 24: 2^p points, primitive points, cycles, the sign
      split, the points in Bad (reproducing the mod-192 note's K5 counts), and the integral points
  F3  Banach in R and Z_3: the inverse word is a real contraction iff x_w < 0 (sign law); Z_3 expands
  F4  Brouwer/IVT for continuous extensions: Chamberland's f = T on Z; integer-cycle multipliers
      3^a/2^p (attracting iff >= 0); Chamberland's cited numbers re-derived; covering chains and
      periodic points of every period on BOTH sides; a monotone divergent orbit
  F5  parity principles: Sperner (1D, 2D) and Redei (n <= 6), and the Banach count "exactly one"

Every check raises on failure.  Runtime about 2-3 min; peak memory about 400 MB (the p = 24 level of F2);
one process.
"""
import sys
import time
import random
import resource
import itertools
from fractions import Fraction as Fr
from math import gcd

import numpy as np
import mpmath as mp

T0 = time.time()
random.seed(924)


def check(cond, msg):
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)


def hdr(s):
    print()
    print("=" * 100)
    print(s)
    print("=" * 100)


def Tb(x, b=1):
    return x // 2 if x % 2 == 0 else (3 * x + b) // 2


def T_rat(x):
    """T on rationals with odd denominator (2-adic integers in Q)"""
    num, den = x.numerator, x.denominator
    check(den % 2 == 1, "odd denominator")
    return x / 2 if num % 2 == 0 else (3 * x + 1) / 2


def c_word(w):
    R = 0
    for i, e in enumerate(w):
        R = 3 ** e * R + e * 2 ** i
    return R


def v2(n):
    n = abs(n)
    if n == 0:
        return 10 ** 9
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def v3(n):
    n = abs(n)
    if n == 0:
        return 10 ** 9
    k = 0
    while n % 3 == 0:
        n //= 3
        k += 1
    return k


def mobius(n):
    res, m, p = 1, n, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            res = -res
        p += 1
    if m > 1:
        res = -res
    return res


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


# ----------------------------------------------------------------------------------------------
hdr("F1  Banach in Z_2: every inverse word has exactly one fixed point, the rational periodic point")
# ----------------------------------------------------------------------------------------------
MOD = 2 ** 64
INV3 = pow(3, -1, MOD)


def W_apply_mod(w, y):
    for e in reversed(w):
        y = (2 * y) % MOD if e == 0 else ((2 * y - 1) * INV3) % MOD
    return y


# contraction factor 1/2 of D and E in Z_2
for _ in range(2000):
    x, y = random.randrange(-10 ** 12, 10 ** 12), random.randrange(-10 ** 12, 10 ** 12)
    if x == y:
        continue
    check(v2(2 * x - 2 * y) == v2(x - y) + 1, "D contracts by 1/2")
    # E(x) - E(y) = 2(x - y)/3 : v_2 increases by exactly 1
    check(v2(2 * (x - y)) == v2(x - y) + 1, "E contracts by 1/2 (3 is a 2-adic unit)")
check(W_apply_mod((0,), 0) == 0 and W_apply_mod((1,), MOD - 1) == MOD - 1, "fixed points 0 (D), -1 (E)")
nwords = 0
pts_by_p = {}
for p in range(1, 13):
    seen = set()
    for bits in range(2 ** p):
        w = tuple((bits >> t) & 1 for t in range(p))
        a = sum(w)
        cw = c_word(w)
        xw = Fr(cw, 2 ** p - 3 ** a)
        # exact periodicity and parity word
        y = xw
        for t in range(p):
            check((y.numerator % 2) == w[t], "parity word of x_w is w")
            y = T_rat(y)
        check(y == xw, "T^p(x_w) = x_w")
        # Banach iteration from 0 converges 2-adically (mod 2^64 after ceil(64/p) rounds)
        z = 0
        for _ in range(64 // p + 2):
            z = W_apply_mod(w, z)
        target = (xw.numerator * pow(xw.denominator, -1, MOD)) % MOD
        check(z == target, "Banach iterate = x_w mod 2^64")
        seen.add(xw)
        nwords += 1
    check(len(seen) == 2 ** p, "distinct words give distinct fixed points")
    pts_by_p[p] = seen
print(f"  D(x) = 2x and E(x) = (2x-1)/3 are 2-adic contractions with factor 1/2 (3 is a unit in Z_2); their")
print("  fixed points are 0 and -1.  For every word w of length p <= 12 (all 8190 words):")
print("    the inverse word W has exactly one fixed point in Z_2 (Banach), namely x_w = c_w/(2^p - 3^a);")
print("    Banach iteration from 0 agrees with x_w modulo 2^64 after ceil(64/p)+1 rounds;")
print("    T^p(x_w) = x_w exactly, with parity word w; distinct words give distinct points.")
print("  So Per_p(T) in Z_2 = {Banach fixed points of the 2^p words of length p}: the 'fixed point chain'.")

# ----------------------------------------------------------------------------------------------
hdr("F2  Growth of the fixed-point chain by period p <= 24; integral points; the points in Bad")
# ----------------------------------------------------------------------------------------------
PMAX = 24
LANE_I = [1, 0, 1, 2, 3, 6, 12, 16, 36, 60, 127, 216, 366, 721, 1290, 2095, 4227, 7451]
amin = {j: next(a for a in range(j + 2) if 3 ** a > 2 ** j) for j in range(1, PMAX + 1)}
c = np.zeros(1, dtype=np.int64)
a_arr = np.zeros(1, dtype=np.int8)
good = np.ones(1, dtype=bool)
stats = {}
int_values = set()
for p in range(1, PMAX + 1):
    c = np.concatenate([c, 3 * c + (1 << (p - 1))])
    a_new = np.concatenate([a_arr, a_arr + 1])
    del a_arr
    a_arr = a_new
    good = np.concatenate([good, good]) & (a_arr >= amin[p])
    n_all = c.shape[0]
    n_bad = int(np.count_nonzero(good))
    n_pos = n_neg = n_zero = n_int = 0
    for av in range(0, p + 1):
        den = 2 ** p - 3 ** av
        mask = a_arr == av
        cnt = int(np.count_nonzero(mask))
        if av == 0:
            n_zero += cnt
        elif den > 0:
            n_pos += cnt
        else:
            n_neg += cnt
        cs = c[mask]
        hit = cs[cs % abs(den) == 0]
        n_int += int(hit.shape[0])
        for v in np.unique(hit):
            int_values.add(int(v) // den)
        del cs, mask, hit
    stats[p] = dict(all=n_all, bad=n_bad, pos=n_pos, neg=n_neg, zero=n_zero, integral=n_int)


def prim(key, p):
    return sum(mobius(p // d) * stats[d][key] for d in divisors(p))


print("   p | points 2^p | primitive | cycles | primitive x_w>0 | primitive x_w<0 | in Bad (primitive) | integral")
rows = []
for p in range(1, PMAX + 1):
    pa = prim("all", p)
    pp_, pn_, pb_, pi_ = prim("pos", p), prim("neg", p), prim("bad", p), prim("integral", p)
    pz_ = prim("zero", p)
    check(pa == pp_ + pn_ + pz_ and pa % p == 0, "primitive counts consistent")
    rows.append((p, pa, pb_, pi_))
    print(f"  {p:2d} | {2 ** p:10d} | {pa:9d} | {pa // p:6d} | {pp_:15d} | {pn_:15d} | {pb_:18d} | {pi_:8d}")
check([prim("bad", p) for p in range(1, 19)] == LANE_I, "Bad counts reproduce the mod-192 note K5 (p <= 18)")
expected_int = {0, 1, 2, -1, -5, -7, -10, -17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34}
check(int_values == expected_int, "integral periodic points for p <= 24 = the five known cycles")
check([prim("integral", p) for p in range(1, PMAX + 1)] ==
      [2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 11] + [0] * (PMAX - 11), "integral primitive points by period")
print(f"  * the in-Bad column reproduces the mod-192 note's K5 counts for p <= 18 (independent code) and")
print(f"    extends them to p = {PMAX}: {', '.join(str(prim('bad', p)) for p in range(19, PMAX + 1))}.")
print(f"  * integral points, all p <= {PMAX}: {sorted(int_values)} = the five known cycles")
print("    {0}, {1,2}, {-1}, {-5,-7,-10}, {-17,...,-136} (periods 1, 2, 1, 3, 11).  Positive integral points occur")
print("    only on contracting words (2^p > 3^a), negative ones only on expanding words (sign law).")
bad_growth = [np.log2(prim("bad", p)) / p for p in (16, 20, 24)]
print(f"  * growth: all points 2^p (exponent 1); points in Bad 2^(e p) with log2(count)/p = "
      f"{', '.join(f'{g:.3f}' for g in bad_growth)} at p = 16, 20, 24")
print("    (slowly approaching h(log_3 2) = 0.950, the PROVED exponent of the word count |Bad_p| in the choice")
print("    ladder; the polynomial factor is still visible at p = 24).  Integral points: 18 in total, no growth.")
del c, a_arr, good

# ----------------------------------------------------------------------------------------------
hdr("F3  Banach in R and in Z_3: the sign law decides which words are real contractions")
# ----------------------------------------------------------------------------------------------
cnt_contr = cnt_exp = 0
for p in range(1, 11):
    for bits in range(2 ** p):
        w = tuple((bits >> t) & 1 for t in range(p))
        a = sum(w)
        if a == 0:
            continue
        xw = Fr(c_word(w), 2 ** p - 3 ** a)
        slope = Fr(2 ** p, 3 ** a)          # W(y) = (2^p y - c_w)/3^a
        check((slope < 1) == (xw < 0), "real contraction iff x_w < 0")

        def W_real(y, w=w):
            for e in reversed(w):
                y = 2 * y if e == 0 else (2 * y - 1) / 3
            return y
        y0 = Fr(7, 3)
        y1 = W_real(y0)
        check(y1 - xw == slope * (y0 - xw), "W is affine with slope 2^p/3^a about x_w")
        if slope < 1:
            cnt_contr += 1
        else:
            cnt_exp += 1
print(f"  words of length <= 10 with a >= 1: {cnt_contr} inverse maps are real contractions (slope 2^p/3^a < 1),")
print(f"  exactly those with x_w < 0; {cnt_exp} are real expansions, exactly those with x_w > 0.")
for w, target in (((1,), -1), ((1, 1, 0), -5), ((1, 0), 1)):
    y = 0.0
    for _ in range(400):
        for e in reversed(w):
            y = 2 * y if e == 0 else (2 * y - 1) / 3
    a = sum(w)
    if 3 ** a > 2 ** len(w):
        check(abs(y - target) < 1e-9, "backward iteration converges to the negative cycle point")
    else:
        check(abs(y) > 1e10 or abs(y - target) > 1, "backward iteration escapes from a positive cycle point")
w17 = None
for bits in range(2 ** 11):
    w = tuple((bits >> t) & 1 for t in range(11))
    if sum(w) == 7 and Fr(c_word(w), 2 ** 11 - 3 ** 7) == -17:
        w17 = w
check(w17 is not None, "word of -17")
y = 0.0
for _ in range(3000):
    for e in reversed(w17):
        y = 2 * y if e == 0 else (2 * y - 1) / 3
check(abs(y + 17) < 1e-6, "backward iteration converges to -17 (slope 2048/2187)")
print("  backward (inverse-word) iteration in R converges to -1 (slope 2/3), -5 (8/9) and -17 (2048/2187),")
print("  and is repelled by 1 (word 10, slope 4/3).  So Banach in R sees exactly the NEGATIVE cycles; the")
print("  positive ones are forward-attracting.  In Z_2 every word contracts (factor 2^-p): Banach sees all.")
for _ in range(3000):
    x, y = random.randrange(-10 ** 9, 10 ** 9), random.randrange(-10 ** 9, 10 ** 9)
    if x == y or (2 * x - 1) % 3 or (2 * y - 1) % 3:
        continue
    check(v3((2 * x - 1) // 3 - (2 * y - 1) // 3) == v3(x - y) - 1, "E expands by 3 in Z_3")
    check(v3(2 * x - 2 * y) == v3(x - y), "D is a 3-adic isometry")
print("  In Z_3, D is an isometry and E expands by 3: no contraction principle at the prime 3, where")
print("  integrality (the gate (2^p - 3^a) | c_w and the legality x = 2 mod 3) is decided.")

# ----------------------------------------------------------------------------------------------
hdr("F4  Brouwer (IVT) for continuous extensions: Chamberland's f, Sharkovskii, covering chains")
# ----------------------------------------------------------------------------------------------
for n in range(-10000, 10001):
    fn = n + Fr(1, 4) - Fr(2 * n + 1, 4) * (1 if n % 2 == 0 else -1)
    check(fn == Tb(n), "f(n) = T(n) exactly (cos(pi n) = (-1)^n)")
print("  f(x) = (x/2) cos^2(pi x/2) + ((3x+1)/2) sin^2(pi x/2) = x + 1/4 - ((2x+1)/4) cos(pi x)  [Chamberland 1996]:")
print("  f(n) = T(n) exactly for |n| <= 10^4, and f'(n) = 1 - (-1)^n/2 in {1/2, 3/2}.")
cycles = {"{0}": [0], "{1,2}": [1, 2], "{-1}": [-1], "{-5,-7,-10}": [-5, -7, -10],
          "{-17,...}": [-17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34]}
for name, cyc in cycles.items():
    for i, x in enumerate(cyc):
        check(Tb(x) == cyc[(i + 1) % len(cyc)], "cycle")
    mult = Fr(1)
    for x in cyc:
        mult *= Fr(1, 2) if x % 2 == 0 else Fr(3, 2)
    a = sum(1 for x in cyc if x % 2)
    check(mult == Fr(3 ** a, 2 ** len(cyc)), "multiplier 3^a/2^p")
    check((mult < 1) == (cyc[0] >= 0), "attracting iff nonnegative (sign law)")
    print(f"    {name:14s} period {len(cyc):2d}  multiplier {str(mult):10s} {'attracting' if mult < 1 else 'repelling'}")
print("  An integer cycle of T is an f-cycle with multiplier 3^a/2^p; by the sign law it is attracting iff its")
print("  points are >= 0.  So in Chamberland's extension the SIGN decides the STABILITY of integer cycles.")

mp.mp.dps = 60
PI = mp.pi


def f(x):
    return x + mp.mpf(1) / 4 - (2 * x + 1) / 4 * mp.cos(PI * x)


def fnp(x):
    return x + 0.25 - (2 * x + 1) / 4 * np.cos(np.pi * x)


def fixed_points(lo, hi, N=20000):
    g = lambda x: f(x) - x
    xs = np.linspace(lo, hi, N + 1)
    gv = fnp(xs) - xs
    out = []
    for i in range(N):
        if gv[i] == 0:
            out.append(mp.mpf(xs[i]))
        elif gv[i] * gv[i + 1] < 0:
            out.append(mp.findroot(g, (mp.mpf(xs[i]), mp.mpf(xs[i + 1])), solver="bisect"))
    return out


fps = fixed_points(0.0001, 3.0)
mu1, mu2, mu3 = fps[0], fps[1], fps[2]
check(abs(mu1 - mp.mpf("0.27773")) < 1e-5 and abs(mu3 - mp.mpf("2.44570")) < 1e-5, "mu1, mu3 as cited")
r2 = mp.findroot(lambda x: f(f(x)) - x, mp.mpf("1.1925"))
check(abs(r2 - mp.mpf("1.192531907")) < 1e-9 and abs(f(r2) - mp.mpf("2.138656335")) < 1e-9, "A2 as cited")
m2 = mp.diff(f, r2) * mp.diff(f, f(r2))
check(abs(m2) < 1, "A2 attracting")
xs = np.linspace(float(mu1), float(mu3), 400001)
fv = fnp(xs)
check(fv.min() >= float(mu1) - 1e-9 and fv.max() <= float(mu3) + 1e-9, "[mu1, mu3] invariant (grid)")
print(f"  Chamberland's numbers re-derived: mu1 = {mp.nstr(mu1, 12)}, mu2 = {mp.nstr(mu2, 12)}, mu3 = {mp.nstr(mu3, 12)};")
print(f"  second attracting 2-cycle A2 = {{{mp.nstr(r2, 12)}, {mp.nstr(f(r2), 12)}}} (multiplier {mp.nstr(m2, 6)});")
print(f"  f([mu1, mu3]) inside [mu1, mu3] on a 4*10^5-point grid (range [{fv.min():.4f}, {fv.max():.4f}]).")
print("  These match the values cited from Chamberland 1996 via Lagarias's bibliography (entry 35) and")
print("  Chamberland's 2003 survey (section 6.5).")


def roots_level(J, c_, N=6000):
    lo, hi = J
    xs_ = np.linspace(float(lo), float(hi), N + 1)
    gv = fnp(xs_) - float(c_)
    out = []
    for i in range(N):
        if gv[i] * gv[i + 1] < 0:
            a_, b_ = mp.mpf(xs_[i]), mp.mpf(xs_[i + 1])
            fa = f(a_) - c_
            for _ in range(230):
                m_ = (a_ + b_) / 2
                fm = f(m_) - c_
                if fm == 0:
                    a_ = b_ = m_
                    break
                if (fm > 0) == (fa > 0):
                    a_, fa = m_, fm
                else:
                    b_ = m_
            out.append((a_ + b_) / 2)
        elif gv[i] == 0:
            out.append(mp.mpf(xs_[i]))
    return out


def pullback(J, K):
    """compact L inside J with f(L) = K = [c,d] (given f(J) contains K): adjacent c- and d-preimages"""
    cc, dd = K
    rc = [(r_, 0) for r_ in roots_level(J, cc)]
    rd = [(r_, 1) for r_ in roots_level(J, dd)]
    allr = sorted(rc + rd)
    for (u, tu), (v, tv) in zip(allr[:-1], allr[1:]):
        if tu != tv:
            return (u, v)
    raise AssertionError("CHECK FAILED: pullback interval not found")


def periodic_point(itin, Js):
    """IVT chain: K_n = J_(u_0); K_k inside J_(u_k) with f(K_k) = K_(k+1); then a fixed point of f^n in K_0"""
    n = len(itin)
    K = (mp.mpf(Js[itin[0]][0]), mp.mpf(Js[itin[0]][1]))
    for k in range(n - 1, -1, -1):
        K = pullback(Js[itin[k]], K)
    lo, hi = K

    def g(x):
        y = x
        for _ in range(n):
            y = f(y)
        return y - x
    glo, ghi = g(lo), g(hi)
    check(glo * ghi <= 0, "f^n - id changes sign on K_0 (IVT)")
    for _ in range(260):
        m_ = (lo + hi) / 2
        gm = g(m_)
        if (gm > 0) == (glo > 0):
            lo, glo = m_, gm
        else:
            hi = m_
    x = (lo + hi) / 2
    orb = [x]
    for _ in range(n - 1):
        orb.append(f(orb[-1]))
    check(abs(f(orb[-1]) - x) < mp.mpf(10) ** -30, "periodic")
    for k, y in enumerate(orb):
        lo_, hi_ = Js[itin[k]]
        check(lo_ - 1e-30 <= y <= hi_ + 1e-30, "itinerary respected")
    return orb


def lyndon(n, alphabet, allowed):
    out = []
    for wd in itertools.product(alphabet, repeat=n):
        rots = [wd[i:] + wd[:i] for i in range(n)]
        if wd != min(rots) or len(set(rots)) != n:
            continue
        if all((wd[i], wd[(i + 1) % n]) in allowed for i in range(n)):
            out.append(wd)
    return out


# the covering relations are forced by integer values only (valid for EVERY continuous extension of T)
check(Tb(2) == 1 and Tb(3) == 5 and Tb(4) == 2, "positive turbulence data")
check(Tb(-10) == -5 and Tb(-7) == -10 and Tb(-5) == -7, "negative 3-cycle data")
check(Tb(-4) == -2 and Tb(-3) == -4 and Tb(-2) == -1, "negative turbulence data")
sides = [
    ("positive side, J0=[2,3], J1=[3,4] (full 2-shift; f(J0), f(J1) contain J0 u J1)",
     {0: (2, 3), 1: (3, 4)}, {(0, 0), (0, 1), (1, 0), (1, 1)}),
    ("negative side, I0=[-10,-7], I1=[-7,-5] from the 3-cycle (f(I0) contains I0 u I1, f(I1) contains I0)",
     {0: (-10, -7), 1: (-7, -5)}, {(0, 0), (0, 1), (1, 0)}),
    ("negative side, [-4,-3], [-3,-2] (full 2-shift, like the positive side)",
     {0: (-4, -3), 1: (-3, -2)}, {(0, 0), (0, 1), (1, 0), (1, 1)}),
]
NMAX = 7
for label, Js, allowed in sides:
    counts = []
    allpts = []
    for n in range(1, NMAX + 1):
        words = lyndon(n, (0, 1), allowed)
        for wd in words:
            orb = periodic_point(wd, Js)
            # exact period n: points pairwise distinct
            for i in range(n):
                for j in range(i + 1, n):
                    check(abs(orb[i] - orb[j]) > mp.mpf(10) ** -12, "exact period")
            allpts.append((n, orb[0]))
        counts.append(len(words))
    print(f"  {label}:")
    print(f"    periodic orbits constructed by the IVT chain, one per admissible Lyndon itinerary, periods 1..{NMAX}: "
          f"{counts}")
print("  So every continuous extension of T has periodic points of EVERY period on the positive side (from the")
print("  three values T(2)=1, T(3)=5, T(4)=2: Block-Coppel turbulence) and on the negative side (from the 3-cycle")
print("  {-5,-7,-10} via Sharkovskii/Li-Yorke, and also from T(-4)=-2, T(-3)=-4, T(-2)=-1).  The sign does NOT")
print("  decide the Sharkovskii type of f; it decides only which INTEGER cycles exist (positive: periods 1, 2;")
print("  negative: 1, 3, 11, for p <= 24 by F2) and, in Chamberland's f, their stability.")
print(f"  Growth of the real chain: at least 2^n points of period dividing n in [2,4] (one per itinerary), the")
print("  same 2^n as the Banach count in Z_2 (F1).")
# monotone divergent orbit: [2j+2, 2j+3] -> [2j+4, 2j+5]
Jm = {j: (2 * j + 2, 2 * j + 3) for j in range(14)}
for j in range(13):
    check(Tb(2 * j + 2) == j + 1 and Tb(2 * j + 3) == 3 * j + 5 and j + 1 <= 2 * j + 4 and 3 * j + 5 >= 2 * j + 5,
          "f([2j+2, 2j+3]) contains [2j+4, 2j+5]")
K = (mp.mpf(Jm[13][0]), mp.mpf(Jm[13][1]))
for j in range(12, -1, -1):
    K = pullback(Jm[j], K)
x0 = (K[0] + K[1]) / 2
orb = [x0]
for _ in range(13):
    orb.append(f(orb[-1]))
for j, y in enumerate(orb):
    check(Jm[j][0] <= y <= Jm[j][1], "monotone divergent itinerary")
check(all(orb[i] < orb[i + 1] for i in range(13)), "strictly increasing")
print(f"  monotone divergent orbit (Chamberland; LSW: a Cantor set of them): f([2j+2,2j+3]) contains [2j+4,2j+5]")
print(f"  for all j >= 0, so nested compact pullbacks give x0 = {mp.nstr(x0, 15)} with f^j(x0) in [2j+2, 2j+3],")
print("  j = 0..13 (checked); by compactness some point does this for all j.  Divergence of a real extension")
print("  is not Collatz divergence: it lives between the integers.")

# ----------------------------------------------------------------------------------------------
hdr("F5  Parity principles: Sperner, Redei, and the Banach count")
# ----------------------------------------------------------------------------------------------
for _ in range(3000):
    N = random.randint(1, 40)
    lab = [0] + [random.randint(0, 1) for _ in range(N - 1)] + [1]
    full = sum(1 for i in range(N) if {lab[i], lab[i + 1]} == {0, 1})
    check(full % 2 == 1, "1D Sperner: odd number of fully labelled edges")
nsp = 0
for _ in range(400):
    n = random.randint(2, 14)
    lab = {}
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            allowed_l = [t for t, coord in enumerate((i, j, k)) if coord > 0]
            lab[(i, j)] = random.choice(allowed_l)
    full = 0
    for i in range(n):
        for j in range(n - i):
            tri = [(i, j), (i + 1, j), (i, j + 1)]
            if {lab[v] for v in tri} == {0, 1, 2}:
                full += 1
            if i + j < n - 1:
                tri2 = [(i + 1, j), (i, j + 1), (i + 1, j + 1)]
                if {lab[v] for v in tri2} == {0, 1, 2}:
                    full += 1
    check(full % 2 == 1, "2D Sperner: odd number of fully labelled triangles")
    nsp += 1


def ham_paths(n, adj):
    out = [sum(1 << u for u in range(n) if adj[v][u]) for v in range(n)]
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1 << n):
        row = dp[S]
        for v in range(n):
            cv = row[v]
            if not cv:
                continue
            nb = out[v] & ~S
            while nb:
                low = nb & -nb
                u = low.bit_length() - 1
                dp[S | low][u] += cv
                nb ^= low
    return sum(dp[(1 << n) - 1])


ntour = 0
hist = {}
for n in range(1, 7):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for mask in range(1 << len(pairs)):
        adj = [[False] * n for _ in range(n)]
        for t, (i, j) in enumerate(pairs):
            if (mask >> t) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
        H = ham_paths(n, adj)
        check(H % 2 == 1, "Redei: odd number of Hamiltonian paths")
        ntour += 1
        if n == 6:
            hist[H] = hist.get(H, 0) + 1
print(f"  Sperner: 3000 random 1D and {nsp} random 2D Sperner labellings -- the number of fully labelled cells is")
print("  always odd (the combinatorial core of Brouwer).")
print(f"  Redei: all {ntour} labelled tournaments on n <= 6 vertices have an odd number of Hamiltonian paths;")
print(f"  n = 6 values: {sorted(hist)}.")
print("  Common mechanism: |S| = |Fix(i)| mod 2 for an involution i (door-to-door pairing for Sperner; arc-flip")
print("  or OCF parity for Redei, repo THM-001/THM-002).  For Collatz the corresponding count is not a parity:")
print("  each word has EXACTLY ONE 2-adic fixed point (F1), and the whole difficulty is integrality at the gate")
print("  (2^p - 3^a) | c_w, a congruence modulo an odd number prime to 6 -- invisible to every mod-2 count.")

rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
rss_mb = rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024
print()
print(f"[fixedpt_chains] ALL CHECKS PASSED  ({time.time() - T0:.1f} s, peak RSS {rss_mb:.0f} MB)", file=sys.stderr)
print("[fixedpt_chains] ALL CHECKS PASSED")

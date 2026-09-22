#!/usr/bin/env python3
"""Adversarial audit of lane grand_circuit_typing (session collatz-mod6-20260917, wave 2026-09-21, mac-mini).

Independent recomputation of the key numbers of
  04-computation/experiments/collatz_mod6_20260921_grand_circuit_typing.py
using different code paths: pure-Python modular power sums (no numpy), sympy.factorint
instead of the spf sieve, a Python continued-fraction Pell solver instead of PARI quadunit,
mpmath instead of PARI for the continued fraction of log_2 3, a linear Mobius sieve, and
direct trajectory iteration for the parities.  Every audit item prints A<k> with a verdict.

Audit items
  A1  primes p = 1 mod 4 below 3000 (211) and the AAC data (u mod p, v2, v3, digits, mod-8 census)
      from the negative Pell equation x^2 - p y^2 = -1 plus the cube-root test for the half-integral unit.
  A2  Giuga direct sum, pure Python, n <= 5000; BBBG / user criteria n <= 10^5 by sympy.factorint;
      Carmichael and Giuga lists.
  A3  the trunk t_j = (4^j-1)/3: the explorer wrote R^j(1); the truth is R^j(0) = R^(j-1)(1) (off-by-one).
  A4  continued fraction of log_2 3 by mpmath; convergents and q||q alpha||; records; sum 388.
  A5  the S12 float64 minimum re-checked at high precision (n = 10946 = F_21) and the float64 artefact
      0.017734 vs exact 0.017732 at q = 190537.
  A6  S15 index slip: with parity(T^k x) = [mu(k+1) = 1] the correlation is sum mu(n)[mu(n+1)=1], NOT
      #{mu(n)=1}; corrected point with parity(T^k x) = [mu(k)=1] for k >= 1; both truncations to 32 and
      2000 bits, and the two partial correlations at N = 1999.
  A7  S16 correlations c_j, j <= 6, N = 10^6, by an independent sieve and direct iteration; periodicity
      mod 2^(j+1) and NON-periodicity mod 2^j (sharpness).
  A8  S20 / S21 coprimality counts.
  A9  repository greps: 182 tournament titles, 10 with '17', body hits (THM-868, THM-871 -- NOT 0),
      S^6 / S_2 x S_3 hits with the own-lane files excluded and each hit classified; Wythoff / Zeckendorf /
      golden-ratio hits in the audited notes (the note's 'none' is false).
  A10 S25 table 2^K - 3^L along the convergents.

Explicit raise only; RAM << 1 GB; runtime ~1 min.
Reproduce:
  python3 04-computation/experiments/collatz_mod6_20260921_grand_circuit_typing_audit.py \
      > 05-knowledge/results/collatz_mod6_20260921_grand_circuit_typing_audit.out
"""
import math
import os
import re
import sys
from fractions import Fraction

import sympy
from mpmath import mp, mpf, log, floor, nint, fabs, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))


def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def v_p(n, p):
    if n == 0:
        return 10 ** 9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


print("=" * 78)
print("AUDIT of lane grand_circuit_typing (independent recomputation)")
print("=" * 78)

# ---------------------------------------------------------------------------
print("\n### A1  AAC: negative Pell x^2 - p y^2 = -1 by continued fractions, then the half-integral unit")


def neg_pell(p):
    """Fundamental solution (x, y) of x^2 - p y^2 = -1 via the continued fraction of sqrt p."""
    a0 = math.isqrt(p)
    m, d, a = 0, 1, a0
    h1, h0, k1, k0 = a0, 1, 1, 0
    if h1 * h1 - p * k1 * k1 == -1:
        return h1, k1
    for _ in range(100000):
        m = d * a - m
        d = (p - m * m) // d
        a = (a0 + m) // d
        h1, h0 = a * h1 + h0, h1
        k1, k0 = a * k1 + k0, k1
        if h1 * h1 - p * k1 * k1 == -1:
            return h1, k1
    raise RuntimeError("no negative Pell solution found for p=%d" % p)


def fundamental_unit_max_order(p):
    """(t, u) with t^2 - p u^2 = -4 giving the fundamental unit (t + u sqrt p)/2 of the maximal order."""
    x, y = neg_pell(p)
    # if eps = (t+u sqrt p)/2 with t,u odd and eps^3 = x + y sqrt p, then x = t(t^2+3)/2
    r, exact = sympy.integer_nthroot(2 * x, 3)
    for t in range(max(1, r - 2), r + 3):
        if t % 2 == 1 and t * (t * t + 3) == 2 * x:
            uu = (t * t + 4)
            if uu % p == 0:
                u, ex = sympy.integer_nthroot(uu // p, 2)
                if ex and u % 2 == 1:
                    # confirm the cube
                    # (t + u s)^3 / 8 = (t^3 + 3 t u^2 p)/8 + (3 t^2 u + u^3 p)/8 s
                    xx = (t ** 3 + 3 * t * u * u * p)
                    yy = (3 * t * t * u + u ** 3 * p)
                    check(xx % 8 == 0 and yy % 8 == 0 and xx // 8 == x and yy // 8 == y, "cube check p=%d" % p)
                    return t, u
    return 2 * x, 2 * y


primes14 = [p for p in sympy.primerange(5, 3000) if p % 4 == 1]
print("A1 primes p = 1 mod 4 in [5, 3000): %d  (explorer: 211)" % len(primes14))
check(len(primes14) == 211, "prime count")
aac = {}
for p in primes14:
    t, u = fundamental_unit_max_order(p)
    check(t * t - p * u * u == -4, "norm at p=%d" % p)
    check(u % p != 0, "AAC violation at p=%d" % p)
    aac[p] = (t, u)
v3h, v2h = {}, {}
for p, (t, u) in aac.items():
    v3h[v_p(u, 3)] = v3h.get(v_p(u, 3), 0) + 1
    v2h[v_p(u, 2)] = v2h.get(v_p(u, 2), 0) + 1
print("    all 211 satisfy t^2 - p u^2 = -4 and p !| u; histogram v3(u): %s; v2(u): %s  (explorer: {0: 211}, {0: 79, 1: 132})"
      % (dict(sorted(v3h.items())), dict(sorted(v2h.items()))))
check(v3h == {0: 211} and v2h == {0: 79, 1: 132}, "valuation histograms")
p1 = [p for p in primes14 if p % 8 == 1]
p5 = [p for p in primes14 if p % 8 == 5]
print("    census: p=1 mod 8: %d (all v2(u)=1: %s); p=5 mod 8: %d, v2(u)=0: %d, v2(u)=1: %d  (explorer: 101, True, 110, 79, 31)"
      % (len(p1), all(v_p(aac[p][1], 2) == 1 for p in p1), len(p5),
         sum(1 for p in p5 if v_p(aac[p][1], 2) == 0), sum(1 for p in p5 if v_p(aac[p][1], 2) == 1)))
check(len(p1) == 101 and len(p5) == 110 and sum(1 for p in p5 if v_p(aac[p][1], 2) == 0) == 79, "mod-8 census")
longest = max(aac, key=lambda p: len(str(aac[p][1])))
print("    longest u: p=%d, %d digits, u mod p = %d  (explorer: p=2689, 38 digits, 2646)"
      % (longest, len(str(aac[longest][1])), aac[longest][1] % longest))
check(longest == 2689 and len(str(aac[2689][1])) == 38 and aac[2689][1] % 2689 == 2646, "longest u")
print("    small cases (p, t, u): %s" % [(p, aac[p][0], aac[p][1]) for p in primes14[:6]])
check(aac[5] == (1, 1) and aac[13] == (3, 1) and aac[17] == (8, 2), "small units (1+sqrt5)/2, (3+sqrt13)/2, 4+sqrt17")
# boundary: the proof of S7b uses t = u mod 2 (maximal order) -- verified on all 211
check(all((t - u) % 2 == 0 for t, u in aac.values()), "t = u mod 2")
print("    S7 / S7b: CONFIRMED (independent solver agrees on every number; t = u mod 2 on all 211).")

# ---------------------------------------------------------------------------
print("\n### A2  GIUGA: direct sum (pure Python) n <= 5000; criteria by sympy.factorint n <= 10^5")
ND = 5000
bad = 0
comp_hits = []
for n in range(2, ND + 1):
    s = sum(pow(i, n - 1, n) for i in range(1, n)) % n
    if sympy.isprime(n):
        if s != n - 1:
            bad += 1
    elif s == n - 1:
        comp_hits.append(n)
print("A2 direct sum n <= %d: primes failing Fermat: %d; composites with sum = -1: %s  (explorer: none up to 30000)" % (ND, bad, comp_hits or "none"))
check(bad == 0 and not comp_hits, "direct Giuga sum")
N1 = 100000
carm, giuga_nums, bbbg, user, disagree = [], [], [], [], 0
for n in range(2, N1 + 1):
    f = sympy.factorint(n)
    ps = list(f)
    sqf = all(a == 1 for a in f.values())
    isp = (len(ps) == 1 and f[ps[0]] == 1)
    kor = sqf and all((n - 1) % (p - 1) == 0 for p in ps)
    giu = sqf and all((n // p - 1) % p == 0 for p in ps)
    bb = sqf and all(((n // p - 1) % (p - 1) == 0) and ((n // p - 1) % p == 0) for p in ps)
    us = all((n - p) % (p * p * (p - 1)) == 0 for p in ps)
    if bb != us:
        disagree += 1
    if not isp:
        if kor and len(ps) >= 3:
            carm.append(n)
        if giu:
            giuga_nums.append(n)
        if bb:
            bbbg.append(n)
        if us:
            user.append(n)
print("    Carmichael numbers <= 10^5: %d %s" % (len(carm), carm))
print("    Giuga numbers <= 10^5: %s; composites satisfying BBBG: %s; user's criterion: %s; disagreements (all n, primes included): %d"
      % (giuga_nums, bbbg or "none", user or "none", disagree))
check(carm == [561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341, 41041, 46657, 52633, 62745, 63973, 75361], "Carmichael list")
check(giuga_nums == [30, 858, 1722, 66198] and not bbbg and not user and disagree == 0, "Giuga / BBBG / user lists")
for g in giuga_nums:
    val = sum(Fraction(1, p) for p in sympy.factorint(g)) - Fraction(1, g)
    check(val == 1, "Giuga sum at %d" % g)
print("    sum 1/p - 1/n = 1 for each Giuga number: True")
# S2 lemma boundary cases: p = 2 with a >= 2 and odd k; a >= 2 odd p with k = -1 mod p (the actual case k = n-1)
lem_bad = 0
for p in list(sympy.primerange(2, 40)):
    for a in range(2, 7):
        q = p ** a
        if q > 5000:
            break
        for k in (q - 1, 2 * q - 1, 3 * q - 1):   # k = -1 mod q, which is what n-1 is when q | n
            s = sum(pow(i, k, q) for i in range(1, q)) % q
            if s % (p ** (a - 1)) != 0 or s == q - 1:
                lem_bad += 1
print("    S2 prime-power lemma at k = -1 mod q (the case that occurs), q = p^a <= 5000, a >= 2: failures %d" % lem_bad)
check(lem_bad == 0, "prime-power lemma")
print("    S1-S4: CONFIRMED.  Note: BBBG's equivalence is already in Giuga (1950); BBBG restate and extend it.")

# ---------------------------------------------------------------------------
print("\n### A3  the trunk: R^j(1) vs (4^j-1)/3")


def R(x):
    return 4 * x + 1


orb1 = [1]
orb0 = [0]
for _ in range(8):
    orb1.append(R(orb1[-1]))
    orb0.append(R(orb0[-1]))
tj = [(4 ** j - 1) // 3 for j in range(9)]
print("A3 R^j(1), j=0..8: %s" % orb1)
print("    R^j(0), j=0..8: %s" % orb0)
print("    t_j=(4^j-1)/3, j=0..8: %s" % tj)
check(orb0 == tj and orb1[:-1] == tj[1:], "trunk identification")
print("    => t_j = R^j(0) = R^(j-1)(1), NOT R^j(1) (explorer's S5 header is off by one; wild_typing lane has N_j = R^(j-1)(1)).")
trunk_rows = {}
for j in range(3, 41):
    t = (4 ** j - 1) // 3
    f = sympy.factorint(t)
    kor = all((t - 1) % (p - 1) == 0 for p in f)
    giu = all((t // p - 1) % p == 0 for p in f)
    fail = next((p for p in sorted(f) if (t - 1) % (p - 1) != 0 or (t // p - 1) % p != 0), None)
    trunk_rows[j] = (all(a == 1 for a in f.values()), kor, giu, fail)
    check(not kor and not giu, "trunk Korselt/Giuga at j=%d" % j)
print("    j: (squarefree, Korselt, Giuga, first failing prime) for j in 3..12, 20, 30, 40:")
for j in list(range(3, 13)) + [20, 30, 40]:
    print("      %2d: %s" % (j, trunk_rows[j]))
expect_fail = {3: 7, 4: 5, 5: 11, 6: 3, 7: 43, 8: 5, 9: 3, 10: 5, 11: 23, 12: 7, 20: 5, 30: 5, 40: 5}
check(all(trunk_rows[j][3] == v for j, v in expect_fail.items()), "first failing primes")
check(trunk_rows[9][0] is False and trunk_rows[10][0] is False and trunk_rows[11][0] is True, "squarefree flags")
print("    S5 table: CONFIRMED (0 Korselt, 0 Giuga for 3 <= j <= 40; all 13 first-failing primes agree).")

# ---------------------------------------------------------------------------
print("\n### A4  continued fraction of log_2 3 by mpmath (independent of PARI)")
mp.dps = 220
alpha = log(mpf(3)) / log(mpf(2))
x = alpha
cf = []
for _ in range(60):
    a = int(floor(x))
    cf.append(a)
    x = 1 / (x - a)
print("A4 mpmath at %d digits; first 60 partial quotients: %s" % (mp.dps, cf))
expl_cf = [1, 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, 15, 1, 9, 2, 5, 7, 1, 1, 4, 8, 1, 11, 1, 20, 2, 1, 10, 1, 4, 1, 1, 1, 1, 1, 37, 4, 55, 1, 1, 49, 1, 1, 1, 4, 1, 3, 2, 3, 3, 1]
check(cf == expl_cf, "continued fraction")
rec, m = [], 0
for i, a in enumerate(cf):
    if a > m:
        m, rec = a, rec + [(i, a)]
print("    records: %s; max %d; sum %d; 55 recurs at index %d" % (rec, max(cf), sum(cf), cf.index(55, 15)))
check(sum(cf) == 388 and rec == [(0, 1), (3, 2), (5, 3), (7, 5), (9, 23), (14, 55)] and cf.index(55, 15) == 46, "records/sum")
h1, h0, k1, k0 = 1, 0, 0, 1
conv = []
for a in cf[:40]:
    h1, h0 = a * h1 + h0, h1
    k1, k0 = a * k1 + k0, k1
    conv.append((h1, k1))
tab = [(h, k, float(k * fabs(k * alpha - h))) for h, k in conv if k <= 10 ** 12]
print("    convergents (q <= 10^12), q||q alpha||: %s" % ["%d/%d:%.6f" % r for r in tab])
expl_q = {1: 0.584963, 5: 0.375937, 665: 0.041881, 190537: 0.017732, 10590737: 0.799110, 753110839881: 0.129252}
for h, k, val in tab:
    if k in expl_q and k != 1:
        check(abs(val - expl_q[k]) < 2e-6, "q||q alpha|| at q=%d" % k)
check(len(tab) == 25, "25 convergents with q <= 10^12")
print("    S10: CONFIRMED (60 quotients, records, sum, all 25 convergent qualities agree to 1e-6).")

# ---------------------------------------------------------------------------
print("\n### A5  S12 float64 minimum re-checked at high precision")
phi = (1 + sqrt(mpf(5))) / 2


def dist(z):
    return fabs(z - nint(z))


n = 10946
val = n * dist(n * alpha) * dist(n * phi)
print("A5 n=10946=F_21: n||n alpha|| ||n phi|| = %.6f (explorer float64: 0.000209); n||n alpha|| = %.6f, n||n phi|| = %.6f"
      % (float(val), float(n * dist(n * alpha)), float(n * dist(n * phi))))
check(abs(float(val) - 0.000209) < 2e-6, "S12 minimum")
q = 190537
print("    q=190537: exact q||q alpha|| = %.6f; the float64 scan reported 0.017734 (error ~ q^2 * 2^-53 = %.1e): a float64 artefact,"
      % (float(q * dist(q * alpha)), q * q * 2 ** -53))
print("    consistent with the HEURISTIC label.  min_n n||n phi|| over n<=10^6 is 1/phi^2 = %.6f at n=1; liminf is 1/sqrt5 = %.6f."
      % (float(1 / phi ** 2), float(1 / sqrt(mpf(5)))))
fibs = [1, 1]
while fibs[-1] < 20000:
    fibs.append(fibs[-1] + fibs[-2])
check(fibs[20] == 10946, "F_21 = 10946 (1-based)")

# ---------------------------------------------------------------------------
print("\n### A6  S15: the index slip in the Mobius-correlated 2-adic point")


def T(x):
    return (3 * x + 1) // 2 if x & 1 else x // 2


def parities(x, k):
    out = []
    for _ in range(k):
        out.append(x & 1)
        x = T(x)
    return out


def lift(target):
    """2-adic truncation x mod 2^K whose first K parities are `target` (unique)."""
    x = 0
    for k in range(1, len(target) + 1):
        if parities(x, k)[k - 1] != target[k - 1]:
            x += 2 ** (k - 1)
        check(parities(x, k) == target[:k], "lift at k=%d" % k)
    return x


mu = [0] * 2101
for n in range(1, 2101):
    mu[n] = int(sympy.mobius(n))
K = 32
target_expl = [1 if mu[k + 1] == 1 else 0 for k in range(K)]            # parity(T^k x) = [mu(k+1)=1], k>=0
target_corr = [0] + [1 if mu[k] == 1 else 0 for k in range(1, K)]       # parity(x)=0, parity(T^k x) = [mu(k)=1], k>=1
x_expl = lift(target_expl)
x_corr = lift(target_corr)
print("A6 explorer's target [mu(k+1)=1], k=0..31: %s -> x = %d mod 2^32 (explorer: 890931253)" % ("".join(map(str, target_expl)), x_expl))
check(x_expl == 890931253, "explorer's point")
print("    corrected target (parity(x)=0; [mu(k)=1] for k=1..31): %s -> x = %d mod 2^32" % ("".join(map(str, target_corr)), x_corr))
print("    With the explorer's indexing f(T^n x) = parity(T^n x) = [mu(n+1)=1], so (1/N) sum mu(n) f(T^n x) = (1/N) sum mu(n)[mu(n+1)=1],")
print("    NOT (1/N)#{mu(n)=1}.  Demonstration on 2000-bit truncations, N = 1999:")
K2 = 2000
te = [1 if mu[k + 1] == 1 else 0 for k in range(K2)]
tc = [0] + [1 if mu[k] == 1 else 0 for k in range(1, K2)]
xe = lift(te)
xc = lift(tc)
pe = parities(xe, K2)
pc = parities(xc, K2)
N = K2 - 1
corr_e = sum(mu[n] * pe[n] for n in range(1, N + 1)) / N
corr_c = sum(mu[n] * pc[n] for n in range(1, N + 1)) / N
frac1 = sum(1 for n in range(1, N + 1) if mu[n] == 1) / N
print("    explorer's point: (1/N) sum mu(n) parity(T^n x) = %+.6f;  corrected point: %+.6f;  #{mu(n)=1}/N = %.6f;  3/pi^2 = %.6f"
      % (corr_e, corr_c, frac1, 3 / math.pi ** 2))
check(abs(corr_c - frac1) < 1e-12, "corrected point gives #{mu=1}/N exactly")
check(abs(corr_e - frac1) > 0.1, "explorer's indexing does not give #{mu=1}/N")
print("    (Trivial repair, and the conclusion 'disjointness is FALSE on Z_2' survives with the corrected point; the note's")
print("    formula (1/N)#{mu(n)=1} was wrong for the point as defined.  The correlation with the explorer's point is")
print("    (1/2) sum mu(n) mu(n+1)^2 + (1/2) sum mu(n) mu(n+1), whose second term is the 2-point Chowla sum: OPEN for natural density.)")
print("    A 2-adic point with parity vector [mu=1] is not in Q: rational x have eventually periodic parity vectors (Lagarias 1985, CITED)")
print("    and [mu=1] is not eventually periodic (every class r mod P contains infinitely many n with mu(n)=1 and with mu(n)=0).")

# ---------------------------------------------------------------------------
print("\n### A7  S16 correlations by an independent linear sieve and direct iteration, N = 10^6")
N4 = 10 ** 6
mu6 = [1] * (N4 + 1)
lp = [0] * (N4 + 1)
primes = []
mu6[0] = 0
for i in range(2, N4 + 1):
    if lp[i] == 0:
        lp[i] = i
        primes.append(i)
        mu6[i] = -1
    for p in primes:
        if p > lp[i] or i * p > N4:
            break
        lp[i * p] = p
        mu6[i * p] = 0 if p == lp[i] else -mu6[i]
sqf_count = sum(1 for n in range(1, N4 + 1) if mu6[n] != 0)
print("A7 linear sieve: squarefree count to 10^6 = %d (expected 607926)" % sqf_count)
check(sqf_count == 607926, "squarefree count to 10^6 (607926, OEIS-known)")
cur = list(range(N4 + 1))
expl_c = [-0.000068, 0.000344, 0.000306, 0.000086, -0.000054, 0.000600, -0.000760]
for j in range(7):
    par = [c & 1 for c in cur]
    c = sum(mu6[n] * (1 - 2 * par[n]) for n in range(1, N4 + 1)) / N4
    mod = 2 ** (j + 1)
    periodic = all(par[r] == par[r + mod * s] for r in range(mod) for s in range(1, 8))
    periodic_full = all(par[n] == par[n % mod] for n in range(N4 + 1))
    half = mod // 2
    periodic_half = all(par[n] == par[n % half] for n in range(N4 + 1)) if half >= 1 else True
    mean_par = sum(par[1:]) / N4
    print("    j=%d  c_j = %+.6f (explorer %+.6f)  parity mean %.6f  periodic mod %d: %s  periodic mod %d: %s"
          % (j, c, expl_c[j], mean_par, mod, periodic_full, half, periodic_half))
    check(abs(c - expl_c[j]) < 1.5e-6, "c_%d" % j)
    check(periodic_full and not periodic_half, "periodicity sharpness at j=%d" % j)
    cur = [T(v) for v in cur]
print("    S16: CONFIRMED; sharpness: parity(T^j n) is periodic mod 2^(j+1) and NOT mod 2^j (Terras 1976 / Everett 1977, CITED).")

# ---------------------------------------------------------------------------
print("\n### A8  S20 / S21 coprimality")
g1 = max(math.gcd(4 ** j * 7 + (4 ** j - 1) // 3, 4 ** (j + 1) * 7 + (4 ** (j + 1) - 1) // 3) for j in range(30))
gmax, pairs = 0, 0
for n0 in range(1, 20001, 2):
    xx = n0
    for _ in range(200):
        y = 3 * xx + 1
        while y % 2 == 0:
            y //= 2
        if y == xx:
            break
        gmax = max(gmax, math.gcd(xx, y))
        pairs += 1
        xx = y
        if xx == 1:
            break
print("A8 max gcd(x_j, x_{j+1}) on the R-orbit of 7 (j<30): %d; consecutive odd iterates: %d pairs, max gcd %d (explorer: 326623, 1)" % (g1, pairs, gmax))
check(g1 == 1 and pairs == 326623 and gmax == 1, "coprimality counts")
print("    S20/S21: CONFIRMED (the pair loop stops at the fixed point 1 = (3*1+1)/4 and at x=1; the fixed point is excluded from the count).")

# ---------------------------------------------------------------------------
print("\n### A9  repository greps re-done")
thm_dir = os.path.join(REPO, "01-canon", "theorems")
tourn = []
for fn in sorted(os.listdir(thm_dir)):
    if not fn.startswith("THM-"):
        continue
    with open(os.path.join(thm_dir, fn), encoding="utf-8", errors="replace") as fh:
        body = fh.read()
    m = re.search(r"^title:\s*(.*)$", body[:4000], re.M)
    title = m.group(1) if m else fn
    if re.search(r"tournament", title, re.I):
        tourn.append((fn, title, body))
with17 = [fn for fn, t, b in tourn if re.search(r"\b17\b|seventeen", t, re.I)]
body17 = [fn for fn, t, b in tourn if re.search(r"17[- ]vert|\bn\s*=\s*17\b|order[- ]17\b|17-tournament|17 vertices", b, re.I)]
print("A9 tournament-titled canon theorems: %d (explorer 182); with 17 in title: %d (explorer 10); body hits for n=17 tournaments: %d %s"
      % (len(tourn), len(with17), len(body17), body17))
check(len(tourn) == 182 and len(with17) == 10 and body17 == ["THM-868-e8-bridge-score-lattice.md", "THM-871-fermat-rung-rigidity.md"], "grep counts")
print("    => the explorer's NOTE said 0 body hits; its own .out said 2.  THM-871 is ABOUT 17-vertex (rotational) tournaments:")
for fn, t, b in tourn:
    if fn.startswith("THM-871"):
        for ln in b.splitlines():
            if re.search(r"n = 17", ln) and not ln.startswith("title:"):
                print("      THM-871 body: " + ln.strip()[:150])
                break
        print("      THM-871 title: " + t[:260])
res_dir = os.path.join(REPO, "05-knowledge", "results")
pat = re.compile(r"S\^6\b|\bS6 monodromy|S_2 x S_3|S2 x S3")
hits = []
for fn in sorted(os.listdir(res_dir)):
    if not fn.endswith(".md") or fn.startswith("collatz_mod6_20260921_grand_circuit_typing"):
        continue
    with open(os.path.join(res_dir, fn), encoding="utf-8", errors="replace") as fh:
        txt = fh.read()
    if pat.search(txt):
        hits.append(fn)
print("    results notes (own lane excluded) matching S^6 / S6 monodromy / S_2 x S_3: %d: %s" % (len(hits), hits))
print("    classification: planar_jc48_sep06_* : 'S^6' is the monomial S^6 in projective coordinates [S:T] (not a sphere, not a paste audit);")
print("      arithmetic_seams_20260921_{synthesis,dynamics}: an explicit S_2 x S_3 action on a six-vector lift of a rational three-point cycle")
print("      (a concurrent opus lane on the same paste family); collatz_blueprint_20260921_synthesis: the S6-monodromy refutation.")
check(set(hits) == {"arithmetic_seams_20260921_dynamics.md", "arithmetic_seams_20260921_synthesis.md", "collatz_blueprint_20260921_synthesis.md",
                    "planar_jc48_sep06_global_curve.md", "planar_jc48_sep06_resolution_budget.md"}, "S^6 hit set")
print("    => the explorer's 'S_2 x S_3 appears once (catalan_elliptic)' is wrong: the catalan_elliptic synthesis writes 'C2 x S3' (not matched")
print("      by its own regex) and the arithmetic_seams lanes carry a defined S_2 x S_3 action.  Also the count 5/6 is unstable (it counted itself).")
audited = ["collatz_mod6_20260917_%s.md" % s for s in ("synthesis", "row_braid_typing", "extended_collatz_scc", "three_adic_g_map", "reverse_tree_pieces",
           "zsigmondy_triad", "pythagorean_semicircle", "sandwich_bias", "cell_ordering_scale", "divisor_balance_family", "floor_sums_odd_functions",
           "wild_typing", "scaffolding_audit")]
audited += [fn for fn in os.listdir(res_dir) if re.match(r"(collatz_blueprint|collatz_guards|catalan_elliptic|glued_xor|odd_square|prime_shells|arithmetic_seams|arithmetic_braids2?)_2026", fn) and fn.endswith(".md")]
wy = {}
for fn in sorted(set(audited)):
    path = os.path.join(res_dir, fn)
    if not os.path.exists(path):
        continue
    with open(path, encoding="utf-8", errors="replace") as fh:
        txt = fh.read()
    found = [w for w in ("wythoff", "zeckendorf", "golden ratio") if w in txt.lower()]
    if found:
        wy[fn] = found
print("    Wythoff / Zeckendorf / 'golden ratio' in the audited notes: %s" % wy)
check("collatz_blueprint_20260921_energy.md" in wy and "collatz_mod6_20260917_scaffolding_audit.md" in wy, "Wythoff hits")
print("    => the explorer's 'Wythoff, Zeckendorf and the golden ratio occur in none of the audited notes' is FALSE as stated:")
print("      the blueprint energy note defines the lower Wythoff sequence W(k)=floor(k phi) and the scaffolding audit / mod6 synthesis")
print("      carry 'Wythoff diffraction' as a SCOPE token; the guards valves note mentions the golden ratio.  The correct statement is that")
print("      no audited note attaches a golden-ratio object to a Collatz orbit (all occurrences are SCOPE typings of the same paste family).")

# ---------------------------------------------------------------------------
print("\n### A10  S25 table: 2^K - 3^L along the convergents")
rows = []
for h, k in conv[:14]:
    d = 2 ** h - 3 ** k
    rows.append((h, k, 1 if d > 0 else -1, abs(d) / 3 ** k, len(str(abs(d)))))
    print("    K=%6d L=%6d sign %+d  |2^K-3^L|/3^L = %.3e  (%d digits)" % rows[-1])
expl = [(1, 1, -1, 1), (2, 1, 1, 1), (3, 2, -1, 1), (8, 5, 1, 2), (19, 12, -1, 4), (65, 41, 1, 18), (84, 53, -1, 23), (485, 306, 1, 144),
        (1054, 665, -1, 313), (24727, 15601, 1, 7439), (50508, 31867, -1, 15200), (125743, 79335, 1, 37847), (176251, 111202, -1, 53052),
        (301994, 190537, 1, 90903)]
check([(r[0], r[1], r[2], r[4]) for r in rows] == expl, "S25 table")
check(abs(rows[3][3] - 13 / 243) < 1e-12 and abs(rows[13][3] - 6.451e-8) < 1e-10, "S25 ratios")
print("    S25: CONFIRMED (signs alternate; 2^8-3^5 = 13 = 5.350e-02 * 3^5).")

print("\n### AUDIT SUMMARY")
print("CONFIRMED: S1-S4 (Giuga, BBBG, user's criterion), S7/S7b (AAC data and forced valuations), S10-S12 (cf of log_2 3; S12 float artefact noted),")
print("           S16 (c_j and periodicity, sharp), S18-S21, S25, the 182/10 title greps.")
print("REFUTED-AS-WRITTEN and repaired: S5 header t_j = R^j(1) (is R^(j-1)(1) = R^j(0)); S15 correlation formula (index slip mu(n) vs mu(n+1));")
print("           S22 '0 body hits' (2: THM-868, THM-871 -- THM-871 is a theorem about 17-vertex rotational tournaments); S24 'S_2 x S_3 appears once'")
print("           (arithmetic_seams lanes carry a defined action; the grep count was unstable and included non-paste planar_jc48 monomials);")
print("           S13 'Wythoff in none of the audited notes' (blueprint energy defines W(k)=floor(k phi)); S23 h-spectrum cited to THM-1745 (it is THM-1370,")
print("           whose completeness 'odds minus {7,21}' is a CONJECTURE there).")
print("DOWNGRADES: 'Stroeker-Tijdeman 1982' exponent not asserted (shape only); the norm -1 fact is CITED classical (Legendre); Terras/Everett added for S16.")
print("ALL AUDIT CHECKS PASSED")

#!/usr/bin/env python3
"""Audit D: constants of THM-4654 (a)-(c), the closed forms in its Setting, and the arithmetic of the Chowla ledger."""
import math
from itertools import product
from fractions import Fraction as F
import mpmath

mpmath.mp.dps = 50
OUT = '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_D_20261014/audit_constants.out'
L = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    L.append(s)


lg3 = mpmath.log(3, 2)
c = 1 - lg3 / 3
say("[THM-4654 constants]")
say("  lambda = 3^(1/3)/2 =", mpmath.nstr(mpmath.cbrt(3) / 2, 12), "; c = 1 - log2(3)/3 =", mpmath.nstr(c, 12),
    "; 2^(-c) =", mpmath.nstr(mpmath.power(2, -c), 12))
say("  1/c =", mpmath.nstr(1 / c, 12), " (THM: 'about 2.12'; the proof of (c) uses 2.12 as an UPPER bound - it is not)")
say("  log2(3)/c =", mpmath.nstr(lg3 / c, 12), " (THM: 3.36, used as an upper bound in (a) and (c) - rounded DOWN)")
say("  log2(3)*2.1 =", mpmath.nstr(lg3 * 2.1, 12), " (THM: 3.33, rounded up: fine)")
say("  2*log2(3) =", mpmath.nstr(2 * lg3, 12), " (additive constant of (a); THM: 3.2: fine)")
# (a): does log2 3 * (n1 + floor(2.1v) + 1) <= 3.33 v + 3.36 log2(6|u|/v + 1) + 3.2 always hold?  v = 1, |u| = 2^k
first_fail = None
for k in range(1, 20000):
    X = mpmath.mpf(k) + mpmath.log(6, 2)            # log2(6|u|) with |u| = 2^k, v = 1
    n1 = int(mpmath.ceil(X / c))
    lhs = lg3 * (n1 + 3)
    rhs = mpmath.mpf('3.33') + mpmath.mpf('3.36') * mpmath.log(mpmath.power(2, X) + 1, 2) + mpmath.mpf('3.2')
    if lhs > rhs:
        first_fail = k
        break
say(f"  (a) as stated is NOT implied by the proved bound once |u| >= 2^{first_fail} (v = 1): the proved term "
    f"log2(3)*ceil(log2(6|u|)/c) exceeds 3.36*log2(6|u|+1) + slack")
# (c): 3.36 (K+2) vs (log2 3 / c)(K+1)
Kx = (2 * mpmath.mpf('3.36') - lg3 / c) / (lg3 / c - mpmath.mpf('3.36'))
say(f"  (c) derivable exponent mu <= (log2 3/c)(K+1) = {mpmath.nstr(lg3 / c, 8)}(K+1); the stated 3.36(K+2) is >= it only "
    f"for K <= {mpmath.nstr(Kx, 8)}")
say(f"  Terras-memory version: mu <= log2(3)(K+1) = {mpmath.nstr(lg3, 8)}(K+1) for K >= 1/c - 1 = {mpmath.nstr(1 / c - 1, 6)}"
    f"  (e.g. Sturmian K = 2: {mpmath.nstr(3 * lg3, 6)} instead of {mpmath.nstr(mpmath.mpf('3.36') * 4, 6)})")
say("  statement 6 equality-family gap log2(27/19) =", mpmath.nstr(mpmath.log(mpmath.mpf(27) / 19, 2), 10))
say("  sup offset 1/2 + (3/8)/(1-3/8) =", F(1, 2) + F(3, 8) / (1 - F(3, 8)))


# closed forms in the Setting: A = -(1/3) sum (2/3)^(n_m) (2-adic) and z_A = -sum 2^(n_m) 3^(-m)
def mahler_code(A, n):
    out, a = [], A
    for _ in range(n):
        out.append(a & 1)
        a = (3 * a + (a & 1)) // 2
    return out


def tcode_2adic(z_num, z_den, n):
    """T-code of the 2-adic rational z_num/z_den (odd den)."""
    out = []
    j, v = z_num, z_den
    for _ in range(n):
        b = j & 1
        out.append(b)
        j = (3 * j + v) // 2 if b else j // 2
    return out


ok = True
for A in (1, 2, 3, 13, 100, 8101, 3532045):
    H = 60
    e = mahler_code(A, H)
    mod = 2 ** H
    inv3 = pow(3, -1, mod)
    # A == -(1/3) sum_m (2/3)^(n_m)  mod 2^H  (terms with n_m >= H vanish mod 2^H)
    s = 0
    for n, b in enumerate(e):
        if b:
            s = (s + pow(2, n, mod) * pow(inv3, n, mod)) % mod
    ok &= ((-s * inv3) % mod == A % mod)
    # z_A = -sum_m 2^(n_m) 3^(-m) mod 2^H, and its T-code equals e on [0, H)
    z = 0
    m = 0
    for n, b in enumerate(e):
        if b:
            m += 1
            z = (z + pow(2, n, mod) * pow(inv3, m, mod)) % mod
    z = (-z) % mod
    ok &= (tcode_2adic(z, 1, H) == e)
say(f"[Setting] A = -(1/3) sum (2/3)^(n_m) and z_A = -sum 2^(n_m) 3^(-m) reproduce A and the code (mod 2^60) "
    f"for 7 test integers: {ok}")

# ---------------------------------------------------------------- Chowla ledger arithmetic
say("[Ledger]")
delta = mpmath.mpf(1) / 200
for W in (10,):
    lhs = delta * (2 + mpmath.log(2)) / (6 * W)
    say(f"  delta(2+log2)/(6W) at W={W}: {mpmath.nstr(lhs, 8)}; (1/20)/that = {mpmath.nstr(mpmath.mpf(1) / 20 / lhs, 8)}")
    say(f"  delta(1+log2)/(6W) = {mpmath.nstr(delta * (1 + mpmath.log(2)) / (6 * W), 8)} (< 100 trivially)")
say(f"  c = delta/(6WA) at W=10, A=1000: 1/{mpmath.nstr(1 / (delta / (6 * 10 * 1000)), 10)}")
say(f"  A/700 >= 0.8 iff A >= {0.8 * 700:.0f};  1000/700 = {F(1000, 700)};  theta <= 1000/0.8 = {1000 / 0.8:.0f};"
    f"  2^(1/700) = {mpmath.nstr(mpmath.power(2, mpmath.mpf(1) / 700), 10)}")
say(f"  large-|B| exponent 4 tau - (3.96 + 0.99 + 0.8): at tau = 21/20 -> {4 * 1.05 - 5.75:.2f} (= -31/20 = {-31 / 20})")
say(f"  balance 4 tau - 5.75 = -tau -> tau* = {F(575, 500)} = 23/20; output exponent of the lemma = min(tau, 5.75 - 4 tau);"
    f" usable iff min(...) > 1 + 0.000224: tau in (1.000224, {(5.75 - 1.000224) / 4:.6f})")
say(f"  6/5 = (8/5)(3/4): {F(8, 5) * F(3, 4)};  need (8/5) alpha_Q > 1 <=> alpha_Q > 5/8")
from sympy import totient
good = [m for m in range(2, 40) if 1 - F(1, int(totient(m))) > F(5, 8)]
say(f"  moduli m < 40 with density 1 - 1/phi(m) of non-1 classes > 5/8: {good}  (smallest = {good[0]})")
# U_(3,4): direct count of points of F_p^3 avoiding x1, x2, x3, x1 + x2 - x3
res = []
for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
    cnt = 0
    for x1, x2, x3 in product(range(p), repeat=3):
        if x1 % p and x2 % p and x3 % p and (x1 + x2 - x3) % p:
            cnt += 1
    chi = p ** 3 - 4 * p ** 2 + 6 * p - 3
    res.append((p, cnt == chi, cnt == (p - 1) * (p * p - 3 * p + 3)))
    # nu_4(p) = 1 - cnt/p^3 = 4/p - 6/p^2 + 3/p^3
    assert F(p ** 3 - cnt, p ** 3) == F(4, p) - F(6, p ** 2) + F(3, p ** 3)
say(f"  U(3,4) avoiding count = t^3-4t^2+6t-3 = (t-1)(t^2-3t+3) at p in {[r[0] for r in res]}: "
    f"{all(r[1] and r[2] for r in res)}; nu_4(p) = 4/p - 6/p^2 + 3/p^3 exactly; chi(2) = {8 - 16 + 12 - 3}")
# every 3 of the 4 forms independent over every F_p (det of each 3x3 = +-1)
import itertools
forms = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, -1)]
dets = []
for tri in itertools.combinations(forms, 3):
    a, b, cc = tri
    det = (a[0] * (b[1] * cc[2] - b[2] * cc[1]) - a[1] * (b[0] * cc[2] - b[2] * cc[0]) + a[2] * (b[0] * cc[1] - b[1] * cc[0]))
    dets.append(det)
say(f"  3x3 minors of the four forms: {dets} (all +-1, so U(3,4) over every prime field incl. F_2)")
# Part II: rho constraints found in the paper (eps = 1e-4)
eps = mpmath.mpf('1e-4')
say("[Ledger Part II, rho]")
say(f"  ledger: needs only rho > 4 eps = {mpmath.nstr(4 * eps, 4)}, 'margin 125x' = {mpmath.nstr(mpmath.mpf(1) / 20 / (4 * eps), 6)}")
say(f"  Lemma 16.3(ii): 1 - rho + 4 eps < 1 - rho/2  <=>  rho > 8 eps = {mpmath.nstr(8 * eps, 4)}  (margin "
    f"{mpmath.nstr(mpmath.mpf(1) / 20 / (8 * eps), 6)}x)")
say(f"  Lemma 14.6: 1 - 8 rho > 2 eps  <=>  rho < (1 - 2 eps)/8 = {mpmath.nstr((1 - 2 * eps) / 8, 8)}  (margin "
    f"{mpmath.nstr((1 - 2 * eps) / 8 * 20, 6)}x above 1/20); at rho = 1/20: 1 - 8 rho = {1 - 8 * F(1, 20)}")
say(f"  Lemma 14.5: (1 - 4 rho)/8 > 2 eps  <=>  rho < (1 - 16 eps)/4 = {mpmath.nstr((1 - 16 * eps) / 4, 8)}; at rho = 1/20: "
    f"(1 - 4 rho)/8 = {(1 - 4 * F(1, 20)) / 8}")
with open(OUT, 'w') as fh:
    fh.write("\n".join(L) + "\n")

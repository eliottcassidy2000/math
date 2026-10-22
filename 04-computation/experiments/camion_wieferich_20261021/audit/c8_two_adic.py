# c8_two_adic.py -- second-pass audit of THM-4750 6(b), 6(c) ("PROVED in outline").
#  (1) Newton polygon of T(x) = sum 2^C(n,2) x^n/n! over Q_2 (exact valuations, n <= 60): segments, zero valuations.
#  (2) The unit root rho = zeta*w of T in Z_2[zeta_3] (Newton iteration mod 2^K), C = 1/(rho T'(rho)), r = sigma(w)/w,
#      t = -zeta sigma(C)/C, and the valuations v(w^3 - 1), v(C), v(r - 1), v(t - 1) used by the outline.
#  (3) Remainder decay: v_2(f(n) - Tr(C rho^-n)) - 3n/2 for n <= 120 (the "O(1)" in v(r_n) >= 3n/2 - O(1)).
#  (4) Range of the quasi-period 6(b): for j <= JMAX list the n <= NMAX with f(n + 3*2^j) != f(n) (mod 2^(j+2)), and on
#      n = 1 (mod 3) the n with failure mod 2^(j+3).  The outline proves the congruence only once 3n/2 - O(1) >= j + 2.
# f(n) = s(n)/n! = -[x^n] 1/T(x) (n >= 1), computed mod 2^K from T mod 2^K (a polynomial: v(tau_n) = C(n,2) - v_2(n!)).
# Second-pass auditor, 2026-10-21.
import sys
from fractions import Fraction
K = 160; MOD = 1 << K
JMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16; NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 60
def v2(x):
    if x == 0: return 10 ** 9
    x = abs(x); v = 0
    while x % 2 == 0: x //= 2; v += 1
    return v
def v2fact(n): return n - bin(n).count('1')
# ---- (1) Newton polygon
pts = [(n, n * (n - 1) // 2 - v2fact(n)) for n in range(0, 61)]
hull = [pts[0]]
i = 0
while i < len(pts) - 1:
    best = None
    for j in range(i + 1, len(pts)):
        sl = Fraction(pts[j][1] - pts[i][1], pts[j][0] - pts[i][0])
        if best is None or sl < best[0] or (sl == best[0] and j > best[1]): best = (sl, j)
    hull.append(pts[best[1]]); i = best[1]
segs = [(Fraction(b[1] - a[1], b[0] - a[0]), b[0] - a[0]) for a, b in zip(hull, hull[1:])]
print("(1) Newton polygon of T over Q_2 (first segments): (slope s, number of zeros; these zeros have v_2(x) = -s, |x| = 2^s):", [(str(s), l) for s, l in segs[:6]])
print("    zeros in the closed unit disc (slope <= 0):", sum(l for s, l in segs if s <= 0), "; next slope:", str(segs[1][0]), "-> |x| =", f"2^{segs[1][0]}")
# ---- arithmetic in Z/2^K [zeta], zeta^2 = -1 - zeta ; element (a, b) = a + b zeta
def mul(x, y):
    a, b = x; c, d = y                  # (a + b z)(c + d z) = ac + (ad + bc) z + bd z^2 = (ac - bd) + (ad + bc - bd) z
    return ((a * c - b * d) % MOD, (a * d + b * c - b * d) % MOD)
def add(x, y): return ((x[0] + y[0]) % MOD, (x[1] + y[1]) % MOD)
def sub(x, y): return ((x[0] - y[0]) % MOD, (x[1] - y[1]) % MOD)
def sc(c, x): return (c * x[0] % MOD, c * x[1] % MOD)
def sigma(x): a, b = x; return ((a - b) % MOD, (-b) % MOD)          # zeta -> zeta^2 = -1 - zeta
def norm(x): y = mul(x, sigma(x)); assert y[1] == 0; return y[0]
def inv(x):
    nx = norm(x); assert nx % 2 == 1
    return sc(pow(nx, -1, MOD), sigma(x))
def v(x): return min(v2(x[0]), v2(x[1]))
ONE = (1, 0); ZETA = (0, 1)
# T mod 2^K as a polynomial (2-integral coefficients)
def tau(n):
    e = n * (n - 1) // 2 - v2fact(n)
    if e >= K: return 0
    f = 1
    for k in range(2, n + 1):
        while k % 2 == 0: k //= 2
        f = f * k % MOD
    return (1 << e) * pow(f, -1, MOD) % MOD
deg = 0
while tau(deg + 1) or (deg + 1) * deg // 2 - v2fact(deg + 1) < K: deg += 1
tc = [tau(n) for n in range(deg + 1)]
def T_at(x):
    s = (0, 0); pw = ONE
    for c in tc: s = add(s, sc(c, pw)); pw = mul(pw, x)
    return s
def Tp_at(x):
    s = (0, 0); pw = ONE
    for n in range(1, len(tc)): s = add(s, sc(n * tc[n] % MOD, pw)); pw = mul(pw, x)
    return s
rho = ZETA
for _ in range(12): rho = sub(rho, mul(T_at(rho), inv(Tp_at(rho))))
assert v(T_at(rho)) >= K - 2
w = mul(rho, mul(ZETA, ZETA))                 # rho / zeta = rho * zeta^2
C = inv(mul(rho, Tp_at(rho)))
rr = mul(sigma(w), inv(w)); t = sub((0, 0), mul(ZETA, mul(sigma(C), inv(C))))
w3 = mul(w, mul(w, w))
print(f"(2) rho = zeta*w with v(T(rho)) >= {v(T_at(rho))}: v(w - 1) = {v(sub(w, ONE))}, v(w^3 - 1) = {v(sub(w3, ONE))}, v(C) = {v(C)}, "
      f"v(r - 1) = {v(sub(rr, ONE))}, v(t - 1) = {v(sub(t, ONE))}, N(r) = 1: {norm(rr) == 1}, N(t) = 1: {norm(t) == 1}")
# ---- f(n) mod 2^K
NF = 3 * (1 << JMAX) + NMAX + 5
W = [0] * (NF + 1); W[0] = 1
for n in range(1, NF + 1):
    acc = 0
    for k in range(1, min(n, deg) + 1): acc += tc[k] * W[n - k]
    W[n] = (-acc) % MOD
f = [0] + [(-W[n]) % MOD for n in range(1, NF + 1)]
# cross-check f(n) with exact s(n)/n! for n <= 25
from math import comb, factorial
s = [0] * 26; Tn = [2 ** (n * (n - 1) // 2) for n in range(26)]
for n in range(1, 26): s[n] = Tn[n] - sum(comb(n, k) * s[k] * Tn[n - k] for k in range(1, n))
okf = all((Fraction(s[n], factorial(n)).numerator * pow(Fraction(s[n], factorial(n)).denominator, -1, MOD) - f[n]) % MOD == 0 for n in range(1, 26))
print(f"    f(n) mod 2^{K} from the truncated series agrees with exact s(n)/n! for n <= 25: {okf}")
# ---- (3) remainder decay
rhoinv = inv(rho); sC = sigma(C); srinv = sigma(rhoinv)
x1 = ONE; x2 = ONE; rem = []
for n in range(1, 121):
    x1 = mul(x1, rhoinv); x2 = mul(x2, srinv)
    tr = add(mul(C, x1), mul(sC, x2)); assert tr[1] == 0
    rem.append((n, v2((f[n] - tr[0]) % MOD) if (f[n] - tr[0]) % MOD else K))
print("(3) v_2(r_n) - 3n/2, r_n = f(n) - Tr(C rho^-n), n = 1..40:", [str(Fraction(2 * vv - 3 * n, 2)) for n, vv in rem[:40]])
per = {n % 4: set() for n in range(4)}
for n, vv in rem:
    if vv < K: per[n % 4].add(str(Fraction(2 * vv - 3 * n, 2)))
print("    v_2(r_n) - 3n/2 by n mod 4 (n <= 120):", {k: sorted(x) for k, x in per.items()})
pred_ok = True
for j in range(0, JMAX + 1):
    sh = 3 * (1 << j)
    for n in range(1, NMAX + 1):
        fails = (f[n + sh] - f[n]) % (1 << (j + 2)) != 0
        if fails != (rem[n - 1][1] < j + 2): pred_ok = False
print(f"    rule 'f(n+3*2^j) = f(n) mod 2^(j+2) iff v_2(r_n) >= j+2' on all j <= {JMAX}, n <= {NMAX}: {pred_ok}")
print("    min over 5 <= n <= 100 of v_2(r_n) - 3n/2 =", str(min(Fraction(2 * vv - 3 * n, 2) for n, vv in rem[4:100] if vv < K)))
# ---- (4) range of the quasi-period
print(f"(4) f(n + 3*2^j) = f(n) mod 2^(j+2): failing n in [1, {NMAX}], and (class n = 1 mod 3) failures mod 2^(j+3):")
for j in range(0, JMAX + 1):
    sh = 3 * (1 << j)
    bad = [n for n in range(1, NMAX + 1) if (f[n + sh] - f[n]) % (1 << (j + 2))]
    bad1 = [n for n in range(1, NMAX + 1, 3) if (f[n + sh] - f[n]) % (1 << (j + 3))]
    print(f"    j={j:2d}: fails mod 2^(j+2) at n = {bad};  n = 1 (mod 3) fails mod 2^(j+3) at n = {bad1}")
print("(5) explicit counterexamples to the unqualified 6(b): v_2(f(n + 3*2^j) - f(n)) against j + 2, and v_2(r_n):")
for j, n in ((11, 5), (11, 6), (13, 7), (14, 8), (17, 9), (17, 10)):
    d = (f[n + 3 * (1 << j)] - f[n]) % MOD
    print(f"    j={j}, n={n}: v_2(f(n+3*2^j) - f(n)) = {v2(d)} < j+2 = {j + 2};  v_2(r_n) = {rem[n - 1][1]}")

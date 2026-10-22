# Test of the heuristic in THM-4750 (6c): the 2-adic zero n* of f(n) = s(n)/n! on n = 2 (mod 3) is a ratio of logarithms in the
# norm-one torus of Q_2(zeta_3).  Mod 2^K, 1/T(x) = 1/Q(x) with Q a polynomial, Q = Phi_3 (mod 2).  Hensel-lift the root rho of Q
# near the cube root of unity zeta in W(F_4) = Z_2[zeta]; then f(n) = Tr(C rho^(-n)) + (terms from the roots with |rho| > 1, which
# vanish for large n), with C = 1/(rho Q'(rho)).  Writing rho = zeta*w (w = 1 mod 2), the zero on n = 2 (mod 3) solves
#   r^n = t,   r = sigma(w)/w,   t = -zeta*sigma(C)/C   (both of norm 1).
# The script checks f(n) = Tr(C rho^-n) for large n, and r^N = t (mod 2^60) for N = n* (mod 2^64) from two_adic_zero_fast.out.
# mac-mini-2026-10-21-fanopsp.
K = 120; M = 1 << K
def mul(x, y):                       # (a + b z)(c + d z), z^2 = -1 - z
    a, b = x; c, d = y
    return ((a * c - b * d) % M, (a * d + b * c - b * d) % M)
def add(x, y): return ((x[0] + y[0]) % M, (x[1] + y[1]) % M)
def sigma(x): a, b = x; return ((a - b) % M, (-b) % M)
def power(x, e):
    r = (1, 0)
    while e:
        if e & 1: r = mul(r, x)
        x = mul(x, x); e >>= 1
    return r
def norm(x): return mul(x, sigma(x))[0]
def inv(x):                           # x^-1 = sigma(x)/N(x)
    nx = norm(x); ninv = pow(nx, -1, M); s = sigma(x); return (s[0] * ninv % M, s[1] * ninv % M)
def tau(n):
    f = 1; e = n * (n - 1) // 2
    for i in range(2, n + 1):
        while i % 2 == 0: i //= 2; e -= 1
        f = f * i % M
    return (pow(f, -1, M) << e) % M if e < K else 0
d = 1
while any(tau(j) for j in range(d + 1, d + 6)): d += 1
Q = [tau(j) for j in range(d + 1)]
def Qev(x):
    r = (0, 0)
    for c in reversed(Q): r = add(mul(r, x), (c, 0))
    return r
def Qder(x):
    r = (0, 0)
    for j in range(d, 0, -1): r = add(mul(r, x), (j * Q[j] % M, 0))
    return r
rho = (0, 1)                           # zeta
for _ in range(12):
    num = Qev(rho); den = Qder(rho); rho = add(rho, mul(((-num[0]) % M, (-num[1]) % M), inv(den)))
print("Q(rho) = 0 mod 2^%d:" % K, Qev(rho) == (0, 0))
C = inv(mul(rho, Qder(rho)))
zeta = (0, 1); w = mul(rho, inv(zeta))
# check f(n) = Tr(C rho^-n) for large n against the companion-matrix value from two_adic_zero_fast
import importlib.util, sys
spec = importlib.util.spec_from_file_location("tz", "two_adic_zero_fast.py")
sys.argv = ["x", "72"]
import io, contextlib
buf = io.StringIO()
with contextlib.redirect_stdout(buf): tz = importlib.util.module_from_spec(spec); spec.loader.exec_module(tz)
M72 = 1 << 72
ok = True
for n in (10**5, 10**5 + 1, 10**5 + 2, 3 * 10**6 + 5):
    val = mul(C, power(inv(rho), n)); tr = (2 * val[0] - val[1]) % M     # Tr(a + b zeta) = 2a - b
    ok &= (tr % M72 == tz.f(n) % M72)
print("f(n) = Tr(C rho^-n) mod 2^72 at large n:", ok)
r = mul(sigma(w), inv(w)); t = mul(((-1) % M, 0), mul(zeta, mul(sigma(C), inv(C))))
print("norms: N(r) =", norm(r) % M == 1, " N(t) =", norm(t) % M == 1)
nstar = 16221783116972227762
rt = power(r, nstar)
diff = ((rt[0] - t[0]) % M, (rt[1] - t[1]) % M)
def v2(x): return (x & -x).bit_length() - 1 if x else K
print("v_2(r^N - t) for N = n* mod 2^64:", min(v2(diff[0]), v2(diff[1])), "(>= 64 + v_2(r - 1) expected if n* = log t / log r)")
print("v_2(r - 1) =", min(v2((r[0] - 1) % M), v2(r[1])))
# quasi-period exponent: f(n + 3*2^j) - f(n) = Tr(C rho^-n (rho^(-3*2^j) - 1)), rho^3 = w^3, so the exponent is j + v_2(w^3 - 1)
w3 = power(w, 3)
print("v_2(w^3 - 1) =", min(v2((w3[0] - 1) % M), v2(w3[1])), " (predicts f(n + 3*2^j) = f(n) mod 2^(j + that))")
print("v_2(C) =", min(v2(C[0]), v2(C[1])))

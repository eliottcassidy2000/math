"""aud_mrb.py -- independent audit of the MRB numerics in the note's section 8 (fermat_self_duality_20261104_audit_S10).
B = sum_{n>=1} (-1)^n (r_n - 1), r_1 = 1, r_2 = 4, r_n = exp(-W0(-ln n / n)) (n >= 3).
Own methods (not mpmath.nsum):
  * r(x) by a bracketed Newton solve of ln r / r = ln x / x on (1, e), cross-checked against lambertw;
  * B by Euler-Boole summation: sum_{j>=0} (-1)^j f(N + j) ~ (1/2) sum_k E_k(0) f^(k)(N) / k!, f = r - 1, with the
    first N - 1 terms summed directly (two values of N);
  * B by the Cohen-Villegas-Zagier alternating-series algorithm (own implementation);
  * length(E) = sum_{k>=1} (r_(2k+1) - r_(2k+2)) by Euler-Maclaurin with an exact tail integral
    (1/2) int_{2K+1}^{2K+2} (r(t) - 1) dt, independently of the alternating series;
  * the MRB constant C by the same Euler-Boole code;
  * r_n = exp(T(ln n / n)) with the Cayley tree series T(z) = sum m^(m-1) z^m / m!.
Then mpmath's default nsum on length(E) and on B, to test the note's caution."""
import mpmath as mp

mp.mp.dps = 50

def r_newton(x):
    x = mp.mpf(x)
    c = mp.log(x) / x                      # target value of ln r / r, r in (1, e)
    lo, hi = mp.mpf(1), mp.e               # ln r / r is increasing on (1, e): bisection, then Newton polish
    for _ in range(60):
        mid = (lo + hi) / 2
        if mp.log(mid) / mid - c > 0: hi = mid
        else: lo = mid
    r = (lo + hi) / 2
    for _ in range(12):
        r = r - (mp.log(r) / r - c) / ((1 - mp.log(r)) / r ** 2)
    return r

def r_w(x):
    x = mp.mpf(x)
    return mp.exp(-mp.lambertw(-mp.log(x) / x, 0).real)

def rr(n):
    if n == 1: return mp.mpf(1)
    if n == 2: return mp.mpf(4)
    return r_w(n)

print("== mirror base")
for n in (3, 4, 5, 10, 1000):
    a, b = r_newton(n), r_w(n)
    print(f"  r_{n}: Newton {mp.nstr(a, 30)}  lambertw {mp.nstr(b, 30)}  |diff| = {mp.nstr(abs(a - b), 3)}  "
          f"r^(1/r) - n^(1/n) = {mp.nstr(a ** (1 / a) - mp.mpf(n) ** (mp.mpf(1) / n), 3)}")

def boole_tail(f, N, K):
    """sum_{j>=0} (-1)^j f(N + j) by Euler-Boole summation with K derivative terms"""
    ders = mp.taylor(f, N, K)          # ders[k] = f^(k)(N) / k!
    return mp.mpf(1) / 2 * sum(mp.eulerpoly(k, 0) * ders[k] for k in range(K + 1))

def B_boole(N, K):
    s = sum((-1) ** n * (rr(n) - 1) for n in range(1, N))
    return s + (-1) ** N * boole_tail(lambda x: r_w(x) - 1, N, K)

def cvz(a, nterms):
    """Cohen-Villegas-Zagier: sum_{k>=0} (-1)^k a(k)"""
    d = (3 + mp.sqrt(8)) ** nterms
    d = (d + 1 / d) / 2
    b, c, s = mp.mpf(-1), -d, mp.mpf(0)
    for k in range(nterms):
        c = b - c
        s += c * a(k)
        b = (k + nterms) * (k - nterms) * b / ((k + mp.mpf(1) / 2) * (k + 1))
    return s / d

print("== B")
with mp.workdps(60):
    B1 = B_boole(200, 24)
    B2 = B_boole(400, 24)
print(f"  Euler-Boole N=200, K=24: {mp.nstr(B1, 32)}")
print(f"  Euler-Boole N=400, K=24: {mp.nstr(B2, 32)}   |difference| = {mp.nstr(abs(B1 - B2), 3)}")
for nt in (40, 80, 120):
    # B = (r_2 - 1) - (r_1 - 1) + sum_{n>=3} (-1)^n (r_n - 1) = 3 - sum_{j>=0} (-1)^j (r_(3+j) - 1)
    Bc = 3 - cvz(lambda j: r_w(3 + j) - 1, nt)
    print(f"  CVZ ({nt} terms):           {mp.nstr(Bc, 32)}   |CVZ - Boole| = {mp.nstr(abs(Bc - B2), 3)}")
PDF = mp.mpf("2.096791124391526125725030675")
print(f"  PDF value 2.096791124391526125725030675: |Boole - PDF| = {mp.nstr(abs(B2 - PDF), 3)}")
print(f"  note value 2.0967911243915261257 (20 digits): |Boole - note| = {mp.nstr(abs(B2 - mp.mpf('2.0967911243915261257')), 3)}")

print("== length(E), independently of the alternating series")
def lenE(K, M):
    g = lambda k: r_w(2 * k + 1) - r_w(2 * k + 2)
    head = sum(g(k) for k in range(1, K))
    integral = mp.quad(lambda t: r_w(t) - 1, [2 * K + 1, 2 * K + 2]) / 2
    ders = mp.taylor(g, K, 2 * M)
    em = g(K) / 2 - sum(mp.bernoulli(2 * j) / mp.factorial(2 * j) * mp.factorial(2 * j - 1) * ders[2 * j - 1] for j in range(1, M + 1))
    return head + integral + em
with mp.workdps(60):
    L1 = lenE(150, 12)
    L2 = lenE(300, 12)
print(f"  Euler-Maclaurin K=150: length(E) = {mp.nstr(L1, 30)};  K=300: {mp.nstr(L2, 30)}")
print(f"  3 - length(E) = {mp.nstr(3 - L2, 30)};  |3 - length(E) - B| = {mp.nstr(abs(3 - L2 - B2), 3)}")

print("== MRB constant C = sum (-1)^n (n^(1/n) - 1)")
with mp.workdps(60):
    f = lambda x: x ** (1 / x) - 1
    N = 300
    C = sum((-1) ** n * (mp.mpf(n) ** (mp.mpf(1) / n) - 1) for n in range(1, N)) + (-1) ** N * boole_tail(f, N, 24)
KNOWN = mp.mpf("0.18785964246206712024851793405427323005590309490013878617200468408947723156466021370329665443310749690384")
print(f"  C (Euler-Boole) = {mp.nstr(C, 40)}")
print(f"  |C - OEIS A037077 digits| = {mp.nstr(abs(C - KNOWN), 3)}")
print(f"  C = sqrt2 - 1 - sum_k [f(r_(2k+1)) - f(r_(2k+2))] is the FTC form of the PDF's integral over E: check at 30 digits:")
with mp.workdps(60):
    # sum_k [ (2k+1)^(1/(2k+1)) - (2k+2)^(1/(2k+2)) ]  =  sqrt2 - 1 - C  (grouped MRB series), via the mirror values
    chk = mp.sqrt(2) - 1 - C
    s_mirror = sum(r_w(2 * k + 1) ** (1 / r_w(2 * k + 1)) - r_w(2 * k + 2) ** (1 / r_w(2 * k + 2)) for k in range(1, 6))
    s_direct = sum(mp.mpf(2 * k + 1) ** (mp.mpf(1) / (2 * k + 1)) - mp.mpf(2 * k + 2) ** (mp.mpf(1) / (2 * k + 2)) for k in range(1, 6))
print(f"  first five interval terms via f(r_n) and via f(n) agree: |diff| = {mp.nstr(abs(s_mirror - s_direct), 3)}")

print("== tree function: r_n = exp(T(ln n / n)), T(z) = sum m^(m-1) z^m / m!")
for n in (3, 5, 10):
    z = mp.log(n) / n
    s, term_log = mp.mpf(0), None
    M = 20000 if n == 3 else 4000
    for m in range(1, M + 1):
        s += mp.exp((m - 1) * mp.log(m) + m * mp.log(z) - mp.loggamma(m + 1))
    # tail bound: terms ~ (e z)^m / (sqrt(2 pi) m^1.5)
    ez = mp.e * z
    tail = ez ** (M + 1) / (mp.sqrt(2 * mp.pi) * (M + 1) ** 1.5 * (1 - ez))
    print(f"  n={n}: exp(series, {M} terms) = {mp.nstr(mp.exp(s), 25)}; r_n = {mp.nstr(r_w(n), 25)}; "
          f"|diff| = {mp.nstr(abs(mp.exp(s) - r_w(n)), 3)} (tail bound {mp.nstr(tail, 3)})")

print("== the note's caution about mpmath's default nsum")
mp.mp.dps = 30
dflt = mp.nsum(lambda k: r_w(2 * k + 1) - r_w(2 * k + 2), [1, mp.inf])
print(f"  default nsum length(E) = {mp.nstr(dflt, 20)}  -> 3 - length(E) = {mp.nstr(3 - dflt, 20)}  (true B = {mp.nstr(B2, 20)})")
dfltB = mp.nsum(lambda n: (-1) ** int(n) * (rr(int(n)) - 1), [1, mp.inf])
print(f"  default nsum on the alternating series B = {mp.nstr(dfltB, 20)}")
for dps in (15, 30, 40):
    mp.mp.dps = dps
    L = mp.nsum(lambda k: r_w(2 * k + 1) - r_w(2 * k + 2), [1, mp.inf])
    print(f"  default nsum at dps={dps}: 3 - length(E) = {mp.nstr(3 - L, 12)}")

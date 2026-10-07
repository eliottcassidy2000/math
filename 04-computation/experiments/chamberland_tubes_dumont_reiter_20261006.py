"""Tubes, shadows and limit cycles of two real extensions of the Collatz map (session opus-2026-10-06-S19).
Note: 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md, section 2.
  C(x) = x + 1/4 - (2x+1)/4 cos(pi x)          Chamberland 1996 (C(n) = T(n) on Z)
  D(x) = (3^s x + s)/2,  s = sin^2(pi x/2)      Dumont-Reiter 2003, the 3-power extension
PRIOR ART for C: Lygeros-Rozier, Ratio Math. 26 (2014) 77-94, arXiv:1402.1979 -- Lemma 2.4 (tube J_n = [n, n + a/(pi^2 n)],
a = 7/2, all n > 0; small odd n checked in floating point), Theorem 3.3 / Corollary 3.4 (odd c_n -> A1, the 3x+1 reformulation),
(5.2)-(5.3) (mirror law, fixed points).  For C this script is a fully interval-certified RE-PROOF; the D part is new.
Rigorous parts use mpmath interval arithmetic (iv); numerical parts are labelled NUMERICAL.
  A. C: displacement even about -1/2; mirror multipliers sum to 2; the attracting fixed points are exactly 0 and -1.27773...
     (interval root isolation: each candidate cell cluster has an interval derivative excluding 0; interval multipliers).
  B. C: negative Schwarzian on [0, oo) (2 C'C''' - 3 C''^2 = pi^2 Q(w), w = pi(x + 1/2)).
  C. Tube theorem: C maps tau(m) = [m, m + a/(2m+1)] into tau(T(m)) for every m >= 1 (a = 0.8, also 0.9);
     D likewise (a = 0.6, also 0.7); the odd critical point lies in tau(n); D'' < 0 on odd tubes (a = 0.6 and 0.7).
  D. Dumont-Reiter: D^2 contracts tau(1) to the cycle (1,2); hence their Odd Critical Point Conjecture (ii), (iii) (real sense)
     hold for every odd n, and (i) holds iff the Collatz orbit of n reaches 1.
  E. Chamberland: in tau(1) the scaled return map W has a fixed point in (0.0705, 0.3) (IVT, rigorous; value 0.0710583 and
     repulsion W' = 1.2148 NUMERICAL) and A2 (0.5775957, NUMERICAL); every odd n >= 7 with a Collatz orbit reaching 1 sends c_n to A1 (tails through 13, 21, 40, 64);
     c_1, c_3, c_5 -> A2 (interval enclosures).
  F. Census (NUMERICAL): odd shadow deviations; the even critical point c_54 leaves the shadow; c_382, c_496, c_502 -> A2.
  G. Flip lemma (rigorous): a point at scaled position xi in [-1.3, -0.45] left of an even integer k >= 2 is mapped into the
     right tube tau(T(k)), hence shadows the Collatz orbit of T(k) forever.
Runtime: about 30 seconds.  Prints ALL CHECKS PASSED."""
import sys, time
from mpmath import mp, iv, mpf, cos, sin, pi, findroot, sqrt, log

iv.dps = 25
mp.dps = 40
PI = iv.pi
LN3 = iv.log(3)
FAIL = []


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        FAIL.append(msg)


# ------------------------------------------------------------------ helpers (rigorous enclosures)
def sinc_iv(u):  # sinc decreasing on [0, pi]; u is an interval inside [0, 1]
    lo, hi = u.a, u.b
    assert lo >= 0 and hi <= 1.0
    s_hi = mpf(1) if lo == 0 else (iv.sin(iv.mpf(lo)) / iv.mpf(lo)).b
    s_lo = (iv.sin(iv.mpf(hi)) / iv.mpf(hi)).a if hi > 0 else mpf(1)
    return iv.mpf([s_lo, s_hi])


def phi_iv(y):  # (1 - e^-y)/y, decreasing for y >= 0
    lo, hi = y.a, y.b
    up = mpf(1) if lo == 0 else ((1 - iv.exp(-iv.mpf(lo))) / iv.mpf(lo)).b
    dn = ((1 - iv.exp(-iv.mpf(hi))) / iv.mpf(hi)).a if hi > 0 else mpf(1)
    return iv.mpf([dn, up])


def psi_iv(y):  # (e^y - 1)/y, increasing for y >= 0
    lo, hi = y.a, y.b
    dn = mpf(1) if lo == 0 else ((iv.exp(iv.mpf(lo)) - 1) / iv.mpf(lo)).a
    up = ((iv.exp(iv.mpf(hi)) - 1) / iv.mpf(hi)).b if hi > 0 else mpf(1)
    return iv.mpf([dn, up])


def parts(xi, h):
    u = PI * xi * h / 2
    return sinc_iv(u), iv.sin(u) ** 2


# exact scaled maps, h = 1/(2m+1), xi = (2m+1)(x - m); image coordinate (2T(m)+1)(F(x) - T(m))
def C_odd_over_xi(xi, h):
    S, _ = parts(xi, h)
    return (3 + h) / 2 * (1 + iv.cos(PI * xi * h) / 2 - PI ** 2 * xi / 8 * S ** 2)


def C_odd(xi, h):
    return xi * C_odd_over_xi(xi, h)


def C_even_over_xi(xi, h):
    S, _ = parts(xi, h)
    return (1 + h) / 2 * ((1 - iv.cos(PI * xi * h) / 2) + PI ** 2 * xi / 8 * S ** 2)


def C_even(xi, h):
    return xi * C_even_over_xi(xi, h)


def C_fprime_odd(xi, h):  # C'(m + xi h), m odd
    w = PI * xi * h
    return 1 + iv.cos(w) / 2 - PI ** 2 * xi / 4 * (1 + 2 * xi * h ** 2) * sinc_iv(w)


def D_odd_over_xi(xi, h):
    S, sig = parts(xi, h)
    return (3 + h) / 4 * (-(3 * (1 - h) / 2) * LN3 * (PI ** 2 * xi / 4) * S ** 2 * phi_iv(sig * LN3)
                          + iv.exp((1 - sig) * LN3) - h * (PI ** 2 * xi / 4) * S ** 2)


def D_odd(xi, h):
    return xi * D_odd_over_xi(xi, h)


def D_even_over_xi(xi, h):
    S, sig = parts(xi, h)
    return (1 + h) / 4 * (((1 - h) / 2) * LN3 * (PI ** 2 * xi / 4) * S ** 2 * psi_iv(sig * LN3)
                          + iv.exp(sig * LN3) + h * (PI ** 2 * xi / 4) * S ** 2)


def D_even(xi, h):
    return xi * D_even_over_xi(xi, h)


def D_prime_odd2(xi, h):  # 2 D'(m + xi h), m odd
    S, sig = parts(xi, h)
    w = PI * xi * h
    sw = sinc_iv(w)
    return (iv.exp((1 - sig) * LN3) * (1 - (PI ** 2 * LN3 / 2) * xi * ((1 - h) / 2 + xi * h ** 2) * sw)
            - (PI / 2) * PI * xi * h * sw)


def hD2_odd(xi, h):  # h D''(m + xi h), m odd
    w = PI * xi * h
    sig = iv.sin(w / 2) ** 2
    xh = (1 - h) / 2 + xi * h ** 2
    sp = -(PI / 2) * iv.sin(w)
    spp = -(PI ** 2 / 2) * iv.cos(w)
    return (iv.exp((1 - sig) * LN3) * LN3 * (LN3 * sp ** 2 * xh + spp * xh + 2 * sp * h) + spp * h) / 2


def cover(pred, x0, x1, h0, h1, nx, nh, maxdepth=10):
    """True iff pred holds on every box of a subdivision of [x0,x1] x [h0,h1] (adaptive)."""
    stack = []
    dx, dh = (x1 - x0) / nx, (h1 - h0) / nh
    for i in range(nx):
        for j in range(nh):
            stack.append((x0 + i * dx, x0 + (i + 1) * dx, h0 + j * dh, h0 + (j + 1) * dh, 0))
    nbox = 0
    while stack:
        a, b, c, d, dep = stack.pop()
        nbox += 1
        if pred(iv.mpf([a, b]), iv.mpf([c, d])):
            continue
        if dep >= maxdepth:
            return False, nbox
        m, n = (a + b) / 2, (c + d) / 2
        stack += [(a, m, c, n, dep + 1), (m, b, c, n, dep + 1), (a, m, n, d, dep + 1), (m, b, n, d, dep + 1)]
    return True, nbox


def Tint(m):
    return m // 2 if m % 2 == 0 else (3 * m + 1) // 2


Cpt = lambda x: x + mpf(1) / 4 - (2 * x + 1) / 4 * cos(pi * x)
Cpp = lambda x: 1 - cos(pi * x) / 2 + (2 * x + 1) * pi / 4 * sin(pi * x)


def Dpt(x):
    s = sin(pi * x / 2) ** 2
    return (mpf(3) ** s * x + s) / 2


def Dpp(x):
    s = sin(pi * x / 2) ** 2
    return (mpf(3) ** s * (1 + x * log(3) * pi / 2 * sin(pi * x)) + pi / 2 * sin(pi * x)) / 2


def Cp_iv_global(X):  # interval C'(x)
    return 1 - iv.cos(PI * X) / 2 + (2 * X + 1) * PI / 4 * iv.sin(PI * X)


t0 = time.time()
# ================================================================== A
print('A. Chamberland: mirror law and attracting fixed points')
g = lambda x: Cpt(x) - x
import random
random.seed(19)
err = max(abs(g(mpf(-0.5) + u) - g(mpf(-0.5) - u)) for u in [mpf(random.uniform(-40, 40)) for _ in range(200)])
check(err < mpf(10) ** -30, f'displacement g(x) = C(x) - x is even about -1/2: max |g(-1/2+u) - g(-1/2-u)| = {float(err):.1e} (identity: g(-1/2+u) = 1/4 - (u/2) sin(pi u))')
fps = []
for k in range(-40, 41):
    for x0 in (k + mpf('0.25'), k + mpf('0.75'), mpf(k)):
        try:
            r = findroot(lambda x: (2 * x + 1) * cos(pi * x) - 1, x0)
        except Exception:
            continue
        if abs(r) < 20 and all(abs(r - q) > mpf(10) ** -15 for q in fps):
            fps.append(r)
fps.sort()
pairs_ok = all(abs(Cpp(r) + Cpp(-1 - r) - 2) < mpf(10) ** -25 for r in fps)
check(pairs_ok, f'mirror multipliers sum to 2 at all {len(fps)} fixed points in [-20, 20]  (C\'(-1-x) = 2 - C\'(x))')
# rigorous classification: |C'| < 1 at a fixed point forces |2x+1| < sqrt(1 + 100/pi^2) = 3.3359
wmax = sqrt(1 + 100 / pi ** 2)
print(f'     at a fixed point C\' = 1 - 1/(2w) +- (pi/4) sqrt(w^2-1), w = 2x+1; |C\'| < 1 forces |w| < sqrt(1 + 100/pi^2) = {float(wmax):.5f}, x in ({float((-wmax-1)/2):.4f}, {float((wmax-1)/2):.4f})')
# all zeros of g(x) = (2x+1)cos(pi x) - 1 in that window, isolated rigorously: interval grid, then on each cluster of
# candidate cells the interval derivative g'(x) = 2 cos(pi x) - pi (2x+1) sin(pi x) excludes 0 and g changes sign at the ends
lo, hi = (-wmax - 1) / 2, (wmax - 1) / 2
N = 4000
cells = []
for k in range(N):
    a = lo + (hi - lo) * k / N
    b = lo + (hi - lo) * (k + 1) / N
    X = iv.mpf([a, b])
    val = (2 * X + 1) * iv.cos(PI * X) - 1
    if val.a > 0 or val.b < 0:
        continue
    cells.append((a, b))
clusters = []
for a, b in cells:
    if clusters and abs(clusters[-1][1] - a) < mpf(10) ** -30:
        clusters[-1] = (clusters[-1][0], b)
    else:
        clusters.append((a, b))
gfun = lambda X: (2 * X + 1) * iv.cos(PI * X) - 1
iso_ok = True
mults = []
for a, b in clusters:
    X = iv.mpf([a, b])
    dg = 2 * iv.cos(PI * X) - PI * (2 * X + 1) * iv.sin(PI * X)
    ga, gb = gfun(iv.mpf(a)), gfun(iv.mpf(b))
    if not ((dg.a > 0 or dg.b < 0) and ((ga.b < 0 and gb.a > 0) or (ga.a > 0 and gb.b < 0))):
        iso_ok = False
    mults.append((float(a), float(b), Cp_iv_global(X)))
check(iso_ok and len(clusters) == 4, f'fixed points of C with |2x+1| < 3.33645: exactly {len(clusters)} (each candidate cluster has g\' excluding 0 and a sign change)')
attr_iv = [(a, b, M) for a, b, M in mults if M.a > -1 and M.b < 1]
rep_iv = [(a, b, M) for a, b, M in mults if M.a > 1 or M.b < -1]
check(len(attr_iv) == 2 and len(rep_iv) == 2 and attr_iv[0][0] < -1.2777337 < attr_iv[0][1] and attr_iv[1][0] <= 0 <= attr_iv[1][1],
      'interval multipliers: attracting exactly at 0 (C\' in [' + f'{float(attr_iv[1][2].a):.4f}, {float(attr_iv[1][2].b):.4f}]) and -1.27773 (C\' in [{float(attr_iv[0][2].a):.4f}, {float(attr_iv[0][2].b):.4f}]); '
      f'0.27773 and -1 repelling')
attr = [r for r in fps if abs(Cpp(r)) < 1]
print(f'     attracting fixed points (40 digits): {[(round(float(r), 10), round(float(Cpp(r)), 9)) for r in attr]}; 0.2777338 has C\' = 1.614292 (repelling)')

# ================================================================== B
print('B. Chamberland: negative Schwarzian on [0, oo)')


def Qw(w):
    s, c = iv.sin(w), iv.cos(w)
    return -(iv.mpf(1) / 2 + s * s / 4) * w * w + c * (1 + s) * w + iv.mpf(3) / 2 * (s * s + 2 * s - 2)


w0 = 2 * sqrt(3)
okB, nb = cover(lambda X, H: Qw(X).b < 0, mp.pi / 2, w0 + mpf('0.01'), 0, 1, 4000, 1)
check(okB, f'Q(w) < 0 on [pi/2, 2 sqrt3] (x in [0, 0.6027]) by {nb} interval cells; for w > 2 sqrt3, Q <= -w^2/2 + (3 sqrt3/4) w + 3/2 < 0')
Qpt = lambda w: -(mpf(1) / 2 + sin(w) ** 2 / 4) * w * w + cos(w) * (1 + sin(w)) * w + mpf(3) / 2 * (sin(w) ** 2 + 2 * sin(w) - 2)
wstar = findroot(Qpt, mpf('1.50'))
check(Qpt(pi / 2) < 0 and abs(wstar / pi - mpf(1) / 2 + mpf('0.0217160')) < 1e-6,
      f'Q(pi/2) = {float(Qpt(pi / 2)):.4f} (x = 0); the sign change of Q nearest 0 is at x = {float(wstar / pi - mpf(1) / 2):.7f} (Q > 0 just left of it)')
from mpmath import diff
f1 = lambda x: diff(Cpt, x, 1); f2 = lambda x: diff(Cpt, x, 2); f3 = lambda x: diff(Cpt, x, 3)
idt = max(abs(2 * f1(x) * f3(x) - 3 * f2(x) ** 2 - pi ** 2 * (lambda w: -(mpf(1)/2 + sin(w)**2/4)*w*w + cos(w)*(1+sin(w))*w + mpf(3)/2*(sin(w)**2 + 2*sin(w) - 2))(pi * (x + mpf(1) / 2))) for x in [mpf('0.13'), mpf('1.7'), mpf('9.31')])
check(idt < mpf(10) ** -20, f'identity 2C\'C\'\'\' - 3C\'\'^2 = pi^2 Q(pi(x+1/2)) at sample points (max err {float(idt):.1e})')

# ================================================================== C
print('C. Tube theorem (both maps, every integer m >= 1)')
for m in (1, 2, 3, 4, 9, 10, 101, 1000):  # exactness of the scaled maps
    for xi in (mpf('0.11'), mpf('0.43'), mpf('0.67')):
        h = mpf(1) / (2 * m + 1)
        for F, Fo, Fe in ((Cpt, C_odd, C_even), (Dpt, D_odd, D_even)):
            direct = (2 * Tint(m) + 1) * (F(m + xi * h) - Tint(m))
            enc = (Fo if m % 2 else Fe)(iv.mpf(xi), iv.mpf(h))
            if not (enc.a - mpf(10) ** -18 <= direct <= enc.b + mpf(10) ** -18):
                FAIL.append(f'scaled map mismatch {F} m={m} xi={xi}')
check(not [f for f in FAIL if 'mismatch' in f], 'scaled maps O, E agree with direct evaluation (48 cases)')
# tube constants are exact decimals: boxes cover xi in [0, A.b] and images are compared with A.a (A = interval enclosing a);
# the h-range [0, 1/3] (odd m >= 1) and [0, 1/5] (even m >= 2) is covered up to mpf(1)/3, mpf(1)/5 >= the exact endpoints
H3 = (iv.mpf(1) / 3).b
H5 = (iv.mpf(1) / 5).b
for s_a in ('0.8', '0.9'):
    A = iv.mpf(s_a)
    r1 = cover(lambda X, H: C_odd_over_xi(X, H).a > 0, 0, A.b, 0, H3, 40, 10)
    r2 = cover(lambda X, H: C_odd(X, H).b <= A.a, 0, A.b, 0, H3, 40, 10)
    r3 = cover(lambda X, H: C_even(X, H).b <= A.a, 0, A.b, 0, H5, 40, 10)
    check(r1[0] and r2[0] and r3[0], f'C, a = {s_a}: odd m >= 1: 0 < O(xi) <= a on (0,a]; even m >= 2: 0 <= E(xi) <= a  ({r1[1]+r2[1]+r3[1]} boxes)')
    rc = cover(lambda X, H: C_fprime_odd(A, H).b < 0, 0, 1, 0, H3, 1, 200)
    check(rc[0], f"C'(m) = 3/2 > 0 and C'(m + a/(2m+1)) < 0 for all odd m (a = {s_a}); C'' < 0 on [m, m+1/2) (analytic): unique critical point c_m in tau(m)")
for s_a in ('0.6', '0.7'):
    A = iv.mpf(s_a)
    r1 = cover(lambda X, H: D_odd_over_xi(X, H).a > 0, 0, A.b, 0, H3, 40, 10)
    r2 = cover(lambda X, H: D_odd(X, H).b <= A.a, 0, A.b, 0, H3, 40, 10)
    r3 = cover(lambda X, H: D_even(X, H).b <= A.a, 0, A.b, 0, H5, 40, 10)
    check(r1[0] and r2[0] and r3[0], f'D, a = {s_a}: same three inequalities ({r1[1]+r2[1]+r3[1]} boxes)')
    rd1 = cover(lambda X, H: D_prime_odd2(A, H).b < 0, 0, 1, 0, H3, 1, 200)
    rd2 = cover(lambda X, H: hD2_odd(X, H).b < 0, 0, A.b, 0, H3, 40, 10)
    check(rd1[0] and rd2[0], f"D (a = {s_a}): D'(m) = 3/2 > 0, D'(m + a/(2m+1)) < 0 and D'' < 0 on tau(m) for all odd m: unique critical point c_m in tau(m)")
aC, aD = mpf('0.8'), mpf('0.6')


def max_enclosure(fun, a_b, hmax, nx=80, nh=20):
    m = mpf(-10)
    for i in range(nx):
        for j in range(nh):
            X = iv.mpf([a_b * i / nx, a_b * (i + 1) / nx])
            H = iv.mpf([hmax * j / nh, hmax * (j + 1) / nh])
            m = max(m, fun(X, H).b)
    return m


print('     upper bounds of the images (interval enclosures): ' + ', '.join(
    f'{nm} {float(max_enclosure(fn, iv.mpf(sa).b, hm)):.4f}' for nm, fn, sa, hm in
    (('C O a=0.8', C_odd, '0.8', H3), ('C E a=0.8', C_even, '0.8', H5), ('C E a=0.9', C_even, '0.9', H5),
     ('D O a=0.6', D_odd, '0.6', H3), ('D E a=0.6', D_even, '0.6', H5), ('D E a=0.7', D_even, '0.7', H5))))
E095 = C_even(iv.mpf('0.95'), iv.mpf(1) / 5)
check(E095.a > iv.mpf('0.95').b, f'the constant a = 0.95 fails for C at m = 2: E(0.95, 1/5) = [{float(E095.a):.4f}, {float(E095.b):.4f}] > 0.95')

# ================================================================== D
print('D. Dumont-Reiter: tau(1) is contracted to (1,2); the Odd Critical Point Conjecture')


def W_D_ratio(X):
    h1, h2 = iv.mpf(1) / 3, iv.mpf(1) / 5
    Oxi = D_odd_over_xi(X, h1)
    return D_even_over_xi(X * Oxi, h2) * Oxi


rW = cover(lambda X, H: W_D_ratio(X).b < 1, 0, iv.mpf('0.6').b, 0, 1, 300, 1)
check(rW[0], 'W_D(xi)/xi < 1 on (0, 0.6] (W_D = E(., 1/5) o O(., 1/3) is D^2 on tau(1) in scaled form): every point of tau(1) -> (1,2)')
mu = sorted(set(round(float(findroot(lambda x: Dpt(x) - x, mpf(x0))), 12) for x0 in ('0.3', '1.5', '2.5')))
check(abs(mu[0] - 0.3158162033) < 1e-8 and abs(mu[1] - 1.5155526112) < 1e-8,
      f'D fixed points mu1, mu2 = {mu[0]:.10f}, {mu[1]:.10f}: tau(1) = [1, 1.2] lies in (mu1, mu2), tau(m) lies in [2, oo) for m >= 2')
print('     => (ii) total stopping time of c_n equals that of n (finite or infinite) for EVERY odd n >= 1;')
print('        (iii) (real sense) tau(n) is a connected set containing n and c_n on which the total stopping time is constant;')
print('        (i) c_n -> (1,2) iff the Collatz orbit of n reaches 1 (otherwise c_n stays in tubes of integers >= 3).')

# ================================================================== E
print('E. Chamberland: the return map on tau(1); A1, the separatrix and A2; where odd critical orbits go')
Wpt = lambda xi: 3 * (Cpt(Cpt(1 + xi / 3)) - 1)
xiR = findroot(lambda t: Wpt(t) - t, mpf('0.07'))
xiA2 = findroot(lambda t: Wpt(t) - t, mpf('0.578'))
xic = 3 * (findroot(Cpp, mpf('1.18')) - 1)
dW = lambda t: diff(Wpt, t)
print(f'     NUMERICAL (40 digits): repelling fixed point of W at xi_R = {float(xiR):.10f} (W\' = {float(dW(xiR)):.6f}); '
      f'A2 at xi = {float(xiA2):.10f} (W\' = {float(dW(xiA2)):.6f}); critical point of W at xi_c1 = {float(xic):.7f}')
check(abs(xiR - mpf('0.0710583592786')) < 1e-10 and abs(xiA2 - 3 * (mpf('1.1925319070466') - 1)) < 1e-10,
      'the two nonzero fixed points of W located (NUMERICAL); the rigorous statements follow')


def W_C(X):
    return C_even(C_odd(X, iv.mpf(1) / 3), iv.mpf(1) / 5)


def W_C_ratio(X):
    Oxi = C_odd_over_xi(X, iv.mpf(1) / 3)
    return C_even_over_xi(X * Oxi, iv.mpf(1) / 5) * Oxi


cut = mpf('0.0705')
rA1 = cover(lambda X, H: W_C_ratio(X).b < 1, 0, iv.mpf('0.0705').b, 0, 1, 400, 1)
W03 = W_C(iv.mpf('0.3'))
check(rA1[0] and W03.a > mpf('0.3'), f'W_C(xi) < xi on (0, {cut}] (and W_C >= 0): [0, {cut}] lies in the basin of A1; W_C(0.3) > 0.3, '
      f'so a fixed point (the separatrix) lies in ({cut}, 0.3) by the IVT')
# backward tree of 1: every odd n >= 7 whose orbit reaches 1 passes through 13, 21, 40 or 64
tails = {13: [13, 20, 10, 5, 8, 4, 2, 1], 40: [40, 20, 10, 5, 8, 4, 2, 1], 21: [21, 32, 16, 8, 4, 2, 1], 64: [64, 32, 16, 8, 4, 2, 1]}


def push(X, path):
    for m in path[:-1]:
        h = iv.mpf(1) / (2 * m + 1)
        X = C_odd(X, h) if m % 2 else C_even(X, h)
        X = iv.mpf([max(mpf(0), X.a), X.b])
    return X


worst = {}
for node, path in tails.items():
    sup = mpf(0)
    K = 400
    for k in range(K):
        top = iv.mpf('0.8').b
        X = iv.mpf([top * k / K, top * (k + 1) / K])
        sup = max(sup, push(X, path).b)
    worst[node] = sup
check(all(v < iv.mpf('0.0705').a for v in worst.values()), 'tube points at 13, 21, 40, 64 enter tau(1) at xi <= ' + ', '.join(f'{k}: {float(v):.5f}' for k, v in worst.items()) + f' < {cut}')


# backward-tree claim (finite, exact): odd n >= 7 never meets 3, 6, 12, ...; the tree of 1 to depth 7 has exactly these frontier nodes
def preds(m):
    out = [2 * m]
    if (2 * m - 1) % 3 == 0 and ((2 * m - 1) // 3) % 2 == 1:
        out.append((2 * m - 1) // 3)
    return out


tree_ok = True
for n in range(7, 200001, 2):
    x, seen = n, False
    while x != 1:
        if x in (13, 21, 40, 64):
            seen = True
            break
        if x in (5, 8, 16, 32, 10, 20):
            break
        x = Tint(x)
    if not seen:
        tree_ok = False
        print('     tail exception', n)
        break
check(tree_ok and sorted(preds(20)) == [13, 40] and sorted(preds(32)) == [21, 64] and preds(5) == [10, 3],
      'every odd 7 <= n <= 200001 passes through 13, 21, 40 or 64 (and the backward tree: 20 <- 13, 40; 32 <- 21, 64; 5 <- 10, 3; 3 is reached only from 3*2^k)')
print('     => for every odd n >= 7 whose Collatz orbit reaches 1, the critical point c_n of C is attracted to A1 = {1,2}.')
# c1, c3, c5 -> A2 : trap J around A2 for W, with |W'| < 1 on J (rigorous, subdivided)
def C_iv(X):
    return X + iv.mpf(1) / 4 - (2 * X + 1) / 4 * iv.cos(PI * X)


def Cp_iv(X):
    return 1 - iv.cos(PI * X) / 2 + (2 * X + 1) * PI / 4 * iv.sin(PI * X)


Ja, Jb = mpf('0.55'), mpf('0.60')
img_lo, img_hi, dmax = mpf(10), mpf(-10), mpf(0)
K = 400
for k in range(K):
    Xi = iv.mpf([Ja + (Jb - Ja) * k / K, Ja + (Jb - Ja) * (k + 1) / K])
    X = 1 + Xi / 3
    Y = C_iv(X)
    Z = C_iv(Y)
    Wk = 3 * (Z - 1)
    img_lo, img_hi = min(img_lo, Wk.a), max(img_hi, Wk.b)
    d = abs(Cp_iv(Y) * Cp_iv(X))
    dmax = max(dmax, d.b)
check(img_lo >= Ja and img_hi <= Jb and dmax < 1,
      f'W_C maps J = [0.55, 0.60] into [{float(img_lo):.5f}, {float(img_hi):.5f}] with |W\'| <= {float(dmax):.4f} < 1 on J: J lies in the basin of A2')
entries = {}
for n in (1, 3, 5):
    c = findroot(Cpp, mpf(n) + mpf('0.6') / (2 * n + 1))
    eps = mpf(10) ** -12
    ok_enc = Cp_iv(iv.mpf(c - eps)).a > 0 and Cp_iv(iv.mpf(c + eps)).b < 0
    X = iv.mpf([c - eps, c + eps])     # rigorous enclosure of the critical point (C' changes sign, C'' < 0)
    m = n
    while m != 1:
        X, m = C_iv(X), Tint(m)
    Xi = 3 * (X - 1)
    k = 0
    while not (Xi.a >= Ja and Xi.b <= Jb) and k < 60:
        Xi = 3 * (C_iv(C_iv(1 + Xi / 3)) - 1)
        k += 1
    entries[n] = (ok_enc, k, float(Xi.b - Xi.a))
check(all(v[0] and v[1] < 60 for v in entries.values()),
      f'c_1, c_3, c_5 (interval enclosures) enter J after {[(n, v[1]) for n, v in entries.items()]} returns: attracted to A2 (rigorous)')

# ================================================================== F
print('F. Census (NUMERICAL)')


def odd_shadow(F, Fp, n, a):
    c = findroot(Fp, mpf(n) + mpf('0.35') / (2 * n + 1) * (1.7 if F is Cpt else 1))
    x, m, mx = c, n, mpf(0)
    while m != 1:
        xi = (2 * m + 1) * (x - m)
        if xi < -mpf(10) ** -30 or xi > a:
            return None
        mx = max(mx, xi)
        x, m = F(x), Tint(m)
    return mx


mxC = max(odd_shadow(Cpt, Cpp, n, mpf('0.8')) for n in range(1, 402, 2))
mxD = max(odd_shadow(Dpt, Dpp, n, mpf('0.6')) for n in range(1, 402, 2))
check(mxC < 0.8 and mxD < 0.6, f'odd n <= 401: scaled shadow deviation stays in tau along the whole orbit (max {float(mxC):.4f} for C, {float(mxD):.4f} for D)')
mp.dps = 60
c54 = findroot(Cpp, mpf(54) - mpf('0.2') / 109)
x, m, esc = c54, 54, None
for t in range(40):
    xi = (2 * m + 1) * (x - m)
    if esc is None and abs(xi) > mpf('1.2'):
        esc = (t, m, float(x))
    x, m = Cpt(x), Tint(m)
check(esc is not None, f'C: the even critical point c_54 = {float(c54):.6f} leaves the near-integer regime at step {esc[0]} (integer {esc[1]}, real {esc[2]:.4f}): even critical points need not shadow')
# even critical points are not A1-only (Lygeros-Rozier 2014, Chamberland 1996 census): c_382, c_496, c_502 -> A2
mp.dps = 1500
A2pts = (mpf('1.192531907046640255'), mpf('2.138656335516705548'))
fates = {}
for m0 in (382, 496, 502):
    c = findroot(Cpp, mpf(m0) - mpf('0.2') / (2 * m0 + 1))
    x, fate = c, None
    for it in range(20000):
        if abs(x - 1) < mpf(10) ** -12 or abs(x - 2) < mpf(10) ** -12:
            fate = 'A1'; break
        if abs(x - A2pts[0]) < mpf(10) ** -12 or abs(x - A2pts[1]) < mpf(10) ** -12:
            fate = 'A2'; break
        x = Cpt(x)
    fates[m0] = fate
check(all(v == 'A2' for v in fates.values()), f'C: even critical points c_382, c_496, c_502 -> {fates} (NUMERICAL, 1500 digits; as in Lygeros-Rozier 2014)')
mp.dps = 40

# ================================================================== G
print('G. Flip lemma for even integers (rigorous): left deviations in [-1.3, -0.45] are thrown into the right tube')


def sinc_abs(u):  # sinc is even and decreasing in |u| on [0, pi]
    if u.a <= 0 <= u.b:
        lo = mpf(0)
    else:
        lo = min(abs(u.a), abs(u.b))
    hi = max(abs(u.a), abs(u.b))
    s_hi = mpf(1) if lo == 0 else (iv.sin(iv.mpf(lo)) / iv.mpf(lo)).b
    s_lo = (iv.sin(iv.mpf(hi)) / iv.mpf(hi)).a
    return iv.mpf([s_lo, s_hi])


def C_even_signed(xi, h):  # exact scaled even map, valid for xi of either sign
    S = sinc_abs(PI * xi * h / 2)
    return (1 + h) / 2 * (xi * (1 - iv.cos(PI * xi * h) / 2) + PI ** 2 * xi ** 2 / 8 * S ** 2)


A08 = iv.mpf('0.8')
rG = cover(lambda X, H: (lambda v: v.a >= 0 and v.b <= A08.a)(C_even_signed(X, H)), iv.mpf('-1.3').a, iv.mpf('-0.45').b, 0, H5, 200, 20)
check(rG[0], f'for every even k >= 2 and xi in [-1.3, -0.45]: 0 <= E(xi, 1/(2k+1)) <= 0.8 ({rG[1]} boxes): the point k + xi/(2k+1) maps into tau(T(k)) and shadows forever (THM-4563)')
print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')

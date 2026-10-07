"""Tubes, shadows and limit cycles of two real extensions of the Collatz map (session opus-2026-10-06-S19).
Note: 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md, section 2.
  C(x) = x + 1/4 - (2x+1)/4 cos(pi x)          Chamberland 1996 (C(n) = T(n) on Z)
  D(x) = (3^s x + s)/2,  s = sin^2(pi x/2)      Dumont-Reiter 2003, the 3-power extension
Rigorous parts use mpmath interval arithmetic (iv); numerical parts are labelled.
  A. C: displacement even about -1/2; mirror multipliers sum to 2; the attracting fixed points are exactly 0 and -1.27773...
  B. C: negative Schwarzian on [0, oo) (2 C'C''' - 3 C''^2 = pi^2 Q(w), w = pi(x + 1/2)).
  C. Tube theorem: C maps tau(m) = [m, m + a/(2m+1)] into tau(T(m)) for every m >= 1 (a = 0.8, also 0.9);
     D likewise (a = 0.6, also 0.7); the odd critical point lies in tau(n); D'' < 0 on odd tubes.
  D. Dumont-Reiter: D^2 contracts tau(1) to the cycle (1,2); hence their Odd Critical Point Conjecture (ii), (iii) hold for every
     odd n, and (i) holds iff the Collatz orbit of n reaches 1.
  E. Chamberland: in tau(1) the scaled return map W has fixed points 0 (A1), 0.0710583 (repelling) and 0.5775957 (A2);
     every odd n >= 7 with a Collatz orbit reaching 1 sends c_n to A1 (backward-tree tails through 13, 21, 40, 64); c_1, c_3, c_5 -> A2.
  F. Census (NUMERICAL): critical-orbit fates; the even critical point near 54 leaves the shadow.
Runtime: about 1-2 minutes.  Prints ALL CHECKS PASSED."""
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
print(f'     at a fixed point C\' = 1 - 1/(2w) +- (pi/4) sqrt(w^2-1), w = 2x+1; |C\'| < 1 forces |w| < {float(wmax):.4f}, x in ({float((-wmax-1)/2):.4f}, {float((wmax-1)/2):.4f})')
# all zeros of (2x+1)cos(pi x) - 1 in that window, by interval sign checks on a fine grid
lo, hi = (-wmax - 1) / 2, (wmax - 1) / 2
N = 4000
signs = []
zero_cells = 0
for k in range(N):
    a = lo + (hi - lo) * k / N
    b = lo + (hi - lo) * (k + 1) / N
    X = iv.mpf([a, b])
    val = (2 * X + 1) * iv.cos(PI * X) - 1
    if val.a > 0 or val.b < 0:
        continue
    zero_cells += 1
window_fps = [r for r in fps if lo < r < hi]
check(zero_cells <= 2 * len(window_fps) + 2, f'interval grid: {zero_cells} cells may contain zeros; known fixed points in window: {[round(float(r), 6) for r in window_fps]}')
attr = [r for r in fps if abs(Cpp(r)) < 1]
check(len(attr) == 2 and abs(attr[0] - mpf('-1.2777337661')) < 1e-9 and abs(attr[1]) < 1e-20,
      f'attracting fixed points: {[(round(float(r), 9), round(float(Cpp(r)), 9)) for r in attr]}  (0.2777... is REPELLING, C\' = {float(Cpp(fps[[abs(r - mpf("0.2777")) < 0.01 for r in fps].index(True)])):.6f})')

# ================================================================== B
print('B. Chamberland: negative Schwarzian on [0, oo)')


def Qw(w):
    s, c = iv.sin(w), iv.cos(w)
    return -(iv.mpf(1) / 2 + s * s / 4) * w * w + c * (1 + s) * w + iv.mpf(3) / 2 * (s * s + 2 * s - 2)


w0 = 2 * sqrt(3)
okB, nb = cover(lambda X, H: Qw(X).b < 0, mp.pi / 2, w0 + mpf('0.01'), 0, 1, 4000, 1)
check(okB, f'Q(w) < 0 on [pi/2, 2 sqrt3] (x in [0, 0.6027]) by {nb} interval cells; for w > 2 sqrt3, Q <= -w^2/2 + (3 sqrt3/4) w + 3/2 < 0')
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
for aC in (mpf('0.8'), mpf('0.9')):
    r1 = cover(lambda X, H: C_odd_over_xi(X, H).a > 0, 0, aC, 0, mpf(1) / 3, 40, 10)
    r2 = cover(lambda X, H: C_odd(X, H).b <= aC, 0, aC, 0, mpf(1) / 3, 40, 10)
    r3 = cover(lambda X, H: C_even(X, H).b <= aC, 0, aC, 0, mpf(1) / 5, 40, 10)
    check(r1[0] and r2[0] and r3[0], f'C, a = {aC}: odd m >= 1: 0 < O(xi) <= a on (0,a]; even m >= 2: 0 <= E(xi) <= a  ({r1[1]+r2[1]+r3[1]} boxes)')
aC = mpf('0.8')
rc = cover(lambda X, H: C_fprime_odd(iv.mpf(aC), H).b < 0, 0, 1, 0, mpf(1) / 3, 1, 200)
check(rc[0], f"C'(m) = 3/2 > 0 and C'(m + a/(2m+1)) < 0 for all odd m (a = {aC}); C'' < 0 on [m, m+1/2) (analytic): unique critical point c_m in tau(m)")
for aD in (mpf('0.6'), mpf('0.7')):
    r1 = cover(lambda X, H: D_odd_over_xi(X, H).a > 0, 0, aD, 0, mpf(1) / 3, 40, 10)
    r2 = cover(lambda X, H: D_odd(X, H).b <= aD, 0, aD, 0, mpf(1) / 3, 40, 10)
    r3 = cover(lambda X, H: D_even(X, H).b <= aD, 0, aD, 0, mpf(1) / 5, 40, 10)
    check(r1[0] and r2[0] and r3[0], f'D, a = {aD}: same three inequalities ({r1[1]+r2[1]+r3[1]} boxes)')
aD = mpf('0.6')
rd1 = cover(lambda X, H: D_prime_odd2(iv.mpf(aD), H).b < 0, 0, 1, 0, mpf(1) / 3, 1, 200)
rd2 = cover(lambda X, H: hD2_odd(X, H).b < 0, 0, aD, 0, mpf(1) / 3, 40, 10)
check(rd1[0] and rd2[0], "D'(m) = 3/2 > 0, D'(m + a/(2m+1)) < 0 and D'' < 0 on tau(m) for all odd m: unique critical point c_m in tau(m)")

# ================================================================== D
print('D. Dumont-Reiter: tau(1) is contracted to (1,2); the Odd Critical Point Conjecture')


def W_D_ratio(X):
    h1, h2 = iv.mpf(1) / 3, iv.mpf(1) / 5
    Oxi = D_odd_over_xi(X, h1)
    return D_even_over_xi(X * Oxi, h2) * Oxi


rW = cover(lambda X, H: W_D_ratio(X).b < 1, 0, aD, 0, 1, 300, 1)
check(rW[0], 'W_D(xi)/xi < 1 on (0, 0.6] (W_D = E(., 1/5) o O(., 1/3) is D^2 on tau(1) in scaled form): every point of tau(1) -> (1,2)')
mu = sorted(set(round(float(findroot(lambda x: Dpt(x) - x, mpf(x0))), 12) for x0 in ('0.3', '1.5', '2.5')))
check(abs(mu[0] - 0.3158162033) < 1e-8 and abs(mu[1] - 1.5155526112) < 1e-8,
      f'D fixed points mu1, mu2 = {mu[0]:.10f}, {mu[1]:.10f}: tau(1) = [1, 1.2] lies in (mu1, mu2), tau(m) lies in [2, oo) for m >= 2')
print('     => (ii) total stopping time of c_n equals that of n (finite or infinite) for EVERY odd n >= 1;')
print('        (iii) tau(n) is a connected set containing n and c_n on which the total stopping time is constant;')
print('        (i) c_n -> (1,2) iff the Collatz orbit of n reaches 1.')

# ================================================================== E
print('E. Chamberland: the return map on tau(1); A1, the separatrix and A2; where odd critical orbits go')
Wpt = lambda xi: 3 * (Cpt(Cpt(1 + xi / 3)) - 1)
xiR = findroot(lambda t: Wpt(t) - t, mpf('0.07'))
xiA2 = findroot(lambda t: Wpt(t) - t, mpf('0.578'))
xic = 3 * (findroot(Cpp, mpf('1.18')) - 1)
dW = lambda t: diff(Wpt, t)
check(abs(xiR - mpf('0.0710583592786')) < 1e-10 and dW(xiR) > 1, f'repelling fixed point of W at xi_R = {float(xiR):.10f} (W\' = {float(dW(xiR)):.6f})')
check(abs(xiA2 - 3 * (mpf('1.1925319070466') - 1)) < 1e-10 and abs(dW(xiA2)) < 1, f'A2 at xi = {float(xiA2):.10f} (W\' = {float(dW(xiA2)):.6f}); critical point of W at xi_c1 = {float(xic):.7f}')


def W_C(X):
    return C_even(C_odd(X, iv.mpf(1) / 3), iv.mpf(1) / 5)


def W_C_ratio(X):
    Oxi = C_odd_over_xi(X, iv.mpf(1) / 3)
    return C_even_over_xi(X * Oxi, iv.mpf(1) / 5) * Oxi


cut = mpf('0.0705')
rA1 = cover(lambda X, H: W_C_ratio(X).b < 1, 0, cut, 0, 1, 400, 1)
check(rA1[0], f'W_C(xi) < xi on (0, {cut}] and W_C increasing there (xi_c1 > cut): [0, {cut}] lies in the basin of A1')
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
        X = iv.mpf([aC * k / K, aC * (k + 1) / K])
        sup = max(sup, push(X, path).b)
    worst[node] = sup
check(all(v < cut for v in worst.values()), 'tube points at 13, 21, 40, 64 enter tau(1) at xi <= ' + ', '.join(f'{k}: {float(v):.5f}' for k, v in worst.items()) + f' < {cut}')


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
print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')

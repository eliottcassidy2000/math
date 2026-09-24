#!/usr/bin/env python3
"""collatz_procgen_20260924_fixedpt_inequality.py

Lane: fixed points / the owner's inequality (session collatz-procgen-20260922, 2026-09-24).

The owner's sandwich (Lemma 2.1 of arXiv:2502.20642v1):
    |d(x,z) - d(z,y)| <= d(x,y) <= d(x,z) + d(z,y),
    d(x,z)^2 - 2 d(x,z) d(z,y) + d(z,y)^2 <= d(x,y)^2 <= d(x,z)^2 + 2 d(x,z) d(z,y) + d(z,y)^2.

Sections (printed):
  I1  Lemma 2.1 is the triangle inequality followed by AM-QM; what it discards (the cross term's sign)
  I2  on the line every consecutive orbit triple is an EQUALITY case; the cross-term sign is parity
      agreement, the same on both half-lines (side-blind); the strict middle needs a second dimension
  I3  the sign law is the pair of equality cases, with z = 0 as the middle point; transport x -> -x
  I4  the central sheet b = 0: midpoint and odd part of the two corners; odd branch = central map
      centred at -b; the shadowing identity T^n(x) = (3^a_n/2^n) xi_n; the THM-4469 instance
  I5  orthogonality in C: Pythagorean locus of affine branches; Chamberland's f vs the
      Letherman-Schleicher-Wood f_0 (integers critical: the fold turns orthogonal into opposite)
  I6  |T(x)| = T_sgn(x)(|x|): the correction is sgn x = d|x|/dx; pair sums (linear) conserved, pair
      products (quadratic) shrink, the cross term flows into the squares
  I7  the fixed points {0,-b} are the centres of the two branches; the self-mapped consecutive pairs

Every check raises on failure.  Runtime about 1 min; peak memory under 300 MB; one process.
"""
import sys
import time
import random
import resource
from fractions import Fraction as Fr
from math import sqrt

import mpmath as mp
import numpy as np
import sympy as sp

T0 = time.time()
random.seed(20260924)


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


def sgn(x):
    return (x > 0) - (x < 0)


# ----------------------------------------------------------------------------------------------
hdr("I1  Lemma 2.1 = triangle inequality + AM-QM; what it discards")
# ----------------------------------------------------------------------------------------------
a, b, c, th = sp.symbols("a b c theta", real=True)
check(sp.expand(2 * (a ** 2 + b ** 2) - (a + b) ** 2 - (a - b) ** 2) == 0, "AM-QM identity")
check(sp.expand((a + b) ** 2 - (a ** 2 + 2 * a * b + b ** 2)) == 0 and
      sp.expand((a - b) ** 2 - (a ** 2 - 2 * a * b + b ** 2)) == 0, "cross terms")
print("  Lemma 2.1: 2 min(theta,0) (a^2 + b^2) <= theta c^2, with a = d(x,z), b = d(z,y), c = d(x,y).")
print("   * theta >= 0: it reads 0 <= theta c^2 -- the LOWER half of the sandwich, (a-b)^2 <= c^2, is thrown away;")
print("   * theta < 0: it reads c^2 <= 2(a^2 + b^2), i.e. c^2 <= (a+b)^2 (triangle) and (a+b)^2 <= 2(a^2+b^2)")
print("     (AM-QM, equivalently 2ab <= a^2 + b^2).  The cross term +-2ab is replaced by its sign-blind bound.")
print("   * equality (theta < 0) iff c = a + b and a = b: z is a metric midpoint of x and y.")


def lemma21_ok(a_, b_, c_, t_):
    return 2 * min(t_, 0) * (a_ * a_ + b_ * b_) <= t_ * c_ * c_ + 1e-9


def random_metric(n):
    # shortest-path metric of a random complete weighted graph
    W = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            W[i][j] = W[j][i] = random.uniform(0.1, 3.0)
    for k_ in range(n):
        for i in range(n):
            for j in range(n):
                if W[i][k_] + W[k_][j] < W[i][j]:
                    W[i][j] = W[i][k_] + W[k_][j]
    return W


ntest = 0
thetas = [-3.0, -1.0, -0.25, 0.0, 0.5, 2.0]
for trial in range(60):
    W = random_metric(7)
    for i in range(7):
        for j in range(7):
            for k_ in range(7):
                for t_ in thetas:
                    check(lemma21_ok(W[i][k_], W[k_][j], W[i][j], t_), "Lemma 2.1 on a random metric")
                    ntest += 1
for trial in range(20000):
    P = [(random.uniform(-5, 5), random.uniform(-5, 5)) for _ in range(3)]
    for norm in ("l2", "l1", "linf"):
        def dist(u, v):
            dx, dy = abs(u[0] - v[0]), abs(u[1] - v[1])
            return sqrt(dx * dx + dy * dy) if norm == "l2" else (dx + dy if norm == "l1" else max(dx, dy))
        for t_ in thetas:
            check(lemma21_ok(dist(P[0], P[2]), dist(P[2], P[1]), dist(P[0], P[1]), t_), "Lemma 2.1 in R^2")
            ntest += 1
print(f"  Lemma 2.1 checked on {ntest} (metric, triple, theta) instances: random 7-point path metrics, "
      "R^2 with l1, l2, l_inf.")
# the loss along Collatz triples (p, Tp, T^2 p): exact
mono = turn = 0
for p in range(2, 100001):
    q, r = Tb(p), Tb(Tb(p))
    d1, d2, d13 = abs(q - p), abs(r - q), abs(r - p)
    loss = 2 * (d1 * d1 + d2 * d2) - d13 * d13
    same = (p % 2) == (q % 2)
    check(not (same and d1 == d2), "Lemma 2.1 never sharp on positive orbits (needs monotone and d1 = d2)")
    if d1 == d2:
        check(p == 2, "d1 = d2 only at p = 2 (the turning triple 2, 1, 2 of the 2-cycle)")
    if same:
        check(d13 == d1 + d2 and loss == (d1 - d2) ** 2, "monotone triple: upper equality, loss (d1-d2)^2")
        mono += 1
    else:
        check(d13 == abs(d1 - d2) and loss == (d1 + d2) ** 2, "turning triple: lower equality, loss (d1+d2)^2")
        turn += 1
print(f"  Collatz triples (p, Tp, T^2 p), 2 <= p <= 10^5: {mono} monotone (upper equality, Lemma 2.1 loses "
      f"(d1-d2)^2), {turn} turning (lower equality, loses (d1+d2)^2).")
print("  Sharpness needs a monotone triple with d1 = d2; d1 = d2 happens only at p = 2 (the turning triple")
print("  2, 1, 2 of the 2-cycle), so Lemma 2.1 is never sharp along a positive Collatz orbit.  What it discards")
print("  is the orientation (turn / no turn) of the orbit together with the size mismatch d1 - d2.")

# ----------------------------------------------------------------------------------------------
hdr("I2  On the line every orbit triple is an equality case; the cross-term sign is side-blind")
# ----------------------------------------------------------------------------------------------
ntrip = 0
zero_steps = set()
for x0 in list(range(-20000, 0)) + list(range(1, 20001)):
    o = [x0]
    for _ in range(120):
        o.append(Tb(o[-1]))
    for n in range(1, len(o) - 1):
        s1, s2 = o[n] - o[n - 1], o[n + 1] - o[n]
        if s1 == 0 or s2 == 0:
            zero_steps.add(o[n] if s2 == 0 else o[n - 1])
            continue
        d13 = abs(o[n + 1] - o[n - 1])
        same_par = (o[n - 1] % 2) == (o[n] % 2)
        if same_par:
            check(d13 == abs(s1) + abs(s2), "same parity -> upper equality (+2ab)")
        else:
            check(d13 == abs(abs(s1) - abs(s2)), "parity change -> lower equality (-2ab)")
        check((s1 * s2 > 0) == same_par, "cross-term sign = parity agreement, both half-lines")
        ntrip += 1
check(zero_steps == {-1}, "zero steps only at the fixed point -1 (0 excluded from starts)")
print(f"  {ntrip} consecutive triples of orbits from [-2*10^4, 2*10^4] \\ {{0}} (120 steps each):")
print("  d(x_(n-1), x_(n+1)) = d1 + d2 exactly when x_(n-1), x_n have the same parity (cross term +2 d1 d2),")
print("  and |d1 - d2| exactly when the parity changes (cross term -2 d1 d2) -- on BOTH half-lines.")
print("  The rule depends only on the parity word, hence (Theorem 6 of the mod-192 note) it is side-blind.")
print("  The 'central part' (cross term 0 with both steps nonzero) never occurs on the line: a zero step")
print(f"  occurs only at a fixed point (here {sorted(zero_steps)}; 0 is the other one).  A strictly-inside")
print("  triangle needs a second dimension, where the middle case is orthogonality (section I5).")

# ----------------------------------------------------------------------------------------------
hdr("I3  The sign law is the pair of equality cases of the sandwich, with z = 0 in the middle")
# ----------------------------------------------------------------------------------------------


def c_word(w):
    """T^p(x) = (3^a x + c_w)/2^p on the cylinder of w (plus sheet); R_() = 0, R_(z e) = 3^e R_z + e 2^|z|"""
    R = 0
    for i, e in enumerate(w):
        R = 3 ** e * R + e * 2 ** i
    return R


cyl = {}
for p in range(1, 13):
    for r in range(2 ** p):
        x, w = r, []
        for _ in range(p):
            w.append(x % 2)
            x = Tb(x)
        cyl[(p, tuple(w))] = r
ncheck = 0
for (p, w), r in cyl.items():
    aa = sum(w)
    cw = c_word(w)
    for t in (-3, -2, -1, 0, 1, 2, 3):
        x = r + 2 ** p * t
        if x == 0:
            continue
        y = x
        for _ in range(p):
            y = Tb(y)
        check(2 ** p * y == 3 ** aa * x + cw, "affine branch")
        ym = -x
        for _ in range(p):
            ym = Tb(ym, -1)
        check(ym == -y, "transport T_+^p(x) = -T_-^p(-x)")
        X, Y = 3 ** aa * x, -cw           # points on the line; z = 0 in the middle
        if x > 0:
            check(abs(X - Y) == abs(X) + abs(Y), "x>0: upper equality, 0 between 3^a x and -c_w")
            check(abs(2 ** p * y) == 3 ** aa * abs(x) + cw, "x>0: |2^p T^p x| = 3^a|x| + c_w")
        else:
            check(abs(X - Y) == abs(abs(X) - abs(Y)), "x<0: lower equality")
            check(abs(2 ** p * y) == 3 ** aa * abs(x) - cw and 3 ** aa * abs(x) >= cw, "x<0: 3^a|x| - c_w >= 0")
        ncheck += 1
print(f"  all words of length p <= 12, 6-7 lifts each (both signs): {ncheck} checks of")
print("    2^p T^p(x) = 3^a x + c_w  and  |2^p T^p(x)| = 3^a|x| + sgn(x) c_w  (c_w >= 0),")
print("  i.e. with X = 3^a x, Y = -c_w and z = 0:  x > 0 <=> 0 lies between X and Y <=> UPPER equality")
print("  d(X,Y) = d(X,0) + d(0,Y) (cross term +2 3^a|x| c_w);  x < 0 <=> LOWER equality d(X,Y) = |d(X,0) - d(0,Y)|.")
print("  Transport T_+^p(-x) = -T_-^p(x) turns the negative half-line of 3n+1 into the positive half-line of")
print("  3n-1: the lower side of the sandwich IS the 3n-1 sheet.  The middle point of the owner's inequality")
print("  is literally 0, and the sign law says on which side of it the orbit's two ingredients sit.")

# ----------------------------------------------------------------------------------------------
hdr("I4  The central sheet b = 0 (Mahler's multiplicative map) and the Mahler bridge")
# ----------------------------------------------------------------------------------------------
xs_ = sp.symbols("x")
nw = 0
for (p, w), r in cyl.items():
    if p > 10:
        continue
    aa, cw = sum(w), c_word(w)
    Tp_ = (3 ** aa * xs_ + cw) / sp.Integer(2 ** p)
    Tm_ = (3 ** aa * xs_ - cw) / sp.Integer(2 ** p)
    T0_ = 3 ** aa * xs_ / sp.Integer(2 ** p)
    check(sp.simplify((Tp_ + Tm_) / 2 - T0_) == 0, "central = midpoint")
    check(sp.simplify((Tp_ - Tp_.subs(xs_, -xs_)) / 2 - T0_) == 0, "central = odd part")
    check(sp.expand(Tp_ ** 2 - (9 ** aa * xs_ ** 2 + 2 * 3 ** aa * cw * xs_ + cw ** 2) / sp.Integer(4 ** p)) == 0,
          "cross term")
    nw += 1
print(f"  for all {nw} words with p <= 10:  T_0^w = (T_+^w + T_-^w)/2 = odd part of T_+^w, and")
print("  (T_b^w x)^2 = (9^a x^2 + 2b 3^a c_w x + b^2 c_w^2)/4^p: the cross term has sign b*sgn(x).  So the three")
print("  sheets b = -1, 0, +1 are the left side, the centre and the right side of the owner's sandwich (x > 0).")
# odd branch = central map centred at -b
for b_ in (1, -1):
    for x in range(-999, 1000, 2):
        check(2 * (Tb(x, b_) + b_) == 3 * (x + b_), "odd branch: T_b(x) + b = (3/2)(x + b)")
print("  odd branch: T_b(x) + b = (3/2)(x + b) -- the central map x -> 3x/2 centred at the fixed point -b;")
print("  the even branch x -> x/2 is centred at 0.  The b = 0 sheet is where both centres coincide.")
# shadowing identity
for b_ in (1, -1):
    for x in range(-500, 501):
        y, aa, cn = x, 0, 0
        xi_prev = Fr(x)
        for n in range(1, 61):
            e_ = y % 2
            cn = 3 ** e_ * cn + e_ * 2 ** (n - 1)
            aa += e_
            y = Tb(y, b_)
            xi = x + b_ * Fr(cn, 3 ** aa)
            check(Fr(3 ** aa, 2 ** n) * xi == y, "shadowing identity T^n x = T_0^(w_n)(xi_n)")
            check((xi >= xi_prev) if b_ == 1 else (xi <= xi_prev), "xi_n monotone")
            xi_prev = xi
print("  shadowing identity (exact, x in [-500,500], n <= 60, both sheets): T_b^n(x) = (3^a_n/2^n) xi_n with")
print("  xi_n = x + b c_n/3^a_n; xi_n increases on the plus sheet and decreases on the minus sheet.  Every")
print("  orbit point is the CENTRAL-sheet image of a shifted start; the shift is the Bernstein partial sum.")
# THM-4469 instance
Bw = tuple(int(ch) for ch in "0111101110")
Bpw = tuple(int(ch) for ch in "1101100111")
RB, RBp = c_word(Bw), c_word(Bpw)
P3, Q2 = 3 ** 7, 2 ** 10
check(RB == 4726 and RBp == 4727 and sum(Bw) == 7 and sum(Bpw) == 7 and P3 - Q2 == 1163, "THM-4469 carries")


def realize(word):
    """least positive integer whose parity vector starts with word (Terras class via inverse branches)"""
    L = len(word)
    mod = 2 ** L
    inv3 = pow(3, -1, mod)
    y = 0
    for e_ in reversed(word):
        y = (2 * y) % mod if e_ == 0 else ((2 * y - 1) * inv3) % mod
    return y if y > 0 else mod


check(realize(Bw) == 990 and realize(Bpw) == 187, "THM-4469 cylinder classes 990, 187 mod 1024")
alpha, rho = Fr(P3, Q2), Fr(Q2, P3)
for trial in range(40):
    J = 6
    blocks = [random.choice((Bw, Bpw)) for _ in range(J)]
    word = tuple(e_ for blk in blocks for e_ in blk)
    x0 = realize(word)
    xs_orbit = [x0]
    y = x0
    for j in range(J):
        for _ in range(10):
            y = Tb(y)
        xs_orbit.append(y)
    Rs = [RB if blk == Bw else RBp for blk in blocks]
    t = [sum(Fr(Rs[j + i], P3) * rho ** i for i in range(J - j)) for j in range(J)] + [Fr(0)]
    xi = x0 + t[0]
    for j in range(J + 1):
        check(xi * alpha ** j == xs_orbit[j] + t[j], "Mahler shadowing x_j = xi alpha^j - t_j")
        check(Fr(RB, P3 - Q2) * (1 - rho ** (J - j)) <= t[j] <= Fr(RBp, P3 - Q2) * (1 - rho ** (J - j)),
              "carry confined to I = [4726,4727]/1163")
    # the minus sheet: same blocks, mirrored start, x_j^- = -x_j, xi_- = -xi
    ym = -x0
    for j in range(J):
        for _ in range(10):
            ym = Tb(ym, -1)
        check(ym == -xs_orbit[j + 1] and (-xi) * alpha ** (j + 1) == ym - t[j + 1], "minus-sheet shadow")
print("  THM-4469 instance (blocks 0111101110 / 1101100111, R = 4726 / 4727, cylinder classes 990 / 187 mod 1024,")
print("  alpha = 2187/1024): for 40 random block sequences of length 6, realized by the least positive integer")
print("  of their Terras class,")
print("  x_j = xi alpha^j - t_j with t_j in [4726,4727]/1163 * (1 - rho^(6-j)) exactly: the +1 orbit is the")
print("  central-sheet (Mahler) orbit xi alpha^j pushed DOWN by a confined carry; the 3n-1 orbit of -x is the")
print("  central orbit of -xi pushed UP by the same carry.  This is the precise sense in which the Mahler bridge")
print("  says '+-1-corner orbits are confined like central-sheet orbits'.")

# ----------------------------------------------------------------------------------------------
hdr("I5  Orthogonality in C: the Pythagorean middle, Chamberland's f and the LSW f_0")
# ----------------------------------------------------------------------------------------------
for trial in range(20000):
    z = complex(random.uniform(-9, 9), random.uniform(-9, 9))
    w_ = complex(random.uniform(-9, 9), random.uniform(-9, 9))
    lhs = abs(z + w_) ** 2
    rhs = abs(z) ** 2 + abs(w_) ** 2 + 2 * (z * w_.conjugate()).real
    check(abs(lhs - rhs) < 1e-9 * (1 + lhs), "|z+w|^2 identity")
    y_ = random.uniform(-9, 9)
    cw_ = random.randint(1, 1000)
    a3 = 3 ** random.randint(0, 6)
    zz = complex(0, y_)
    check(abs(abs(a3 * zz + cw_) ** 2 - (a3 ** 2 * y_ ** 2 + cw_ ** 2)) < 1e-6 * (1 + abs(a3 * zz + cw_) ** 2),
          "Pythagorean on the imaginary axis")
print("  |z + w|^2 = |z|^2 + |w|^2 + 2 Re(z conj w): the middle (Pythagorean) case is Re(z conj w) = 0.  For an")
print("  affine branch L_w(z) = (3^a z + c_w)/2^p the cross term is 2 3^a c_w Re(z): right half-plane = '+2ab',")
print("  left half-plane = '-2ab', imaginary axis = the orthogonal middle.  (20000 random checks.)")
mp.mp.dps = 40
pi = mp.pi


def f_ch(x):
    return x + mp.mpf(1) / 4 - (2 * x + 1) / 4 * mp.cos(pi * x)


def f_lsw(z):
    return z / 2 + (1 - mp.cos(pi * z)) * (z + mp.mpf(1) / 2) / 2 + (mp.mpf(1) / 2 - mp.cos(pi * z)) * mp.sin(pi * z) / pi


for n in range(-30, 31):
    check(abs(f_ch(n) - Tb(n)) < mp.mpf(10) ** -30, "Chamberland f = T on Z")
    check(abs(mp.diff(f_ch, n) - (mp.mpf(1) / 2 if n % 2 == 0 else mp.mpf(3) / 2)) < mp.mpf(10) ** -25,
          "Chamberland f'(n) = branch slope")
    check(abs(f_lsw(n) - Tb(n)) < mp.mpf(10) ** -30, "LSW f_0 = T on Z")
    check(abs(mp.diff(f_lsw, n)) < mp.mpf(10) ** -25, "LSW f_0'(n) = 0")
print("  Chamberland f(x) = x + 1/4 - ((2x+1)/4) cos(pi x): f = T on Z and f'(n) = 1/2 (even), 3/2 (odd) --")
print("  f keeps the branch slopes, so it is conformal at integers (orthogonal stays orthogonal).")
print("  LSW f_0(z) = z/2 + (1 - cos pi z)(z + 1/2)/2 + (1/2 - cos pi z) sin(pi z)/pi (h = 0, as printed in")
print("  Chamberland's survey): f_0 = T on Z and f_0'(n) = 0 -- every integer is a critical point.  (n in [-30,30].)")
# f_0'(x) = sin(pi x) (2 sin(pi x) + (pi/2)(x + 1/2))
for trial in range(200):
    xv = mp.mpf(random.uniform(-6, 6))
    check(abs(mp.diff(f_lsw, xv) - mp.sin(pi * xv) * (2 * mp.sin(pi * xv) + pi / 2 * (xv + mp.mpf(1) / 2)))
          < mp.mpf(10) ** -20, "closed form of f_0'")
q = lambda x: 2 * mp.sin(pi * x) + pi / 2 * (x + mp.mpf(1) / 2)
roots = []
grid = [mp.mpf(i) / 1000 for i in range(-4000, 4001)]
for u, v in zip(grid[:-1], grid[1:]):
    if q(u) * q(v) < 0:
        roots.append(mp.findroot(q, (u, v), solver="bisect"))
check(len(roots) == 3 and all(-2 < r_ < 0 for r_ in roots), "three non-integer real critical points")
check(abs(mp.diff(f_lsw, mp.mpf(-1) / 2) - 2) < mp.mpf(10) ** -25, "f_0'(-1/2) = 2")
print("  f_0'(x) = sin(pi x) (2 sin(pi x) + (pi/2)(x + 1/2)) (checked); its real zeros are Z together with")
print(f"  {', '.join(mp.nstr(r_, 10) for r_ in roots)} (only |x + 1/2| <= 4/pi can solve the second factor).")
print("  NOTE: Lagarias's annotated bibliography (entry 115) summarizes the real critical points of f_0 as")
print("  'the integers together with -1/2'; for the formula as printed, f_0'(-1/2) = 2, so this could not be")
print("  reconciled (primary source not read: publisher bot-check).  Nothing below depends on it.")
# the fold: orthogonal -> opposite at a critical integer
for n in (-5, -2, 1, 2, 7):
    c2 = mp.diff(f_lsw, n, 2) / 2
    eps = mp.mpf(10) ** -8
    re_ = (f_lsw(n + eps) - Tb(n)) / eps ** 2
    im_ = (f_lsw(n + 1j * eps) - Tb(n)) / eps ** 2
    check(abs(re_ - c2) < mp.mpf(10) ** -6 and abs(im_ + c2) < mp.mpf(10) ** -6 and abs(c2) > 0.1,
          "fold: real direction -> +c y^2, imaginary -> -c y^2")
print("  At each critical integer n, f_0(n + y) - T(n) ~ c_n y^2 and f_0(n + iy) - T(n) ~ -c_n y^2 (c_n real,")
print("  checked at n = -5,-2,1,2,7): the fold sends the ORTHOGONAL direction onto the OPPOSITE real ray.  This is")
print("  the exact content of 'the corner at 0 introducing orthogonality': squaring (the smooth version of |x|)")
print("  identifies +-x and turns the orthogonal axis into the sign.")
print("  Stability consequence: in Chamberland's f an integer cycle has multiplier 3^a/2^p, so (sign law) it is")
print("  attracting iff its points are >= 0; in LSW's f_h every integer cycle is superattracting on BOTH sides.")

# ----------------------------------------------------------------------------------------------
hdr("I6  |T(x)| = T_sgn(x)(|x|); linear pair sums, quadratic pair products, and the cross term")
# ----------------------------------------------------------------------------------------------
xs_all = np.arange(-10 ** 6, 10 ** 6 + 1, dtype=np.int64)
Tplus = np.where(xs_all % 2 == 0, xs_all // 2, (3 * xs_all + 1) // 2)
absx = np.abs(xs_all)
sg = np.sign(xs_all)
Tsg = np.where(absx % 2 == 0, absx // 2, (3 * absx + sg) // 2)
check(np.array_equal(np.abs(Tplus), Tsg), "|T_+(x)| = T_sgn(x)(|x|)")
del xs_all, Tplus, absx, sg, Tsg
for x in range(-10000, 10001):
    y, z_ = x, abs(x)
    for _ in range(50):
        y, z_ = Tb(y), Tb(z_, sgn(x)) if x != 0 else 0
        check(abs(y) == z_, "iterated |T_+^p(x)| = T_sgn(x)^p(|x|)")
print("  |T_+(x)| = T_sgn(x)(|x|) for all |x| <= 10^6 (T_0(0) = 0), and iterated for p <= 50, |x| <= 10^4.")
print("  On the size coordinate |x| the single map T_+ becomes the PAIR T_+1, T_-1, selected by sgn x = d|x|/dx:")
print("  on odd x, |T(x)| = (3|x| + sgn x)/2.  The corner of |x| at 0 is exactly where the selector jumps.")
for i in range(1, 10 ** 6 + 1):
    aa_, bb_ = 3 * i - 1, i
    x1, x2 = 2 * i - 1, 2 * i
    if i <= 20000:
        check(Fr(aa_ * bb_, x1 * x2) == Fr(3, 4) + Fr(1, 8 * i - 4), "product ratio 3/4 + 1/(8i-4)")
    check(aa_ + bb_ == x1 + x2, "pair sum conserved")
    check(2 * aa_ * bb_ - 2 * x1 * x2 == -2 * i * (i - 1), "cross term drops by 2i(i-1)")
    check(aa_ ** 2 + bb_ ** 2 - x1 ** 2 - x2 ** 2 == 2 * i * (i - 1), "squares gain 2i(i-1)")
print("  THM-4470 pairs {2i-1, 2i} -> {3i-1, i}, i <= 10^6:  (a+b) conserved (LINEAR invariant),")
print("  ab scaled by 3/4 + 1/(8i-4) (QUADRATIC contraction), and exactly  D(2ab) = -2i(i-1) = -D(a^2+b^2):")
print("  in (a+b)^2 = a^2 + 2ab + b^2 the left side is fixed and T moves mass from the cross term to the squares.")
print("  Equivalently QM^2 + GM^2 = 2 AM^2 with AM fixed: T lowers GM and raises QM by the same amount.  Lemma 2.1")
print("  keeps only QM >= AM; the Collatz drift log(2/sqrt3) (the AM-GM gap) lives in the discarded cross term.")
print("  Even observables (|x|, x^2, |x - y|) are invariant under x -> -x, hence side-blind (Theorem 6 of the")
print("  mod-192 note); only odd observables (x, sgn x, ...) can carry the sign law.")

# ----------------------------------------------------------------------------------------------
hdr("I7  The fixed points {0, -b}: the centres of the two branches; the self-mapped pairs")
# ----------------------------------------------------------------------------------------------
for b_ in (1, -1):
    fx = [x for x in range(-10 ** 5, 10 ** 5 + 1) if Tb(x, b_) == x]
    check(fx == sorted({0, -b_}), f"Fix(T_b) = {{0,-b}} for b = {b_}")
    for x in range(-3000, 3001):
        if (2 * x - b_) % 3 == 0:
            check(3 * ((2 * x - b_) // 3 + b_) == 2 * (x + b_), "E_b(x) + b = (2/3)(x + b)")
print("  Fix(T_b) on Z is {0, -b} (checked on [-10^5, 10^5]); D(0) = 0 and E_b(-b) = -b are the Banach fixed")
print("  points of the one-letter inverse words (chains script F1); E_b(x) + b = (2/3)(x + b).")
selfp = [i for i in range(-10 ** 5, 10 ** 5 + 1) if {Tb(2 * i - 1), Tb(2 * i)} == {2 * i - 1, 2 * i}]
check(selfp == [0, 1], "self-mapped plus pairs: i = 0, 1")
selfm = [i for i in range(-10 ** 5, 10 ** 5 + 1) if {Tb(2 * i, -1), Tb(2 * i + 1, -1)} == {2 * i, 2 * i + 1}]
check(selfm == [-1, 0], "self-mapped minus pairs: i = -1, 0")
print("  Consecutive pairs of 3n+1, {2i-1, 2i}, mapped onto themselves (|i| <= 10^5): i = 0, the pair {-1, 0}")
print("  of the two fixed points (pointwise), and i = 1, {1, 2} (swapped: the 2-cycle).  For 3n-1, {2i, 2i+1}:")
print("  {0, 1} (pointwise) and {-2, -1} (swapped) -- the mirror images.  The 'central pair' straddling the")
print("  sign is {-1, 0} for 3n+1 and {0, 1} for 3n-1: b decides on which side of 0 the odd branch is centred.")

rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
rss_mb = rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024
print()
print(f"[fixedpt_inequality] ALL CHECKS PASSED  ({time.time() - T0:.1f} s, peak RSS {rss_mb:.0f} MB)",
      file=sys.stderr)
print("[fixedpt_inequality] ALL CHECKS PASSED")

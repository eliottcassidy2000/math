#!/usr/bin/env python3
"""Part 4: Collatz's own repunits -- the dictionary, the 2-adic macrocosm, the Bernstein--Lagarias
fixed points, and the test of "the length of all-one primes".

 (4a) Exact identities (checked exhaustively in the stated ranges):
        odd branch = dilation about -1:  T(x) + 1 = (3/2)(x + 1) for odd x; the initial run of odd steps
        of n has length v_2(n+1);  T^j(2^k-1) = 3^j 2^(k-j) - 1 (j <= k), T^k = 3^k - 1,
        T^(k+1)(2^k-1) = (3^k-1)/2;  trunk (4^k-1)/3: T = 2^(2k-1), sigma = 2, tau = 2k;
        Wagstaff/minus-trunk (2^k+1)/3, k odd: T = 2^(k-1)+1, T^(1+2i) = 3^i 2^(k-1-2i)+1, sigma = 5 (k >= 5),
        minus sheet T_-(W) = 2^(k-1); base-3 repunits: sigma = 1 (k even), 2 (k odd);
        base-10 repunits: T^(2+2i)(R_k) = 3^i 5^k 2^(k-2-2i) + 1;
        the SHADOWING LEMMA T^j(R_k^(b)) = T^j(xi_b) + 3^(a_j(xi_b)) b^k / ((b-1) 2^j), j <= k v_2(b),
        xi_b = 1/(1-b), for even b.
 (4b) 2-adic orbits of the repunit limits xi_b = 1/(1-b), even b in [-64, 64], on both sheets.
 (4c) The Bernstein--Lagarias conjugacy Phi = Q^(-1) (Q = parity-vector map): Q(-1) = -1, Q(1/3) = 1/3,
      Q(1) = -1/3, Q(-1/3) = 1; search of small-height odd rationals for fixed points and 2-cycles of Q.
 (4d) Rotation <-> Collatz cycles: rational cycles are shift orbits of parity words; integral cycles of
      period <= 18 (the "circular primes" of Collatz); the finiteness heuristics side by side.
 (4e) The 3-adic side: base-3 repunits -> -1/2 in Z_3, the fixed point of x -> 3x+1.
 (4f) "The length of all-one primes": Collatz statistics of 2^k-1 (and of the other repunit families)
      as functions of k, exponent families vs other lengths (permutation tests), 2-adic continuity in k.
 (4g) The convergent-numerator coincidence (2, 19, 317 are best-approximation numerators of log2 3).
Needs the GMP helper procgen_repunit_20260924_collatz_repunits.c.  Runtime ~ 3 min, memory < 300 MB.
Session collatz-procgen-20260922, lane procgen_repunit_20260924.
"""
import os, sys, subprocess, tempfile, math, time, random
from fractions import Fraction as F
from itertools import product
import numpy as np
from scipy import stats
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="procgen_repunit_dict_")
BIN = os.path.join(TMP, "crep")
KM = int(os.environ.get("PR_DICT_KM", "12000"))       # Mersenne family range
KW = int(os.environ.get("PR_DICT_KW", "12000"))       # Wagstaff family range
K3 = int(os.environ.get("PR_DICT_K3", "6000"))        # base-3 repunit family range
K10 = int(os.environ.get("PR_DICT_K10", "3000"))      # base-10 repunit family range

# exponent lists (OEIS, read 2026-09-24)
A000043 = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941,
           11213, 19937, 21701, 23209, 44497]
A000978 = [3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101, 127, 167, 191, 199, 313, 347, 701, 1709, 2617, 3539, 5807,
           10501, 10691, 11279, 12391, 14479, 42737]
A028491 = [3, 7, 13, 71, 103, 541, 1091, 1367, 1627, 4177, 9011, 9551, 36913]
A004023 = [2, 19, 23, 317, 1031, 49081, 86453, 109297, 270343, 5794777, 8177207]


def v2(n):
    n = abs(n)
    return (n & -n).bit_length() - 1 if n else 10**9


def T(x, sheet=1):
    """shortcut map on Z_(2) (rationals with odd denominator); sheet=+1: 3x+1, sheet=-1: 3x-1"""
    x = F(x)
    if x.numerator % 2 == 0:
        return x / 2
    return (3 * x + sheet) / 2


def parity(x):
    return F(x).numerator % 2


def orbit_until_periodic(x, sheet=1, cap=20000):
    seen = {}
    seq = []
    x = F(x)
    while x not in seen:
        if len(seq) > cap:
            return None
        seen[x] = len(seq)
        seq.append(x)
        x = T(x, sheet)
    m = seen[x]
    return seq[:m], seq[m:]                    # preperiodic part, cycle


def Qmap(x, sheet=1):
    """Q(x) = sum_i x_i 2^i, x_i = parity of T^i(x), for x with an eventually periodic orbit (exact)"""
    r = orbit_until_periodic(x, sheet)
    if r is None:
        return None
    pre, cyc = r
    u = [parity(y) for y in pre]
    w = [parity(y) for y in cyc]
    val = F(sum(b << i for i, b in enumerate(u)))
    cw = sum(b << i for i, b in enumerate(w))
    val += F(2 ** len(u) * cw, 1 - 2 ** len(w))
    return val


def word(x, n, sheet=1):
    out = []
    x = F(x)
    for _ in range(n):
        out.append(parity(x))
        x = T(x, sheet)
    return "".join(map(str, out))


# ------------------------------------------------------------------------------------------ (4a)
def part4a():
    print("\n(4a) Exact identities of the Collatz repunits")
    # dilation about -1 and run length
    for n in range(1, 200001):
        r = v2(n + 1)
        x = n
        for j in range(r):
            assert x % 2 == 1
            assert x + 1 == F(3, 2) ** j * (n + 1)
            x = (3 * x + 1) // 2
        assert x % 2 == 0 or r == 0
    print("    for every n <= 2*10^5: the initial run of odd T-steps has length exactly v_2(n+1), and along it"
          " T^j(n) + 1 = (3/2)^j (n+1)  [the odd branch is the dilation by 3/2 about the fixed point -1].")
    for k in range(1, 401):
        n = 2 ** k - 1
        x = n
        for j in range(k + 1):
            assert x == 3 ** j * 2 ** (k - j) - 1
            if j < k:
                x = int(T(x))
        assert x == 3 ** k - 1 and int(T(x)) == (3 ** k - 1) // 2
        assert v2(3 ** k - 1) == (1 if k % 2 else v2(k) + 2)
        R3 = (3 ** k - 1) // 2
        s = next(j for j in range(1, 10) if int(_Tpow(R3, j)) < R3) if k >= 2 else None
        if k >= 2:
            assert s == (1 if k % 2 == 0 else 2)
    print("    k <= 400: T^j(2^k-1) = 3^j 2^(k-j) - 1 for j <= k, T^(k+1)(2^k-1) = (3^k-1)/2 (the base-3 repunit);"
          " v_2(3^k-1) = 1 (k odd), v_2(k)+2 (k even) [lifting the exponent]; sigma((3^k-1)/2) = 1 (k even), 2 (k odd).")
    for k in range(2, 301):
        R4 = (4 ** k - 1) // 3
        assert int(T(R4)) == 2 ** (2 * k - 1)
        xs = _orbit(R4)
        assert len(xs) - 1 == 2 * k and max(xs) == 2 ** (2 * k - 1)
        if k >= 2:
            assert next(j for j in range(1, len(xs)) if xs[j] < R4) == 2
    print("    2 <= k <= 300: trunk (4^k-1)/3: T = 2^(2k-1), total stopping time 2k, height 2^(2k-1), stopping time 2.")
    for k in range(3, 302, 2):
        W = (2 ** k + 1) // 3
        assert int(T(W)) == 2 ** (k - 1) + 1
        x = int(T(W))
        for i in range(0, (k - 1) // 2 + 1):
            assert x == 3 ** i * 2 ** (k - 1 - 2 * i) + 1
            if i < (k - 1) // 2:
                x = int(T(T(x)))
        xs = _orbit(W)
        s = next(j for j in range(1, len(xs)) if xs[j] < W)
        assert s == (4 if k == 3 else 5)
        assert int(T(W, -1)) == 2 ** (k - 1)
    print("    odd k <= 301: W = (2^k+1)/3: T(W) = 2^(k-1)+1, T^(1+2i)(W) = 3^i 2^(k-1-2i) + 1 (so T^k(W) = 3^((k-1)/2)+1),"
          " stopping time 5 (k >= 5; 4 at k = 3); on the 3n-1 sheet T_-(W) = 2^(k-1): the minus trunk.")
    for j in range(1, 150):
        x = -(4 ** j - 1) // 3
        assert T(x, -1) == -(2 ** (2 * j - 1))
    print("    even length: R_(2j)^(-2) = -(4^j-1)/3 and T_-(R) = -2^(2j-1): the negative trunk of the 3n-1 sheet"
          " (the nu-image of the 3n+1 trunk).")
    for k in range(2, 301):
        R = (10 ** k - 1) // 9
        x = R
        for j in range(1, k + 1):
            x = int(T(x))
            if j >= 2 and j % 2 == 0 and j <= k:
                i = (j - 2) // 2
                assert x == 3 ** i * 5 ** k * 2 ** (k - 2 - 2 * i) + 1
    print("    k <= 300: base-10 repunits T^(2+2i)(R_k) = 3^i 5^k 2^(k-2-2i) + 1 (2+2i <= k): the (10)-pattern that shadows"
          " the cycle {1,2}.")
    # shadowing lemma for even bases
    nb = 0
    for b in (2, -2, 4, -4, 6, -6, 8, -8, 10, -10, 12, 16, -16, 18, 20, 24, 32, 100):
        xi = F(1, 1 - b)
        for k in range(1, 61):
            R = F(b ** k - 1, b - 1)
            J = k * v2(b)
            x, y, a = R, xi, 0
            for j in range(J + 1):
                assert x == y + F(3 ** a * b ** k, (b - 1) * 2 ** j), (b, k, j)
                if j < J:
                    assert parity(x) == parity(y)
                    a += parity(y)
                    x, y = T(x), T(y)
            nb += 1
    print(f"    SHADOWING LEMMA checked for 18 even bases b and k <= 60 ({nb} cases): the first k*v_2(b) parity bits of"
          " R_k^(b) are those of xi_b = 1/(1-b), and T^j(R_k^(b)) = T^j(xi_b) + 3^(a_j) b^k/((b-1) 2^j).")


def _Tpow(x, j):
    for _ in range(j):
        x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
    return x


def _orbit(n):
    xs = [n]
    while xs[-1] != 1:
        x = xs[-1]
        xs.append(x // 2 if x % 2 == 0 else (3 * x + 1) // 2)
    return xs


# ------------------------------------------------------------------------------------------ (4b)
def part4b():
    print("\n(4b) The 2-adic limits xi_b = 1/(1-b) of the base-b repunits (b even) and their Collatz orbits")
    print("     b   xi_b     sheet  preperiod period  cycle (as rationals)                     parity word (pre|cycle)")
    periodic_b = []
    for b in range(-64, 65, 2):
        if b == 0:
            continue
        xi = F(1, 1 - b)
        for sheet in (1, -1):
            pre, cyc = orbit_until_periodic(xi, sheet)
            w = "".join(str(parity(y)) for y in pre) + "|" + "".join(str(parity(y)) for y in cyc)
            if abs(b) <= 16 or (sheet == 1 and len(pre) == 0):
                cs = ", ".join(str(c) for c in cyc[:6]) + (" ..." if len(cyc) > 6 else "")
                print(f"    {b:>4}  {str(xi):>7}   {'+' if sheet == 1 else '-'}   {len(pre):>5} {len(cyc):>6}   {cs:<40} {w[:40]}")
            if sheet == 1 and len(pre) == 0:
                periodic_b.append(b)
    print(f"    even b in [-64, 64] with xi_b T_+-periodic: {periodic_b}")
    # b = 4 - 2^p gives xi = 1/(2^p - 3), the cycle of the word 1 0^(p-1)
    for p in range(1, 12):
        b = 4 - 2 ** p
        if b == 0:
            continue
        xi = F(1, 1 - b)
        pre, cyc = orbit_until_periodic(xi)
        assert len(pre) == 0 and "".join(str(parity(y)) for y in cyc) == "1" + "0" * (p - 1)
    print("    b = 4 - 2^p (p >= 1: 2, -4, -12, -28, ...) makes xi_b = 1/(2^p - 3) the periodic point of the word 1 0^(p-1)"
          " (checked p <= 11); b = 2 is the fixed point -1 (word 1), and b = 0 would be the cycle {1, 2}.")
    print("    Every orbit above is eventually periodic (FINITE-EXACT), as the Periodicity Conjecture predicts for all"
          " rationals with odd denominator (open in general).")
    # the Moebius action on bases: T(xi_b) = xi_{M(b)}, M(b) = -(b+2)/(b-4)
    for b in range(-40, 41, 2):
        if b in (0, 4):
            continue
        xi = F(1, 1 - b)
        Mb = F(-(b + 2), b - 4)
        assert T(xi) == 1 / (1 - Mb)
    print("    T acts on the base: T(xi_b) = xi_M(b) with M(b) = -(b+2)/(b-4): 10 -> -2 -> 0 (then the cycle {1,2}),"
          " 6 -> -4 (periodic), 2 -> 2 (fixed), 4 -> infinity (xi = 0).")


# ------------------------------------------------------------------------------------------ (4c)
def part4c():
    print("\n(4c) The Bernstein--Lagarias conjugacy Phi = Q^(-1): its known odd fixed points are 2-adic repunits")
    checks = [(F(-1), F(-1)), (F(1, 3), F(1, 3)), (F(1), F(-1, 3)), (F(-1, 3), F(1))]
    for x, y in checks:
        assert Qmap(x) == y, (x, Qmap(x))
    print("    Q(-1) = -1, Q(1/3) = 1/3, Q(1) = -1/3, Q(-1/3) = 1 (exact).  -1 = ...1111 (base 2), 1/3 = ...1111 (base -2),")
    print("    -1/3 = ...1111 (base 4) and 1 = R_1^(4): the two known odd fixed points and the known odd 2-cycle of Phi")
    print("    [Bernstein--Lagarias 1996, Fixed Point Conjecture: exactly two odd fixed points; 2-cycle {1, -1/3}]")
    print("    are the 2-adic repunits of bases 2 and -2, and the two ends of the base-4 trunk.")
    print("    parity word of 1/3: " + word(F(1, 3), 24) + "   2-adic digits of 1/3 (low first): " +
          "".join(str((F(1, 3) * 1).numerator * pow(3, -1, 2 ** 24) % 2 ** 24 >> i & 1) for i in range(24)))
    for m in range(1, 6):
        for x in (F(-1), F(1, 3)):
            assert Qmap(2 ** m * x) == 2 ** m * x
    print("    Q(2x) = 2Q(x), so the even fixed points are 0, -2^m, 2^m/3 (checked m <= 5).")
    # search small-height odd rationals
    t0 = time.time()
    fixed, two = [], []
    H = int(os.environ.get("PR_DICT_QH", "201"))
    for q in range(1, 100, 2):
        for p in range(-H, H + 1, 2):
            if math.gcd(p, q) != 1:
                continue
            x = F(p, q)
            y = Qmap(x)
            if y is None:
                continue
            if y == x:
                fixed.append(x)
            elif y.denominator % 2 == 1 and y.denominator < 10 ** 6:
                z = Qmap(y)
                if z == x and x < y:
                    two.append((x, y))
    print(f"    search: odd x = p/q, |p| <= {H}, odd q < 100 ({time.time() - t0:.0f} s): Q-fixed points {fixed};"
          f" Q-2-cycles {two}")
    assert sorted(fixed) == [F(-1), F(1, 3)]
    assert sorted(two) == [(F(-1, 3), F(1)), (F(-1, 5), F(5, 7))]
    # independent check of every cycle found: first 80 parity bits from x mod 2^160 (Terras), compared mod 2^80
    NB = 80

    def Qmod(x):
        Mo = 1 << (2 * NB)
        r = x.numerator * pow(x.denominator, -1, Mo) % Mo
        v = 0
        for i in range(NB):
            b = r & 1
            v |= b << i
            r = (r >> 1) if b == 0 else ((3 * r + 1) >> 1)
        return v

    def res(x):
        return x.numerator * pow(x.denominator, -1, 1 << NB) % (1 << NB)
    for x in fixed:
        assert Qmod(x) == res(x)
    for x, y in two:
        assert Qmod(x) == res(y) and Qmod(y) == res(x)
    print("    independent check (first 80 parity bits from residues mod 2^160): every fixed point and 2-cycle confirmed.")
    print("    The fixed points reproduce Bernstein--Lagarias's search.  The 2-cycle {-1/5, 5/7} is NOT among the cycles")
    print("    listed by B--L 1996 ('We know of one odd periodic cycle of Phi of length 2, namely {1, -1/3}'), so it gives")
    print("    F_1 >= 2 in their notation; whether later work (e.g. Hotzel 2003, Ch. 7, periodic points of the conjugacy")
    print("    map -- not read) records it is UNVERIFIED.  Its point -1/5 = ...1111 (base 6) is again a repunit limit:")
    print("    Q(-1/5) = 5/7 because -1/5 -> 1/5 -> 4/5 -> 2/5 -> 1/5 (word 1(100)^inf), and 5/7 -> 11/7 -> 20/7 -> 10/7 -> 5/7"
          " (word (1100)^inf) gives Q(5/7) = 3/(1-16) = -1/5.")
    # which repunit limits are Q-fixed / 2-periodic
    res = []
    for b in range(-64, 65, 2):
        xi = F(1, 1 - b)
        y = Qmap(xi)
        tag = "fixed" if y == xi else ("2-cycle" if Qmap(y) == xi else "")
        if tag:
            res.append((b, str(xi), tag))
    print(f"    among the repunit limits xi_b (even b in [-64,64], and b = 0 -> xi = 1): {res}")


# ------------------------------------------------------------------------------------------ (4d)
def part4d():
    print("\n(4d) Rotation <-> Collatz: rational cycles are shift orbits of parity words")
    print("    The shift S on Z_2 (drop the last binary digit) acts on the period-k points -B/(2^k-1) (B = the k-bit")
    print("    block) as rotation of the block: base-2 digit rotation.  T = Phi S Phi^-1, so the T-cycle of a primitive")
    print("    word w of length p is the image of the rotation orbit of w; its points are x_w = c_w/(2^p - 3^a).")
    print("    Fixed points of rotation (blocks 0^k, 1^k) <-> the two fixed points 0, -1 of T.  For prime p there are")
    print("    (2^p - 2)/p cycles of exact period p (Fermat, necklace form).  'All rotations prime' is not rotation-")
    print("    invariant; 'integral' is (T maps Z to Z along the cycle), so integral cycles are the Collatz analogue.")
    t0 = time.time()
    integral = []
    PMAX = int(os.environ.get("PR_DICT_PMAX", "18"))
    for p in range(1, PMAX + 1):
        seen = set()
        for bits in range(2 ** p):
            w = [(bits >> i) & 1 for i in range(p)]
            # primitive words only, one per rotation class
            rots = [tuple(w[i:] + w[:i]) for i in range(p)]
            key = min(rots)
            if key in seen:
                continue
            seen.add(key)
            if len(set(rots)) != p:
                continue
            a = sum(w)
            # T^p(x) = (3^a x + c_w)/2^p, c_w = sum over odd steps i of 3^(#odd steps after i) 2^i
            c = 0
            ones = 0
            for i, bit in enumerate(w):
                if bit:
                    ones += 1
                    c += 3 ** (a - ones) * 2 ** i
            den = 2 ** p - 3 ** a
            if den != 0 and c % den == 0:
                x = F(c, den)
                assert _Tfrac_pow(x, p) == x
                integral.append((p, "".join(map(str, w)), int(x)))
    print(f"    integral periodic points with a primitive word of length <= {PMAX} ({time.time() - t0:.0f} s):")
    for p, w, x in integral:
        print(f"      period {p:>2}: word {w:<20} point {x}")
    assert sorted(x for _, _, x in integral) == ([-17, -5, -1, 0, 1] if PMAX >= 11 else [-5, -1, 0, 1]), integral
    print("    = the five known integer cycles 0, -1, {1,2}, {-5,-7,-10}, {-17,...,-91} (one point per cycle), and no other"
          " integral cycle of period <= %d (FINITE-EXACT; the known cycle-length bounds are far larger)." % PMAX)
    print("    Heuristic count of integral cycles of length p: sum_a C(p,a)/(p |2^p - 3^a|)  vs circular primes:")
    tot = 0.0
    for p in range(2, 61):
        e = sum(math.comb(p, a) / abs(2 ** p - 3 ** a) for a in range(1, p) if 2 ** p != 3 ** a) / p
        tot += e
        if p in (5, 10, 20, 30, 40, 50, 60):
            print(f"      p = {p:>2}: expected integral p-cycles {e:.3e}")
    print("    Both expectations decay (Collatz like 2^(-(1-h(log_3 2))p) = 2^(-0.05 p) up to Diophantine spikes at")
    print("    convergents of log_2 3; circular primes like (6.5/k)^k), so both 'finitely many' statements are")
    print("    heuristically natural; neither is proved (the Collatz one is a Diophantine problem at 2^p - 3^a).")


def _Tfrac_pow(x, p):
    for _ in range(p):
        x = T(x)
    return x


# ------------------------------------------------------------------------------------------ (4e)
def part4e():
    print("\n(4e) The 3-adic side")
    for k in range(1, 200):
        R3 = F(3 ** k - 1, 2)
        d = R3 + F(1, 2)
        assert d == F(3 ** k, 2)
        assert 3 * ((3 ** (k - 1) - 1) // 2) + 1 == (3 ** k - 1) // 2 if k >= 1 else True
    print("    (3^k-1)/2 = 0 -> 1 -> 4 -> 13 -> 40 -> 121 -> ... under x -> 3x+1 (the unaccelerated odd step, the")
    print("    E-graph climb), and (3^k-1)/2 + 1/2 = 3^k/2: the base-3 repunits converge in Z_3 to -1/2, the fixed")
    print("    point of x -> 3x+1 (a 3-adic contraction).  Since T^(k+1)(2^k-1) = (3^k-1)/2, the k+1 steps of the")
    print("    longest rising run carry the k-th 2-adic approximant of the fixed point -1 (distance 2^k in Z_2)")
    print("    to the k-th 3-adic approximant of -1/2 (distance 3^k/2 in Z_3).  Every base-b repunit is the orbit")
    print("    of 0 under the affine map x -> bx + 1 with fixed point -1/(b-1): b = 4 is the sibling ladder")
    print("    S(p) = 4p+1 of the inverse tree (fixed point -1/3 = E(0)), b = 3 the E-graph climb, b = 2 the")
    print("    all-ones binary tail (fixed point -1 = T's hostile fixed point).")


# ------------------------------------------------------------------------------------------ (4f)
def run_crep(fam, K):
    out = subprocess.run([BIN, fam, str(K)], capture_output=True, text=True, check=True).stdout.splitlines()
    rows = {}
    for l in out:
        t = l.split()
        rows[int(t[1])] = t[2:]
    return rows


def zscores(d, keys, trend, fam=()):
    """standardized residuals.  trend 'lin': remove a least-squares line in k (hitting times grow linearly);
    'const': no trend.  Then standardize inside each octave [2^j, 2^(j+1)) by the mean and sd of the
    NON-family lengths of that octave, so every octave's reference values have mean 0 and sd 1."""
    ks = np.array(keys, dtype=float)
    ys = np.array([d[k] for k in keys], dtype=float)
    if trend == "lin":
        ys = ys - np.polyval(np.polyfit(ks, ys, 1), ks)
    r = {k: ys[i] for i, k in enumerate(keys)}
    fs = set(fam)
    octs = {}
    for k in keys:
        if k not in fs:
            octs.setdefault(int(math.log2(k)), []).append(r[k])
    out = {}
    for k in keys:
        pool = octs.get(int(math.log2(k)))
        if pool is None or len(pool) < 5:
            continue
        mu, sd = float(np.mean(pool)), float(np.std(pool)) or 1.0
        out[k] = (r[k] - mu) / sd
    return out


def strat_test(z, fam, ref, n_perm=20000, seed=1):
    """Stratified randomization test: each family member is replaced by a random non-family reference length
    from the same octave [2^j, 2^(j+1)); effect = mean z over the family; two-sided p."""
    rng = np.random.default_rng(seed)
    fs = set(fam)
    oct_ = lambda k: int(math.log2(k))
    pools = {}
    for k in ref:
        if k not in fs and k in z:
            pools.setdefault(oct_(k), []).append(z[k])
    fam = [k for k in fam if k in z and oct_(k) in pools]
    if not fam:
        return float("nan"), float("nan"), 0
    obs = float(np.mean([z[k] for k in fam]))
    sims = np.mean([np.array(pools[oct_(k)])[rng.integers(0, len(pools[oct_(k)]), n_perm)] for k in fam], axis=0)
    p = (np.sum(np.abs(sims - sims.mean()) >= abs(obs - sims.mean())) + 1) / (n_perm + 1)
    return obs, p, len(fam)


def local_rank(d, k, ref, w=20):
    """percentile rank of d[k] among the w nearest reference lengths on each side"""
    import bisect
    i = bisect.bisect_left(ref, k)
    win = [x for x in ref[max(0, i - w): i + w + 1] if x != k]
    return sum(1 for x in win if d[x] < d[k]) / len(win)


def part4f():
    print("\n(4f) 'The length of all-one primes': does the primality of the exponent mean anything for Collatz?")
    subprocess.run(["cc", "-O2", "-I/opt/homebrew/include", "-L/opt/homebrew/lib", "-o", BIN,
                    os.path.join(HERE, "procgen_repunit_20260924_collatz_repunits.c"), "-lgmp", "-lm"], check=True)
    t0 = time.time()
    M = run_crep("M", KM)
    W = run_crep("W", KW)
    R3 = run_crep("R3", K3)
    R10 = run_crep("R10", K10)
    print(f"    GMP runs: 2^k-1 for k <= {KM}, (2^k+1)/3 for odd k <= {KW}, (3^k-1)/2 for k <= {K3},"
          f" (10^k-1)/9 for k <= {K10}: {time.time() - t0:.0f} s")
    # identity flags
    bad = [k for k, r in M.items() if (r[6] != "1" or r[7] != "1") and r[6] != "-1"]
    assert not bad, bad[:5]
    badW = [k for k, r in W.items() if (r[6] != "1" or r[7] != "1") and r[6] != "-1"]
    assert not badW
    print(f"    identities T^k(2^k-1) = 3^k-1 and T^(k+1) = (3^k-1)/2 re-checked by direct GMP iteration for k <= 4000;"
          f" W identities for k <= 4000.")
    eq = sum(1 for r in M.values() if r[8] == "1")
    print(f"    height of 2^k-1 equals 3^k-1 (the rising run's peak is the maximum) for {eq} of {len(M)} lengths k.")
    # heights of the other families: the maximum of the forced prefix
    eqW = [k for k, r in W.items() if k >= 5 and abs(float(r[3]) - math.log2(3 * 2 ** (k - 2) + 2)) < 2e-4]
    nW = sum(1 for k in W if k >= 5)
    eq10 = [k for k, r in R10.items() if k >= 4 and abs(float(r[3]) - math.log2(3 * 5 ** k * 2 ** (k - 3) + 2)) < 2e-4]
    n10 = sum(1 for k in R10 if k >= 4)
    sW, s10 = set(eqW), set(eq10)
    lastW = max([k for k in W if k >= 5 and k not in sW] or [0])
    last10 = max([k for k in R10 if k >= 4 and k not in s10] or [0])
    print(f"    height of (2^k+1)/3 is the prefix peak T^2 = 3*2^(k-2)+2 for {len(eqW)} of {nW} odd k in [5, {KW}]"
          f" (last exception k = {lastW}); height of (10^k-1)/9 is the prefix peak T^3 = 3*5^k*2^(k-3)+2 for"
          f" {len(eq10)} of {n10} k in [4, {K10}] (last exception k = {last10}).")
    v2ok = all(int(r[9]) == (1 if k % 2 else v2(k) + 2) for k, r in M.items())
    s3ok = all(int(r[10]) == (1 if k % 2 == 0 else 2) for k, r in M.items())
    assert v2ok and s3ok
    # statistics as functions of k
    ks = np.array(sorted(M))

    def col(rows, i):
        return {k: float(r[i]) for k, r in rows.items()}
    tau = col(M, 0); a = col(M, 1); sig = col(M, 2); lh = col(M, 3); dw = col(M, 4); db = col(M, 5)
    kk = np.array(sorted(tau))
    tt = np.array([tau[k] for k in kk])
    beta = np.polyfit(kk, tt, 1)
    print(f"    total stopping time tau(2^k-1) = {beta[0]:.4f} k {beta[1]:+.1f} (least squares, k <= {KM}); the random-walk"
          f" model predicts 1 + log(3)/log(2/sqrt3)... = 1 + 7.638 = {1 + math.log(3) / math.log(2 / math.sqrt(3)) / 1:.3f}"
          " (k+1 prefix steps, then log(3^k/2)/0.1438).")
    ss = np.array([sig[k] for k in kk])
    bs = np.polyfit(kk, ss, 1)
    print(f"    stopping time sigma(2^k-1) = {bs[0]:.4f} k {bs[1]:+.1f}; model 1 + log(3/2)/0.1438 = "
          f"{1 + math.log(1.5) / math.log(2 / math.sqrt(3)):.3f} (the tail must fall by (3/2)^k/2).")

    sd_step = math.log(3) / 2
    K0 = 64
    # statistics of 2^k - 1: the stopping time, and the tail = orbit of (3^k-1)/2 = T^(k+1)(2^k-1)
    st = {
        "sigma(2^k-1)": ({k: float(r[2]) for k, r in M.items()}, "lin"),
        "tau of tail": ({k: float(r[11]) for k, r in M.items()}, "lin"),
        "height excess (bits)": ({k: 0.0 if r[8] == "1" else float(r[3]) - k * math.log2(3) for k, r in M.items()},
                                 "const"),
        "odd fraction of tail": ({k: float(r[12]) / float(r[11]) for k, r in M.items() if float(r[11]) > 0}, "const"),
        "Dword(tail)/sqrt(tau)": ({k: float(r[13]) / math.sqrt(float(r[11])) for k, r in M.items() if float(r[11]) > 0},
                                  "const"),
        "Dbridge(tail)/(sd sqrt)": ({k: float(r[14]) / (sd_step * math.sqrt(float(r[11]))) for k, r in M.items()
                                     if float(r[11]) > 0}, "const"),
    }
    assert all(abs(float(M[k][11]) - float(R3[k][0])) < 0.5 and abs(float(M[k][13]) - float(R3[k][4])) < 1e-3
               for k in R3 if k in M), "tail statistics != those of (3^k-1)/2"
    print("    tail check: the orbit statistics of T^(k+1)(2^k-1) equal those of (3^k-1)/2 computed separately (k <= %d)." % K3)
    pairs = sum(1 for k in range(3, KM, 2) if k + 1 in M and int(M[k + 1][0]) == int(M[k][0]) + 1)
    print(f"    pairing (PROVED: T(x) = T(3x+1) for odd x, and (3^(2j)-1)/2 = 3 (3^(2j-1)-1)/2 + 1): the orbits of 2^(2j-1)-1 and"
          f" 2^(2j)-1 merge, tau(2^(2j)-1) = tau(2^(2j-1)-1) + 1; checked for {pairs} of {len(range(3, KM, 2))} odd k.")
    oddprimes = [p for p in sympy.primerange(K0, KM + 1)]
    families = {"Mersenne (2^p-1 prime)": [p for p in A000043 if K0 <= p <= KM],
                "Wagstaff ((2^p+1)/3 prime)": [p for p in A000978 if K0 <= p <= KM],
                "base-3 ((3^p-1)/2 prime)": [p for p in A028491 if K0 <= p <= KM],
                "base-10 (R_p prime)": [p for p in A004023 if K0 <= p <= KM]}
    print(f"    reference: odd primes {K0} <= k <= {KM}.  z-scores: a linear trend in k is removed from the hitting times;")
    print("    every statistic is then standardized inside its octave [2^j, 2^(j+1)) by the non-family primes of that")
    print("    octave; each family member is compared with random non-family primes of the same octave (stratified")
    print("    randomization test, 20000 draws).  Statistics of the TAIL use the orbit of (3^k-1)/2, since the first")
    print("    k+1 steps of 2^k-1 are forced (identity above).")
    print("    Statistic of n = 2^k - 1         family                       n   mean z    p")
    pvals = []
    for sname, (dd, trend) in st.items():
        keys = [k for k in oddprimes if k in dd]
        allfam = sorted({k for f in families.values() for k in f})
        z = zscores(dd, keys, trend, allfam)
        for fname, fam in families.items():
            obs, p, nf = strat_test(z, [k for k in fam if k in z], keys)
            pvals.append((p, sname, fname))
            print(f"      {sname:<25} {fname:<26} {nf:>3}   {obs:+.3f}   {p:.3f}")
    ps = sorted(pvals)
    print(f"    {len(pvals)} tests; smallest p = {ps[0][0]:.3f} ({ps[0][1]}, {ps[0][2]}); Bonferroni threshold"
          f" 0.05/{len(pvals)} = {0.05 / len(pvals):.4f}: {'no test survives' if ps[0][0] > 0.05 / len(pvals) else 'SOME TEST SURVIVES'};"
          f" {sum(1 for p in ps if p[0] < 0.05)} of {len(pvals)} below 0.05 (about {0.05 * len(pvals):.1f} expected by chance).")
    # what a naive comparison would have shown
    fam = families["Mersenne (2^p-1 prime)"]
    dd = {k: float(r[4]) / math.sqrt(float(r[0])) for k, r in M.items()}
    naive = (np.mean([dd[k] for k in fam]) - np.mean([dd[k] for k in oddprimes])) / np.std([dd[k] for k in oddprimes])
    print(f"    (for the record: the WHOLE-orbit discrepancy Dword/sqrt(tau) of 2^k-1 grows like sqrt(k) because of the forced"
          f" prefix 1^k 0; a naive unstratified comparison gives Mersenne {naive:+.2f} sd -- an artefact of the family's"
          " small k, which the octave stratification removes.)")
    print("    Each family on the total stopping time of its own repunit (z: linear trend, sqrt(k) scale; stratified):")
    own = [("Wagstaff", W, A000978, KW), ("base-3", R3, A028491, K3), ("base-10", R10, A004023, K10)]
    for fname, rows, lst, KK in own:
        ref = [p for p in sympy.primerange(K0, KK + 1) if p in rows]
        dd = {k: float(r[0]) for k, r in rows.items()}
        fam = [p for p in lst if K0 <= p <= KK]
        z = zscores(dd, ref, "lin", fam)
        obs, p, nf = strat_test(z, fam, ref)
        pvals.append((p, "own tau", fname))
        print(f"      {fname:<9} exponents on tau((their repunit)): n = {nf:>2}, mean z {obs:+.3f}, p = {p:.3f}")
    print(f"    exponents below {K0} are not tested: the octaves below {K0} hold too few non-family primes, and there the")
    print("    statistics are dominated by the smallness of the numbers (a naive comparison is biased, see above).")
    # power: minimal detectable shift
    n = len([k for k in families["Mersenne (2^p-1 prime)"]])
    n10 = max(1, len([p for p in A004023 if K0 <= p <= K10]))
    print(f"    power: with n = {n} (Mersenne) a mean shift of {2.8 / math.sqrt(n):.2f} sd would be detected at alpha = 0.05"
          f" with 80% power; for the {n10} base-10 exponents in [{K0}, {K10}] only shifts > {2.8 / math.sqrt(n10):.1f} sd.")
    # 2-adic continuity in k: first m parity bits of (3^k-1)/2 depend only on k mod 2^(m-1)
    for m in range(2, 14):
        mod = 2 ** (m - 1)
        table = {}
        for k in range(2, 3000):
            wbits = word(F((3 ** k - 1) // 2), m)
            key = k % mod
            if key in table:
                assert table[key] == wbits, (m, k)
            else:
                table[key] = wbits
    print("    2-adic continuity (PROVED + checked m <= 13, k < 3000): the first m parity bits after the prefix 1^k 0 of")
    print("    2^k-1 (the bits of (3^k-1)/2) depend only on k mod 2^(m-1), since 3 has order 2^(m-1) mod 2^(m+1).")
    print("    So the orbit data of 2^k-1 beyond the forced prefix are (at every finite resolution) a function of k in Z_2,")
    print("    and primality of k -- equidistributed over the odd classes mod 2^m (Dirichlet) -- cannot bias them.")
    # the odd/even k effect is exact and trivial
    ev = [int(M[k][9]) for k in M if k % 2 == 0]
    print(f"    the one exact k-effect: after the prefix, 2^k-1 halves v_2(3^k-1)-1 = v_2(k)+1 more times when k is even"
          f" (0 when k is odd); mean over even k <= {KM}: {np.mean(ev) - 1:.3f} (-> 1 + E[v_2(k) | k even] = 3).")


# ------------------------------------------------------------------------------------------ (4g)
def part4g():
    print("\n(4g) A tested coincidence: best-approximation numerators of log_2 3")
    cf = []
    import mpmath as mp
    mp.mp.dps = 60
    x = mp.log(3, 2)
    for _ in range(25):
        a = int(mp.floor(x))
        cf.append(a)
        x = 1 / (x - a)
    nums = set()
    p0, q0, p1, q1 = 1, 0, cf[0], 1
    nums.add(p1)
    for a in cf[1:]:
        for j in range(1, a + 1):             # semiconvergents (j < a) and the convergent (j = a)
            nums.add(j * p1 + p0)
        p0, q0, p1, q1 = p1, q1, a * p1 + p0, a * q1 + q0
    print(f"    continued fraction of log_2 3: {cf[:15]} ...")
    lim = 1100
    S = sorted(n for n in nums if n <= lim)
    print(f"    numerators of convergents and semiconvergents <= {lim}: {S}")
    primes = list(sympy.primerange(2, lim + 1))
    special = [p for p in primes if p in nums]
    fam = [p for p in A004023 if p <= lim]
    hits = [p for p in fam if p in nums]
    N, K, n = len(primes), len(special), len(fam)
    pv = sum(math.comb(K, i) * math.comb(N - K, n - i) for i in range(len(hits), min(K, n) + 1)) / math.comb(N, n)
    print(f"    base-10 repunit exponents <= {lim}: {fam}; hits {hits}; {K} of the {N} primes <= {lim} are such numerators;"
          f" hypergeometric P(>= {len(hits)} hits) = {pv:.2e} (a posteriori choice of the target set!)")
    later = [p for p in A004023 if p > lim]
    big = sorted(n for n in nums if n <= 10 ** 7)
    print(f"    out-of-sample: the 6 later exponents {later} hit {[p for p in later if p in nums]} of the numerators"
          f" {[m for m in big if m > lim][:12]} ...")
    for name, lst in (("Mersenne", A000043), ("Wagstaff", A000978), ("base-3", A028491)):
        f2 = [p for p in lst if p <= lim]
        h2 = [p for p in f2 if p in nums]
        pv2 = sum(math.comb(K, i) * math.comb(N - K, len(f2) - i) for i in range(len(h2), min(K, len(f2)) + 1)) / math.comb(N, len(f2))
        print(f"    {name:<9} exponents <= {lim}: {len(f2)} with hits {h2}: P(>= {len(h2)}) = {pv2:.3f}")
    for name, lst in (("base-10", A004023), ("Mersenne", A000043), ("Wagstaff", A000978), ("base-3", A028491)):
        f3 = [p for p in lst if 20 < p <= lim]
        print(f"    {name:<9} exponents in (20, {lim}]: {f3} -> hits {[p for p in f3 if p in nums]}")
    print("    Verdict: NUMEROLOGY.  2 and 19 are small numbers in a set that is dense among small numbers; 317/200 is")
    print("    a semiconvergent; no later base-10 exponent is one, and the other three families show no excess.")


def main():
    t0 = time.time()
    print("=" * 100)
    print("PART 4. Collatz's repunits: the dictionary, the 2-adic macrocosm, and the exponent test")
    print("=" * 100)
    part4a()
    part4b()
    part4c()
    part4d()
    part4e()
    part4f()
    part4g()
    print(f"\nPart 4 done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()

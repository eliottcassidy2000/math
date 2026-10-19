#!/usr/bin/env python3
"""Independent finite controls for the cusp arithmetic transfer theorem.

No author imports, input files, network, or writes. Requires SymPy only.
Actual recurrence certificates use an independently implemented Decimal
interval arithmetic, a Machin-series enclosure of pi, and the inherited C_2
enclosure. They certify finite samples, not an asymptotic onset.
"""

from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as Q
from math import ceil, floor
import platform
import sympy as sp


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def symbolic_controls():
    e, t, p = sp.symbols("e t p", positive=True)
    c, s = sp.symbols("c s", real=True)
    phase = p / sp.sqrt(e) + c + p * sp.sqrt(e) / 3
    density = p * (1 - e / 3) / (2 * e ** sp.Rational(3, 2))
    vector = -2 * e ** sp.Rational(3, 2) / (p * (1 - e / 3))
    check(sp.simplify(sp.diff(phase, e) + density) == 0, "Calibration density")
    check(sp.simplify(sp.diff(phase, e) * vector - 1) == 0, "Parameter flow")
    check(sp.diff(density, c) == 0, "Calibration independence")
    H = t + p ** 2 / (3 * t)
    check(sp.simplify((H + c + s) - c - (H + s)) == 0, "Translation C cancellation")
    w = 2 * p ** 2 / (t ** 3 * (1 - p ** 2 / (3 * t ** 2)))
    check(sp.simplify(w * density.subs(e, p ** 2 / t ** 2) - 1) == 0,
          "Reciprocal Jacobians")
    print("SYMBOLIC: 5 exact phase, density, flow, and constant-cancellation identities")


def domain_controls():
    u0 = Q(17, 4)
    shift = lambda u, q: u + q
    reflect = lambda u, m: m - 1 - u
    in_h = lambda u, m: u0 < u < m - 1 - u0
    count = boundary = 0
    for j in range(-5, 161):
        u = Q(j, 4)
        for q in range(-9, 10):
            v = shift(u, q)
            source = u > u0 and v > u0
            image = v > u0 and v - q > u0
            check(source == image, "D_q to D_-q")
            if source:
                check(shift(v, -q) == u, "Partial translation inverse")
                for p in (-20, -1, 0, 3, 27):
                    if v + p > u0:
                        check(shift(v, p) == shift(u, p + q), "Partial composition")
                        count += 1
            elif u > u0 and v <= u0:
                boundary += 1
        for m in (10, 11, 20, 61):
            if in_h(u, m):
                v = reflect(u, m)
                check(in_h(v, m) and reflect(v, m) == u, "Reflection domain/involution")
                count += 1
                for q in (-7, -1, 0, 2, 10):
                    if in_h(u + q, m) and v - q > u0:
                        check(reflect(u + q, m) == v - q, "Conjugacy")
                        count += 1
                for ell in (10, 11, 20, 61):
                    if in_h(u, ell) and in_h(reflect(u, ell), m):
                        check(reflect(reflect(u, ell), m) == u + m - ell,
                              "Reflection-composition sign")
                        count += 1
    for m in range(2, 500):
        I = [a for a in range(m + 1) if Q(m, 3) < a <= Q(2 * m, 3)]
        check(I == list(range(m // 3 + 1, 2 * m // 3 + 1)), "Full-cell endpoint range")
        for a in I:
            check(1 < 6*a-1 < 6*m and 1 < 6*a+1 < 6*m
                  and 1 < 6*(m-a)-1 < 6*m and 1 < 6*(m-a)+1 < 6*m,
                  "Four forms below 6m")
    # A=3, T=2 gives a genuine two-positive-root hostile: t=4 and 3/4.
    check(Q(4) + Q(3, 4) == Q(3, 4) + Q(3) / Q(3, 4), "Two inverse roots")
    check(Q(3, 4) < 2 < 4, "Small root is outside the domain")
    u = u0 + 1
    check(u > u0 and u-2 <= u0 and u-2+2 == u,
          "An admissible final translation need not admit the intermediate composition")
    check(boundary > 0, "No domain-exit hostile was exercised")
    print(f"EXACT: {count} affine domain/group controls; {boundary} rejected translation exits")
    print("EXACT: every middle-third endpoint class and four-form bound for 2 <= m < 500")
    print("EXACT HOSTILE: the small inverse root leaves t>T; negative translations can leave D")


def interval_length(*intervals):
    return max(Q(0), min(b for a, b in intervals) - max(a for a, b in intervals))


def pair_integral(lo, hi, slope, offset, statistic, threshold, eta=None):
    """Intersect each pair of staircase cells directly; no breakpoint sampling."""
    lo, hi, offset = Q(lo), Q(hi), Q(offset)
    vlo, vhi = sorted((slope * lo + offset, slope * hi + offset))
    good = [(lo, hi)] if eta is None else [
        (Q(j) + eta, Q(j+1) - eta) for j in range(floor(lo), ceil(hi))]
    value = Q(0)
    for n in range(floor(lo) - 2, ceil(hi) + 3):
        first = (threshold(n-1), threshold(n))
        if not interval_length((lo, hi), first):
            continue
        for k in range(floor(vlo) - 2, ceil(vhi) + 3):
            second = tuple(sorted(((threshold(k-1) - offset) / slope,
                                   (threshold(k) - offset) / slope)))
            joint = (max(lo, first[0], second[0]), min(hi, first[1], second[1]))
            if joint[0] < joint[1]:
                length = sum(interval_length(joint, part) for part in good)
                value += length * statistic(n, k)
    return value


def integral_controls():
    rho, eta = Q(1, 1000), Q(1, 100)
    threshold = lambda n: Q(n) + (-1) ** n * rho
    statistic = lambda a, b: ((2*a + b*b) % 11) - 5
    count = 0
    for M in (20, 21, 33):
        for q in (-10, -1, 0, 1, 2, 25):
            total = sum(statistic(n, n+q) for n in range(M+1, 2*M+1))
            good = pair_integral(M, 2*M, 1, q, statistic, threshold, eta)
            full = pair_integral(M, 2*M, 1, q, statistic, threshold)
            check(good == (1-2*eta)*total, "Exact guarded translation count")
            check(abs(full-total) <= 4*5*M*eta, "Full translation error")
            count += 1
    for m in (60, 61, 62, 120, 121, 122):
        a, b = m//3, (2*m)//3
        total = sum(statistic(n, m-n) for n in range(a+1, b+1))
        good = pair_integral(a, b, -1, m-1, statistic, threshold, eta)
        full = pair_integral(a, b, -1, m-1, statistic, threshold)
        check(good == (1-2*eta)*total, "Exact guarded reflection count")
        check(abs(full-total) <= 4*5*(b-a)*eta, "Full reflection error")
        count += 1
    wrong = lambda a, b: int(b-a == 2)
    hostile = pair_integral(20, 40, 1, 1, wrong, threshold)
    check(hostile == Q(1, 50), "Constant-error positivity hostile")
    check(pair_integral(20, 40, 1, 1, wrong, threshold, eta) == 0, "Hostile removed by guard")
    shrinking = lambda n: Q(n) + Q((-1) ** n, 100*n*n)
    thin_hostile = pair_integral(20, 40, 1, 1, wrong, shrinking)
    check(thin_hostile > 0, "Shrinking-error positivity hostile")
    check(pair_integral(20, 40, 1, 1, wrong, shrinking, eta) == 0, "Shrinking hostile guarded")
    for n in range(10, 100):
        check(ceil(Q(n)) + ceil(200 - 1 - Q(n)) == 199, "Integer-phase exclusion")
        for theta in (Q(1, 1000), Q(1, 3), Q(999, 1000)):
            u = Q(n) + theta
            check(ceil(u) + ceil(200 - 1 - u) == 200, "Correct reflection offset")
            check(ceil(u) + ceil(200-u) == 201, "Wrong reflection offset")
    print(f"EXACT: {count} direct-cell-intersection integral identities and error controls")
    print(f"EXACT HOSTILE: unguarded false positive = {hostile}; guarded value and discrete sum = 0")
    print(f"EXACT HOSTILE: O(n^-2) threshold perturbations also give a false positive ({float(thin_hostile):.12g})")


PREC = 75
DOWN = Context(prec=PREC, rounding=ROUND_FLOOR)
UP = Context(prec=PREC, rounding=ROUND_CEILING)


class Interval:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = D(lo), D(lo if hi is None else hi)
        check(self.lo <= self.hi, "Invalid interval")

    @staticmethod
    def rational(q):
        q = Q(q)
        return Interval(DOWN.divide(D(q.numerator), D(q.denominator)),
                        UP.divide(D(q.numerator), D(q.denominator)))

    def __add__(self, b):
        return Interval(DOWN.add(self.lo, b.lo), UP.add(self.hi, b.hi))

    def __sub__(self, b):
        return Interval(DOWN.subtract(self.lo, b.hi), UP.subtract(self.hi, b.lo))

    def __mul__(self, b):
        lows = [DOWN.multiply(x, y) for x in (self.lo, self.hi) for y in (b.lo, b.hi)]
        highs = [UP.multiply(x, y) for x in (self.lo, self.hi) for y in (b.lo, b.hi)]
        return Interval(min(lows), max(highs))

    def __truediv__(self, b):
        check(b.lo > 0, "Division interval must be positive")
        return self * Interval(DOWN.divide(D(1), b.hi), UP.divide(D(1), b.lo))

    def sqrt(self):
        check(self.lo > 0, "Square-root interval must be positive")
        # Decimal square roots are correctly rounded to nearest. One adjacent
        # representable number on either side makes a rigorous enclosure.
        return Interval(DOWN.next_minus(DOWN.sqrt(self.lo)), UP.next_plus(UP.sqrt(self.hi)))


def atan_bounds(inv, terms):
    partial = sum((Q((-1)**k, (2*k+1)*inv**(2*k+1)) for k in range(terms)), Q(0))
    next_partial = partial + Q((-1)**terms, (2*terms+1)*inv**(2*terms+1))
    return min(partial, next_partial), max(partial, next_partial)


def pi_interval():
    a, b = atan_bounds(5, 80), atan_bounds(239, 25)
    lo, hi = 16*a[0] - 4*b[1], 16*a[1] - 4*b[0]
    check(hi-lo < Q(1, 10**110), "Machin enclosure width")
    return Interval(Interval.rational(lo).lo, Interval.rational(hi).hi)


PI = pi_interval()
C2 = Interval("-2.463836664", "-2.463836643")
PI2 = PI * PI
A = PI2 / Interval(3)
CACHE = {}


def epsilon_interval(u):
    v = Interval.rational(u) - C2
    t = (v + (v*v - Interval(4)*A).sqrt()) / Interval(2)
    return PI2 / (t*t)


def escape(u):
    u = Q(u)
    if u in CACHE:
        return CACHE[u]
    eps = epsilon_interval(u)
    x = Interval(0)
    constant = Interval("0.25") + eps
    for n in range(1, ceil(u)+50):
        x = x*x + constant
        if x.lo > 2:
            CACHE[u] = n
            return n
        check(x.hi <= 2, f"Inconclusive strict crossing for phase {u}, step {n}")
    raise RuntimeError("Finite-test cap reached")


def rational_escape(eps, radius):
    x = Q(0)
    for n in range(1, 10):
        x = x*x + Q(1, 4) + eps
        if x > radius:
            return n
    raise RuntimeError("Rational endpoint cap reached")


def dynamics_controls():
    translations = reflections = 0
    for n in (31, 75, 250, 1000, 2500):
        for theta in (Q(1, 20), Q(1, 3), Q(19, 20)):
            u = n + theta
            check(escape(u) == ceil(u), "Finite actual rounding sample")
            for q in (-(n//2), -7, 1, 2, 53):
                check(escape(u+q) == escape(u)+q, "Finite actual translation pair")
                translations += 1
    for m in (93, 94, 95, 301, 1001, 3002):
        a, b = m//3, 2*m//3
        for n in (a+1, (a+b+1)//2, b):
            for theta in (Q(1, 20), Q(1, 2), Q(19, 20)):
                u = n-1+theta
                check(escape(u) == n and escape(m-1-u) == m-n, "Finite reflected pair")
                reflections += 1
    tests = [(Q(1,4), 3), (Q(6,25), 3), (Q(13,50), 2),
             (Q(1,2), 2), (Q(49,100), 2), (Q(51,100), 1)]
    for eps, expected in tests:
        check(rational_escape(eps, Q(3,4)) == expected, "Strict rational endpoint")
    print(f"DECIMAL INTERVAL: {translations} actual translation pairs and {reflections} actual reflection pairs")
    print(f"DECIMAL INTERVAL: {len(CACHE)} distinct x-variable escape evaluations; 75-digit outward arithmetic")
    print("PI: independent alternating-series Machin certificate; inherited C_2 enclosure used unchanged")
    print("EXACT: 6 recurrence endpoint/interior tests at R=3/4; equality is not escape")


if __name__ == "__main__":
    print(f"Python {platform.python_version()}, SymPy {sp.__version__}")
    symbolic_controls()
    domain_controls()
    integral_controls()
    dynamics_controls()
    print("PASS: all independent controls; no effective guard/onset or prime lower bound is claimed")

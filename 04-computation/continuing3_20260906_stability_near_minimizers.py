"""Exact identities and finite controls for THM-4454 near-minimizer rigidity.

This standalone verifier uses SymPy only for polynomial/radical identities;
all numerical comparisons use rational outward intervals.  No producer or
inherited computation is imported.  The inherited universal inequalities
remain analytical dependencies, not deductions from this finite universe.
Run with Python 3 and SymPy, normally and with -O.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from math import isqrt
import sys
import sympy as sp

sys.stdout.reconfigure(newline="\n")
GATES = 0


def gate(ok, label):
    global GATES
    if not bool(ok):
        raise RuntimeError(label)
    GATES += 1


def identity(expr, label):
    gate(sp.simplify(sp.expand(expr)) == 0, label)


@dataclass(frozen=True)
class I:
    lo: Q
    hi: Q

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @staticmethod
    def point(x):
        return I(Q(x), Q(x))

    @staticmethod
    def lift(x):
        return x if isinstance(x, I) else I.point(x)

    def __add__(self, other):
        other = I.lift(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-I.lift(other))

    def __rsub__(self, other):
        return I.lift(other) - self

    def __mul__(self, other):
        other = I.lift(other)
        values = [x * y for x in (self.lo, self.hi)
                  for y in (other.lo, other.hi)]
        return I(min(values), max(values))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = I.lift(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval contains zero")
        return self * I(1 / other.hi, 1 / other.lo)

    def __rtruediv__(self, other):
        return I.lift(other) / self

    def __pow__(self, k):
        if k < 0:
            return 1 / (self ** (-k))
        if k == 0:
            return I.point(1)
        result = I.point(1)
        for _ in range(k):
            result = result * self
        return result

    def contains(self, x):
        return self.lo <= x <= self.hi


BITS = 256


def sqrtq(x):
    x = Q(x)
    if x < 0:
        raise ValueError("negative square root")
    k = isqrt((x.numerator << (2 * BITS)) // x.denominator)
    lo = Q(k, 1 << BITS)
    hi = Q(k + 1, 1 << BITS)
    if lo * lo == x:
        hi = lo
    gate(lo * lo <= x <= hi * hi, "outward square root")
    return I(lo, hi)


def sqrt_i(x):
    return I(sqrtq(x.lo).lo, sqrtq(x.hi).hi)


def abs_upper(x):
    return max(abs(x.lo), abs(x.hi))


def display(x, digits=10):
    # The printed decimal is only a display of a proved rational enclosure.
    gate(x.hi - x.lo < Q(1, 10 ** (digits + 2)), "display enclosure width")
    return f"{float((x.lo + x.hi) / 2):.{digits}f}"


u, v = sp.sqrt(2), sp.sqrt(3)
z, h = 1 / v, 1 / u
K3 = 4 * v / (3 * (1 + u) * (1 + v))
cstar = (13 - 8 * u) / 3
Kt = (64 - 44 * u) / 3
Kd = (28 * u - 32) / 3
Ko = u - sp.Rational(2, 3)
A, B = 2 - u, u - 1
gamma = 3 * K3 / (4 * u)
C0 = A + u * gamma
a, b, p3, p4 = sp.symbols("a b p3 p4", real=True)
C = C0 - gamma * (a + b)
E = (1 - p4) / 2
D = (5 - 8 * p3 + 3 * p4) / 6
d2 = 2 - u * (a + b)
identity(1 - C - p3 + C * p4 - sp.Rational(3, 4) * E *
         (D / E - cstar - K3 * d2), "objective identity")
identity(gamma - (2 - u) * (3 - v) / 4, "gamma")
identity(C0 - (1 + u + v - u * v) / 2, "C0")
identity(C.subs({a: z, b: z}) - (3 - v) / 2, "C triple limit")
identity(z - C.subs({a: z, b: z}) * z ** 2 - (v - 1) / 2,
         "strict secant value")
identity(1 - 2 * C.subs({a: z, b: z}) * z - (2 - v),
         "strict secant derivative")

# Independently reconstruct the two lower-bound boundary polynomials used
# to extract zero sets.  Their universal sign proofs are inherited.
Fsec = 1 - C - a ** 3 - (1 - a ** 2) * b + C * (
    a ** 4 + (1 - a ** 2) * b ** 2)
identity(Fsec.subs(b, a) - 2 * gamma * (1 - a) * (a - z) * (a - h),
         "equal-top secant factor")
P = a ** 3 - (1 + u) * a ** 2 + sp.Rational(2, 3) * (a + 1)
identity(Fsec.subs(b, z) - gamma * (1 - a) * (a - z) * P,
         "threshold-top secant factor")
for aa, bb in ((1, 0), (z, z)):
    identity(Fsec.subs({a: aa, b: bb}), "first region zero")
for aa, bb, cc in ((z, z, z), (h, h, 0)):
    identity((1 - C - (aa ** 3 + bb ** 3 + cc ** 3) +
              C * (aa ** 4 + bb ** 4 + cc ** 4)).subs({a: aa, b: bb}),
             "second region zero")

# Both degenerating quotient profiles, differentiated independently in the
# square-mass and squared-imbalance coordinates.
q, w2 = sp.symbols("q w2", nonnegative=True)
t = sp.sqrt(2 - 2 * q - w2)
top3 = t * (1 - q + w2) / 2
top4 = ((1 - q) ** 2 + (2 - 2 * q - w2) * w2) / 2
aa, bb = (t + sp.sqrt(w2)) / 2, (t - sp.sqrt(w2)) / 2
identity(aa ** 3 + bb ** 3 - top3, "two-atom third moment")
identity(aa ** 4 + bb ** 4 - top4, "two-atom fourth moment")
g = B - top3 + A * top4
identity(g.subs({q: 0, w2: 0}), "two-atom objective zero")
identity(sp.diff(g, q).subs({q: 0, w2: 0}) - (7 * u / 4 - 2),
         "dust-direction derivative")
identity(sp.diff(g, w2).subs({q: 0, w2: 0}) - (2 - 11 * u / 8),
         "imbalance-direction derivative")
identity((16 * sp.Rational(1, 3)) * (7 * u / 4 - 2) - Kd,
         "dust quotient coefficient")
identity((32 * sp.Rational(1, 3)) * (2 - 11 * u / 8) - Kt,
         "imbalance quotient coefficient")
identity(sp.diff(2 - u * t, q).subs({q: 0, w2: 0}) - 1,
         "distance dust derivative")
identity(sp.diff(2 - u * t, w2).subs({q: 0, w2: 0}) - sp.Rational(1, 2),
         "distance imbalance derivative")
identity((1 - cstar) / (2 - u) - Ko, "one-atom quotient constant")
one_num = 5 - 8 * a ** 3 + 3 * a ** 4
identity(one_num - (a - 1) * (3 * a ** 3 - 5 * a ** 2 - 5 * a - 5),
         "one-atom numerator factor")
identity(sp.limit(one_num / (1 - a ** 2), a, 1) - 6,
         "one-atom numerator limit")
identity(sp.limit(3 * (1 - a ** 4) / (1 - a ** 2), a, 1) - 6,
         "one-atom denominator limit")
identity(((5 - 8 * z + 1) / 6) / sp.Rational(1, 3) - cstar -
         K3 * (2 - 2 * u * z), "three-atom quotient value")

# Moment-band constants use only rational inequalities.  This is an explicit
# count of roots, not a mass-measure argument forgetting multiplicity.
gate(Q(4 * 49, 64 * 3) > 1, "at most three band roots")
gate(Q(2 * 81, 64 * 3) == Q(27, 32), "two-root band capacity")
gate(Q(192 * 5, 6144) == Q(5, 32), "outside mass threshold")
gate(Q(192, 49) + 192 == Q(9600, 49), "distance modulus coefficient")
identity(4 * z ** 4 - 2 * z * (3 * z ** 3) + sp.Rational(1, 3) - z ** 4,
         "three-atom moment residual zero")

# Newton identities recover moments from the locally uniform product limit.
s = sp.symbols("s")
target = sp.series((1 + z * s) ** 3 * sp.exp((1 - v) * s), s, 0, 5).removeO()
eco = [sp.expand(target).coeff(s, k) for k in range(5)]
identity(eco[1] - 1, "target p1")
identity(1 - 2 * eco[2] - 1, "target p2")
identity(1 + 3 * eco[3] - z, "target p3")
identity((4 * z - 1 - 12 * eco[4]) / 3 - sp.Rational(1, 3), "target p4")

# Exact mixed-dust normalization in a symbolic quadratic extension.
L, cp, dd = sp.symbols("L cp dd", positive=True)
T = 3 + cp
Delta = 3 * L ** 2 + L * (T ** 2 - 3 + cp ** 2) - cp ** 2
rel = sp.expand(L * (T - dd) ** 2 - 3 * L - cp ** 2 - dd ** 2)
identity(Delta - T ** 2 - (L - 1) * (3 * L + T ** 2 + cp ** 2),
         "positive normalization radical")
identity(sp.discriminant(rel, dd) / 4 - Delta, "quadratic discriminant")
identity(L * (T ** 2 - 3) - cp ** 2 - ((L - 1) * cp ** 2 + 6 * L * cp + 6 * L),
         "positive negative-dust total")
for cvalue in range(4):
    for lvalue in (2, 3, 5):
        tv = sp.Rational(3 + cvalue)
        rv = sp.expand(rel.subs({L: lvalue, cp: cvalue}))
        raw = [sp.Integer(1)] * 3 + [sp.Rational(cvalue, lvalue)] * lvalue + [
            -dd / lvalue] * lvalue
        sums = [sum(x ** k for x in raw) for k in range(1, 5)]
        identity(sums[0] - (tv - dd), "literal p1 normalization")
        gate(sp.rem(sp.expand(sums[1] - (tv - dd) ** 2), rv, dd) == 0,
             "literal p2 normalization")
        # Multiply every original factor, to degree four, then square it.
        co = [sp.Integer(1)] + [sp.Integer(0)] * 4
        for x in raw:
            for k in range(4, 0, -1):
                co[k] = sp.expand(co[k] + x * co[k - 1])
        coeff = sum(co[k] * co[4 - k] for k in range(5))
        discrepancy = sp.expand(-6 * coeff - (5 * (tv - dd) ** 4 -
                                8 * sums[2] * (tv - dd) + 3 * sums[3]))
        gate(sp.rem(discrepancy, rv, dd) == 0,
             "original doubled-product quartic formula")
        # Normalized moment residual from the literal list as an identity.
        literal_M_scaled = sum(x ** 2 * (x - z * (tv - dd)) ** 2 for x in raw)
        expr = sp.expand(literal_M_scaled - (sums[3] -
            2 * z * (tv - dd) * sums[2] + (tv - dd) ** 4 / 3))
        gate(sp.simplify(sp.rem(expr, rv, dd)) == 0, "literal moment residual")

ui, vi = sqrtq(2), sqrtq(3)
zi = 1 / vi
ki = 4 * vi / (3 * (1 + ui) * (1 + vi))
ct = (13 - 8 * ui) / 3
kti, kdi, koi = (64 - 44 * ui) / 3, (28 * ui - 32) / 3, ui - Q(2, 3)
gate(kdi.lo > kti.hi > Q(1, 2) > ki.hi, "all boundary gaps")
gate(koi.lo > ki.hi, "one-atom boundary gap")
print("Boundary constants (rational enclosures, displayed):")
print("K3", display(ki), "Ktwo", display(kti), "Kdust", display(kdi),
      "Kone", display(koi))


def mixed(lvalue, cvalue):
    cv, lv = Q(cvalue), Q(lvalue)
    tv = 3 + cv
    delta = 3 * lv ** 2 + lv * (tv ** 2 - 3 + cv ** 2) - cv ** 2
    # Rationalized S avoids large-L subtraction of near-equal terms.
    si = (3 * lv + tv ** 2 + cv ** 2) / (sqrtq(delta) + tv)
    di = tv - si
    gate(si.lo > 0 and di.lo > 0, "positive mixed normalization")
    main, plus, minus = 1 / si, cv / (lv * si), -di / (lv * si)
    gate(main.lo > plus.hi and main.lo > -minus.lo, "three largest roots")
    moments = [3 * main ** k + lv * plus ** k + lv * minus ** k
               for k in range(1, 5)]
    gate(moments[0].contains(1) and moments[1].contains(1), "enclosed normalization")
    energy = (1 - moments[3]) / 2
    distance2 = 2 - 2 * ui * main
    distance3 = 2 - 6 * zi * main
    residual = moments[3] - 2 * zi * moments[2] + Q(1, 3)
    quotient = ((5 - 8 * moments[2] + 3 * moments[3]) / (6 * energy) - ct) / distance2
    tail = lv * (plus ** 2 + minus ** 2)
    gate(energy.lo > 0 and distance2.lo > 0 and tail.lo > 0, "eligible mixed family")
    gate(quotient.lo > ki.hi, "strict finite quotient control")
    gate(distance3.hi < ((1 + zi) ** 2 * 3).hi, "distance sanity")
    gate(residual.hi <= ((1 + zi) ** 2 * distance3).lo, "global moment-distance bound")
    if residual.hi < Q(5, 6144):
        gate(distance3.hi <= (Q(9600, 49) * residual).lo, "local integer-band bound")
    return quotient, distance3, tail, cv / si, di / si


print("Mixed-dust exact controls: L, c, quotient, square dust, positive/negative first mass")
for cv in (0, 1, 3):
    for lv in (8, 16, 64, 256, 4096, 65536):
        row = mixed(lv, cv)
        if lv == 65536:
            print(lv, cv, *(display(row[i]) for i in (0, 2, 3, 4)))
            gate(abs_upper(row[0] - ki) < Q(1, 50), "bounded-dust approach")
for nn in (2, 4, 8, 16, 32, 64):
    row = mixed(nn ** 4, nn)
    if nn == 64:
        print(nn ** 4, nn, *(display(row[i]) for i in (0, 2, 3, 4)))
        gate(row[3].lo > 30 and row[4].lo > 30, "unbounded-dust hostile control")
        gate(row[2].hi < Q(1, 1000), "unbounded-dust small square mass")
        gate(abs_upper(row[0] - ki) < Q(1, 1000), "unbounded-dust near minimizer")

# A literal rational three-root family approaches the one-atom boundary.
for nn in (10, 20, 40, 80, 160, 320, 640, 4096):
    small = Q(1, nn)
    raw = [Q(1), small, -small / (1 + small)]
    ss = sum(raw)
    roots = [x / ss for x in raw]
    pp = [sum(x ** k for x in roots) for k in range(1, 5)]
    gate(pp[0] == pp[1] == 1, "one-atom exact normalization")
    ee = (1 - pp[3]) / 2
    jj = (5 - 8 * pp[2] + 3 * pp[3]) / (6 * ee)
    gate(jj == 1, "one-atom rational family exact J")
    rr = (jj - ct) / (2 - ui * (roots[0] + roots[1]))
    gate(rr.lo > ki.hi, "one-atom hostile quotient")
    if nn == 4096:
        gate(abs_upper(rr - koi) < Q(1, 1000), "one-atom limit control")
        print("One-atom n=4096 quotient", display(rr))


def two_family(lvalue, eta):
    lv, eta = Q(lvalue), Q(eta)
    rawsquare = 2 + 2 * eta ** 2
    radical = lv ** 2 * rawsquare + lv * (4 - rawsquare)
    ss = (lv * rawsquare + 4) / (sqrtq(radical) + 2)
    dd = 2 - ss
    aa, bb, dust = (1 + eta) / ss, (1 - eta) / ss, -dd / (lv * ss)
    gate(bb.lo > -dust.lo > 0, "two-atom largest positives")
    pp = [aa ** k + bb ** k + lv * dust ** k for k in range(1, 5)]
    gate(pp[0].contains(1) and pp[1].contains(1), "two-atom normalization")
    ee = (1 - pp[3]) / 2
    rr = ((5 - 8 * pp[2] + 3 * pp[3]) / (6 * ee) - ct) / (2 - ui * (aa + bb))
    tail = lv * dust ** 2
    imbalance = (aa - bb) ** 2
    mixing = tail / (tail + imbalance / 2)
    predicted = kti + (kdi - kti) * mixing
    gate(rr.lo > ki.hi, "two-atom hostile quotient")
    gate(abs_upper(rr - predicted) < 100 * sqrt_i(tail + imbalance).lo,
         "uniform two-scale profile control")
    return rr, mixing


print("Two-atom boundary controls: n, regime, quotient, mixing fraction")
for nn in (8, 16, 32, 64, 128, 256):
    for label, lv, eta in (("dust", nn * nn, Q(0)),
                           ("imbalance", nn ** 4, Q(1, nn)),
                           ("mixed", nn * nn, Q(1, nn))):
        rr, mix = two_family(lv, eta)
        if nn == 256:
            print(nn, label, display(rr), display(mix))
            if label == "dust":
                gate(abs_upper(rr - kdi) < Q(1, 1000), "dust endpoint approach")
            if label == "imbalance":
                gate(abs_upper(rr - kti) < Q(1, 1000), "imbalance endpoint approach")

print("PASS", GATES, "always-active exact gates; no census or optimum inference.")
print("Universal conclusions depend on the inherited regional sign proofs plus the new analytic limit arguments.")

#!/usr/bin/env python3
"""Independent coordinates over Q(sqrt(5)); does not import producer code."""
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
import argparse
import json

CHECKS = 0


def need(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise ValueError(message)


@dataclass(frozen=True)
class R:
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __post_init__(self):
        need(type(self.a) in (int, Fraction) and type(self.b) in (int, Fraction),
             "exact rational coefficients required")
        object.__setattr__(self, "a", Fraction(self.a))
        object.__setattr__(self, "b", Fraction(self.b))

    def __add__(self, y):
        y = field(y)
        return R(self.a + y.a, self.b + y.b)

    __radd__ = __add__

    def __neg__(self):
        return R(-self.a, -self.b)

    def __sub__(self, y):
        return self + -field(y)

    def __rsub__(self, y):
        return field(y) + -self

    def __mul__(self, y):
        y = field(y)
        return R(self.a * y.a + 5 * self.b * y.b,
                 self.a * y.b + self.b * y.a)

    __rmul__ = __mul__

    def __truediv__(self, y):
        y = field(y)
        norm = y.a * y.a - 5 * y.b * y.b
        need(norm != 0, "zero divisor")
        z = self * R(y.a, -y.b)
        return R(z.a / norm, z.b / norm)

    def sign(self):
        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        if self.a == 0:
            return (self.b > 0) - (self.b < 0)
        sa = (self.a > 0) - (self.a < 0)
        sb = (self.b > 0) - (self.b < 0)
        if sa == sb:
            return sa
        norm = self.a * self.a - 5 * self.b * self.b
        need(norm != 0, "irrationality violated")
        return sa if norm > 0 else sb

    def serial(self):
        return [str(self.a), str(self.b)]


def field(x):
    return x if isinstance(x, R) else R(x)


ZERO, ONE = R(0), R(1)
PHI = R(Fraction(1, 2), Fraction(1, 2))
Q = 2 - PHI
Z, U, V, W = (ZERO, ZERO), (ONE, ZERO), (ZERO, ONE), (ONE, ONE)


def plus(a, b):
    return (a[0] + b[0], a[1] + b[1])


def minus(a, b):
    return (a[0] - b[0], a[1] - b[1])


def determinant(a, b):
    return a[0] * b[1] - a[1] * b[0]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + (PHI - 1) / 2 * \
        (a[0] * b[1] + a[1] * b[0])


def area2(poly):
    return sum((determinant(a, b) for a, b in zip(poly, poly[1:] + poly[:1])), ZERO)


def rotate(p):
    return ((1 - PHI) * p[0] - p[1], p[0])


def signs(poly):
    return [determinant(minus(poly[i], poly[i - 1]),
                        minus(poly[(i + 1) % 4], poly[i])).sign() for i in range(4)]


def prove_shape(poly, angles, reflex, lengths):
    orient = area2(poly).sign()
    need(orient != 0, "zero area")
    turns = signs(poly)
    for i in range(4):
        prev = minus(poly[i - 1], poly[i])
        nxt = minus(poly[(i + 1) % 4], poly[i])
        need(dot(prev, prev) == lengths[i - 1] * lengths[i - 1], "previous length")
        need(dot(nxt, nxt) == lengths[i] * lengths[i], "next length")
        need(dot(prev, nxt) / (lengths[i - 1] * lengths[i]) == angles[i], "corner cosine")
        need((turns[i] != orient) == (i in reflex), "convex/reflex flag")
    # Non-adjacent edges cannot properly cross: their endpoint-side tests
    # have a same-side pair. All involved orientations are nonzero here.
    for i, j in ((0, 2), (1, 3)):
        a, b = poly[i], poly[(i + 1) % 4]
        c, d = poly[j], poly[(j + 1) % 4]
        orient_values = [determinant(minus(b, a), minus(c, a)).sign(),
                         determinant(minus(b, a), minus(d, a)).sign(),
                         determinant(minus(d, c), minus(a, c)).sign(),
                         determinant(minus(d, c), minus(b, c)).sign()]
        need(all(orient_values), "degenerate shape")
        need(not (orient_values[0] != orient_values[1] and orient_values[2] != orient_values[3]),
             "self-intersection")


def coordinates():
    C = (Q, Q)
    D, K = [Z, U, C, V], [W, U, C, V]
    dc, kc = [0, 1, 0, 1], [1, 0, 1, 0]
    lengths = [ONE, PHI - 1, PHI - 1, ONE]
    cos36, cos72 = PHI / 2, (PHI - 1) / 2
    prove_shape(D, [cos72, cos36, -cos36, cos36], {2}, lengths)
    prove_shape(K, [cos72, cos72, -cos36, cos72], set(), lengths)
    need(area2(D) == 2 * Q, "dart oriented area coefficient")
    need(area2(K) == -2 * (PHI - 1), "kite oriented area coefficient")
    need(-area2(K) / area2(D) == PHI, "area ratio")
    need(determinant(rotate(U), rotate(V)) == ONE, "not a direct map")
    for a, b in product((U, V), repeat=2):
        need(dot(rotate(a), rotate(b)) == dot(a, b), "Gram matrix changed")
    good = [minus(p, V) for p in K]
    bad = [plus(rotate(p), (PHI, -ONE)) for p in K]
    expected_good = [U, (ONE, -ONE), (Q, 1 - PHI), Z]
    expected_bad = [Z, U, (ONE, 1 - PHI), (PHI - 1, -ONE)]
    need(good == expected_good and bad == expected_bad, "coordinate reduction mismatch")
    need(all(p[1].sign() >= 0 for p in D), "dart halfplane")
    for P in (D, good, bad):
        need(all(p[0].sign() >= 0 and (1 - p[0]).sign() >= 0 for p in P), "strip bound")
    for P in (good, bad):
        need(all(p[1].sign() <= 0 and (p[1] + 1).sign() >= 0 for p in P), "kite lower band")
        need({p for p in P if p[1] == ZERO} == {Z, U}, "shared seam wrong")
        need(area2(P) == area2(K), "direct area changed")
    colors = dict(zip(D, dc))
    endpoints = []
    for P in (bad, good):
        marked = dict(zip(P, kc))
        endpoints.append([marked[Z] == colors[Z], marked[U] == colors[U]])
    need(endpoints == [[False, False], [True, True]], "wrong endpoint colors")
    return D, (bad, good), endpoints


def guards(D, choices, k):
    count = 0
    state_count = 1 << k
    moments = [[0] * state_count for _ in range(2)]
    for state in range(state_count):
        for subset in range(state_count):
            all_match = True
            for j in range(k):
                bit = (state >> j) & 1
                selected = (subset >> j) & 1
                # A named move translates the initially staged tile by +10V.
                offset = (R(3 * j), R(0 if selected else -10))
                final = [plus(p, offset) for p in choices[bit]]
                need(all((p[0] - 3 * j).sign() >= 0 and (3 * j + 1 - p[0]).sign() >= 0
                         for p in final), "gadget column collision")
                upper = 0 if selected else -10
                need(all((upper - p[1]).sign() >= 0 and (p[1] - upper + 1).sign() >= 0
                         for p in final), "staging band collision")
                if selected:
                    marks = dict(zip(final, (1, 0, 1, 0)))
                    all_match &= marks[(R(3 * j), ZERO)] == 0
                    all_match &= marks[(R(3 * j + 1), ZERO)] == 1
            expected = (state & subset) == subset
            need(all_match == expected, "native geometric guard != conjunction")
            moments[state.bit_count() % 2][subset] += int(all_match)
            count += 1
    for subset in range(state_count - 1):
        need(moments[0][subset] == moments[1][subset], "proper marginal distinguishes parity")
    need(abs(moments[0][-1] - moments[1][-1]) == 1, "top conjunction fails to distinguish")
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutate", choices=("color", "rotation", "stage"))
    args = parser.parse_args()
    D, choices, endpoints = coordinates()
    census = [{"k": k, "state_subset_placements": guards(D, choices, k)} for k in range(1, 6)]
    if args.mutate == "color":
        need(endpoints[0] == endpoints[1], "CORRUPTION: endpoint marking erased")
    if args.mutate == "rotation":
        need(determinant(rotate(U), rotate(V)) == -ONE, "CORRUPTION: direct map called reflection")
    if args.mutate == "stage":
        need(R(-10).sign() >= 0, "CORRUPTION: staged band overlaps dart band")
    print(json.dumps({"status": "PASS", "field": "independent Q(sqrt(5)) implementation",
                      "active_checks": CHECKS, "placement_census": census,
                      "total_state_subset_placements": sum(x["state_subset_placements"] for x in census),
                      "dart_angles_degrees": [72, 36, 216, 36],
                      "kite_angles_degrees": [72, 72, 144, 72],
                      "bad_good_endpoint_matches": endpoints,
                      "all_k_proof": "disjoint oblique columns [3j,3j+1], dart band [0,1], staged band [-11,-10], attached band [-1,0]",
                      "scope": "finite marked local-matching guards; no infinite tiling extendibility assertion"},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

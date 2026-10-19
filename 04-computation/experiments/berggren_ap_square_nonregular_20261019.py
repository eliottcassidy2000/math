#!/usr/bin/env python3
"""Exact finite controls for the quartic-wall nonregularity proof.

Python standard library only. Every mathematical check remains active under -O.
The infinite-language conclusion is proved in berggren_ap_square_nonregular_20261019.md, not inferred
from this finite computation.
"""
from fractions import Fraction as Q
from itertools import product
from math import gcd, isqrt

ALPHABET = "LMR"
IDENTITY = (1, 0, 0, 1)
INVERSE = {"L": (0, 1, -1, 2), "M": (0, 1, 1, -2), "R": (1, -2, 0, 1)}


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def polynomial(t):
    return t**4 - 8*t**3 + 2*t**2 + 8*t + 1


def branch(letter, t):
    if letter == "L":
        return 2 - 1/t
    if letter == "M":
        return 2 + 1/t
    if letter == "R":
        return t + 2
    raise ValueError(letter)


def root_eval(word):
    t = Q(2)
    for letter in word:
        t = branch(letter, t)
    return t


def outer_eval(word):
    t = Q(2)
    for letter in reversed(word):
        t = branch(letter, t)
    return t


def words_through(depth):
    for n in range(depth + 1):
        for letters in product(ALPHABET, repeat=n):
            yield "".join(letters)


def multiply(A, B):
    a, b, c, d = A
    e, f, g, h = B
    return a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h


def canonical(A):
    first = next(x for x in A if x)
    return A if first > 0 else tuple(-x for x in A)


class Alpha:
    """A shared exact isolating interval for the root alpha in (1,2).

    P'(t)<0 on (1,2), and P(1)>0>P(2). Rational bisection therefore
    preserves this unique root and provides certified comparisons.
    """
    def __init__(self):
        self.lo, self.hi = Q(1), Q(2)
        self.bisections = 0

    def refine(self):
        mid = (self.lo+self.hi)/2
        value = polynomial(mid)
        check(value != 0, "rational quartic root unexpectedly encountered")
        if value > 0:
            self.lo = mid
        else:
            self.hi = mid
        self.bisections += 1

    def bounds(self, A):
        a, b, c, d = A
        while (c*self.lo+d)*(c*self.hi+d) <= 0:
            self.refine()
        u = (a*self.lo+b)/(c*self.lo+d)
        v = (a*self.hi+b)/(c*self.hi+d)
        return min(u, v), max(u, v)

    def compare(self, A, rational):
        """Return sign(A(alpha)-rational), certified by exact enclosure."""
        for _ in range(10000):
            lo, hi = self.bounds(A)
            if lo > rational:
                return 1
            if hi < rational:
                return -1
            self.refine()
        raise RuntimeError("comparison did not isolate within safety bound")

    def itinerary_step(self, A):
        check(self.compare(A, Q(1)) > 0, "inverse orbit escaped J")
        if self.compare(A, Q(2)) < 0:
            letter = "L"
        elif self.compare(A, Q(3)) < 0:
            letter = "M"
        else:
            letter = "R"
        return letter, multiply(INVERSE[letter], A)


def descent(t):
    """Return the unique root-to-leaf word by exact numerator descent."""
    backwards = []
    while t != 2:
        old_height = t.numerator+t.denominator
        if 1 < t < 2:
            letter, parent = "L", 1/(2-t)
        elif 2 < t < 3:
            letter, parent = "M", 1/(t-2)
        elif t > 3:
            letter, parent = "R", t-2
        else:
            raise RuntimeError("invalid opposite-parity descent seam")
        check(parent > 1, "descent outside J")
        check((parent.numerator+parent.denominator) % 2 == 1, "parity lost")
        check(parent.numerator+parent.denominator < old_height, "height did not decrease")
        backwards.append(letter)
        t = parent
    return "".join(reversed(backwards))


def quadratic_lower(t):
    # On J, t < 1+sqrt(2) iff t^2-2t-1 < 0.
    return t*t-2*t-1 < 0


def quadratic_dfa(word, wrong_middle_sign=False):
    # States -1,+1 mean lower/upper cut; 0,2 are rejecting/accepting sinks.
    state = -1
    for letter in word:
        if state in (0, 2):
            continue
        if letter == "M":
            if not wrong_middle_sign:
                state = -state
        elif letter == "L":
            state = 2 if state == -1 else 0
        else:
            state = 0 if state == -1 else 2
    return state in (-1, 2)  # Root 2 lies below the quadratic wall.


def main():
    alpha = Alpha()
    check(polynomial(Q(1)) == polynomial(Q(-1)) == 4, "rational-root control")
    check(polynomial(Q(2)) < 0, "alpha interval")
    check(polynomial(Q(3)) < 0 and polynomial(Q(8)) > 0, "upper-wall interval")
    check(isqrt(48)**2 != 48, "quadratic-factor obstruction")

    height_count = 0
    max_descent = 0
    for n in range(2, 241):
        for m in range(1, n):
            if gcd(n, m) != 1 or (n+m) % 2 != 1:
                continue
            t = Q(n, m)
            word = descent(t)
            check(root_eval(word) == t, "descent reconstruction")
            height_count += 1
            max_descent = max(max_descent, len(word))

    counts = []
    total = 0
    quadratic_checks = 0
    wrong_sign_witness = None
    for depth in range(9):
        accepted = 0
        for letters in product(ALPHABET, repeat=depth):
            word = "".join(letters)
            t = root_eval(word)
            check(t == outer_eval(word[::-1]), "word reversal mismatch")
            check((t.numerator+t.denominator) % 2 == 1, "forward parity failure")
            n, m = t.numerator, t.denominator
            a, b, C = n*n-m*m, 2*n*m, n*n+m*m
            check((C*C > 4*a*b) == (polynomial(t) > 0), "cone polynomial mismatch")
            check(polynomial(t) != 0, "rational wall equality")
            accepted += polynomial(t) > 0
            expected = quadratic_lower(outer_eval(word))
            check(quadratic_dfa(word) == expected, "quadratic automaton failed")
            if wrong_sign_witness is None and quadratic_dfa(word, True) != expected:
                wrong_sign_witness = word
            quadratic_checks += 1
            total += 1
        counts.append(accepted)

    orbit = [IDENTITY]
    itinerary = ""
    for _ in range(128):
        letter, A = alpha.itinerary_step(orbit[-1])
        itinerary += letter
        orbit.append(A)
    check(len({canonical(A) for A in orbit}) == len(orbit), "repeated finite inverse state")
    for A in orbit:
        a, b, c, d = A
        check(abs(a*d-b*c) == 1, "inverse matrix determinant")

    residual_checks = 0
    for depth in range(17):
        prefix = itinerary[:depth]
        lower = prefix.count("M") % 2 == 0
        for suffix in words_through(4):
            actual = alpha.compare(IDENTITY, outer_eval(prefix+suffix)) > 0
            comparison = alpha.compare(orbit[depth], outer_eval(suffix))
            expected = (comparison > 0) if lower else (comparison < 0)
            check(actual == expected, "exact residual recurrence failed")
            residual_checks += 1

    isolate_checks = 0
    _, lower_residual = alpha.itinerary_step(IDENTITY)
    for suffix in words_through(6):
        image = outer_eval("L"+suffix)
        check(1 < image < 2, "lower branch range failed")
        actual = polynomial(image) > 0
        expected = alpha.compare(lower_residual, outer_eval(suffix)) > 0
        check(actual == expected, "two-wall isolation quotient failed")
        isolate_checks += 1

    check(wrong_sign_witness == "M", "wrong-orientation hostile did not trigger at M")
    print("STATUS: PASS (FINITE-EXACT controls; infinite nonregularity is proved in note)")
    print(f"Opposite-parity slopes with numerator <=240: {height_count}; max descent length: {max_descent}")
    print(f"Root/outer reversal and cone checks through depth 8: {total}")
    print(f"Accepted cone counts at depths 0..8: {counts}")
    print(f"Quadratic positive-control DFA checks: {quadratic_checks}")
    print(f"Wrong-orientation hostile minimal witness: {wrong_sign_witness}")
    print(f"Distinct exact inverse-wall matrices: {len(orbit)}")
    print(f"Quartic inverse itinerary, first 128 symbols: {itinerary}")
    print(f"Residual checks (17 prefixes x 121 suffixes): {residual_checks}")
    print(f"Two-wall isolation quotient checks through suffix depth 6: {isolate_checks}")
    print(f"Exact alpha bisections used: {alpha.bisections}")
    print("No floating arithmetic, assertions, external modules, or empirical periodicity inference.")


if __name__ == "__main__":
    main()

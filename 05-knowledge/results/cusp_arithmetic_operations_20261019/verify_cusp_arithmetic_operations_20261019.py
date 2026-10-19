#!/usr/bin/env python3
"""Exact, numerical, and directed-interval controls for cusp arithmetic.

Requires Python 3 and mpmath. No input files, network, or filesystem writes.
Finite controls are not a numerical calibration of an asymptotic guard.
The C_2 enclosure is inherited from the initial cusp-phase certificate.
"""

from fractions import Fraction as Q
from math import ceil, floor
import mpmath as mp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dist_integer(u):
    return min(u - floor(u), ceil(u) - u)


def exact_rounding_controls():
    count = 0
    for n in range(10, 81):
        for theta in (Q(1, 97), Q(1, 4), Q(1, 2), Q(96, 97)):
            u = n - 1 + theta
            for q in (-7, -1, 0, 1, 2, 19):
                require(ceil(u + q) == n + q, "Integer translation rounding")
                require(dist_integer(u + q) == dist_integer(u), "Translation guard")
                count += 1
            for m in (200, 201, 202):
                v = m - 1 - u
                require(ceil(u) + ceil(v) == m, "Complement rounding")
                require(dist_integer(v) == dist_integer(u), "Reflection guard")
                require(ceil(u) + ceil(m - u) == m + 1, "Wrong-offset hostile")
                count += 1
    for m in range(20, 81):
        for k in range(1, m):
            require(ceil(Q(k)) + ceil(m - 1 - Q(k)) == m - 1,
                    "Integer-phase endpoint hostile")
    print(f"EXACT: {count} rational rounding/guard controls, plus integer endpoint hostiles")


RHO = Q(1, 1000)
ETA = Q(1, 100)


def model_threshold(n):
    return Q(n) + (RHO if n % 2 == 0 else -RHO)


def model_escape(u):
    n = ceil(u)
    while u < model_threshold(n - 1):
        n -= 1
    while u >= model_threshold(n):
        n += 1
    return n


def exact_model_integral(lo, hi, slope, offset, statistic, guarded):
    """Integrate a finite rational step function, splitting every breakpoint."""
    lo, hi, offset = Q(lo), Q(hi), Q(offset)
    endpoints = {lo, hi}
    for n in range(floor(lo) - 3, ceil(hi) + 4):
        endpoints.add(model_threshold(n))
        if guarded:
            endpoints.add(Q(n) + ETA)
            endpoints.add(Q(n) - ETA)
    vlo, vhi = sorted((slope * lo + offset, slope * hi + offset))
    for n in range(floor(vlo) - 3, ceil(vhi) + 4):
        endpoints.add((model_threshold(n) - offset) / slope)
    endpoints = sorted(x for x in endpoints if lo <= x <= hi)
    value = Q(0)
    for left, right in zip(endpoints, endpoints[1:]):
        u = (left + right) / 2
        if guarded and dist_integer(u) <= ETA:
            continue
        value += (right - left) * statistic(
            model_escape(u), model_escape(slope * u + offset))
    return value


def exact_integral_controls():
    statistic = lambda a, b: (a * a + 2 * a * b + 3 * b) % 7 - 3
    norm = 3
    controls = 0
    for M in (20, 21, 50):
        for q in (-9, -1, 0, 1, 2, 7, 100):
            total = sum(statistic(n, n + q) for n in range(M + 1, 2 * M + 1))
            guarded = exact_model_integral(M, 2 * M, 1, q, statistic, True)
            unguarded = exact_model_integral(M, 2 * M, 1, q, statistic, False)
            require(guarded == (1 - 2 * ETA) * total, "Guarded translation integral")
            require(abs(unguarded - total) <= 4 * norm * M * ETA,
                    "Unguarded translation error bound")
            controls += 1
    for m in (120, 121, 122, 123, 250):
        lo, hi = m // 3, (2 * m) // 3
        total = sum(statistic(a, m - a) for a in range(lo + 1, hi + 1))
        guarded = exact_model_integral(lo, hi, -1, m - 1, statistic, True)
        unguarded = exact_model_integral(lo, hi, -1, m - 1, statistic, False)
        require(guarded == (1 - 2 * ETA) * total, "Guarded reflection integral")
        require(abs(unguarded - total) <= 4 * norm * (hi - lo) * ETA,
                "Unguarded reflection error bound")
        controls += 1
    wrong_difference = lambda a, b: int(b - a == 2)
    false_positive = exact_model_integral(20, 40, 1, 1, wrong_difference, False)
    guarded = exact_model_integral(20, 40, 1, 1, wrong_difference, True)
    require(false_positive > 0 and guarded == 0, "Unguarded positivity hostile")
    print(f"EXACT: {controls} piecewise-rational guarded sum identities and unguarded bounds")
    print(f"EXACT HOSTILE: wrong-difference unguarded integral = {false_positive}; "
          "arithmetic sum and guarded integral = 0")


mp.mp.dps = 80
C_MID = mp.mpf("-2.4638366535")
A = mp.pi ** 2 / 3


def tau(t):
    return t + C_MID + A / t


def large_inverse(u):
    v = u - C_MID
    return (v + mp.sqrt(v * v - 4 * A)) / 2


def xi(u):
    return mp.pi ** 2 / large_inverse(u) ** 2


def phase(e):
    return tau(mp.pi / mp.sqrt(e))


def translate(e, q):
    return xi(phase(e) + q)


def universal_translate(e, s):
    t = mp.pi / mp.sqrt(e)
    h = t + A / t + s
    t_new = (h + mp.sqrt(h * h - 4 * A)) / 2
    return mp.pi ** 2 / t_new ** 2


def reflect(e, m):
    return xi(m - 1 - phase(e))


def near(a, b, label):
    require(abs(a - b) < mp.mpf("1e-65") * max(1, abs(a), abs(b)), label)


def analytic_controls():
    count = 0
    for u in map(mp.mpf, ("40.125", "100.25", "333.75", "1000.5")):
        e = xi(u)
        near(tau(large_inverse(u)), u, "Large inverse identity")
        near(phase(e), u, "Xi inverse identity")
        t = large_inverse(u)
        w = 2 * mp.pi ** 2 / (t ** 3 * (1 - A / t ** 2))
        near(-mp.diff(xi, u), w, "Parameter density identity")
        density = mp.pi / (2 * e ** mp.mpf("1.5")) * (1 - e / 3)
        near(-mp.diff(phase, e), density, "Calibration density identity")
        near(w * density, 1, "Change-of-variable reciprocity")
        near(mp.diff(lambda s: universal_translate(e, s), 0),
             -2 * e ** mp.mpf("1.5") / (mp.pi * (1 - e / 3)),
             "Universal parameter-flow vector field")
        for p, q in ((1, 2), (-7, 9), (31, -20), (0, 0)):
            near(translate(e, q), universal_translate(e, q),
                 "Endpoint phase cancellation")
            near(translate(translate(e, q), p), translate(e, p + q),
                 "Partial action composition")
            near(translate(translate(e, q), -q), e, "Partial action inverse")
            count += 3
        m = 3000
        near(reflect(reflect(e, m), m), e, "Reflection involution")
        near(reflect(translate(e, 7), m), translate(reflect(e, m), -7),
             "Reflection translation conjugacy")
        near(reflect(reflect(e, 2990), 3000), translate(e, 10),
             "Two-reflection composition")
        count += 9
    print(f"NUMERICAL: {count} analytic map/density controls at 80 digits, tolerance 1e-65")
    print("NUMERICAL scope: algebra is checked at a phase midpoint; this does not certify C_2")


mp.iv.dps = 60
C_IV = mp.iv.mpf(["-2.463836664", "-2.463836643"])


def iv_rational(q):
    q = Q(q)
    return mp.iv.mpf(q.numerator) / q.denominator


def interval_xi(u):
    v = iv_rational(u) - C_IV
    t = (v + mp.iv.sqrt(v * v - 4 * mp.iv.pi ** 2 / 3)) / 2
    return mp.iv.pi ** 2 / (t * t)


ESCAPE_CACHE = {}


def interval_escape(u):
    """Certify the first strict crossing for every C in the inherited interval."""
    u = Q(u)
    if u in ESCAPE_CACHE:
        return ESCAPE_CACHE[u]
    epsilon = interval_xi(u)
    y = mp.iv.mpf("-0.5")
    b = mp.iv.mpf("1.5")
    # N is near ceil(u); the loose cap is only a finite-test safety stop.
    for k in range(1, ceil(u) + 20):
        y = y + y * y + epsilon
        if y.a > b:
            ESCAPE_CACHE[u] = k
            return k
        require(y.b <= b, f"Interval overlap makes first crossing inconclusive at {u}, step {k}")
    raise RuntimeError(f"No crossing within finite-test cap for phase {u}")


def interval_dynamics_controls():
    translations = reflections = 0
    for n in (40, 100, 300, 1000):
        for theta in (Q(1, 10), Q(1, 4), Q(1, 2), Q(3, 4), Q(9, 10)):
            u = n + theta
            Nu = interval_escape(u)
            require(Nu == ceil(u), "Actual escape ceiling sample")
            for q in (-7, -1, 0, 1, 2, 19):
                require(interval_escape(u + q) == Nu + q, "Actual escape translation sample")
                translations += 1
    for m in (120, 121, 122, 300, 1000):
        lo, hi = m // 3, (2 * m) // 3
        for a in (lo + 1, (lo + hi + 1) // 2, hi):
            for theta in (Q(1, 10), Q(1, 2), Q(9, 10)):
                u = a - 1 + theta
                v = m - 1 - u
                require(interval_escape(u) == a, "Actual reflected first center")
                require(interval_escape(v) == m - a, "Actual reflected complement")
                require(interval_escape(u) + interval_escape(v) == m,
                        "Actual exact complementary sum sample")
                reflections += 1
    print(f"DIRECTED-INTERVAL: {translations} actual translation pairs and {reflections} "
          "actual complementary pairs certified")
    print(f"DIRECTED-INTERVAL: {len(ESCAPE_CACHE)} distinct actual escape evaluations; R=2")
    print("Inherited phase certificate: -2.463836664 < C_2 < -2.463836643")
    print("Scope: finite samples only; no numerical K_R, asymptotic onset, or prime lower bound")


if __name__ == "__main__":
    exact_rounding_controls()
    exact_integral_controls()
    analytic_controls()
    interval_dynamics_controls()
    print("PASS: all exact, numerical, and directed-interval controls")

#!/usr/bin/env python3
"""Finite checks for the independently proved cusp clock and phase theorem.

The mpmath.iv computation is a directed-interval phase enclosure with analytic
tail bounds. Other printed experiments are high-precision NUMERICAL controls,
not proofs of exact integer equality or a next-order asymptotic.
Requires mpmath. Run from any directory. No files are modified.
"""
import mpmath as mp


def require(predicate, message):
    if not predicate:
        raise RuntimeError(message)


def phase_interval(m=10_000):
    mp.iv.dps = 50
    I = mp.iv.mpf
    u, v = I("-0.5"), I("1.5")
    for _ in range(m):
        u = u + u * u
        v = 2 * v / (1 + mp.iv.sqrt(1 + 4 * v))
    h = lambda x: -1 / x + mp.iv.log(abs(x))
    cm = h(v) - h(u) + 2 * m
    lower_expr = cm + u / (2 * (1 + u) ** 2) - v / 2
    upper_expr = cm + u / 2 - v / (2 * (1 + v) ** 2)
    require(I("-2.463836664") < lower_expr,
            "Outward lower certificate failed")
    require(upper_expr < I("-2.463836643"),
            "Outward upper certificate failed")
    print("DIRECTED-INTERVAL PHASE ENCLOSURE, m =", m)
    print("lower expression:", lower_expr)
    print("upper expression:", upper_expr)
    print("Certified outward enclosure: -2.463836664 < C_2 < -2.463836643")


def defect(y, delta):
    return (mp.atan(delta / (1 + y)) / delta
            + mp.log((1 + y) ** 2 + delta * delta) / 2 - 1)


def clock(y, delta):
    return mp.atan(y / delta) / delta + mp.log(y * y + delta * delta) / 2


def numerical_identity_controls():
    mp.mp.dps = 70
    for delta in map(mp.mpf, ["0.0001", "0.01", "0.3", "1.0"]):
        for y in map(mp.mpf, ["-0.5", "-0.25", "0", "0.2", "1.5", "10"]):
            q = y * y + delta * delta
            d = defect(y, delta)
            require(abs(clock(y + q, delta) - clock(y, delta) - 1 - d) < mp.mpf("1e-60"),
                    "Exact-clock numerical consistency failed")
            require(0 <= d <= 2 * q, "Defect bound numerical check failed")
    print("NUMERICAL: clock identity and defect bounds passed 24 controls")


def endpoint_above(epsilon, n, b=mp.mpf("1.5")):
    y = mp.mpf("-0.5")
    for _ in range(n):
        y += y * y + epsilon
        if y > b:
            return True
    return y > b


def numerical_threshold(n):
    # Monotonicity is proved in the note. The bisection evaluations here use
    # high precision rather than interval arithmetic; label their output so.
    lo, hi = mp.mpf(0), mp.mpf(2)
    for _ in range(210):
        mid = (lo + hi) / 2
        if endpoint_above(mid, n):
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def numerical_threshold_controls():
    # An approximation inside the certified interval, used only to display
    # the proposed next coefficient. It is not a high-precision certificate.
    phase_hint = mp.mpf("-2.4638366534")
    print("NUMERICAL: threshold phase, adjacent-gap ratio, proposed next coefficient")
    print("n | epsilon_n | n-pi/sqrt(epsilon_n) | gap*n^3/(2*pi^2) | phase remainder/delta")
    for n in [10, 50, 100, 200, 500, 1000]:
        en = numerical_threshold(n)
        en1 = numerical_threshold(n + 1)
        delta = mp.sqrt(en)
        phase = n - mp.pi / delta
        ratio = (en - en1) * n ** 3 / (2 * mp.pi ** 2)
        coeff = (phase - phase_hint) / delta
        require(en > en1 > 0, "Threshold monotonicity numerical check failed")
        print(n, mp.nstr(en, 18), mp.nstr(phase, 16),
              mp.nstr(ratio, 14), mp.nstr(coeff, 14), sep=" | ")
    print("Proposed coefficient pi/3 =", mp.nstr(mp.pi / 3, 16))
    print("The last column is evidence for hypothesis (H), not its proof.")


def numerical_escape(epsilon, threshold=mp.mpf("1.5")):
    y, n = mp.mpf("-0.5"), 0
    while y <= threshold:
        y += y * y + epsilon
        n += 1
    return n


def numerical_base_controls():
    print("NUMERICAL: base-independent scaling controls")
    print("base | m | epsilon | escape N | N/base^m | N-pi*base^m | floor(pi*base^m)-N")
    for base, m in [(2, 3), (2, 7), (2, 12), (3, 4), (10, 1), (10, 2), (10, 3), (10, 4)]:
        scale = mp.mpf(base) ** m
        e = scale ** -2
        n = numerical_escape(e)
        error = n - mp.pi * scale
        require(-11 < error < 3, "Uniform error numerical check failed")
        print(base, m, mp.nstr(e, 12), n, mp.nstr(n / scale, 14),
              mp.nstr(error, 14), int(mp.floor(mp.pi * scale)) - n, sep=" | ")


if __name__ == "__main__":
    phase_interval()
    numerical_identity_controls()
    numerical_threshold_controls()
    numerical_base_controls()
    print("PASS: interval certificate and all finite numerical controls")

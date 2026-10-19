#!/usr/bin/env python3
"""Modest NUMERICAL controls for a separately proved next-order theorem.

No finite computation is used as proof of the asymptotic remainder.
Dependencies: Python 3, mpmath 1.3.0. No input files or network required.
"""
import mpmath as mp

mp.mp.dps = 70


def require(predicate, message):
    if not predicate:
        raise RuntimeError(message)


def improved_clock(y, delta):
    e = delta * delta
    return ((1 + e / 3) * mp.atan(y / delta) / delta
            + mp.log(y * y + e) / 2 - y / 2)


def defect(y, e):
    w = 1 + y
    if not e:
        return 1 / w + mp.log(w) - 1 - y * y / 2
    delta = mp.sqrt(e)
    return ((1 + e / 3) * mp.atan(delta / w) / delta
            + mp.log(w * w + e) / 2 - 1 - (y * y + e) / 2)


def identity_controls():
    for delta in map(mp.mpf, ["0.0001", "0.01", "0.2"]):
        for y in map(mp.mpf, ["-0.5", "-0.1", "0", "0.1", "1.5"]):
            e = delta * delta
            q = y * y + e
            require(abs(improved_clock(y + q, delta) - improved_clock(y, delta)
                        - 1 - defect(y, e)) < mp.mpf("1e-60"),
                    "Improved clock identity numerical failure")
            derivative = (y - e / 3) / ((1 + y) ** 2 + e) - y
            require(abs(mp.diff(lambda z: defect(z, e), y) - derivative)
                    < mp.mpf("1e-60"), "Derivative identity numerical failure")
    print("NUMERICAL: 15 improved-clock and exact derivative controls passed")


def endpoint_controls():
    tested = 0
    tol = mp.mpf("1e-60")  # Permit rounding at exact endpoint equalities.
    for delta in map(mp.mpf, ["0.01", "0.001"]):
        e = delta * delta
        m = int(mp.floor(1 / (10 * delta)))
        u = y = mp.mpf("-0.5")
        for k in range(1, m + 1):
            u += u * u
            y += y * y + e
            require(u - tol <= y <= u + k * e + tol and u + k * e <= u / 2 < 0,
                    "Forward side/perturbation control failed")
        for b in map(mp.mpf, ["0.1", "1.5", "10"]):
            v = z = b
            for j in range(1, m + 1):
                v = 2 * v / (1 + mp.sqrt(1 + 4 * v))
                z = 2 * (z - e) / (1 + mp.sqrt(1 + 4 * (z - e)))
                require(0 < v / 2 <= v - j * e and v - j * e - tol <= z <= v + tol,
                        "Reverse side/perturbation control failed")
            tested += 1
    print("NUMERICAL:", tested, "endpoint-side regimes passed")


def endpoint_above(e, n):
    y = mp.mpf("-0.5")
    for _ in range(n):
        y += y * y + e
        if y > mp.mpf("1.5"):
            return True
    return False


def threshold(n):
    lo, hi = mp.mpf(0), mp.mpf(2)
    for _ in range(210):
        mid = (lo + hi) / 2
        if endpoint_above(mid, n):
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def next_order_controls():
    # This interval comes from the initial packet's independent directed-
    # interval recurrence computation and analytic tails. The roots here are
    # merely high-precision numerical bisections, so this whole table is
    # expressly NUMERICAL, not a new interval certificate.
    c_lo, c_hi = mp.mpf("-2.463836664"), mp.mpf("-2.463836643")
    print("NUMERICAL: threshold residual after universal pi*delta/3 correction")
    print("n | delta_n | residual lower | residual upper | residual/(delta^2 log(1/delta)) midpoint")
    for n in [20, 100, 500, 1000]:
        e = threshold(n)
        delta = mp.sqrt(e)
        r0 = n - mp.pi / delta - mp.pi * delta / 3
        lo, hi = r0 - c_hi, r0 - c_lo
        scale = e * mp.log(1 / delta)
        print(n, mp.nstr(delta, 12), mp.nstr(lo, 14), mp.nstr(hi, 14),
              mp.nstr((lo + hi) / (2 * scale), 12), sep=" | ")


if __name__ == "__main__":
    identity_controls()
    endpoint_controls()
    next_order_controls()
    print("PASS: all finite numerical controls; the theorem rests on its analytic proof")

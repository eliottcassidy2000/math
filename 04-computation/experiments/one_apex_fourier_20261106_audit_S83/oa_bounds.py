# Bound formulas for the S83 audit.
import math
import numpy as np

def thmA(n, q):
    """THM-4886 (A): |S_n| / 2^(n-1) <= min over splits n = a + L + l (a >= 2, L, l >= 1) of
    2^(-L/2) (1 + 2 (2^L/q + 1) sqrt((3 + 2^(3-a) q)(3 + 2^(2-l) q)))^(1/2).  Exact minimum (convexity in a, all L)."""
    best = math.inf; arg = None
    for L in range(1, n - 2):
        N = n - L
        amid = (N + 1) // 2
        for a in range(amid - 2, amid + 3):
            a = min(max(a, 2), N - 1); l = N - a
            if l < 1:
                continue
            R = 2 * (2.0 ** L / q + 1) * math.sqrt((3 + 2.0 ** (3 - a) * q) * (3 + 2.0 ** (2 - l) * q))
            b = 2.0 ** (-L / 2) * math.sqrt(1 + R)
            if b < best:
                best, arg = b, (a, L, l)
    return best, arg

def twoblock(n, q):
    """Audit upgrade (two-block split, Session 193 eq. (1)): |S_n| / 2^(n-1) <= min over 2 <= a <= n-1 of
    (6/q + 7 2^(a-n) + 8 (1 + ln q) 2^(-a) + 8 q (1 + ln q) 2^(-n))^(1/2)."""
    lq = 1 + math.log(q)
    best = math.inf; arg = None
    for a in range(2, n):
        v = 6 / q + 7 * 2.0 ** (a - n) + 8 * lq * 2.0 ** (-a) + 8 * q * lq * 2.0 ** (-n)
        if v < best:
            best, arg = v, a
    return math.sqrt(best), arg

def admissible_fracs(beta, qmax):
    """All reduced r/q, 1 <= q <= qmax, with |beta - r/q| <= 1/q^2 (beta a Fraction)."""
    from fractions import Fraction
    out = []
    for q in range(1, qmax + 1):
        r = round(beta * q)
        if abs(beta - Fraction(r, q)) <= Fraction(1, q * q) and math.gcd(r, q) == 1:
            out.append((r, q))
        # also r +- 1 can qualify only when q <= 2; handle generally
        for rr in (r - 1, r + 1):
            if abs(beta - Fraction(rr, q)) <= Fraction(1, q * q) and math.gcd(rr, q) == 1:
                out.append((rr, q))
    return sorted(set(out), key=lambda t: t[1])

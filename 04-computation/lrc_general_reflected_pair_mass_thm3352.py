#!/usr/bin/env python3
"""Exact O(log min(p,q)) reflected pair mass on the LRC ruler domain.

The bulk is an unbounded periodized trapezoid evaluated through Euclidean
floor moments.  Any complete turns of the tent sum to a constant, leaving
only one residue prefix and one suffix.  For the reflected atlas domain
`L>=168`, `14|L`, and `1<=e,f<=14`, at most the four endpoint teeth
`{-1,0,p-1,p}` need correction; those are replaced by an exact periodic
antiderivative after clipping to `[0,1]`.
"""
from fractions import Fraction as F


def floor_moments(n, m, a, b):
    if n == 0:
        return 0, 0, 0
    s1 = n * (n - 1) // 2
    s2 = n * (n - 1) * (2 * n - 1) // 6
    qa, a0 = divmod(a, m)
    qb, b0 = divmod(b, m)
    base0 = qa * s1 + qb * n
    base1 = qa * s2 + qb * s1
    base2 = qa * qa * s2 + 2 * qa * qb * s1 + qb * qb * n
    if a0 == 0:
        return base0, base1, base2
    height = (a0 * (n - 1) + b0) // m
    if height == 0:
        return base0, base1, base2
    u0, u1, u2 = floor_moments(height, a0, m, m - b0 + a0 - 1)
    r0 = n * height - u0
    r1 = height * s1 - (u2 - u0) // 2
    r2 = n * height * height - 2 * u1 - u0
    return base0 + r0, base1 + r1, base2 + 2 * qa * r1 + 2 * qb * r0 + r2


def residue_prefix(n, m, a, b, threshold, base, total):
    """Count and sum (a*t+b mod m), 0<=t<n, strictly below threshold."""
    if threshold <= 0:
        return 0, 0
    if threshold >= m:
        return n, total
    shifted = floor_moments(n, m, a, b + m - threshold)
    d0, d1 = shifted[0] - base[0], shifted[1] - base[1]
    y0d = (shifted[2] - base[2] - d0) // 2
    high_sum = a * d1 + b * d0 - m * y0d
    return n - d0, total - high_sum


def triangle_sum(n, m, a, b, peak, L, base, total):
    """Sum_t sum_k [peak-L*abs((a*t+b mod m)+k*m)]_+ in O(log m)."""
    if peak <= 0:
        return 0
    radius = (peak - 1) // L
    turns, tail = divmod(radius, m)
    # The 2*turns full lifts have cancelling residue slopes.
    answer = n * (2 * turns * peak - L * m * turns * turns)
    low_count, low_sum = residue_prefix(n, m, a, b, tail + 1, base, total)
    answer += (peak - L * turns * m) * low_count - L * low_sum
    before_count, before_sum = residue_prefix(n, m, a, b, m - tail, base, total)
    high_count, high_sum = n - before_count, total - before_sum
    answer += (peak - L * (turns + 1) * m) * high_count + L * high_sum
    return answer


def one_triangle(residue, m, peak, L):
    if peak <= 0:
        return 0
    radius = (peak - 1) // L
    turns, tail = divmod(radius, m)
    answer = 2 * turns * peak - L * m * turns * turns
    if residue <= tail:
        answer += peak - L * (turns * m + residue)
    if residue >= m - tail:
        answer += peak - L * ((turns + 1) * m - residue)
    return answer


def periodic_primitive(y):
    """Primitive of the 1-periodic radius-1/14 indicator centred on Z."""
    n = y.numerator // y.denominator
    t = y - n
    return F(n, 7) + min(t, F(1, 14)) + max(F(0), t - F(13, 14))


def clipped_tooth_mass(left, right, phase, modulus, L):
    """Intersection with every second-clause tooth, already clipped to [0,1]."""
    lo = F(modulus * left - phase, L)
    hi = F(modulus * right - phase, L)
    return F(L, modulus) * (periodic_primitive(hi) - periodic_primitive(lo))


def mass(L, cell, e, p, f, q):
    if L < 168 or L % 14 or not (1 <= e <= 14 and 1 <= f <= 14 and p >= 1 and q >= 1):
        raise RuntimeError(("outside LRC ruler domain", L, cell, e, p, f, q))
    z, w = L * p - e, L * q - f
    if z > w:
        return mass(L, cell, f, q, e, p)
    r, s = e * cell % L, f * cell % L
    determinant = r * w - s * z
    if determinant % L:
        raise RuntimeError(("nonintegral phase", L, cell, e, p, f, q, determinant))
    b, a = (determinant // L) % z, w % z
    base = floor_moments(p, z, a, b)
    total = a * p * (p - 1) // 2 + b * p - z * base[0]
    unit = L // 14
    outer, inner = unit * (z + w), unit * (w - z)
    numerator = triangle_sum(p, z, a, b, outer, L, base, total)
    numerator -= triangle_sum(p, z, a, b, inner, L, base, total)

    # On the declared domain only -1,0,p-1,p can be boundary/extra teeth.
    radius = F(L, 14 * z)
    for k in {-1, 0, p - 1, p}:
        centre = F(r + k * L, z)
        actual = centre + radius > 0 and centre - radius < 1
        nominal = 0 <= k < p
        boundary = actual and (centre - radius < 0 or centre + radius > 1)
        residue = (a * k + b) % z
        if nominal and (not actual or boundary):
            numerator -= one_triangle(residue, z, outer, L)
            numerator += one_triangle(residue, z, inner, L)
        if actual and (not nominal or boundary):
            left, right = max(F(0), centre - radius), min(F(1), centre + radius)
            numerator += clipped_tooth_mass(left, right, s, w, L) * z * w
    return F(numerator, z * w)


if __name__ == "__main__":
    print("exact general reflected pair mass; import mass(L,cell,e,p,f,q)")

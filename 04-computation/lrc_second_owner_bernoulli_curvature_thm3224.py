#!/usr/bin/env python3
"""Self-contained exact audit for the second reflected-cell LRC corrector.

This is a scratch companion, not a repository truth source.  It checks:

* the exact complete-orbit Bernoulli-B2 formula and its residue-periodic
  g^-2 coefficient;
* the sharp reflected-lane coefficient bound and the absence of a zero
  complete-orbit coefficient;
* the arithmetic-carry sign flip on one fixed orbit;
* the cellwise residue coefficient forced by a quadratic cleared numerator;
* two infinite residue-tail branch certificates, including a cell with zero
  first correction but mixed-sign nonzero second correction; and
* one full 168-owner second-order word and its cyclic Hodge sum.

Only Fraction/integer arithmetic is used.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import gcd


L0 = 168
THETA = F(1, 14)
H = (1, 2, 3, 4, 6, 12)
EDGES = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
         (1, 2), (1, 3), (1, 4), (1, 5))
LANES = tuple((H[i], H[j]) for i, j in EDGES
              for i, j in ((i, j), (j, i)))


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def frac(x):
    return x - x.numerator // x.denominator


def b2bar(x):
    x = frac(x)
    return x * x - x + F(1, 6)


def b3bar(x):
    x = frac(x)
    return x ** 3 - F(3, 2) * x ** 2 + F(1, 2) * x


def comb(frequency, phase=F(0)):
    """Reduced representatives of chi(frequency*x-phase) on [0,1]."""
    phase = frac(phase)
    ans = []
    for k in range(-2, frequency + 3):
        lo = max(F(0), (F(k) + phase - THETA) / frequency)
        hi = min(F(1), (F(k) + phase + THETA) / frequency)
        if lo < hi:
            ans.append((lo, hi))
    ans.sort()
    return ans


def interval_overlap(left, right, centered_moment=False):
    i = j = 0
    total = F(0)
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            if centered_moment:
                total += ((hi - F(1, 2)) ** 2
                          - (lo - F(1, 2)) ** 2) / 2
            else:
                total += hi - lo
        if left[i][1] < right[j][1]:
            i += 1
        elif right[j][1] < left[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return total


def barycenter(P, alpha, Q, beta):
    return interval_overlap(comb(P, alpha), comb(Q, beta), True)


def correction_cell(P, Q, e, f, cell):
    return (barycenter(P, F((cell + 1) * e, L0),
                       Q, F((cell + 1) * f, L0))
            - barycenter(P, F(cell * e, L0),
                         Q, F(cell * f, L0)))


def limit_cell(P, Q, e, f, cell):
    cross = Q * e - P * f
    require(cross != 0, ("rank-one cell limit", P, Q, e, f))
    R, S = cell * e % L0, cell * f % L0
    D = Q * R - P * S
    a, b = F(P, 14), F(Q, 14)
    u, v = F(D + cross, L0), F(-D, L0)
    psi = (b3bar(u + a - b) + b3bar(u - a + b)
           + b3bar(v + a - b) + b3bar(v - a + b)
           - b3bar(u + a + b) - b3bar(u - a - b)
           - b3bar(v + a + b) - b3bar(v - a - b))
    return F(1, 49) + F(28, P * Q * cross) * psi


def located_comb(p, label, cell):
    z, phase = L0 * p - label, cell * label % L0
    ans = []
    for k in range(-2, p + 3):
        lo = max(F(0), F(phase + L0 * k - 12, z))
        hi = min(F(1), F(phase + L0 * k + 12, z))
        if lo < hi:
            ans.append((lo, hi))
    ans.sort()
    return ans


def cell_mass(P, Q, e, f, cell, g):
    return interval_overlap(located_comb(g * P, e, cell),
                            located_comb(g * Q, f, cell))


def direct_cleared(P, Q, e, f, cell, g):
    return (cell_mass(P, Q, e, f, cell, g)
            * (L0 * g * P - e) * (L0 * g * Q - f))


def orbit_data(P, Q, e, f, g):
    d = gcd(L0, gcd(e, f))
    T, a, b = L0 // d, e // d, f // d
    K = abs(Q * e - P * f) // d
    require(K > 0, ("rank-one orbit", P, Q, e, f))
    n, m = g * P * T - a, g * Q * T - b
    h = gcd(n, m)
    require(K % h == 0, ("determinant divisor", P, Q, e, f, g, h, K))
    require(gcd(h, 14) == 1, ("nonunit carry", P, Q, e, f, g, h))
    A, B = n // h, m // h
    inverse = pow(h, -1, 14)
    xi_direct = (b2bar(F(A - B, 14))
                 - b2bar(F(A + B, 14)))
    xi_owner = (b2bar(F((a - b) * inverse, 14))
                - b2bar(F((a + b) * inverse, 14)))
    require(xi_direct == xi_owner,
            ("owner Bernoulli mismatch", P, Q, e, f, g))
    error = F(T * h * h, n * m) * xi_direct
    omega = F(h * h, P * Q * T) * xi_direct
    exact_remainder = (F(h * h) * xi_direct
                       * (T * (P * b + Q * a) * g - a * b)
                       / (P * Q * T * g * g * n * m))
    require(error == omega / (g * g) + exact_remainder,
            ("exact second-order remainder", P, Q, e, f, g))
    return {
        "d": d, "T": T, "a": a, "b": b, "K": K,
        "n": n, "m": m, "h": h, "A": A, "B": B,
        "xi": xi_direct, "error": error, "omega": omega,
        "exact_remainder": exact_remainder,
    }


def direct_closed_orbit(P, Q, e, f, g):
    row = orbit_data(P, Q, e, f, g)
    reduced = row["T"] * interval_overlap(comb(row["A"]), comb(row["B"]))
    formula = F(row["T"], 49) + row["error"]
    require(reduced == formula, ("closed Bernoulli formula", P, Q, e, f, g))
    return reduced


# Exact floor-moment engine and affine-branch stability test.  This is kept
# local so the two declared infinite cell controls do not depend on a repo
# import or on a finite interpolation guess.
def ceildiv(a, b):
    return -((-a) // b)


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
    return (base0 + r0, base1 + r1,
            base2 + 2 * qa * r1 + 2 * qb * r0 + r2)


def residue_prefix(n, m, a, b, threshold, base):
    shifted = floor_moments(n, m, a, b + m - threshold)
    d0 = shifted[0] - base[0]
    d1 = shifted[1] - base[1]
    y0d = (shifted[2] - base[2] - d0) // 2
    high_sum = a * d1 + b * d0 - m * y0d
    total = a * n * (n - 1) // 2 + b * n - m * base[0]
    return n - d0, total - high_sum


def triangle_sum(n, m, a, b, peak, base, total):
    if peak <= 0:
        return 0
    radius = (peak - 1) // L0
    require(2 * radius < m, ("overlapping triangle tails", n, m, peak))
    low_count, low_sum = residue_prefix(n, m, a, b, radius + 1, base)
    before_count, before_sum = residue_prefix(n, m, a, b, m - radius, base)
    high_count, high_sum = n - before_count, total - before_sum
    return (peak * low_count - L0 * low_sum
            + (peak - L0 * m) * high_count + L0 * high_sum)


def fast_cleared(cell, e, p, f, q):
    require(p <= q, ("fast engine order", p, q))
    z, w = L0 * p - e, L0 * q - f
    r, s = e * cell % L0, f * cell % L0
    determinant = r * w - s * z
    require(determinant % L0 == 0,
            ("nonintegral phase", cell, e, p, f, q))
    b = (determinant // L0) % z
    a = w % z
    base = floor_moments(p, z, a, b)
    total = a * p * (p - 1) // 2 + b * p - z * base[0]
    return (triangle_sum(p, z, a, b, 12 * (z + w), base, total)
            - triangle_sum(p, z, a, b, 12 * (w - z), base, total))


def term_candidates(cell, P, Q, e, f, residue_index, shift_index, kernel, g):
    z, w = L0 * g * P - e, L0 * g * Q - f
    r0, s0 = e * cell % L0, f * cell % L0
    base = (r0 * w - s0 * z) // L0
    cross = P * w - Q * z
    lo = max(0, ceildiv(-shift_index, Q))
    hi = min(g - 1, (g * Q - 1 - shift_index) // Q)
    affine = base + residue_index * w - shift_index * z
    peak = 12 * ((z + w) if kernel == 0 else (w - z))
    radius = (peak - 1) // L0
    if cross < 0:
        affine, cross = -affine, -cross
    if cross:
        positive_lo = ceildiv(-radius - affine, cross)
        positive_hi = (radius - affine) // cross
        negative_hi = ceildiv(-affine, cross) - 1
        return (lo, positive_lo, hi, positive_hi,
                negative_hi, negative_hi + 1), cross
    return (lo, hi, affine, peak - L0 * abs(affine)), 0


def stable_term(cell, P, Q, e, f, residue_index, shift_index, kernel,
                residue, period, h):
    here, cross = term_candidates(
        cell, P, Q, e, f, residue_index, shift_index, kernel,
        residue + period * h)
    nxt, cross2 = term_candidates(
        cell, P, Q, e, f, residue_index, shift_index, kernel,
        residue + period * (h + 1))
    require(cross == cross2, ("cross instability", cell, residue, period, h))
    slopes = tuple(y - x for x, y in zip(here, nxt))
    nxt2, _ = term_candidates(
        cell, P, Q, e, f, residue_index, shift_index, kernel,
        residue + period * (h + 2))
    require(tuple(y - x for x, y in zip(nxt, nxt2)) == slopes,
            ("endpoint nonaffinity", cell, residue, period, h))
    if cross:
        lo, positive_lo, hi, positive_hi, negative_hi, _ = here
        slo, splo, shi, sphi, sneg, _ = slopes
        if positive_lo > lo or (positive_lo == lo and splo > slo):
            left, sleft = positive_lo, splo
        else:
            left, sleft = lo, slo
        if not (sleft >= slo and sleft >= splo):
            return False
        if positive_hi < hi or (positive_hi == hi and sphi < shi):
            right, sright = positive_hi, sphi
        else:
            right, sright = hi, shi
        if not (sright <= shi and sright <= sphi):
            return False
        if left > right:
            return sleft >= sright
        if negative_hi < left:
            return sneg <= sleft
        if negative_hi >= right:
            return sneg >= sright
        return sleft <= sneg <= sright
    lo, hi, affine, value = here
    slo, shi, saffine, svalue = slopes
    if hi < lo or shi < slo:
        return False
    if affine > 0 and saffine < 0:
        return False
    if affine < 0 and saffine > 0:
        return False
    if affine == 0 and saffine != 0:
        return False
    if value > 0:
        return svalue >= 0
    if value < 0:
        return svalue <= 0
    return svalue == 0


def ray_is_stable(cell, P, Q, e, f, residue, period, h):
    for residue_index in range(P):
        for shift_index in range(-2, Q + 3):
            for kernel in (0, 1):
                if not stable_term(cell, P, Q, e, f, residue_index,
                                   shift_index, kernel, residue, period, h):
                    return False
    return True


def interpolate(points):
    (x0, y0), (x1, y1), (x2, y2) = points
    den = (x0 - x1) * (x0 - x2)
    aa = (F(y0, den) + F(y1, (x1 - x0) * (x1 - x2))
          + F(y2, (x2 - x0) * (x2 - x1)))
    bb = (-F(y0 * (x1 + x2), den)
          - F(y1 * (x0 + x2), (x1 - x0) * (x1 - x2))
          - F(y2 * (x0 + x1), (x2 - x0) * (x2 - x1)))
    cc = (F(y0 * x1 * x2, den)
          + F(y1 * x0 * x2, (x1 - x0) * (x1 - x2))
          + F(y2 * x0 * x1, (x2 - x0) * (x2 - x1)))
    return aa, bb, cc


def polynomial_value(poly, g):
    return poly[0] * g * g + poly[1] * g + poly[2]


def compile_infinite_cell(cell, P, Q, e, f):
    period = abs(Q * e - P * f) or 1
    polynomials = []
    for residue in range(1, period + 1):
        require(ray_is_stable(cell, P, Q, e, f, residue, period, 1),
                ("unstable declared cell", cell, P, Q, e, f, residue))
        points = []
        for h in range(1, 5):
            g = residue + period * h
            points.append((g, fast_cleared(cell, e, g * P, f, g * Q)))
        poly = interpolate(points[:3])
        require(polynomial_value(poly, points[3][0]) == points[3][1],
                ("fourth tail point", cell, residue, poly))
        # The one omitted head in this affine branch agrees too, so the
        # formula covers every positive g in the residue class.
        head = fast_cleared(cell, e, residue * P, f, residue * Q)
        require(polynomial_value(poly, residue) == head,
                ("residue head", cell, residue, poly, head))
        require(F(head) == direct_cleared(P, Q, e, f, cell, residue),
                ("independent geometry head", cell, residue))
        polynomials.append(poly)
    return tuple(polynomials)


def q_from_kappa(P, Q, e, f, limit, correction, kappa):
    d2 = L0 * L0 * P * Q
    d1 = -L0 * (P * f + Q * e)
    d0 = e * f
    return (kappa - d0 * limit - d1 * correction) / d2


# 1. The Bernoulli owner table.  h is always a unit modulo 14; this table is
# a complete finite audit over the 18 normalized reflected lanes.
units14 = (1, 3, 5, 9, 11, 13)
lane_xi = set()
for e, f in LANES:
    d = gcd(L0, gcd(e, f))
    a, b = e // d, f // d
    for h in units14:
        inv = pow(h, -1, 14)
        lane_xi.add(b2bar(F((a - b) * inv, 14))
                    - b2bar(F((a + b) * inv, 14)))
expected_lane_xi = {F(x, 49) for x in (-8, -6, -5, -3, 1, 2, 3, 4, 5, 8, 9)}
require(lane_xi == expected_lane_xi, ("lane Xi table", lane_xi))
require(F(0) not in lane_xi and max(map(abs, lane_xi)) == F(9, 49),
        ("lane zero/bound", lane_xi))


# 2. Exact direct controls for the complete-orbit formula.
closed_controls = ((3, 5, 1, 2, 1),
                   (3, 5, 1, 12, 2),
                   (3, 5, 1, 12, 4),
                   (4, 5, 1, 6, 11))
for control in closed_controls:
    direct_closed_orbit(*control)

negative = orbit_data(3, 5, 1, 12, 2)
require(negative["error"] == -F(10, 979811),
        ("negative aggregate hostile", negative))
sharp = orbit_data(4, 5, 1, 6, 11)
require(sharp["h"] == 19 and sharp["xi"] == F(9, 49)
        and sharp["omega"] == F(1083, 54880),
        ("sharp 9/49 control", sharp))

# The same orbit has opposite second-order signs according to one gcd carry.
switch_rows = []
for g in range(1, 63):
    row = orbit_data(3, 5, 1, 12, g)
    if g % 31 == 4:
        require((row["h"], row["xi"], row["omega"])
                == (31, F(8, 49), F(961, 15435)),
                ("exceptional carry", g, row))
    else:
        require((row["h"], row["xi"], row["omega"])
                == (1, -F(5, 49), -F(1, 24696)),
                ("generic carry", g, row))
    switch_rows.append((g, row["h"], row["xi"], row["omega"]))

# Direct concatenation check for the old negative hostile: cell owner sum is
# the closed geodesic, but the individual owner distribution is retained here.
cell_sum = sum(cell_mass(3, 5, 1, 12, j, 2) for j in range(168))
require(cell_sum == direct_closed_orbit(3, 5, 1, 12, 2),
        ("cell concatenation", cell_sum))


# 3. Broad exact finite arithmetic-carry census.  The theorem is analytic;
# this census is a hostile control for period, signs, nonzero, and the lane
# bound rather than the source of the universal quantifiers.
census = []
signs = {"positive": 0, "negative": 0, "zero": 0}
for Q in range(2, 26):
    for P in range((Q + 1) // 2, Q):
        if gcd(P, Q) != 1 or P + Q < 8:
            continue
        for e, f in LANES:
            d = gcd(L0, gcd(e, f))
            K = abs(Q * e - P * f) // d
            for g in range(1, K + 1):
                row = orbit_data(P, Q, e, f, g)
                periodic = orbit_data(P, Q, e, f, g + K)
                require((row["h"], row["xi"], row["omega"])
                        == (periodic["h"], periodic["xi"], periodic["omega"]),
                        ("K-period failure", P, Q, e, f, g, K))
                require(row["xi"] in lane_xi,
                        ("Xi outside lane table", P, Q, e, f, g, row))
                sign = ("positive" if row["omega"] > 0
                        else "negative" if row["omega"] < 0 else "zero")
                signs[sign] += 1
                census.append((P, Q, e, f, g, row["K"], row["h"],
                               row["xi"], row["omega"]))
require(signs["positive"] and signs["negative"] and signs["zero"] == 0,
        ("aggregate sign census", signs))


# 4. Minimal endpoint-only obstruction at one fixed cell.  The static phase
# endpoints, L, and c do not change with g, but q has three residue values.
cell90_polys = compile_infinite_cell(90, 3, 5, 6, 1)
expected90_kappa = tuple((F(-336), F(84), F(0))[(r - 1) % 3]
                         for r in range(1, 28))
require(tuple(poly[2] for poly in cell90_polys) == expected90_kappa,
        ("cell90 kappa word", cell90_polys))
L90 = limit_cell(3, 5, 6, 1, 90)
c90 = correction_cell(3, 5, 6, 1, 90)
require((L90, c90) == (F(17, 1680), -F(8213, 1411200)),
        ("cell90 limit/correction", L90, c90))
q90 = tuple(q_from_kappa(3, 5, 6, 1, L90, c90, k)
            for k in expected90_kappa[:3])
require(q90 == (-F(343771, 395136000), F(48229, 395136000),
                -F(10057, 131712000)),
        ("cell90 q word", q90))


# 5. Zero first correction does not kill the second owner word.  Exact branch
# stability at h=1 plus the checked omitted heads proves this nine-periodic
# formula for every positive g, not just for a sampled prefix.
zero_polys = compile_infinite_cell(56, 5, 6, 1, 12)
zero_kappa9 = (F(-64, 3), F(2188, 3), F(920), F(1076, 3),
               F(1312, 3), F(-44), F(1820, 3), F(3460, 3), F(0))
require(tuple(poly[2] for poly in zero_polys)
        == tuple(zero_kappa9[(r - 1) % 9] for r in range(1, 55)),
        ("zero-correction kappa word", zero_polys))
Lzero = limit_cell(5, 6, 1, 12, 56)
czero = correction_cell(5, 6, 1, 12, 56)
require((Lzero, czero) == (F(4, 189), F(0)),
        ("zero first correction", Lzero, czero))
qzero9 = tuple(q_from_kappa(5, 6, 1, 12, Lzero, czero, k)
               for k in zero_kappa9)
expected_qzero9 = (-F(17, 666792), F(11483, 13335840),
                   F(7243, 6667920), F(1129, 2667168),
                   F(1721, 3333960), -F(697, 13335840),
                   F(9551, 13335840), F(18161, 13335840),
                   -F(1, 3333960))
require(qzero9 == expected_qzero9,
        ("zero-correction q word", qzero9))
require(any(x > 0 for x in qzero9) and any(x < 0 for x in qzero9)
        and all(x != 0 for x in qzero9),
        ("zero-correction sign hostile", qzero9))


# 6. FINITE-EXACT 168-owner scout.  Five separated values fit and hostile-
# check one period-one quadratic per owner, but no stabilization bound proves
# that these samples are already on THM-3200's eventual branch.  This scout is
# therefore a control only and is not an input to the general Hodge theorem.
owner_q = []
P, Q, e, f = 3, 5, 1, 2
d2 = L0 * L0 * P * Q
d1 = -L0 * (P * f + Q * e)
for cell in range(168):
    limit = limit_cell(P, Q, e, f, cell)
    correction = correction_cell(P, Q, e, f, cell)
    points = tuple((g, direct_cleared(P, Q, e, f, cell, g))
                   for g in (80, 81, 82))
    poly = interpolate(points)
    require(polynomial_value(poly, 83)
            == direct_cleared(P, Q, e, f, cell, 83),
            ("owner fourth point", cell, poly))
    require(polynomial_value(poly, 160)
            == direct_cleared(P, Q, e, f, cell, 160),
            ("owner separated point", cell, poly))
    require(poly[0] == d2 * limit
            and poly[1] == d2 * correction + d1 * limit,
            ("owner top coefficients", cell, poly, limit, correction))
    owner_q.append(q_from_kappa(P, Q, e, f, limit, correction, poly[2]))
owner_q = tuple(owner_q)
owner_omega = orbit_data(P, Q, e, f, 1)["omega"]
require(sum(owner_q) == owner_omega == F(1, 24696),
        ("owner holonomy", sum(owner_q), owner_omega))
require(all(owner_q[167 - j] == owner_q[j] for j in range(168)),
        ("owner reflection palindrome", owner_q))
owner_signs = (sum(x > 0 for x in owner_q),
               sum(x < 0 for x in owner_q),
               sum(x == 0 for x in owner_q))
require(owner_signs == (156, 12, 0), ("owner signs", owner_signs))

# Finite-scout cyclic potential after removing the nonzero harmonic mean.
potential = [F(0)]
for value in owner_q:
    potential.append(potential[-1] + value - owner_omega / 168)
require(potential[-1] == 0, ("Hodge closure", potential[-1]))
require(all(potential[168 - j] == -potential[j] for j in range(169)),
        ("Hodge oddness", potential))


census_digest = sha256("\n".join(map(str, census)).encode()).hexdigest()
owner_digest = sha256("\n".join(map(str, owner_q)).encode()).hexdigest()
print("LRC SECOND-OWNER BERNOULLI CURVATURE EXACT AUDIT")
print(f"lane_xi={tuple(sorted(lane_xi))};sharp_abs=9/49;zero=NONE")
print(f"complete_orbit_controls={closed_controls}")
print("orbit_sign_switch=(3,5;1,12): generic -1/24696; "
      "g=4 mod 31 gives 961/15435")
print(f"census_states={len(census)};Q<=25;signs={signs};sha256={census_digest}")
print(f"cell90_q_mod3={q90};static_endpoint_functional=REFUTED")
print(f"zero_first_correction_cell=(5,6;1,12;j=56);q_mod9={qzero9}")
print(f"owner_word_FINITE_SCOUT=(3,5;1,2);signs={owner_signs};sum={sum(owner_q)};"
      f"sha256={owner_digest}")
print("SECOND_ORDER_HODGE=Bernoulli_gcd_carry_harmonic_plus_residue_owner_coboundary")
print("CELLWISE_SAFETY=NOT_CLAIMED;FINITE_HEAD_REPLACEMENT=NOT_CLAIMED;")
print("PHYSICAL_LRC14=NOT_CLAIMED")
print("FAILED_CHECKS=NONE")

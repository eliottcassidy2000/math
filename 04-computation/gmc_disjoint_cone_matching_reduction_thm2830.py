#!/usr/bin/env python3
"""Exact referee for the THM-2830 disjoint-cone theorem.

The universal proof bypasses the stronger matching inequality.  This
companion verifies the factorial identities, bounded theorem/equality
controls, and the exact nonnegative-coefficient certificates used in the
cyclic-split proof.  Every truth-bearing check uses an explicit exception.
"""

from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import comb, factorial

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def h(a, b):
    return comb(a + b, a)


def multinomial3(a, b, c):
    return factorial(a + b + c) // (factorial(a) * factorial(b) * factorial(c))


def triple(a, b, c):
    return (
        multinomial3(a + 1, b, c)
        + multinomial3(a, b + 1, c)
        + multinomial3(a, b, c + 1)
        - multinomial3(a, b, c)
    )


def tau(a, b, c):
    total = a + b + c
    return (
        Fraction(total + 1, a + 1)
        + Fraction(total + 1, b + 1)
        + Fraction(total + 1, c + 1)
        - 1
    )


def p_sum(i, values):
    ans = 0
    for slot in range(4):
        rest = values[:slot] + values[slot + 1 :]
        ans += h(i, values[slot]) * triple(*rest)
    return ans


def matching_gap(i, a, b, c, d):
    values = (a, b, c, d)
    return p_sum(i, values) - 3 * (
        triple(i, a, b) * h(c, d) + triple(i, c, d) * h(a, b)
    )


def normalized_matching_gap(i, a, b, c, d):
    return (
        factorial(i)
        * factorial(a)
        * factorial(b)
        * factorial(c)
        * factorial(d)
        * matching_gap(i, a, b, c, d)
    )


def all_pair_sum(i, values):
    total = 0
    for mask in range(16):
        if mask.bit_count() != 2:
            continue
        left = [values[j] for j in range(4) if (mask >> j) & 1]
        right = [values[j] for j in range(4) if not ((mask >> j) & 1)]
        total += triple(i, *left) * h(*right)
    return total


def polarized(i, values):
    return Fraction(p_sum(i, values) - all_pair_sum(i, values), 2)


def poly_mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            out[j + k] += x * y
    return out


def poly_add(a, b):
    out = [Fraction(0)] * max(len(a), len(b))
    for j in range(len(out)):
        out[j] = (a[j] if j < len(a) else 0) + (b[j] if j < len(b) else 0)
    return out


def poly_scale(scale, poly):
    return [scale * value for value in poly]


def moment(poly):
    return sum(x * factorial(j) for j, x in enumerate(poly))


def d_poly(index):
    out = [Fraction(0)] * (index + 2)
    out[index] = -Fraction(1, factorial(index))
    out[index + 1] = Fraction(1, factorial(index + 1))
    return out


def direct_triple(a, b, c):
    return moment(poly_mul(poly_mul(d_poly(a), d_poly(b)), d_poly(c)))


def direct_orientation(i, upper_atoms):
    v = [Fraction(0)]
    for atom in upper_atoms:
        v = poly_add(v, d_poly(atom))
    v2 = poly_mul(v, v)
    v3 = poly_mul(v2, v)
    di = d_poly(i)
    return (
        2 * moment(v3) * moment(poly_mul(di, v))
        - 3 * moment(poly_mul(di, v2)) * moment(v2)
    )


def direct_polarized(i, values):
    total = Fraction(0)
    for mask in range(16):
        chosen = [values[j] for j in range(4) if (mask >> j) & 1]
        sign = -1 if (4 - len(chosen)) % 2 else 1
        total += sign * direct_orientation(i, chosen)
    return total / factorial(4)


def pair_product_control(p, q):
    left = poly_mul(d_poly(p), d_poly(q))
    degree = p + q
    right = [Fraction(0)] * len(left)
    right[degree] += Fraction(comb(degree, p), factorial(degree))
    atom = d_poly(degree + 1)
    scale = comb(degree + 2, p + 1)
    for j, value in enumerate(atom):
        right[j] += scale * value
    return left == right


def falling(x, order):
    if x < order:
        return 0
    return factorial(x) // factorial(x - order)


def pair_factorial_hostile():
    i, y, z, order = 0, 3, 44, 4
    total = y + z
    rho = Fraction(196, 705)
    low = Fraction(falling(total - 1, order), 5)
    spike = falling(total + 1, order)
    phi = (low + rho * spike) / (1 + rho)
    endpoint = Fraction(falling(y, order) + falling(z, order), 2)
    return phi, endpoint, phi - endpoint


def laguerre_hostile(t):
    return 2 * t * (5640 * t * t + 371 * t - 3)


def kernel(index, argument):
    return comb(index + argument, index)


def ratio_a(index, upper):
    return sum(kernel(index, atom) for atom in upper)


def ratio_b(index, upper):
    return sum(triple(index, left, right) for left in upper for right in upper)


def strict_ratio_control(i, m, upper):
    ai = ratio_a(i, upper)
    am = ratio_a(m, upper)
    bi = ratio_b(i, upper)
    bm = ratio_b(m, upper)
    return ai * bm - bi * am


def rising(x, order):
    value = 1
    for offset in range(order):
        value *= x + offset
    return value


def newton_source_term(i, order, source, left, right):
    total = left + right
    c_order = Fraction(i + 1, i + order + 1)
    f_source = factorial(i + source) * factorial(i + total)
    t_source = Fraction(
        (i + 1) * (total + 2) * (i + total + 1),
        (left + 1) * (right + 1),
    )
    low = total * (
        c_order * falling(total - 1, order) - falling(source, order)
    )
    high = t_source * (
        falling(total + 1, order) - falling(source, order)
    )
    return f_source * (low + high)


def newton_bridge(i, m, source, left, right):
    common = Fraction(
        1,
        factorial(i)
        * factorial(i + 1)
        * factorial(source)
        * factorial(left)
        * factorial(right),
    )
    expansion = sum(
        Fraction(comb(m - i, order), rising(i + 1, order))
        * newton_source_term(i, order, source, left, right)
        for order in range(1, m - i + 1)
    )
    return common * expansion


def symbolic_certificate(expression, generators, term_count, minimum, constant=None):
    polynomial = sp.Poly(sp.expand(expression), *generators)
    raw_coefficients = [value for _, value in polynomial.terms()]
    require(
        all(value.is_Integer for value in raw_coefficients),
        "symbolic integral coefficients",
    )
    coefficients = [int(value) for value in raw_coefficients]
    require(len(coefficients) == term_count, "symbolic term count")
    require(min(coefficients) == minimum, "symbolic minimum coefficient")
    require(all(value >= 0 for value in coefficients), "symbolic coefficient sign")
    if constant is not None:
        require(polynomial.coeff_monomial(1) == constant, "symbolic constant")


def symbolic_cyclic_split_certificates():
    """Return the exact certificate statistics used in the proof."""

    i, w, e, d, h0, u, big_d, sigma, tail = sp.symbols(
        "i w e d h0 u big_d sigma tail", nonnegative=True, integer=True
    )
    p_var, z_var, x_var = sp.symbols(
        "p_var z_var x_var", nonnegative=True, integer=True
    )

    def symbolic_falling(value, order):
        return sp.prod(value - offset for offset in range(order))

    c1 = (i + 1) / (i + 2)
    c2 = (i + 1) / (i + 3)
    y_var = p_var + z_var - 1
    high_var = x_var + z_var + 1
    base_cross = sp.factor(
        (symbolic_falling(high_var, 2) - symbolic_falling(p_var, 2))
        * (x_var - c1 * y_var)
        - (high_var - p_var)
        * (symbolic_falling(x_var, 2) - c2 * symbolic_falling(y_var, 2))
    )

    p0 = i + 1 + w
    z0 = p0 + e
    y0 = p0 + z0 - 1
    x0 = y0 + d
    base_ge = (
        d * e * i
        + 3 * d * e
        + 4 * d * i**2
        + 2 * d * i * w
        + 17 * d * i
        + 6 * d * w
        + 15 * d
        + 2 * e**2
        + 11 * e * i
        + 8 * e * w
        + 11 * e
        + 14 * i**2
        + 22 * i * w
        + 25 * i
        + 8 * w**2
        + 22 * w
        + 9
    )
    require(
        sp.factor(
            base_cross.subs({p_var: p0, z_var: z0, x_var: x0})
            - (d + 2 * e + 2 * i + 2 * w + 2)
            * base_ge
            / ((i + 2) * (i + 3))
        )
        == 0,
        "ratio base x>=Y",
    )
    symbolic_certificate(base_ge, (i, w, e, d), 17, 1)

    x_low = (i + 1) * h0 + u
    z_low = x_low + h0 - p_var + 1
    high_low = x_low + z_low + 1
    base_lt = (
        h0**2 * i**2
        + 3 * h0**2 * i
        + 2 * h0**2
        + h0 * i**2
        + 3 * h0 * i * u
        + 3 * h0 * i
        + 5 * h0 * u
        + 2 * h0
        + 3 * i * u
        + 2 * u**2
        + 7 * u
    )
    require(
        sp.factor(
            base_cross.subs({x_var: x_low, z_var: z_low})
            - (high_low - p_var) * base_lt / ((i + 2) * (i + 3))
        )
        == 0,
        "ratio base x<Y",
    )
    symbolic_certificate(base_lt, (i, h0, u), 11, 1)

    p0 = i + 1 + w
    z0 = p0 + e

    def t_prefactor(x_value):
        return (
            (i + 1)
            * (x_value + z0 + 2)
            * (i + x_value + z0 + 1)
            / ((x_value + 1) * (z0 + 1))
        )

    def m1_expression(delta, ratio_lower):
        x_value = p0 + delta
        return (
            3
            * ratio_lower
            * t_prefactor(x_value)
            * (x_value + z0 + 1 - p0)
            - 8 * (p0 + z0) * (x_value - c1 * (p0 + z0 - 1))
        )

    numerator, denominator = sp.together(
        m1_expression(sp.Integer(0), sp.Integer(1))
    ).as_numer_denom()
    require(
        sp.factor(denominator - (i + 2) * (i + w + 2)) == 0,
        "M1 d=0 denominator",
    )
    symbolic_certificate(numerator, (i, w, e), 27, 8, 40)

    delta = 1 + big_d
    x_value = p0 + delta
    quotient_increment = z0 / (i + x_value)
    ratio_lower = (
        1
        + delta * quotient_increment
        + delta * (delta - 1) * quotient_increment**2 / 2
    )
    numerator, denominator = sp.together(
        m1_expression(delta, ratio_lower)
    ).as_numer_denom()
    require(
        sp.factor(
            denominator
            - 2
            * (i + 2)
            * (big_d + i + w + 3)
            * (big_d + 2 * i + w + 2) ** 2
            * (e + i + w + 2)
        )
        == 0,
        "M1 d>=1 denominator",
    )
    symbolic_certificate(numerator, (i, w, e, big_d), 424, 3, 2016)

    delta = z0 + 2 + sigma
    x_value = p0 + delta
    quotient_increment = z0 / (i + x_value)
    ratio_lower = (
        1
        + delta * quotient_increment
        + delta * (delta - 1) * quotient_increment**2 / 2
    )
    t_p = (x_value + z0 + 2) * (i + x_value + z0 + 1) / (x_value + 1)
    t_x = (p0 + z0 + 2) * (i + p0 + z0 + 1) / (p0 + 1)
    m2_expression = (
        ratio_lower * t_p * (delta + z0 + 1)
        - 8 * t_x * (delta - z0 - 1)
    )
    numerator, denominator = sp.together(m2_expression).as_numer_denom()
    require(
        sp.factor(
            denominator
            - 2
            * (i + w + 2)
            * (sigma + e + 2 * i + 2 * w + 5)
            * (sigma + e + 3 * i + 2 * w + 4) ** 2
        )
        == 0,
        "M2 denominator",
    )
    for value, constant in ((0, 10680), (1, 10944), (2, 8064), (3, 4224)):
        symbolic_certificate(
            numerator.subs(sigma, value), (i, w, e), 164, 8, constant
        )
    symbolic_certificate(
        numerator.subs(sigma, tail + 4), (i, w, e, tail), 474, 1, 5400
    )

    return {
        "ratio_base_ge_terms": 17,
        "ratio_base_lt_terms": 11,
        "m1_d0_terms": 27,
        "m1_positive_d_terms": 424,
        "m2_small_terms_each": 164,
        "m2_tail_terms": 474,
    }


def symbolic_global_quarter_certificates():
    """Certify the global adjacent-row beta/Abel prefix polynomials."""

    n, m, a, b = sp.symbols(
        "n m a b", nonnegative=True, integer=True
    )
    x = m + 2
    y = x + a
    z = y + b
    total = x + y + z

    def prefix_atom(value):
        return (
            (n + 1)
            * (n + 2)
            * value
            * (total - value)
            * (total - 2 * value)
            * (n + total - value - 1)
            + 4
            * x
            * y
            * z
            * (total - value - 2)
            * (
                (n + 1) * total
                - (2 * n + 3) * value
                - 2 * n
                - 1
            )
        )

    px = sp.expand(prefix_atom(x))
    pxy = sp.expand(prefix_atom(x) + prefix_atom(y))
    pxyz = sp.expand(prefix_atom(x) + prefix_atom(y) + prefix_atom(z))
    symbolic_certificate(px, (n, m, a, b), 140, 1, 32)
    symbolic_certificate(pxy, (n, m, a, b), 145, 1, 64)
    symbolic_certificate(pxyz, (n, m, a, b), 136, 1, 96)

    def high_atom(value):
        return (
            value
            * (total - value)
            * (total - 2 * value)
            * (n + total - value - 1)
        )

    qx = sp.expand(high_atom(x))
    qxy = sp.expand(high_atom(x) + high_atom(y))
    qxyz = sp.expand(high_atom(x) + high_atom(y) + high_atom(z))
    symbolic_certificate(qx, (n, m, a, b), 46, 1, 48)
    symbolic_certificate(qxy, (n, m, a, b), 51, 1, 96)
    symbolic_certificate(qxyz, (n, m, a, b), 46, 1, 144)
    return {
        "global_high_qx_terms": 46,
        "global_high_qxy_terms": 51,
        "global_high_qxyz_terms": 46,
        "global_quarter_px_terms": 140,
        "global_quarter_pxy_terms": 145,
        "global_quarter_pxyz_terms": 136,
    }


def adjacent_cyclic_parts(n, p, q, r):
    total = p + q
    high = (
        Fraction(comb(total + 2, p + 1) * kernel(n, total + 1) * kernel(n, r), n + 1)
        * (total + 1 - r)
    )
    low = (
        Fraction(comb(total, p) * kernel(n, r) * comb(n + total, n + 1), n + 1)
        * (Fraction((n + 1) * (total - 1), n + 2) - r)
    )
    return high, low


def global_quarter_scaled(n, p, q, r):
    x, y, z = p + 1, q + 1, r + 1
    total = x + y + z
    high = 0
    low = Fraction(0)
    for value in (x, y, z):
        high += (
            value
            * (total - value)
            * (total - 2 * value)
            * factorial(n + value - 1)
            * factorial(n + total - value - 1)
        )
        low += (
            Fraction(x * y * z, (n + 1) * (n + 2))
            * (total - value - 2)
            * factorial(n + value - 1)
            * factorial(n + total - value - 2)
            * (
                (n + 1) * total
                - (2 * n + 3) * value
                - 2 * n
                - 1
            )
        )
    common = Fraction(
        1,
        factorial(n) ** 2
        * factorial(x)
        * factorial(y)
        * factorial(z)
        * (n + 1),
    )
    return common * high, common * low


for a in range(7):
    for b in range(7):
        for c in range(7):
            require(triple(a, b, c) == direct_triple(a, b, c), "triple tensor")

for p in range(8):
    for q in range(8):
        require(pair_product_control(p, q), "pair product")

require(
    tau(3, 4, 5)
    == Fraction(
        triple(3, 4, 5) * factorial(3) * factorial(4) * factorial(5),
        factorial(12),
    ),
    "tau normalization",
)

polarization_controls = 0
for i in range(3):
    for values in combinations_with_replacement(range(i + 1, i + 5), 4):
        require(polarized(i, values) == direct_polarized(i, values), "polarization")
        polarization_controls += 1

matching_cells = 0
strict_cells = 0
equality_cells = 0
forward_cells = 0
minimum_positive = None

for i in range(5):
    values_range = range(i + 1, i + 9)
    for values in combinations_with_replacement(values_range, 4):
        a, b, c, d = values
        matchings = (
            (a, b, c, d),
            (a, c, b, d),
            (a, d, b, c),
        )
        gaps = []
        for matching in matchings:
            gap = matching_gap(i, *matching)
            matching_cells += 1
            require(gap >= 0, "matching gap")
            gaps.append(gap)
            if gap == 0:
                equality_cells += 1
                require(values == (i + 1,) * 4, "unexpected equality")
            else:
                strict_cells += 1
                if minimum_positive is None or gap < minimum_positive:
                    minimum_positive = gap
            before = normalized_matching_gap(i, *matching)
            for slot in range(4):
                raised = list(matching)
                raised[slot] += 1
                after = normalized_matching_gap(i, *raised)
                forward_cells += 1
                require(after >= before, "normalized forward difference")
        require(polarized(i, values) >= 0, "polarized coefficient")
        require(sum(gaps) == 3 * (2 * polarized(i, values)), "matching average")

phi, endpoint, hostile_gap = pair_factorial_hostile()
require(phi == Fraction(1467522360, 901), "pair hostile phi")
require(endpoint == 1629012, "pair hostile endpoint")
require(hostile_gap == -Fraction(217452, 901), "pair hostile gap")

laguerre_value = laguerre_hostile(Fraction(1, 1000))
require(laguerre_value < 0, "Laguerre hostile")

cone_cells = 0
cone_equalities = 0
for boundary in range(1, 5):
    lower_vectors = product(range(3), repeat=boundary)
    for lower in lower_vectors:
        if not any(lower):
            continue
        for upper in product(range(3), repeat=3):
            if not any(upper):
                continue
            u = [Fraction(0)]
            v = [Fraction(0)]
            for index, coefficient in enumerate(lower):
                u = poly_add(u, poly_scale(coefficient, d_poly(index)))
            for offset, coefficient in enumerate(upper):
                v = poly_add(v, poly_scale(coefficient, d_poly(boundary + offset)))
            v2 = poly_mul(v, v)
            value = (
                2 * moment(poly_mul(v2, v)) * moment(poly_mul(u, v))
                - 3 * moment(poly_mul(u, v2)) * moment(v2)
            )
            cone_cells += 1
            require(value >= 0, "cone orientation")
            if value == 0:
                cone_equalities += 1
                require(
                    all(coefficient == 0 for coefficient in lower[:-1])
                    and lower[-1] > 0
                    and upper[0] > 0
                    and upper[1:] == (0, 0),
                    "cone equality",
                )

pascal_identity_cells = 0
strict_ratio_cells = 0
for i in range(5):
    for length in range(1, 4):
        for upper in combinations_with_replacement(range(i + 1, i + 6), length):
            quadratic = sum(h(left, right) for left in upper for right in upper)
            cubic = sum(
                triple(left, middle, right)
                for left in upper
                for middle in upper
                for right in upper
            )
            require(
                quadratic == 2 * sum(ratio_a(atom - 1, upper) for atom in upper),
                "quadratic Pascal identity",
            )
            require(
                cubic == 3 * sum(ratio_b(atom - 1, upper) for atom in upper),
                "cubic Pascal identity",
            )
            pascal_identity_cells += 1
            for m in range(i + 1, i + 6):
                require(strict_ratio_control(i, m, upper) > 0, "strict R ratio")
                strict_ratio_cells += 1

newton_bridge_cells = 0
cyclic_aggregation_cells = 0
for i in range(4):
    for m in range(i + 1, i + 5):
        for source in range(i + 1, i + 6):
            for left in range(i + 1, i + 6):
                for right in range(i + 1, i + 6):
                    direct = (
                        h(i, source) * triple(m, left, right)
                        - h(m, source) * triple(i, left, right)
                    )
                    require(
                        newton_bridge(i, m, source, left, right) == direct,
                        "Newton source bridge",
                    )
                    newton_bridge_cells += 1

        for values in combinations_with_replacement(range(i + 1, i + 6), 3):
            p, q, r = values
            for order in range(1, m - i + 1):
                if p == r:
                    aggregate = newton_source_term(i, order, p, p, p)
                elif p == q:
                    aggregate = (
                        2 * newton_source_term(i, order, p, p, r)
                        + newton_source_term(i, order, r, p, p)
                    )
                elif q == r:
                    aggregate = (
                        newton_source_term(i, order, p, q, q)
                        + 2 * newton_source_term(i, order, q, p, q)
                    )
                else:
                    aggregate = 2 * (
                        newton_source_term(i, order, p, q, r)
                        + newton_source_term(i, order, q, p, r)
                        + newton_source_term(i, order, r, p, q)
                    )
                require(aggregate >= 0, "cyclic Newton aggregation")
                if order == 1:
                    require(aggregate > 0, "strict first Newton aggregation")
                cyclic_aggregation_cells += 1

adjacent_part_cells = 0
global_quarter_cells = 0
for n in range(5):
    for p in range(1, 7):
        for q in range(1, 7):
            for r in range(1, 7):
                high, low = adjacent_cyclic_parts(n, p, q, r)
                direct = (
                    triple(n + 1, p, q) * h(n, r)
                    - triple(n, p, q) * h(n + 1, r)
                )
                require(high + low == direct, "adjacent high/low split")
                adjacent_part_cells += 1

                labels = ((p, q, r), (q, r, p), (r, p, q))
                cyclic_high = sum(adjacent_cyclic_parts(n, *entry)[0] for entry in labels)
                cyclic_low = sum(adjacent_cyclic_parts(n, *entry)[1] for entry in labels)
                scaled_high, scaled_low = global_quarter_scaled(n, p, q, r)
                require(cyclic_high == scaled_high, "global high normalization")
                require(cyclic_low == scaled_low, "global low normalization")
                require(scaled_high > 0, "global high positivity")
                require(4 * scaled_low + scaled_high > 0, "global quarter sign")
                global_quarter_cells += 1

symbolic_stats = symbolic_cyclic_split_certificates()
global_quarter_stats = symbolic_global_quarter_certificates()

print("THM-2830 DISJOINT-CONE CYCLIC-SPLIT THEOREM -- exact referee")
print("status=PROVED ANALYTICALLY; VERIFIED-EXACT CERTIFICATES")
print("triple_tensor_cells=343")
print("pair_product_cells=64")
print(f"polarization_controls={polarization_controls}")
print(f"matching_cells={matching_cells}")
print(f"matching_strict={strict_cells}")
print(f"matching_equalities={equality_cells}")
print(f"normalized_forward_cells={forward_cells}")
print(f"minimum_positive_matching_gap={minimum_positive}")
print(f"pair_factorial_phi={phi}")
print(f"pair_factorial_endpoint={endpoint}")
print(f"pair_factorial_gap={hostile_gap}")
print(f"laguerre_t=1/1000_gap={laguerre_value}")
print(f"cone_cells={cone_cells}")
print(f"cone_equalities={cone_equalities}")
print(f"pascal_identity_cells={pascal_identity_cells}")
print(f"strict_ratio_cells={strict_ratio_cells}")
print(f"newton_bridge_cells={newton_bridge_cells}")
print(f"cyclic_aggregation_cells={cyclic_aggregation_cells}")
print(f"adjacent_part_cells={adjacent_part_cells}")
print(f"global_quarter_cells={global_quarter_cells}")
for key, value in symbolic_stats.items():
    print(f"{key}={value}")
for key, value in global_quarter_stats.items():
    print(f"{key}={value}")
print("stronger_matching_inequality_status=OPEN")
print("all_exact_controls=PASS")

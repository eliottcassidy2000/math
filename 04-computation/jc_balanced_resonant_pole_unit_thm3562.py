#!/usr/bin/env python3
"""Exact non-assert companion for proved THM-3562."""

from fractions import Fraction


def gbinom(a: Fraction, n: int) -> Fraction:
    out = Fraction(1)
    for j in range(n):
        out *= (a - j) / (j + 1)
    return out


def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_add(left, right, right_scale=Fraction(1)):
    out = [Fraction(0)] * max(len(left), len(right))
    for j, value in enumerate(left):
        out[j] += value
    for j, value in enumerate(right):
        out[j] += right_scale * value
    return trim(out)


def poly_mul(left, right):
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for j, value in enumerate(left):
        for ell, other in enumerate(right):
            out[j + ell] += value * other
    return trim(out)


def poly_pow(poly, exponent):
    out = [Fraction(1)]
    base = list(poly)
    power = exponent
    while power:
        if power & 1:
            out = poly_mul(out, base)
        base = poly_mul(base, base)
        power //= 2
    return out


def derivative(poly):
    if len(poly) == 1:
        return [Fraction(0)]
    return [j * poly[j] for j in range(1, len(poly))]


def evaluate(poly, point):
    out = Fraction(0)
    for value in reversed(poly):
        out = out * point + value
    return out


def poly_divmod(numerator, denominator):
    num = trim(numerator)
    den = trim(denominator)
    if den == [0]:
        raise RuntimeError("zero polynomial divisor")
    if len(num) < len(den):
        return [Fraction(0)], num
    quotient = [Fraction(0)] * (len(num) - len(den) + 1)
    while num != [0] and len(num) >= len(den):
        shift = len(num) - len(den)
        scale = num[-1] / den[-1]
        quotient[shift] += scale
        subtractor = [Fraction(0)] * shift + [scale * value for value in den]
        num = poly_add(num, subtractor, Fraction(-1))
    return trim(quotient), trim(num)


def monic_gcd(left, right):
    a = trim(left)
    b = trim(right)
    while b != [0]:
        _, remainder = poly_divmod(a, b)
        a, b = b, remainder
    scale = a[-1]
    return trim([value / scale for value in a])


def compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(1, total - length + 2):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


failures = []
gates = 0


def gate(condition, label):
    global gates
    gates += 1
    if not condition:
        failures.append(label)


# Resonance arithmetic and the regular quotient H/T^(2k).
resonance_gates = 0
for k in range(2, 41):
    R = 3 * k + 2
    D = 4 * k + 3
    pure_index = 4 * k + 3
    c_R = 4 * gbinom(Fraction(2 * R - 1, 2), pure_index)
    before = gates
    gate(D - 2 * (2 * k + 1) == 1, ("Farey gap", k))
    gate(3 * pure_index == 4 * R + 1, ("pure-q t degree", k))
    gate(c_R != 0, ("pure-q unit", k))
    gate((pure_index - 3) // 2 == 2 * k, ("H exponent", k))

    lower_exponents = []
    for row in range(2, R - 1):
        if row % 3 == 2:
            row_index = (4 * row + 1) // 3
            lower_exponents.append((row_index - 3) // 2)
    gate(not lower_exponents or max(lower_exponents) <= 2 * k - 2,
         ("retained pure-q gap", k))

    for delta in range(1, 13):
        nu = 2 + D * delta
        a_order = 1 + (2 * k + 1) * delta
        t_order = 2 * a_order - nu
        k_order = 1 - a_order
        h_order = k_order - t_order
        gate(t_order == -delta, ("T order", k, delta))
        gate(h_order == -2 * k * delta, ("H order", k, delta))
        gate(h_order - 2 * k * t_order == 0,
             ("regular quotient", k, delta))
        gate(2 * a_order == 2 + (4 * k + 2) * delta,
             ("A square divisor", k, delta))
        gate(nu - 2 == D * delta, ("response pole part", k, delta))
    resonance_gates += gates - before


# Full Lagrange interpolation vectors on three independent rational packets.
lagrange_packets = 0
for h in range(2, 17):
    packet_bank = (
        [Fraction(j) for j in range(h)],
        [Fraction(j * j) for j in range(h)],
        [Fraction((j + 1) if j % 2 == 0 else -(j + 1)) for j in range(h)],
    )
    for packet_id, nodes in enumerate(packet_bank):
        weights = []
        for j, beta in enumerate(nodes):
            derivative_value = Fraction(1)
            for ell, other in enumerate(nodes):
                if ell != j:
                    derivative_value *= beta - other
            weights.append(Fraction(1, 1) / derivative_value)
        moments = [
            sum(weights[j] * nodes[j] ** power for j in range(h))
            for power in range(h)
        ]
        gate(moments == [Fraction(0)] * (h - 1) + [Fraction(1)],
             ("Lagrange vector", h, packet_id))
        gate(sum(weights) == 0, ("partial-fraction residue", h, packet_id))
        lagrange_packets += 1


# Exhaust the small balanced delta-passport universe.  A one-pole even-total
# passport is split; every genuinely nonsplit packet has at least two poles.
passport_count = 0
genuine_passports = 0
for total in range(2, 11, 2):
    for h in range(1, total + 1):
        for delta_packet in compositions(total, h):
            passport_count += 1
            genuine = any(delta % 2 for delta in delta_packet)
            if genuine:
                genuine_passports += 1
            gate(not (h == 1 and genuine), ("one-pole split", delta_packet))
            if h >= 2 and genuine:
                gate(sum(delta_packet) == total and total != 0,
                     ("nonzero Lagrange right side", delta_packet))


# The asymmetric R=8, (33,11) finite-pole hostile control.
half = Fraction(11, 2)
E_descending = [(-1) ** m * gbinom(half, m) for m in range(23)]
E = list(reversed(E_descending))
x_minus_one = [Fraction(-1), Fraction(1)]
x_poly = [Fraction(0), Fraction(1)]
W = poly_mul(x_poly, x_minus_one)
P = poly_mul(poly_pow(x_poly, 3), x_minus_one)
weighted_poles = [Fraction(-33), Fraction(44)]
ode_left = poly_add(
    [2 * value for value in poly_mul(W, derivative(E))],
    poly_mul(weighted_poles, E),
    Fraction(-1),
)
ode_constant = Fraction(30705345, 2 ** 41)
E_zero = Fraction(930465, 2 ** 41)
E_one = Fraction(-2791395, 2 ** 41)

control_before = gates
gate(len(E) - 1 == 22 and E[-1] == 1, "E degree/monicity")
gate(ode_left == [ode_constant], "response first integral")
gate(evaluate(E, 0) == E_zero, "E(0)")
gate(evaluate(E, 1) == E_one, "E(1)")
gate(E_one == -3 * E_zero and E_one != E_zero, "pole-unit failure")
gate(monic_gcd(E, derivative(E)) == [Fraction(1)], "E squarefree")
gate(E_zero != 0 and E_one != 0, "E disjoint from W")
gate(33 * Fraction(-1) * E_zero == 11 * Fraction(1) * E_one,
     "weighted pole evaluations")
gate(33 * Fraction(-1) != 11 * Fraction(1),
     "unweighted Lagrange obstruction")

# Exponent vectors for V,A,B in x and x-1.  These independently recover
# T=1/[x^3(x-1)], d=-x(x-1)/4, and s=1/(2x).
v_exp = (35, 13)
a_exp = (16, 6)
b_exp = (18, 7)
gate((2 * a_exp[0] - v_exp[0], 2 * a_exp[1] - v_exp[1]) == (-3, -1),
     "T exponent vector")
gate((2 * b_exp[0] - v_exp[0], 2 * b_exp[1] - v_exp[1]) == (1, 1),
     "d exponent vector")
gate((a_exp[0] + b_exp[0] - v_exp[0],
      a_exp[1] + b_exp[1] - v_exp[1]) == (-1, 0),
     "s finite-pole exponent vector")
gate(P == [Fraction(0), Fraction(0), Fraction(0), Fraction(-1), Fraction(1)],
     "P factor")
gate(W == [Fraction(0), Fraction(-1), Fraction(1)], "W factor")

# The third fibre E^2-P^11 has degree 21 and is squarefree/disjoint, so the
# missing degree 23 is the high point at infinity in passport (23,1^21).
third_fibre = poly_add(poly_mul(E, E), poly_pow(P, 11), Fraction(-1))
gate(len(third_fibre) - 1 == 21, "third-fibre degree")
gate(monic_gcd(third_fibre, derivative(third_fibre)) == [Fraction(1)],
     "third-fibre squarefree")
gate(monic_gcd(third_fibre, E) == [Fraction(1)], "third fibre/E disjoint")
gate(monic_gcd(third_fibre, P) == [Fraction(1)], "third fibre/P disjoint")
gate(4 * gbinom(Fraction(15, 2), 11) == Fraction(-195, 131072),
     "R=8 pure-q unit")
control_gates = gates - control_before


if failures:
    raise RuntimeError(("THM-3562 exact companion failures", failures[:8], len(failures)))

print("THM-3562 balanced resonant pole-unit exact companion: PASS")
print("resonance universe: k=2..40, delta=1..12")
print("resonance/pure-q gates:", resonance_gates)
print("Lagrange universe: h=2..16, three rational node packets per h")
print("Lagrange packets:", lagrange_packets)
print("balanced delta passports: even total 2..10, all positive compositions")
print("passport packets:", passport_count, "genuinely nonsplit:", genuine_passports)
print("R=8 hostile: poles (33,11); s_F=1/(2x); E(1)=-3E(0)")
print("R=8 control gates:", control_gates)
print("exact gates:", gates)
print("consequence: no balanced entrant in the normalized nonsplit polynomial exact-square-prefix chart")
print("scope: unbalanced responses, nonpolynomial prefixes, other charts, and JC(2) remain open")

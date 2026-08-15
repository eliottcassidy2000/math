#!/usr/bin/env python3
"""Exact audit for THM-3383's terminal monomial-cone initial ring.

The theorem is algebraic and all-parameter.  This standard-library companion
checks a normalized hostile grid, the residue/rational-decoding gate, the
cyclic invariant and target-polynomial intersection generators, the complete
normalized C3 atlas on that grid, and the two multistage controls inherited
from the toric-key tower.  Every truth-bearing gate remains active under
``python -O``.
"""

from fractions import Fraction
from hashlib import sha256
import json
from math import gcd


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def clean(poly):
    return {key: value for key, value in poly.items() if value}


def add(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, Fraction()) + value
    return clean(answer)


def scale(poly, scalar):
    return clean({key: scalar * value for key, value in poly.items()})


def mul(left, right):
    answer = {}
    for (lu, ls), lv in left.items():
        for (ru, rs), rv in right.items():
            key = (lu + ru, ls + rs)
            answer[key] = answer.get(key, Fraction()) + lv * rv
    return clean(answer)


def monomial(coefficient=1, u_power=0, s_power=0):
    coefficient = Fraction(coefficient)
    return {} if not coefficient else {(u_power, s_power): coefficient}


ONE = monomial()


def monomial_power(poly, exponent):
    require(len(poly) == 1, ("not monomial", poly, exponent))
    ((u_power, s_power), coefficient), = poly.items()
    return monomial(
        coefficient ** exponent,
        u_power * exponent,
        s_power * exponent,
    )


def power(poly, exponent):
    require(exponent >= 0, ("negative polynomial power", exponent))
    answer = ONE
    base = poly
    current = exponent
    while current:
        if current & 1:
            answer = mul(answer, base)
        base = mul(base, base)
        current //= 2
    return answer


def derivative(poly, axis):
    answer = {}
    coordinate = 0 if axis == "u" else 1
    for powers, coefficient in poly.items():
        exponent = powers[coordinate]
        if not exponent:
            continue
        new_powers = list(powers)
        new_powers[coordinate] -= 1
        answer[tuple(new_powers)] = coefficient * exponent
    return clean(answer)


def wedge(left, right):
    return add(
        mul(derivative(left, "u"), derivative(right, "s")),
        scale(mul(derivative(left, "s"), derivative(right, "u")), -1),
    )


def terminal_packet(E, g, a):
    require(gcd(E, g) == 1, (E, g, a, "unnormalized packet"))
    determinant = g - a * E
    require(determinant, (E, g, a, "zero determinant"))
    coefficient = Fraction(E, determinant)
    R = monomial(1, a, g)
    Z = monomial(coefficient, 1, E)
    M = add(ONE, Z)
    x = monomial_power(R, -1)
    y = mul(M, R)
    T = add(mul(x, y), scale(ONE, -1))
    require(T == Z, (E, g, a, "T"))
    require(wedge(x, y) == monomial(E, 0, E - 1),
            (E, g, a, "Jacobian", wedge(x, y)))

    divisor = gcd(g, E)
    require(divisor == 1, (E, g, divisor, "terminal primitivity"))
    theta = mul(
        monomial_power(Z, g // divisor),
        monomial_power(R, -E // divisor),
    )
    expected_theta = monomial(
        coefficient ** (g // divisor),
        determinant // divisor,
        0,
    )
    require(theta == expected_theta, (E, g, a, "theta"))

    if abs(determinant) != divisor:
        return {
            "rational_decode": False,
            "determinant": determinant,
            "divisor": divisor,
        }

    require(E % determinant == 0 and g % determinant == 0,
            (E, g, a, determinant, "integrality"))
    u_from_ring = mul(
        monomial_power(x, E // determinant),
        monomial_power(scale(T, 1 / coefficient), g // determinant),
    )
    t_from_ring = mul(
        monomial_power(x, -E // determinant),
        monomial_power(
            scale(T, 1 / coefficient), -a * E // determinant
        ),
    )
    require(u_from_ring == monomial(1, 1, 0),
            (E, g, a, "u decoder", u_from_ring))
    require(t_from_ring == monomial(1, 0, E),
            (E, g, a, "t decoder", t_from_ring))

    u_exponents = (E // determinant, g // determinant)
    t_exponents = (-E // determinant, -a * E // determinant)
    v = scale(T, 1 / coefficient)
    L = M
    u_variable = monomial(1, 1, 0)
    t_variable = monomial(1, 0, E)
    X = monomial_power(x, E)
    Y = power(y, E)
    require(mul(X, Y) == power(L, E),
            (E, g, a, "cyclic invariant relation"))

    if determinant > 0:
        require(min(u_exponents) >= 0 and t_exponents[0] < 0,
                (E, g, a, u_exponents, t_exponents))
        polynomial_target = "u"
        reciprocal_target = "t"
        G = mul(
            mul(
                monomial_power(u_variable, g - 1),
                monomial_power(t_variable, g),
            ),
            power(L, E),
        )
        require(G == Y, (E, g, a, "positive generator", G, Y))
        require(mul(u_variable, G)
                == mul(monomial_power(v, g), power(L, E)),
                (E, g, a, "positive hypersurface relation"))
        require(monomial_power(x, E)
                == mul(u_variable, monomial_power(v, -g)),
                (E, g, a, "positive target-field generator"))
    else:
        require(max(u_exponents) < 0 and min(t_exponents) >= 0,
                (E, g, a, u_exponents, t_exponents))
        polynomial_target = "t"
        reciprocal_target = "u"
        ae = a * E
        G = mul(
            mul(
                monomial_power(u_variable, ae),
                monomial_power(t_variable, ae - 1),
            ),
            power(L, E),
        )
        require(G == Y, (E, g, a, "negative generator", G, Y))
        require(mul(t_variable, G)
                == mul(monomial_power(v, ae), power(L, E)),
                (E, g, a, "negative hypersurface relation"))
        require(monomial_power(x, E)
                == mul(t_variable, monomial_power(v, -ae)),
                (E, g, a, "negative target-field generator"))

    return {
        "rational_decode": True,
        "determinant": determinant,
        "divisor": divisor,
        "polynomial_target": polynomial_target,
        "reciprocal_target": reciprocal_target,
        "u_exponents": u_exponents,
        "t_exponents": t_exponents,
    }


def multistage_controls():
    u = monomial(1, 1, 0)
    s2 = monomial(1, 0, 2)
    s3 = monomial(1, 0, 3)

    # THM-3080's 7=4+3 packet.  The triangular coordinate Y=y-x is the
    # terminal (g,a)=(2,1) monomial-cone coordinate.
    R = mul(u, s2)
    M1 = add(ONE, scale(mul(u, s3), -3))
    M0 = add(ONE, mul(monomial_power(R, 2), M1))
    x = monomial_power(R, -1)
    y = mul(M0, x)
    Y = add(y, scale(x, -1))
    require(Y == mul(M1, R), ("4+3 shear", Y, mul(M1, R)))
    standard_43 = terminal_packet(3, 2, 1)
    require(standard_43["polynomial_target"] == "t", standard_43)

    # THM-3080's 11=4+4+3 packet.  Y=y-x-1 is the terminal
    # (g,a)=(4,1) coordinate.
    R = monomial(1, 1, 4)
    M2 = add(ONE, scale(mul(u, s3), 3))
    M1 = add(ONE, mul(R, M2))
    M0 = add(ONE, mul(R, M1))
    x = monomial_power(R, -1)
    y = mul(M0, x)
    Y = add(add(y, scale(x, -1)), scale(ONE, -1))
    require(Y == mul(M2, R), ("4+4+3 shear", Y, mul(M2, R)))
    standard_443 = terminal_packet(3, 4, 1)
    require(standard_443["polynomial_target"] == "u", standard_443)

    return {
        "4+3": {
            "terminal": (2, 1),
            "polynomial_target": standard_43["polynomial_target"],
        },
        "4+4+3": {
            "terminal": (4, 1),
            "polynomial_target": standard_443["polynomial_target"],
        },
    }


def c3_expected(g):
    residue = g % 3
    if residue == 1:
        return ((g - 1) // 3, 1, "u"),
    if residue == 2:
        return ((g + 1) // 3, -1, "t"),
    return ()


def main():
    checked = 0
    rational_decode = []
    for E in range(1, 25):
        for g in range(1, 49):
            if gcd(E, g) != 1:
                continue
            for a in range(0, 49):
                if g == a * E:
                    continue
                packet = terminal_packet(E, g, a)
                checked += 1
                if packet["rational_decode"]:
                    rational_decode.append((E, g, a, packet))

    c3_atlas = []
    for g in range(1, 25):
        actual = []
        for E, current_g, a, packet in rational_decode:
            if E == 3 and current_g == g:
                actual.append((a, packet["determinant"],
                               packet["polynomial_target"]))
        actual = tuple(sorted(actual))
        expected = tuple(sorted(c3_expected(g)))
        require(actual == expected, ("C3 atlas", g, actual, expected))
        c3_atlas.append((g, actual))

    controls = multistage_controls()
    semantic = {
        "checked": checked,
        "rational_decode_cells": len(rational_decode),
        "c3_atlas": c3_atlas,
        "controls": controls,
    }
    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    print("THM-3383 TERMINAL MONOMIAL-CONE INITIAL-RING AUDIT")
    print(f"normalized_parameter_grid=PASS checked={checked} "
          f"rational_decode_cells={len(rational_decode)}")
    print("initial_ring=PASS cyclic_invariant_hypersurface; "
          "positive_and_negative_target_polynomial_slices")
    print("fork=PASS n>0 gives polynomial u / reciprocal t; "
          "n<0 gives reciprocal u / polynomial t")
    print("normalized_C3_atlas=PASS " + ";".join(
        f"g{g}:" + (",".join(
            f"a{a}/n{det}/{target}" for a, det, target in rows
        ) or "excluded_nonprimitive")
        for g, rows in c3_atlas
    ))
    print("multistage_controls=PASS " + json.dumps(controls, sort_keys=True))
    print(f"semantic_sha256={digest}")


if __name__ == "__main__":
    main()

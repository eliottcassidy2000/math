#!/usr/bin/env python3
"""Bounded exact root-sign probe plus one unbounded symbolic family.

No repository theorem or producer import. Requires SymPy. The fourteen named
controls are a declared probe, not a census; thirteen have >=3 first channels.
All gates survive python -O. See the matching result note for scope.
"""
from __future__ import annotations

import hashlib
import json
from math import factorial, gcd

import sympy as S


T, G, Z = S.symbols("tau g s")
GATES = 0


def require(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def raw_channels(a, b, c, mass):
    """Solve the original charge equation, scanning the negative count."""
    answer = []
    for x in range(mass + 1):
        numerator = a * x - b * (mass - x)
        if numerator % (c - b):
            continue
        z = numerator // (c - b)
        y = mass - x - z
        if min(y, z) >= 0:
            weight = factorial(mass) // (factorial(x) * factorial(y) * factorial(z))
            answer.append(((x, y, z), weight))
    return sorted(answer, key=lambda item: item[0][2])


def native_rows(a, b, c, final_mass):
    """Literal multiplication of (u^-a+u^b+v*u^c), preserving v degree."""
    current = {(0, 0): 1}
    rows = []
    for _ in range(final_mass):
        following = {}
        for (charge, z), value in current.items():
            for delta, dz in ((-a, 0), (b, 0), (c, 1)):
                key = charge + delta, z + dz
                following[key] = following.get(key, 0) + value
        current = following
        rows.append({z: value for (charge, z), value in current.items() if charge == 0})
    return rows


def interval_polynomial(poly, left, right):
    """Exact Horner enclosure; a nonzero enclosure certifies the root sign."""
    lo = hi = S.Rational(0)
    for coefficient in poly.all_coeffs():
        products = (lo * left, lo * right, hi * left, hi * right)
        lo, hi = min(products) + coefficient, max(products) + coefficient
    return lo, hi


def certified_root_signs(p, remainder):
    for precision in (4, 8, 16, 32, 64, 128):
        intervals = p.intervals(eps=S.Rational(1, 10) ** precision)
        require(len(intervals) == p.degree(), "all first roots real and distinct")
        data = []
        for (left, right), multiplicity in intervals:
            require(multiplicity == 1 and right < 0, "simple negative root")
            require(p.count_roots(left, right) == 1, "exact Sturm isolation")
            lo, hi = interval_polynomial(remainder, left, right)
            if lo > 0:
                sign = 1
            elif hi < 0:
                sign = -1
            else:
                break
            data.append((str(left), str(right), str(lo), str(hi), sign))
        else:
            return [entry[-1] for entry in data], data, precision
    raise RuntimeError("bounded interval refinement failed")


def falling(value, length):
    return S.prod(value - j for j in range(length))


def symbolic_family():
    # For g integer >=5 these are the exact mass-g and mass-2g rows.
    # Only odd g makes g the first support return for primitive charges.
    original_p = S.Poly(sum(falling(G, 4-j) * T**j /
                           (factorial(4-2*j) * factorial(j)) for j in range(3)), T)
    q = S.Poly(sum(falling(2*G, 8-j) * T**j /
                  (factorial(8-2*j) * factorial(j)) for j in range(5)), T)
    p = S.Poly(12*T**2 + 12*(G-2)*T + (G-2)*(G-3), T)
    require(S.cancel(original_p.as_expr() - G*(G-1)*p.as_expr()/24) == 0,
            "positive first-row content")
    constant = G*(G-3)*(G-2)*(G-1)*(2*G-3)*(2*G-1) * \
        (503*G**2-1723*G+1470)/30240
    slope = G*(G-2)*(G-1)*(2*G-3)*(2*G-1)*(3*G-5)*(11*G-18)/180
    remainder = q.rem(p)
    require(S.cancel(remainder.nth(0)-constant) == 0 and
            S.cancel(remainder.nth(1)-slope) == 0 and remainder.degree() == 1,
            "complete symbolic Euclidean remainder")
    rho = -(G-3)*(503*G**2-1723*G+1470)/(168*(3*G-5)*(11*G-18))
    require(S.cancel(constant + slope*rho) == 0, "exact remainder zero")
    h = 5141*G**4-31762*G**3+69537*G**2-62820*G+18900
    j = 2269*G**3-11468*G**2+19233*G-10710
    require(S.cancel(p.eval(rho) - (G-3)*(5*G-7)*h /
                     (2352*(3*G-5)**2*(11*G-18)**2)) == 0,
            "P at remainder zero")
    require(S.cancel(rho+(G-2)/2-j/(168*(3*G-5)*(11*G-18))) == 0,
            "remainder zero to right of first vertex")
    require(S.cancel(S.discriminant(p.as_expr(), T)-48*(G-2)*(2*G-3)) == 0,
            "two distinct negative first roots for g>=5")
    shifted = {}
    for name, expression in (("H", h), ("J", j),
                             ("constant_factor", 503*G**2-1723*G+1470)):
        shifted[name] = [str(value) for value in
                         S.Poly(S.expand(expression.subs(G, Z+5)), Z).all_coeffs()]
        require(all(S.Rational(value) > 0 for value in shifted[name]),
                "strict shifted coefficient positivity: " + name)
    # Independent fixed-support channel recovery at three primitive parameters.
    for g in (5, 7, 9):
        a, b, c = 4, g-4, 2*g-4
        for mass, expression in ((g, original_p), (2*g, q)):
            raw = raw_channels(a, b, c, mass)
            require([value for _, value in raw] ==
                    list(reversed(S.Poly(expression.as_expr().subs(G, g), T).all_coeffs())),
                    "symbolic versus original charge equation")
    print("symbolic family: (-4,g-4,2g-4), odd integers g>=5: Q negative at both P roots")
    print("factorial extension: real g>=5; first-return interpretation requires odd integer g")
    print("shifted positive factors:", json.dumps(shifted, sort_keys=True))
    return {"p": str(p.as_expr()), "constant": str(S.factor(constant)),
            "slope": str(S.factor(slope)), "rho": str(rho), "shifted": shifted}


def main():
    # Explicit named supports only; no support-height/gcd search is performed.
    controls = [(4, 1, 6), (4, 3, 10), (4, 5, 14), (5, 1, 7),
                (6, 1, 8), (6, 5, 16), (7, 1, 9), (7, 2, 11),
                (8, 1, 10), (8, 3, 14), (9, 1, 11), (10, 1, 12),
                (15, 1, 9), (13, 1, 8)]
    manifest = []
    for a, b, c in controls:
        require(gcd(a, gcd(b, c)) == 1 and 0 < b < c, "primitive positive charges")
        g = gcd(a+b, a+c)
        step = (a+b)//g
        native = native_rows(a, b, c, 2*g)
        require(all(not row for row in native[:g-1]), "first support mass")
        polynomials, anchors = [], []
        for mass in (g, 2*g):
            raw = raw_channels(a, b, c, mass)
            require(bool(raw), "nonempty complete row")
            require(native[mass-1] == {channel[2]: value for channel, value in raw},
                    "literal Laurent multiplication equals charge-equation weights")
            require([channel[2] for channel, _ in raw] ==
                    [raw[0][0][2]+step*j for j in range(len(raw))], "full affine index")
            polynomials.append(S.Poly(sum(value*T**j for j, (_, value) in enumerate(raw)), T))
            anchors.append(raw[0][0])
        p, q = polynomials
        remainder = q.rem(p)
        signs, intervals, precision = certified_root_signs(p, remainder)
        epsilon = 2*anchors[0][2]//step
        other_step = (a+c)//g
        affine_step = (other_step-step, -other_step, step)
        require(anchors[1] == tuple(2*v-epsilon*d for v, d in zip(anchors[0], affine_step)),
                "second anchor retains the lower carry")
        expected = (-1)**(epsilon+1)
        require(all(sign == expected for sign in signs), "one strict sign on first roots")
        require(all((-1)**epsilon*sign == -1 for sign in signs), "negative first-anchor normalization")
        require(S.gcd(p, q).degree() == 0, "independent Euclidean coprimality")
        print(f"support={(-a,b,c)} first_mass={g} channels={p.degree()+1}/{q.degree()+1} "
              f"Q_root_signs={signs} lower_carry={epsilon} interval_epsilon=10^-{precision}")
        if (a, b, c) == (4, 1, 6):
            require(p.as_expr() == 5+30*T+10*T*T, "specified first polynomial")
            require(q.as_expr() == 45+840*T+3150*T*T+2520*T**3+210*T**4,
                    "specified complete second polynomial")
            require(remainder.as_expr() == S.Rational(2715,2)+7770*T, "literal remainder")
            rho = -S.Rational(181,1036)
            require(p.eval(rho) == S.Rational(34305,536648), "positive quadratic obstruction test")
            print("specified Q mod P =", remainder.as_expr(),
                  "; rho =", rho, "; P(rho) =", p.eval(rho))
        manifest.append({"support": [-a,b,c], "first_mass": g, "lower_carry": epsilon,
                         "P": [str(v) for v in reversed(p.all_coeffs())],
                         "Q": [str(v) for v in reversed(q.all_coeffs())],
                         "remainder": str(remainder.as_expr()), "intervals": intervals})
    family = symbolic_family()
    # This is deliberately not asserted to be an actual factorial-row pair.
    toy_p = S.Poly((T+1)*(T+3), T)
    toy_q = S.Poly((2*T+1)*(T+2)*(T+4)*(T+5), T)
    require(toy_q.rem(toy_p).as_expr() == -23-11*T, "logical hostile remainder")
    require(toy_q.eval(-3) == 10 and toy_q.eval(-1) == -12, "logical hostile opposite signs")
    require(toy_p.eval(-S.Rational(23,11)) == -S.Rational(120,121),
            "logical hostile strict root separation")
    require(S.gcd(toy_p, toy_q).degree() == 0, "coprimality does not imply one root sign")
    print("logical hostile only: P=(tau+1)(tau+3), Q=(2tau+1)(tau+2)(tau+4)(tau+5)")
    print("Q mod P=-23-11tau; Q(-3)=10; Q(-1)=-12; P(-23/11)=-120/121")
    print("scope: 13 named >=3-channel supports PASS; one two-carry normalization control PASS")
    print("general actual higher-channel root-sign and two-rung problems remain OPEN")
    encoded = json.dumps({"controls": manifest, "family": family}, sort_keys=True,
                         separators=(",", ":")).encode()
    print("semantic_sha256", hashlib.sha256(encoded).hexdigest())
    print("explicit_gates", GATES)


if __name__ == "__main__":
    main()

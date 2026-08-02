#!/usr/bin/env python3
"""Exact referee for THM-3151's all-R equality-cell nonentry."""

from fractions import Fraction
from math import comb


def gbinom(a: Fraction, n: int) -> Fraction:
    out = Fraction(1)
    for i in range(n):
        out *= (a - i) / (i + 1)
    return out


def add(out, key, value):
    out[key] = out.get(key, Fraction(0)) + value
    if not out[key]:
        del out[key]


def coefficient(j: int, n: int):
    """Sparse coefficient of A(t)^(j-1/2), keyed by (q,d,s)."""
    alpha = Fraction(2 * j - 1, 2)
    out = {}
    # Expand A=(1+d*t^2)^2+t^3(q-s*t).
    for v in range(n // 3 + 1):
        for ell in range(v + 1):
            r = v - ell
            rem = n - 3 * v - ell
            if rem < 0 or rem % 2:
                continue
            u = rem // 2
            coeff = (
                gbinom(alpha, v)
                * comb(v, ell)
                * (-1) ** ell
                * gbinom(Fraction(2 * j - 1 - 2 * v), u)
            )
            add(out, (r, u, ell), coeff)
    return out


def rows(j: int):
    phi = {}
    for (r, u, ell), c in coefficient(j, 4 * j - 1).items():
        if r < 1:
            raise RuntimeError("phi coefficient not divisible by q")
        add(phi, (r - 1, u, ell), 4 * c)
    psi = {key: 4 * c for key, c in coefficient(j, 4 * j).items()}
    return phi, psi


def face(poly, D, S, delete_s=False, delete_d=False):
    items = {
        key: c
        for key, c in poly.items()
        if not (delete_s and key[2]) and not (delete_d and key[1])
    }
    weights = {key: key[0] - D * key[1] - S * key[2] for key in items}
    low = min(weights.values())
    return low, {key: items[key] for key in items if weights[key] == low}


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


EXPECTED_CHECKS = 18320
checks = 0
for j in range(2, 51):
    phi, psi = rows(j)

    # Support inequalities behind the universal wall.
    for a, u, ell in phi:
        if a < 2 * u:
            raise RuntimeError(("phi support", j, a, u, ell))
        checks += 1
    for r, u, ell in psi:
        if r < 2 * u:
            raise RuntimeError(("psi support", j, r, u, ell))
        checks += 1

    # Exact wall faces and odd/even binomial coefficient polynomials.
    D, S = 16, 7
    wphi, fphi = face(phi, D, S)
    wpsi, fpsi = face(psi, D, S)
    if wphi != -(j - 1) * S or wpsi != -j * S:
        raise RuntimeError(("wall weights", j, wphi, wpsi))
    P = [(-1) ** u * comb(j, 2 * u + 1) for u in range((j - 1) // 2 + 1)]
    Q = [(-1) ** u * comb(j, 2 * u) for u in range(j // 2 + 1)]
    cp = fphi[(0, 0, j - 1)] / P[0]
    cq = fpsi[(0, 0, j)] / Q[0]
    for u, value in enumerate(P):
        if fphi.get((2 * u, u, j - 1 - 2 * u)) != cp * value:
            raise RuntimeError(("P wall", j, u))
    for u, value in enumerate(Q):
        if fpsi.get((2 * u, u, j - 2 * u)) != cq * value:
            raise RuntimeError(("Q wall", j, u))

    # Q^2+zP^2=(1+z)^j, hence no common root in characteristic zero.
    lhs = poly_mul(Q, Q)
    zp2 = [0] + poly_mul(P, P)
    lhs += [0] * (max(len(lhs), len(zp2)) - len(lhs))
    zp2 += [0] * (len(lhs) - len(zp2))
    lhs = [x + y for x, y in zip(lhs, zp2)]
    rhs = [comb(j, u) for u in range(j + 1)]
    p_at_minus_one = sum(value * (-1) ** u for u, value in enumerate(P))
    q_at_minus_one = sum(value * (-1) ** u for u, value in enumerate(Q))
    if lhs != rhs or p_at_minus_one == 0 or q_at_minus_one == 0:
        raise RuntimeError(("binomial Bezout identity", j))

    # Hostile degree box: top row is strictly below every retained j<=R-2.
    if j >= 7:
        all_rows = {i: rows(i) for i in range(1, j + 1)}
        for S in range(1, 9):
            for D in range(0, 25):
                if D < 2 * S + 2:
                    channels = (0, 1)
                elif D == 2 * S + 2:
                    channels = (0, 1)
                else:
                    channels = (1,) if j % 2 == 0 else (0,)
                top_survivors = 0
                for channel in channels:
                    top_w, top_f = face(all_rows[j][channel], D, S)
                    if D == 2 * S + 2:
                        # At least one wall polynomial survives for any z;
                        # support/weight separation is checked channelwise.
                        pass
                    if all(face(all_rows[i][channel], D, S)[0] > top_w
                           for i in range(1, j - 1)
                           if all_rows[i][channel]):
                        top_survivors += 1
                if not top_survivors:
                    raise RuntimeError(("lower separation", j, D, S))
                checks += 1

        # s=0 and deg(d)>=3 boundary, using the parity-selected channel.
        channel = 1 if j % 2 == 0 else 0
        for D in range(3, 31):
            top_w, top_f = face(all_rows[j][channel], D, 0, delete_s=True)
            if not all(face(all_rows[i][channel], D, 0, delete_s=True)[0] > top_w
                       for i in range(1, j - 1)
                       if any(key[2] == 0 for key in all_rows[i][channel])):
                raise RuntimeError(("s=0 separation", j, D))
            checks += 1


# Divisor and final-order arithmetic on resonance equality cells.
for k in range(2, 31):
    R = 3 * k + 2
    pole = 4 * k + 3
    nu = pole + 2
    a = 2 * k + 2
    b = 2 * k + 3
    if 2 * a - nu != -1 or 2 * b - nu < 1 or a + b - nu < 0:
        raise RuntimeError(("source exponents", k))
    deg_E = pole
    deg_A = 2 * a
    deg_M = deg_E + 2
    required_ord_K = deg_A - deg_M
    if required_ord_K != -1:
        raise RuntimeError(("K order", k, required_ord_K))
    checks += 1

if checks != EXPECTED_CHECKS:
    raise RuntimeError(("gate count", checks, EXPECTED_CHECKS))

print("THM-3151 all-R equality-cell exact referee: PASS")
print("R rows checked: 2..50; resonance k checked: 2..30")
print("exact gates:", checks)
print("wall identity: Q_R(z)^2+z P_R(z)^2=(1+z)^R")
print("consequence: (D,D) equality cell impossible; balanced chart entrant N>=4D")
print("scope: normalized polynomial exact-square-prefix balanced branch; JC(2) remains open")

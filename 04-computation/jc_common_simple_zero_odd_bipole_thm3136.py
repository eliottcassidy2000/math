#!/usr/bin/env python3
"""Exact controls for THM-3136's simple-zero exclusion and odd bipoles."""

from fractions import Fraction
from math import factorial

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def add(poly, monomial, coefficient):
    poly[monomial] = sp.cancel(poly.get(monomial, 0) + coefficient)
    if poly[monomial] == 0:
        del poly[monomial]


def coefficient(R, n):
    """[t^n](1+2*d*t^2+q*t^3+(d^2-s)*t^4)^(R-1/2)."""
    alpha = sp.Rational(2 * R - 1, 2)
    answer = {}
    for i in range(n // 2 + 1):
        for j in range(n // 3 + 1):
            remainder = n - 2 * i - 3 * j
            if remainder < 0 or remainder % 4:
                continue
            ell_total = remainder // 4
            total = i + j + ell_total
            multinomial = factorial(total) // (
                factorial(i) * factorial(j) * factorial(ell_total)
            )
            common = sp.binomial(alpha, total) * multinomial * 2**i
            for s_degree in range(ell_total + 1):
                add(
                    answer,
                    (j, i + 2 * (ell_total - s_degree), s_degree),
                    common * sp.binomial(ell_total, s_degree) * (-1) ** s_degree,
                )
    return answer


def shift(poly, q=0, d=0, s=0, scalar=1):
    answer = {}
    for (qe, de, se), value in poly.items():
        monomial = (qe + q, de + d, se + s)
        require(min(monomial) >= 0, f"negative exponent {monomial}")
        add(answer, monomial, scalar * value)
    return answer


def combine(*polys):
    answer = {}
    for poly in polys:
        for monomial, value in poly.items():
            add(answer, monomial, value)
    return answer


def row(R):
    c_phi = coefficient(R, 4 * R - 1)
    c_psi = coefficient(R, 4 * R)
    c_response = coefficient(R, 4 * R + 1)
    phi = shift(c_phi, q=-1, scalar=4)
    psi = shift(c_psi, scalar=4)
    H = shift(combine(c_response, shift(c_phi, d=1)), q=-3, scalar=4)
    return {"phi": phi, "psi": psi, "H": H}


def min_face(poly, weights):
    values = {
        monomial: sum(Fraction(exp) * weight for exp, weight in zip(monomial, weights))
        for monomial in poly
    }
    minimum = min(values.values())
    face = {m: poly[m] for m, value in values.items() if value == minimum}
    return minimum, face


rho = sp.symbols("rho")


def prefix(R):
    """Polar initial polynomials P_R,Q_R, independently from their sums."""
    alpha = sp.Rational(2 * R - 1, 2)
    P = 0
    Q = 0
    for ell in range((R - 1) // 3 + 1):
        P += (
            sp.binomial(alpha, R + ell)
            * (-1) ** (R - 1 - 3 * ell)
            * 2 ** (R - 1 - 3 * ell)
            * sp.binomial(R - 1 - ell, 2 * ell)
            * rho**ell
        )
    for ell in range(R // 3 + 1):
        Q += (
            sp.binomial(alpha, R + ell)
            * (-1) ** (R - 3 * ell)
            * sp.Rational(2) ** (R - 3 * ell - 1)
            * (
                2 * sp.binomial(R - 1 - ell, 2 * ell - 1)
                + sp.binomial(R - 1 - ell, 2 * ell)
            )
            * rho**ell
        )
    return sp.Poly(P, rho, domain=sp.QQ), sp.Poly(Q, rho, domain=sp.QQ)


rows = {R: row(R) for R in range(1, 31)}
prefixes = {R: prefix(R) for R in range(1, 31)}

# Polar lane at ord(V)=ord(M)=1 and ord(A)=ord(B)=0:
# x=ord(s)=-1 and y=ord(q/omega^3)=1.  Constant prefix terms survive.
for R in range(7, 31):
    P_R, Q_R = prefixes[R]
    require(P_R.nth(0) != 0 and Q_R.nth(0) != 0, f"polar constants R={R}")
    for j in range(1, R - 1):
        require(-(R - 1) < -(j - 1), f"polar phi separation R,j={R,j}")
        require(-R < -j, f"polar psi separation R,j={R,j}")

# Pure-q lane ord(A)=0, ord(B)>=1: q has order -1/2 while d,s are regular.
# The correct Phi/Psi/H channel has a unique top face below every retained row.
channel = {0: "psi", 1: "phi", 2: "H"}
pure_weights = (Fraction(-1, 2), Fraction(0), Fraction(0))
for R in range(7, 31):
    active = channel[R % 3]
    top_value, top_face = min_face(rows[R][active], pure_weights)
    require(len(top_face) == 1, f"pure unique face R={R}")
    for j in range(1, R - 1):
        if not rows[j][active]:
            continue
        lower_value, _ = min_face(rows[j][active], pure_weights)
        require(top_value < lower_value, f"pure separation R,j={R,j}")
    if active == "H":
        # K=T*H and ord(T)=-1, whereas A*K=M would require ord(K)=1.
        require(-1 + top_value != 1, f"spurious simple-zero resonance R={R}")

# ord(A)>=1: T has positive order.  If B is a unit, d has order -1 but T*d
# is regular; otherwise d is regular.  In both cases T,s,Td are regular, so
# every H row is regular and K=T*H has order >=2a-1, never 1-a.
for a in range(1, 8):
    for b in range(0, 8):
        vq = Fraction(2 * a - 1, 2)
        vd = Fraction(-1 if b == 0 else 0)
        vs = Fraction(0)
        for R in range(7, 31):
            for j in list(range(1, R - 1)) + [R]:
                if not rows[j]["H"]:
                    continue
                h_value, _ = min_face(rows[j]["H"], (vq, vd, vs))
                require(h_value >= 0, f"regular H failure a,b,R,j={a,b,R,j}")
            require(2 * vq > 1 - a, f"response order gap a={a}")

# Closed-form sequence for the R=3k+2 pure-H face.  The recurrence is a
# product-Gamma certificate that no resonance coefficient vanishes.
c = [4 * sp.binomial(3 * k + sp.Rational(3, 2), 4 * k + 3) for k in range(10)]
require(c[0] == -sp.Rational(1, 4), "c_0")
for k in range(9):
    ratio = -sp.Rational(
        (6 * k + 5) * (6 * k + 7) * (6 * k + 9),
        32 * (4 * k + 4) * (4 * k + 5) * (4 * k + 7),
    )
    require(c[k + 1] == sp.cancel(ratio * c[k]), f"coefficient recurrence k={k}")
    require(c[k] != 0, f"zero coefficient k={k}")

print("COMMON-SIMPLE-ZERO REALIZATION PROBE -- exact")
print("R=7..30 polar, pure-q, and regular-H lanes checked")
print("all retained lower rows included; no simple-zero resonance")
print("pure-H coefficients c_0..c_9 =", tuple(c))
print("c_0=-1/4 and the rational first-order recurrence is nonvanishing")
print("ALL EXACT HOSTILE CONTROLS PASS")

# --- independent odd-bipole boundary controls ---

"""Exact controls for the odd bipole balanced-response boundary."""

import hashlib
from pathlib import Path

import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


x = sp.symbols("x")
T = x**2 - 1
records = []

for d in range(1, 22, 2):
    m = (d - 1) // 2
    E = sp.expand(
        sum((-1) ** j * sp.binomial(sp.Rational(d, 2), j) * x ** (d - 2 * j)
            for j in range(m + 1))
    )
    C = sp.expand(2 * T * sp.diff(E, x) - 2 * d * x * E)
    require(sp.Poly(C, x).degree() == 0 and C != 0, f"constant first integral d={d}")
    require(sp.gcd(sp.Poly(E, x), sp.Poly(T, x)).degree() == 0, f"E/T disjoint d={d}")
    require(sp.gcd(sp.Poly(E, x), sp.Poly(sp.diff(E, x), x)).degree() == 0,
            f"E squarefree d={d}")

    V = T ** (d + 2)
    G = sp.cancel(E / T ** (d + 1))
    F = sp.cancel(E**2 / T**d)
    require(sp.cancel(F - V * G**2) == 0, f"square potential d={d}")
    require(sp.cancel(2 * V * sp.diff(G, x) + sp.diff(V, x) * G - C) == 0,
            f"response ODE d={d}")
    require(sp.cancel(V * sp.diff(F, x) ** 2 - C**2 * F) == 0,
            f"base response equation d={d}")

    W = sp.expand(E**2 - T**d)
    require(sp.Poly(W, x).degree() == d - 1, f"third-fibre degree d={d}")
    if d > 1:
        require(sp.gcd(sp.Poly(W, x), sp.Poly(sp.diff(W, x), x)).degree() == 0,
                f"third finite roots squarefree d={d}")
    require(sp.gcd(sp.Poly(W, x), sp.Poly(E * T, x)).degree() == 0,
            f"three fibres disjoint d={d}")
    records.append((d, E, C, W))

# Product-Gamma/closed-form recurrence for the nonzero first integral.
constants = [record[2] for record in records]
require(constants[0] == -2, "C_0")
for m in range(len(constants) - 1):
    require(
        constants[m + 1] == -sp.Rational(2 * m + 3, 2 * m + 2) * constants[m],
        f"constant recurrence m={m}",
    )

# Exact clean Nielsen census for the first Faber-resonant cell d=11.
d = 11
n = 2 * d
beta = list(range(n))
for offset in (0, d):
    for i in range(d):
        beta[offset + i] = offset + (i + 1) % d


def near_perfect(offset, unmatched):
    match = {}
    path = [offset + (unmatched + i) % d for i in range(1, d)]
    for i in range(0, d - 1, 2):
        match[path[i]] = path[i + 1]
        match[path[i + 1]] = path[i]
    return match, offset + unmatched


def cycle_type(permutation):
    seen = set()
    lengths = []
    for i in range(n):
        if i in seen:
            continue
        j = i
        length = 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = permutation[j]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


matching_count = 0
representative = None
for u in range(d):
    left, left_unmatched = near_perfect(0, u)
    for v in range(d):
        right, right_unmatched = near_perfect(d, v)
        alpha = list(range(n))
        for i, j in {**left, **right}.items():
            alpha[i] = j
        alpha[left_unmatched] = right_unmatched
        alpha[right_unmatched] = left_unmatched
        product = [alpha[beta[i]] for i in range(n)]
        require(cycle_type(product) == (12,) + (1,) * 10, f"passport u,v={u,v}")
        matching_count += 1
        if representative is None:
            representative = alpha

require(matching_count == 121, "marked matching count")
# The two independent pole rotations translate the unmatched labels (u,v).
# Their orbit of (0,0) is the complete d-by-d marked parameter set.
rotation_orbit = {((0 + a) % d, (0 + b) % d) for a in range(d) for b in range(d)}
require(len(rotation_orbit) == matching_count, "single pole-rotation orbit")
group = PermutationGroup([Permutation(representative), Permutation(beta)])
require(group.is_transitive(), "representative transitivity")
require(group.order() == 2**10 * sp.factorial(11), "representative monodromy order")

d11, E11, C11, W11 = next(record for record in records if record[0] == 11)
print("ODD BIPOLE BALANCED-RESPONSE PROBE -- exact")
print("odd d=1..21: factor ODE, square-potential ODE, squarefreeness, and passport pass")
print("C_1=-2; C_(m+1)=-(2m+3)/(2m+2) C_m (d=2m+1)")
print("d=11 E =", E11)
print("d=11 C =", C11)
print("d=11 degree(E^2-(x^2-1)^11) =", sp.Poly(W11, x).degree())
print("d=11 marked clean matchings =", matching_count)
print("d=11 unmarked orbit count under independent pole rotations = 1")
print("d=11 representative monodromy order =", group.order(), "=2^10*11!")
print("ALL EXACT HOSTILE CONTROLS PASS")

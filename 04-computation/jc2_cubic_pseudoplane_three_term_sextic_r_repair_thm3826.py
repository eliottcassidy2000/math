#!/usr/bin/env python3
"""Exact companion for THM-3826's at-most-three-term sextic closure."""

from __future__ import annotations

import ast
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


r, z, e, u = sp.symbols("r z e u")
b0, b1, b2, b3, b4, b5, b6 = sp.symbols("b0:7")
G, D = sp.symbols("G D")
coefficients = (b0, b1, b2, b3, b4, b5, b6)
g = sum(coefficient * e**degree for degree, coefficient in enumerate(coefficients))
gp = sp.diff(g, e)

# THM-3820's coefficient-free residual, independently replayed here.
H_universal = (
    -729 * D**4 * e**2 + 1458 * D**3 * G * e + 8748 * D**3 * e**4
    + 2125764 * D**2 * G**3 * e**4 - 729 * D**2 * G**2
    - 17496 * D**2 * G * e**3 - 34992 * D**2 * e**6
    - 4251528 * D * G**4 * e**3 - 8503056 * D * G**3 * e**6
    + 11664 * D * G**2 * e**2 + 52488 * D * G * e**5
    + 46656 * D * e**8 - 2 * D * e + 2125764 * G**5 * e**2
    + 8503056 * G**4 * e**5 + 8503056 * G**3 * e**8
    - 2916 * G**3 * e - 17496 * G**2 * e**4 - 23328 * G * e**7 + G
)
H = sp.Poly(sp.expand(H_universal.subs({G: g, D: gp})), e)
gate(H.degree() == 32, "sextic residual degree thirty-two")
gate(H.LC() == 53144100 * b6**5, "sextic residual leading coefficient")
gate(len(H.terms()) == 33, "all sextic residual e-layers present generically")

# Exact logarithmic quotient.  The integer numerator is primitive and has no
# b6 factor, so the displayed b6^5 localization is not hiding a support seam.
denominator = 15625 * b6**5
q_clear = (
    256250 * b0 * b6**5 - 25000 * b1 * b5 * b6**4
    + 265625 * b1 * b6**5 * e - 40000 * b2 * b4 * b6**4
    + 21000 * b2 * b5**2 * b6**3 - 25000 * b2 * b5 * b6**4 * e
    + 287500 * b2 * b6**5 * e**2 - 22500 * b3**2 * b6**4
    + 50500 * b3 * b4 * b5 * b6**3 - 37500 * b3 * b4 * b6**4 * e
    - 17600 * b3 * b5**3 * b6**2 + 20000 * b3 * b5**2 * b6**3 * e
    - 22500 * b3 * b5 * b6**4 * e**2 + 321875 * b3 * b6**5 * e**3
    + 7500 * b3 * b6**4 + 9000 * b4**3 * b6**3
    - 28200 * b4**2 * b5**2 * b6**2 + 22500 * b4**2 * b5 * b6**3 * e
    - 15000 * b4**2 * b6**4 * e**2 + 14720 * b4 * b5**4 * b6
    - 16000 * b4 * b5**3 * b6**2 * e
    + 17000 * b4 * b5**2 * b6**3 * e**2
    - 17500 * b4 * b5 * b6**4 * e**3 - 5500 * b4 * b5 * b6**3
    + 368750 * b4 * b6**5 * e**4 - 2048 * b5**6
    + 2560 * b5**5 * b6 * e - 3200 * b5**4 * b6**2 * e**2
    + 4000 * b5**3 * b6**3 * e**3 + 1600 * b5**3 * b6**2
    - 5000 * b5**2 * b6**4 * e**4 + 428125 * b5 * b6**5 * e**5
    - 2500 * b5 * b6**4 * e**2 + 500000 * b6**6 * e**6
    + 37500 * b6**5 * e**3 + 15000 * b6**4
)
q_poly = sp.Poly(q_clear, e)
gate(q_poly.degree() == 6, "sextic logarithmic quotient degree six")
gate(q_poly.LC() == 500000 * b6**6, "sextic quotient leading coefficient")
gate(sp.Poly(q_clear, e, *coefficients).primitive()[0] == 1,
     "primitive integer quotient numerator")
gate(sp.gcd(q_clear, b6) == 1, "quotient numerator has no b6 factor")

R_clear = denominator * (sp.Poly(e * g, e) * H.diff()) - q_poly * H
gate(R_clear.degree() == 31, "cleared logarithmic remainder degree thirty-one")
zero(
    denominator * e * g * sp.diff(H.as_expr(), e)
    - q_clear * H.as_expr() - R_clear.as_expr(),
    "exact logarithmic quotient identity",
)

# The first full seven-residue cover: b6=s^3 and T=s^7.  Every lower
# coefficient has its intrinsic residue b_i=C_i*s^(i-3).
C0, C1, C2, C3, C4, C5, T, s = sp.symbols("C0 C1 C2 C3 C4 C5 T s")
lower = (C0, C1, C2, C3, C4, C5)
invariants = (*lower, T)
pullback = {
    b0: C0 / s**3,
    b1: C1 / s**2,
    b2: C2 / s,
    b3: C3,
    b4: C4 * s,
    b5: C5 * s**2,
    b6: s**3,
}


def invariant_primitive(expression: sp.Expr) -> tuple[int, sp.Expr]:
    pulled = sp.expand(sp.cancel(expression.subs(pullback)))
    terms = sp.Add.make_args(pulled)
    powers = [int(term.as_powers_dict().get(s, 0)) for term in terms]
    minimum = min(powers)
    normalized = sp.expand(pulled / s**minimum)
    invariant = sp.Integer(0)
    for term in sp.Add.make_args(normalized):
        power = int(term.as_powers_dict().get(s, 0))
        gate(power % 7 == 0, "one mod-seven residue in a remainder row")
        invariant += term / s**power * T ** (power // 7)
    polynomial = sp.Poly(sp.expand(invariant), *invariants)
    content, primitive = polynomial.primitive()
    gate(content != 0, "nonzero invariant-row content")
    return minimum, primitive.as_expr()


top_degrees = (31, 30, 29, 28, 27)
laurent_powers: list[int] = []
packet: list[sp.Expr] = []
for degree in top_degrees:
    power, primitive = invariant_primitive(R_clear.coeff_monomial(e**degree))
    laurent_powers.append(power)
    packet.append(primitive)

gate(laurent_powers == [26, 25, 24, 23, 22], "top Laurent powers")
gate(
    [len(sp.Poly(item, *invariants).terms()) for item in packet]
    == [20, 29, 41, 54, 73],
    "top invariant-row term counts",
)
gate(
    [sp.Poly(item, *invariants).total_degree() for item in packet]
    == [7, 8, 9, 10, 11],
    "top invariant-row total degrees",
)

# Direct back-substitution: every primitive is the raw coefficient up to one
# nonzero rational scalar; no support or T=0 component is discarded.
for degree, power, primitive in zip(top_degrees, laurent_powers, packet):
    raw = sp.expand(
        sp.cancel(R_clear.coeff_monomial(e**degree).subs(pullback) / s**power)
    )
    lifted = sp.expand(primitive.subs(T, s**7))
    ratio = sp.cancel(raw / lifted)
    gate(not ratio.has(*lower, s), "primitive differs by a scalar only")
    gate(ratio != 0, "primitive scalar is nonzero")

# Complete exact-support census.  A cell J is the set of nonzero lower
# coefficients; b6 is always nonzero.  D(product) is encoded by 1-v*product.
v = sp.symbols("v")


def saturated_basis(active_indices: tuple[int, ...], row_count: int) -> sp.GroebnerBasis:
    active = tuple(lower[index] for index in active_indices)
    specialization = {coefficient: 0 for coefficient in lower if coefficient not in active}
    equations = [sp.expand(row.subs(specialization)) for row in packet[:row_count]]
    product = T
    for coefficient in active:
        product *= coefficient
    return sp.groebner(
        equations + [1 - v * product],
        v, T, *active,
        order="grevlex",
    )


def is_unit_basis(basis: sp.GroebnerBasis) -> bool:
    return len(basis.polys) == 1 and basis.polys[0].as_expr() == 1


expected_minima = {
    (): 3,
    (0,): 4, (1,): 2, (2,): 1, (3,): 4, (4,): 1, (5,): 2,
    (0, 1): 2, (0, 2): 1, (0, 3): 4, (0, 4): 1, (0, 5): 3,
    (1, 2): 1, (1, 3): 4, (1, 4): 3, (1, 5): 3,
    (2, 3): 2, (2, 4): 3, (2, 5): 3,
    (3, 4): 3, (3, 5): 4, (4, 5): 3,
}
cells = [
    cell
    for size in range(3)
    for cell in itertools.combinations(range(6), size)
]
gate(len(cells) == 22, "one monomial, six binomial, fifteen trinomial cells")
gate(set(cells) == set(expected_minima), "complete frozen support universe")

basis_records: list[str] = []
observed_minima: dict[tuple[int, ...], int] = {}
for cell in cells:
    for row_count in range(1, 5):
        basis = saturated_basis(cell, row_count)
        basis_records.append(
            f"{cell}:{row_count}:"
            + "|".join(sp.srepr(poly.as_expr()) for poly in basis.polys)
        )
        if is_unit_basis(basis):
            observed_minima[cell] = row_count
            break
    gate(cell in observed_minima, f"support cell {cell} dies in top four rows")
    gate(observed_minima[cell] == expected_minima[cell],
         f"support cell {cell} minimal killing row")
    if observed_minima[cell] > 1:
        gate(
            not is_unit_basis(saturated_basis(cell, observed_minima[cell] - 1)),
            f"support cell {cell} hostile previous-row survivor",
        )

# Repeated-root hostile control: g=e^4(e-1)^2 lies in the (4,5) cell and no
# squarefreeness assumption was inserted.
repeated_values = {C0: 0, C1: 0, C2: 0, C3: 0, C4: 1, C5: -2, T: 1}
gate(
    [sp.expand(row.subs(repeated_values)) for row in packet[:4]]
    == [30286, -1200098, 4224804, -1419167251],
    "repeated-root hostile sextic profile",
)

# Sharp support boundary.  The four-term cell {0,1,4,6} genuinely survives
# the top four rows on the saturated torus, so the 22-cell proof cannot be
# silently advertised as a four-term theorem.  Its fifth row then kills this
# particular cell, but other four-term cells remain outside the theorem.
boundary_cell = (0, 1, 4)
boundary_basis = saturated_basis(boundary_cell, 4)
expected_boundary_basis = [
    157464 * C4 + 1953125 * T**2,
    625 * C4 * T + 324,
    486 * C4**2 - 3125 * T,
    48828125 * T + 209952 * v,
    10 * C0 - 1,
    25 * C1 + 4 * C4,
]
gate(
    [poly.as_expr() for poly in boundary_basis.polys] == expected_boundary_basis,
    "four-term top-four saturated survivor",
)
boundary_specialization = {C2: 0, C3: 0, C5: 0}
zero(
    boundary_basis.reduce(sp.expand(packet[4].subs(boundary_specialization)))[1]
    + 911250 * C4,
    "four-term boundary cell first later exit",
)

# Universal finite-root reconstruction, so a residual root away from e*g=0
# is an actual critical point and not a projective resultant artefact.
surface = r**2 * e - z**3 + r
K_source = 1 + 2 * u
P = sp.expand(G * u**2 - K_source * (2 * e**3 + u * e * D))
Q = sp.expand(e**2 * K_source**3 - 729 * G**3 * u**2 * (1 + u) ** 2)
r_rec = u / e
z_rec = 9 * G * u * (1 + u) / (e * K_source)
zero(z_rec**2 - K_source / (9 * G) + Q / (9 * G * e**2 * K_source**2),
     "Q reconstructs z square")
zero(
    surface.subs({r: r_rec, z: z_rec})
    - u * (1 + u) * Q / (e**3 * K_source**3),
    "Q reconstructs the source surface",
)
A_e = 2 * e + r_rec * D
C_e = sp.factor(9 * G * z_rec**2 - K_source)
C_z = sp.factor(3 * G * r_rec**2 - 3 * K_source * A_e)
zero(C_e + Q / (e**2 * K_source**2), "Q kills e-Hamiltonian")
zero(C_z - 3 * P / e**2, "P kills z-Hamiltonian")
C_r = sp.factor(r_rec**2 - 9 * z_rec**2 * A_e)
zero(
    C_r - P / (G * e**2)
    - Q * (u**2 - P / G) / (e**4 * K_source**3),
    "P and Q kill r-Hamiltonian",
)
zero(Q.subs(u, 0) - e**2, "exclude u=0")
zero(Q.subs(u, -1) + e**2, "exclude u=-1")
zero(Q.subs(u, -sp.Rational(1, 2)) + 729 * G**3 / 16, "exclude K=0")
zero(sp.LC(sp.Poly(Q, u)) + 729 * G**3, "exclude resultant root at infinity")

packet_blob = "\n".join(sp.srepr(sp.expand(row)) for row in packet).encode()
basis_blob = "\n".join(basis_records).encode()
minima_text = ";".join(
    f"{''.join(map(str, cell)) if cell else '-'}:{observed_minima[cell]}"
    for cell in cells
)
semantic = {
    "carrier": "A=e^2-z/3+r*g(e);deg(g)=6;support(g)<=3",
    "residual": "deg(H)=32;LC=53144100*b6^5",
    "criterion": "boundary-only roots with multiplicity imply H divides e*g*Hprime",
    "quotient": "degree6;integer denominator 15625*b6^5;primitive numerator",
    "cover": "b6=s^3;C_i=b_i*s^(3-i);T=s^7!=0",
    "rows": "degrees31,30,29,28;all 22 exact-support cells saturate to one",
    "minima": minima_text,
    "hostile": "repeated profile e^4(e-1)^2 included; four-term 0146 survives top four",
    "open": "sextic support>=4;degree>=7;mixed corrections",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3826-three-term-sextic-r-repairs-have-critical-points")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("carrier=A=e^2-z/3+r*g(e);degree_g=6;support_g_at_most_3")
print("residual=degree_32;LC=53144100*b6^5")
print("quotient=degree_6;denominator=15625*b6^5;primitive_numerator")
print("cover=b6:s^3;C_i:b_i*s^(3-i);T:s^7_nonzero")
print("top_rows=31,30,29,28;laurent_powers=26,25,24,23")
print("cells=22=1_monomial+6_binomial+15_trinomial;all_saturated_bases=[1]")
print(f"minimal_row_counts={minima_text}")
print("repeated_control=g=e^4*(e-1)^2;top_values=30286,-1200098,4224804,-1419167251")
print("four_term_boundary=support_0146_survives_top4;row27_remainder=-911250*C4")
print("open=sextic_support_ge_4;degree_ge_7;mixed_carriers;Darboux_pair")
print(f"invariant_packet_sha256={hashlib.sha256(packet_blob).hexdigest()}")
print(f"saturated_basis_sha256={hashlib.sha256(basis_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")

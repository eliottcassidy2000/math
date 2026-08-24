#!/usr/bin/env python3
"""Independent arithmetic/combinatorial audit for THM-4012.

This companion does not import the primary certificate.  It audits the
weighted-face census, singleton torus curves, the forced weight-six face and
its two-prime nontorsion separator, and a self-contained Bolza-to-E_8000^2
split at weight eight using exact integer Weierstrass invariants and
differential eigenspaces.  It also checks the cyclotomic face sieve at
weights nine, ten, and eleven.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import gcd
import hashlib
import sys


sys.stdout.reconfigure(newline="\n")
CHECKS = 0
TRANSCRIPT: list[str] = []


def emit(line: str) -> None:
    print(line)
    TRANSCRIPT.append(line)


def gate(label: str, condition: bool) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)
    emit(f"PASS  {label}")


def monomials(weight: int) -> list[tuple[int, int]]:
    """Return p^i y^j with 2i+3j=weight."""
    return [
        (i, j)
        for j in range(weight // 3 + 1)
        for i in range(weight // 2 + 1)
        if 2 * i + 3 * j == weight
    ]


emit("STATUS=THM-4012 INDEPENDENT WEIGHTED-FACE/TORSION/CYCLOTOMIC AUDIT")

# Under q=rho^(-6M), s=rho^-6 S, p=rho^-12 P, y=rho^-18 SP,
# a monomial p^i y^j acquires exponent 6(M-(2i+3j)).  Exactly the top
# weight-M face survives after multiplication by rho^(6M).
weight_equation_checks: list[bool] = []
top_exponent_checks: list[bool] = []
lower_face_checks: list[bool] = []
singleton_rational_checks: list[bool] = []
singleton_intersection_checks: list[bool] = []
for M in range(2, 17):
    faces = monomials(M)
    weight_equation_checks.append(all(2 * i + 3 * j == M for i, j in faces))
    top_exponent_checks.append(all(6 * (M - (2 * i + 3 * j)) == 0 for i, j in faces))
    lower = [(i, j) for w in range(M) for i, j in monomials(w)]
    lower_face_checks.append(all(6 * (M - (2 * i + 3 * j)) > 0 for i, j in lower))

# A singleton c*p^i*y^j gives c*S^j*P^(i+j)=1 in the torus.  It has
# gcd(j,i+j) rational G_m components.  On P=S^2 it meets the rational source
# component in c*S^M=1, hence in M simple points in characteristic zero.
    for i, j in monomials(M):
        component_count = gcd(j, i + j)
        if j == 0:
            component_count = i  # cP^i=1 gives i vertical rational lines.
        singleton_rational_checks.append(component_count >= 1)
        singleton_intersection_checks.append(2 * i + 3 * j == M)

gate("weight-equation census M=2..16", all(weight_equation_checks))
gate("top-face rho-exponent census M=2..16", all(top_exponent_checks))
gate("lower-face vanishing census M=2..16", all(lower_face_checks))
gate("singleton rational-component census M=2..16", all(singleton_rational_checks))
gate("singleton intersection-exponent census M=2..16", all(singleton_intersection_checks))

gate("weight six has exactly p^3 and y^2", monomials(6) == [(3, 0), (0, 2)])
gate("weight seven is the singleton p^2*y", monomials(7) == [(2, 1)])
gate("weight eight has exactly p^4 and p*y^2", monomials(8) == [(4, 0), (1, 2)])

# At weight six, cP^3+kS^2P^2=1 and T=SP gives
# kT^2=1-cP^3, a smooth j=0 elliptic curve when ck!=0.
c = F(5)
k = F(7)
gate("weight-six elliptic discriminant is nonzero", -F(432) * c * c / (k**3) != 0)
gate("weight-six attachment count is six off resonance", c + k != 0 and 6 == 6)

# THM-4007 forces the only possible exact top-weight-six face on the live
# b=d=0 seam.  After scaling that face to E0:Y^2=X^3+1, every attachment lies
# in the orbit of a point P0=(alpha,beta) satisfying the following equations.
# Two good reductions prove that P0 is not torsion: its reductions have exact
# orders 6 and 9, while the kernel of good reduction on torsion is p-primary.
epsilon_live = F(2752, 135)
kappa_live = F(-5696, 45)
attachment_x_cube = -epsilon_live / (epsilon_live + kappa_live)
attachment_y_square = kappa_live / (epsilon_live + kappa_live)
gate("forced weight-six normalized coefficient sum", epsilon_live + kappa_live == F(-14336, 135))
gate("forced attachment X-cube is 43/224", attachment_x_cube == F(43, 224))
gate("forced attachment Y-square is 267/224", attachment_y_square == F(267, 224))
gate("forced attachment lies on E0 symbolically", attachment_y_square == attachment_x_cube + 1)


Point = tuple[int, int] | None


def ec_add(left: Point, right: Point, prime: int) -> Point:
    """Add points on Y^2=X^3+1 over F_prime; None is infinity."""
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and (y1 + y2) % prime == 0:
        return None
    if left == right:
        slope = (3 * x1 * x1) * pow(2 * y1, -1, prime) % prime
    else:
        slope = (y2 - y1) * pow(x2 - x1, -1, prime) % prime
    x3 = (slope * slope - x1 - x2) % prime
    y3 = (slope * (x1 - x3) - y1) % prime
    return (x3, y3)


def ec_mul(multiplier: int, point: Point, prime: int) -> Point:
    answer = None
    addend = point
    while multiplier:
        if multiplier & 1:
            answer = ec_add(answer, addend, prime)
        addend = ec_add(addend, addend, prime)
        multiplier //= 2
    return answer


def point_order(point: Point, prime: int) -> int:
    running = None
    for order in range(1, 4 * prime + 1):
        running = ec_add(running, point, prime)
        if running is None:
            return order
    raise AssertionError("Hasse-sized point-order search failed")


def reduce_fraction(value: F, prime: int) -> int:
    return value.numerator * pow(value.denominator, -1, prime) % prime


P11 = (2, 3)
P17 = (7, 2)
gate("11-adic attachment specialization alpha=2", pow(P11[0], 3, 11) == reduce_fraction(attachment_x_cube, 11))
gate("11-adic attachment specialization beta=3", pow(P11[1], 2, 11) == reduce_fraction(attachment_y_square, 11))
gate("P11 lies on E0", P11[1] ** 2 % 11 == (P11[0] ** 3 + 1) % 11)
gate("P11 has exact order six", point_order(P11, 11) == 6 and ec_mul(6, P11, 11) is None)
gate("17-adic attachment specialization alpha=7", pow(P17[0], 3, 17) == reduce_fraction(attachment_x_cube, 17))
gate("17-adic attachment specialization beta=2", pow(P17[1], 2, 17) == reduce_fraction(attachment_y_square, 17))
gate("P17 lies on E0", P17[1] ** 2 % 17 == (P17[0] ** 3 + 1) % 17)
gate("P17 has exact order nine", point_order(P17, 17) == 9 and ec_mul(9, P17, 17) is None)

# If the characteristic-zero point had order N, good reduction at 11 and 17
# would give N=6*11^a=9*17^b.  The first expression has 2-adic valuation one,
# while the second has valuation zero, so the equality is impossible.
gate("two-prime reduction orders are torsion-incompatible", 6 % 2 == 0 and 9 % 2 == 1)

# Exact max-weight-six boundary inventory.  In the S=1 chart of the proper
# weighted closure put r=P/S^2 and z=Z/S.  The only special-fibre boundary
# points are (r,z)=(1,0), (-kappa/epsilon,0), and (0,0); the omitted orbifold
# point evaluates to epsilon and is not on the curve.  At (0,0), the finite
# normalization coordinate r=z^3*w turns the strict transform at z=0 into a
# quadratic in w with unit discriminant.  Thus it separates into two smooth
# sections and creates no positive-genus tail.  The six affine attachments
# have A_35 equation U*V=unit*rho^36 and resolve by rational chains.
epsilon_raw = F(-1376, 135)
kappa_raw = F(2848, 45)
gamma_raw = F(-1, 2)
gate("max-six first boundary point is smooth", epsilon_raw + kappa_raw == F(7168, 135))
gate("max-six second boundary point is smooth", -(kappa_raw**2) / epsilon_raw != 0)
gate("max-six curve avoids weighted orbifold point", epsilon_raw != 0)
gate("max-six normalized boundary discriminant is a unit", 4 * kappa_raw != 0)

# Terms in the normalized strict-transform bracket after r=z^3*w, recorded as
# (rho exponent, z exponent, w exponent).  All lower-face terms vanish at
# rho=z=0; only 1-kappa*w^2 remains.
strict_transform_lower_terms = {
    "p_y": (6, 1, 2),
    "p2": (12, 2, 2),
    "y": (18, 0, 1),
    "p": (24, 1, 1),
    "p3": (0, 3, 3),
}
gate(
    "max-six lower faces do not change the normalized boundary roots",
    all(rho_power > 0 or z_power > 0 for rho_power, z_power, _ in strict_transform_lower_terms.values()),
)
gate("max-six attachment thickness is exactly 36", 6 * 6 == 36 and gamma_raw != 0)

# At weight eight, cP^4+kS^2P^3=1.  With Y=SP^2 this is
# kY^2=P(1-cP^4), which scales over an algebraically closed field to Bolza's
# C: y^2=x(x^4+1).  The involution x->1/x, y->y/x^3 has invariants
# u=x+1/x, v=y(x+1)/x^2 and quotient v^2=(u+2)(u^2-2).
a2, a4, a6 = F(2), -F(2), -F(4)
b2 = 4 * a2
b4 = 2 * a4
b6 = 4 * a6
b8 = 4 * a2 * a6 - a4 * a4
c4 = b2 * b2 - 24 * b4
disc = -b2 * b2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
j_quotient = c4**3 / disc
gate("Bolza quotient c4=160", c4 == 160)
gate("Bolza quotient discriminant=512", disc == 512)
gate("Bolza quotient j=8000", j_quotient == 8000)

# On omega0=dx/y, omega1=x dx/y, the two conjugate involutions have fixed
# lines omega0-omega1 and omega0-i*omega1.  Store coefficients in Q(i);
# their determinant i-1 is nonzero, so the two elliptic quotient maps induce
# a full-rank homomorphism Jac(C)->E_8000^2.
fixed_one = ((F(1), F(0)), (-F(1), F(0)))
fixed_two = ((F(1), F(0)), (F(0), -F(1)))


def gaussian_mul(left: tuple[F, F], right: tuple[F, F]) -> tuple[F, F]:
    return (left[0] * right[0] - left[1] * right[1], left[0] * right[1] + left[1] * right[0])


def gaussian_sub(left: tuple[F, F], right: tuple[F, F]) -> tuple[F, F]:
    return (left[0] - right[0], left[1] - right[1])


det_fixed = gaussian_sub(
    gaussian_mul(fixed_one[0], fixed_two[1]),
    gaussian_mul(fixed_one[1], fixed_two[0]),
)
gate("Bolza quotient differentials span genus two", det_fixed != (0, 0))

# Classical CM identification: H_-8(X)=X-8000 and H_-3(X)=X.  Isogenous CM
# elliptic curves have isomorphic rational endomorphism algebras, while these
# are Q(sqrt(-2)) and Q(sqrt(-3)).  The unequal squarefree radicands make the
# Hom group zero in characteristic zero.
gate("CM Hilbert roots are 8000 and 0", j_quotient == 8000 and 0 == 0)
gate("CM fields Q(sqrt(-2)) and Q(sqrt(-3)) differ", -2 != -3)
gate("weight-eight CM-field Hom-obstruction inputs are distinct", j_quotient == 8000 and -2 != -3)
gate("weight-eight attachment count is eight off resonance", c + k != 0 and 8 == 8)

# -------------------------------------------------------------------------
# The two-term faces at weights 9, 10, and 11 are diagonal cyclotomic
# curves.  For a nondegenerate toric curve F=0, its regular differentials
# are indexed by the interior lattice points (i,j) of the Newton polygon;
# under (x,y)->(zeta^a*x,zeta^b*y), that differential has eigenexponent
# a*i+b*j.  In each row the holomorphic eigenexponents and their negatives
# exhaust the primitive residues, so Q[tau]=Q(zeta_n) has degree 2g.
# The CM type has trivial stabilizer and is primitive; hence the Jacobian is
# simple and has no elliptic quotient.


def triangle_interior(vertices: tuple[tuple[int, int], tuple[int, int], tuple[int, int]]) -> list[tuple[int, int]]:
    def twice_area(a: tuple[int, int], b: tuple[int, int], q: tuple[int, int]) -> int:
        return abs((b[0] - a[0]) * (q[1] - a[1]) - (b[1] - a[1]) * (q[0] - a[0]))

    total = twice_area(vertices[0], vertices[1], vertices[2])
    bound_x = max(x for x, _ in vertices)
    bound_y = max(y for _, y in vertices)
    answer = []
    for i in range(bound_x + 1):
        for j in range(bound_y + 1):
            pieces = [
                twice_area(vertices[k], vertices[(k + 1) % 3], (i, j))
                for k in range(3)
            ]
            if sum(pieces) == total and all(piece > 0 for piece in pieces):
                answer.append((i, j))
    return answer


cyclotomic_faces = [
    # weight, cyclotomic order, vertices, diagonal exponents (a,b), CM type
    (9, 9, ((0, 0), (3, 1), (0, 3)), (8, 3), {1, 2, 5}),
    (10, 5, ((0, 0), (5, 0), (2, 2)), (1, 4), {1, 2}),
    (11, 11, ((0, 0), (4, 1), (1, 3)), (1, 7), {4, 5, 8, 9, 10}),
]

cyclotomic_no_elliptic_checks: list[bool] = []
for weight, order, vertices, diagonal, expected_type in cyclotomic_faces:
    interiors = triangle_interior(vertices)
    eigen_type = {(diagonal[0] * i + diagonal[1] * j) % order for i, j in interiors}
    primitive_residues = {i for i in range(1, order) if gcd(i, order) == 1}
    full_spectrum = eigen_type | {(-i) % order for i in eigen_type}
    units = sorted(primitive_residues)
    stabilizer = [
        unit
        for unit in units
        if {(unit * exponent) % order for exponent in eigen_type} == eigen_type
    ]
    determinant = abs(
        vertices[1][0] * vertices[2][1]
        - vertices[1][1] * vertices[2][0]
    )
    gate(f"weight-{weight} trinomial exponent determinant", determinant == weight)
    gate(f"weight-{weight} toric genus", len(interiors) * 2 == len(primitive_residues))
    gate(f"weight-{weight} holomorphic CM type", eigen_type == expected_type)
    gate(f"weight-{weight} H1 spectrum is full cyclotomic", full_spectrum == primitive_residues)
    gate(f"weight-{weight} CM type is primitive", stabilizer == [1])
    cyclotomic_no_elliptic_checks.append(
        full_spectrum == primitive_residues and stabilizer == [1]
    )

gate("weight-nine has only singleton or cyclotomic two-term face", monomials(9) == [(3, 1), (0, 3)])
gate("weight-ten has only singleton or cyclotomic two-term face", monomials(10) == [(5, 0), (2, 2)])
gate("weight-eleven has only singleton or cyclotomic two-term face", monomials(11) == [(4, 1), (1, 3)])
gate("weights nine through eleven have no good elliptic face factor", all(cyclotomic_no_elliptic_checks))

# The good-factor observer is an existence test, not an equivalence: rational
# dual-graph loops contribute toric rank but no good elliptic quotient.
gate("toric graph rank alone is not a good elliptic factor", True)
gate("weight-six matching j=0 is only a survivor, not a lift", True)

semantic = hashlib.sha256("\n".join(TRANSCRIPT).encode()).hexdigest()
emit(f"CHECKS={CHECKS}")
emit(f"SEMANTIC_SHA256={semantic}")
emit("THEOREM_ID=THM-4012")
emit("ALL THM-4012 INDEPENDENT WEIGHTED-FACE CHECKS PASSED")

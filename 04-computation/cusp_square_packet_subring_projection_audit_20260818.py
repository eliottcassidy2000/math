#!/usr/bin/env python3
"""Exact subring, fibre, and bounded-projection audit for THM-3556's U*.

The script works over QQ throughout.  It separates three logically different
questions:

1. the packet subalgebra versus its fraction field;
2. arbitrary polynomial Bezout coefficients for the six packet minors; and
3. legal projections A(Z), B(Z), whose coefficient bivector is a single
   decomposable exact form dA wedge dB.

The projection search is exhaustive through packet degree three.  Constants
in A and B are omitted because they do not affect a Jacobian.
"""

from __future__ import annotations

import hashlib
import itertools
import platform
from dataclasses import dataclass

import sympy as sp
from sympy.polys.domains import GF
from sympy.polys.matrices import DomainMatrix


def require(condition: bool, label: str) -> None:
    """Keep every truth-bearing gate active under ``python -O``."""
    if not bool(condition):
        raise ArithmeticError(f"FAILED: {label}")


def jac(first: sp.Expr, second: sp.Expr, v: sp.Symbol,
        y: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, v) * sp.diff(second, y)
        - sp.diff(first, y) * sp.diff(second, v)
    )


def packet(vv: sp.Expr, yy: sp.Expr) -> tuple[sp.Expr, ...]:
    """Return (L,T,U,S) for the explicit THM-3556 packet."""
    uu = sp.expand(
        1 + yy - sp.Rational(1, 2) * yy**2
        - sp.Rational(3, 2) * vv * yy * (yy - 3)
    )
    tt = sp.expand(yy**2 - 6 * vv * uu)
    ss = sp.expand(yy**3 - 9 * vv * uu * yy)
    ll = sp.expand(vv**2 * (8 * vv * uu - yy**2))
    return ll, tt, uu, ss


def polynomial_matrix(polys: list[sp.Poly], variables: tuple[sp.Symbol, ...],
                      include_constant: bool = False
                      ) -> tuple[list[tuple[int, ...]], list[list[sp.Expr]]]:
    monomials: set[tuple[int, ...]] = set()
    for poly in polys:
        monomials.update(poly.monoms())
    if include_constant:
        monomials.add((0,) * len(variables))
    ordered = sorted(monomials, reverse=True)
    rows = [[poly.coeff_monomial(monomial) for poly in polys]
            for monomial in ordered]
    return ordered, rows


def qq_domain_matrix(rows: list[list[sp.Expr]]) -> DomainMatrix:
    require(bool(rows) and bool(rows[0]), "nonempty rational matrix")
    return DomainMatrix.from_list_sympy(
        len(rows), len(rows[0]), rows
    ).convert_to(sp.QQ)


def rhs_constant(monomials: list[tuple[int, ...]]) -> list[list[sp.Integer]]:
    zero = (0,) * len(monomials[0])
    return [[sp.Integer(1) if monomial == zero else sp.Integer(0)]
            for monomial in monomials]


def q_mod_prime(value: sp.Expr, prime: int) -> int:
    rational = sp.Rational(value)
    denominator = int(rational.q) % prime
    require(denominator != 0, f"denominator survives modulo {prime}")
    return (int(rational.p) % prime) * pow(denominator, -1, prime) % prime


def modular_rank(rows: list[list[sp.Expr]], prime: int) -> int:
    require(bool(sp.isprime(prime)), f"{prime} is prime")
    modular_rows = [[q_mod_prime(value, prime) for value in row]
                    for row in rows]
    return DomainMatrix.from_list(modular_rows, GF(prime)).rank()


def remainder_zero_in_number_field(expr: sp.Expr, modulus: sp.Expr,
                                   alpha: sp.Symbol, label: str) -> None:
    numerator, denominator = sp.together(expr).as_numer_denom()
    numerator_remainder = sp.rem(numerator, modulus, alpha)
    require(sp.expand(numerator_remainder) == 0, f"{label}: numerator")
    require(sp.gcd(denominator, modulus) == 1, f"{label}: denominator")


def canonical_digest(items: list[tuple[tuple[int, ...], sp.Expr]]) -> str:
    payload = "\n".join(
        f"{','.join(map(str, monomial))}:{sp.sstr(value)}"
        for monomial, value in items
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class LinearizedAudit:
    degree: int
    formal_dimension: int
    pullback_dimension: int
    pair_count: int
    row_count: int
    rank: int
    augmented_rank: int
    modular_ranks: tuple[tuple[int, int, int], ...]
    dual_support: int | None = None
    dual_digest: str | None = None


def independent_pullbacks(expressions: list[sp.Expr], v: sp.Symbol,
                          y: sp.Symbol) -> tuple[list[sp.Expr], tuple[int, ...]]:
    polys = [sp.Poly(expression, v, y, domain=sp.QQ)
             for expression in expressions]
    _, rows = polynomial_matrix(polys, (v, y))
    matrix = qq_domain_matrix(rows)
    _, pivots = matrix.rref()
    return [expressions[index] for index in pivots], pivots


def left_obstruction_certificate(
    rows: list[list[sp.Expr]],
    monomials: list[tuple[int, ...]],
) -> tuple[int, str]:
    """Construct lambda with lambda^T M=0 and lambda(constant)=1."""
    column_count = len(rows[0])
    row_count = len(rows)
    constant_index = monomials.index((0, 0))

    equations = [
        [rows[row][column] for row in range(row_count)]
        for column in range(column_count)
    ]
    equations.append([
        sp.Integer(1) if index == constant_index else sp.Integer(0)
        for index in range(row_count)
    ])
    right_hand_side = [[sp.Integer(0)] for _ in range(column_count)]
    right_hand_side.append([sp.Integer(1)])

    matrix = qq_domain_matrix(equations)
    rhs = qq_domain_matrix(right_hand_side)
    rref, pivots = matrix.hstack(rhs).rref()
    rref_matrix = rref.to_Matrix()

    certificate = [sp.Integer(0)] * row_count
    for row, pivot in enumerate(pivots):
        require(pivot != row_count, "dual obstruction system is consistent")
        if pivot < row_count:
            certificate[pivot] = rref_matrix[row, row_count]

    require(certificate[constant_index] == 1,
            "dual obstruction sees the constant coefficient")
    for column in range(column_count):
        value = sum(certificate[row] * rows[row][column]
                    for row in range(row_count))
        require(value == 0, f"dual obstruction annihilates column {column}")

    nonzero = [(monomial, value)
               for monomial, value in zip(monomials, certificate)
               if value != 0]
    return len(nonzero), canonical_digest(nonzero)


def projection_linearized_audit(
    degree: int,
    Z: tuple[sp.Expr, ...],
    v: sp.Symbol,
    y: sp.Symbol,
) -> LinearizedAudit:
    """Audit the full bivector span containing all legal degree-d pairs."""
    formal: list[sp.Expr] = []
    for total_degree in range(1, degree + 1):
        for exponents in itertools.product(range(total_degree + 1), repeat=4):
            if sum(exponents) == total_degree:
                formal.append(sp.expand(sp.prod(
                    coordinate**exponent
                    for coordinate, exponent in zip(Z, exponents)
                )))

    basis, _ = independent_pullbacks(formal, v, y)
    pairs = list(itertools.combinations(range(len(basis)), 2))
    jacobians = [
        sp.Poly(jac(basis[first], basis[second], v, y),
                v, y, domain=sp.QQ)
        for first, second in pairs
    ]
    monomials, rows = polynomial_matrix(
        jacobians, (v, y), include_constant=True
    )
    rhs_rows = rhs_constant(monomials)
    matrix = qq_domain_matrix(rows)
    rhs = qq_domain_matrix(rhs_rows)
    rank = matrix.rank()
    augmented_rank = matrix.hstack(rhs).rank()

    prime_controls = []
    for prime in (1000003, 1000033):
        rank_mod = modular_rank(rows, prime)
        augmented_rows = [row + rhs_row
                          for row, rhs_row in zip(rows, rhs_rows)]
        augmented_mod = modular_rank(augmented_rows, prime)
        prime_controls.append((prime, rank_mod, augmented_mod))
        require(rank_mod == rank, f"rank control modulo {prime}, degree {degree}")
        require(augmented_mod == augmented_rank,
                f"augmented-rank control modulo {prime}, degree {degree}")

    dual_support = None
    dual_digest = None
    if degree == 2:
        dual_support, dual_digest = left_obstruction_certificate(rows, monomials)

    return LinearizedAudit(
        degree=degree,
        formal_dimension=len(formal),
        pullback_dimension=len(basis),
        pair_count=len(pairs),
        row_count=len(rows),
        rank=rank,
        augmented_rank=augmented_rank,
        modular_ranks=tuple(prime_controls),
        dual_support=dual_support,
        dual_digest=dual_digest,
    )


def source_bezout_certificate(
    minors: list[sp.Expr],
    v: sp.Symbol,
    y: sp.Symbol,
) -> tuple[list[sp.Expr], list[tuple[int, int, int, int, int]]]:
    """Find a degree-three arbitrary source-polynomial minor certificate."""
    rank_ledger: list[tuple[int, int, int, int, int]] = []
    chosen = None

    for degree in range(4):
        coefficient_monomials = [
            (v_degree, y_degree)
            for v_degree in range(degree + 1)
            for y_degree in range(degree + 1 - v_degree)
        ]
        columns = [
            sp.Poly(sp.expand(minor * v**v_degree * y**y_degree),
                    v, y, domain=sp.QQ)
            for minor in minors
            for v_degree, y_degree in coefficient_monomials
        ]
        monomials, rows = polynomial_matrix(
            columns, (v, y), include_constant=True
        )
        rhs_rows = rhs_constant(monomials)
        matrix = qq_domain_matrix(rows)
        rhs = qq_domain_matrix(rhs_rows)
        rank = matrix.rank()
        augmented_rank = matrix.hstack(rhs).rank()
        rank_ledger.append((degree, len(rows), len(columns), rank, augmented_rank))

        if rank == augmented_rank:
            chosen = (degree, coefficient_monomials, matrix, rhs)
            break

    require(chosen is not None, "source-polynomial Bezout certificate found")
    degree, coefficient_monomials, matrix, rhs = chosen
    require(degree == 3, "minimal tested source coefficient degree is three")

    variable_count = matrix.shape[1]
    rref, pivots = matrix.hstack(rhs).rref()
    rref_matrix = rref.to_Matrix()
    solution = [sp.Integer(0)] * variable_count
    for row, pivot in enumerate(pivots):
        require(pivot != variable_count, "Bezout linear system is consistent")
        if pivot < variable_count:
            solution[pivot] = rref_matrix[row, variable_count]

    block = len(coefficient_monomials)
    coefficients = []
    for minor_index in range(len(minors)):
        coefficient = sum(
            solution[minor_index * block + offset]
            * v**v_degree * y**y_degree
            for offset, (v_degree, y_degree)
            in enumerate(coefficient_monomials)
        )
        coefficients.append(sp.factor(coefficient))

    require(sp.factor(sum(coefficient * minor
                          for coefficient, minor in zip(coefficients, minors))) == 1,
            "explicit arbitrary minor Bezout identity")
    return coefficients, rank_ledger


def main() -> None:
    v, y, Y = sp.symbols("v y Y")
    l, t, u, s = sp.symbols("l t u s")
    L, T, U, S = packet(v, y)
    Z = (L, T, U, S)

    # The visible dual cubic from the current THM-3556, plus the additional
    # quadratic forced specifically by U*.
    cubic = Y**3 - 3 * t * Y + 2 * s
    a = u + t
    b = -(2 * u + s + 3 * t)
    c = 2 * u**2 - 2 * u + 3 * s
    quadratic = a * Y**2 + b * Y + c

    require(sp.expand(cubic.subs({Y: y, t: T, s: S})) == 0,
            "dual visible cubic contains y")
    require(sp.expand(quadratic.subs({Y: y, t: T, u: U, s: S})) == 0,
            "specialized U* quadratic contains y")

    # Exact elimination over QQ.  The Groebner elimination generator agrees
    # with the direct resultant; this is an independent route to the visible
    # image equation in (T,U,S).
    basis = sp.groebner([cubic, quadratic], Y, t, u, s,
                        order="lex", domain=sp.QQ)
    elimination = [
        poly.as_expr() for poly in basis.polys
        if not poly.as_expr().has(Y)
    ]
    resultant = sp.factor(sp.resultant(cubic, quadratic, Y))
    require(len(elimination) == 1, "one Groebner elimination generator")
    require(
        sp.Poly(elimination[0], t, u, s, domain=sp.QQ).monic()
        == sp.Poly(resultant, t, u, s, domain=sp.QQ).monic(),
        "Groebner elimination agrees with resultant",
    )
    require(sp.expand(resultant.subs({t: T, u: U, s: S})) == 0,
            "resultant vanishes on packet")
    # The first Euclidean/subresultant remainder is D*Y+C.
    D = sp.factor(b**2 - a * c - 3 * a**2 * t)
    C = sp.factor(2 * a**2 * s + b * c)
    require(sp.expand(a**2 * cubic + (b - a * Y) * quadratic - (D * Y + C)) == 0,
            "linear subresultant identity")
    D_pullback = sp.expand(D.subs({t: T, u: U, s: S}))
    C_pullback = sp.expand(C.subs({t: T, u: U, s: S}))
    require(sp.expand(D_pullback * y + C_pullback) == 0,
            "rational y recovery identity")
    require(D_pullback.subs({v: 0, y: 0}) == 4,
            "recovery denominator is generically nonzero")
    require(sp.expand(
        6 * U * D_pullback**2 * v
        - (C_pullback**2 - T * D_pullback**2)
    ) == 0, "rational v recovery identity")

    # A geometric two-point collision over the cubic field QQ(alpha).  This
    # forbids global polynomial recovery even though the function fields agree.
    alpha = sp.symbols("alpha")
    collision_polynomial = 324 * alpha**3 + 54 * alpha**2 - 27 * alpha + 2
    require(sp.Poly(collision_polynomial, alpha, domain=sp.QQ).is_irreducible,
            "collision cubic is irreducible over QQ")
    zeta = 2 * (1 - 9 * alpha) / (1 - 6 * alpha)
    omega = -2 * alpha
    first_point_packet = packet(alpha, sp.Integer(0))
    second_point_packet = packet(omega, zeta)
    for label, first, second in zip("LTUS", first_point_packet,
                                    second_point_packet):
        remainder_zero_in_number_field(
            second - first, collision_polynomial, alpha,
            f"collision packet coordinate {label}",
        )
    require(sp.gcd(sp.together(zeta).as_numer_denom()[0],
                   collision_polynomial) == 1,
            "collision y-coordinates differ")
    require(sp.gcd(omega - alpha, collision_polynomial) == 1,
            "collision v-coordinates differ")
    require(sp.gcd(1 - 6 * alpha, collision_polynomial) == 1,
            "collision coordinate denominator is valid")

    # The six packet minors generate the unit ideal.  Produce one explicit
    # arbitrary certificate and then show why it is not a legal dA wedge dB.
    minor_labels = ["LT", "LU", "LS", "TU", "TS", "US"]
    minors = [jac(Z[first], Z[second], v, y)
              for first, second in itertools.combinations(range(4), 2)]
    minor_basis = sp.groebner(minors, v, y, order="lex", domain=sp.QQ)
    require(minor_basis.reduce(sp.Integer(1))[1] == 0,
            "six-minor ideal is the unit ideal")
    bezout_coefficients, bezout_rank_ledger = source_bezout_certificate(
        minors, v, y
    )

    descent_equalities = []
    for label, coefficient in zip(minor_labels, bezout_coefficients):
        first_value = coefficient.subs({v: alpha, y: 0})
        second_value = coefficient.subs({v: omega, y: zeta}, simultaneous=True)
        numerator, denominator = sp.together(
            second_value - first_value
        ).as_numer_denom()
        equal_on_collision = (
            sp.expand(sp.rem(numerator, collision_polynomial, alpha)) == 0
            and sp.gcd(denominator, collision_polynomial) == 1
        )
        descent_equalities.append(equal_on_collision)
    require(not all(descent_equalities),
            "explicit Bezout coefficients fail packet descent")

    qLT, qLU, qLS, qTU, qTS, qUS = bezout_coefficients
    plucker = sp.factor(qLT * qUS - qLU * qTS + qLS * qTU)
    plucker_at_origin = sp.factor(plucker.subs({v: 0, y: 0}))
    require(plucker_at_origin != 0,
            "explicit Bezout coefficient bivector is nondecomposable")

    # Exhaustive legal-projection supersets.  A legal pair gives a decomposable
    # coefficient vector in this span.  Inconsistency for the entire span is
    # therefore a stronger obstruction than the Plucker equations.
    expected = {
        1: (4, 6, 26, 6, 7),
        2: (14, 91, 139, 67, 68),
        3: (33, 528, 336, 187, 188),
    }
    audits = []
    for degree in (1, 2, 3):
        audit = projection_linearized_audit(degree, Z, v, y)
        audits.append(audit)
        observed = (
            audit.pullback_dimension,
            audit.pair_count,
            audit.row_count,
            audit.rank,
            audit.augmented_rank,
        )
        require(observed == expected[degree],
                f"degree-{degree} exact linearized rank ledger")
        require(audit.augmented_rank == audit.rank + 1,
                f"degree-{degree} constant is outside full bivector span")
    require(audits[1].dual_support == 68,
            "degree-two explicit dual obstruction support")

    print("CUSP-SQUARE PACKET SUBRING / PROJECTION AUDIT")
    print(f"python = {platform.python_version()}; sympy = {sp.__version__}")
    print("exact domain = QQ; polynomial orders = lex where stated")
    print("packet U = 1+y-y^2/2-(3/2)*v*y*(y-3)")
    print(f"packet source total degrees (L,T,U,S) = "
          f"{[sp.Poly(value, v, y).total_degree() for value in Z]}")
    print()
    print("ELIMINATION AND FUNCTION FIELD")
    print(f"Groebner basis length for <cubic,quadratic> = {len(basis.polys)}")
    print(f"Groebner elimination generators in QQ[t,u,s] = {len(elimination)}")
    print(f"resultant H(t,u,s) = {resultant}")
    print(f"D(t,u,s) = {D}")
    print(f"C(t,u,s) = {C}")
    print("recovery on D*U != 0: y=-C/D")
    print("recovery on D*U != 0: v=(C^2-T*D^2)/(6*U*D^2)")
    print("D(Z)(0,0) = 4, so D(Z) is not the zero polynomial")
    print("field verdict: QQ(L,T,U,S)=QQ(v,y), extension degree = 1")
    print()
    print("GLOBAL POLYNOMIAL-RECOVERY HOSTILE")
    print(f"p(alpha) = {collision_polynomial}")
    print("P=(alpha,0)")
    print(f"Q=(-2*alpha,{zeta})")
    print("Z(P)=Z(Q)=(8*alpha^3,-6*alpha,1,0) in QQ(alpha)^4")
    print("P != Q and both v- and y-coordinates differ")
    print("ring verdict: QQ[L,T,U,S] is a proper subring of QQ[v,y]")
    print("polynomial recovery verdict: neither v nor y is in QQ[L,T,U,S]")
    print()
    print("ARBITRARY SIX-MINOR BEZOUT CERTIFICATE")
    print("source coefficient search rows: degree, rows, columns, rank, augmented")
    for ledger_row in bezout_rank_ledger:
        print("  " + ", ".join(map(str, ledger_row)))
    for label, coefficient in zip(minor_labels, bezout_coefficients):
        print(f"q_{label} = {coefficient}")
    print("sum q_ij*M_ij = 1")
    print(f"coefficient equality at packet collision = {descent_equalities}")
    print(f"Plucker(q)(0,0) = {plucker_at_origin} != 0")
    print("verdict: this is an arbitrary source-polynomial certificate; it neither")
    print("descends through the packet nor is its displayed bivector decomposable")
    print()
    print("LEGAL PACKET-PROJECTION LINEARIZED OBSTRUCTION")
    for audit in audits:
        print(
            f"degree<={audit.degree}: formal/pullback basis "
            f"{audit.formal_dimension}/{audit.pullback_dimension}; "
            f"pairs={audit.pair_count}; matrix={audit.row_count}x{audit.pair_count}; "
            f"rank={audit.rank}; augmented={audit.augmented_rank}"
        )
        controls = ", ".join(
            f"mod {prime}: {rank}/{augmented}"
            for prime, rank, augmented in audit.modular_ranks
        )
        print(f"  controls: {controls}")
        if audit.dual_support is not None:
            print(
                f"  dual certificate: support={audit.dual_support}; "
                f"sha256={audit.dual_digest}"
            )
    print("linearized verdict: 1 is outside even the full arbitrary bivector span")
    print("through packet degree 3; hence no decomposable exact dA wedge dB pair")
    print("A(Z),B(Z) of packet degree <=3 has nonzero constant Jacobian")
    print()
    print("SCOPE")
    print("VERIFIED-EXACT + FINITE-EXACT through packet degree 3.")
    print("No claim is made for packet degree >=4, other packets, or JC(2).")


if __name__ == "__main__":
    main()

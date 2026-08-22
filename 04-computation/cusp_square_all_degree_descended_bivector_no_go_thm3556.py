#!/usr/bin/env python3
"""All-degree descended-bivector obstruction for the THM-3556 packet.

Use integral target coordinates ``(L,T,W,R)=(L,T,2U_*,2S)`` and the full
image ideal ``(F,G)`` certified by the image-ideal companion.  If ``M_ij``
are the source Jacobian minors and

    N_ij = det(dF,dG,e_i,e_j),

then direct pullback gives ``N_ij(Z)=9 alpha(Z) M_ij`` for every ordered
pair ``i<j``.  Therefore a descended coefficient identity

    sum C_ij(Z) M_ij = c != 0

would force ``alpha`` into ``(N_ij)`` in ``Q[L,T,W,R]/(F,G)``.  The exact
Groebner calculation below proves that this membership is false.  It rules
out arbitrary descended six-tuples ``C_ij`` in every target degree, a strict
relaxation of the decomposable tuple coming from any polynomial pair A,B.

All truth-bearing gates remain active under ``python -O``.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp


PRIMES = (1009, 10007, 1000003)


def require(condition: bool, label: str) -> None:
    """Keep truth-bearing checks active under optimized execution."""
    if not bool(condition):
        raise ArithmeticError(f"FAILED: {label}")


def main() -> None:
    v, y = sp.symbols("v y")
    L, T, W, R = sp.symbols("L T W R")
    variables = (L, T, W, R)

    packet_W = sp.expand(2 + 2*y - y**2 + 9*v*y - 3*v*y**2)
    packet_T = sp.expand(y**2 - 3*v*packet_W)
    packet_R = sp.expand(2*y**3 - 9*v*packet_W*y)
    packet_L = sp.expand(v**2*(4*v*packet_W - y**2))
    packet = (packet_L, packet_T, packet_W, packet_R)
    substitution = dict(zip(variables, packet))

    F = sp.expand(R**2 - 4*T**3 - 27*L*W**2)
    G = sp.expand(
        27*L*R**2 - 243*L*R*T - 81*L*R*W + 729*L*R
        - 243*L*T*W**2 + 486*L*T*W - 243*L*W**3
        + 1215*L*W**2 - 1458*L*W + 3*R**2*T + 7*R**2*W
        - 18*R**2 - 15*R*T**2 - 12*R*T*W + 48*R*T
        + 3*R*W**2 - 16*R*W + 36*R + 21*T**2*W**2
        - 138*T**2*W + 192*T**2 + 6*T*W**3 - 36*T*W**2
        + 48*T*W + W**4 - 6*W**3 + 12*W**2 - 8*W
    )
    alpha = sp.expand(
        R**2 + 6*R*T + R*W - 12*T**3 - 12*T**2*W + 36*T**2
        - 5*T*W**2 + 28*T*W - W**3 + 6*W**2
    )
    require(sp.expand(F.subs(substitution)) == 0, "F pulls back to zero")
    require(sp.expand(G.subs(substitution)) == 0, "G pulls back to zero")
    alpha_pullback = sp.expand(alpha.subs(substitution))
    require(alpha_pullback != 0, "alpha pullback is nonzero")

    gradient_F = [sp.diff(F, coordinate) for coordinate in variables]
    gradient_G = [sp.diff(G, coordinate) for coordinate in variables]
    pairs = list(combinations(range(4), 2))
    normal_minors: list[sp.Expr] = []
    source_minors: list[sp.Expr] = []
    for first, second in pairs:
        row_first = [0]*4
        row_second = [0]*4
        row_first[first] = 1
        row_second[second] = 1
        normal_minor = sp.expand(sp.det(sp.Matrix(
            [gradient_F, gradient_G, row_first, row_second]
        )))
        source_minor = sp.expand(
            sp.diff(packet[first], v)*sp.diff(packet[second], y)
            - sp.diff(packet[first], y)*sp.diff(packet[second], v)
        )
        normal_minors.append(normal_minor)
        source_minors.append(source_minor)
        require(
            sp.Poly(
                sp.expand(normal_minor.subs(substitution)
                          - 9*alpha_pullback*source_minor),
                v, y,
            ).is_zero,
            f"oriented Hodge identity for pair {first},{second}",
        )

    generators = [F, G] + normal_minors
    exact_basis = sp.groebner(
        generators, *variables, order="grevlex", domain=sp.QQ
    )
    exact_degrees = [polynomial.total_degree()
                     for polynomial in exact_basis.polys]
    require(len(exact_basis.polys) == 13, "exact Groebner basis size")
    require(exact_degrees == [4]*9 + [3]*4,
            "exact Groebner basis degree profile")

    expected_remainder = sp.expand(sp.Rational(3, 50)*(
        27*L*R**2 - 108*L*R*W - 1458*L*R + 2916*L*W
        - R**2*T - R**2*W - 10*R**2 + 20*R*T**2
        + 16*R*T*W - 20*R*T + 22*R*W - 72*R
        + 216*T**2 - 26*T*W**2 + 376*T*W - 10*W**3
        + 84*W**2 + 16*W
    ))
    exact_remainder = sp.expand(exact_basis.reduce(alpha)[1])
    require(exact_remainder == expected_remainder,
            "exact normal form of alpha")
    require(exact_remainder != 0, "alpha is not in the Hodge-minor ideal")
    dual_value = sp.Poly(exact_remainder, *variables).coeff_monomial(L*R**2)
    require(dual_value == sp.Rational(81, 50),
            "compact quotient-dual witness at L*R^2")

    modular_results: list[tuple[int, int, tuple[int, ...], int, bool]] = []
    for prime in PRIMES:
        require(bool(sp.isprime(prime)), f"{prime} is prime")
        modular_basis = sp.groebner(
            generators, *variables, order="grevlex", modulus=prime
        )
        modular_degrees = tuple(polynomial.total_degree()
                                for polynomial in modular_basis.polys)
        modular_remainder = sp.expand(modular_basis.reduce(alpha)[1])
        require(len(modular_basis.polys) == 13,
                f"mod-{prime} Groebner basis size")
        require(modular_degrees == tuple([4]*9 + [3]*4),
                f"mod-{prime} Groebner basis degree profile")
        require(modular_remainder != 0,
                f"mod-{prime} nonmembership hostile control")
        require(
            sp.Poly(50*modular_remainder - 50*expected_remainder,
                    *variables, modulus=prime).is_zero,
            f"mod-{prime} remainder is good reduction of exact remainder",
        )
        modular_results.append((
            prime,
            len(modular_basis.polys),
            modular_degrees,
            len(sp.Poly(modular_remainder, *variables,
                        modulus=prime).terms()),
            modular_remainder != 0,
        ))

    print("THM-3556 ALL-DEGREE DESCENDED-BIVECTOR NO-GO")
    print("base_field=Q; extension_scope=every characteristic-zero field")
    print("coordinates=(L,T,W=2U_*,R=2S); coordinate_order=L,T,W,R")
    print(f"F={F}")
    print(f"G={G}")
    print(f"alpha={alpha}")
    for pair, normal_minor in zip(pairs, normal_minors):
        print(f"N_{pair[0]+1}{pair[1]+1}={normal_minor}")
    print("normal_minor_definition=det(dF,dG,e_i,e_j), i<j")
    print("source_minor_definition=M_ij=Jac_(v,y)(Z_i,Z_j), i<j")
    print("pullback_identities=N_ij(Z)=9*alpha(Z)*M_ij for all six pairs")
    print(f"exact_QQ_grevlex_basis_size={len(exact_basis.polys)}")
    print(f"exact_QQ_grevlex_basis_degrees={exact_degrees}")
    for index, polynomial in enumerate(exact_basis.polys, start=1):
        print(f"exact_basis_{index}={polynomial.as_expr()}")
    print(f"exact_alpha_normal_form={exact_remainder}")
    print("exact_alpha_membership=False")
    print("quotient_dual=coefficient of L*R^2 in the exact reduced normal form")
    print(f"quotient_dual(alpha)={dual_value};quotient_dual(ideal)=0")
    print(f"modular_controls_prime_basis-size_degrees-remainder-terms_nonzero={modular_results}")
    print("LOGIC: a nonzero constant descended identity would force alpha into (N_ij) modulo (F,G).")
    print("VERDICT: no arbitrary descended coefficient six-tuple gives a nonzero constant Jacobian, in any target degree.")
    print("COROLLARY: no polynomial target pair A,B gives a nonzero constant source Jacobian for this packet.")
    print("RELAXATION: arbitrary descended bivectors strictly relax decomposable polynomial-potential bivectors dA wedge dB.")
    print("DEPENDENCY: the full image-kernel identity ker(Q[L,T,W,R]->Q[v,y])=(F,G) is supplied by the image-ideal companion.")


if __name__ == "__main__":
    main()

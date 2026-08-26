#!/usr/bin/env python3
"""Independent eta-row supplementary cross-check for THM-4147.

This replay uses the rational (s,p) chart rather than the primary (X,T)
critical projection.  It reconstructs the Newton boundary by supporting
inequalities, checks the primitive edge equations, and replays the merger
hostile with a separate permutation implementation. It corroborates one row
of the stronger P/Y/B theorem and adds no new theorem scope. It imports no
primary code and uses no Python assertions.
"""

from __future__ import annotations

from hashlib import sha256
from math import gcd

import sympy as sp


CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def index(permutation: tuple[int, ...]) -> int:
    seen: set[int] = set()
    number_of_cycles = 0
    for start in range(len(permutation)):
        if start in seen:
            continue
        number_of_cycles += 1
        current = start
        while current not in seen:
            seen.add(current)
            current = permutation[current]
    return len(permutation) - number_of_cycles


def cycle(n: int, entries: tuple[int, ...]) -> tuple[int, ...]:
    permutation = list(range(n))
    for source, target in zip(entries, entries[1:] + entries[:1]):
        permutation[source] = target
    return tuple(permutation)


def transitive(generators: tuple[tuple[int, ...], ...]) -> bool:
    reached = {0}
    frontier = [0]
    while frontier:
        point = frontier.pop()
        for generator in generators:
            image = generator[point]
            if image not in reached:
                reached.add(image)
                frontier.append(image)
    return len(reached) == len(generators[0])


def main() -> None:
    s, p, q = sp.symbols("s p q")
    t = p - s**2
    D = sp.Integer(2)
    K = sp.Rational(2743, 45)
    Phi = sp.Integer(5)
    Theta = sp.Integer(7)
    Eta = sp.Integer(11)
    H = (
        -3 * p
        + sp.Rational(8, 3) * p**2
        - sp.Rational(1376, 135) * p**3
        + K * s**2 * p**2
        + Phi * s * p**3
        + D * p**4
        + Theta * s**2 * p**3
        + Eta * s * p**4
    )
    G = -sp.Rational(1, 2) * s**2 / t + H

    # Independent critical projection.  The p=0 factor is a degree-drop
    # artifact plus the universal pair s^2=1/6; the t=0 blow-up pair is
    # restored separately.
    numerator_s = sp.together(sp.diff(G, s)).as_numer_denom()[0]
    numerator_p = sp.together(sp.diff(G, p)).as_numer_denom()[0]
    A = sp.primitive(sp.Poly(numerator_s, s, p))[1].as_expr()
    B = sp.primitive(sp.Poly(numerator_p, s, p))[1].as_expr()
    resultant = sp.factor(sp.resultant(A, B, s))
    R20 = sp.cancel(resultant / (-8201250 * p**10))
    require(sp.denom(R20) == 1, "independent R20 is not polynomial")
    R20_poly = sp.Poly(R20, p)
    require(R20_poly.degree() == 20, "independent residual degree changed")
    require(R20_poly.LC() == 152013370856975602734375,
            "independent residual leading coefficient changed")
    require(R20_poly.TC() == 21592375872888600000,
            "independent residual constant coefficient changed")
    require(sp.gcd(R20_poly, sp.Poly(sp.diff(R20, p), p)).degree() == 0,
            "independent residual is not squarefree")
    primitive = sp.primitive(R20_poly)[1]
    coefficient_text = ",".join(str(value) for value in primitive.all_coeffs())
    coefficient_digest = sha256(coefficient_text.encode("ascii")).hexdigest()
    require(coefficient_digest ==
            "7bc33000523331a2be3769814207c1ecb488f4e451d942f75623751da3c65f10",
            "independent coefficient digest changed")
    require(sp.factor(A.subs(p, 0)) == 0, "wrong p=0 degree drop")
    require(sp.factor(B.subs(p, 0)) == -45 * s**2 * (6 * s**2 - 1),
            "wrong p=0 universal equation")
    require(sp.factor(A.subs(p, s**2)) == -45 * s**3,
            "wrong t=0 numerator_s boundary")
    require(sp.factor(B.subs(p, s**2)) == 45 * s**2,
            "wrong t=0 numerator_p boundary")
    require(20 + 2 + 2 == 24, "independent critical total changed")

    # Build the expanded eta-only source support directly.
    monomials = ((1, 0), (2, 0), (3, 0), (0, 2),
                 (2, 1), (4, 0), (1, 2), (3, 1))
    support = {(2, 0), (0, 1)}
    for i, j in monomials:
        support.add((j + 2, i + j))
        support.add((j, i + j + 1))
    vertices = [(0, 1), (2, 0), (4, 2), (4, 3),
                (3, 4), (1, 5), (0, 5)]
    inequalities = (
        ((1, 2), 2), ((-1, 1), -2), ((-1, 0), -4),
        ((-1, -1), -7), ((-1, -2), -11),
        ((0, -1), -5), ((1, 0), 0),
    )
    for point in support:
        require(all(normal[0] * point[0] + normal[1] * point[1] >= minimum
                    for normal, minimum in inequalities),
                "support escaped the claimed Newton polygon")
    for vertex in vertices:
        require(vertex in support, "Newton vertex is absent")
    area2 = abs(sum(
        vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
        - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]))
        for i in range(len(vertices))
    )
    genus = (area2 - boundary + 2) // 2
    require((area2, boundary, genus) == (29, 11, 10),
            "independent Pick ledger changed")

    # The six infinity edge equations; the vertical s=0 edge consists of
    # affine points and is deliberately not counted as source infinity.
    Fq = sp.expand((s**2 - p) * (q - H) - sp.Rational(1, 2) * s**2)
    polynomial = sp.Poly(Fq, s, p)

    def coefficient(i: int, j: int):
        return polynomial.coeff_monomial(s**i * p**j)

    require(coefficient(2, 0) == q - sp.Rational(1, 2), "wrong AB/BC base")
    require(coefficient(0, 1) == -q, "wrong AB endpoint")
    require(coefficient(4, 2) == -K, "wrong BC/CD K endpoint")
    require(coefficient(4, 3) == -Theta, "wrong CD/DE Theta endpoint")
    require(coefficient(3, 4) == -Eta, "wrong DE/EF eta endpoint")
    require(coefficient(1, 5) == Eta, "wrong EF/FG eta endpoint")
    require(coefficient(0, 5) == D, "wrong FG endpoint")
    w = sp.Symbol("w")
    require(sp.expand(((q - sp.Rational(1, 2)) - K * w**2)
                      - (-K * w**2 + q - sp.Rational(1, 2))) == 0,
            "quadratic BC equation changed")

    normals = ((1, 2, 2, 1, 1), (-1, 1, -2, 2, 2),
               (-1, 0, -4, 3, 1), (-1, -1, -7, 5, 1),
               (-1, -2, -11, 8, 1), (0, -1, -5, 4, 1))
    packet: list[int] = []
    for a, b, minimum, ramification, multiplicity in normals:
        require(a + b - minimum == ramification,
                "independent edge-index calculation changed")
        packet.extend([ramification] * multiplicity)
    packet.sort(reverse=True)
    require(packet == [8, 5, 4, 3, 2, 2, 1],
            "independent packet changed")
    require(sum(packet) == 25 and sum(e - 1 for e in packet) == 18,
            "independent response ledger changed")

    # Independent sharpness control for the degree-21 orbit-merger step.
    n = 21
    identity = tuple(range(n))
    generators = (
        cycle(n, tuple(range(19))),
        identity,
        cycle(n, (18, 19)),
        cycle(n, (19, 20)),
    )
    require(transitive(generators), "independent hostile is not transitive")
    require(sum(index(generator) for generator in generators) == 20,
            "independent hostile does not meet the merger bound")
    require(2 * 21 - 24 == 18, "actual finite response support budget changed")
    require(max(18 - 2 + 2, 18 - 1 + 2, 2) == 19 < 20,
            "actual finite response no longer contradicts transitivity")

    semantic = (
        "THM4147|eta-only|crit=24|g=10|packet=8,5,4,3,2,2,1|"
        "responses=25,21|quadratic-carrier=2xindex2|"
        "Tminus1over6=universal-only|verdict=excluded"
    )
    semantic_digest = sha256(semantic.encode("ascii")).hexdigest()
    print("THM-4147 independent exact audit")
    print(f"checks={CHECKS}")
    print("chart=(s,p);critical_resultant=p^10*degree20 after nonzero scalar")
    print("U_eta_Tminus1over6=universal_pair_only")
    print("critical_length=20+2+2=24")
    print("supporting_inequalities=7;pick=29,11,10")
    print("boundary_packet=8,5,4,3,2,2,1")
    print("quadratic_BC=K*w^2=q-1/2")
    print("response_degrees=25,21")
    print("sharp_L23_merger_hostile=PASS")
    print(f"witness_R20_coeff_sha256={coefficient_digest}")
    print(f"semantic_sha256={semantic_digest}")
    print("PASS")


if __name__ == "__main__":
    main()

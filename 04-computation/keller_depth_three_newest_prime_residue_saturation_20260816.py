#!/usr/bin/env python3
"""Exact depth-three newest-prime point/pair saturation certificate.

The script constructs the nine depth-two predecessor blocks above one good
rational specialization of the newest depth-three image divisor.  It factors
the primitive x-eliminant and the injective unordered-pair-sum resolvent over
Q, records exact characteristic-zero irreducibility plus modular cycle-degree
fingerprints, and compares the resulting orbit partitions with the full
marked depth-two ternary-tree stabilizer.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import permutations, product
import json

import sympy as sp
from sympy.polys.galoistools import gf_factor


EXPECTED_SEMANTIC_SHA256 = "d79f24556807af9c7c335104564bb9f6baf047d42ca46a0123da669a228e6835"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def l_polynomial(a: sp.Expr, b: sp.Expr, c: sp.Expr) -> sp.Expr:
    return sp.expand(27 * a**2 * c**2 - 18 * a * b * c
                     + 16 * a + b**3 * c - b**2)


def sporadic_map(x: sp.Expr, y: sp.Expr, z: sp.Expr) -> tuple[sp.Expr, ...]:
    u = 1 + x * y
    return (
        sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y)),
        sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)),
        sp.expand(2 * x - 3 * x**2 * y - x**3 * z),
    )


def inverse_core(target: tuple[sp.Expr, ...], variable: sp.Symbol) -> sp.Expr:
    a, b, c = target
    return sp.expand(l_polynomial(a, b, c) * variable**3
                     + (4 - 3 * b * c) * variable - 2 * c)


def inverse_section(target: tuple[sp.Expr, ...], w: sp.Symbol
                    ) -> tuple[sp.Expr, ...]:
    a, b, c = target
    denominator = sp.expand((12 * a - b**2) * w**2 + b * w + 2)
    y = sp.cancel(b - 3 * a * w * ((9 * a * c - b) * w + 2)
                  / denominator)
    z = sp.cancel((2 * w - c - 3 * w**2 * y) / w**3)
    return w, y, z


def primitive_polynomial(expression: sp.Expr, variable: sp.Symbol) -> sp.Poly:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    require(variable not in denominator.free_symbols,
            ("nonconstant univariate denominator", denominator))
    polynomial = sp.Poly(numerator, variable, domain=sp.QQ)
    _, primitive = polynomial.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive


def polynomial_ledger_hash(polynomial: sp.Poly) -> str:
    ledger = tuple(
        (tuple(int(exponent) for exponent in monomial), str(coefficient))
        for monomial, coefficient in polynomial.terms()
    )
    return sha256(
        json.dumps(ledger, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def normalized_factors(expression: sp.Expr, variable: sp.Symbol
                       ) -> tuple[tuple[sp.Poly, int], ...]:
    _, raw_factors = sp.factor_list(expression, variable)
    factors = tuple(
        (sp.Poly(factor, variable, domain=sp.QQ).monic(), exponent)
        for factor, exponent in raw_factors
    )
    return tuple(sorted(
        factors,
        key=lambda item: (item[0].degree(), polynomial_ledger_hash(item[0]),
                          item[1]),
    ))


def subset_sums(degrees: tuple[int, ...]) -> set[int]:
    sums = {0}
    for degree in degrees:
        sums |= {value + degree for value in tuple(sums)}
    return sums


def irreducibility_certificate(polynomial: sp.Poly
                               ) -> tuple[tuple[tuple[int, tuple[int, ...]], ...],
                                          tuple[int, ...]]:
    """Record modular factor-degree restrictions on possible Q factors.

    A rational factor of degree d reduces, at every good prime, to a union of
    irreducible modular factors.  Their subset-degree intersections can prove
    irreducibility without an n-cycle.  An imprimitive transitive action can
    leave a proper degree forever possible, so exact QQ irreducibility is also
    checked by the caller rather than misreporting that situation as failure.
    """
    _, integral = polynomial.clear_denoms(convert=True)
    _, integral = integral.primitive()
    degree = integral.degree()
    possible_degrees = set(range(degree + 1))
    certificate = []
    for prime in sp.primerange(2, 2000):
        reduced = sp.Poly(integral.as_expr(), *integral.gens, modulus=prime)
        if reduced.degree() != degree or sp.gcd(reduced, reduced.diff()).degree() != 0:
            continue
        coefficients = [int(coefficient) % prime
                        for coefficient in integral.all_coeffs()]
        _, modular_factors = gf_factor(coefficients, prime, sp.ZZ)
        if any(exponent != 1 for _, exponent in modular_factors):
            continue
        factor_degrees = tuple(sorted(
            len(factor) - 1
            for factor, _ in modular_factors
        ))
        refined = possible_degrees & subset_sums(factor_degrees)
        if refined != possible_degrees:
            certificate.append((int(prime), factor_degrees))
            possible_degrees = refined
        if possible_degrees == {0, degree}:
            break
    return tuple(certificate), tuple(sorted(possible_degrees))


def unit_mod_outer(rational: sp.Expr, outer: sp.Poly, w: sp.Symbol,
                   label: str) -> None:
    numerator, denominator = sp.fraction(sp.cancel(rational))
    for part_name, part in (("numerator", numerator),
                            ("denominator", denominator)):
        part_polynomial = sp.Poly(part, w, domain=sp.QQ)
        require(sp.gcd(outer, part_polynomial).degree() == 0,
                (label, part_name, part))


def zero_mod_outer(rational: sp.Expr, outer: sp.Poly, w: sp.Symbol,
                   label: str) -> None:
    numerator, denominator = sp.fraction(sp.cancel(rational))
    require(sp.gcd(outer, sp.Poly(denominator, w,
                                  domain=sp.QQ)).degree() == 0,
            (label, "denominator", denominator))
    remainder = sp.rem(sp.Poly(numerator, w, domain=sp.QQ), outer)
    require(remainder.is_zero, (label, remainder))


def orbit_partition(items: tuple[object, ...], actions: tuple[dict, ...]
                    ) -> tuple[tuple[object, ...], ...]:
    universe = set(items)
    unseen = set(items)
    orbits = []
    while unseen:
        seed = min(unseen, key=repr)
        orbit = {action[seed] for action in actions}
        require(orbit <= universe, ("orbit escaped universe", seed, orbit))
        unseen -= orbit
        orbits.append(tuple(sorted(orbit, key=repr)))
    return tuple(sorted(orbits, key=lambda orbit: (len(orbit), repr(orbit))))


def marked_depth_two_actions() -> tuple[tuple[dict, ...], tuple[dict, ...]]:
    letters = (0, 1, 2)
    words = tuple(product(letters, repeat=2))
    letter_permutations = tuple(permutations(letters))
    word_actions = []
    pair_actions = []
    pairs = tuple((left, right) for index, left in enumerate(words)
                  for right in words[index + 1:])
    for root_permutation in letter_permutations:
        if root_permutation[0] != 0:
            continue
        for local_zero in letter_permutations:
            if local_zero[0] != 0:
                continue
            for local_one in letter_permutations:
                for local_two in letter_permutations:
                    local = (local_zero, local_one, local_two)
                    action = {
                        word: (root_permutation[word[0]],
                               local[word[0]][word[1]])
                        for word in words
                    }
                    require(action[(0, 0)] == (0, 0), action[(0, 0)])
                    word_actions.append(action)
                    pair_actions.append({
                        pair: tuple(sorted((action[pair[0]], action[pair[1]])))
                        for pair in pairs
                    })
    require(len(word_actions) == 144, len(word_actions))
    return tuple(word_actions), tuple(pair_actions)


def main() -> None:
    w, x, t = sp.symbols("w x t")

    # THM-2570 normalization at (tau,lambda)=(0,1).
    q0 = (sp.Rational(1, 9), sp.Rational(4, 3), sp.Rational(0))
    require(l_polynomial(*q0) == 0, ("q0 is not on old L", q0))
    q1 = tuple(sp.factor(value) for value in sporadic_map(*q0))
    q2 = tuple(sp.factor(value) for value in sporadic_map(*q1))
    require(l_polynomial(*q1) != 0, "q1 unexpectedly lies on old L")
    require(l_polynomial(*q2) != 0, "q2 unexpectedly lies on old L")

    # First inverse step: three actual finite middle points.
    outer = sp.Poly(inverse_core(q2, w), w, domain=sp.QQ)
    require(outer.degree() == 3, outer.degree())
    require(sp.gcd(outer, outer.diff()).degree() == 0,
            "outer inverse cubic is not squarefree")
    middle = inverse_section(q2, w)
    require(tuple(sp.factor(value.subs(w, q1[0]) - expected)
                  for value, expected in zip(middle, q1)) == (0, 0, 0),
            "inverse section misses the marked middle point")
    for index, coordinate in enumerate(middle):
        _, denominator = sp.fraction(sp.cancel(coordinate))
        require(sp.gcd(outer, sp.Poly(denominator, w,
                                      domain=sp.QQ)).degree() == 0,
                ("middle section denominator", index, denominator))
    for index, (image_coordinate, target_coordinate) in enumerate(
            zip(sporadic_map(*middle), q2)):
        zero_mod_outer(image_coordinate - target_coordinate, outer, w,
                       f"middle reconstruction coordinate {index}")

    # Second inverse step: the child core remains cubic and reconstructible
    # over all three middle points.
    middle_l = l_polynomial(*middle)
    unit_mod_outer(middle_l, outer, w, "middle L")
    child_rational = inverse_core(middle, x)
    child_numerator, child_denominator = sp.fraction(sp.cancel(child_rational))
    require(sp.gcd(outer, sp.Poly(child_denominator, w,
                                  domain=sp.QQ)).degree() == 0,
            "child-core denominator meets outer roots")
    child = sp.Poly(child_numerator, w, domain=sp.QQ.frac_field(x))

    child_section_denominator = ((12 * middle[0] - middle[1]**2) * x**2
                                 + middle[1] * x + 2)
    section_denominator_numerator, _ = sp.fraction(
        sp.cancel(child_section_denominator)
    )
    reconstruction_obstruction = sp.resultant(
        child_numerator, section_denominator_numerator, x
    )
    unit_mod_outer(reconstruction_obstruction, outer, w,
                   "child reconstruction denominator resultant")

    coordinate_resultant = sp.resultant(outer.as_expr(), child.as_expr(), w)
    point_polynomial = primitive_polynomial(coordinate_resultant, x)
    require(point_polynomial.degree() == 9, point_polynomial.degree())
    require(point_polynomial.TC() != 0, "a specialized child x-coordinate is zero")
    require(sp.gcd(point_polynomial, point_polynomial.diff()).degree() == 0,
            "point eliminant is not squarefree")
    require(point_polynomial.eval(q0[0]) == 0,
            "marked ancestry x-coordinate is absent")

    point_factors = normalized_factors(point_polynomial.as_expr(), x)
    point_factor_data = tuple((factor.degree(), exponent)
                              for factor, exponent in point_factors)
    require(point_factor_data == ((1, 1), (2, 1), (6, 1)),
            ("point factorization", point_factor_data))
    point_product = sp.Poly(sp.prod(
        factor.as_expr()**exponent for factor, exponent in point_factors
    ), x, domain=sp.QQ)
    require(point_product == point_polynomial.monic(),
            "point factors do not reconstruct the monic eliminant")

    marked_child = primitive_polynomial(inverse_core(q1, x), x).monic()
    require(marked_child.degree() == 3, marked_child.degree())
    marked_quotient, marked_remainder = sp.div(point_polynomial.monic(),
                                               marked_child)
    require(marked_remainder.is_zero and marked_quotient.degree() == 6,
            ("marked child division", marked_remainder,
             marked_quotient.degree()))
    marked_factor_data = tuple((factor.degree(), exponent) for factor, exponent
                               in normalized_factors(marked_child.as_expr(), x))
    require(marked_factor_data == ((1, 1), (2, 1)), marked_factor_data)

    # For monic P with roots r_i, Res_x(P(x),P(t-x)) differs only by a sign
    # from prod_(i,j)(t-r_i-r_j).  Removing the diagonal prod_i(t-2r_i)
    # leaves the square of the 36 unordered pair sums.
    monic_point = point_polynomial.monic()
    symmetric_resultant = sp.Poly(
        sp.resultant(monic_point.as_expr(),
                     monic_point.as_expr().subs(x, t - x), x),
        t,
        domain=sp.QQ,
    )
    diagonal = sp.Poly(2**9 * monic_point.as_expr().subs(x, t / 2),
                       t, domain=sp.QQ)
    pair_square, remainder = sp.div(symmetric_resultant, diagonal)
    require(remainder.is_zero, ("diagonal division", remainder))
    require(pair_square.degree() == 72, pair_square.degree())
    pair_square_factors = normalized_factors(pair_square.as_expr(), t)
    pair_square_data = tuple((factor.degree(), exponent)
                             for factor, exponent in pair_square_factors)
    require(pair_square_data == ((1, 2), (2, 2), (6, 2), (6, 2),
                                 (9, 2), (12, 2)),
            ("pair-square factorization", pair_square_data))
    pair_factors = tuple((factor, exponent // 2)
                         for factor, exponent in pair_square_factors)
    pair_polynomial = sp.Poly(sp.prod(
        factor.as_expr()**exponent for factor, exponent in pair_factors
    ), t, domain=sp.QQ)
    require(pair_polynomial.degree() == 36, pair_polynomial.degree())
    require(pair_square.monic() == pair_polynomial**2,
            "symmetric resultant quotient is not the pair-resolvent square")
    require(sp.gcd(pair_polynomial, pair_polynomial.diff()).degree() == 0,
            "unordered pair sums collide")
    pair_factor_data = tuple((factor.degree(), exponent)
                             for factor, exponent in pair_factors)
    require(pair_factor_data == ((1, 1), (2, 1), (6, 1), (6, 1),
                                 (9, 1), (12, 1)), pair_factor_data)

    # Exact QQ irreducibility is paired with modular subset-degree fingerprints.
    # The latter deliberately need not close for an imprimitive transitive
    # action; this avoids falsely demanding one full cycle.
    require(all(factor.is_irreducible for factor, _ in point_factors),
            "a reported point factor is reducible over QQ")
    require(all(factor.is_irreducible for factor, _ in pair_factors),
            "a reported pair factor is reducible over QQ")
    point_certificates = tuple(
        (factor.degree(), polynomial_ledger_hash(factor),
         irreducibility_certificate(factor))
        for factor, _ in point_factors
    )
    pair_certificates = tuple(
        (factor.degree(), polynomial_ledger_hash(factor),
         irreducibility_certificate(factor))
        for factor, _ in pair_factors
    )

    # Full H_2=Stab_W2(00) census, independently enumerated from its 144
    # rooted-tree automorphisms.
    words = tuple(product((0, 1, 2), repeat=2))
    pairs = tuple((left, right) for index, left in enumerate(words)
                  for right in words[index + 1:])
    word_actions, pair_actions = marked_depth_two_actions()
    point_orbits = orbit_partition(words, word_actions)
    pair_orbits = orbit_partition(pairs, pair_actions)
    point_orbit_sizes = tuple(sorted(map(len, point_orbits)))
    pair_orbit_sizes = tuple(sorted(map(len, pair_orbits)))
    require(point_orbit_sizes == (1, 2, 6), point_orbit_sizes)
    require(pair_orbit_sizes == (1, 2, 6, 6, 9, 12), pair_orbit_sizes)
    require(sum(point_orbit_sizes) == 9, point_orbit_sizes)
    require(sum(pair_orbit_sizes) == 36, pair_orbit_sizes)
    require(len(point_orbits) + len(pair_orbits) == 9,
            (point_orbits, pair_orbits))

    ledger = {
        "specialization": {
            "tau": 0,
            "lambda": 1,
            "q0": tuple(map(str, q0)),
            "q1": tuple(map(str, q1)),
            "q2": tuple(map(str, q2)),
            "L_q1": str(sp.factor(l_polynomial(*q1))),
            "L_q2": str(sp.factor(l_polynomial(*q2))),
        },
        "point_eliminant": {
            "degree": point_polynomial.degree(),
            "terms": len(point_polynomial.terms()),
            "ledger_sha256": polynomial_ledger_hash(point_polynomial),
            "factors": point_certificates,
            "marked_child_factors": marked_factor_data,
        },
        "pair_resolvent": {
            "degree": pair_polynomial.degree(),
            "terms": len(pair_polynomial.terms()),
            "ledger_sha256": polynomial_ledger_hash(pair_polynomial),
            "factors": pair_certificates,
            "squarefree": True,
        },
        "marked_tree_stabilizer": {
            "order": len(word_actions),
            "point_orbit_sizes": point_orbit_sizes,
            "pair_orbit_sizes": pair_orbit_sizes,
            "packet_count": len(point_orbits) + len(pair_orbits),
        },
    }
    semantic_sha256 = sha256(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    print("== THM-3542 depth-three newest-prime point/pair saturation ==")
    print("specialization=(tau,lambda)=(0,1);q0=(1/9,4/3,0);q2=F^2(q0)")
    print(f"L(q1)={sp.factor(l_polynomial(*q1))};"
          f"L(q2)={sp.factor(l_polynomial(*q2))}")
    print("fibre_gates=outer_squarefree;all reconstruction denominators units;"
          "degree9;point_x_squarefree;marked_child=1+2;PASS")
    print(f"point_factor_degrees={point_factor_data};"
          f"point_modular_certificates={tuple(row[2] for row in point_certificates)}")
    print(f"pair_factor_degrees={pair_factor_data};pair_sums_squarefree=PASS;"
          f"pair_modular_certificates={tuple(row[2] for row in pair_certificates)}")
    print(f"H2_order={len(word_actions)};point_orbits={point_orbit_sizes};"
          f"pair_orbits={pair_orbit_sizes};packets=9=3^2")
    print("specialized_orbits=H2_orbits;therefore generic residue point/pair "
          "orbits=H2_orbits by good-specialization refinement")
    print(f"point_ledger_sha256={polynomial_ledger_hash(point_polynomial)}")
    print(f"pair_ledger_sha256={polynomial_ledger_hash(pair_polynomial)}")
    print(f"semantic_sha256={semantic_sha256}")
    print("status=PROVED FINITE-EXACT DEPTH-THREE POINT/PAIR SATURATION;"
          "no full decomposition group, n>=4, coordinate-index, or LRC claim")


if __name__ == "__main__":
    main()

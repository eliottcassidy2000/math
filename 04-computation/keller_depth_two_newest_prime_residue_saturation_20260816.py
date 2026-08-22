#!/usr/bin/env python3
"""Exact depth-two newest-prime residue-action certificate for THM-3540.

The script pulls the inverse cubic back to the normalization of V(L), removes
the birational ancestry root, and proves that the residual quadratic has
square class -L(target).  Its restriction to the lambda=0 axis has odd tau
valuation, which is the decisive nonsquare certificate.
"""

from __future__ import annotations

from hashlib import sha256
import json

import sympy as sp


EXPECTED_SEMANTIC_SHA256 = "81e36ca0659eb75dc7a0f660f4c8a40ae70bf0d4fba56c08f4ef17241349488e"


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


def primitive_numerator(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Poly:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    require(not any(variable in denominator.free_symbols for variable in variables),
            ("nonconstant denominator", expression, denominator))
    polynomial = sp.Poly(numerator, *variables, domain=sp.QQ)
    _, primitive = polynomial.primitive()
    return primitive


def polynomial_ledger_hash(polynomial: sp.Poly) -> str:
    ledger = tuple(
        (tuple(int(exponent) for exponent in monomial), str(coefficient))
        for monomial, coefficient in polynomial.terms()
    )
    return sha256(
        json.dumps(ledger, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def orbit_partition(items: tuple[object, ...], permutations: tuple[dict[int, int], ...],
                    action) -> tuple[tuple[object, ...], ...]:
    unseen = set(items)
    orbits = []
    while unseen:
        seed = min(unseen, key=repr)
        orbit = {action(seed, permutation) for permutation in permutations}
        require(orbit <= set(items), ("orbit escaped universe", seed, orbit))
        unseen -= orbit
        orbits.append(tuple(sorted(orbit, key=repr)))
    return tuple(sorted(orbits, key=lambda orbit: (len(orbit), repr(orbit))))


def main() -> None:
    tau, lam, w = sp.symbols("tau lambda w")

    # THM-2570 normalization of the old Jelonek component.
    source_x = lam**2 * (3 - tau * lam) / 27
    source_y = lam * (4 - tau * lam) / 3
    source_z = tau
    require(sp.factor(l_polynomial(source_x, source_y, source_z)) == 0,
            "normalization does not lie on V(L)")

    target_a, target_b, target_c = sporadic_map(source_x, source_y, source_z)
    target_l = sp.factor(l_polynomial(target_a, target_b, target_c))
    target_t = sp.factor(4 - 3 * target_b * target_c)
    target_s = sp.factor(27 * target_a * target_c**2
                         - 9 * target_b * target_c + 8)
    require(target_l != 0, "image H is not generically contained in old L")

    # THM-2473 inverse x-core over the generic point of H.
    core = sp.Poly(target_l * w**3 + target_t * w - 2 * target_c,
                   w, domain=sp.QQ.frac_field(tau, lam))
    require(sp.factor(core.eval(source_x)) == 0,
            "birational ancestry x-coordinate is not an inverse root")
    quotient, remainder = sp.div(core, sp.Poly(w - source_x, w,
                                               domain=core.domain))
    require(remainder.is_zero, ("inverse division remainder", remainder))
    expected_quotient = sp.Poly(
        target_l * w**2 + target_l * source_x * w
        + target_l * source_x**2 + target_t,
        w,
        domain=core.domain,
    )
    require(quotient == expected_quotient,
            ("residual quadratic", quotient, expected_quotient))

    quadratic_disc = sp.factor(sp.discriminant(quotient.as_expr(), w))
    core_disc = sp.factor(sp.discriminant(core.as_expr(), w))
    derivative_at_root = sp.factor(sp.diff(core.as_expr(), w).subs(w, source_x))
    require(sp.factor(core_disc + 4 * target_s**2 * target_l) == 0,
            "cubic discriminant identity")
    require(sp.factor(core_disc - quadratic_disc * derivative_at_root**2) == 0,
            "linear-times-quadratic discriminant product")
    require(sp.factor(
        quadratic_disc * derivative_at_root**2
        + 4 * target_s**2 * target_l
    ) == 0, "residual square-class identity")
    square_ratio = sp.cancel(quadratic_disc / (-target_l))
    require(sp.factor(square_ratio - (2 * target_s / derivative_at_root)**2) == 0,
            "quadratic discriminant / -L is not the advertised square")

    # The lambda=0 divisor is the axis q=(0,0,tau), F(q)=(tau,0,0).
    axis_substitution = {lam: 0}
    axis_source = tuple(sp.factor(value.subs(axis_substitution))
                        for value in (source_x, source_y, source_z))
    axis_target = tuple(sp.factor(value.subs(axis_substitution))
                        for value in (target_a, target_b, target_c))
    axis_l = sp.factor(target_l.subs(axis_substitution))
    axis_s = sp.factor(target_s.subs(axis_substitution))
    axis_derivative = sp.factor(derivative_at_root.subs(axis_substitution))
    axis_core = sp.factor(core.as_expr().subs(axis_substitution))
    axis_quotient = sp.factor(quotient.as_expr().subs(axis_substitution))
    axis_quadratic_disc = sp.factor(quadratic_disc.subs(axis_substitution))
    require(axis_source == (0, 0, tau), axis_source)
    require(axis_target == (tau, 0, 0), axis_target)
    require(axis_l == 16 * tau, axis_l)
    require(axis_s == 8, axis_s)
    require(axis_derivative == 4, axis_derivative)
    require(sp.expand(axis_core - 4 * w * (4 * tau * w**2 + 1)) == 0,
            axis_core)
    require(sp.expand(axis_quotient - 4 * (4 * tau * w**2 + 1)) == 0,
            axis_quotient)
    require(axis_quadratic_disc == -256 * tau, axis_quadratic_disc)

    # In the lambda-adic DVR, a hypothetical square root is a unit because
    # -L(F nu) has lambda-valuation zero.  Its residue would square to
    # -16*tau, impossible because the tau-valuation is one.
    minus_l_axis = sp.Poly(-axis_l, tau, domain=sp.QQ)
    tau_valuation = min(monomial[0] for monomial, _ in minus_l_axis.terms())
    require(tau_valuation == 1, ("axis odd valuation", minus_l_axis))

    # Full H_1 point/pair orbit census: marked block 0 and one off-ray pair.
    permutations = (
        {0: 0, 1: 1, 2: 2},
        {0: 0, 1: 2, 2: 1},
    )
    blocks = (0, 1, 2)
    pairs = ((0, 1), (0, 2), (1, 2))
    block_orbits = orbit_partition(
        blocks, permutations, lambda vertex, permutation: permutation[vertex]
    )
    pair_orbits = orbit_partition(
        pairs,
        permutations,
        lambda pair, permutation: tuple(sorted((permutation[pair[0]],
                                                 permutation[pair[1]]))),
    )
    require(tuple(sorted(map(len, block_orbits))) == (1, 2), block_orbits)
    require(tuple(sorted(map(len, pair_orbits))) == (1, 2), pair_orbits)
    require(len(block_orbits) + len(pair_orbits) == 4,
            (block_orbits, pair_orbits))

    l_pullback = primitive_numerator(target_l, (tau, lam))
    qdisc_pullback = primitive_numerator(quadratic_disc, (tau, lam))
    ledger = {
        "normalization": {
            "x": str(source_x),
            "y": str(source_y),
            "z": str(source_z),
        },
        "image_L": {
            "degrees_tau_lambda": [l_pullback.degree(tau), l_pullback.degree(lam)],
            "terms": len(l_pullback.terms()),
            "ledger_sha256": polynomial_ledger_hash(l_pullback),
        },
        "quadratic_discriminant": {
            "degrees_tau_lambda": [qdisc_pullback.degree(tau),
                                    qdisc_pullback.degree(lam)],
            "terms": len(qdisc_pullback.terms()),
            "ledger_sha256": polynomial_ledger_hash(qdisc_pullback),
            "square_class": "-L(target)",
        },
        "axis": {
            "source": tuple(map(str, axis_source)),
            "target": tuple(map(str, axis_target)),
            "L": str(axis_l),
            "quadratic_discriminant": str(axis_quadratic_disc),
            "tau_valuation_of_minus_L": tau_valuation,
        },
        "orbits": {
            "blocks": block_orbits,
            "pairs": pair_orbits,
            "packet_count": 4,
        },
    }
    semantic_sha256 = sha256(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    print("== THM-3540 depth-two newest-prime residue saturation ==")
    print("inverse_factorization=one birational V(L) root + residual quadratic")
    print("quadratic_disc_square_class=[-L(target)]")
    print("axis=lambda=0;source=(0,0,tau);target=(tau,0,0);"
          "L=16*tau;quadratic_disc=-256*tau;tau_valuation=1")
    print("residual_quadratic=irreducible over kappa(H);"
          "residue_action_on_predecessor_blocks=1+2=S2-stabilizer")
    print(f"block_orbits={block_orbits};pair_orbits={pair_orbits};"
          "THM3538_packets=4=2^2")
    print("index_packet_formula=v(A_marked)+2v(A_off)+"
          "2v(R_marked_off)+v(R_off_off)")
    print(f"image_L_ledger_sha256={polynomial_ledger_hash(l_pullback)}")
    print(f"quadratic_disc_ledger_sha256={polynomial_ledger_hash(qdisc_pullback)}")
    print(f"semantic_sha256={semantic_sha256}")
    print("status=PROVED SYMBOLIC DEPTH-TWO POINT/PAIR SATURATION;"
          "no full-centralizer, n>=3, all-level-index, or LRC claim")


if __name__ == "__main__":
    main()

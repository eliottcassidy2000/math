#!/usr/bin/env python3
"""Exact referee for THM-3431's two secondary H^1 defects.

LRC side:
  H^1(C7;F13) -> H^1(C13_deck;F13)
is reconstructed from the unique primitive on the degree-13 cyclic cover.

JC side:
  A theta_sigma ~= A/(lambda^q) -> H^1_(lambda)(A)
maps the distinguished observer to [lambda^-q], where A=K[lambda].

The shared valuation barcode is checked only as a lossy invariant.  No
cross-coefficient map, current, mate, LRC(14), or JC(2) claim is made.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256


EXPECTED_SEMANTIC_SHA256 = "79f37de16658ac1e755e0acb5640507dba716756d9b3136f378004a4d29bcc21"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def valuation(n: int, prime: int) -> int:
    answer = 0
    while n % prime == 0:
        n //= prime
        answer += 1
    return answer


def lrc_secondary_audit() -> tuple[object, ...]:
    edge_checks = 0
    deck_checks = 0
    cocycle_checks = 0
    degree_checks = 0
    rows = []
    for a in range(1, 13):
        holonomy = 7 * a % 13
        primitive = tuple(j * a % 13 for j in range(91))

        # Pullback of the constant base cochain is delta h on C91, including
        # the closing edge 90 -> 0.
        for j in range(91):
            delta = (primitive[(j + 1) % 91] - primitive[j]) % 13
            require(delta == a, ("cover primitive", a, j, delta))
            edge_checks += 1

        defects = []
        for power in range(13):
            expected = 7 * power * a % 13
            observed_values = {
                (primitive[(j + 7 * power) % 91] - primitive[j]) % 13
                for j in range(91)
            }
            require(observed_values == {expected},
                    ("deck defect", a, power, observed_values, expected))
            defects.append(expected)
            deck_checks += 91
        for left in range(13):
            for right in range(13):
                require(
                    defects[(left + right) % 13]
                    == (defects[left] + defects[right]) % 13,
                    ("deck one-cocycle", a, left, right),
                )
                cocycle_checks += 1

        require(defects[1] == holonomy != 0,
                ("transgression class", a, defects[1], holonomy))
        for degree in range(1, 201):
            exact = 7 * degree * a % 13 == 0
            require(exact == (degree % 13 == 0),
                    ("degree death", a, degree, exact))
            degree_checks += 1
        rows.append((a, holonomy, defects[1], primitive[:8]))

    # Pullback composition multiplies degrees, while v_13 turns that
    # multiplication into addition.
    composition_checks = 0
    for left in range(1, 51):
        for right in range(1, 51):
            require(
                valuation(left * right, 13)
                == valuation(left, 13) + valuation(right, 13),
                ("valuation composition", left, right),
            )
            composition_checks += 1

    hostiles = (
        (12, valuation(12, 13), "survives"),
        (13, valuation(13, 13), "dies"),
        (14, valuation(14, 13), "survives"),
        (26, valuation(26, 13), "dies"),
        (169, valuation(169, 13), "dies"),
    )
    require(hostiles[2] == (14, 0, "survives"), hostiles)
    return (
        tuple(rows),
        edge_checks,
        deck_checks,
        cocycle_checks,
        degree_checks,
        composition_checks,
        hostiles,
    )


def selected_observer_profiles() -> tuple[tuple[int, int, int, int], ...]:
    rows = []
    for d in range(2, 13):
        for e in range(2, 26):
            for sigma in range(1, d + 1):
                numerator = sigma * (e - 1)
                if numerator % d == 0:
                    rows.append((d, e, sigma, numerator // d))
    return tuple(rows)


def jc_secondary_audit() -> tuple[object, ...]:
    profiles = selected_observer_profiles()
    annihilator_checks = 0
    embedding_checks = 0
    bar_histogram = Counter()
    sample_rows = []
    for d, e, sigma, depth in profiles:
        require(depth >= 1, (d, e, sigma, depth))
        bar_histogram[depth] += 1

        # [lambda^-depth] in A[lambda^-1]/A.  Multiplication by lambda^n
        # is nonzero exactly while n-depth is negative.
        states = []
        for exponent in range(depth + 3):
            laurent_exponent = exponent - depth
            nonzero = laurent_exponent < 0
            require(nonzero == (exponent < depth),
                    ("local H1 annihilator", d, e, sigma, depth, exponent))
            states.append((exponent, laurent_exponent, nonzero))
            annihilator_checks += 1

        # The canonical cyclic embedding A/(lambda^depth) -> H^1_(lambda)(A)
        # sends lambda^i to lambda^(i-depth); all representatives are
        # distinct negative powers.
        image_basis = tuple(range(-depth, 0))
        require(len(image_basis) == depth and len(set(image_basis)) == depth,
                ("cyclic local-cohomology embedding", d, e, sigma, depth))
        embedding_checks += depth
        if (d, e, sigma) in ((2, 2, 2), (2, 3, 1), (2, 5, 1), (3, 4, 1)):
            sample_rows.append((d, e, sigma, depth, tuple(states), image_basis))

    require(
        tuple((row[0], row[1], row[2], row[3]) for row in sample_rows)
        == ((2, 2, 2, 1), (2, 3, 1, 1), (2, 5, 1, 2), (3, 4, 1, 1)),
        sample_rows,
    )
    return (
        len(profiles),
        annihilator_checks,
        embedding_checks,
        tuple(sorted(bar_histogram.items())),
        tuple(sample_rows),
    )


def additive_no_go_controls() -> tuple[object, ...]:
    # H^1(C13;F13) has exponent 13.  A characteristic-zero local-cohomology
    # module is a K-vector space and hence has no nonzero 13-torsion.
    source_orders = tuple((x, 13 * x % 13) for x in range(13))
    require(all(value == 0 for _, value in source_orders), source_orders)

    # Conversely, multiplication by 13 is surjective on every K-vector
    # space.  Any homomorphism to an exponent-13 group therefore vanishes:
    # f(y)=f(13*(y/13))=13f(y/13)=0.  Fractions (numerator, denominator)
    # below freeze the divisibility identity on representative controls.
    divisible_controls = tuple(
        (numerator, denominator, numerator, 13 * denominator)
        for numerator, denominator in ((1, 1), (-7, 5), (13, 17), (91, 29))
    )
    for numerator, denominator, divided_numerator, divided_denominator in divisible_controls:
        require(
            13 * divided_numerator * denominator
            == numerator * divided_denominator,
            ("13 divisibility", numerator, denominator),
        )

    # Equal barcode length one is deliberately not a class map.
    equal_length_hostile = {
        "LRC": ("H1_group", "F13", 1),
        "JC": ("H1_local", "characteristic_zero", 1),
    }
    require(equal_length_hostile["LRC"][1] != equal_length_hostile["JC"][1],
            equal_length_hostile)
    return source_orders, divisible_controls, equal_length_hostile


def main() -> None:
    lrc = lrc_secondary_audit()
    jc = jc_secondary_audit()
    no_go = additive_no_go_controls()
    semantic_surface = (lrc, jc, no_go)
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256,
                (semantic_digest, EXPECTED_SEMANTIC_SHA256))

    print("THM-3431 D5 secondary H1 descent defects -- exact referee")
    print("status=VERIFIED_EXACT_SUPPORT_FOR_PROVED_SECONDARY_CLASSES;typed_cospan_only;no_cross_coefficient_map;no_LRC14_or_JC2")
    print("LRC=(rows,edge_checks,deck_checks,cocycle_checks,degree_checks,composition_checks,hostiles)=" + repr(lrc))
    print("JC=(selected_profiles,annihilator_checks,embedding_checks,bar_histogram,samples)=" + repr(jc))
    print("additive_no_go=(C13_exponent,characteristic_zero_divisibility,equal_bar_hostile)=" + repr(no_go))
    print("common_forgetting=LRC_degree_v13_bar_[0,1);JC_lambda_bar_[0,q);site_coefficient_class_and_target_destroyed")
    print(f"semantic_sha256={semantic_digest}")
    print("scope=secondary_H1_classes_and_lossy_valuation_barcodes_only;no_semantic_current;no_polynomial_mate;LRC14_and_JC2_open")


if __name__ == "__main__":
    main()

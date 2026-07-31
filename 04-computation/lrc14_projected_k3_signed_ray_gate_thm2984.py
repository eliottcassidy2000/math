#!/usr/bin/env python3
"""Independent finite hostile audit for the new THM-2984 residue claims.

The proofs in THM-2984 are uniform.  This script is only a counterexample
search and boundary-control companion; it is not a dependency of the theorem.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import random
from fractions import Fraction


def require(condition: bool, witness: object) -> None:
    """Optimization-stable proof gate (`python -O` must not erase checks)."""
    if not condition:
        raise RuntimeError(f"hostile audit failed: {witness!r}")


def units(d: int) -> tuple[int, ...]:
    return tuple(u for u in range(d) if math.gcd(u, d) == 1)


def band(d: int) -> frozenset[int]:
    return frozenset(r for r in range(d) if 14 * min(r, d - r) < d)


def beta(d: int) -> int:
    return 2 * ((d - 1) // 14) + 1


def order_mod(x: int, d: int) -> int:
    return d // math.gcd(x % d, d)


def transporter(d: int, S: frozenset[int], B: frozenset[int]) -> frozenset[int]:
    return frozenset(u for u in units(d) if {(u * x) % d for x in S} <= B)


def translated_band(d: int, delta: Fraction) -> frozenset[int]:
    """Strict band after scaling theta by d: delta=d*theta on R/dZ."""
    answer = set()
    for r in range(d):
        x = (Fraction(r) + delta) % d
        if 14 * min(x, d - x) < d:
            answer.add(r)
    return frozenset(answer)


def translated_phase_samples(d: int) -> tuple[Fraction, ...]:
    """One exact point in every membership chamber, plus every boundary.

    Integer shifts of delta only rotate residues.  Modulo one, membership can
    change only at the endpoint fractions +-d/14, so this is exhaustive over
    all real translations for cardinality purposes.
    """
    a = Fraction(d % 14, 14)
    events = sorted({Fraction(0), a % 1, (-a) % 1})
    samples = set(events)
    cyclic = events + [events[0] + 1]
    for left, right in zip(cyclic, cyclic[1:]):
        samples.add(((left + right) / 2) % 1)
    return tuple(sorted(samples))


def cyclic_block(d: int, a: int, kappa: int) -> frozenset[int]:
    return frozenset((a + j) % d for j in range(kappa))


def affine_transporter(d: int, S: frozenset[int]) -> frozenset[tuple[int, int]]:
    kappa = (d + 6) // 7
    blocks = tuple(cyclic_block(d, a, kappa) for a in range(d))
    return frozenset(
        (u, a)
        for u in units(d)
        for a, C in enumerate(blocks)
        if {u * x % d for x in S} <= C
    )


def cyclic_gaps(d: int, S: frozenset[int]) -> tuple[int, ...]:
    """Positive cyclic step gaps of a nonempty residue set.

    `S` is a set, so repeated input representatives have already been
    discarded.  For a singleton the unique positive gap is the full turn
    `d`, not zero.
    """
    require(bool(S), ("cyclic-gaps-requires-nonempty-set", d, S))
    xs = sorted(x % d for x in S)
    return tuple(
        (xs[(i + 1) % len(xs)] - xs[i]) % d or d
        for i in range(len(xs))
    )


def cyclic_gap_cover_count(d: int, S: frozenset[int], kappa: int) -> int:
    """Predicted number of kappa-block starts containing nonempty `S`."""
    return sum(max(0, gap - (d - kappa)) for gap in cyclic_gaps(d, S))


def composition_orbit_face(d: int, S: frozenset[int]) -> bool:
    """Literal diagonal-unit/positive-composition normal form (13o)--(13p)."""
    kappa = (d + 6) // 7
    if not S:
        return True
    if len(S) == 1:
        return True
    if len(S) > kappa:
        return False
    for ordering in itertools.permutations(S):
        for u in units(d):
            # Since kappa-1<d, a positive target a_i in the allowed simplex
            # is the unique least-positive residue of this unit-scaled step.
            steps = tuple(
                u * (ordering[i + 1] - ordering[i]) % d
                for i in range(len(ordering) - 1)
            )
            if sum(steps) <= kappa - 1:
                return True
    return False


def prime_triple_ratio_face(p: int, S: frozenset[int]) -> bool:
    """Six-ordering projective ratio test (13r), for a three-set mod prime p."""
    require(len(S) == 3, ("prime-ratio-requires-triple", p, S))
    kappa = (p + 6) // 7
    ratios = frozenset(
        b * pow(a, -1, p) % p
        for a in range(1, kappa)
        for b in range(1, kappa)
        if a + b <= kappa - 1
    )
    return any(
        (z - y) * pow((y - x) % p, -1, p) % p in ratios
        for x, y, z in itertools.permutations(S)
    )


def check_modulus(d: int) -> tuple[object, ...]:
    B = band(d)
    U = units(d)
    b = (d - 1) // 14
    require(len(B) == beta(d), ("band-cardinality", d, B, beta(d)))
    ceil_seventh = (d + 6) // 7
    require(beta(d) <= ceil_seventh, ("beta-vs-ceil", d, beta(d), ceil_seventh))
    R = max(r for r in range(1, 8) if d % r == 0)
    require(ceil_seventh <= d // R, ("ceil-vs-pair", d, R, ceil_seventh))

    stabilizer = tuple(u for u in U if frozenset(u * x % d for x in B) == B)
    expected_stabilizer = U if d <= 14 else (1, d - 1)
    require(stabilizer == expected_stabilizer, ("stabilizer", d, stabilizer))

    divisors = tuple(g for g in range(1, d + 1) if d % g == 0)
    capacities = []
    for g in divisors:
        actual = sum(math.gcd(r, d) == g for r in B)
        expected = (
            1
            if g == d
            else 2
            * sum(
                math.gcd(a, d // g) == 1
                for a in range(1, b // g + 1)
            )
        )
        require(actual == expected, ("stratum-capacity", d, g, actual, expected))
        capacities.append((g, actual))

    for x, y in itertools.combinations(B, 2):
        require(order_mod(x - y, d) > 7, ("band-independence", d, x, y))

    # Every short-order pair escapes the band for every unit.
    for delta in range(1, d):
        r = order_mod(delta, d)
        if 2 <= r <= 7:
            require(
                r * min(delta, d - delta) >= d,
                ("order-separation", d, delta, r),
            )
            if d <= 100:
                for x in range(d):
                    y = (x + delta) % d
                    for u in U:
                        require(
                            not (u * x % d in B and u * y % d in B),
                            ("pair-to-unit", d, x, y, r, u),
                        )

    return d, beta(d), R, stabilizer, tuple(capacities)


def main() -> None:
    rows = [check_modulus(d) for d in range(2, 501)]
    require(len(rows) == 499, ("formula-modulus-row-count", len(rows)))

    # Exact continuum reduction for the arbitrary translated band.  The
    # sample list contains every boundary and one point per open chamber.
    translated_rows = []
    translated_phase_checks = 0
    affine_containment_checks = 0
    affine_realization_checks = 0
    centered_subcomplex_checks = 0
    for d in range(2, 501):
        samples = translated_phase_samples(d)
        translated_bands = tuple(translated_band(d, delta) for delta in samples)
        sizes = tuple(len(B) for B in translated_bands)
        kappa = (d + 6) // 7
        require(max(sizes) == kappa, ("translated-capacity", d, samples, sizes))
        blocks = tuple(cyclic_block(d, a, kappa) for a in range(d))
        require(all(len(C) == kappa for C in blocks), ("affine-facet-size", d, kappa))
        contain_witnesses = []
        for delta, B in zip(samples, translated_bands):
            witness = next((a for a, C in enumerate(blocks) if B <= C), None)
            require(witness is not None, ("translated-to-block", d, delta, B))
            contain_witnesses.append(witness)
            affine_containment_checks += 1

        # The centered band is a consecutive cyclic block of length beta and
        # therefore lies in an affine facet.
        centered_witness = next((a for a, C in enumerate(blocks) if band(d) <= C), None)
        require(centered_witness is not None, ("centered-subcomplex", d, band(d)))
        centered_subcomplex_checks += 1

        # A centered open interval around the midpoint of C_0 realizes C_0
        # exactly; integer translation realizes every C_a.  Check all starts
        # through d=100 and the first eight thereafter.
        base_delta = -Fraction(kappa - 1, 2)
        starts = range(d) if d <= 100 else range(min(8, d))
        for a in starts:
            actual = translated_band(d, base_delta - a)
            require(actual == blocks[a], ("block-realization", d, a, actual, blocks[a]))
            affine_realization_checks += 1

        sample_keys = tuple((q.numerator, q.denominator) for q in samples)
        translated_rows.append(
            (d, sample_keys, sizes, kappa, tuple(contain_witnesses), centered_witness)
        )
        translated_phase_checks += len(samples)

    # Direct continuous-phase versus finite affine-transporter equivalence on
    # every residue subset for small moduli.  Integer shifts of delta restore
    # every cyclic rotation omitted by the modulo-one chamber quotient.
    affine_subset_equivalence_checks = 0
    cyclic_gap_formula_checks = 0
    for d in range(2, 11):
        phase_bands = frozenset(
            translated_band(d, delta + shift)
            for delta in translated_phase_samples(d)
            for shift in range(d)
        )
        U = units(d)
        for mask in range(1 << d):
            S = frozenset(x for x in range(d) if mask >> x & 1)
            continuous = any(
                {u * x % d for x in S} <= B for u in U for B in phase_bands
            )
            finite = bool(affine_transporter(d, S))
            require(continuous == finite, ("affine-equivalence", d, S))
            affine_subset_equivalence_checks += 1

            # For each fixed unit, compare the proposed cyclic-gap functional
            # with the literal number of translated block starts.  The empty
            # set is deliberately excluded: every start contains it, whereas
            # a nonempty cyclic gap list is the theorem's stated domain.
            if S:
                kappa = (d + 6) // 7
                blocks = tuple(cyclic_block(d, a, kappa) for a in range(d))
                for u in U:
                    uS = frozenset(u * x % d for x in S)
                    literal = sum(uS <= C for C in blocks)
                    predicted = cyclic_gap_cover_count(d, uS, kappa)
                    require(
                        literal == predicted,
                        ("cyclic-gap-cover-count", d, S, u, uS, literal, predicted),
                    )
                    require(
                        bool(literal)
                        == (max(cyclic_gaps(d, uS)) >= d - kappa + 1),
                        ("cyclic-gap-existence", d, S, u, uS, literal),
                    )
                    cyclic_gap_formula_checks += 1

    # Exact affine pair relation and the first affine flag failure.
    affine_pair_checks = 0
    affine_flag_subset_checks = 0
    affine_small_flag_witnesses: dict[int, tuple[int, ...]] = {}
    expected_small_flag_moduli = frozenset({2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 15})
    for d in range(2, 101):
        kappa = (d + 6) // 7
        facets = tuple(cyclic_block(d, a, kappa) for a in range(d))
        # Unit images of the consecutive facets give the full affine complex.
        affine_facets = tuple(
            frozenset(pow(u, -1, d) * x % d for x in C)
            for u in units(d)
            for C in facets
        )
        affine_edge_set = frozenset(
            pair
            for C in affine_facets
            for pair in (frozenset(p) for p in itertools.combinations(C, 2))
        )
        for x, y in itertools.combinations(range(d), 2):
            pair = frozenset((x, y))
            expected = order_mod(x - y, d) > 7
            require(
                (pair in affine_edge_set) == expected,
                ("affine-pair-order", d, x, y, order_mod(x - y, d)),
            )
            require(
                (pair not in affine_edge_set)
                == (2 <= order_mod(x - y, d) <= 7),
                ("affine-pair-nonface", d, x, y, order_mod(x - y, d)),
            )
            affine_pair_checks += 1

        if d <= 15:
            for mask in range(1 << d):
                S = frozenset(x for x in range(d) if mask >> x & 1)
                face = any(S <= C for C in affine_facets)
                clique = all(
                    frozenset((x, y)) in affine_edge_set
                    for x, y in itertools.combinations(S, 2)
                )
                require(not face or clique, ("affine-face-must-be-clique", d, S))
                if clique and not face and d not in affine_small_flag_witnesses:
                    affine_small_flag_witnesses[d] = tuple(sorted(S))
                if d in expected_small_flag_moduli:
                    require(
                        face == clique,
                        ("affine-flag-small", d, S, face, clique),
                    )
                affine_flag_subset_checks += 1

            require(
                (d in expected_small_flag_moduli)
                == (d not in affine_small_flag_witnesses),
                ("small-flag-classification", d, affine_small_flag_witnesses.get(d)),
            )

    d9_triple = frozenset({0, 1, 2})
    d9_facets = tuple(
        frozenset(pow(u, -1, 9) * x % 9 for x in cyclic_block(9, a, 2))
        for u in units(9)
        for a in range(9)
    )
    require(
        all(
            any(pair <= C for C in d9_facets)
            for pair in (frozenset(p) for p in itertools.combinations(d9_triple, 2))
        ),
        "d=9 affine triple must be a one-skeleton clique",
    )
    require(
        not any(d9_triple <= C for C in d9_facets),
        "d=9 affine triple must not be a face",
    )

    # All-cardinality positive-composition/diagonal-unit normal form.  Exhaust
    # every subset up to one beyond facet size through d=18; this includes the
    # empty, singleton, pair, triple, and capacity-overflow boundaries.
    composition_orbit_checks = 0
    for d in range(2, 19):
        kappa = (d + 6) // 7
        for size in range(0, min(d, kappa + 1) + 1):
            for xs in itertools.combinations(range(d), size):
                S = frozenset(xs)
                direct = bool(affine_transporter(d, S))
                normal_form = composition_orbit_face(d, S)
                require(
                    direct == normal_form,
                    ("composition-orbit-normal-form", d, S, direct, normal_form),
                )
                composition_orbit_checks += 1

    # Higher-cardinality structured and hostile controls once kappa exceeds
    # three.  Every short positive composition is a positive control; seeded
    # sets also probe composite gcd strata and capacity overflow.
    rng_composition = random.Random(130013)
    for d in (22, 23, 29, 36):
        kappa = (d + 6) // 7
        controls: set[frozenset[int]] = {frozenset(), frozenset({0})}
        for length in range(1, kappa):
            for cuts in itertools.combinations(range(1, kappa), length):
                steps = tuple(b - a for a, b in zip((0,) + cuts, cuts + (kappa,)))
                # Drop the final slack part: the retained positive steps have
                # total at most kappa-1 and therefore define a face.
                prefixes = tuple(itertools.accumulate(steps[:-1]))
                controls.add(frozenset((0,) + prefixes))
        for _ in range(512):
            size = rng_composition.randrange(0, min(d, kappa + 2) + 1)
            controls.add(frozenset(rng_composition.sample(range(d), size)))
        for S in controls:
            direct = bool(affine_transporter(d, S))
            normal_form = composition_orbit_face(d, S)
            require(
                direct == normal_form,
                ("composition-orbit-higher", d, S, direct, normal_form),
            )
            composition_orbit_checks += 1

    # Prime triple projectivization.  The small primes with kappa<=2 verify
    # that the ratio fan is empty; later primes exercise nontrivial ratio
    # orbits, including possible repetitions among the six orderings.
    prime_ratio_checks = 0
    for p in (5, 7, 11, 13, 17, 23, 43):
        for xs in itertools.combinations(range(p), 3):
            S = frozenset(xs)
            direct = bool(affine_transporter(p, S))
            ratio = prime_triple_ratio_face(p, S)
            require(direct == ratio, ("prime-triple-ratio", p, S, direct, ratio))
            prime_ratio_checks += 1

    # Exact sign-paired unit tables for the six finite witnesses used at the
    # equality boundary d/R=kappa.  Each displayed tuple is the cyclic gap
    # word for the positive representative u; -u reverses the cyclic word and
    # therefore has the same multiset and maximum.
    finite_flag_gap_expected = {
        18: (
            frozenset({0, 1, 5}),
            {1: (1, 4, 13), 5: (5, 2, 11), 7: (7, 10, 1)},
        ),
        21: (
            frozenset({0, 1, 5}),
            {
                1: (1, 4, 16),
                2: (2, 8, 11),
                4: (4, 16, 1),
                5: (4, 1, 16),
                8: (8, 11, 2),
                10: (8, 2, 11),
            },
        ),
        24: (
            frozenset({0, 1, 10}),
            {1: (1, 9, 14), 5: (2, 3, 19), 7: (7, 15, 2), 11: (11, 3, 10)},
        ),
        30: (
            frozenset({0, 1, 8}),
            {1: (1, 7, 22), 7: (7, 19, 4), 11: (11, 17, 2), 13: (13, 1, 16)},
        ),
        36: (
            frozenset({0, 1, 11}),
            {
                1: (1, 10, 25),
                5: (5, 14, 17),
                7: (5, 2, 29),
                11: (11, 2, 23),
                13: (13, 22, 1),
                17: (7, 10, 19),
            },
        ),
        42: (
            frozenset({0, 1, 10}),
            {
                1: (1, 9, 32),
                5: (5, 3, 34),
                11: (11, 15, 16),
                13: (4, 9, 29),
                17: (2, 15, 25),
                19: (19, 3, 20),
            },
        ),
    }
    finite_flag_gap_rows = []
    finite_flag_gap_checks = 0
    for d, (S, expected) in finite_flag_gap_expected.items():
        kappa = (d + 6) // 7
        threshold = d - kappa + 1
        require(
            all(order_mod(x - y, d) > 7 for x, y in itertools.combinations(S, 2)),
            ("finite-flag-witness-clique", d, S),
        )
        require(
            frozenset(units(d))
            == frozenset({v for u in expected for v in (u, d - u)}),
            ("finite-gap-sign-classes", d, expected, units(d)),
        )
        unit_rows = []
        for u, expected_gaps in expected.items():
            positive = cyclic_gaps(d, frozenset(u * x % d for x in S))
            negative = cyclic_gaps(d, frozenset((-u) * x % d for x in S))
            require(positive == expected_gaps, ("finite-gap-word", d, S, u, positive))
            require(
                negative == tuple(reversed(expected_gaps)),
                ("finite-gap-negative-word", d, S, u, negative, expected_gaps),
            )
            require(
                max(positive) < threshold and max(negative) < threshold,
                ("finite-gap-nonface", d, S, u, threshold, positive, negative),
            )
            unit_rows.append((u, positive, max(positive)))
            finite_flag_gap_checks += 2
        require(not composition_orbit_face(d, S), ("finite-composition-nonface", d, S))
        finite_flag_gap_rows.append((d, tuple(sorted(S)), threshold, tuple(unit_rows)))
    require(
        len(finite_flag_gap_rows) == 6,
        ("finite-flag-gap-row-count", len(finite_flag_gap_rows)),
    )

    # The sole infinite equality family d=7m has a symbolic clique/nonfacet
    # witness.  A broad exact replay checks the formulas, but the argument is
    # algebraic: old/new differences have order >7, while N_1=m-3 and
    # N_2=m-2 would force 3v=+-1 and 2v=+-2 in any unit-step facet.
    symbolic_m_checks = 0
    for m in range(4, 1001):
        if m == 6:
            continue
        d = 7 * m
        S = (frozenset(range(m)) - {1}) | {m + 1}
        require(len(S) == m, ("symbolic-m-cardinality", m, S))
        require(
            all(order_mod(h, d) > 7 for h in range(1, m)),
            ("symbolic-m-old-difference", m),
        )
        new_differences = tuple(range(2, m)) + (m + 1,)
        require(
            all(order_mod(h, d) > 7 for h in new_differences),
            ("symbolic-m-new-difference", m, new_differences),
        )
        expected_new_order = m if (m + 1) % 7 == 0 else 7 * m
        require(
            order_mod(m + 1, d) == expected_new_order > 7,
            ("symbolic-m-outlier-order", m, order_mod(m + 1, d)),
        )
        N1 = sum(((x + 1) % d) in S for x in S)
        N2 = sum(((x + 2) % d) in S for x in S)
        require((N1, N2) == (m - 3, m - 2), ("symbolic-m-pair-counts", m, N1, N2))
        for eps1 in (-1, 1):
            for eps2 in (-1, 1):
                obstruction = 2 * eps1 - 6 * eps2
                require(
                    obstruction % d != 0,
                    ("symbolic-m-facet-congruence", m, eps1, eps2, obstruction),
                )
        symbolic_m_checks += 1

    # Check the arithmetic split behind the full classification: outside the
    # displayed equality set, THM-2979 supplies a clique strictly larger than
    # every affine facet.
    equality_modulus_checks = 0
    finite_equalities = frozenset({2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 24, 30, 36})
    for d in range(2, 5001):
        R = max(r for r in range(1, 8) if d % r == 0)
        kappa = (d + 6) // 7
        equality = d // R == kappa
        expected_equality = d % 7 == 0 or d in finite_equalities
        require(
            equality == expected_equality,
            ("flag-equality-modulus-split", d, R, kappa, equality),
        )
        if not equality:
            require(d // R > kappa, ("flag-cardinality-obstruction", d, R, kappa))
        equality_modulus_checks += 1

    # Exact transporter equivalence on every subset through d=12.
    subset_checks = 0
    active_checks = 0
    for d in range(2, 13):
        B = band(d)
        U = units(d)
        for mask in range(1 << d):
            S = frozenset(x for x in range(d) if mask >> x & 1)
            T = transporter(d, S, B)
            universal_escape = all(any(u * x % d not in B for x in S) for u in U)
            require(
                universal_escape == (not T),
                ("transporter-equivalence", d, S, T),
            )
            subset_checks += 1
            for active_mask in range(1 << len(U)):
                active = frozenset(u for j, u in enumerate(U) if active_mask >> j & 1)
                unresolved = bool(T & active)
                direct = any(all(u * x % d in B for x in S) for u in active)
                require(
                    unresolved == direct,
                    ("active-transporter", d, S, T, active),
                )
                active_checks += 1

    # Equality boundary: the first min(8, phi(d)) canonical unit preimages
    # have transporter exactly the stabilizer coset, across d=2..500.
    equality_checks = 0
    for d in range(2, 501):
        B = band(d)
        U = units(d)
        stab = frozenset(u for u in U if frozenset(u * x % d for x in B) == B)
        for u0 in U[:8]:
            inv = pow(u0, -1, d)
            S = frozenset(inv * x % d for x in B)
            T = transporter(d, S, B)
            expected = frozenset(v * u0 % d for v in stab)
            require(T == expected, ("equality-coset", d, u0, T, expected))
            equality_checks += 1

    # Gcd-stratum implication on deterministic hostile random sets.
    rng = random.Random(2984)
    stratum_checks = 0
    for d in range(2, 201):
        B = band(d)
        U = units(d)
        for _ in range(64):
            S = frozenset(r for r in range(d) if rng.randrange(4) == 0)
            for g in (g for g in range(1, d + 1) if d % g == 0):
                O = frozenset(r for r in range(d) if math.gcd(r, d) == g)
                if len(S & O) > len(B & O):
                    require(
                        all(any(u * c % d not in B for c in S) for u in U),
                        ("gcd-stratum-gate", d, g, S),
                    )
                stratum_checks += 1

    # Boundary and hostile controls.
    require(band(14) == frozenset({0}), "strict-open d=14 boundary")
    require(band(15) == frozenset({0, 1, 14}), "first positive band d=15")
    require(bool(transporter(29, band(29), band(29))), "equality hostile control")
    require(
        all(u % 8 not in band(8) for u in units(8)),
        "d=8 singleton strictness control",
    )
    require(
        math.gcd(2, 6) != 1
        and math.gcd(2 * 1 % 6, 6) != math.gcd(1, 6),
        "nonunit gcd-preservation hostile control",
    )
    d28_shifted = translated_band(28, Fraction(-3, 2))
    require(beta(28) == 3, ("d28-beta", beta(28)))
    require((28 + 6) // 7 == 4, "d28-kappa")
    require(
        d28_shifted == frozenset({0, 1, 2, 3}),
        ("d28-translated-hostile", d28_shifted),
    )
    S28 = frozenset({0, 1, 2, 3})
    require(bool(affine_transporter(28, S28)), "d28 affine transporter")
    require(
        not any({u * x % 28 for x in S28} <= band(28) for u in units(28)),
        "d28 centered-vs-affine strictness",
    )
    require(
        len({math.gcd(x, 28) for x in cyclic_block(28, 0, 4)}) > 1,
        "affine blocks must cross centered gcd strata",
    )

    # The d=43 non-flag control: all proper nonempty subsets of S are faces,
    # but the triple is a minimal nonface.
    B43 = band(43)
    S43 = frozenset({1, 2, 15})
    require(
        B43 == frozenset({0, 1, 2, 3, 40, 41, 42}),
        ("d43-band", B43),
    )
    expected_pair_transporters = {
        frozenset({1, 2}): frozenset({1, 42}),
        frozenset({1, 15}): frozenset({3, 40}),
        frozenset({2, 15}): frozenset({20, 23}),
    }
    for pair, expected in expected_pair_transporters.items():
        actual = transporter(43, pair, B43)
        require(actual == expected, ("d43-pair", pair, actual, expected))
    require(not transporter(43, S43, B43), "d43 triple must be a nonface")
    require(
        all(
            order_mod(x - y, 43) == 43
            for x, y in itertools.combinations(S43, 2)
        ),
        "d43 pair orders",
    )

    profile = {
        "formula_stabilizer_capacity_band_moduli": [2, 500],
        "pair_to_unit_full_moduli": [2, 100],
        "translated_capacity_moduli": [2, 500],
        "translated_capacity_reduction": "all boundaries and one point per chamber modulo integer rotation",
        "affine_block_containment_moduli": [2, 500],
        "affine_block_realization_all_starts_moduli": [2, 100],
        "affine_block_realization_later_starts": "first min(8,d)",
        "centered_subcomplex_moduli": [2, 500],
        "affine_all_subset_equivalence_moduli": [2, 10],
        "cyclic_gap_formula_all_nonempty_subsets_and_units_moduli": [2, 10],
        "affine_pair_relation_moduli": [2, 100],
        "affine_flag_all_subset_moduli": [2, 15],
        "affine_flag_expected_small_moduli": sorted(expected_small_flag_moduli),
        "affine_small_nonflag_witnesses": affine_small_flag_witnesses,
        "affine_first_nonflag_control": [9, [0, 1, 2]],
        "composition_orbit_exhaustive_through_modulus": 18,
        "composition_orbit_exhaustive_cardinalities": "0..kappa+1",
        "composition_orbit_higher_control_moduli": [22, 23, 29, 36],
        "composition_orbit_random_sets_per_higher_modulus": 512,
        "composition_orbit_random_seed": 130013,
        "prime_triple_ratio_exhaustive_primes": [5, 7, 11, 13, 17, 23, 43],
        "finite_flag_gap_table_moduli": sorted(finite_flag_gap_expected),
        "symbolic_multiple_seven_m_range": [4, 1000],
        "symbolic_multiple_seven_exception": 6,
        "flag_equality_modulus_split_range": [2, 5000],
        "all_subset_transporter_moduli": [2, 12],
        "all_active_subsets_per_subset": True,
        "equality_coset_moduli": [2, 500],
        "equality_units_per_modulus": "first min(8,phi(d))",
        "random_stratum_moduli": [2, 200],
        "random_stratum_sets_per_modulus": 64,
        "random_seed": 2984,
    }
    payload = {
        "profile": profile,
        "modulus_rows": rows,
        "translated_rows": translated_rows,
        "translated_phase_checks": translated_phase_checks,
        "affine_containment_checks": affine_containment_checks,
        "affine_realization_checks": affine_realization_checks,
        "centered_subcomplex_checks": centered_subcomplex_checks,
        "affine_subset_equivalence_checks": affine_subset_equivalence_checks,
        "cyclic_gap_formula_checks": cyclic_gap_formula_checks,
        "affine_pair_checks": affine_pair_checks,
        "affine_flag_subset_checks": affine_flag_subset_checks,
        "composition_orbit_checks": composition_orbit_checks,
        "prime_ratio_checks": prime_ratio_checks,
        "finite_flag_gap_rows": finite_flag_gap_rows,
        "finite_flag_gap_checks": finite_flag_gap_checks,
        "symbolic_m_checks": symbolic_m_checks,
        "equality_modulus_checks": equality_modulus_checks,
        "all_subset_transporter_checks": subset_checks,
        "all_active_subset_checks": active_checks,
        "equality_coset_checks": equality_checks,
        "random_stratum_checks": stratum_checks,
    }
    semantic = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    print("THM2984_NEW_CLAIMS_HOSTILE_AUDIT=PASS")
    print(f"FORMULA_MODULI=2..500 ({len(rows)})")
    print(f"TRANSLATED_CAPACITY_PHASE_CHECKS={translated_phase_checks};MODULI=2..500")
    print(f"AFFINE_CONTAINMENT_CHECKS={affine_containment_checks}")
    print(f"AFFINE_REALIZATION_CHECKS={affine_realization_checks}")
    print(f"CENTERED_SUBCOMPLEX_CHECKS={centered_subcomplex_checks}")
    print(f"AFFINE_SUBSET_EQUIVALENCE_CHECKS={affine_subset_equivalence_checks}")
    print(f"CYCLIC_GAP_FORMULA_CHECKS={cyclic_gap_formula_checks}")
    print(f"AFFINE_PAIR_CHECKS={affine_pair_checks};MODULI=2..100")
    print(f"AFFINE_FLAG_SUBSET_CHECKS={affine_flag_subset_checks};MODULI=2..15")
    print(f"AFFINE_SMALL_NONFLAG_WITNESSES={affine_small_flag_witnesses}")
    print("AFFINE_FIRST_NONFLAG_CONTROL=D9:{0,1,2}")
    print(f"COMPOSITION_ORBIT_CHECKS={composition_orbit_checks}")
    print(f"PRIME_TRIPLE_RATIO_CHECKS={prime_ratio_checks}")
    print(f"FINITE_FLAG_GAP_CHECKS={finite_flag_gap_checks};MODULI=18,21,24,30,36,42")
    print(f"SYMBOLIC_MULTIPLE_SEVEN_CHECKS={symbolic_m_checks};M=4..1000_EXCEPT_6")
    print(f"FLAG_EQUALITY_MODULUS_CHECKS={equality_modulus_checks};D=2..5000")
    print("PAIR_TO_UNIT_FULL_MODULI=2..100")
    print("ALL_SUBSET_AND_ACTIVE_MODULI=2..12")
    print("EQUALITY_COSET_MODULI=2..500;UNITS=FIRST_MIN(8,PHI(D))")
    print("RANDOM_STRATUM_MODULI=2..200;SETS_PER_D=64;SEED=2984")
    print(f"ALL_SUBSET_TRANSPORTER_CHECKS={subset_checks}")
    print(f"ALL_ACTIVE_SUBSET_CHECKS={active_checks}")
    print(f"EQUALITY_COSET_CHECKS={equality_checks}")
    print(f"RANDOM_STRATUM_CHECKS={stratum_checks}")
    print(f"SEMANTIC_SHA256={semantic}")


if __name__ == "__main__":
    main()

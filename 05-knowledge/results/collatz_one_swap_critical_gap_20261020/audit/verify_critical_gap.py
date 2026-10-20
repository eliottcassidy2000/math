#!/usr/bin/env python3
"""Independent finite checks for the critical-shape exponential-content gap.

The uniform quantifiers are proved in AUDIT_CRITICAL_GAP.md, not inferred from
this finite audit. No repository mathematical code is imported. All comparisons
below are exact integer or rational comparisons. This script changes no files.
"""

import argparse
from collections import Counter
from fractions import Fraction
from math import gcd
import json

import sympy as sp


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def rho_monomial(p, exponent, coefficient=1):
    """Normal form in Q[rho]/(rho**p-2), represented by a sparse dictionary."""
    quotient, remainder = divmod(exponent, p)
    return {remainder: coefficient * 2**quotient}


def add_term(polynomial, p, exponent, coefficient):
    for power, value in rho_monomial(p, exponent, coefficient).items():
        polynomial[power] = polynomial.get(power, 0) + value
        if polynomial[power] == 0:
            del polynomial[power]


def direct_fourier(p, a):
    """Expand Q times conjugate Q directly, without the proposed Y formula."""
    q = [(0, 1), (a, -1), (a + 1, 1)]
    result = [{} for _ in range(p)]
    for exponent1, sign1 in q:
        for exponent2, sign2 in q:
            add_term(
                result[(exponent1 - exponent2) % p],
                p,
                exponent1 + exponent2,
                sign1 * sign2,
            )
    return result


def audit_aliases(max_p):
    cases = Counter()
    for p in range(3, max_p + 1):
        for a in range(1, p - 1):
            direct = direct_fourier(p, a)
            proposed = [{} for _ in range(p)]
            for exponent in (0, 2 * a, 2 * a + 2):
                add_term(proposed[0], p, exponent, 1)
            for frequency, exponent, sign in (
                (a + 1, a + 1, 1),
                (a, a, -1),
                (1, 2 * a + 1, -1),
            ):
                add_term(proposed[frequency % p], p, exponent, sign)
                add_term(proposed[-frequency % p], p, exponent, sign)
            require(direct == proposed, ("Fourier expansion", p, a))
            require(
                all(direct[k] == direct[-k % p] for k in range(p)),
                ("real conjugate symmetry", p, a),
            )
            expected = {}
            if p == 3:
                tag = "p3_a1"
                for exponent, coefficient in ((2, 1), (1, -1), (0, -2)):
                    add_term(expected, p, exponent, coefficient)
            elif p == 4 and a == 2:
                tag = "p4_a2_self_conjugate"
                expected = rho_monomial(p, 2, -2)
                require(direct[2] == expected, ("self-conjugate frequency", p, a))
                cases[tag] += 1
                continue
            elif a == 1:
                tag = "a1"
                for exponent in (1, 3):
                    add_term(expected, p, exponent, -1)
            elif a == p - 2:
                tag = "a_pminus2"
                add_term(expected, p, p - 1, 1)
                add_term(expected, p, 2 * p - 3, -1)
            else:
                tag = "generic"
                expected = rho_monomial(p, 2 * a + 1, -1)
            require(direct[1] == expected, ("frequency 1 classification", p, a))
            cases[tag] += 1
    return dict(sorted(cases.items()))


def audit_constants():
    lower_c = Fraction(7497, 2500)
    require(lower_c > Fraction(299, 100), "c rational lower bound")
    require(3**5 < 2**8, "log_2(3)<8/5 and p>=5 alias constants")
    require(3**3 < 2**5, "3^(3/5)<2")
    require(299**6 > 9 * 200**6, "endpoint absorption from p=6")
    require(7 < Fraction(299, 100) ** 2, "minus-interior mean bound")
    require(Fraction(3, 2) ** 3 > 2, "p=3 rho<3/2")
    require(Fraction(3, 2) ** 2 - Fraction(3, 2) < 1, "p=3 alias coefficient")
    require(4**2 * 2 < 6**2, "1+4*sqrt(2)<7")
    require(4 < 2**3, "4^(1/3)<2")
    require(4**3 < 3**4, "4^(3/4)<3")
    critical_small_r = {}
    for p in range(3, 7):
        r = (3**p).bit_length() - p
        require(4 * r <= 3 * p, ("small-p r/p<=3/4", p, r))
        critical_small_r[str(p)] = r
    small = {}
    for p in range(2, 6):
        require(p * 3**p < 5000, ("small-p exponential lower tangent", p))
        require(
            (3**p - 2) * 2500**p < 7497**p,
            ("stronger exact finite small-p rational comparison", p),
        )
        small[str(p)] = p * 3**p
    return {
        "c_lower_bound": str(lower_c),
        "endpoint_absorption_integer_margin": 299**6 - 9 * 200**6,
        "small_p_times_3_to_p": small,
        "small_critical_r": critical_small_r,
        "2_to_8_minus_3_to_5": 2**8 - 3**5,
        "2_to_5_minus_3_cubed": 2**5 - 3**3,
    }


def audit_endpoint_admissibility(max_p):
    """Check the two large endpoints are forbidden at k=1 for every r."""
    coprime_shapes = 0
    for p in range(2, max_p + 1):
        for r in range(1, p):
            if gcd(p, r) != 1:
                continue
            coprime_shapes += 1
            # a=0: the previous residue is r, so the preceding gap is 1.
            a = 0
            previous = (a + r) % p
            preceding_delta = (r + a - previous) // p
            require(preceding_delta == 0, ("large minus endpoint allowed", p, r))
            # a=p-1: the next residue is p-1-r, so the following gap is 1.
            a = p - 1
            following = (a - r) % p
            following_delta = (r + following - a) // p
            require(following_delta == 0, ("large plus endpoint allowed", p, r))
    return {
        "coprime_k1_shapes": coprime_shapes,
        "large_endpoint_forbidden_checks": 2 * coprime_shapes,
    }


def audit_norms(max_p):
    z = sp.Symbol("z")
    counts = Counter()
    for p in range(2, max_p + 1):
        f = z**p - 2
        for a in range(p):
            for direction in (-1, 1):
                q = (
                    1 + z ** (a + 1) - z**a
                    if direction == 1
                    else z ** (p - a) - z + 1
                )
                norm = abs(int(sp.resultant(f, q, z)))
                require(norm > 0, ("zero norm", p, a, direction))
                if 1 <= a <= p - 2:
                    require(
                        norm * 2500**p < 7497**p,
                        ("interior norm versus rational lower bound for c", p, a, direction),
                    )
                    counts["interior_norms"] += 1
                else:
                    expected = {
                        (0, 1): 2,
                        (p - 1, 1): 3**p - 2 ** (p - 1),
                        (0, -1): 3**p - 2,
                        (p - 1, -1): 1,
                    }[(a, direction)]
                    require(norm == expected, ("endpoint norm", p, a, direction))
                    counts["endpoint_norms"] += 1
                critical_q = (3**p).bit_length()
                critical_r = critical_q - p
                if gcd(p, critical_q) == 1:
                    legal = (
                        a < critical_r
                        if direction == 1
                        else a >= p - critical_r
                    )
                    if legal:
                        require(
                            norm**2 < 7**p,
                            ("critical admissible norm versus sqrt(7)", p, a, direction),
                        )
                        counts["critical_admissible_norms_below_sqrt7"] += 1
    return dict(sorted(counts.items()))


def carry(word):
    """Affine composition, without a positional carry sum."""
    numerator = 0
    power_two = 1
    for bit in word:
        numerator = 3 * numerator + power_two if bit else numerator
        power_two *= 2
    return numerator


def audit_critical_words(max_p):
    counts = Counter()
    counts["minus_large_endpoint_occurrences"] = 0
    counts["plus_large_endpoint_occurrences"] = 0
    largest_content = {"g": 0}
    for p in range(2, max_p + 1):
        three_p = 3**p
        q = three_p.bit_length()
        require(2 ** (q - 1) < three_p < 2**q, ("critical q", p, q))
        if gcd(p, q) != 1:
            counts["noncoprime_critical_shapes_omitted"] += 1
            continue
        counts["coprime_critical_shapes"] += 1
        r = q - p
        require(1 <= r < p, ("critical k=1", p, q))
        require(4 * r <= 3 * p, ("critical r/p<=3/4", p, q))
        # Exact r < 3p/5+1, needed by the endpoint proof.
        require(5 * r < 3 * p + 5, ("endpoint r bound", p, q))
        word = [((i + 1) * p) // q - (i * p) // q for i in range(q)]
        positions = {
            (j + (j * r + p - 1) // p - 1) % q: j
            for j in range(p)
        }
        require(set(positions) == {i for i, bit in enumerate(word) if bit}, "indexing")
        for position, index in positions.items():
            a = (-index * r) % p
            plus_legal = word[(position + 1) % q] == 0
            minus_legal = word[(position - 1) % q] == 0
            require(plus_legal == (a < r), ("plus legality iff", p, q, a))
            require(minus_legal == (a >= p - r), ("minus legality iff", p, q, a))
            counts["direction_legality_iff_checks"] += 2
        d = 2**q - three_p
        for i in range(q):
            j = (i + 1) % q
            if word[i] == word[j]:
                continue
            direction = 1 if word[i] else -1
            moved_position = i if direction == 1 else j
            a = (-positions[moved_position] * r) % p
            require(
                (a, direction) not in ((0, -1), (p - 1, 1)),
                ("large endpoint must be forbidden at k=1", p, q, i),
            )
            edited = word.copy()
            edited[i], edited[j] = edited[j], edited[i]
            g = gcd(d, carry(edited))
            counts["actual_cyclic_swap_occurrences"] += 1
            counts["cut_crossing_occurrences"] += int(i == q - 1)
            require(
                g * 2500**p < 7497**p,
                ("critical content versus rational lower bound for c", p, q, i, g),
            )
            require(g**2 < 7**p, ("critical content versus sqrt(7)", p, q, i, g))
            if (a, direction) in ((0, 1), (p - 1, -1)):
                require(g == 1, ("free endpoint", p, q, i, g))
                counts["free_endpoint_occurrences"] += 1
            if g > largest_content["g"]:
                largest_content = {"g": g, "p": p, "q": q, "a": a, "direction": direction}
    return {"counts": dict(sorted(counts.items())), "largest_content": largest_content}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-alias-p", type=int, default=400)
    parser.add_argument("--max-norm-p", type=int, default=20)
    parser.add_argument("--max-word-p", type=int, default=200)
    args = parser.parse_args()
    result = {
        "status": "PASS",
        "scope": "Finite exact corroboration, not a proof of uniform quantifiers.",
        "parameters": vars(args),
        "constants": audit_constants(),
        "fourier_alias_classes": audit_aliases(args.max_alias_p),
        "k1_endpoint_admissibility": audit_endpoint_admissibility(args.max_alias_p),
        "exact_norms": audit_norms(args.max_norm_p),
        "critical_words": audit_critical_words(args.max_word_p),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

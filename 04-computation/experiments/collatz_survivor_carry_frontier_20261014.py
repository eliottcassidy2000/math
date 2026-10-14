#!/usr/bin/env python3
"""Bounded carry-frontier exploration, exact all-height source inequalities."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import comb, log2
import hashlib
import json

checks = 0


def require(ok, message):
    global checks
    checks += 1
    if not ok:
        raise ArithmeticError(message)


def carrier(word):
    A, B, P = 0, 0, 1
    for a in word:
        B, A, P = 3 * B + (1 << A), A + a, 3 * P
    return A, B, P


def value(n):
    x = 3 * n + 1
    a = (x & -x).bit_length() - 1
    return a, x // (1 << a)


def encode(word):
    A, B, P = carrier(word)
    return B * pow(1 << A, -1, P) % P


minimum_cost, polynomials = {0: 0}, {0: Counter({0: 1})}
rows = []
for L in range(1, 7):
    cap_bound = 4 * Fraction(3, 2) ** (L - 1) - Fraction(2, 3)
    require(cap_bound < 32, "all-height per-letter valuation cap four")
    cells = defaultdict(list)
    finite_only = []
    excluded = 0
    for word in product(range(1, 5), repeat=L):
        A, B, P, bound = 0, 0, 1, None
        for a in word:
            B, A, P = 3 * B + (1 << A), A + a, 3 * P
            if P < 1 << A:
                threshold = Fraction(B, (1 << A) - P)
                bound = threshold if bound is None else min(bound, threshold)
        modulus = 1 << (A + 1)
        root = ((1 << A) - B) * pow(P, -1, modulus) % modulus
        first = root if root >= 3 else root + modulus
        # Independently authenticate this exact word at its two first sources >1.
        for n in (first, first + modulus):
            x = n
            actual_survival = True
            for a in word:
                got, x = value(x)
                require(got == a, "native cylinder actual valuation")
                actual_survival &= x >= n
            require(actual_survival == (bound is None or n <= bound), "source inequalities exact")
        if bound is not None and first > bound:
            excluded += 1
            continue
        residue = B * pow(1 << A, -1, P) % P
        item = {"word": word, "A": A, "B": B, "least_source_above_one": first,
                "native_modulus": modulus, "relative_odd_mass": str(Fraction(1, 1 << A)),
                "upper_source_bound": None if bound is None else str(bound)}
        cells[residue].append(item)
        if bound is not None:
            finite_only.append(item)
    require(not finite_only, "no finite-only exceptions through length six")

    # Independent minimum-cost and weighted-polynomial projected compilers.
    budget, new_min, new_poly = (3 ** L).bit_length() - 1, {}, defaultdict(Counter)
    for old_r, A in minimum_cost.items():
        for a in range(1, budget - A + 1):
            r = (3 * old_r + 1) * pow(1 << a, -1, 3 ** L) % (3 ** L)
            new_min[r] = min(new_min.get(r, budget + 1), A + a)
    for old_r, coefficients in polynomials.items():
        for A, multiplicity in coefficients.items():
            for a in range(1, budget - A + 1):
                r = (3 * old_r + 1) * pow(1 << a, -1, 3 ** L) % (3 ** L)
                new_poly[r][A + a] += multiplicity
    require(new_min == {r: min(w["A"] for w in ws) for r, ws in cells.items()}, "minimum-cost compiler")
    require(dict(new_poly) == {r: Counter(w["A"] for w in ws) for r, ws in cells.items()},
            "weighted-polynomial compiler")
    minimum_cost, polynomials = new_min, dict(new_poly)

    # Binary parity-decoding control and the all-length central-binomial lower family.
    binary_words = list(product((1, 2), repeat=L))
    binary_codes = [encode(w) for w in binary_words]
    require(len(set(binary_codes)) == 2 ** L, "binary alphabet injectivity")
    ballot_words = [w for w in binary_words if all(sum(w[:i]) * 2 <= 3 * i for i in range(1, L + 1))]
    require(len(ballot_words) == comb(L, L // 2), "ballot central-binomial count")
    require(all(encode(w) in cells for w in ballot_words), "ballot family is slope admissible")

    mass = sum((Fraction(1, 1 << w["A"]) for ws in cells.values() for w in ws), Fraction())
    naive_mass = sum((Fraction(1, 1 << min(w["A"] for w in ws)) for ws in cells.values()), Fraction())
    collisions = []
    for r, ws in sorted(cells.items()):
        if len(ws) > 1:
            require(len(ws) == 2, "small collision fibres are pairs")
            high, low = sorted(ws, key=lambda w: w["A"], reverse=True)
            wh, wl = tuple(high["word"]), tuple(low["word"])
            rewritten = False
            for m in range(L - 1):
                if wh[:m + 1] == (1,) * (m + 1) and wh[m + 1] >= 3:
                    candidate = (1,) * m + (2, wh[m + 1] - 2) + wh[m + 2:]
                    rewritten |= candidate == wl
            require(rewritten, "all observed collisions explained by suffix-stable rewrite")
            collisions.append({"carry": r, "words": ws,
                               "correct_fibre_mass": str(sum((Fraction(1, 1 << w["A"]) for w in ws), Fraction()))})
    rows.append({"length": L, "valuation_universe": [1, 2, 3, 4],
                 "all_height_letter_bound": str(cap_bound), "examined_words": 4 ** L,
                 "excluded_words": excluded, "admissible_words": sum(map(len, cells.values())),
                 "infinite_source_families": sum(map(len, cells.values())), "finite_only_families": 0,
                 "carry_classes": len(cells), "unit_classes": 2 * 3 ** (L - 1),
                 "relative_odd_survival_density": str(mass),
                 "wrong_minimum_representative_density": str(naive_mass),
                 "mass_lost_by_wrong_replacement": str(mass - naive_mass),
                 "binary_ballot_lower_bound": comb(L, L // 2), "collision_fibres": collisions,
                 "minimum_cost_frontier": {str(r): A for r, A in sorted(new_min.items())}})

rewrite_tests = 0
for m in range(5):
    for b in range(3, 6):
        left, right = (1,) * (m + 1) + (b,), (1,) * m + (2, b - 2)
        for suffix in ((), (1,), (2,), (1, 2)):
            if len(left + suffix) > 6:
                continue
            A1, B1, P1 = carrier(left + suffix)
            A2, B2, P2 = carrier(right + suffix)
            require(A1 == A2 + 1 and P1 == P2, "rewrite length and cost")
            require(2 * B2 - B1 == P1, "rewrite exact affine difference")
            require(encode(left + suffix) == encode(right + suffix), "suffix-stable carry collision")
            rewrite_tests += 1

alpha = log2(3)
beta = alpha ** alpha / (alpha - 1) ** (alpha - 1)
result = {"status": "FINITE-EXACT all-height classification through six; general identities proved in note",
          "checks": checks, "rows": rows, "rewrite_examples_checked": rewrite_tests,
          "entropy_composition_rate_approx": beta,
          "open_hypothesis": "limsup(carry_class_count(L)^(1/L)) is strictly below the slope-word rate beta",
          "limits": ["No proof of a strict entropy gap", "Carry compression is not sourcewise descent",
                     "Finite-only no-descent source families are excluded only through length six"]}
payload = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
result["payload_sha256"] = hashlib.sha256(payload).hexdigest()
print(json.dumps(result, sort_keys=True, indent=2))

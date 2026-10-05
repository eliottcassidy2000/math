#!/usr/bin/env python3
"""Exact finite audit; universal weak coverage uses Applegate--Lagarias, not this census.

Run from the repository root. No third-party packages. Assertions are deliberately
not used as the verification mechanism, so -O must produce identical JSON.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
import json


CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def v2(n):
    n = abs(n)
    if not n:
        raise ValueError("valuation of zero")
    return (n & -n).bit_length() - 1


def step(n, sign=1):
    return (3*n + sign) // (1 << v2(3*n + sign))


def root_path(n, sign=1, limit=10000):
    path, seen = [], set()
    while n != 1:
        if n in seen or len(path) >= limit:
            return None
        seen.add(n)
        path.append(n)
        n = step(n, sign)
    return path


def tidy(vector):
    return {x: a for x, a in sorted(vector.items()) if a}


def boundary(edges, sign=1):
    result = Counter()
    for u, multiplicity in edges.items():
        check(u > 0 and u % 2 and multiplicity >= 0, "edge type")
        result[u] += multiplicity
        result[step(u, sign)] -= multiplicity
    return tidy(result)


def value(edges, sign=1):
    result = Fraction(1)
    for u, multiplicity in edges.items():
        result *= Fraction(u, step(u, sign))**multiplicity
    return result


def factors(n):
    result = Counter()
    p = 2
    while p*p <= n:
        while n % p == 0:
            result[p] += 1
            n //= p
        p += 1
    if n > 1:
        result[n] += 1
    return result


def prime_image(vector):
    result = Counter()
    for u, multiplicity in vector.items():
        for p, exponent in factors(u).items():
            result[p] += multiplicity * exponent
    return tidy(result)


def defect(n, edges):
    result = Counter(boundary(edges))
    result[n] -= 1
    result[1] += 1
    return tidy(result)


def factor_relation(a, b):
    result = Counter({a*b: 1})
    result[1] += 1
    result[a] -= 1
    result[b] -= 1
    return tidy(result)


def normal_form(vector):
    """Reduce [u] to prime labels minus (Omega(u)-1)[1]."""
    result = Counter()
    for u, multiplicity in vector.items():
        fs = factors(u)
        for p, exponent in fs.items():
            result[p] += multiplicity * exponent
        result[1] += multiplicity * (1-sum(fs.values()))
    return tidy(result)


def grounded_extract(n, edges, grounded, sign=1):
    """Check the sufficient source/deficit rule; then extract its literal path."""
    b = boundary(edges, sign)
    if n in grounded:
        return []
    if b.get(n, 0) <= 0 or any(a < 0 and u not in grounded for u, a in b.items()):
        return None
    route, seen = [], set()
    u = n
    while u not in grounded:
        check(u not in seen and edges.get(u, 0) > 0, "grounded-sink extraction")
        seen.add(u)
        route.append(u)
        u = step(u, sign)
    return route


def factor_certificate(n):
    """Finite test constructor: multiply independently checked prime root paths."""
    result = Counter()
    for p, exponent in factors(n).items():
        path = root_path(p)
        check(path is not None, "finite prime certificate")
        for u in path:
            result[u] += exponent
    return result


def prefix_data(n, m):
    p, q, b = 1, 1, 0
    word = []
    d = n
    for _ in range(m):
        a = v2(3*d+1)
        word.append(a)
        b = 3*b + q
        p *= 3
        q <<= a
        d = step(d)
    check(p*n+b == q*d, "prefix carry")
    return p, q, b, d, word


def log4_mod3(target, level):
    """Unique e mod 3^(level-1), for target=1 mod 3."""
    e, period = 0, 1
    check(target % 3 == 1 and level >= 1, "principal unit")
    for ell in range(2, level+1):
        modulus = 3**ell
        choices = [e+j*period for j in range(3)
                   if pow(4, e+j*period, modulus) == target % modulus]
        check(len(choices) == 1, "ternary logarithm lift")
        e = choices[0]
        period *= 3
    return e, period


def rooted_shadow(n, m, h):
    p, q, b, d, word = prefix_data(n, m)
    next_a = v2(3*d+1)
    e, period = log4_mod3(3*d+1, m+h+1)
    lower = max(2, next_a//2+1)
    if e < lower:
        e += ((lower-e+period-1)//period)*period
    while q*(4**e-1) <= 3*b:
        e += period
    numerator = q*(4**e-1)-3*b
    check(numerator % (3*p) == 0, "shadow integrality")
    r = numerator // (3*p)
    check(r > 0 and r % 2, "positive odd shadow")
    check((r-n) % 3**h == 0, "ternary source agreement")
    check(v2(r-n) == (q.bit_length()-1)+next_a, "exact first binary difference")
    # Independent literal execution, including exact valuations and first hit.
    x = r
    for a in word + [2*e]:
        check(x > 1, "no root padding")
        check(v2(3*x+1) == a, "literal shadow valuation")
        x = step(x)
    check(x == 1, "shadow terminal")
    check(p*r+b == q*((4**e-1)//3), "shadow affine endpoint")
    return {"n": n, "m": m, "h": h, "e": e,
            "source_bits": r.bit_length(), "word": word + [2*e],
            "source": str(r) if r.bit_length() < 256 else None}


def main():
    report = {}
    # Direct certificates and factorization defects throughout a declared census.
    accepted_factor = []
    anchored_3 = 0
    for n in range(1, 1002, 2):
        c = factor_certificate(n)
        check(value(c) == n, "finite weak certificate")
        d = defect(n, c)
        check(sum(d.values()) == 0 and not prime_image(d), "defect kernel")
        check(not normal_form(d), "factorization normal form")
        literal = root_path(n)
        check(literal is not None, "finite target control")
        check(not defect(n, Counter(literal)), "literal zero defect")
        check(grounded_extract(n, Counter(literal), {1}) == literal, "path extraction")
        extracted = grounded_extract(n, c, {1})
        if extracted is not None:
            accepted_factor.append(n)
            check(extracted == literal, "factor extraction correctness")
        # Source anchoring preserves multiplicative value; v3 pins source n.
        if n % 3 == 0:
            a = factor_certificate(step(n))
            a[n] += 1
            check(value(a) == n, "anchored value")
            check(boundary(a).get(n) == 1, "distinguished source port")
            check(all(u == n or u % 3 for u in a), "all other sources 3-units")
            anchored_3 += 1
    report["finite_factor_census"] = {"positive_odds_through": 1001,
        "count": 501, "grounded_factor_certificates": len(accepted_factor),
        "source_anchored_multiples_of_3": anchored_3,
        "composite_grounded_factor_certificates":
            [n for n in accepted_factor if n > 1 and sum(factors(n).values()) > 1]}
    examples = {}
    for n in (9, 27, 65, 81, 243):
        c = factor_certificate(n)
        examples[str(n)] = {"edges": dict(c), "boundary": boundary(c),
                           "defect": defect(n, c),
                           "grounded_from_1": grounded_extract(n, c, {1}) is not None}
    report["factor_examples"] = examples
    for k in range(1, 41):
        n = 3**k
        c = Counter({3: k, 5: k})
        expected = Counter({3: k})
        expected[n] -= 1
        expected[1] -= k-1
        check(value(c) == n and defect(n, c) == tidy(expected), "power-3 fusion")
    # The graph lemma also works without a multiplicative certificate assumption.
    extra = Counter(root_path(3)) + Counter(root_path(7))
    check(value(extra) == 21 and defect(3, extra), "extra flow is not weak at target3")
    check(grounded_extract(3, extra, {1}) == [3, 5], "extra sources are harmless")
    # Applegate--Lagarias Table 3's reciprocal certificate for 5, accelerated.
    inverse5 = Counter({7: 2, 11: 1, 17: 1, 55: 1, 65: 1, 83: 1})
    check(value(inverse5) == Fraction(1, 5), "published wild-5 certificate")
    neutral = inverse5 + Counter({5: 1})
    check(value(neutral) == 1 and boundary(neutral), "nontrivial neutral packet")
    check(all(u % 3 for u in neutral), "neutral packet contains no 3-source")
    adorned7 = Counter(root_path(7)) + neutral
    grounded = {u for u, coefficient in boundary(adorned7).items() if coefficient < 0}
    check(all(root_path(u) is not None for u in grounded), "checked deficit roots")
    check(value(adorned7) == 7 and defect(7, adorned7), "weak with nonzero defect")
    route = grounded_extract(7, adorned7, grounded)
    check(route is not None, "grounded nonzero-defect receipt")
    nine = Counter({3: 2, 5: 2}) + neutral
    check(boundary(nine).get(3) == 2 and boundary(nine).get(9, 0) == 0,
          "neutral addition cannot repair distinguished 3-source")
    report["neutral_packet"] = {"edges": dict(neutral), "boundary": boundary(neutral),
        "grounded7_deficit_set": sorted(grounded), "grounded7_route": route}
    # A positive-sheet 3x-1 cycle is a hostile to weak => rooted.
    minus = Counter({25: 1, 37: 1, 11: 1})
    check(value(minus, -1) == 5, "3x-1 weak cycle certificate")
    check(root_path(5, -1) is None, "3x-1 cycle hostile")
    check(grounded_extract(5, minus, {1}, -1) is None, "cycle hostile rejected")
    cycle_added = minus + Counter({5: 1, 7: 1})
    check(value(cycle_added, -1) == 5, "cycle addition is neutral")
    check(boundary(cycle_added, -1).get(5, 0) == 0, "cycle is not source surplus")
    check(grounded_extract(5, cycle_added, {1}, -1) is None, "anchor alone insufficient")
    report["minus_cycle_hostile"] = {"weak_value": 5, "edges": dict(minus),
                                     "boundary": boundary(minus, -1),
                                     "cycle": [5, 7, 5]}
    # Independent finite functional-graph audit of the grounded deficit criterion.
    graph_cases = 0
    for mapping in product(range(4), repeat=4):
        for stock in product(range(2), repeat=4):
            b = [0]*4
            for u in range(4):
                b[u] += stock[u]
                b[mapping[u]] -= stock[u]
            if any(b[u] < 0 for u in (1, 2, 3)):
                continue
            for source in (1, 2, 3):
                if b[source] <= 0:
                    continue
                u, seen = source, set()
                while u != 0 and u not in seen and stock[u]:
                    seen.add(u)
                    u = mapping[u]
                check(u == 0, "independent four-vertex grounded flow")
                graph_cases += 1
    report["four_vertex_grounded_flow_cases"] = graph_cases
    # Exact fusion transport and obstruction to one-step multiplicativity.
    pairs, below, above = 0, 0, 0
    for a in range(3, 202, 2):
        for b in range(3, 202, 2):
            ratio = Fraction((3*a+1)*(3*b+1), 3*a*b+1)
            check(3 < ratio < 4, "one-step fusion ratio")
            k = v2(3*a+1)+v2(3*b+1)-v2(3*a*b+1)
            check((step(a*b) < step(a)*step(b)) == (k <= 1), "fusion orientation")
            check(step(a*b) != step(a)*step(b), "nonmultiplicativity")
            left = Counter(factor_relation(a, b))
            right = Counter(boundary(Counter({a*b: 1})))
            for x, coef in boundary(Counter([a, b])).items():
                right[x] -= coef
            for x, coef in factor_relation(step(a), step(b)).items():
                right[x] += coef
            right[step(a*b)] += 1
            right[step(a)*step(b)] -= 1
            check(tidy(left) == tidy(right), "fusion transport")
            pairs += 1
            below += k <= 1
            above += k >= 2
    report["fusion_pairs"] = {"pairs": pairs, "product_image_smaller": below,
                              "product_image_larger": above}
    cocycles = 0
    for a, b, c in product(range(1, 20, 2), repeat=3):
        left, right = Counter(), Counter()
        for relation in (factor_relation(a, b), factor_relation(a*b, c)):
            for u, coefficient in relation.items():
                left[u] += coefficient
        for relation in (factor_relation(b, c), factor_relation(a, b*c)):
            for u, coefficient in relation.items():
                right[u] += coefficient
        check(tidy(left) == tidy(right), "associative fusion cocycle")
        cocycles += 1
    report["fusion_cocycle_triples"] = cocycles
    # Re-express the inherited 4u+1 inverse ray as a checked infinite fusion family.
    fusion_rays = []
    for a in range(3, 52, 2):
        k = 1
        while pow(4, k, 3*a) != 1:
            k += 1
        for t in (1, 2, 3):
            exponent = k*t
            b = 4**exponent + (4**exponent-1)//(3*a)
            check((4**exponent-1) % (3*a) == 0, "fusion ray integrality")
            check(a*b == 4**exponent*a + (4**exponent-1)//3, "fusion ray source")
            check(step(a*b) == step(a), "fusion ray common future")
            check(v2(3*a*b+1) == 2*exponent+v2(3*a+1), "fusion ray price")
            first = root_path(a)
            check(first is not None, "finite fusion seed")
            route = [a*b] + first[1:]
            check(not defect(a*b, Counter(route)), "fusion ray zero defect")
            fusion_rays.append({"a": a, "k": exponent, "b": str(b),
                                "product_bits": (a*b).bit_length()})
    report["fusion_rays"] = {"count": len(fusion_rays),
        "seeds": "3,5,...,51", "order_multiples": [1, 2, 3],
        "examples": [x for x in fusion_rays if x["a"] == 3]}
    # Rooted shadows even of negative-cycle trajectories.
    sources = list(range(1, 100, 2)) + [-1, -5, -7, -17, -27]
    shadows = [rooted_shadow(n, m, h) for n in sources
               for m in range(6) for h in range(3)]
    report["rooted_shadows"] = {"count": len(shadows),
        "positive_sources": "1,3,...,99", "negative_sources": sources[-5:],
        "prefix_lengths": "0..5", "ternary_depths": "0..2",
        "max_source_bits": max(s["source_bits"] for s in shadows),
        "examples": [s for s in shadows if s["n"] in (-1, -5, -17, 27)
                     and s["m"] == 2 and s["h"] == 1]}
    report["checks"] = CHECKS
    report["status"] = "PASS: finite exact audit only; Collatz remains OPEN"
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

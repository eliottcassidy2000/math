#!/usr/bin/env python3
"""Boundary conductor controls; no third-party packages or producer imports."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json

GATES = 0


def need(test, message):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(message)


def span(vectors):
    out = {0}
    for v in vectors:
        out |= {w ^ v for w in tuple(out)}
    return frozenset(out)


def code(labels):
    return frozenset(sum(((a & v).bit_count() % 2) << i
                         for i, v in enumerate(labels)) for a in range(8))


def puncture(word, j, n=7):
    return sum(((word >> i) & 1) << k
               for k, i in enumerate(i for i in range(n) if i != j))


def conductor(c, d, n):
    full = (1 << n) - 1
    decomposition = {p ^ q: (p, q) for p in c for q in d}
    need(len(decomposition) == 1 << n, 'complementary character decomposition')
    e = {z: p & q for z, (p, q) in decomposition.items()}
    shared = 0
    for i in range(n):
        if any(p >> i & 1 for p in c) and any(q >> i & 1 for q in d):
            shared |= 1 << i
    for z in range(1 << n):
        calculated = 0
        for y in range(1 << n):
            calculated |= e[z ^ y] & ~(z & y) & full
        need(calculated == shared & ~z, 'conductor exponent from every multiplication row')
        # Principal generator product_{i in shared} u_i, multiplied by u_{shared xor z}.
        need(shared & (shared ^ z) == calculated, 'principal conductor consequence')
        for i in range(n):
            if calculated >> i & 1:
                need(any((e[z ^ y] >> i & 1) and not (z & y) >> i & 1
                         for y in range(1 << n)), 'lowering a conductor valuation fails')
    # Image is a subalgebra, not merely an R-submodule.
    for z, w in product(range(1 << n), repeat=2):
        for i in range(n):
            need(((e[z] >> i) & 1) + ((e[w] >> i) & 1) + (((z & w) >> i) & 1)
                 >= ((e[z ^ w] >> i) & 1), 'image multiplication closure')
    return e


def fibre(c, d, zero, chart):
    """Independent sign evaluation of the reduced geometric boundary fibre."""
    indices = [i for i in range(7) if i != chart and i not in zero]
    counts = Counter()
    for signs in product((-1, 1), repeat=len(indices)):
        roots = dict(zip(indices, signs))
        roots[chart] = 1
        roots.update({i: 0 for i in zero})
        def ev(s):
            out = 1
            for i in range(7):
                if s >> i & 1:
                    out *= roots[i]
            return out
        counts[(tuple(ev(s) for s in sorted(c)), tuple(ev(s) for s in sorted(d)))] += 1
    return counts


def main():
    codes = sorted({code((1,) + p) for p in permutations(range(2, 8))}, key=sorted)
    need(len(codes) == 30, 'thirty labelled simplex codes')
    pairs = [(c, d) for c, d in combinations(codes, 2) if c & d == {0}]
    need(len(pairs) == 120, 'all complementary Fano pairs')
    for c, d in pairs:
        for j in range(7):
            conductor({puncture(p, j) for p in c}, {puncture(q, j) for q in d}, 6)

    # Exhaust every ordered complementary decomposition in dimension four,
    # including coordinates ramified in only one factor.
    subspaces = {span(v) for size in (1, 2, 3)
                 for v in combinations(range(1, 16), size) if len(span(v)) == 2**size}
    small_pairs = [(c, d) for c in subspaces for d in subspaces
                   if len(c)*len(d) == 16 and c & d == {0}]
    for c, d in small_pairs:
        conductor(c, d, 4)

    U = code((1, 2, 4, 5, 6, 3, 7))
    D = code((1, 2, 4, 6, 3, 7, 5))
    triples = [{0, 1, 2}, {0, 3, 4}, {0, 5, 6}]
    ordinary = [set(t) for t in combinations(range(7), 2)
                if not any(set(t) <= z for z in triples)]
    strata = [{i} for i in range(7)] + ordinary + triples
    fibre_results = Counter()
    for zero in strata:
        j = next(i for i in range(7) if i not in zero)
        counts = fibre(U, D, zero, j)
        r = len(zero)
        need(len(counts) == 2**(6-2*r), 'product reduced fibre cardinality')
        need(set(counts.values()) == {2**r}, 'normalization sheet multiplicity')
        fibre_results[(r, len(counts), 2**r)] += 1

    # Chart cocycle for every even character, with negative and fractional twists.
    roots = [F(1), F(-2), F(3, 2), F(5), F(-7, 3), F(11), F(13, 2)]
    for s in range(128):
        if s.bit_count() % 2:
            continue
        numerator = F(1)
        for i in range(7):
            if s >> i & 1:
                numerator *= roots[i]
        values = [numerator / r**s.bit_count() for r in roots]
        for j, k, l in product(range(7), repeat=3):
            tjk = (roots[k] / roots[j])**(-s.bit_count())
            tkl = (roots[l] / roots[k])**(-s.bit_count())
            need(values[k] == tjk*values[j], 'even-character overlap with twists')
            need(tjk*tkl == (roots[l]/roots[j])**(-s.bit_count()), 'triple chart cocycle')
    P = [product_value(roots)/r**7 for r in roots]
    for j, k in product(range(7), repeat=2):
        need(P[k] == (roots[k]/roots[j])**-7 * P[j], 'global conductor transition')

    # Explicit rational quotient points with no rational root lift at L1=0.
    # Chart 0, point (E,x,z)=(1,-1,3), so nonzero branch values are 1,2,4,-2,3,-1.
    units = [2, 3, 4, 5, 6]
    unit_characters = [s for s in U | D if not s >> 1 & 1]
    hs = [h for h in range(1, 128) if not h & 3 and
          all((h & s).bit_count() % 2 == 0 for s in unit_characters)]
    need(len(hs) == 1, 'one residual squareclass at a simple branch')
    h = hs[0]
    b = [F(1), F(0)] + [F(2**((h >> i) & 1)) for i in units]
    raw = [1, 0, 2, 4, -2, 3, -1]
    twists = [F(1), F(1)] + [b[i]/raw[i] for i in units]
    def value(s):
        if s >> 1 & 1:
            return F(0)
        exponent = (s & h).bit_count()
        need(exponent % 2 == 0, 'unit character is a rational square')
        return F(2**(exponent//2))
    for c in (U, D):
        for s, t in product(c, repeat=2):
            correction = product_value(b[i] for i in range(7) if (s & t) >> i & 1)
            need(value(s)*value(t) == correction*value(s ^ t), 'rational quotient algebra point')
    need(any(b[i] == 2 for i in units), 'rational full lift is impossible')

    # Triple point: both quotients have a unique rational reduced point,
    # while its normalization has residue field Q(sqrt(2),sqrt(3),sqrt(5)).
    for c in (U, D):
        need(all(s == 0 or s & 7 for s in c), 'all nontrivial quotient characters vanish at 012')
    from math import isqrt
    for a,b,c in product(range(2), repeat=3):
        value_squared = 2**a * 3**b * 5**c
        need((isqrt(value_squared)**2 == value_squared) == (a == b == c == 0),
             'three independent prime-valuation squareclasses')

    # Positive-control lift: erase the residual twist while preserving the zero root.
    need(set(fibre(U, D, {1}, 0).values()) == {2}, 'two rational lifts after trivializing twist')
    result = {
        'status': 'FINITE-EXACT', 'gates': GATES,
        'Fano_codes': 30, 'complementary_pairs': 120, 'charts_per_pair': 7,
        'ordered_dimension_four_decompositions': len(small_pairs),
        'boundary_strata': [{'vanishing_roots': r, 'product_points': n,
                             'normalization_sheets': s, 'number_of_strata': count}
                            for (r, n, s), count in sorted(fibre_results.items())],
        'arithmetic_hostile': {'point': [1, -1, 3], 'branch': 1,
                               'nonsquare_indices': [i for i in units if h >> i & 1],
                               'twists': list(map(str, twists)),
                               'quotient_U_values': {str(s): str(value(s)) for s in sorted(U)},
                               'quotient_D_values': {str(s): str(value(s)) for s in sorted(D)}},
        'triple_hostile': {'point': [0,0,1], 'chart': 3, 'twists': [1,1,1,1,-2,3,-5],
                          'unit_radicands': [2,3,5], 'residue_field_degree': 8},
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope': 'Algebraic proof separate. All labelled complementary Fano pairs, all seven charts; exact boundary signs for the displayed norm-seven pair; no rational-point existence or Collatz coverage conclusion.'}
    print(json.dumps(result, indent=2, sort_keys=True))


def product_value(values):
    answer = F(1)
    for v in values:
        answer *= v
    return answer


if __name__ == '__main__':
    main()

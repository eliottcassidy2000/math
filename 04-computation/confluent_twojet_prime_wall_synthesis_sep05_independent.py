#!/usr/bin/env python3
"""Independent SymPy Smith and symbolic-tariff audit; no primary imports."""

from hashlib import sha256
from itertools import combinations
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def vp(n, p):
    n = abs(int(n))
    require(n != 0, "nonzero valuation")
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def smith(nodes, p):
    n = 2*len(nodes)
    a = s.Matrix([row for x in nodes for row in (
        [x**q for q in range(n)],
        [q*x**(q-1) if q else 0 for q in range(n)])])
    diagonal = smith_normal_form(a, domain=s.ZZ)
    return tuple(vp(diagonal[i, i], p) for i in range(n))


def main():
    walls = []
    for p in (2, 3, 5, 7):
        for e in (1, 2, 3):
            nodes = tuple(7-2*p+p**e*(i+p*(i*i+2)) for i in range(p))
            actual = smith(nodes, p)
            expected = ((0, 0)+tuple(e*j for j in range(1, p-1))
                        +((p-1)*e+1,)+tuple(e*j for j in range(p+1, 2*p-1))
                        +((2*p-1)*e-1,))
            require(actual == expected, ("wall", p, e, actual, expected))
            walls.append((p, e, actual))
    triples = []
    for e in range(9):
        actual = smith((0, 2**e, 2**(e+1)), 2)
        delta = (0, min(e+1, 2*e), min(3*e+2, 4*e), 7*e+2, 12*e+4)
        expected = (0, 0)+tuple(y-x for x, y in zip(delta, delta[1:]))
        require(actual == expected, ("dyadic triple", e, actual, expected))
        triples.append((e, actual))
    rows = ((1, 0), (1, 1), (2, 0), (2, 1))
    table = {}
    minor_count = 0
    for h in range(1, 5):
        costs = {}
        for rr in combinations(range(4), h):
            for cc in combinations(range(2, 6), h):
                a = s.Matrix([[q**r*u**(q-r) for q in cc]
                              for u, r in (rows[i] for i in rr)])
                determinant = a.det(method="domain-ge")
                minor_count += 1
                if determinant:
                    w = sum(cc)-sum(rows[i][1] for i in rr)
                    costs[w] = min(costs.get(w, 10**9), vp(determinant, 2))
        table[h] = costs
    expected_table = {
        1: {1: 1, 2: 0, 3: 0, 4: 0, 5: 0},
        2: {3: 2, 4: 0, 5: 1, 6: 0, 7: 1, 8: 0, 9: 4},
        3: {7: 2, 8: 2, 9: 2, 10: 2, 11: 3},
        4: {12: 4},
    }
    require(table == expected_table, "complete symbolic minor tariff table")
    semantic = json.dumps((walls, triples, table), sort_keys=True, separators=(",", ":"))
    print("CONFLUENT_TWOJET_PRIME_WALL_INDEPENDENT_SYNTHESIS_SEP05")
    print(f"sympy_version={s.__version__}")
    print(f"independent_wall_Smith_cases={len(walls)}")
    print(f"independent_dyadic_triple_Smith_cases={len(triples)}")
    print(f"complete_symbolic_residual_minors={minor_count}")
    print(f"dyadic_triple_tariff_table={json.dumps(table, sort_keys=True)}")
    print(f"semantic_sha256={sha256(semantic.encode()).hexdigest()}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()

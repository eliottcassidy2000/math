#!/usr/bin/env python3
"""Independent exhaustive controls for unit-flow Hodge and chart information.

Universe: every labelled tournament through order five, and every fixed-path
chart through order five. All arithmetic proof gates are integral/rational.
Canonicalization here uses all permutations, independent of the augmentation
producer's refinement/deletion algorithm. No asserts disappear under -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def run():
    result = []
    expected = [1, 1, 1, 2, 4, 12]
    collision = None
    for n in range(6):
        pairs = list(combinations(range(n), 2))
        perms = list(permutations(range(n)))
        chart_fibres = Counter()
        classes = {}
        score_h = defaultdict(list)
        for code in range(1 << len(pairs)):
            a = [[0] * n for _ in range(n)]
            for k, (i, j) in enumerate(pairs):
                a[i][j] = 1 if (code >> k) & 1 else -1
                a[j][i] = -a[i][j]
            scores = [sum(v == 1 for v in row) for row in a]
            d = [sum(row) for row in a]
            x = sum(v * v for v in d)
            c3 = sum(abs(a[i][j] + a[j][k] + a[k][i]) == 3
                     for i, j, k in combinations(range(n), 3))
            need(24 * c3 == n * (n * n - 1) - 3 * x, "triangle-score identity")
            if n:
                residual = [[n * a[i][j] - (d[i] - d[j])
                             for j in range(n)] for i in range(n)]
                need(all(sum(row) == 0 for row in residual), "normal equations")
                need(sum(residual[i][j] * (d[i] - d[j]) for i, j in pairs) == 0,
                     "orthogonality")
                grad = Fraction(sum((d[i] - d[j]) ** 2 for i, j in pairs), n * n)
                cyc = Fraction(sum(residual[i][j] ** 2 for i, j in pairs), n * n)
                need(grad == Fraction(x, n), "correct gradient energy")
                need(grad + cyc == len(pairs), "Pythagoras")
                need(cyc == Fraction((n - 1) * (n - 2), 6) + Fraction(8 * c3, n),
                     "transitive baseline and cyclic excess")
                curl2 = sum((a[i][j] + a[j][k] + a[k][i]) ** 2
                            for i, j, k in combinations(range(n), 3))
                need(curl2 == n * cyc, "complete-complex curl norm")
            codes = [sum((a[p[i]][p[j]] == 1) << k for k, (i, j) in enumerate(pairs))
                     for p in perms]
            canonical = min(codes)
            aut = codes.count(canonical)
            hamilton = sum(all(a[p[i]][p[i + 1]] == 1 for i in range(n - 1))
                           for p in perms)
            need(hamilton % aut == 0, "free action on paths")
            h = hamilton // aut
            data = (hamilton, aut, h, tuple(sorted(scores)), c3)
            if canonical in classes:
                need(classes[canonical] == data, "class invariant transport")
            classes[canonical] = data
            if all(a[i][i + 1] == 1 for i in range(n - 1)):
                chart_fibres[canonical] += 1
        need(len(classes) == expected[n], "independent small-order class census")
        m = comb(n - 1, 2) if n >= 2 else 0
        need(sum(chart_fibres.values()) == 1 << m, "chart universe")
        need(chart_fibres == Counter({c: data[2] for c, data in classes.items()}),
             "chart fibre equals H/Aut")
        need(sum(Fraction(count, classes[c][2]) for c, count in chart_fibres.items())
             == len(classes), "reciprocal chart mass")
        need(sum(factorial(n) // data[1] for data in classes.values()) == 1 << len(pairs),
             "labelled class masses")
        for c, data in classes.items():
            score_h[(data[3], data[0])].append(c)
        if collision is None:
            bad = [(key, cs) for key, cs in score_h.items() if len(cs) > 1]
            if bad:
                key, cs = min(bad)
                collision = {"n": n, "scores": key[0], "H": key[1], "codes": cs,
                             "automorphisms": [classes[c][1] for c in cs],
                             "cyclic_triangles": [classes[c][4] for c in cs]}
        result.append({"n": n, "labelled": 1 << len(pairs), "classes": len(classes),
                       "charts": 1 << m, "fibre_histogram": dict(sorted(Counter(chart_fibres.values()).items()))})
    need(collision == {"n": 5, "scores": (1, 1, 2, 3, 3), "H": 9,
                       "codes": [8, 10], "automorphisms": [1, 1],
                       "cyclic_triangles": [3, 3]}, "full scalar collision witness")
    # Hostile K2: the former half-size potential is not the least-squares fit.
    wrong_gradient = Fraction(1, 2)
    need((1 - wrong_gradient) != 0, "old K2 residual divergence must be nonzero")
    need(wrong_gradient ** 2 == Fraction(1, 4), "old tautological energy check")
    # This one is a counting obstruction, not an asymptotic/log float test.
    n = 32
    lower_numerator = 1 << comb(n, 2)
    coarse_capacity = factorial(n) ** 3 * comb(2 * n - 1, n)
    need(lower_numerator > coarse_capacity,
         "scores, H and Aut cannot identify every order-32 class")
    # The planar reflection exchanges two declared port pairs.
    n = 6
    pairs = list(combinations(range(n), 2))
    g = (1, 0, 3, 2, 4, 5)
    transport = [(pairs.index(tuple(sorted((g[i], g[j])))), g[i] > g[j])
                 for i, j in pairs]
    minimum = len(pairs)
    defect_counts = Counter()
    for code in range(1 << len(pairs)):
        defect = sum(((code >> k) & 1) != (((code >> index) & 1) ^ flip)
                     for k, (index, flip) in enumerate(transport))
        minimum = min(minimum, defect)
        defect_counts[defect] += 1
    need(minimum == 2, "two-pair reflection orientation defect")
    need(defect_counts == Counter({2 + 2*j: 512 * comb(6, j) for j in range(7)}),
         "full reflection defect polynomial")
    return {"status": "FINITE-EXACT", "universe": "all labelled tournaments n=0..5",
            "rows": result, "first_score_H_collision": collision,
            "unit_flow_normalization_hostile": "K2 rejects phi=d/(2n)",
            "order32_capacity_test": True,
            "reflection_defect_n6": {"tournaments": 1 << len(pairs), "minimum": minimum,
                                     "histogram": dict(sorted(defect_counts.items()))}}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))

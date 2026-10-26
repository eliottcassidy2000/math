"""Independent binary coset and finite-field incidence audit; standard library."""
import sys as _output_sys
_output_sys.stdout.reconfigure(newline="\n")
import itertools
import json
from collections import Counter
from fractions import Fraction


def need(test, message):
    if not test:
        raise RuntimeError(message)


def span(gens):
    out = {0}
    for g in gens:
        out |= {x ^ g for x in out}
    return frozenset(out)


def dim(s):
    n = len(s).bit_length() - 1
    need(len(s) == 2 ** n, "subspace size")
    return n


def add(a, b):
    return frozenset(x ^ y for x in a for y in b)


G = frozenset(range(64))
INERTIA = [63] + [1 << i for i in range(6)]
MARKINGS = {
    "U2": (1, 2, 4, 5, 6, 3, 7),
    "D2": (1, 2, 4, 6, 3, 7, 5),
    "U4": (1, 2, 4, 5, 7, 3, 6),
    "D4": (1, 2, 4, 7, 3, 6, 5),
}


def kernel(labels):
    out = []
    for x in G:
        syn = 0
        for i in range(6):
            if x >> i & 1:
                syn ^= labels[i + 1]
        if syn == 0:
            out.append(x)
    return frozenset(out)


def reps(s):
    return sorted({min(x ^ y for y in s) for x in G})


def count_fixed(s, twist):
    # Literally act on cosets; no rank/cardinality formula here.
    return sum(min(x ^ twist ^ y for y in s) == x for x in reps(s))


def audit_orbit(h, k, indices):
    inertia = span(INERTIA[i] for i in indices)
    r, a, b = dim(inertia), dim(inertia & h), dim(inertia & k)
    d = r - a - b
    hi, ki = add(h, inertia), add(k, inertia)
    image = Counter()
    for x in reps(inertia):
        image[(min(x ^ y for y in hi), min(x ^ y for y in ki))] += 1
    need(len(image) == len(reps(hi)) * len(reps(ki)), "surjectivity")
    need(set(image.values()) == {2 ** d}, "missing sheets")
    distribution = Counter()
    for twist in G:
        x = count_fixed(inertia, twist)
        y = count_fixed(hi, twist) * count_fixed(ki, twist)
        distribution[x, y] += 1
    fx, fy = 2 ** (6 - r), 2 ** (6 - r - d)
    expected = Counter()
    expected[fx, fy] += 2 ** r
    expected[0, fy] += 2 ** (r + d) - 2 ** r
    expected[0, 0] += 64 - 2 ** (r + d)
    expected = +expected
    need(distribution == expected, "twist distribution")
    mean = lambda f: sum(Fraction(n, 64) * f(x, y)
                         for (x, y), n in distribution.items())
    need(mean(lambda x, y: x) == mean(lambda x, y: y) == 1, "mass")
    need(mean(lambda x, y: abs(x-y)) == 2 * (1-Fraction(1, 2 ** d)), "absolute defect")
    need(mean(lambda x, y: (x-y)**2) == fx-fy, "square defect")
    return {"zero_labels": list(indices), "r": r, "a": a, "b": b,
            "missing_bits": d, "geometric_fibres": [fx, fy],
            "twist_distribution": [dict(X=x, Y=y, multiplicity=n)
                                   for (x, y), n in sorted(distribution.items())],
            "mean_absolute_defect": str(mean(lambda x, y: abs(x-y))),
            "mean_square_defect": str(mean(lambda x, y: (x-y)**2))}


def projective_points(p):
    for x in range(p):
        for z in range(p):
            yield (1, x, z)
    for z in range(p):
        yield (0, 1, z)
    yield (0, 0, 1)


def zero_labels(point, p):
    e, x, z = point
    return tuple(i for i, l in enumerate((e, e+x, e-x, e+z, e-z, e+x+z, e-x-z)) if l % p == 0)


def main():
    kernels = sorted({kernel(row) for row in itertools.permutations(range(1, 8))}, key=lambda s: tuple(sorted(s)))
    need(len(kernels) == 30, "all Fano kernels")
    pairs = [(h, k) for h, k in itertools.combinations(kernels, 2) if h & k == {0}]
    need(len(pairs) == 120, "complementary Fano pairs")
    named = {name: kernel(row) for name, row in MARKINGS.items()}
    cases, census = [], []
    for p in (3, 5, 7, 11, 13):
        strata = Counter(zero_labels(pt, p) for pt in projective_points(p))
        triples = sorted(s for s in strata if len(s) == 3)
        need(triples == ([(0,1,2),(0,3,4),(0,5,6),(1,3,6),(2,4,5)] if p == 3 else [(0,1,2),(0,3,4),(0,5,6)]), "triple census")
        bad = lambda h: any(h & span(INERTIA[i] for i in s) != {0} for s in strata)
        census.append({"prime": p, "stratum_counts": dict(sorted(Counter({r: sum(v for s,v in strata.items() if len(s)==r) for r in range(4)}).items())),
                       "free_kernels": sum(not bad(h) for h in kernels),
                       "both_free_complementary_pairs": sum(not bad(h) and not bad(k) for h,k in pairs)})
        for root in (2,4):
            h,k = named[f"U{root}"], named[f"D{root}"]
            rows = [audit_orbit(h,k,s) for s in sorted(strata)]
            abs_total = sum(Fraction(row["mean_absolute_defect"]) * strata[tuple(row["zero_labels"])] for row in rows)
            sq_total = sum(Fraction(row["mean_square_defect"]) * strata[tuple(row["zero_labels"])] for row in rows)
            need(abs_total == Fraction(28*p-11,4), "global pointwise absolute defect")
            need(sq_total == (112*p-215 if p==3 else 112*p-251), "global pointwise square defect")
            cases.append({"prime":p,"norm_seven_root":root,"strata":rows,
                          "sum_pointwise_mean_absolute_defect":str(abs_total),
                          "sum_pointwise_mean_square_defect":str(sq_total)})
    # Every actual stratum, every complementary Fano pair, independently of freeness.
    universal_cases = 0
    all_strata = set(zero_labels(pt,p) for p in (3,5) for pt in projective_points(p))
    for h,k in pairs:
        for s in all_strata:
            audit_orbit(h,k,s)
            universal_cases += 1
    return {"status":"FINITE-EXACT", "universe":{"markings":5040,"Fano_kernels":30,
            "complementary_pairs":120,"primes":[3,5,7,11,13],
            "all_pair_stratum_twist_checks":universal_cases*64},
            "census":census,"named_pairs":cases}


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True))

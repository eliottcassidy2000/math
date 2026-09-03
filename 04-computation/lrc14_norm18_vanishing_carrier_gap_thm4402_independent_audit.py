#!/usr/bin/env python3
"""Independent exact audit of the norm-18 vanishing carrier-gap family.

This scratch referee imports no repository mathematics.  It checks the
definitions used in THM-4386/4393 directly with exact Fraction arithmetic and
keeps every theorem gate live under ``python -O``.
"""

from fractions import Fraction
from math import gcd


R = Fraction(3, 14)
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def l1(a):
    return sum(abs(x) for x in a)


def lin_add(*polys):
    """Add affine polynomials encoded as (constant, m-coefficient)."""
    return (sum(p[0] for p in polys), sum(p[1] for p in polys))


def lin_scale(k, poly):
    return (k * poly[0], k * poly[1])


def lin_eval(poly, m):
    return poly[0] + poly[1] * m


def data(m):
    w = (1, m, 16 * m - 1)
    c = (1, -16, 1)
    n = (0, 1, 16)
    return w, c, n


def owners(w, n):
    return tuple((-pow(wi, -1, 3) * ni) % 3 for wi, ni in zip(w, n))


def carrier_length(w, C):
    w1, w2, w3 = w
    C1, C2, C3 = C
    terms = (
        2 * R / w1,
        2 * R / w2,
        2 * R / w3,
        R / w1 + R / w2 - Fraction(abs(C3), w1 * w2),
        R / w1 + R / w3 - Fraction(abs(C2), w1 * w3),
        R / w2 + R / w3 - Fraction(abs(C1), w2 * w3),
    )
    return max(Fraction(0), min(terms)), terms


def literal_interval_length(w, n):
    intervals = [
        (Fraction(ni, wi) - R / wi, Fraction(ni, wi) + R / wi)
        for wi, ni in zip(w, n)
    ]
    left = max(a for a, _ in intervals)
    right = min(b for _, b in intervals)
    return max(Fraction(0), right - left), intervals, left, right


def shortest_relations(w, cutoff=18):
    """Enumerate every nonzero relation of l1 norm <= cutoff."""
    found = []
    w1, w2, w3 = w
    for a in range(-cutoff, cutoff + 1):
        for b in range(-cutoff, cutoff + 1):
            numerator = -(a * w1 + b * w2)
            if numerator % w3:
                continue
            d = numerator // w3
            u = (a, b, d)
            if u == (0, 0, 0) or l1(u) > cutoff:
                continue
            found.append(u)
    minimum = min(map(l1, found)) if found else None
    return minimum, sorted(u for u in found if l1(u) == minimum)


def main():
    # Dependency-free symbolic audit in Z[m].  Pairs encode constant and
    # m-coefficients.  This proves the relation and explicit lift identities
    # for all m, before the hostile finite controls below.
    w_poly = ((1, 0), (0, 1), (-1, 16))
    c = (1, -16, 1)
    n = (0, 1, 16)
    relation_poly = lin_add(*(lin_scale(ci, wi) for ci, wi in zip(c, w_poly)))
    cross_poly = (
        lin_add(lin_scale(n[2], w_poly[1]), lin_scale(-n[1], w_poly[2])),
        lin_add(lin_scale(n[0], w_poly[2]), lin_scale(-n[2], w_poly[0])),
        lin_add(lin_scale(n[1], w_poly[0]), lin_scale(-n[0], w_poly[1])),
    )
    require(relation_poly == (0, 0), "symbolic relation in Z[m]")
    require(cross_poly == tuple((ci, 0) for ci in c), "symbolic raw lift in Z[m]")

    # Uniform finite coefficient audit for every possible relation of norm
    # at most 17.  Its affine equation is A+B*m=0 with
    # A=a-d and B=b+16d.  There is no root m>=17 congruent to 5 mod 6.
    shorter_roots = []
    for a in range(-17, 18):
        for b in range(-17, 18):
            for d in range(-17, 18):
                u = (a, b, d)
                if u == (0, 0, 0) or l1(u) > 17:
                    continue
                A, B = a - d, b + 16 * d
                if B == 0:
                    require(A != 0, f"unexpected constant zero relation {u}")
                elif (-A) % B == 0:
                    root = (-A) // B
                    if root >= 17 and root % 6 == 5:
                        shorter_roots.append((u, root))
    require(shorter_roots == [], f"uniform shorter relations: {shorter_roots}")

    # The literal intersection is the third interval for every m>=17.  These
    # are the three cleared-denominator containment numerators, all increasing.
    containment_numerators = ((11, 45), (-17, 45), (-230, 48))
    for poly in containment_numerators:
        require(poly[1] > 0 and lin_eval(poly, 17) > 0,
                f"uniform interval containment: {poly}")

    test_ms = [6 * t + 5 for t in range(2, 1002)]
    test_ms += [6_000_005, 60_000_005]
    sample_rows = []
    for m in test_ms:
        w, c, n = data(m)
        W = w[2]

        # Exact admissibility and relation predicates of THM-4393.
        require(1 < m < W, f"positive distinct sorted speeds at m={m}")
        require(all(x % 2 == 1 for x in w), f"odd speeds at m={m}")
        require(all(x % 3 != 0 for x in w), f"ternary-unit speeds at m={m}")
        require(gcd(gcd(*w[:2]), w[2]) == 1, f"primitive speed triple at m={m}")
        require(dot(c, w) == 0, f"relation at m={m}")
        require(gcd(gcd(abs(c[0]), abs(c[1])), abs(c[2])) == 1,
                f"primitive relation at m={m}")
        require(all(ci % 3 != 0 for ci in c), f"full ternary support at m={m}")
        require(l1(c) == 18, f"relation norm at m={m}")

        # Explicit raw-carrier lift, zero defect, and distinct owners.
        require(cross(w, n) == c, f"raw carrier lift at m={m}")
        require(dot(c, n) == 0, f"zero defect at m={m}")
        require(owners(w, n) == (0, 1, 2), f"owner permutation at m={m}")
        require(all(Ci % 3 != 0 for Ci in c), f"live carrier residue at m={m}")

        # Independent raw formula and literal interval intersection agree.
        formula_length, terms = carrier_length(w, c)
        literal_length, intervals, left, right = literal_interval_length(w, n)
        target = Fraction(3, 7 * W)
        require(formula_length == target, f"carrier formula at m={m}")
        require(literal_length == target, f"literal component at m={m}")
        require(left == Fraction(221, 14 * W), f"left endpoint at m={m}")
        require(right == Fraction(227, 14 * W), f"right endpoint at m={m}")
        require(terms[2] == target, f"third-width active at m={m}")
        require(all(term >= target for term in terms), f"active minimum at m={m}")
        require(Fraction(0) < left < right < Fraction(1),
                f"physical-circle component at m={m}")

        # Brute-force hostile control: the entire integer relation lattice,
        # not only the ternary-unit/full-support sublattice, has l1 minimum 18.
        minimum, minimizers = shortest_relations(w)
        require(minimum == 18, f"shortest relation norm at m={m}: {minimum}")
        require(c in minimizers and tuple(-x for x in c) in minimizers,
                f"fixed minimizer at m={m}")

        if m in (17, 101, 1001, 6_000_005, 60_000_005):
            sample_rows.append((m, W, target, len(minimizers)))

    # A compact paper proof matching the uniform coefficient audit is:
    #   (a-d)+(b+16d)m=0.
    # If b+16d != 0 and ||(a,b,d)||_1 <= 16, then m <= |a-d| <= 16.
    # If b+16d = 0, d != 0 already gives |b|+|d| >= 17; the relation
    # further forces a=d and gives norm at least 18.  If d=0 it is trivial.
    # Every relation has even norm because all three speeds are odd.
    print("LRC14 NORM-18 VANISHING RAW-CARRIER GAP -- INDEPENDENT EXACT AUDIT")
    print("family: m=6t+5, t>=2; w=(1,m,16m-1); c=C=(1,-16,1)")
    print("symbolic relation: 1-16m+(16m-1)=0; ||c||_1=18")
    print("explicit lift: n=(0,1,16); w cross n=C; owners=(0,1,2); defect=0")
    print("literal component: (221/[14W],227/[14W]), W=16m-1")
    print("exact component length: 3/[7(16m-1)]")
    print("analytic minimality: no nonzero relation has l1<=16; odd l1 is impossible")
    print(f"hostile scan: {len(test_ms)} admissible m values; full relation minimum=18 throughout")
    for m, W, length, multiplicity in sample_rows:
        print(f"sample m={m} W={W} L={length} shortest_relation_vectors={multiplicity}")
    print("conclusion: infimum of positive live-carrier lengths in this minimal norm-18 shell is 0")
    print("scope: refutes only a uniform per-component lower quantum; not a total-comb bound or LRC(14)")
    print(f"checks={CHECKS}")
    print("PASS")


if __name__ == "__main__":
    main()

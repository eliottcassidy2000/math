#!/usr/bin/env python3
"""Exact verifier for the THM-4387 coefficient-one defect-three boundary.

The script uses only Python's standard library and exact Fraction arithmetic.
Every correctness check is explicit, so ``python -O`` retains the audit.
"""

from fractions import Fraction as F
from itertools import permutations, product
from math import gcd


# Exactly the coefficient-one |coefficient|_1=16 slice.  The additional
# primitive full-support shapes (2,7,7) and (4,5,7) are deliberately outside
# this verifier.
R = F(3, 14)
PATTERNS = ((1, 14), (2, 13), (4, 11), (5, 10), (7, 8))
EXPECTED = {
    (1, 14): {
        "target": F(4, 133), "winner": (1, 5, 19), "cutoff": 205,
        "presentations": 298, "triples": 149,
        "parts": {-3: F(2, 133), 0: F(0), 3: F(2, 133)},
    },
    (2, 13): {
        "target": F(228, 11165), "winner": (1, 11, 145), "cutoff": 409,
        "presentations": 1295, "triples": 1295,
        "parts": {-3: F(48, 11165), 0: F(12, 1015), 3: F(48, 11165)},
    },
    (4, 11): {
        "target": F(218, 10465), "winner": (7, 13, 115), "cutoff": 372,
        "presentations": 1270, "triples": 1270,
        "parts": {-3: F(31, 10465), 0: F(12, 805), 3: F(31, 10465)},
    },
    (5, 10): {
        "target": F(16, 715), "winner": (1, 11, 65), "cutoff": 335,
        "presentations": 1105, "triples": 1104,
        "parts": {-3: F(23, 5005), 0: F(6, 455), 3: F(23, 5005)},
    },
    (7, 8): {
        "target": F(36, 1309), "winner": (1, 11, 85), "cutoff": 257,
        "presentations": 792, "triples": 792,
        "parts": {-3: F(5, 1309), 0: F(26, 1309), 3: F(5, 1309)},
    },
}

CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def ceil_fraction(x):
    return -((-x.numerator) // x.denominator)


def nint(x):
    z = x + F(1, 2)
    return z.numerator // z.denominator


def component_length(p, b, q, m, t, delta, k):
    """Length of the (delta,k) lift class in z=y-n_p/p coordinates."""
    centers = (
        F(0),
        -F(k, p * b),
        -F(t * (m * k + p * delta), p * q),
    )
    radii = (R / p, R / b, R / q)
    lower = max(c - a for c, a in zip(centers, radii))
    upper = min(c + a for c, a in zip(centers, radii))
    return max(F(0), upper - lower)


def formula_parts(p, b, q, h, m, s, t, audit=False):
    del h, s  # They select the presentation; geometry then uses its m,t,p,b,q.
    result = {}
    rows = {}
    for delta in (-3, 0, 3):
        # Endpoint overlap is necessary: |m k+p delta| < R(p+q).
        lo = F(-p * delta, m) - R * F(p + q, m)
        hi = F(-p * delta, m) + R * F(p + q, m)
        first = lo.numerator // lo.denominator + 1
        last = ceil_fraction(hi) - 1
        terms = []
        for k in range(first, last + 1):
            if k % 3 == 0:
                continue
            value = component_length(p, b, q, m, t, delta, k)
            if value:
                terms.append((k, value))
                if audit:
                    check(
                        component_length(p, b, q, m, t, -delta, -k) == value,
                        "reflection symmetry failed",
                    )
        result[delta] = sum((v for _, v in terms), F(0))
        rows[delta] = tuple(terms)
    check(result[-3] == result[3], "defect reflection totals disagree")
    return result, rows


def formula_measure(p, b, q, h, m, s, t, audit=False):
    parts, rows = formula_parts(p, b, q, h, m, s, t, audit)
    return sum(parts.values(), F(0)), parts, rows


def physical_walls(speeds):
    walls = {F(0), F(1)}
    for w in speeds:
        for n in range(-1, w + 2):
            for sign in (-1, 1):
                y = (F(n) + sign * R) / w
                if 0 <= y <= 1:
                    walls.add(y)
    return sorted(walls)


def direct_physical_measure(speeds):
    """Definition-level y-circle cell decomposition, with no relation input."""
    total = F(0)
    walls = physical_walls(speeds)
    for left, right in zip(walls, walls[1:]):
        y = (left + right) / 2
        owners = []
        eligible = True
        for w in speeds:
            n = nint(w * y)
            e = w * y - n
            if abs(e) >= R:
                eligible = False
                break
            owners.append((-pow(w, -1, 3) * n) % 3)
        if eligible and set(owners) == {0, 1, 2}:
            total += right - left
    return total, len(walls) - 1


def presentations_below_width(h, m, width):
    pool = [w for w in range(1, width, 2) if w % 3]
    for p in pool:
        for q in pool:
            for s, t in product((-1, 1), repeat=2):
                numerator = s * h * p + t * q
                if numerator <= 0 or numerator % m:
                    continue
                b = numerator // m
                if b >= width or b % 2 == 0 or b % 3 == 0:
                    continue
                if len({p, b, q}) != 3 or gcd(gcd(p, b), q) != 1:
                    continue
                yield p, b, q, s, t


def bulk(h, m):
    # (2/3) times [4R^2/m + 2*(2R^2/(mh))].
    return F(6 * (h + 1), 49 * m * h)


def width_cutoff(h, m, target):
    width = 1
    while bulk(h, m) + F(18, 7 * width) >= target:
        width += 1
    return width


def egcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


def endpoint_coordinates(p, b, q, h, m, s, t, delta, k):
    # b*u-p*v=1 gives a section of k=b*n_p-p*n_b.
    g0, u, x = egcd(b, p)
    check(g0 == 1, "p,b section is not primitive")
    v = -x
    np = u * k
    nb = v * k
    nq = t * (m * nb - s * h * np - delta)
    check(m * nb - s * h * np - t * nq == delta, "defect section failed")
    check(b * np - p * nb == k, "k section failed")

    g = gcd(p, q)
    P, Q = p // g, q // g
    gg, alpha, beta = egcd(P, Q)
    check(gg == 1, "reduced endpoints not coprime")
    determinant = q * np - p * nq
    j = determinant // g
    torsion = (alpha * np + beta * nq) % g
    check(determinant % g == 0, "endpoint determinant seam failed")
    check(determinant == t * (m * k + p * delta), "determinant identity failed")
    return j, torsion, (np, nb, nq)


def audit_seams():
    # M=1: all three defects collide in the determinant coordinate, while the
    # endpoint torsion coordinate separates them.
    all_state_examples = (
        (2, 13, 13, 7, 65, 1, 1),
        (4, 11, 11, 1, 55, -1, 1),
    )
    seam_rows = []
    for h, m, p, b, q, s, t in all_state_examples:
        g = gcd(p, q)
        check(g == m, "M=1 seam example changed")
        coordinates = []
        for delta in (-3, 0, 3):
            # t*(m/g*k+p/g*delta)=1.
            k = t - (p // g) * delta
            check(k % 3 != 0, "owner gate lost in M=1 seam")
            coordinates.append(endpoint_coordinates(p, b, q, h, m, s, t, delta, k))
        check(len({row[0] for row in coordinates}) == 1, "M=1 determinant did not collide")
        check(len({row[1] for row in coordinates}) == 3, "M=1 torsion did not separate")
        seam_rows.append(((h, m), g, tuple((r[0], r[1]) for r in coordinates)))

    # M=2: delta=-3 and +3 collide in determinant, but not in torsion.
    for h, m, p, b, q, s, t in (
        (1, 14, 7, 5, 77, -1, 1),
        (5, 10, 25, 13, 5, 1, 1),
    ):
        g = gcd(p, q)
        check(m // g == 2, "M=2 seam example changed")
        coordinates = []
        for delta in (-3, 3):
            numerator = t - (p // g) * delta
            check(numerator % (m // g) == 0, "M=2 collision has no integral k")
            k = numerator // (m // g)
            check(t * ((m // g) * k + (p // g) * delta) == 1,
                  "M=2 collision solve failed")
            check(k % 3 != 0, "owner gate lost in M=2 seam")
            coordinates.append(endpoint_coordinates(p, b, q, h, m, s, t, delta, k))
        check(coordinates[0][0] == coordinates[1][0], "M=2 determinant did not collide")
        check(coordinates[0][1] != coordinates[1][1], "M=2 torsion did not separate")
        seam_rows.append(((h, m), g, tuple((r[0], r[1]) for r in coordinates)))
    return seam_rows


def audit_boundary_witnesses():
    rows = (
        (1, 14, 19, 1, 5, 1, -1, F(0)),
        (2, 13, 5, 1, 23, -1, 1, F(34, 161)),
        (4, 11, 19, 7, 1, 1, 1, F(213, 1862)),
        (5, 10, 13, 7, 5, 1, 1, F(509, 1274)),
        (7, 8, 23, 17, 25, 1, -1, F(2189, 5474)),
    )
    out = []
    for h, m, p, b, q, s, t, y in rows:
        if (h, m) == (1, 14):
            # Use the midpoint of the winner's positive delta=3 component.
            parts, terms = formula_parts(p, b, q, h, m, s, t)
            check(parts[3] > 0 and terms[3], "(1,14) has no defect-three component")
            delta, k = 3, terms[3][0][0]
            # Reconstruct a point from the common intersection midpoint.
            g0, u, x = egcd(b, p)
            check(g0 == 1, "winner section failed")
            np, nb = u * k, -x * k
            nq = t * (m * nb - s * h * np - delta)
            lows = (F(np, p) - R / p, F(nb, b) - R / b, F(nq, q) - R / q)
            highs = (F(np, p) + R / p, F(nb, b) + R / b, F(nq, q) + R / q)
            y = (max(lows) + min(highs)) / 2
            y -= y.numerator // y.denominator
        ns = tuple(nint(w * y) for w in (p, b, q))
        errors = tuple(w * y - n for w, n in zip((p, b, q), ns))
        owners = tuple((-pow(w, -1, 3) * n) % 3 for w, n in zip((p, b, q), ns))
        delta = m * ns[1] - s * h * ns[0] - t * ns[2]
        check(all(abs(e) < R for e in errors), "boundary witness lost eligibility")
        check(set(owners) == {0, 1, 2}, "boundary witness lost distinct owners")
        check(abs(delta) == 3, "boundary witness lost nonzero defect")
        out.append(((h, m), (p, b, q), y, ns, errors, owners, delta))
    return out


def coefficient_vectors(h, m):
    """All signed coefficient rows of magnitude multiset {m,h,1}."""
    rows = set()
    for magnitudes in set(permutations((m, h, 1))):
        for signs in product((-1, 1), repeat=3):
            row = tuple(a * sign for a, sign in zip(magnitudes, signs))
            if next(z for z in row if z) < 0:
                row = tuple(-z for z in row)
            rows.add(row)
    return rows


def positive_unit_intersection_rays(left, right):
    rays = set()
    for u in coefficient_vectors(*left):
        for v in coefficient_vectors(*right):
            if u == v:
                continue
            cross = (
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            )
            if not all(cross):
                continue
            if not (all(z > 0 for z in cross) or all(z < 0 for z in cross)):
                continue
            cross = tuple(abs(z) for z in cross)
            divisor = gcd(gcd(cross[0], cross[1]), cross[2])
            ray = tuple(sorted(z // divisor for z in cross))
            if len(set(ray)) == 3 and all(z % 2 and z % 3 for z in ray):
                rays.add(ray)
    return tuple(sorted(rays))


def audit_relation_chart_incidence():
    old = ((1, 2), (1, 4), (1, 8), (1, 10), (2, 5),
           (2, 7), (2, 11), (4, 5), (4, 7), (5, 8))
    expected_old = {
        ((1, 14), (1, 4)): ((1, 5, 19),),
        ((1, 14), (5, 8)): ((1, 23, 37),),
        ((2, 13), (1, 10)): ((5, 7, 43),),
        ((2, 13), (2, 5)): ((1, 5, 23),),
        ((2, 13), (4, 7)): ((5, 7, 29), (5, 17, 31)),
        ((4, 11), (1, 8)): ((5, 11, 29),),
        ((4, 11), (1, 10)): ((1, 7, 17), (5, 13, 37)),
        ((4, 11), (2, 5)): ((1, 7, 17), (1, 7, 19)),
        ((4, 11), (4, 7)): ((1, 5, 31), (5, 13, 17)),
        ((5, 10), (1, 4)): ((5, 7, 13),),
        ((5, 10), (2, 7)): ((1, 7, 25), (5, 11, 23)),
        ((5, 10), (2, 11)): ((1, 7, 25), (5, 19, 37)),
        ((5, 10), (5, 8)): ((1, 5, 35), (5, 13, 25), (11, 17, 25)),
    }
    found_old = {}
    for fresh in PATTERNS:
        for previous in old:
            rays = positive_unit_intersection_rays(fresh, previous)
            if rays:
                found_old[(fresh, previous)] = rays
    check(found_old == expected_old, "old/new coefficient-incidence atlas changed")

    expected_new = {
        ((1, 14), (5, 10)): ((5, 23, 47),),
        ((2, 13), (4, 11)): ((1, 17, 47),),
    }
    found_new = {}
    for index, left in enumerate(PATTERNS):
        for right in PATTERNS[index + 1:]:
            rays = positive_unit_intersection_rays(left, right)
            if rays:
                found_new[(left, right)] = rays
    check(found_new == expected_new, "new/new coefficient-incidence atlas changed")

    # At the (1,14) maximizer, the two pure defect-three components are exactly
    # the k'=+/-1 defect-zero components in the older (1,4) chart.
    new_parts, new_terms = formula_parts(19, 1, 5, 1, 14, 1, -1)
    check(new_parts[0] == 0, "new chart unexpectedly has a zero-defect component")
    transported = []
    for delta in (-3, 3):
        for k, value in new_terms[delta]:
            old_k = -k - delta
            old_value = component_length(1, 5, 19, 4, 1, 0, old_k)
            check(old_value == value, "relation-chart component transport failed")
            transported.append((delta, k, old_k, value))
    check(sorted(transported) == [(-3, 4, -1, F(2, 133)), (3, -4, 1, F(2, 133))],
          "relation-chart address map changed")
    return found_old, found_new, tuple(sorted(transported))


def main():
    check(PATTERNS == tuple(
        (h, 15 - h) for h in range(1, 8)
        if h % 3 and (15 - h) % 3 and ((15 - h) - h - 1) % 2 == 0
    ), "coefficient-one first-boundary parity enumeration changed")

    direct_cache = {}
    summaries = []
    all_gcds = {}
    for h, m in PATTERNS:
        expected = EXPECTED[(h, m)]
        cutoff = width_cutoff(h, m, expected["target"])
        check(cutoff == expected["cutoff"], "analytic width cutoff changed")
        check(bulk(h, m) + F(18, 7 * cutoff) < expected["target"], "cutoff does not close tail")
        check(cutoff == 1 or bulk(h, m) + F(18, 7 * (cutoff - 1)) >= expected["target"],
              "cutoff is not minimal for the stated envelope")

        fibres = {}
        presentations = 0
        gcds = set()
        for p, b, q, s, t in presentations_below_width(h, m, cutoff):
            presentations += 1
            check(m * b == s * h * p + t * q, "enumerator emitted false relation")
            check(gcd(p, b) == 1, "primitive relation did not force gcd(p,b)=1")
            g = gcd(p, q)
            gcds.add(g)
            check(m % g == 0 and gcd(p, m) == g, "gcd seam identity failed")

            mu, parts, terms = formula_measure(p, b, q, h, m, s, t, audit=True)
            key = tuple(sorted((p, b, q)))
            if key not in direct_cache:
                direct_cache[key] = direct_physical_measure(key)
            check(mu == direct_cache[key][0], "lift formula disagrees with physical circle")
            if key in fibres:
                check(fibres[key][0] == mu, "multiple presentations disagree")
                fibres[key][2].append((p, b, q, s, t))
            else:
                fibres[key] = [mu, parts, [(p, b, q, s, t)], terms]

            # The normalized endpoint determinant has the same owner gate.
            for delta, rows in terms.items():
                for k, value in rows:
                    del value
                    j = t * (m * k + p * delta) // g
                    check((j % 3 != 0) == (k % 3 != 0), "normalized endpoint owner gate failed")

        check(presentations == expected["presentations"], "finite presentation count changed")
        check(len(fibres) == expected["triples"], "finite triple count changed")
        ranking = sorted(((v[0], key) for key, v in fibres.items()), reverse=True)
        check(ranking[0] == (expected["target"], expected["winner"]), "sharp maximum changed")
        check(sum(mu == expected["target"] for mu, _ in ranking) == 1, "sharp maximum not unique")
        winner_parts = fibres[expected["winner"]][1]
        check(winner_parts == expected["parts"], "winner defect decomposition changed")
        all_gcds[(h, m)] = tuple(sorted(gcds))
        summaries.append((
            (h, m), expected["target"], expected["winner"], bulk(h, m), cutoff,
            presentations, len(fibres), winner_parts,
        ))

    check(all_gcds == {
        (1, 14): (1, 7),
        (2, 13): (1, 13),
        (4, 11): (1, 11),
        (5, 10): (1, 5),
        (7, 8): (1,),
    }, "finite gcd seam census changed")

    # The sole multiple-presentation fibre in the asymmetric rows.
    duplicate = []
    for row in presentations_below_width(5, 10, EXPECTED[(5, 10)]["cutoff"]):
        if tuple(sorted(row[:3])) == (5, 17, 35):
            duplicate.append(row)
    check(sorted(duplicate) == [(17, 5, 35, 1, -1), (35, 17, 5, 1, -1)],
          "(5,10) duplicate fibre changed")

    witnesses = audit_boundary_witnesses()
    seams = audit_seams()
    old_incidence, new_incidence, transported = audit_relation_chart_incidence()

    print("status=PASS")
    print(f"patterns={PATTERNS}")
    for summary in summaries:
        print(
            "pattern=%s max=%s winner=%s bulk=%s width_cutoff=%s "
            "presentations=%s triples=%s parts=%s" % summary
        )
    print(f"gcd_seams={all_gcds}")
    print(f"duplicate_fibre_5_10={sorted(duplicate)}")
    print(f"boundary_witnesses={witnesses}")
    print(f"endpoint_seam_controls={seams}")
    print(f"old_new_relation_rays={old_incidence}")
    print(f"new_new_relation_rays={new_incidence}")
    print(f"winner_chart_transport={transported}")
    print(f"direct_physical_triples={len(direct_cache)}")
    print(f"explicit_checks={CHECKS}")
    print("coefficient_scope=coefficient_one_only; omitted_shapes=((2,7,7),(4,5,7))")
    print("scope=LRC(14)_OPEN")


if __name__ == "__main__":
    main()

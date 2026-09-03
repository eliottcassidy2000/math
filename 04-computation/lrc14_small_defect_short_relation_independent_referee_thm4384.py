from fractions import Fraction as F
from itertools import permutations
from math import gcd, isqrt
import hashlib


R = F(3, 14)
PATTERNS = ((1, 2), (1, 4), (1, 8), (1, 10),
            (2, 5), (2, 7), (2, 11),
            (4, 5), (4, 7), (5, 8))
SIGNS = ((1, 1), (1, -1), (-1, 1))
EXPECTED = {
    (1, 2): (F(6, 77), 13, {(1, 5, 11, -1, 1)}),
    (1, 4): (F(12, 301), 31, {(1, 11, 43, 1, 1)}),
    (1, 8): (F(6, 287), 58, {(1, 5, 41, -1, 1)}),
    (1, 10): (F(6, 343), 65, {(1, 5, 49, 1, 1)}),
    (2, 5): (F(12, 371), 75, {(1, 11, 53, 1, 1)}),
    (2, 7): (F(12, 469), 75, {(67, 19, 1, 1, -1)}),
    (2, 11): (F(6, 371), 132, {(1, 5, 53, 1, 1)}),
    (4, 5): (F(82, 2597), 84,
             {(7, 5, 53, -1, 1), (53, 41, 7, 1, -1)}),
    (4, 7): (F(12, 497), 97, {(5, 13, 71, 1, 1)}),
    (5, 8): (F(58, 2765), 117, {(5, 13, 79, 1, 1)}),
}
EXPECTED_TABLE_COUNTS = {
    (1, 2): (26, 13), (1, 4): (62, 31), (1, 8): (116, 58),
    (1, 10): (130, 65), (2, 5): (75, 75), (2, 7): (75, 74),
    (2, 11): (132, 131), (4, 5): (84, 84), (4, 7): (97, 97),
    (5, 8): (117, 116),
}


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_unit(w):
    return w > 0 and w % 2 == 1 and w % 3 != 0


def e3(t):
    # E(t) = sum_{k>=1, 3 does not divide k} (t-k)_+ - t^2/3.
    rho = t - 3 * (t // 3)
    if rho <= 1:
        return -rho * rho / 3
    if rho <= 2:
        return rho - 1 - rho * rho / 3
    return 2 * rho - 3 - rho * rho / 3


def measure_series(m, p, q):
    a = R * abs(q - p) / m
    bcap = R * (p + q) / m
    total = F(0)
    for k in range(1, (bcap.numerator + bcap.denominator - 1) // bcap.denominator):
        if k % 3:
            total += 2 * F(m, p * q) * (max(F(0), bcap - k) - max(F(0), a - k))
    return total


def measure_quad(m, p, q):
    a = R * abs(q - p) / m
    bcap = R * (p + q) / m
    return F(6, 49 * m) + F(2 * m, p * q) * (e3(bcap) - e3(a))


def add_circle_interval(out, center, radius):
    center -= center // 1
    lo, hi = center - radius, center + radius
    if lo < 0:
        out.append((F(0), hi))
        out.append((lo + 1, F(1)))
    elif hi > 1:
        out.append((lo, F(1)))
        out.append((F(0), hi - 1))
    else:
        out.append((lo, hi))


def merged(intervals):
    intervals = sorted(intervals)
    out = []
    for lo, hi in intervals:
        if not out or lo > out[-1][1]:
            out.append([lo, hi])
        elif hi > out[-1][1]:
            out[-1][1] = hi
    return [(lo, hi) for lo, hi in out]


def danger(w, j):
    out = []
    radius = F(1, 14 * w)
    for n in range(w):
        add_circle_interval(out, F(n, w) - F(j, 3), radius)
    return merged(out)


def intersect_two(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            out.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return out


def physical_measure(speeds):
    cache = {(w, j): danger(w, j) for w in speeds for j in range(3)}
    pieces = []
    for perm in permutations(range(3)):
        cur = cache[(speeds[0], perm[0])]
        cur = intersect_two(cur, cache[(speeds[1], perm[1])])
        cur = intersect_two(cur, cache[(speeds[2], perm[2])])
        pieces.extend(cur)
    return sum((hi - lo for lo, hi in merged(pieces)), F(0))


def nearest_integer(z):
    return (z + F(1, 2)) // 1


def criterion_cell_audit(h, m, p, b, q, s, t):
    walls = {F(0), F(1)}
    for w in (p, b, q):
        for n in range(w):
            for sign in (-1, 1):
                point = (F(n, 1) + sign * R) / w
                point -= point // 1
                walls.add(point)
    walls = sorted(walls)
    checked = 0
    for lo, hi in zip(walls, walls[1:]):
        y = (lo + hi) / 2
        ns = tuple(nearest_integer(w * y) for w in (p, b, q))
        es = tuple(w * y - n for w, n in zip((p, b, q), ns))
        eligible = tuple(abs(e) < R for e in es)
        owners = tuple((-pow(w, -1, 3) * n) % 3
                       for w, n in zip((p, b, q), ns))
        physical = all(eligible) and len(set(owners)) == 3
        delta = s * h * ns[0] + t * ns[2] - m * ns[1]
        k = b * ns[0] - p * ns[1]
        cyclic = eligible[0] and eligible[2] and delta == 0 and k % 3 != 0
        check(physical == cyclic,
              f"cell criterion {(h,m,p,b,q,s,t,y,physical,cyclic,delta,k)}")
        if physical:
            check(delta == 0, f"owner did not kill defect {(h,m,p,b,q,y,delta)}")
        checked += 1
    return checked


def presentations_below(h, m, product_limit):
    # h=1 has symmetric coefficient-one endpoints, so use p<q once.
    rows = []
    for p in range(1, product_limit + 1, 2):
        if not is_unit(p):
            continue
        qmax = product_limit // p
        for q in range(1, qmax + 1, 2):
            if not is_unit(q) or q == p or (h == 1 and q < p):
                continue
            for s, t in SIGNS:
                numer = s * h * p + t * q
                if numer <= 0 or numer % m:
                    continue
                b = numer // m
                if not is_unit(b) or b in (p, q):
                    continue
                if gcd(gcd(p, b), q) != 1:
                    continue
                rows.append((p, b, q, s, t, measure_series(m, p, q)))
    return rows


def canonical_coefficient_vectors(h, m):
    vectors = set()
    for mag in set(permutations((m, h, 1))):
        for signs in ((a, b, c) for a in (-1, 1) for b in (-1, 1) for c in (-1, 1)):
            v = tuple(a * b for a, b in zip(mag, signs))
            if v[0] < 0:
                v = tuple(-x for x in v)
            vectors.add(v)
    return sorted(vectors)


def positive_primitive_cross(u, v):
    w = (u[1] * v[2] - u[2] * v[1],
         u[2] * v[0] - u[0] * v[2],
         u[0] * v[1] - u[1] * v[0])
    if 0 in w or not (all(x > 0 for x in w) or all(x < 0 for x in w)):
        return None
    if w[0] < 0:
        w = tuple(-x for x in w)
    d = gcd(gcd(w[0], w[1]), w[2])
    return tuple(x // d for x in w)


def valid_speed_ray(w):
    return w is not None and len(set(w)) == 3 and all(is_unit(x) for x in w)


def presentation_overlap_audit():
    vectors = {pattern: canonical_coefficient_vectors(*pattern) for pattern in PATTERNS}
    duplicate_rays = {}
    for pattern in PATTERNS:
        vecs = vectors[pattern]
        rays = set()
        for i, u in enumerate(vecs):
            for v in vecs[i + 1:]:
                w = positive_primitive_cross(u, v)
                if valid_speed_ray(w):
                    rays.add(tuple(sorted(w)))
        duplicate_rays[pattern] = rays
    overlaps = {}
    for i, left in enumerate(PATTERNS):
        for right in PATTERNS[i + 1:]:
            rays = set()
            for u in vectors[left]:
                for v in vectors[right]:
                    w = positive_primitive_cross(u, v)
                    if valid_speed_ray(w):
                        rays.add(tuple(sorted(w)))
            if rays:
                overlaps[(left, right)] = rays
    return vectors, duplicate_rays, overlaps


def main():
    derived_patterns = tuple((h, m) for h in range(1, 14) for m in range(1, 14)
                             if m >= h + 1 and m + h <= 13
                             and gcd(3, h * m) == 1 and (m + h) % 2 == 1)
    check(derived_patterns == PATTERNS, f"pattern derivation {derived_patterns}")
    print("PATTERNS", PATTERNS)

    summaries = []
    all_physical_checks = 0
    for h, m in PATTERNS:
        expected_mu, expected_count, expected_maximizers = EXPECTED[(h, m)]
        base = F(6, 49 * m)
        check(expected_mu > base, f"maximum not above base {(h, m)}")
        threshold = F(2 * m, 3) / (expected_mu - base)
        exceptional_limit = threshold.numerator // threshold.denominator
        exceptional = presentations_below(h, m, exceptional_limit)
        max_mu = max(row[-1] for row in exceptional)
        maximizers = [row for row in exceptional if row[-1] == max_mu]
        got_maximizers = {row[:-1] for row in maximizers}
        check(len(exceptional) == expected_count,
              f"exceptional count {(h,m,len(exceptional),expected_count)}")
        check(max_mu == expected_mu, f"maximum {(h,m,max_mu,expected_mu)}")
        check(got_maximizers == expected_maximizers,
              f"maximizers {(h,m,got_maximizers,expected_maximizers)}")
        raw_presentations = len(exceptional) * (2 if h == 1 else 1)
        triple_count = len({tuple(sorted(row[:3])) for row in exceptional})
        check((raw_presentations, triple_count) == EXPECTED_TABLE_COUNTS[(h, m)],
              f"table counts {(h,m,raw_presentations,triple_count)}")
        check(all(measure_series(m, row[0], row[2]) == measure_quad(m, row[0], row[2])
                  for row in exceptional), f"series/quadrature {(h, m)}")
        for p, b, q, s, t, mu in exceptional:
            check(m * b == s * h * p + t * q, "enumerated relation")
            check(gcd(p, b) == 1, "primitive cyclic pair")
            residues = (m * b % 3, (-s * h * p) % 3, (-t * q) % 3)
            check(len(set(residues)) == 1, f"three-term mod3 coincidence {residues}")
            check(mu <= expected_mu, "finite sharp bound")
            direct = physical_measure((p, b, q))
            check(direct == mu,
                  f"all-physical mismatch {(h,m,p,b,q,s,t,mu,direct)}")
            all_physical_checks += 1
        summaries.append((h, m, max_mu, threshold, exceptional_limit,
                          len(exceptional), maximizers))
        print("SECTOR", h, m, "MAX", max_mu, "THRESHOLD", threshold,
              "EXCEPTIONAL_PRODUCT_LE", exceptional_limit,
              "NORMALIZED_COUNT", len(exceptional),
              "RAW_PRESENTATIONS", raw_presentations, "TRIPLES", triple_count,
              "MAXIMIZERS", maximizers)

    # Definition-level physical controls: every maximizer and the next two
    # exceptional presentations by measure in each sector.
    criterion_controls = 0
    criterion_cells = 0
    for h, m, max_mu, threshold, exceptional_limit, count, maximizers in summaries:
        rows = presentations_below(h, m, exceptional_limit)
        controls = sorted(rows, key=lambda row: (-row[-1], row[:-1]))[:3]
        for p, b, q, s, t, mu in controls:
            direct = physical_measure((p, b, q))
            check(direct == mu, f"physical mismatch {(h,m,p,b,q,s,t,mu,direct)}")
            criterion_cells += criterion_cell_audit(h, m, p, b, q, s, t)
            criterion_controls += 1
            print("PHYSICAL", h, m, p, b, q, s, t, mu)

    vectors, duplicate_rays, overlaps = presentation_overlap_audit()
    print("COEFFICIENT_VECTOR_COUNTS", tuple((p, len(vectors[p])) for p in PATTERNS))
    print("DUPLICATE_PRESENTATION_RAYS",
          tuple((p, tuple(sorted(duplicate_rays[p]))) for p in PATTERNS if duplicate_rays[p]))
    print("CROSS_SECTOR_OVERLAPS",
          tuple((pair, tuple(sorted(rays))) for pair, rays in sorted(overlaps.items())))

    # Independently check the requested defect-three hostile boundary.
    h, m, p, b, q, s, t = 2, 13, 5, 1, 23, -1, 1
    y = F(34, 161)
    ns = (1, 0, 5)
    es = (p * y - ns[0], b * y - ns[1], q * y - ns[2])
    delta = s * h * ns[0] + t * ns[2] - m * ns[1]
    owners = tuple((-pow(w, -1, 3) * n) % 3 for w, n in zip((p, b, q), ns))
    check(m * b == s * h * p + t * q, "hostile relation")
    check(all(abs(e) < R for e in es), "hostile eligibility")
    check(owners == (1, 0, 2) and len(set(owners)) == 3, "hostile owners")
    check(delta == 3, f"hostile defect sign/value {delta}")
    check(m + h == 15, "hostile boundary")
    print("HOSTILE", (h, m, p, b, q, s, t), "y", y, "n", ns,
          "errors", es, "owners", owners, "delta", delta)
    semantic_lines = []
    for h, m, mu, threshold, limit, count, maximizers in summaries:
        semantic_lines.append(f"{h},{m}|{mu}|{threshold}|{limit}|{count}|{maximizers}")
    semantic_lines.append(repr(tuple((p, tuple(sorted(duplicate_rays[p])))
                                     for p in PATTERNS if duplicate_rays[p])))
    semantic_lines.append(repr(tuple((pair, tuple(sorted(rays)))
                                     for pair, rays in sorted(overlaps.items()))))
    digest = hashlib.sha256("\n".join(semantic_lines).encode("ascii")).hexdigest()
    print("TOTAL_EXCEPTIONAL_PRESENTATIONS", sum(item[5] for item in summaries))
    print("TOTAL_ALL_PHYSICAL_CHECKS", all_physical_checks)
    print("TOTAL_CRITERION_CONTROLS", criterion_controls)
    print("TOTAL_CRITERION_CELLS", criterion_cells)
    print("SEMANTIC_SHA256", digest)


if __name__ == "__main__":
    main()

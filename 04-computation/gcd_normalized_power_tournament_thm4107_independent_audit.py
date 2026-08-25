from fractions import Fraction
from itertools import combinations
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def edge(a, b):
    g = gcd(a, b)
    x, y = a // g, b // g
    return x**y > y**x


def exceptional(a, b):
    require(a < b, (a, b))
    return b % a == 0 or 3 * a == 2 * b


def cyclic(a, b, c):
    vertices = (a, b, c)
    degrees = [sum(x != y and edge(x, y) for y in vertices) for x in vertices]
    return degrees == [1, 1, 1] or sorted(degrees) == [1, 1, 1]


def c3_formula(a, b, c):
    eab = exceptional(a, b)
    ebc = exceptional(b, c)
    eac = exceptional(a, c)
    return eab == ebc and eab != eac


def is_dth_power(n, d):
    if d == 1:
        return True
    lo, hi = 1, n
    while lo <= hi:
        mid = (lo + hi) // 2
        value = mid**d
        if value == n:
            return True
        if value < n:
            lo = mid + 1
        else:
            hi = mid - 1
    return False


def distance(v, t):
    x = v * t
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def observer(row, t):
    return min(distance(v, t) for v in row)


def exact_max(row):
    """Independent exact max of min_v ||vt|| on [0,1/2]."""
    boundaries = sorted(
        {Fraction(k, 2 * v) for v in row for k in range(v + 1)}
    )
    best = Fraction(0)
    owners = []
    for left, right in zip(boundaries, boundaries[1:]):
        p = (2 * left + right) / 3
        q = (left + 2 * right) / 3
        lines = []
        for v in row:
            dp, dq = distance(v, p), distance(v, q)
            slope = (dq - dp) / (q - p)
            lines.append((slope, dp - slope * p))
        candidates = {left, right}
        for (s1, b1), (s2, b2) in combinations(lines, 2):
            if s1 != s2:
                t = (b2 - b1) / (s1 - s2)
                if left <= t <= right:
                    candidates.add(t)
        for t in candidates:
            value = observer(row, t)
            if value > best:
                best, owners = value, [t]
            elif value == best:
                owners.append(t)
    return best, tuple(sorted(set(owners)))


def labelled_signature(row):
    return tuple(edge(row[i], row[j]) for i, j in combinations(range(len(row)), 2))


def first_target_denominator(row, target, limit):
    for q in range(2, limit + 1):
        for a in range(1, q):
            if observer(row, Fraction(a, q)) >= target:
                return q, a, observer(row, Fraction(a, q))
    return None


# Classification and tie-freeness in a large finite hostile deck.
pair_count = 0
for a in range(1, 201):
    for b in range(a + 1, 201):
        pair_count += 1
        g = gcd(a, b)
        x, y = a // g, b // g
        require(x**y != y**x, (a, b, x, y))
        require(edge(a, b) == (not exceptional(a, b)), (a, b))

# Exact C3 Boolean criterion, including the first maximum shell.
cycles = []
for a, b, c in combinations(range(1, 81), 3):
    seen = cyclic(a, b, c)
    require(seen == c3_formula(a, b, c), (a, b, c))
    if seen:
        cycles.append((a, b, c))
first_cycles = sorted(cycles, key=lambda z: (z[2], z))[:3]
require(first_cycles == [(2, 5, 6), (3, 5, 6), (4, 5, 6)], first_cycles)

# Rational equality ray: the rationality condition is d=m-n=1.
rational_ratio_candidates = []
for n in range(1, 101):
    for m in range(n + 1, 101):
        if gcd(m, n) != 1:
            continue
        d = m - n
        rational = is_dth_power(m, d) and is_dth_power(n, d)
        require(rational == (d == 1), (m, n, d))
        if rational:
            rational_ratio_candidates.append((m, n))
for n in range(1, 50):
    # Cleared equality packet n(n+1)^n + (n+1)^n = (n+1)^(n+1).
    common = (n + 1) ** n
    terms = (n * common, common, (n + 1) * common)
    require(terms[0] + terms[1] == terms[2], n)
    require(gcd(gcd(*terms[:2]), terms[2]) == common, n)
    require((terms[0] // common, terms[1] // common, terms[2] // common) == (n, 1, n + 1), n)
    require((n + 2) * n - (n + 1) ** 2 == -1, n)

# AP13 versus all tested terminal primes: labelled equality, not just isomorphism.
ap13 = tuple(range(1, 14))
primes = (17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97)
for p in primes:
    bp = tuple(range(1, 13)) + (p,)
    require(labelled_signature(bp) == labelled_signature(ap13), p)
    require(observer(bp, Fraction(1, 13)) >= Fraction(1, 13), p)

# AP13/V26: the raw labels differ on one edge, but endpoint swap gives an isomorphism.
v26 = tuple(26 * k for k in range(1, 13)) + (339,)
raw_distance = sum(x != y for x, y in zip(labelled_signature(ap13), labelled_signature(v26)))
require(raw_distance == 1, raw_distance)
mapping = {0: 12, 12: 0, **{i: i for i in range(1, 12)}}
for i, j in combinations(range(13), 2):
    require(edge(ap13[i], ap13[j]) == edge(v26[mapping[i]], v26[mapping[j]]), (i, j))

# Independently exact LRC optima and first-target denominator.
ap_max, ap_times = exact_max(ap13)
v26_max, v26_times = exact_max(v26)
pair12_max, pair12_times = exact_max((1, 2))
pair23_max, pair23_times = exact_max((2, 3))
require(ap_max == Fraction(1, 14), (ap_max, ap_times))
require(v26_max == Fraction(1, 13), (v26_max, v26_times))
require(pair12_max == Fraction(1, 3), (pair12_max, pair12_times))
require(pair23_max == Fraction(2, 5), (pair23_max, pair23_times))
first_v26 = first_target_denominator(v26, Fraction(1, 14), 27)
require(first_v26[:2] == (27, 2), first_v26)

# The prime-row statement deliberately gives bounds, not an exact mixed optimum.
odd = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
mixed = (2, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
require(labelled_signature(odd) == labelled_signature(mixed), "prime tournaments")
require(all(labelled_signature(odd)), "odd transitive orientation")
odd_max, _ = exact_max(odd)
mixed_max, mixed_times = exact_max(mixed)
require(odd_max == Fraction(1, 2), odd_max)
require(Fraction(1, 4) <= mixed_max < Fraction(1, 2), (mixed_max, mixed_times))

print("THM4107 HOSTILE REFEREE")
print(f"classified_pairs={pair_count}")
print(f"classified_triples={len(tuple(combinations(range(1, 81), 3)))} cycles={len(cycles)}")
print(f"first_cycles={first_cycles}")
print(f"rational_ratio_candidates={len(rational_ratio_candidates)} all_d=1")
print(f"AP13_labelled_prime_twins={len(primes)}")
print(f"AP13_V26_raw_label_distance={raw_distance} isomorphic_after_endpoint_swap=True")
print(f"M_AP13={ap_max} M_V26={v26_max} first_V26_target_witness={first_v26}")
print(f"M_pair_12={pair12_max} M_pair_23={pair23_max}")
print(f"M_odd={odd_max} M_mixed={mixed_max}")
print("PASS")



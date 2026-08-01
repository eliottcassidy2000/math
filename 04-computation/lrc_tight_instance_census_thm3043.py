"""Exact corrected referee for THM-3043.

Measure-zero safe sets are tight, not empty covers.  This companion performs
the finite census with exact rational interval arithmetic, explicitly checks
that every zero-measure row has a nonempty safe-point witness set, verifies
the six displayed canonical witness sets, and checks the quantisation law.
Truth-bearing checks use explicit raises, so ``python -O`` runs the same audit.
"""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lcm(a, b):
    return a // gcd(a, b) * b


def lcm_many(values):
    return reduce(lcm, values, 1)


def safe_measure(vs, n1):
    """Lebesgue measure of {t in [0,1): ||v t|| >= 1/n1}."""
    bad = []
    for v in vs:
        width = F(1, n1 * v)
        for k in range(v + 1):
            center = F(k, v)
            lo, hi = max(F(0), center - width), min(F(1), center + width)
            if lo < hi:
                bad.append((lo, hi))
    bad.sort()
    merged = []
    for start, end in bad:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return F(1) - sum(end - start for start, end in merged)


def safe_at(vs, n1, t):
    for v in vs:
        residue = (v * t) % 1
        if min(residue, 1 - residue) < F(1, n1):
            return False
    return True


def safe_grid_points(vs, n1):
    """All safe endpoints when the safe measure is zero.

    Every danger endpoint lies on the N=(n+1)lcm(v) grid.  If safe measure is
    zero, every safe point is an endpoint of the finite open-interval union,
    hence lies on this grid.
    """
    N = n1 * lcm_many(vs)
    return tuple(F(k, N) for k in range(N) if safe_at(vs, n1, F(k, N)))


def normalize(vs):
    common = reduce(gcd, vs)
    return tuple(v // common for v in vs)


def is_dilated_ap(vs):
    normalized = normalize(vs)
    return normalized == tuple(range(1, len(vs) + 1))


census_rows = []
for n, bound, expected_zero, expected_ap, expected_sporadic in (
    (3, 14, 4, 4, ()),
    (4, 14, 5, 3, ((1, 3, 4, 7), (2, 6, 8, 14))),
    (5, 12, 3, 2, ((1, 3, 4, 5, 9),)),
):
    zero_rows = []
    for speeds in combinations(range(1, bound + 1), n):
        if safe_measure(speeds, n + 1) == 0:
            witnesses = safe_grid_points(speeds, n + 1)
            require(witnesses, f"zero-measure row was falsely treated as empty: {speeds}")
            zero_rows.append((speeds, witnesses))
    ap_rows = [speeds for speeds, _ in zero_rows if is_dilated_ap(speeds)]
    sporadic = tuple(speeds for speeds, _ in zero_rows if not is_dilated_ap(speeds))
    require(len(zero_rows) == expected_zero, f"zero-measure count changed at n={n}")
    require(len(ap_rows) == expected_ap, f"AP count changed at n={n}")
    require(sporadic == expected_sporadic, f"sporadic rows changed at n={n}")
    census_rows.append((n, bound, len(zero_rows), len(ap_rows), len(sporadic)))


canonical_witnesses = {
    (1, 2, 3): (F(1, 4), F(3, 4)),
    (1, 2, 3, 4): tuple(F(a, 5) for a in range(1, 5)),
    (1, 3, 4, 7): tuple(F(a, 5) for a in range(1, 5)),
    (1, 2, 3, 4, 5, 6): tuple(F(a, 7) for a in range(1, 7)),
    (1, 3, 4, 5, 9): (F(1, 6), F(5, 6)),
    (1, 2, 3, 4, 5, 6, 7): tuple(F(a, 8) for a in (1, 3, 5, 7)),
}
for speeds, expected in canonical_witnesses.items():
    require(safe_measure(speeds, len(speeds) + 1) == 0,
            f"canonical row lost tightness: {speeds}")
    actual = safe_grid_points(speeds, len(speeds) + 1)
    require(actual == expected, f"canonical witness set changed: {speeds}")


quantization_samples = (
    (1, 2, 3),
    (1, 2, 3, 4),
    (1, 3, 4, 5, 7),
    (1, 2, 3, 4, 5, 6),
    (2, 3, 5, 7),
    (1, 4, 5, 6, 7),
    (1, 2, 3, 4, 5, 6, 7),
    (3, 5, 7, 11, 13),
)
positive_ratios = []
for speeds in quantization_samples:
    n1 = len(speeds) + 1
    N = n1 * lcm_many(speeds)
    measure = safe_measure(speeds, n1)
    quantum_count = measure * N
    require(quantum_count.denominator == 1, f"quantisation failed: {speeds}")
    if measure == 0:
        require(safe_grid_points(speeds, n1), f"tight sample became empty: {speeds}")
    else:
        positive_ratios.append(int(quantum_count))
require(tuple(positive_ratios) == (30, 66, 58, 11020), "positive quantum ratios changed")


print("THM3043 LRC TIGHT-INSTANCE CENSUS: CORRECTED PASS")
for n, bound, zero, ap, sporadic in census_rows:
    print(f"n={n} B={bound} zero_measure_tight={zero} dilated_AP={ap} sporadic={sporadic}")
print("sporadic_rows=(1,3,4,7);(2,6,8,14);(1,3,4,5,9)")
print(f"canonical_witness_sets={len(canonical_witnesses)} all_nonempty=1")
print("witness_denominators=4,5,5,7,6,8 respectively")
print(f"quantization_samples={len(quantization_samples)} all_integral=1")
print("positive_mu_over_quantum=30,66,58,11020")
print("zero_measure_label=TIGHT_NONEMPTY not COVERED")
print("all_runtime_checks_explicit=1")

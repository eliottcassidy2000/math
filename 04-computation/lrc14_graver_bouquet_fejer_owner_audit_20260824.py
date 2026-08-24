#!/usr/bin/env python3
"""Exact Graver-bouquet / Fejer / endpoint-owner laboratory for LRC(14).

This is deliberately a bounded proof-design audit, not a proof of LRC(14).
The rows are

    V_w = {1,...,11,13,w},       12 <= w <= 200, w != 13.

They all contain the primitive l1-three Graver relation

    2*(speed 1) - 1*(speed 2) = 0.

For 1 <= k <= 118, the scalar relation has l1 norm 3k <= 354 < 356.
Reading the closed frequency walk 0 -> k -> 2k -> 0 produces the bounded
"Graver bouquet"

    B = union_{k=1}^{118} {0,k,2k}.

For a row V let N_V(x) count its open radius-1/14 danger combs.  If a
trigonometric polynomial P satisfies

    integral (N_V-1)|P|^2 < 0,

then N_V vanishes on a positive-measure open set.  The implication is
pointwise and therefore exact.  All load-bearing signs below are certified
with python-flint/Arb ball arithmetic; floating point is used only to order
already certified diagnostic margins.

The sharpened THM-4009 counterexample reduction gives one Graver row with
square norm at most 195 and l1 norm at most 50.  This script separately tests
the cap-respecting scalar bank of a most favourable l1-three row: only
1 <= k <= 16 keeps 3k <= 50.  On the strict row V_38, every *single ordering*
prefix bouquet of every Euclidean-shortest relation is rigorously positive
definite.  Combining the two orderings creates a negative cross-order
direction for two of the thirty shortest relations, and an explicit integer
coefficient vector supplies a rigorous trigonometric certificate.  Thus the
new cap has genuine finite-LP value, but only when ordering/phase coupling is
retained; relation shortness alone does not close even this toy family.

The AP and Goddyn--Wong rows are mandatory equality hostiles.  Their open
danger combs cover almost everywhere, while t=1/14 is a weak safe endpoint
with two opposite-slope owners.  Consequently no absolutely continuous
Toeplitz/Fejer quadratic can certify them; the endpoint-owner channel is
logically separate.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd

from flint import arb, ctx
import sympy as sp


ctx.prec = 256

LEGACY_CAP = 356
BASE_RELATION_L1 = 3
HARMONIC_MAX = 118
FEJER_K = HARMONIC_MAX + 1
SHORT_CAP = 50
SHORT_HARMONIC_MAX = 16
SHORT_FEJER_K = SHORT_HARMONIC_MAX + 1
CORE = tuple(range(1, 12)) + (13,)
BOUQUET = tuple(
    sorted({0} | {x for k in range(1, HARMONIC_MAX + 1) for x in (k, 2 * k)})
)
DIFF_MULT = Counter(x - y for x in BOUQUET for y in BOUQUET)

# Exact negative vectors discovered from the lowest numerical eigenvectors and
# then frozen as ordinary integers.  The Arb replay, not the discovery float,
# is load bearing.  Coefficients follow the sorted full two-ordering bouquet
# union_{1<=k<=16}{0,k*x,k*y,k*z} for x+y=z.
SHORT_NEGATIVE_COEFFICIENTS = {
    (2, 6, 8): (
        -3, -3, 0, 5, 6, 0, -5, -6, -1, 3, 5, 3, -2, -6, -2, 2, 4,
        -2, 0, 2, -3, 3, 1, -4, 2, 3, -4, 4, 1, -5, 2, 3, -3, 1, 1,
        -2, 2,
    ),
    (2, 4, 6): (
        -88, -79, 15, 115, 147, 25, -121, -155, -24, 86, 138, 69, -57,
        -157, -67, 64, 94, -37, -27, 36, 40, -91, 30, 59, 20, -113, 68,
        77, -85, 108, -112, 89, -62,
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def A(q: Fraction | int) -> arb:
    if isinstance(q, int):
        return arb(q)
    return arb(q.numerator) / q.denominator


def circle_distance(t: Fraction, v: int) -> Fraction:
    r = (t.numerator * v) % t.denominator
    return Fraction(min(r, t.denominator - r), t.denominator)


def best_pair_sum_center(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    """Exact scheduler only; negativity of the later quadratic is the proof.

    Candidate points p/(u+v) are the opposite-slope lower-envelope vertices.
    We do not use completeness of this menu in any negative conclusion.
    """

    best = Fraction(-1)
    center = Fraction(0)
    for i, u in enumerate(speeds):
        for v in speeds[i:]:
            q = u + v
            for p in range(1, q):
                t = Fraction(p, q)
                height = min(circle_distance(t, speed) for speed in speeds)
                if height > best:
                    best = height
                    center = t
    return best, center


def danger_minus_one_hat(speeds: tuple[int, ...], r: int) -> arb:
    """The exact r-th Fourier coefficient of N_V-1."""

    if r == 0:
        return arb(len(speeds)) / 7 - 1
    rr = abs(r)
    value = arb(0)
    for v in speeds:
        if rr % v == 0:
            k = rr // v
            value += (arb.pi() * k / 7).sin() / (arb.pi() * k)
    return value


def additive_triples(speeds: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    """All coefficient-(1,1,-1) relations x+y=z, with x<y."""

    speed_set = set(speeds)
    return tuple((x, y, x + y) for x, y in combinations(speeds, 2) if x + y in speed_set)


def single_order_bouquet(x: int, z: int) -> tuple[int, ...]:
    """Prefixes 0 -> k*x -> k*z -> 0 for 1<=k<=16."""

    return tuple(
        sorted(
            {0}
            | {
                frequency
                for k in range(1, SHORT_HARMONIC_MAX + 1)
                for frequency in (k * x, k * z)
            }
        )
    )


def two_order_bouquet(x: int, y: int, z: int) -> tuple[int, ...]:
    """Union of both non-equivalent prefix orderings for x+y=z."""

    return tuple(
        sorted(
            {0}
            | {
                frequency
                for k in range(1, SHORT_HARMONIC_MAX + 1)
                for frequency in (k * x, k * y, k * z)
            }
        )
    )


def moment_quadratic(
    speeds: tuple[int, ...], bouquet: tuple[int, ...], coefficients: tuple[int, ...]
) -> arb:
    """Exact Hermitian danger quadratic for a real integer coefficient word."""

    require(len(bouquet) == len(coefficients), "coefficient/bouquet mismatch")
    total = arb(0)
    for i, lam in enumerate(bouquet):
        for j, mu in enumerate(bouquet):
            total += (
                coefficients[i]
                * coefficients[j]
                * danger_minus_one_hat(speeds, lam - mu)
            )
    return total


def certify_positive_definite(
    speeds: tuple[int, ...], bouquet: tuple[int, ...], pivot_floor: Fraction
) -> arb:
    """Rigorous interval LDL^T certificate, returning its smallest pivot.

    No floating-point sign decision occurs.  A midpoint float is used only to
    select which already-certified positive Arb pivot to print.
    """

    size = len(bouquet)
    lower: list[list[arb]] = [
        [arb(0) for _ in range(size)] for _ in range(size)
    ]
    diagonal: list[arb] = []
    smallest: arb | None = None

    for i in range(size):
        lower[i][i] = arb(1)
        pivot = danger_minus_one_hat(speeds, 0)
        for k in range(i):
            pivot -= lower[i][k] * lower[i][k] * diagonal[k]
        require(
            pivot > A(pivot_floor),
            f"LDL pivot {i} did not clear {pivot_floor}: {pivot}",
        )
        diagonal.append(pivot)
        if smallest is None or float(pivot.mid()) < float(smallest.mid()):
            smallest = pivot

        for j in range(i + 1, size):
            entry = danger_minus_one_hat(speeds, bouquet[j] - bouquet[i])
            for k in range(i):
                entry -= lower[j][k] * lower[i][k] * diagonal[k]
            lower[j][i] = entry / pivot

    require(smallest is not None, "empty bouquet has no LDL pivot")
    return smallest


def complex_bouquet_quadratic(speeds: tuple[int, ...], center: Fraction) -> arb:
    """P has coefficient exp(-2*pi*i*lambda*center) on every lambda in B."""

    total = arb(0)
    for d, multiplicity in DIFF_MULT.items():
        phase = (2 * arb.pi() * d * center.numerator / center.denominator).cos()
        total += multiplicity * danger_minus_one_hat(speeds, d) * phase
    return total / len(BOUQUET)


def oriented_bouquet_quadratic(speeds: tuple[int, ...], center: Fraction) -> arb:
    """Reflection-breaking real probe c_lambda=cos(theta)-sin(theta)."""

    coeffs: list[arb] = []
    for lam in BOUQUET:
        theta = 2 * arb.pi() * lam * center.numerator / center.denominator
        coeffs.append(theta.cos() - theta.sin())
    norm = sum((c * c for c in coeffs), arb(0))
    total = arb(0)
    for i, lam in enumerate(BOUQUET):
        for j, mu in enumerate(BOUQUET):
            total += coeffs[i] * coeffs[j] * danger_minus_one_hat(speeds, lam - mu)
    return total / norm


def centered_fejer_quadratic(
    speeds: tuple[int, ...], center: Fraction, order: int
) -> arb:
    """Integral (N_V-1) F_order(x-center), exactly in Arb."""

    total = danger_minus_one_hat(speeds, 0)
    for r in range(1, order):
        weight = arb(order - r) / order
        phase = (2 * arb.pi() * r * center.numerator / center.denominator).cos()
        total += 2 * weight * danger_minus_one_hat(speeds, r) * phase
    return total


def closed_danger_union_length(speeds: tuple[int, ...]) -> Fraction:
    """Exact measure of the union; endpoint convention is measure-immaterial."""

    intervals: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        radius = Fraction(1, 14 * v)
        for j in range(v + 1):
            center = Fraction(j, v)
            left = max(Fraction(0), center - radius)
            right = min(Fraction(1), center + radius)
            if left <= right:
                intervals.append((left, right))
    intervals.sort()
    merged: list[list[Fraction]] = []
    for left, right in intervals:
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    return sum((right - left for left, right in merged), Fraction(0))


def boundary_owners(speeds: tuple[int, ...], t: Fraction) -> tuple[tuple[int, int], ...]:
    owners: list[tuple[int, int]] = []
    for v in speeds:
        if circle_distance(t, v) != Fraction(1, 14):
            continue
        residue = (t.numerator * v) % t.denominator
        # At residue +1 the distance rises with t; at residue -1 it falls.
        sign = 1 if residue * 2 < t.denominator else -1
        owners.append((v, sign))
    return tuple(owners)


def symbolic_s38_certificate() -> tuple[str, int, Fraction, Fraction, Fraction]:
    """Human-readable exact algebraic certificate for the K=119 sign."""

    speeds = CORE + (38,)
    S = sp.S(0)
    for v in speeds:
        for k in range(1, (FEJER_K - 1) // v + 1):
            S += (
                2
                * (1 - sp.Rational(k * v, FEJER_K))
                * sp.sin(sp.pi * k / 7)
                * sp.cos(sp.pi * k * v / 6)
                / k
            )
    S = sp.simplify(S)
    x = sp.symbols("x")
    polynomial = sp.Poly(sp.minpoly(S, x), x)
    intervals = polynomial.intervals(eps=sp.Rational(1, 10**12))

    s_ball = arb(0)
    for v in speeds:
        for k in range(1, (FEJER_K - 1) // v + 1):
            s_ball += (
                2
                * (arb(1) - arb(k * v) / FEJER_K)
                * (arb.pi() * k / 7).sin()
                * (arb.pi() * k * v / 6).cos()
                / k
            )

    chosen: tuple[Fraction, Fraction] | None = None
    for (left, right), multiplicity in intervals:
        lo = Fraction(int(left.p), int(left.q))
        hi = Fraction(int(right.p), int(right.q))
        if s_ball > A(lo) and s_ball < A(hi):
            require(multiplicity == 1, "S38 algebraic root was not simple")
            chosen = (lo, hi)
            break
    require(chosen is not None, "failed to isolate the exact S38 algebraic sum")
    lo, hi = chosen

    # Q=6/7+S/pi.  Since pi<355/113, S<hi<0 gives
    # S+6*pi/7 < hi+6*(355/113)/7.  The latter is rational and negative.
    pi_upper = Fraction(355, 113)
    rational_margin = hi + Fraction(6, 7) * pi_upper
    require(rational_margin < 0, "rational S38 sign margin was not negative")
    require(arb.pi() < A(pi_upper), "Arb did not verify pi<355/113")

    return str(S), polynomial.degree(), lo, hi, rational_margin


def main() -> None:
    require(BASE_RELATION_L1 * HARMONIC_MAX == 354, "relation budget changed")
    require(354 <= LEGACY_CAP < 3 * 119, "legacy harmonic cutoff changed")
    require(
        BASE_RELATION_L1 * SHORT_HARMONIC_MAX == 48
        and 48 <= SHORT_CAP < 3 * 17,
        "short harmonic cutoff changed",
    )
    require(len(BOUQUET) == 178 and max(BOUQUET) == 236, "bouquet census changed")
    require(set(range(FEJER_K)) <= set(BOUQUET), "Fejer bank escaped bouquet")

    pair_atlas = tuple(
        (p, q)
        for p in range(1, 14)
        for q in range(p + 1, 14)
        if gcd(p, q) == 1 and p * p + q * q <= 195
    )
    require(len(pair_atlas) == 47, "THM-4009 support-two atlas count changed")
    require(max(p + q for p, q in pair_atlas) == 19, "pair l1 maximum changed")
    unoriented_pair_packets = 78 * len(pair_atlas)
    oriented_pair_assignments = 2 * unoriented_pair_packets
    require(unoriented_pair_packets == 3666, "unoriented pair packet count changed")
    require(oriented_pair_assignments == 7332, "oriented pair count changed")

    ap = tuple(range(1, 14))
    gw = CORE + (24,)
    for name, row in (("AP", ap), ("GW", gw)):
        require(closed_danger_union_length(row) == 1, f"{name} did not cover a.e.")
        require(
            min(circle_distance(Fraction(1, 14), v) for v in row) == Fraction(1, 14),
            f"{name} lost its weak boundary witness",
        )
        require(boundary_owners(row, Fraction(1, 14)) == ((1, 1), (13, -1)),
                f"{name} owner regression")

    s38 = CORE + (38,)
    shortest_triples = additive_triples(s38)
    require(len(shortest_triples) == 30, "V_38 additive-triple census changed")

    # Distinct speeds rule out norm squares 1 and 2.  Every displayed triple
    # has norm square 3, so these are exactly all Euclidean-shortest relations.
    single_order_minimum: arb | None = None
    single_order_count = 0
    for x, y, z in shortest_triples:
        for first in (x, y):
            pivot = certify_positive_definite(
                s38, single_order_bouquet(first, z), Fraction(2, 5)
            )
            single_order_count += 1
            if (
                single_order_minimum is None
                or float(pivot.mid()) < float(single_order_minimum.mid())
            ):
                single_order_minimum = pivot
    require(single_order_count == 60, "single-order bouquet count changed")
    require(single_order_minimum is not None, "missing single-order pivot")

    short_negative_values: dict[tuple[int, int, int], arb] = {}
    full_order_positive_count = 0
    full_order_minimum: arb | None = None
    for triple in shortest_triples:
        x, y, z = triple
        bouquet = two_order_bouquet(x, y, z)
        if triple in SHORT_NEGATIVE_COEFFICIENTS:
            coefficients = SHORT_NEGATIVE_COEFFICIENTS[triple]
            value = moment_quadratic(s38, bouquet, coefficients)
            require(value < 0, f"short exact negative certificate failed at {triple}")
            short_negative_values[triple] = value
        else:
            pivot = certify_positive_definite(s38, bouquet, Fraction(1, 3))
            full_order_positive_count += 1
            if (
                full_order_minimum is None
                or float(pivot.mid()) < float(full_order_minimum.mid())
            ):
                full_order_minimum = pivot
    require(len(short_negative_values) == 2, "short negative triple count changed")
    require(full_order_positive_count == 28, "short positive triple count changed")
    require(full_order_minimum is not None, "missing full-order positive pivot")

    # The complete two-order bank stays positive on the two equality controls.
    equality_full_order_minimum: dict[str, arb] = {}
    for name, row in (("AP", ap), ("GW", gw)):
        minimum: arb | None = None
        for x, y, z in shortest_triples:
            pivot = certify_positive_definite(
                row, two_order_bouquet(x, y, z), Fraction(2, 5)
            )
            if minimum is None or float(pivot.mid()) < float(minimum.mid()):
                minimum = pivot
        require(minimum is not None, f"missing {name} full-order pivot")
        equality_full_order_minimum[name] = minimum

    short_fejer_values = [
        centered_fejer_quadratic(s38, Fraction(1, 12), order)
        for order in range(1, SHORT_FEJER_K + 1)
    ]
    require(
        all(value > 0 for value in short_fejer_values),
        "a cap-50 V_38 centered Fejer order became negative",
    )

    fejer_values = [
        centered_fejer_quadratic(s38, Fraction(1, 12), order)
        for order in range(1, FEJER_K + 1)
    ]
    require(all(value > 0 for value in fejer_values[:-1]),
            "an earlier S38 centered Fejer order became negative")
    require(fejer_values[-1] < 0, "the order-119 S38 Fejer certificate failed")

    expression, degree, iso_lo, iso_hi, rational_margin = symbolic_s38_certificate()

    complex_negative = 0
    oriented_only = 0
    boundary_controls = 0
    uncertified: list[int] = []
    closest_complex: tuple[float, int, arb] | None = None
    special_rows: dict[int, tuple[Fraction, Fraction, arb, arb | None]] = {}

    for w in range(12, 201):
        if w in CORE:
            continue
        speeds = CORE + (w,)
        candidate_height, center = best_pair_sum_center(speeds)
        complex_q = complex_bouquet_quadratic(speeds, center)
        oriented_q: arb | None = None

        if w in (12, 24):
            boundary_controls += 1
            require(complex_q > 0, "boundary control acquired a false negative")
        elif complex_q < 0:
            complex_negative += 1
            diagnostic = float(complex_q.mid())
            if closest_complex is None or diagnostic > closest_complex[0]:
                closest_complex = (diagnostic, w, complex_q)
        else:
            oriented_q = oriented_bouquet_quadratic(speeds, center)
            if oriented_q < 0:
                oriented_only += 1
            else:
                uncertified.append(w)

        if w in (12, 24, 36, 38):
            special_rows[w] = (candidate_height, center, complex_q, oriented_q)

    require(complex_negative == 185, "complex bouquet negative count changed")
    require(oriented_only == 1, "oriented-only count changed")
    require(boundary_controls == 2, "boundary control count changed")
    require(not uncertified, f"uncertified rows remain: {uncertified}")
    require(closest_complex is not None and closest_complex[1] == 47,
            "closest complex margin row changed")
    require(special_rows[36][3] is not None and special_rows[36][3] < 0,
            "w=36 orientation certificate failed")

    print("LRC14 GRAVER-BOUQUET / FEJER / OWNER AUDIT")
    print("status=FINITE-EXACT bounded laboratory; LRC(14) remains OPEN")
    print("universe=V_w={1,...,11,13,w}, 12<=w<=200, w!=13 (188 rows)")
    print("")
    print("THM-4009 sharpened-cap audit")
    print("counterexample_row_bound=square_norm<=195; l1<=50; height<=13")
    print(
        "support_two_ratios=47; unoriented_ratio_support_packets=3666; "
        "oriented_labelled_assignments=7332; pair_l1_max=19"
    )
    print("cap_respecting_l1_three_harmonics=1..16; largest_l1=48")
    print("V_38 shortest_relation_norm_square=3")
    print("V_38 shortest_relations=30 additive triples x+y=z")
    print(f"single_order_bouquets={single_order_count}; positive_definite=60")
    print(f"single_order_all_LDL_pivots_gt=2/5; minimum={single_order_minimum}")
    print("two_order_bouquets=30; positive_definite=28; exact_negative=2")
    print(f"two_order_positive_all_LDL_pivots_gt=1/3; minimum={full_order_minimum}")
    for triple in sorted(short_negative_values):
        coefficients = SHORT_NEGATIVE_COEFFICIENTS[triple]
        norm_square = sum(coefficient * coefficient for coefficient in coefficients)
        value = short_negative_values[triple]
        print(
            f"two_order_negative={triple}; coefficient_norm_square={norm_square}; "
            f"quadratic={value}; normalized={value / norm_square}"
        )
    print(
        "V_38 short_Fejer_orders_1_through_17_positive=yes; "
        f"order_17={short_fejer_values[-1]}"
    )
    print(
        "AP all_30_two_order_bouquets_positive_definite=yes; "
        f"all_LDL_pivots_gt=2/5; minimum={equality_full_order_minimum['AP']}"
    )
    print(
        "GW all_30_two_order_bouquets_positive_definite=yes; "
        f"all_LDL_pivots_gt=2/5; minimum={equality_full_order_minimum['GW']}"
    )
    print("verdict=shortness alone fails; cross-order phase creates the strict certificate")
    print("endpoint_owner_still_absent=yes")
    print("")
    print("legacy THM-3743 cap-respecting harmonic laboratory")
    print("base_relation=2*(speed 1)-1*(speed 2)=0; primitive_l1=3; Graver=yes")
    print("flatness_cap=356; harmonic_range=1..118; largest_l1=354")
    print(f"bouquet_size={len(BOUQUET)}; max_frequency={max(BOUQUET)}")
    print("bouquet=union_{1<=k<=118}{0,k,2k}")
    print("")
    print("strict spectral census")
    print(f"complex_bouquet_negative={complex_negative}")
    print(f"reflection_breaking_oriented_only={oriented_only} (w=36)")
    print(f"boundary_controls={boundary_controls} (w=12 AP, w=24 GW)")
    print(f"uncertified={len(uncertified)}")
    print(f"closest_certified_complex=w={closest_complex[1]} q={closest_complex[2]}")
    print("")
    for w in (12, 24, 36, 38):
        height, center, complex_q, oriented_q = special_rows[w]
        oriented_text = "not_needed" if oriented_q is None else str(oriented_q)
        print(
            f"w={w} candidate_height={height} center={center} "
            f"complex_q={complex_q} oriented_q={oriented_text}"
        )
    print("")
    print("boundary owner channel")
    print("AP danger_union_measure=1; weak_time=1/14; owners=((1,+),(13,-))")
    print("GW danger_union_measure=1; weak_time=1/14; owners=((1,+),(13,-))")
    print("therefore every absolutely-continuous Toeplitz quadratic is nonnegative")
    print("and the weak equality witness requires the atomic endpoint-owner sidecar")
    print("")
    print("sharp centered Fejer control on V_38")
    print(f"orders_1_through_118_positive=yes; order_118={fejer_values[117]}")
    print(f"order_119={fejer_values[118]} (negative)")
    print(f"algebraic_sum={expression}")
    print(f"algebraic_degree={degree}")
    print(f"isolating_interval=[{iso_lo},{iso_hi}]")
    print("pi_upper=355/113")
    print(f"exact_upper_margin_for_S_plus_6pi_over_7={rational_margin}<0")
    print("")
    print("connection_contract")
    print("source=THM-4009 short Graver row plus cap-respecting scalar relations")
    print("map=closed-walk prefixes -> frequency bouquet -> Hermitian danger matrix")
    print("preserved=sign partition, harmonic scale, prefix ordering, cross-order phase")
    print("destroyed=measure-zero equality, endpoint owner, physical component incidence")
    print("sidecar=exact signed endpoint owners / cross-phase word")
    print("cheapest_next_test=exact SDP on 47 pair types x 17 THM-3910 hostiles")
    print("PASS")


if __name__ == "__main__":
    main()

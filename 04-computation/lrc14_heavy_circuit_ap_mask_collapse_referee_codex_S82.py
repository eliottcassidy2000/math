#!/usr/bin/env python3
"""Exact referee for THM-1218: heavy circuits collapse the beat-mask nerve.

For four ordered integer frequencies ``x0<x1<x2<x3``, let BAD be the
continuum max-gap event from THM-1203.  Put

    H = 60/637.

THM-1203 gives two ingredients.  First, BAD is contained in each of the
four additive-triangle events belonging to the deletions of one point.
Second, after reducing unequal triangle gaps ``p,q`` by their gcd,

    p+q >= 26  => J(p,q) <= H,

while the exact 99-pair core has only the ratio ``1:2`` strictly above H.
Thus ``mu(BAD)>H`` puts all four deletion-gap pairs in the ratio set
``{1,1/2,2}``; the elementary four-obligation system forces all three
successive gaps to agree.

Applied to ``(bi,bj,b5,b6)`` with ``q=b6-b5``, the only possible heavy
quadruple is

    (b5-2q, b5-q, b5, b5+q).

The four speeds are congruent modulo q.  Their direct q-danger masks, and
therefore their restrictions/lifts to the actual master clock, coincide.  The
actual mixed-period clock is always a divisor ``L=q/d0`` of ``q``.  Longer
arbitrary numerator ranges also preserve the pointwise equality, but those
tests are reported separately as restriction replays rather than being
mislabelled as master-clock constructions.  Among
the usual five labels (the defining pair already contributes one label), a
heavy pair leaves at most three distinct masks.  In the common-Q unit case
this changes THM-1216's threshold from ``B5=5A-3`` to

    B3 = 2+3(A-1) = 6*ceil(Q/14)-4.

All computations below use integers and ``fractions.Fraction``.  Checks use
``require`` rather than ``assert`` so optimized Python executes the same
certificate.

Tournament/assumption audit
---------------------------
The six candidate heavy-pair vertices form a lexicographically gauged
transitive tournament, but the orientation contains no proof information.
The faithful object is the four-deletion additive-circuit hypergraph, then
the equality classes in the five-mask nerve.  This challenges the default
choice of runners as vertices: pair candidates, deletion circuits, residue
masks, and master-clock residues preserve the predicates used here; a bare
runner tournament destroys the continuum mass and mask equality.
"""

from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import gcd


CUTOFF = F(60, 637)
TOP = F(2, 21)
CORE_SUM_LIMIT = 25
BAD_GAP_LIMIT = F(2, 7)
THREE_BANDS = ((1, 2), (3, 4), (5, 6))
QUARTET_DIAMETER = 24
SIX_SPEED_HEIGHT = 16


def require(condition: bool, message: object) -> None:
    """Always-on certificate check, including under ``python -O``."""
    if not condition:
        raise RuntimeError(message)


def ceil_div(numerator: int, denominator: int) -> int:
    require(denominator > 0, denominator)
    return -((-numerator) // denominator)


def fractional_part(value: F) -> F:
    return value - value.numerator // value.denominator


def in_three_band(value: F) -> bool:
    value = fractional_part(value)
    return any(F(left, 7) < value < F(right, 7) for left, right in THREE_BANDS)


@lru_cache(maxsize=None)
def three_band_intervals(multiplier: int) -> tuple[tuple[F, F], ...]:
    require(multiplier > 0, multiplier)
    denominator = 7 * multiplier
    return tuple(
        (F(7 * k + left, denominator), F(7 * k + right, denominator))
        for k in range(multiplier)
        for left, right in THREE_BANDS
    )


def intersect_interval_lists(
    first: tuple[tuple[F, F], ...], second: tuple[tuple[F, F], ...]
) -> tuple[tuple[F, F], ...]:
    """Exact intersection of sorted internally disjoint open intervals."""
    i = j = 0
    answer: list[tuple[F, F]] = []
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            answer.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        elif second[j][1] < first[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return tuple(answer)


def interval_measure(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in intervals), F(0))


@lru_cache(maxsize=None)
def triangle_measure(p: int, q: int) -> F:
    """The exact additive three-band mass J(p,q)."""
    require(p > 0 and q > 0, (p, q))
    common = gcd(p, q)
    p //= common
    q //= common
    intersection = intersect_interval_lists(
        three_band_intervals(p), three_band_intervals(q)
    )
    intersection = intersect_interval_lists(
        intersection, three_band_intervals(p + q)
    )
    return interval_measure(intersection)


def triangle_carry_formula(p: int, q: int) -> F:
    """Independent THM-1203 carry formula for a reduced ordered pair."""
    require(0 < p < q and gcd(p, q) == 1, (p, q))
    total = p + q
    residues = (
        (2 * p - q) % 7,
        (4 * p - q) % 7,
        (2 * p - 3 * q) % 7,
    )
    weighted = F(0)
    for k in range(1, total):
        tent = F(k, q) if k <= q else F(total - k, p)
        multiplicity = sum(k % 7 == residue for residue in residues)
        weighted += multiplicity * tent
    return F(2, 7 * total) * weighted


def ratio12(p: int, q: int) -> bool:
    return p == q or p == 2 * q or q == 2 * p


def strict_triangle_spectrum() -> dict[str, object]:
    """Certify that J>60/637 forces the ratio 1:2 (or 2:1)."""
    rows: list[tuple[F, int, int]] = []
    for p in range(1, CORE_SUM_LIMIT):
        for q in range(p + 1, CORE_SUM_LIMIT + 1 - p):
            if gcd(p, q) != 1:
                continue
            value = triangle_measure(p, q)
            require(value == triangle_carry_formula(p, q), (p, q, value))
            rows.append((value, p, q))

    rows.sort(reverse=True)
    heavy = tuple((p, q) for value, p, q in rows if value > CUTOFF)
    require(len(rows) == 99, len(rows))
    require(heavy == ((1, 2),), heavy)
    require(rows[0] == (F(2, 21), 1, 2), rows[:3])
    require(rows[1] == (F(13, 147), 1, 6), rows[:3])
    require(rows[2] == (F(8, 105), 3, 10), rows[:3])

    tail_endpoint = F(3, 49) + F(6, 7 * 26)
    require(tail_endpoint == CUTOFF, tail_endpoint)
    require(TOP - CUTOFF == F(2, 1911), TOP - CUTOFF)
    require(CUTOFF - F(13, 147) == F(11, 1911), CUTOFF - F(13, 147))

    # Guard the weak/strict boundary: N=26 is at the cutoff and is therefore
    # excluded by the theorem's strict hypothesis, while N>26 lies below it.
    for total in range(26, 10_001):
        tail = F(3, 49) + F(6, 7 * total)
        require(tail <= CUTOFF, (total, tail))
        require((tail == CUTOFF) == (total == 26), (total, tail))

    # Common dilation must preserve the event and the heavy classification.
    for p, q, scale in product(range(1, 20), range(1, 20), range(1, 8)):
        require(triangle_measure(scale * p, scale * q) == triangle_measure(p, q),
                (p, q, scale))
        if triangle_measure(p, q) > CUTOFF:
            require(ratio12(p, q), (p, q))

    return {
        "rows": len(rows),
        "heavy": heavy,
        "top": rows[:3],
        "tail_endpoint": tail_endpoint,
        "top_margin": TOP - CUTOFF,
        "second_margin": CUTOFF - F(13, 147),
    }


def four_obligation_rigidity(limit: int = 64) -> dict[str, int]:
    """Exact finite audit of the unbounded elementary omega lemma."""
    rows = survivors = 0
    for a, b, c in product(range(1, limit + 1), repeat=3):
        obligations = ((a, b), (b, c), (a, b + c), (a + b, c))
        if all(ratio12(p, q) for p, q in obligations):
            require(a == b == c, (a, b, c, obligations))
            survivors += 1
        rows += 1
    require(survivors == limit, survivors)
    return {"rows": rows, "survivors": survivors}


def max_cyclic_gap(points: tuple[F, ...]) -> F:
    ordered = sorted(points)
    gaps = [right - left for left, right in zip(ordered, ordered[1:])]
    gaps.append(1 + ordered[0] - ordered[-1])
    return max(gaps)


def arrangement_walls(differences: set[int]) -> tuple[F, ...]:
    walls = {F(0), F(1)}
    for difference in differences:
        d = abs(difference)
        require(d > 0, differences)
        for k in range(d + 1):
            walls.add(F(k, d))
        for residue in (2, 5):
            for k in range(d):
                walls.add(F(7 * k + residue, 7 * d))
    return tuple(sorted(walls))


@lru_cache(maxsize=None)
def bad_measure(spacings: tuple[int, int, int]) -> F:
    """Exact continuum BAD mass for {0,d2,d3,d4}."""
    d2, d3, d4 = spacings
    require(0 < d2 < d3 < d4, spacings)
    values = (0, d2, d3, d4)
    differences = {right - left for left, right in combinations(values, 2)}
    walls = arrangement_walls(differences)
    answer = F(0)
    for left, right in zip(walls, walls[1:]):
        u = (left + right) / 2
        points = tuple(fractional_part(value * u) for value in values)
        if max_cyclic_gap(points) <= BAD_GAP_LIMIT:
            answer += right - left
    return answer


def quartet_census() -> dict[str, object]:
    """Independently test actual BAD mass on all gap triples of bounded diameter."""
    rows = heavy = 0
    non_ap_max = F(0)
    for a in range(1, QUARTET_DIAMETER + 1):
        for b in range(1, QUARTET_DIAMETER + 1 - a):
            for c in range(1, QUARTET_DIAMETER + 1 - a - b):
                measure = bad_measure((a, a + b, a + b + c))
                is_heavy = measure > CUTOFF
                require(is_heavy == (a == b == c), (a, b, c, measure))
                if is_heavy:
                    require(measure == TOP, (a, b, c, measure))
                    heavy += 1
                else:
                    non_ap_max = max(non_ap_max, measure)
                rows += 1
    require(rows == 2024, rows)
    require(heavy == QUARTET_DIAMETER // 3, heavy)
    require(non_ap_max < CUTOFF, non_ap_max)
    return {"rows": rows, "heavy": heavy, "non_ap_max": non_ap_max}


def heavy_pairs(speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    require(len(speeds) == 6 and tuple(sorted(speeds)) == speeds, speeds)
    b5, b6 = speeds[4], speeds[5]
    answer = []
    for i, j in combinations(range(4), 2):
        bi, bj = speeds[i], speeds[j]
        translated = (bj - bi, b5 - bi, b6 - bi)
        if bad_measure(translated) > CUTOFF:
            answer.append((i, j))
    return tuple(answer)


def six_speed_unique_pair_census() -> dict[str, int]:
    """Audit uniqueness among the six candidate pairs in ordered sextuples."""
    rows = rows_with_heavy = 0
    for speeds in combinations(range(1, SIX_SPEED_HEIGHT + 1), 6):
        pairs = heavy_pairs(speeds)
        require(len(pairs) <= 1, (speeds, pairs))
        q = speeds[5] - speeds[4]
        target = (speeds[4] - 2 * q, speeds[4] - q)
        expected = tuple(
            (i, j) for i, j in combinations(range(4), 2)
            if (speeds[i], speeds[j]) == target
        )
        require(pairs == expected, (speeds, pairs, expected))
        rows_with_heavy += bool(pairs)
        rows += 1
    return {"rows": rows, "heavy_rows": rows_with_heavy}


def danger_mask(q: int, speed: int) -> frozenset[int]:
    """Direct strict q-danger mask on beat numerators modulo q."""
    require(q > 0 and speed > 0, (q, speed))
    return frozenset(
        p for p in range(q)
        if 14 * min((speed * p) % q, (-speed * p) % q) < q
    )


def clock_restriction(L: int, q: int, speed: int) -> frozenset[int]:
    """Restrict the direct q-predicate to numerator representatives ``0..L-1``."""
    require(L > 0 and q > 0, (L, q))
    return frozenset(
        p for p in range(L)
        if 14 * min((speed * p) % q, (-speed * p) % q) < q
    )


def divisors(n: int) -> tuple[int, ...]:
    require(n > 0, n)
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def gcd_many(values: tuple[int, ...]) -> int:
    answer = 0
    for value in values:
        answer = gcd(answer, value)
    return answer


def reduced_mask(Q: int, unit: int) -> frozenset[int]:
    """The strict unit mask on the reduced quotient clock Z/Q."""
    require(Q > 0 and gcd(Q, unit) == 1, (Q, unit))
    return frozenset(
        p for p in range(Q)
        if 14 * min((unit * p) % Q, (-unit * p) % Q) < Q
    )


def quotient_lift(L: int, Q: int, unit: int) -> frozenset[int]:
    """Lift a reduced mask along Z/L -> Z/Q; this requires Q|L."""
    require(L > 0 and Q > 0 and L % Q == 0, (L, Q))
    base = reduced_mask(Q, unit)
    return frozenset(p for p in range(L) if p % Q in base)


def window_count(Q: int) -> int:
    return 2 * ceil_div(Q, 14) - 1


def class_threshold(Q: int, classes: int) -> int:
    return 2 + classes * (window_count(Q) - 1)


def mask_collapse_audit() -> dict[str, int]:
    """Separate genuine divisor master clocks from arbitrary-range replays."""
    actual_master_rows = actual_row_masks = direct_rows = 0
    supplementary_replay_rows = 0

    # Every divisor L|q is realized as the *actual* master clock q/d0 by the
    # ordered packet below.  Its last four entries form the AP completion and
    # its last pair has difference q.  The first speed d forces d0=d=q/L.
    for q in range(1, 101):
        for L in divisors(q):
            d = q // L
            speeds = (
                d,
                2 * d,
                d * (8 * L + 1),
                d * (9 * L + 1),
                d * (10 * L + 1),
                d * (11 * L + 1),
            )
            require(tuple(sorted(speeds)) == speeds and len(set(speeds)) == 6,
                    (q, L, speeds))
            require(speeds[5] - speeds[4] == q, (q, L, speeds))
            require(speeds[2] + 2 * q == speeds[4], (q, L, speeds))
            require(speeds[3] + q == speeds[4], (q, L, speeds))

            d0 = gcd_many((q, *speeds))
            require(d0 == d and q // d0 == L, (q, L, d0, speeds))

            lifted: list[frozenset[int]] = []
            for speed in speeds:
                gi = gcd(q, speed)
                Qi = q // gi
                unit = speed // gi
                require(L % Qi == 0 and gcd(unit, Qi) == 1,
                        (q, L, speed, Qi, unit))
                direct_on_master = clock_restriction(L, q, speed)
                quotient = quotient_lift(L, Qi, unit)
                require(direct_on_master == quotient,
                        (q, L, speed, Qi, direct_on_master, quotient))
                lifted.append(quotient)
                actual_row_masks += 1

            quartet_direct = tuple(danger_mask(q, speed) for speed in speeds[2:])
            require(len(set(quartet_direct)) == 1, (q, L, speeds))
            require(len(set(lifted[2:])) == 1, (q, L, speeds, lifted[2:]))

            # Five relevant labels are b1,...,b5.  The heavy AP pair b3,b4
            # and defining label b5 are one class, leaving at most two others.
            require(len(set(lifted[:5])) <= 3, (q, L, speeds, lifted[:5]))
            actual_master_rows += 1
            direct_rows += 1

    # Supplementary guardrail only: equality is pointwise, so it persists on
    # arbitrary longer ranges.  These R do not divide q and are not described
    # as mixed-period master clocks.
    for q in range(1, 101):
        quartet = (8 * q + 1, 9 * q + 1, 10 * q + 1, 11 * q + 1)
        require(len({danger_mask(q, speed) for speed in quartet}) == 1,
                (q, quartet))
        for R in (2 * q, 3 * q, 5 * q):
            require(q % R != 0, (q, R))
            replay = tuple(clock_restriction(R, q, speed) for speed in quartet)
            require(len(set(replay)) == 1, (q, R, quartet))
            supplementary_replay_rows += 1

    for Q in range(2, 10_001):
        h = ceil_div(Q, 14)
        A = window_count(Q)
        B3 = class_threshold(Q, 3)
        B5 = class_threshold(Q, 5)
        require(B3 == 6 * h - 4, (Q, B3, h))
        require(B3 == 3 * A - 1, (Q, B3, A))
        require(2 <= B3 <= B5 <= Q, (Q, B3, B5))
        require(B5 - B3 == 4 * (h - 1), (Q, B3, B5))

    # Exact common-zero union ledger for every choice of at most three unit
    # masks in small periods.  This is diagnostic; the proof is cardinality.
    union_rows = 0
    for Q in range(2, 41):
        masks = sorted(
            {danger_mask(Q, unit) for unit in range(1, Q) if gcd(unit, Q) == 1},
            key=lambda mask: tuple(sorted(mask)),
        )
        A = window_count(Q)
        for r in range(1, min(3, len(masks)) + 1):
            for selected in combinations(masks, r):
                union = set().union(*selected)
                require(all(0 in mask and len(mask) == A for mask in selected),
                        (Q, selected))
                require(len(union) <= 1 + r * (A - 1), (Q, r, union))
                union_rows += 1
    return {
        "actual_master_rows": actual_master_rows,
        "actual_row_masks": actual_row_masks,
        "direct_rows": direct_rows,
        "supplementary_replay_rows": supplementary_replay_rows,
        "union_rows": union_rows,
    }


def main() -> None:
    spectrum = strict_triangle_spectrum()
    rigidity = four_obligation_rigidity()
    quartets = quartet_census()
    sextuples = six_speed_unique_pair_census()
    masks = mask_collapse_audit()

    print("THM-1218 HEAVY-CIRCUIT / AP-MASK-COLLAPSE EXACT REFEREE")
    print("arithmetic=fractions.Fraction + integers; all checks survive python -O")
    print(f"strict cutoff H={CUTOFF}; AP mass={TOP}; AP margin={spectrum['top_margin']}")
    print(
        f"tail guard: 3/49+6/[7N] <= H for N>=26; equality only N=26; "
        f"endpoint={spectrum['tail_endpoint']}"
    )
    print(
        f"exact reduced core: rows={spectrum['rows']}; J>H pairs={spectrum['heavy']}; "
        f"second-stratum margin={spectrum['second_margin']}"
    )
    print("top reduced strata=" + "; ".join(
        f"{value}@({p},{q})" for value, p, q in spectrum["top"]
    ))
    print(
        f"four deletion obligations: rows={rigidity['rows']}; "
        f"survivors={rigidity['survivors']} (exactly a=b=c)"
    )
    print(
        f"actual BAD census, positive gaps with diameter<={QUARTET_DIAMETER}: "
        f"rows={quartets['rows']}; heavy={quartets['heavy']}; "
        f"largest non-AP mass={quartets['non_ap_max']}"
    )
    print(
        f"ordered sextuple census through {SIX_SPEED_HEIGHT}: rows={sextuples['rows']}; "
        f"rows with a heavy pair={sextuples['heavy_rows']}; max heavy-pair count=1"
    )
    print(
        "actual divisor-clock realizations L=q/d0|q: "
        f"rows={masks['actual_master_rows']}; "
        f"row-mask identities={masks['actual_row_masks']}; PASS"
    )
    print(
        f"direct-q AP coincidences on those rows={masks['direct_rows']}; PASS"
    )
    print(
        "supplementary arbitrary-range restrictions (not master clocks): "
        f"rows={masks['supplementary_replay_rows']}; PASS"
    )
    print(
        f"common-zero unions of <=3 unit masks, Q=2..40: "
        f"rows={masks['union_rows']}; PASS"
    )
    print("common-Q threshold: B3=2+3(A-1)=3A-1=6*ceil(Q/14)-4 <= Q")
    print("saving against B5: B5-B3=4*(ceil(Q/14)-1)")
    print("Tournament Analysis:")
    print("  candidate vertices=the 6 pairs among b1,b2,b3,b4")
    print("  gauge=lexicographic pair order; scores=(0,1,2,3,4,5); cycles=0; SCCs=6; H=1")
    print("  faithful quotient=four deletion circuits -> equality classes in five-mask nerve")
    print("  alternate vertices=pair candidates|gaps|residues|masks|proof obligations")
    print("  destroyed by runner tournament=Haar mass|additive compatibility|mask equality")
    print("VERDICT: mu(BAD)>60/637 forces the unique AP completion and a three-mask beat threshold")


if __name__ == "__main__":
    main()

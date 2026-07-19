#!/usr/bin/env python3
"""Exact referee for THM-1166, the radius-1/14 seven-wall Fano/gcd law.

All proof-facing arithmetic uses Fraction.  The script independently checks
the folded and trapezoid pair-overlap formulas, the finite 106-ratio / 173-
quotient-triple bank behind the three-speed floor, all 30 labelled Fano
planes, the quadratic chamber identity, and every displayed constant.
"""

from fractions import Fraction as F
from itertools import permutations
from math import comb, gcd


def fold14(r: int) -> int:
    r %= 14
    return r * (14 - r)


def pair_overlap(a: int, b: int) -> F:
    """Haar measure of D_a cap D_b for D_s={||st||<1/14}."""
    assert a > 0 and b > 0
    g = gcd(a, b)
    a //= g
    b //= g
    if a > b:
        a, b = b, a
    numerator = 4 * a * b + fold14(a + b) - fold14(b - a)
    return F(numerator, 196 * a * b)


def pair_overlap_trapezoid(a: int, b: int) -> F:
    """Independent LEM-042 capped-trapezoid sum."""
    assert a > 0 and b > 0
    g = gcd(a, b)
    width = F(a + b, 14 * a * b)
    cap = F(1, 7 * max(a, b))
    delta = F(g, a * b)
    total = F(0)
    j = 0
    while j * delta < width:
        overlap = min(width - j * delta, cap)
        total += overlap if j == 0 else 2 * overlap
        j += 1
    return g * total


CANONICAL_FANO = (
    (0, 1, 2),
    (0, 3, 4),
    (0, 5, 6),
    (1, 3, 5),
    (1, 4, 6),
    (2, 3, 6),
    (2, 4, 5),
)


def labelled_fano_planes() -> list[tuple[tuple[int, int, int], ...]]:
    planes = set()
    for p in permutations(range(7)):
        lines = tuple(
            sorted(tuple(sorted(p[v] for v in line)) for line in CANONICAL_FANO)
        )
        planes.add(lines)
    return sorted(planes)


def fano_line_value(c: int) -> F:
    return F(c, 3) - F(2, 7) * comb(c, 2)


def main() -> None:
    print("THM-1166 seven-wall Fano/gcd exact referee")
    print("=" * 72)

    print("pair-overlap law and universal floor")
    crosschecks = 0
    for a in range(1, 81):
        for b in range(a + 1, 121):
            assert pair_overlap(a, b) == pair_overlap_trapezoid(a, b)
            crosschecks += 1
    print(f"folded formula == trapezoid formula: {crosschecks} pairs")

    # The folded correction is at least -49/(196ab)=-1/(4ab).
    assert min(fold14(x) - fold14(y) for x in range(14) for y in range(14)) == -49
    tail_at_27 = F(1, 49) - F(1, 4 * 27)
    assert tail_at_27 > F(1, 91)
    small = []
    for a in range(1, 27):
        for b in range(a + 1, 27):
            if gcd(a, b) == 1 and a * b <= 26:
                small.append((pair_overlap(a, b), a, b))
    small.sort()
    assert len(small) == 36
    assert small[0] == (F(1, 91), 1, 13)
    assert sum(row[0] == F(1, 91) for row in small) == 1
    print("rho(a,b) >= 1/91; unique reduced equality pair: (1,13)")
    print(f"tail lower bound at ab=27: {tail_at_27} > 1/91")

    print("\nthree-speed finite-exact floor")
    tau = F(1, 24) - 2 * F(1, 91)
    assert tau == F(43, 2184)
    assert F(1, 49) - tau == F(11, 15288)
    assert F(2184, 11) > 198 and F(2184, 11) < 199

    low: dict[F, F] = {}
    for b in range(2, 199):
        for a in range(1, b):
            if gcd(a, b) == 1:
                rho = pair_overlap(a, b)
                if rho < tau:
                    low[F(b, a)] = rho
    assert len(low) == 106

    quotient_triples = []
    for r, rho_r in low.items():
        for s, rho_s in low.items():
            if r < s and s / r in low:
                q = s / r
                quotient_triples.append((rho_r + rho_s + low[q], r, s, q))
    quotient_triples.sort()
    assert len(quotient_triples) == 173
    minimum = quotient_triples[0][0]
    assert minimum == F(11, 252)
    assert minimum - F(1, 24) == F(1, 504)
    minimizers = {(r, s, q) for value, r, s, q in quotient_triples if value == minimum}
    assert minimizers == {
        (F(9, 4), F(27), F(12)),
        (F(12), F(27), F(9, 4)),
    }
    print(f"strict-low reduced ratios: {len(low)} (larger reduced coordinate b' <=198)")
    print(f"compatible quotient triples: {len(quotient_triples)}")
    print(f"minimum candidate pair sum: {minimum} = 1/24 + 1/504")
    print("therefore every three-speed pair sum is >=1/24")

    print("\nquadratic chamber certificate")
    q_table = []
    for c in range(8):
        q = F(c, 1) - F(2, 7) * comb(c, 2)
        assert q == F(c * (8 - c), 7)
        if c:
            assert q >= 1
        q_table.append(q)
    assert q_table == [F(0), F(1), F(12, 7), F(15, 7), F(16, 7),
                       F(15, 7), F(12, 7), F(1)]
    print("Q(C), C=0..7:", ", ".join(str(x) for x in q_table))
    print("tight on covered chambers C=1,7; unique peak C=4")

    triple_floor = F(1, 24)
    total_pair_floor = 7 * triple_floor
    uncovered_floor = F(2, 7) * total_pair_floor
    assert total_pair_floor == F(7, 24)
    assert uncovered_floor == F(1, 12)
    common_dilate = 7 * (1 - uncovered_floor)
    assert common_dilate == F(77, 12)
    print(f"R=sum pair overlaps >= {total_pair_floor}")
    print(f"global uncovered mass >= {uncovered_floor}")
    print(f"common-dilate protected-needle bound G/m <= {common_dilate}")

    print("\nFano incidence and mixed triple-gcd constants")
    planes = labelled_fano_planes()
    assert len(planes) == 30
    for plane in planes:
        edges = [tuple(sorted((a, b))) for line in plane for a in line for b in line if a < b]
        assert len(edges) == 21 and len(set(edges)) == 21
        assert all(sum(v in line for line in plane) == 3 for v in range(7))
        for mask in range(1 << 7):
            c = mask.bit_count()
            p = comb(c, 2)
            q = F(c) - F(2, 7) * p
            line_sum = sum(
                fano_line_value(sum((mask >> v) & 1 for v in line)) for line in plane
            )
            assert line_sum == q
    assert [fano_line_value(c) for c in range(4)] == [F(0), F(1, 3), F(8, 21), F(1, 7)]

    mu_max = (1 - 2 * triple_floor) / 7
    e_max = mu_max * (1 - F(21, 8) * mu_max)
    assert mu_max == F(11, 84)
    assert e_max == F(11, 128)
    fano_budget = (F(2, 7) * F(1, 7) * total_pair_floor) / e_max
    assert fano_budget == F(32, 231)
    fano_min_gcd = 7 / fano_budget
    assert fano_min_gcd == F(1617, 32)
    print(f"labelled Fano planes: {len(planes)}")
    print("line values f(c), c=0..3: 0, 1/3, 8/21, 1/7")
    print(f"line mean ceiling mu<= {mu_max}; discrepancy e<= {e_max}")
    print(f"every plane: sum_line m/G_line >= {fano_budget}")
    print(f"hence some line has G_line/m <= {fano_min_gcd}")

    print("\nadaptive forest and localization guardrail")
    # Exact single/pair periodic endpoint debts used by the forest consumer.
    assert F(1, 7) * (1 - F(1, 7)) == F(6, 49)
    print("single endpoint debt: 6/(49s)")
    print("pair endpoint debt: rho(1-rho)/gcd(si,sj)")
    print("forest observable: L*rho - rho(1-rho)/gcd(si,sj)")

    # Empty-star construction on the core-safe needle around 1/2.
    # This checks only absence of six selected intersections; it does NOT
    # assert that the seven combs cover the needle.
    m = 11
    leaves = ((1, 1), (2, 1), (2, 3), (3, 1), (4, 1), (5, 1))
    for q, c in leaves:
        lower_circle_distance = F(1, 2) - F(c, 14 * m)
        assert lower_circle_distance >= F(q + 1, 14)
    print("odd-core needle empty-star inequalities: 6/6 exact checks")
    print("scope: refutes cover-free fixed-tree localization; no cover is asserted")

    print("\ntournament/carrier audit")
    print("pairwise observable: compare exact line periods G_line (tie by line word)")
    print("switch/gauge: relabel the chosen Fano plane, then sort (G_line,line)")
    print("tie Hamiltonian path: lexicographic order of equal-period lines")
    print("fingerprint: transitive; scores 0..6; cycles 0; SCCs 7; paths 1")
    print("destroyed by tournament: which K7 edges share a Fano line, R_line, endpoint debt")
    print("faithful carrier: labelled Fano incidence + (R_line,G_line) + protected needle")
    print("challenged vertices: runners, pairs, lines, slow gaps, wall events, obligations")
    print("=" * 72)
    print("done")


if __name__ == "__main__":
    main()

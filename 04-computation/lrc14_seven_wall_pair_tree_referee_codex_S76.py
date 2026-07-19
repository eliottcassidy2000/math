#!/usr/bin/env python3
"""Exact referee for THM-1166 (seven-wall pair/tree and Fano/gcd discrepancy).

The program is intentionally dependency-free.  It checks the folded and
tent forms of the one-fourteenth pair overlap, both finite low-channel
classifications, quotient compatibility of the channels, the intermediate
1/24 and sharp 51/1183 three-speed pair-sum floors, the Fano-line and
spanning-tree consequences,
the quadratic wall algebra, and exact sample instances of the periodic
interval discrepancy lemma.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd
from pathlib import Path


ONE = F(1)
PAIR_FLOOR = F(1, 91)
PAIR_BASE = F(1, 49)
TRIPLE_BAR = F(1, 63)
TRIPLE_SUM_INTERMEDIATE = F(1, 24)
TRIPLE_SUM_FLOOR = F(51, 1183)
TRIPLE_SUM_CHANNEL_CUT = F(43, 2184)
TREE_FLOOR = F(110, 1183)
FANO_TOTAL_FLOOR = F(51, 169)
LOCAL_TREE_COEFFICIENT = F(102, 1183)


def require(condition: bool, message: str) -> None:
    """Optimization-stable certificate check."""
    if not condition:
        raise RuntimeError(message)


def fold(r: int, modulus: int) -> int:
    r %= modulus
    return r * (modulus - r)


def tent_sum(z: F) -> F:
    """sum_(m in Z) (z-|m|)_+, in closed form."""
    q = z.numerator // z.denominator
    return (2 * q + 1) * z - q * (q + 1)


def rho_fold(a: int, b: int) -> F:
    """Exact Haar measure of D_a intersect D_b at radius 1/14."""
    require(a > 0 and b > 0, "pair speeds must be positive")
    if a > b:
        a, b = b, a
    g = gcd(a, b)
    modulus = 14 * g
    numerator = (
        4 * a * b
        + fold(a + b, modulus)
        - fold(b - a, modulus)
    )
    return F(numerator, 196 * a * b)


def rho_tent(a: int, b: int) -> F:
    """Independent reduced sawtooth/tent evaluation of the same measure."""
    if a > b:
        a, b = b, a
    g = gcd(a, b)
    a //= g
    b //= g
    return (
        tent_sum(F(a + b, 14)) - tent_sum(F(b - a, 14))
    ) / (a * b)


def low_reduced_channels() -> list[tuple[int, int, F]]:
    out: list[tuple[int, int, F]] = []
    for a in range(1, 56):
        for b in range(a + 1, 56):
            if a * b <= 55 and gcd(a, b) == 1:
                value = rho_fold(a, b)
                if value < TRIPLE_BAR:
                    out.append((a, b, value))
    return out


def triple_sum_low_channels() -> list[tuple[int, int, F]]:
    """Complete reduced channels below 1/24-2/91.

    The folded lower bound proves that products at least 348 are outside the
    list, so this finite loop is the entire classification rather than a box
    sample.
    """
    out: list[tuple[int, int, F]] = []
    for a in range(1, 348):
        for b in range(a + 1, 348):
            if a * b <= 347 and gcd(a, b) == 1:
                value = rho_fold(a, b)
                if value < TRIPLE_SUM_CHANNEL_CUT:
                    out.append((a, b, value))
    return out


def oriented_reduced_ratios(product_cap: int) -> list[F]:
    """All oriented nonunit reduced ratios with numerator*denominator <= cap."""
    return sorted(
        {
            F(numerator, denominator)
            for numerator in range(1, product_cap + 1)
            for denominator in range(1, product_cap + 1)
            if numerator != denominator
            and numerator * denominator <= product_cap
            and gcd(numerator, denominator) == 1
        }
    )


def mst_weight(speeds: tuple[int, ...]) -> F:
    edges = sorted(
        (
            (rho_fold(speeds[i], speeds[j]), i, j)
            for i in range(7)
            for j in range(i + 1, 7)
        ),
        reverse=True,
    )
    parent = list(range(7))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    total = F(0)
    used = 0
    for weight, i, j in edges:
        ri, rj = root(i), root(j)
        if ri != rj:
            parent[ri] = rj
            total += weight
            used += 1
            if used == 6:
                break
    require(used == 6, "K7 maximum spanning tree did not connect")
    return total


FANO_LINES = (
    (0, 1, 2),
    (0, 3, 4),
    (0, 5, 6),
    (1, 3, 5),
    (1, 4, 6),
    (2, 3, 6),
    (2, 4, 5),
)


def fano_line_sums(speeds: tuple[int, ...]) -> list[F]:
    ans: list[F] = []
    for line in FANO_LINES:
        ans.append(
            sum(
                (
                    rho_fold(speeds[i], speeds[j])
                    for i, j in combinations(line, 2)
                ),
                F(0),
            )
        )
    return ans


def danger_intervals(speed: int) -> list[tuple[F, F]]:
    """Closed representatives in [0,1]; endpoints do not affect measure."""
    h = F(1, 14 * speed)
    pieces: list[tuple[F, F]] = []
    for k in range(speed):
        center = F(k, speed)
        left, right = center - h, center + h
        if left < 0:
            pieces.extend(((F(0), right), (ONE + left, ONE)))
        elif right > 1:
            pieces.extend(((left, ONE), (F(0), right - ONE)))
        else:
            pieces.append((left, right))
    pieces.sort()
    merged: list[tuple[F, F]] = []
    for left, right in pieces:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return merged


def intersect_intervals(
    lefts: list[tuple[F, F]], rights: list[tuple[F, F]]
) -> list[tuple[F, F]]:
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(lefts) and j < len(rights):
        left = max(lefts[i][0], rights[j][0])
        right = min(lefts[i][1], rights[j][1])
        if left < right:
            out.append((left, right))
        if lefts[i][1] < rights[j][1]:
            i += 1
        else:
            j += 1
    return out


def interval_mass(parts: list[tuple[F, F]], interval: tuple[F, F]) -> F:
    left, right = interval
    return sum(
        (max(F(0), min(right, b) - max(left, a)) for a, b in parts),
        F(0),
    )


def tournament_fingerprint(speeds: tuple[int, ...]) -> tuple:
    """Orient high edges up the speed order and strict-low edges down it."""
    edge: set[tuple[int, int]] = set()
    for i in range(7):
        for j in range(i + 1, 7):
            if rho_fold(speeds[i], speeds[j]) >= TRIPLE_BAR:
                edge.add((i, j))
            else:
                edge.add((j, i))
    scores = tuple(sorted(sum((i, j) in edge for j in range(7) if j != i) for i in range(7)))
    triangles = sum(
        (i, j) in edge and (j, k) in edge and (k, i) in edge
        for i, j, k in permutations(range(7), 3)
    ) // 3

    # Standard insertion proof of Redei's Hamiltonian-path theorem.
    path: list[int] = []
    for vertex in range(7):
        pos = 0
        while pos < len(path) and (path[pos], vertex) in edge:
            pos += 1
        path.insert(pos, vertex)
    require(
        all((path[i], path[i + 1]) in edge for i in range(6)),
        "Redei insertion path is not directed",
    )

    ham_count = sum(
        all((p[i], p[i + 1]) in edge for i in range(6))
        for p in permutations(range(7))
    )
    return scores, triangles, tuple(path), ham_count


def fmt(value: F) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    # Independent folded/tent agreement, including nonprimitive pairs.
    checked = 0
    for a in range(1, 121):
        for b in range(a + 1, 121):
            require(
                rho_fold(a, b) == rho_tent(a, b),
                f"folded/tent mismatch at {(a, b)}",
            )
            checked += 1
    print("THM-1166 exact seven-wall pair/tree and Fano/gcd referee")
    print(f"folded-vs-tent exact pair checks: {checked}")

    # The analytic tail 1/49-1/(4ab) handles all coprime ab >= 27.
    require(
        F(1, 49) - F(1, 4 * 27) > PAIR_FLOOR,
        "pair-floor analytic tail does not start at ab=27",
    )
    small_floor = [
        (rho_fold(a, b), a, b)
        for a in range(1, 27)
        for b in range(a + 1, 27)
        if a * b <= 26 and gcd(a, b) == 1
    ]
    floor_value, floor_a, floor_b = min(small_floor)
    require(
        (floor_value, floor_a, floor_b) == (PAIR_FLOOR, 1, 13),
        "wrong reduced pair-floor minimizer",
    )
    require(
        sum(value == PAIR_FLOOR for value, _, _ in small_floor) == 1,
        "reduced pair-floor minimizer is not unique",
    )
    print("global pair floor: rho >= 1/91; unique reduced minimizer (1,13)")
    print("pair-floor tail starts at ab=27; finite cells ab<=26: PASS")

    # The same analytic tail handles ab >= 56 for the strict 1/63 cut.
    require(
        F(1, 49) - F(1, 4 * 56) > TRIPLE_BAR,
        "triple-bar analytic tail does not start at ab=56",
    )
    low = low_reduced_channels()
    expected_low = [
        (1, 10, F(1, 70)),
        (1, 11, F(1, 77)),
        (1, 12, F(1, 84)),
        (1, 13, F(1, 91)),
        (2, 11, F(1, 77)),
        (3, 10, F(1, 70)),
        (3, 11, F(1, 77)),
    ]
    require(low == expected_low, "strict-low channel classification mismatch")
    print("strict-low channels rho<1/63 (complete):")
    print("  " + ", ".join(f"({a},{b}):{fmt(v)}" for a, b, v in low))

    low_pairs = {(a, b) for a, b, _ in low}
    quotient_checks = 0
    for a, b, _ in low:
        for c, d, _ in low:
            quotient_checks += 1
            ratio = F(b, a) * F(d, c)
            reduced = (ratio.denominator, ratio.numerator)
            require(reduced not in low_pairs, f"low quotient/product at {reduced}")
            require(
                rho_fold(*reduced) >= TRIPLE_BAR,
                f"quotient/product misses triple bar at {reduced}",
            )
    print(f"oriented low-channel product/quotient checks: {quotient_checks}; violations: 0")
    print("therefore every speed triple has an edge rho>=1/63")
    sharp_triple = (1, 12, 27)
    sharp_values = tuple(rho_fold(*edge) for edge in combinations(sharp_triple, 2))
    require(max(sharp_values) == TRIPLE_BAR, "sharp triple does not attain 1/63")
    print(
        "triple constant sharp at (1,12,27): "
        + ", ".join(fmt(v) for v in sharp_values)
    )

    # A counterexample to the 1/24 three-speed pair-sum floor would force
    # every edge below 1/24-2/91.  The analytic tail makes that a finite
    # reduced-channel classification, and sorted speed ratios turn triple
    # compatibility into exact multiplication of two rational channels.
    require(
        TRIPLE_SUM_CHANNEL_CUT + 2 * PAIR_FLOOR == TRIPLE_SUM_INTERMEDIATE,
        "three-speed channel cut arithmetic mismatch",
    )
    require(
        F(1, 49) - F(1, 4 * 348) > TRIPLE_SUM_CHANNEL_CUT,
        "three-speed pair-sum tail does not start at ab=348",
    )
    require(
        F(1, 49) - F(1, 4 * 347) <= TRIPLE_SUM_CHANNEL_CUT,
        "reported three-speed pair-sum tail is not the first integral cut",
    )
    sum_low = triple_sum_low_channels()
    require(len(sum_low) == 106, "wrong three-speed sub-threshold channel count")
    require(
        max(a * b for a, b, _ in sum_low) == 327,
        "wrong largest occupied product in sub-threshold channels",
    )
    require(
        max(b for _, b, _ in sum_low) == 167,
        "wrong largest numerator in sub-threshold channels",
    )
    sum_low_by_ratio = {F(b, a): value for a, b, value in sum_low}
    require(
        len(sum_low_by_ratio) == len(sum_low),
        "reduced channel ratios are not unique",
    )
    compatible: list[tuple[F, F, F]] = []
    for first, first_value in sum_low_by_ratio.items():
        for second, second_value in sum_low_by_ratio.items():
            product = first * second
            if product in sum_low_by_ratio:
                compatible.append(
                    (
                        first,
                        second,
                        first_value + second_value + sum_low_by_ratio[product],
                    )
                )
    require(len(compatible) == 173, "wrong quotient-compatible triple count")
    candidate_minimum = min(total for _, _, total in compatible)
    require(
        candidate_minimum == F(11, 252),
        "wrong compatible three-speed pair-sum minimum",
    )
    require(
        candidate_minimum == TRIPLE_SUM_INTERMEDIATE + F(1, 504),
        "three-speed finite margin arithmetic mismatch",
    )
    minimizers = {
        (first, second)
        for first, second, total in compatible
        if total == candidate_minimum
    }
    require(
        minimizers == {(F(12), F(9, 4)), (F(9, 4), F(12))},
        "wrong compatible three-speed minimizers",
    )
    print("three-speed pair-sum finite-exact certificate:")
    print("  counterexample channel cut: rho < 43/2184")
    print("  analytic tail starts at ab=348; complete reduced channels: 106")
    print("  quotient-compatible triples: 173")
    print("  candidate minimum: 11/252=1/24+1/504")
    print("therefore every speed triple has total pair mass >=1/24 (intermediate)")

    # Sharp upgrade.  If p1<=p2<=p3 are the products of the three reduced
    # pair channels, the folded lower bound forces p1<=41 and p2<=57 in any
    # triple whose pair sum is at most 51/1183.  Two edges share a vertex, so
    # enumerating their two oriented ratios covers the full triangle.
    first_product_threshold = F(3, 4) / (3 * PAIR_BASE - TRIPLE_SUM_FLOOR)
    second_product_threshold = F(1, 2) / (
        PAIR_FLOOR + 2 * PAIR_BASE - TRIPLE_SUM_FLOOR
    )
    require(
        first_product_threshold == F(8281, 200),
        "wrong smallest-product threshold for sharp triple sum",
    )
    require(
        second_product_threshold == F(8281, 144),
        "wrong second-product threshold for sharp triple sum",
    )
    require(
        41 < first_product_threshold < 42,
        "smallest-product integral cutoff is not 41",
    )
    require(
        57 < second_product_threshold < 58,
        "second-product integral cutoff is not 57",
    )
    first_ratios = oriented_reduced_ratios(41)
    second_ratios = oriented_reduced_ratios(57)
    require(len(first_ratios) == 124, "wrong oriented ratio count at product 41")
    require(len(second_ratios) == 184, "wrong oriented ratio count at product 57")
    sharp_candidates: list[tuple[F, F, F]] = []
    for first in first_ratios:
        for second in second_ratios:
            if first == second:
                continue
            total = (
                rho_fold(first.numerator, first.denominator)
                + rho_fold(second.numerator, second.denominator)
                + rho_fold(
                    (first / second).numerator,
                    (first / second).denominator,
                )
            )
            sharp_candidates.append((total, first, second))
    require(len(sharp_candidates) == 22692, "wrong sharp triple candidate count")
    sharp_minimum = min(total for total, _, _ in sharp_candidates)
    require(sharp_minimum == TRIPLE_SUM_FLOOR, "wrong sharp triple-sum minimum")
    sharp_minimizers = {
        (first, second)
        for total, first, second in sharp_candidates
        if total == sharp_minimum
    }
    require(
        sharp_minimizers == {(F(1, 13), F(13)), (F(13), F(1, 13))},
        "wrong sharp triple-sum minimizers",
    )
    require(
        sum((rho_fold(*edge) for edge in combinations((1, 13, 169), 2)), F(0))
        == TRIPLE_SUM_FLOOR,
        "displayed sharp triple does not attain 51/1183",
    )
    print("sharp three-speed pair-sum certificate:")
    print("  counterexample forces two smallest reduced products <=41 and <=57")
    print("  oriented ratio banks: 124 and 184; configurations: 22692")
    print("  exact minimum: 51/1183, uniquely at the ratio triple (1,13,169)")

    # If the high graph is disconnected, it has two components.  Let m be
    # the largest crossing weight.  Every one of the five internal tree
    # edges has weight at least max(1/63,51/1183-2m), so the completed tree
    # has weight at least 5*max(1/63,51/1183-2m)+m.  Its two
    # affine branches meet at the exact universal minimum below.
    tree_tradeoff_crossing = (TRIPLE_SUM_FLOOR - TRIPLE_BAR) / 2
    require(
        tree_tradeoff_crossing == F(145, 10647),
        "wrong tree-tradeoff crossing point",
    )
    require(
        5 * TRIPLE_BAR + tree_tradeoff_crossing == TREE_FLOOR,
        "wrong disconnected-high-graph tree floor",
    )
    require(
        6 * TRIPLE_BAR > TREE_FLOOR,
        "connected-high-graph branch does not clear tree floor",
    )
    print("high-component tree tradeoff:")
    print("  min_m [5*max(1/63,51/1183-2m)+m] = 110/1183")
    print("  connected-high branch gives 6/63=2/21 >110/1183")

    # Validate the Fano incidence design itself.
    pair_mult = {(i, j): 0 for i in range(7) for j in range(i + 1, 7)}
    vertex_mult = [0] * 7
    for line in FANO_LINES:
        for i in line:
            vertex_mult[i] += 1
        for i, j in combinations(line, 2):
            pair_mult[min(i, j), max(i, j)] += 1
    require(vertex_mult == [3] * 7, "Fano vertex multiplicities are not three")
    require(set(pair_mult.values()) == {1}, "Fano pair multiplicities are not one")

    # Exhaust a nontrivial seven-packet box.
    packet_checks = 0
    min_tree: tuple[F, tuple[int, ...]] | None = None
    min_total: tuple[F, tuple[int, ...]] | None = None
    for speeds in combinations(range(1, 19), 7):
        packet_checks += 1
        tree = mst_weight(speeds)
        require(tree >= TREE_FLOOR, f"MST floor failure at {speeds}")
        line_sums = fano_line_sums(speeds)
        require(
            all(value >= TRIPLE_SUM_FLOOR for value in line_sums),
            f"Fano line floor failure at {speeds}",
        )
        total = sum(line_sums, F(0))
        require(total >= FANO_TOTAL_FLOOR, f"Fano total floor failure at {speeds}")
        if min_tree is None or tree < min_tree[0]:
            min_tree = (tree, speeds)
        if min_total is None or total < min_total[0]:
            min_total = (total, speeds)
    require(min_tree is not None and min_total is not None, "empty seven-packet audit")
    print(f"seven-packet exact box C(18,7): {packet_checks}; violations: 0")
    print(f"universal spanning-tree floor: {fmt(TREE_FLOOR)}")
    print(f"box minimum MST: {fmt(min_tree[0])} at {min_tree[1]}")
    print(f"Fano total-pair floor: 7*(51/1183)={fmt(FANO_TOTAL_FLOOR)}")
    print(f"box minimum total pair mass: {fmt(min_total[0])} at {min_total[1]}")

    # Pointwise quadratic wall and its Fano decomposition, all 2^7 states.
    for mask in range(1 << 7):
        active = {i for i in range(7) if mask & (1 << i)}
        count = len(active)
        q = F(count * (8 - count), 7)
        if count:
            require(q >= 1, f"quadratic wall failure at mask {mask}")
        fano_numerator = 0
        for line in FANO_LINES:
            c_line = len(active.intersection(line))
            fano_numerator += 7 * c_line - 6 * (c_line * (c_line - 1) // 2)
        require(21 * q == fano_numerator, f"Fano decomposition failure at mask {mask}")
    print("quadratic Q(C)=C(8-C)/7 and all 128 Fano decompositions: PASS")

    # Exact samples of the sharp universal periodic-indicator debt
    # |mu(I cap B)-|I|rho(B)| <= rho(B)(1-rho(B))/gcd(a,b).
    discrepancy_checks = 0
    worst_ratio = F(0)
    intervals = (
        (F(1, 17), F(13, 29)),
        (F(7, 53), F(41, 71)),
        (F(2, 7), F(5, 6)),
        (F(1, 7), F(6, 7)),
    )
    for a in range(1, 31):
        da = danger_intervals(a)
        for b in range(a + 1, 31):
            pair = intersect_intervals(da, danger_intervals(b))
            require(
                sum((right - left for left, right in pair), F(0)) == rho_fold(a, b),
                f"direct interval pair mass mismatch at {(a, b)}",
            )
            for interval in intervals:
                length = interval[1] - interval[0]
                error = abs(interval_mass(pair, interval) - length * rho_fold(a, b))
                density = rho_fold(a, b)
                bound = density * (1 - density) / gcd(a, b)
                require(error <= bound, f"periodic discrepancy failure at {(a, b, interval)}")
                worst_ratio = max(worst_ratio, error / bound)
                discrepancy_checks += 1
    print(f"periodic interval discrepancy exact samples: {discrepancy_checks}; violations: 0")
    print(f"sample worst error/[rho(1-rho)/gcd]: {fmt(worst_ratio)}")
    require(
        F(2, 7) * FANO_TOTAL_FLOOR == LOCAL_TREE_COEFFICIENT,
        "positioned-tree coefficient arithmetic mismatch",
    )
    print("positioned-tree lower bound (every interval):")
    print("  max_T sum_T mu(I cap D_i cap D_j) >= (102/1183)|I|-(2/7)E")
    print("  E=sum_(i<j)rho_ij(1-rho_ij)/gcd(s_i,s_j)")

    harmonic7 = sum((F(1, i) for i in range(1, 8)), F(0))
    require(harmonic7 == F(363, 140), "H7 arithmetic mismatch")
    # Since rho<=1/7, rho(1-rho)<=6/49.  This gives a direct but weak
    # common-dilate consequence from the positioned inequality.  The global
    # Hunter/tree period argument below is dramatically stronger.
    gcd_common_dilate = (F(18, 7) + harmonic7 / 2) / FANO_TOTAL_FLOOR
    require(gcd_common_dilate == F(61009, 4760), "gcd common-dilate mismatch")
    global_safe_floor = TREE_FLOOR
    periodic_common_dilate = 1 - global_safe_floor
    require(
        periodic_common_dilate == F(1073, 1183),
        "periodic common-dilate interval constant mismatch",
    )
    protected_common_dilate = 7 * periodic_common_dilate
    require(
        protected_common_dilate == F(1073, 169),
        "protected common-dilate constant mismatch",
    )
    print("covered-interval necessary inequality:")
    print("  (51/169)|I| <= (1/2)sum_i(1/s_i) + E")
    print(f"raw gcd-error common-dilate consequence: g|I| <= {fmt(gcd_common_dilate)}")
    print(f"global Hunter safe-set floor: {fmt(global_safe_floor)}")
    print(f"periodic common-dilate consequence: g|I| <= {fmt(periodic_common_dilate)}")
    print(f"protected-needle consequence: g/max(core) <= {fmt(protected_common_dilate)}")

    # Localization counterexample: the difference-comb inclusion makes this
    # exact for every N, not merely a finite scan.
    for n in (1, 2, 5, 10, 100, 1000):
        for i in range(6):
            pair = intersect_intervals(
                danger_intervals(n + i), danger_intervals(n + i + 1)
            )
            require(
                interval_mass(pair, (F(1, 7), F(6, 7))) == 0,
                f"consecutive localization counterexample failed at {(n, i)}",
            )
    print("localization obstruction: for every N,")
    print("  mu([1/7,6/7] cap D_N cap D_(N+1)) = 0, interval length 5/7")
    print("  the six-edge order path on {N,...,N+6} therefore has local weight 0")
    print("  hence no positive per-edge local floor with O(1/N) endpoint error")

    # Use a genuinely cyclic gauge sample rather than the box minimizer, whose
    # switched tournament happens to be transitive.
    sample = (1, 2, 3, 4, 5, 6, 10)
    scores, triangles, path, ham_count = tournament_fingerprint(sample)
    print("tournament audit (vertices=speeds; high edge forward, strict-low backward)")
    print(f"sample packet: {sample}")
    print(f"score histogram: {scores}; directed triangles: {triangles}")
    print(f"tie Hamiltonian path: {path}; Hamiltonian-path count: {ham_count}")
    print("preserved: global low/high channel, Fano line debt, MST floor")
    print("destroyed: interval phase, tooth address, owner chronology, gcd cell position")

    source_hash = sha256(Path(__file__).read_bytes()).hexdigest()
    print(f"source_sha256={source_hash}")


if __name__ == "__main__":
    main()

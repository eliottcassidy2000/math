#!/usr/bin/env python3
"""Exact referee for THM-2086's Fourier split and lacunary cone.

The analytic proof uses Fejer convergence and bounded variation.  This
companion independently replays the exact finite arithmetic: every bounded
THM-2081 tree, the hostile channel decomposition, the odd reduced-ratio
spectrum, and a broad exact bank for the BV mixing inequality.  Runtime
checks remain active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import gcd, pi, sin


MAX_CORE_SPEED = 24
FULL_DIVISOR_MASK = (1 << 13) - 1
DIVISOR_MASK = {
    q: sum(1 << (d - 2) for d in range(2, 15) if q % d == 0)
    for q in range(1, MAX_CORE_SPEED + 1)
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def divisor_complete(Q: tuple[int, ...]) -> bool:
    mask = 0
    for q in Q:
        mask |= DIVISOR_MASK[q]
    return mask == FULL_DIVISOR_MASK


def hereditarily_primitive(Q: tuple[int, ...]) -> bool:
    return all(gcd(*(Q[:i] + Q[i + 1 :])) == 1 for i in range(len(Q)))


def mixed_overlap(h: int, q: int) -> F:
    """Exact I(h,q)=measure(E_h intersect D_q), from THM-2080."""
    g = gcd(h, q)
    a, b = h // g, q // g
    x, y = F(a % 14, 14), F(b % 7, 7)
    correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
    return F(2, 49) + F(2, a * b) * correction


def global_pair(p: int, q: int) -> F:
    """Exact rho(p,q)=measure(D_p intersect D_q), from THM-1166."""
    g = gcd(p, q)
    a, b = p // g, q // g

    def fold(r: int) -> int:
        residue = r % 14
        return residue * (14 - residue)

    return F(1, 49) + F(fold(a + b) - fold(b - a), 196 * a * b)


def restricted_pair(
    h: int,
    p: int,
    q: int,
    cache: dict[tuple[int, int, int], F],
) -> F:
    """Exact measure(D_p intersect D_q intersect E_h^c)."""
    key = (h, min(p, q), max(p, q))
    if key in cache:
        return cache[key]

    boundaries = {F(0), F(1)}
    for speed, radius in ((h, F(1, 7)), (p, F(1, 14)), (q, F(1, 14))):
        for k in range(speed + 1):
            center = F(k, speed)
            for sign in (-1, 1):
                point = center + sign * radius / speed
                if 0 <= point <= 1:
                    boundaries.add(point)

    points = sorted(boundaries)
    measure = F(0)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        if (
            circle_distance(h * midpoint) >= F(1, 7)
            and circle_distance(p * midpoint) < F(1, 14)
            and circle_distance(q * midpoint) < F(1, 14)
        ):
            measure += right - left
    cache[key] = measure
    return measure


def genuine_residual(
    h: int,
    p: int,
    q: int,
    restricted: F,
) -> F:
    """The exact all-nonzero Fourier channel R_h(p,q)."""
    epsilon_p = mixed_overlap(h, p) - F(2, 49)
    epsilon_q = mixed_overlap(h, q) - F(2, 49)
    return (
        restricted
        - F(5, 7) * global_pair(p, q)
        + (epsilon_p + epsilon_q) / 7
    )


def maximum_spanning_tree(
    Q: tuple[int, ...],
    h: int,
    cache: dict[tuple[int, int, int], F],
) -> tuple[F, tuple[tuple[F, int, int], ...]]:
    used = {0}
    total = F(0)
    edges: list[tuple[F, int, int]] = []
    while len(used) < len(Q):
        best = None
        for i in used:
            for j in range(len(Q)):
                if j in used:
                    continue
                weight = restricted_pair(h, Q[i], Q[j], cache)
                candidate = (weight, i, j)
                if best is None or candidate > best:
                    best = candidate
        require(best is not None, "Prim search lost a vertex")
        total += best[0]
        edges.append(best)
        used.add(best[2])
    return total, tuple(edges)


def truncated_genuine_channel(h: int, p: int, q: int, cutoff: int) -> float:
    """Symmetric cubic-box truncation of the all-nonzero Fourier channel."""
    danger = {
        n: sin(pi * n / 7) / (pi * n)
        for n in range(-cutoff, cutoff + 1)
        if n != 0
    }
    guard_complement = {
        n: -sin(2 * pi * n / 7) / (pi * n)
        for n in range(-cutoff, cutoff + 1)
        if n != 0
    }
    total = 0.0
    for a, coefficient_a in danger.items():
        for b, coefficient_b in danger.items():
            numerator = -(a * p + b * q)
            if numerator % h != 0:
                continue
            c = numerator // h
            if c == 0 or abs(c) > cutoff:
                continue
            total += coefficient_a * coefficient_b * guard_complement[c]
    return total


def audit_bounded_terminal_bank(
    cache: dict[tuple[int, int, int], F],
) -> tuple[int, int, int, tuple[object, ...]]:
    core_count = 0
    pair_count = 0
    scalar_survivors = 0
    relative_survivors = 0
    worst = None

    for Q in combinations(range(1, MAX_CORE_SPEED + 1), 7):
        if not divisor_complete(Q) or not hereditarily_primitive(Q):
            continue
        core_count += 1
        for h in range(1, (8 * max(Q) - 1) // 3 + 1, 2):
            if 6 * h >= 16 * max(Q):
                continue
            pair_count += 1
            intersections = tuple(mixed_overlap(h, q) for q in Q)
            epsilons = tuple(value - F(2, 49) for value in intersections)
            deficit = -sum(epsilons, F(0))
            if deficit < 0:
                continue
            scalar_survivors += 1

            tree, edges = maximum_spanning_tree(Q, h, cache)
            degrees = [0] * 7
            bulk = F(0)
            genuine = F(0)
            edge_word = []
            for weight, i, j in edges:
                degrees[i] += 1
                degrees[j] += 1
                rho = global_pair(Q[i], Q[j])
                residual = genuine_residual(h, Q[i], Q[j], weight)
                reconstructed = (
                    F(5, 7) * rho
                    - (epsilons[i] + epsilons[j]) / 7
                    + residual
                )
                require(
                    reconstructed == weight,
                    f"edge split failed Q={Q}, h={h}, edge={(i, j)}",
                )
                bulk += F(5, 7) * rho
                genuine += residual
                edge_word.append((Q[i], Q[j], weight, residual))

            epsilon_term = sum(
                (F(7 - degrees[i], 7) * epsilons[i] for i in range(7)),
                F(0),
            )
            margin = tree - deficit
            require(
                margin == bulk + epsilon_term + genuine,
                f"tree split failed Q={Q}, h={h}",
            )
            if margin <= 0:
                relative_survivors += 1
            row = (
                margin,
                Q,
                h,
                tree,
                deficit,
                bulk,
                epsilon_term,
                genuine,
                tuple(degrees),
                tuple(edge_word),
            )
            if worst is None or row < worst:
                worst = row

    require(core_count == 131, f"core count changed: {core_count}")
    require(pair_count == 4120, f"pair count changed: {pair_count}")
    require(scalar_survivors == 1322, "scalar survivor count changed")
    require(relative_survivors == 0, "relative Hunter left a bounded survivor")
    require(worst is not None, "no bounded scalar survivor was audited")
    require(worst[0] == F(561797, 8288280), "hostile margin changed")
    require(worst[1] == (1, 9, 10, 11, 13, 14, 24), "hostile core changed")
    require(worst[2] == 23, "hostile guard changed")
    require(worst[5] == F(100421, 1177176), "hostile bulk term changed")
    require(worst[6] == -F(16117, 4512508), "hostile epsilon term changed")
    require(worst[7] == -F(2833331, 203062860), "hostile genuine term changed")
    require(worst[8] == (1, 1, 1, 1, 1, 2, 5), "hostile degree word changed")
    expected_edges = (
        (1, 14, F(1, 98), -F(11, 2254)),
        (14, 24, F(16, 1127), -F(11, 27048)),
        (24, 10, F(11, 840), -F(71, 118335)),
        (24, 9, F(1, 84), -F(281, 284004)),
        (24, 11, F(1, 88), -F(7663, 2082696)),
        (24, 13, F(71, 6279), -F(4181, 1230684)),
    )
    require(worst[9] == expected_edges, "hostile edge/channel word changed")
    return pair_count, scalar_survivors, len(cache), worst


def audit_fourier_truncations(worst: tuple[object, ...]) -> float:
    exact_channels = [((3, 1, 2), -F(40, 1029))]
    for p, q, _weight, residual in worst[9]:
        exact_channels.append(((23, p, q), residual))

    maximum_error = 0.0
    for (h, p, q), exact in exact_channels:
        approximation = truncated_genuine_channel(h, p, q, 640)
        error = abs(approximation - float(exact))
        maximum_error = max(maximum_error, error)
        require(error < 4e-7, f"Fourier truncation drifted at {(h, p, q)}: {error}")
    return maximum_error


def audit_odd_reduced_spectrum() -> tuple[int, tuple[tuple[int, int, F], ...]]:
    checked = 0
    low_rows = []
    for a in range(1, 21, 2):
        for b in range(1, 21):
            if gcd(a, b) != 1 or a * b > 20:
                continue
            value = mixed_overlap(a, b)
            if value <= F(1, 35):
                low_rows.append((a, b, value))
            checked += 1

    expected = (
        (1, 5, F(1, 35)),
        (1, 6, F(1, 42)),
        (3, 5, F(1, 35)),
        (11, 1, F(2, 77)),
    )
    require(tuple(low_rows) == expected, "odd-a low spectrum changed")
    tail_floor = F(2, 49) - F(1, 4 * 21)
    require(tail_floor > F(1, 35), "ab>=21 tail no longer clears 1/35")
    return checked, tuple(low_rows)


def audit_bv_mixing(
    cache: dict[tuple[int, int, int], F],
) -> tuple[int, tuple[F, int, int, int, F, F]]:
    checks = 0
    sharpest = None
    mixing_speeds = (1, 2, 3, 5, 8, 13, 21, 34, 55, 89)
    for h in range(1, 24, 2):
        for q in range(1, 25):
            outside_mass = F(1, 7) - mixed_overlap(h, q)
            target = outside_mass / 7
            for B in mixing_speeds:
                actual = restricted_pair(h, B, q, cache)
                error = abs(actual - target)
                bound = F(q + h, 3 * B)
                require(
                    error <= bound,
                    f"BV mixing failure h={h}, q={q}, B={B}: {error}>{bound}",
                )
                normalized = F(0) if q + h == 0 else error / bound
                row = (normalized, h, q, B, error, bound)
                if sharpest is None or row > sharpest:
                    sharpest = row
                checks += 1
    require(checks == 2880, f"BV bank count changed: {checks}")
    require(sharpest is not None, "empty BV bank")
    return checks, sharpest


def audit_seven_dividing_guard_cross_edges(
    cache: dict[tuple[int, int, int], F],
) -> int:
    """Exact broad bank for the 7|h low/high cross-edge identity."""
    checks = 0
    for h in (7, 21, 35, 49):
        for p in range(1, 37):
            if p % 7 == 0:
                continue
            epsilon_p = mixed_overlap(h, p) - F(2, 49)
            require(epsilon_p == 0, f"7|h low epsilon did not vanish: {(h, p)}")
            for q in range(7, 71, 7):
                epsilon_q = mixed_overlap(h, q) - F(2, 49)
                rho = global_pair(p, q)
                require(rho == F(1, 49), f"apex pair did not decorrelate: {(p, q)}")
                weight = restricted_pair(h, p, q, cache)
                residual = genuine_residual(h, p, q, weight)
                expected = F(5, 343) - epsilon_q / 7
                require(residual == 0, f"7|h genuine channel survived: {(h, p, q)}")
                require(weight == expected, f"7|h cross-edge formula failed: {(h, p, q)}")
                checks += 1
    require(checks == 1240, f"7|h cross-edge bank count changed: {checks}")
    return checks


def seven_dividing_guard_tree(
    Q: tuple[int, ...],
    h: int,
    cache: dict[tuple[int, int, int], F],
) -> tuple[F, tuple[tuple[int, int, F], ...]]:
    high = tuple(i for i, q in enumerate(Q) if q % 7 == 0)
    low = tuple(i for i, q in enumerate(Q) if q % 7 != 0)
    require(1 <= len(high) <= 5, f"invalid hereditary mod-7 profile: {Q}")
    require(len(low) >= 2, f"hereditary primitivity lost its two low entries: {Q}")

    # A deterministic spanning tree of the complete bipartite graph K_(low,high).
    edges = [(low[0], j) for j in high]
    edges.extend((i, high[0]) for i in low[1:])
    require(len(edges) == 6, f"bipartite tree edge count changed: {Q}")

    total = F(0)
    word = []
    for i, j in edges:
        weight = restricted_pair(h, Q[i], Q[j], cache)
        epsilon_high = mixed_overlap(h, Q[j]) - F(2, 49)
        expected = F(5, 343) - epsilon_high / 7
        require(
            genuine_residual(h, Q[i], Q[j], weight) == 0,
            f"finite terminal 7|h channel survived: {(Q, h, i, j)}",
        )
        require(weight == expected, f"finite terminal cross edge failed: {(Q, h, i, j)}")
        total += weight
        word.append((Q[i], Q[j], weight))
    return total, tuple(word)


def audit_finite_seven_dividing_guard_pairs(
    cache: dict[tuple[int, int, int], F],
) -> tuple[int, int, tuple[object, ...]]:
    pair_count = 0
    scalar_survivors = 0
    worst = None
    for Q in combinations(range(1, MAX_CORE_SPEED + 1), 7):
        if not divisor_complete(Q) or not hereditarily_primitive(Q):
            continue
        for h in range(7, (8 * max(Q) - 1) // 3 + 1, 14):
            if 6 * h >= 16 * max(Q):
                continue
            pair_count += 1
            intersections = tuple(mixed_overlap(h, q) for q in Q)
            deficit = F(2, 7) - sum(intersections, F(0))
            if deficit >= 0:
                scalar_survivors += 1
            tree, edges = seven_dividing_guard_tree(Q, h, cache)
            margin = tree - deficit
            require(margin > 0, f"finite terminal 7|h branch survived: {(Q, h)}")
            row = (margin, Q, h, tree, deficit, edges)
            if worst is None or row < worst:
                worst = row
    require(worst is not None, "no finite 7|h terminal pair was audited")
    return pair_count, scalar_survivors, worst


def audit_seven_dividing_guard_symbolic_margin() -> tuple[tuple[F, ...], F]:
    floors = tuple(
        F(30, 343) - F(5, 294) * (F(m) - F(6, 7))
        for m in range(1, 6)
    )
    require(all(value > 0 for value in floors), "7|h symbolic floor lost positivity")
    require(floors[-1] == F(5, 294), "7|h worst symbolic margin changed")
    require(min(floors) == F(5, 294), "7|h symbolic minimum changed")
    return floors, min(floors)


def maximum_global_tree(speeds: tuple[int, ...]) -> F:
    used = {0}
    total = F(0)
    while len(used) < len(speeds):
        best = None
        for i in used:
            for j in range(len(speeds)):
                if j in used:
                    continue
                candidate = (global_pair(speeds[i], speeds[j]), i, j)
                if best is None or candidate > best:
                    best = candidate
        require(best is not None, "global Prim search lost a vertex")
        total += best[0]
        used.add(best[2])
    return total


def audit_nonseven_guard_five_high(
    cache: dict[tuple[int, int, int], F],
) -> tuple[int, int, int, tuple[F, tuple[int, ...]]]:
    """Audit the 7-not-dividing-h, five 7-divisible speeds branch."""
    guards = (1, 3, 5, 9, 11, 13)
    high_pool = tuple(7 * k for k in range(1, 11))
    edge_checks = 0
    for h in guards:
        require(h % 7 != 0, f"nonseven guard bank contains 7|h: {h}")
        for p, q in combinations(high_pool, 2):
            epsilon_p = mixed_overlap(h, p) - F(2, 49)
            epsilon_q = mixed_overlap(h, q) - F(2, 49)
            require(epsilon_p == 0 and epsilon_q == 0, "high epsilon did not vanish")
            weight = restricted_pair(h, p, q, cache)
            require(
                genuine_residual(h, p, q, weight) == 0,
                f"five-high genuine channel survived: {(h, p, q)}",
            )
            require(
                weight == F(5, 7) * global_pair(p, q),
                f"five-high restricted edge failed: {(h, p, q)}",
            )
            edge_checks += 1
    require(edge_checks == 270, f"five-high edge count changed: {edge_checks}")

    five_set_checks = 0
    worst_tree = None
    for high in combinations(high_pool, 5):
        total_pair_mass = sum(
            (global_pair(p, q) for p, q in combinations(high, 2)), F(0)
        )
        require(total_pair_mass >= F(44, 273), f"THM-1234 replay failed: {high}")
        tree = maximum_global_tree(high)
        require(tree >= F(88, 1365), f"Cayley tree average failed: {high}")
        restricted_tree = F(5, 7) * tree
        require(restricted_tree >= F(88, 1911), f"restricted high tree failed: {high}")
        row = (restricted_tree, high)
        if worst_tree is None or row < worst_tree:
            worst_tree = row
        five_set_checks += 1
    require(five_set_checks == 252, f"five-high set count changed: {five_set_checks}")
    require(worst_tree is not None, "empty five-high set bank")

    finite_terminal_pairs = 0
    for Q in combinations(range(1, MAX_CORE_SPEED + 1), 7):
        if not divisor_complete(Q) or not hereditarily_primitive(Q):
            continue
        if sum(q % 7 == 0 for q in Q) != 5:
            continue
        for h in range(1, (8 * max(Q) - 1) // 3 + 1, 2):
            if h % 7 != 0 and 6 * h < 16 * max(Q):
                finite_terminal_pairs += 1
    require(finite_terminal_pairs == 0, "unexpected bounded five-high terminal pair")
    return edge_checks, five_set_checks, finite_terminal_pairs, worst_tree


def audit_five_high_symbolic_margin() -> tuple[F, F, F, F]:
    pair_floor = F(44, 273)
    global_tree_floor = F(2, 5) * pair_floor
    restricted_tree_floor = F(5, 7) * global_tree_floor
    deficit_ceiling = 2 * F(5, 294)
    margin = restricted_tree_floor - deficit_ceiling
    require(global_tree_floor == F(88, 1365), "five-high Cayley constant changed")
    require(restricted_tree_floor == F(88, 1911), "five-high restricted constant changed")
    require(deficit_ceiling == F(5, 147), "five-high deficit ceiling changed")
    require(margin == F(23, 1911), "five-high final margin changed")
    return global_tree_floor, restricted_tree_floor, deficit_ceiling, margin


def audit_lacunary_constant() -> tuple[F, F]:
    weighted_floor = (
        F(1, 42) + F(6, 7) * (F(2, 77) + 5 * F(1, 35)) - F(8, 49)
    )
    require(weighted_floor == F(17, 3234), "weighted lacunary floor changed")
    threshold = 3 * weighted_floor
    require(threshold == F(17, 1078), "lacunary ratio threshold changed")
    return weighted_floor, threshold


def main() -> None:
    print("THM-2086 RELATIVE HUNTER FOURIER/LACUNARY REFEREE")
    cache: dict[tuple[int, int, int], F] = {}

    pair_count, scalar_survivors, bounded_edge_count, worst = (
        audit_bounded_terminal_bank(cache)
    )
    fourier_error = audit_fourier_truncations(worst)
    spectrum_checks, low_rows = audit_odd_reduced_spectrum()
    bv_checks, sharpest_bv = audit_bv_mixing(cache)
    seven_cross_checks = audit_seven_dividing_guard_cross_edges(cache)
    seven_pair_count, seven_scalar_survivors, seven_worst = (
        audit_finite_seven_dividing_guard_pairs(cache)
    )
    seven_floors, seven_worst_floor = audit_seven_dividing_guard_symbolic_margin()
    five_high_edges, five_high_sets, five_high_finite, five_high_worst = (
        audit_nonseven_guard_five_high(cache)
    )
    five_high_constants = audit_five_high_symbolic_margin()
    weighted_floor, threshold = audit_lacunary_constant()

    print(f"bounded allowed core/guard pairs: {pair_count}")
    print(f"bounded scalar survivors: {scalar_survivors}")
    print(f"bounded distinct triple edges: {bounded_edge_count}")
    print(f"minimum bounded margin: {worst[0]}")
    print(f"hostile core/guard: {worst[1]}, {worst[2]}")
    print(f"hostile bulk term: {worst[5]}")
    print(f"hostile degree-epsilon term: {worst[6]}")
    print(f"hostile genuine rank-two term: {worst[7]}")
    print(f"hostile degree word: {worst[8]}")
    print("hostile edge/channel word:")
    for edge in worst[9]:
        print(f"  {edge}")
    print(f"cutoff-640 maximum Fourier error: {fourier_error:.6e}")
    print(f"odd-a coprime ab<=20 checks: {spectrum_checks}")
    print(f"odd-a rows at or below 1/35: {low_rows}")
    print(f"exact BV mixing checks: {bv_checks}")
    print(
        "largest normalized BV error: "
        f"{sharpest_bv[0]} at h={sharpest_bv[1]}, q={sharpest_bv[2]}, "
        f"B={sharpest_bv[3]}"
    )
    print(f"7|h broad exact cross-edge checks: {seven_cross_checks}")
    print(f"7|h bounded terminal pairs: {seven_pair_count}")
    print(f"7|h bounded scalar survivors: {seven_scalar_survivors}")
    print(
        "7|h minimum bounded bipartite margin: "
        f"{seven_worst[0]} at Q={seven_worst[1]}, h={seven_worst[2]}"
    )
    print(f"7|h symbolic margins m=1..5: {seven_floors}")
    print(f"7|h worst symbolic margin: {seven_worst_floor}")
    print(f"7nmid h five-high exact edge checks: {five_high_edges}")
    print(f"7nmid h five-high five-set checks: {five_high_sets}")
    print(f"7nmid h bounded five-high terminal pairs: {five_high_finite}")
    print(
        "7nmid h five-high minimum sampled restricted tree: "
        f"{five_high_worst[0]} at {five_high_worst[1]}"
    )
    print(
        "7nmid h five-high constants "
        "(global tree, restricted tree, deficit, margin): "
        f"{five_high_constants}"
    )
    print(f"lacunary weighted floor: {weighted_floor}")
    print(f"lacunary side-mass threshold: {threshold}")
    print("PASS")


if __name__ == "__main__":
    main()

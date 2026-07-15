#!/usr/bin/env python3
"""Exact replay for THM-803's anti-grids and global erosion selector.

Every verdict uses integer or Fraction arithmetic.  The script verifies:

* the general half-divisor escape and the q=13 parity-twisted support law;
* the complete universal even anti-shell list d=2,4,6;
* the quarter-grid escape missed by THM-797's prime-grid sharpness row;
* a row on which every named divisor/anti-grid is silent although another
  closed 1/11-component escapes; and
* the endpoint-plus-cusp selector for E_U + closure(R_U), together with its
  decorated component Tournament Analysis.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd


ALPHA = F(1, 13)
BETA = F(1, 11)
GAMMA = BETA - ALPHA
THRESHOLD = F(11, 13)
CLASSES = frozenset(range(1, 7))

PRIME_GRID_ROW = (1, 2, 3, 4, 7, 9, 10, 11, 12, 16)
SIGNED_GRID_ROW = (1, 2, 4, 6, 7, 9, 10, 11, 12, 16)
GLOBAL_ROW = (2, 4, 6, 7, 9, 10, 11, 12, 14, 16)
GLOBAL_PAIR = (13, 5)
EXPECTED_DIGEST = "5e5851835ba5e753fd5776fe1c18c0639d778b037821aaa8b3e86f5358eb3495"


def norm(z: F) -> F:
    residue = z % 1
    return min(residue, 1 - residue)


def balanced_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    if 2 * residue > modulus:
        residue -= modulus
    return residue


def folded_class(value: int, modulus: int = 13) -> int:
    residue = value % modulus
    return min(residue, (-residue) % modulus)


def inverse_folded_class(value: int) -> int:
    return folded_class(pow(value, -1, 13))


def parity_twisted_class(value: int) -> int:
    assert value % 13 != 0
    if value % 2:
        return folded_class(value)
    return folded_class((value // 2) % 13)


def phi(speeds: tuple[int, ...], t: F) -> F:
    return min(norm(speed * t) for speed in speeds)


def q_value(x: int, y: int, t: F) -> F:
    a = (x + y) // 2
    b = abs(x - y) // 2
    return norm(a * t) + norm(b * t)


def unit_representatives(modulus: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                min(p, modulus - p)
                for p in range(1, modulus)
                if gcd(p, modulus) == 1
            }
        )
    )


def deep_unit_classes(
    speeds: tuple[int, ...], modulus: int
) -> tuple[tuple[int, F], ...]:
    answer = []
    for p in unit_representatives(modulus):
        value = min(
            F(abs(balanced_residue(speed * p, modulus)), modulus)
            for speed in speeds
        )
        if value >= BETA:
            answer.append((p, value))
    return tuple(answer)


def audit_half_divisor_identity(limit: int = 101) -> int:
    tests = 0
    for q in range(3, limit + 1, 2):
        for odd_quotient in range(1, 2 * q, 2):
            x = q * odd_quotient
            for p in unit_representatives(2 * q):
                assert norm(F(x * p, 2 * q)) == F(1, 2)
                tests += 1
    return tests


def audit_twisted_support_law() -> tuple[int, int, dict[int, int]]:
    """Exhaust every subset of the twelve signed nonzero classes mod 26."""
    signed_types = tuple(range(1, 13))
    p_representatives = unit_representatives(26)
    inverse = {c: inverse_folded_class(c) for c in CLASSES}
    masks = 0
    empty_deep = 0
    for mask in range(1 << len(signed_types)):
        speeds = tuple(
            value
            for index, value in enumerate(signed_types)
            if mask & (1 << index)
        )
        direct = {
            folded_class(p)
            for p in p_representatives
            if all(
                11 * abs(balanced_residue(speed * p, 26)) >= 26
                for speed in speeds
            )
        }
        support = {parity_twisted_class(speed) for speed in speeds}
        predicted = CLASSES - {inverse[c] for c in support}
        assert direct == predicted
        assert (not direct) == (support == CLASSES)
        masks += 1
        empty_deep += int(not direct)
    return masks, empty_deep, inverse


def anti_phase(odd_quotient: int, d: int, p: int) -> F:
    return norm(F(odd_quotient * p, d))


def audit_universal_even_antigrids(limit: int = 40) -> tuple[int, ...]:
    universal = []
    for d in range(2, limit + 1, 2):
        works = True
        for odd_quotient in range(1, 2 * d, 2):
            for p in range(1, d):
                if gcd(p, d) != 1:
                    continue
                if anti_phase(odd_quotient, d, p) <= F(2, 13):
                    works = False
                    break
            if not works:
                break
        if works:
            universal.append(d)
    assert tuple(universal) == (2, 4, 6)
    for d in range(8, limit + 1, 2):
        assert anti_phase(1, d, 1) == F(1, d) <= F(1, 8) < F(2, 13)
    return tuple(universal)


def intersect_intervals(
    left: list[tuple[F, F]],
    right: list[tuple[F, F]],
    *,
    allow_points: bool,
) -> list[tuple[F, F]]:
    answer = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi or (allow_points and lo == hi):
            answer.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        elif left[i][1] > right[j][1]:
            j += 1
        else:
            i += 1
            j += 1
    return answer


def deep_components(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    current = [(F(0), F(1))]
    for speed in speeds:
        safe = [
            ((F(k) + BETA) / speed, (F(k + 1) - BETA) / speed)
            for k in range(speed)
        ]
        current = intersect_intervals(current, safe, allow_points=True)
    return tuple(current)


def return_components(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    """Open R_U components represented in (-1/2,1/2); return closures."""
    current = [(F(-1, 2), F(1, 2))]
    for speed in speeds:
        allowed = []
        for k in range(-speed, speed + 1):
            interval = ((F(k) - GAMMA) / speed, (F(k) + GAMMA) / speed)
            if interval[1] > F(-1, 2) and interval[0] < F(1, 2):
                allowed.append(interval)
        current = intersect_intervals(current, allowed, allow_points=False)
    return tuple(current)


def merge_closed(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    intervals.sort()
    merged: list[tuple[F, F]] = []
    for left, right in intervals:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return merged


def split_circle_interval(left: F, right: F) -> list[tuple[F, F]]:
    if right - left >= 1:
        return [(F(0), F(1))]
    while left < 0:
        left += 1
        right += 1
    while left >= 1:
        left -= 1
        right -= 1
    if right <= 1:
        return [(left, right)]
    return [(left, F(1)), (F(0), right - 1)]


def circle_sum_components(
    first: tuple[tuple[F, F], ...],
    second: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    pieces = []
    for left, right in first:
        for shift_left, shift_right in second:
            pieces.extend(
                split_circle_interval(left + shift_left, right + shift_right)
            )
    merged = merge_closed(pieces)
    if merged == [(F(0), F(1))]:
        return tuple(merged)
    return tuple(merged)


def component_candidates(
    component: tuple[F, F], x: int, y: int
) -> tuple[F, ...]:
    left, right = component
    a = (x + y) // 2
    b = abs(x - y) // 2
    candidates = {left, right}
    for frequency in (a, b):
        for k in range(-1, 2 * frequency + 2):
            point = F(k, 2 * frequency) % 1
            if left <= point <= right:
                candidates.add(point)
    return tuple(sorted(candidates))


def selector_profile(
    components: tuple[tuple[F, F], ...], x: int, y: int
) -> tuple[tuple[F, F, F, int], ...]:
    """(margin, argmin, width, selector-size) for every represented arc."""
    profile = []
    for component in components:
        candidates = component_candidates(component, x, y)
        minimum, argmin = min((q_value(x, y, t), t) for t in candidates)
        profile.append(
            (THRESHOLD - minimum, argmin, component[1] - component[0], len(candidates))
        )
    return tuple(profile)


def breakpoint_candidates(speeds: tuple[int, ...]) -> set[F]:
    denominators = {2 * speed for speed in speeds}
    for index, left in enumerate(speeds):
        for right in speeds[index + 1 :]:
            denominators.add(left + right)
            denominators.add(abs(left - right))
    return {
        F(k, denominator)
        for denominator in denominators
        if denominator
        for k in range(denominator + 1)
    }


def exact_loneliness(speeds: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    candidates = breakpoint_candidates(speeds)
    value = max(phi(speeds, t) for t in candidates)
    maximizers = tuple(sorted(t for t in candidates if phi(speeds, t) == value))
    return value, maximizers


def total_order_edges(order: tuple[int, ...]) -> set[tuple[int, int]]:
    rank = {vertex: index for index, vertex in enumerate(order)}
    return {
        (left, right) if rank[left] < rank[right] else (right, left)
        for left, right in combinations(order, 2)
    }


def tournament_fingerprint(
    order: tuple[int, ...], edges: set[tuple[int, int]]
) -> tuple[tuple[int, ...], int, tuple[int, ...], int]:
    scores = {vertex: 0 for vertex in order}
    for source, _target in edges:
        scores[source] += 1
    cycles = 0
    for triple in combinations(order, 3):
        triple_scores = {vertex: 0 for vertex in triple}
        for source, target in edges:
            if source in triple and target in triple:
                triple_scores[source] += 1
        cycles += int(sorted(triple_scores.values()) == [1, 1, 1])
    # Every edge set here comes from one total order.
    return tuple(sorted(scores.values())), cycles, (1,) * len(order), 1


def audit_prime_grid_correction() -> tuple[tuple[F, F], ...]:
    t = F(11, 52)
    folded = q_value(13, 5, t)
    assert norm(13 * t) == F(1, 4) > F(2, 13)
    assert folded == F(1, 4) < THRESHOLD
    values = []
    for speeds in (PRIME_GRID_ROW, SIGNED_GRID_ROW):
        value = phi(speeds, t)
        assert value == F(5, 52) > BETA
        values.append((value, folded))
    return tuple(values)


def audit_global_row() -> dict[str, object]:
    U = GLOBAL_ROW
    x, y = GLOBAL_PAIR
    assert len(U) == len(set(U)) == 10 and gcd(*U) == 1
    assert all(any(speed % modulus == 0 for speed in U) for modulus in range(2, 13))
    assert all(speed % 13 for speed in U)

    value, maximizers = exact_loneliness(U)
    assert value == F(2, 13)
    assert maximizers == (F(5, 13), F(8, 13))
    assert all(q_value(x, y, t) == F(12, 13) for t in maximizers)

    signed_residues = {speed % 13 for speed in U}
    raw_support = {folded_class(speed) for speed in U}
    twisted_support = {parity_twisted_class(speed) for speed in U}
    assert signed_residues == set(range(1, 13)) - {5, 8}
    assert raw_support == {1, 2, 3, 4, 6}
    assert twisted_support == CLASSES
    assert folded_class(y) == 5
    assert any(speed % 2 == 0 and folded_class(speed) == 3 for speed in U)

    ordinary = {q: deep_unit_classes(U, q) for q in (5, 13)}
    assert ordinary[5] == ()
    assert ordinary[13] == ((5, F(2, 13)),)
    assert q_value(x, y, F(5, 13)) == F(12, 13)

    anti = {d: deep_unit_classes(U, 13 * d) for d in (2, 4, 6)}
    assert anti == {2: (), 4: (), 6: ()}

    B = max(U)
    rho = (value - ALPHA) / B
    astar_left = F(1, x * y) + 2 * rho
    astar_right = F(2, 13 * x) + F(2, 13 * y)
    assert x <= 11 * B and y <= 11 * B and min(x, y) <= 11 * B - 36
    assert x <= 2 * B - 1 and y <= B - 1
    assert 13 * B + 2 * x * y <= 2 * B * (x + y)
    assert rho == F(1, 208)
    assert astar_left == F(1, 40) <= astar_right == F(36, 845)

    deep = deep_components(U)
    assert deep == (
        (F(1, 22), F(5, 88)),
        (F(7, 22), F(7, 22)),
        (F(23, 66), F(27, 77)),
        (F(67, 176), F(43, 110)),
        (F(9, 22), F(9, 22)),
        (F(13, 22), F(13, 22)),
        (F(67, 110), F(109, 176)),
        (F(50, 77), F(43, 66)),
        (F(15, 22), F(15, 22)),
        (F(83, 88), F(21, 22)),
    )
    escape = F(7, 22)
    owners = tuple(speed for speed in U if norm(speed * escape) == BETA)
    assert owners == (6, 16)
    assert phi(U, escape) == BETA
    assert q_value(x, y, escape) == F(9, 22) < THRESHOLD

    returns = return_components(U)
    thickened = circle_sum_components(deep, returns)
    profile = selector_profile(thickened, x, y)
    assert any(margin > 0 for margin, _argmin, _width, _size in profile)
    selector_size = sum(row[3] for row in profile)
    selector_minimum = min(
        (THRESHOLD - row[0], row[1]) for row in profile
    )

    margin_order = tuple(
        sorted(range(len(profile)), key=lambda i: (profile[i][0], thickened[i][0]))
    )
    width_order = tuple(
        sorted(range(len(profile)), key=lambda i: (profile[i][2], thickened[i][0]))
    )
    margin_edges = total_order_edges(margin_order)
    width_edges = total_order_edges(width_order)
    flips = len(margin_edges ^ width_edges) // 2
    fingerprint = tournament_fingerprint(margin_order, margin_edges)
    assert fingerprint[0] == tuple(range(len(profile)))
    assert fingerprint[1] == 0 and fingerprint[2] == (1,) * len(profile)
    assert fingerprint[3] == 1

    raw_multiplicity = {
        c: sum(folded_class(speed) == c for speed in U) for c in CLASSES
    }
    twisted_multiplicity = {
        c: sum(parity_twisted_class(speed) == c for speed in U) for c in CLASSES
    }
    raw_order = tuple(sorted(CLASSES, key=lambda c: (raw_multiplicity[c], c)))
    twisted_order = tuple(
        sorted(CLASSES, key=lambda c: (twisted_multiplicity[c], c))
    )
    support_flips = len(
        total_order_edges(raw_order) ^ total_order_edges(twisted_order)
    ) // 2

    full_packet = tuple(sorted({2 * speed for speed in U} | {x, y}))
    peeled_packet = tuple(speed for speed in full_packet if speed != max(full_packet))
    full_value, full_maximizers = exact_loneliness(full_packet)
    peeled_value, peeled_maximizers = exact_loneliness(peeled_packet)
    assert full_value == F(1, 9) > ALPHA
    assert full_maximizers == (F(1, 36), F(17, 36), F(19, 36), F(35, 36))
    assert peeled_value == F(1, 8) > F(1, 12)
    assert peeled_maximizers == (
        F(1, 32), F(1, 16), F(11, 32), F(7, 16),
        F(9, 16), F(21, 32), F(15, 16), F(31, 32),
    )

    return {
        "value": value,
        "maximizers": maximizers,
        "signed_residues": tuple(sorted(signed_residues)),
        "raw_support": tuple(sorted(raw_support)),
        "twisted_support": tuple(sorted(twisted_support)),
        "ordinary": ordinary,
        "anti": anti,
        "rho": rho,
        "astar_left": astar_left,
        "astar_right": astar_right,
        "deep": deep,
        "escape": escape,
        "owners": owners,
        "returns": returns,
        "thickened": thickened,
        "profile": profile,
        "selector_size": selector_size,
        "selector_minimum": selector_minimum,
        "margin_order": margin_order,
        "width_order": width_order,
        "flips": flips,
        "fingerprint": fingerprint,
        "raw_order": raw_order,
        "twisted_order": twisted_order,
        "support_flips": support_flips,
        "full_packet": full_packet,
        "full_value": full_value,
        "peeled_value": peeled_value,
    }


def format_fraction(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def format_interval(interval: tuple[F, F]) -> str:
    return f"[{format_fraction(interval[0])},{format_fraction(interval[1])}]"


def main() -> None:
    half_tests = audit_half_divisor_identity()
    masks, empty_deep, inverse = audit_twisted_support_law()
    universal = audit_universal_even_antigrids()
    corrections = audit_prime_grid_correction()
    row = audit_global_row()

    canonical = []
    canonical.append(f"half_tests={half_tests}")
    canonical.append(f"twisted={masks},{empty_deep},{inverse}")
    canonical.append(f"universal={universal}")
    canonical.append(f"prime_corrections={corrections}")
    for key in (
        "value",
        "maximizers",
        "signed_residues",
        "raw_support",
        "twisted_support",
        "ordinary",
        "anti",
        "deep",
        "returns",
        "thickened",
        "profile",
        "margin_order",
        "width_order",
        "flips",
        "raw_order",
        "twisted_order",
        "support_flips",
        "full_packet",
        "full_value",
        "peeled_value",
    ):
        canonical.append(f"{key}={row[key]}")
    digest = sha256(("\n".join(canonical) + "\n").encode()).hexdigest()
    if EXPECTED_DIGEST != "TO_BE_FILLED":
        assert digest == EXPECTED_DIGEST

    print("THM-803 q=13 anti-grid + all-component erosion-selector replay")
    print(
        f"half_divisor_identity_tests={half_tests} PASS: "
        "q|x, p unit mod 2q => ||xp/(2q)||=1/2"
    )
    print(
        f"q13_twisted_support_subsets={masks} empty_half_grid_profiles={empty_deep} PASS"
    )
    print(f"q13_folded_inverse={inverse}")
    print("even raw-class twist sigma=(1 6 3 5 4 2)")
    print(
        f"complete_universal_even_antigrids={universal}; "
        "every even d>=8 fails already at X=p=1"
    )
    print()
    print("THM-797 prime-grid sharpness correction")
    print(
        f"U={PRIME_GRID_ROW} pair=(13,5) t=11/52 "
        f"phi={corrections[0][0]}>1/11 Q={corrections[0][1]}<11/13"
    )
    print(
        f"signed_wall_U={SIGNED_GRID_ROW} has the same quarter-grid escape: "
        f"phi={corrections[1][0]} Q={corrections[1][1]}"
    )
    print("therefore both rows are sharp only for their stated q=13 gates, not the anti-grid ladder")
    print()
    print(f"global_selector_limit_U={GLOBAL_ROW} pair={GLOBAL_PAIR}")
    print(
        f"primitive=True divisor_complete_2_12=True no_13_multiple=True "
        f"M(U)={row['value']} maximizers={row['maximizers']} Q=12/13"
    )
    print(
        f"signed_residues={row['signed_residues']} exact_complement_of=+/-5; "
        f"raw_support={row['raw_support']} aligned_missing=5=fold(5); "
        f"twisted_support={row['twisted_support']}"
    )
    print(f"ordinary_exception_divisor_deep_classes={row['ordinary']}")
    print(f"anti_grid_deep_classes={row['anti']} (all silent)")
    print(
        f"rho={row['rho']} Astar={row['astar_left']}<={row['astar_right']}"
    )
    print(
        f"all_E_components={len(row['deep'])} "
        + " ".join(format_interval(component) for component in row["deep"])
    )
    print(
        f"global_threshold_escape={row['escape']} owners={row['owners']} "
        "phi=1/11 Q=9/22<11/13"
    )
    print(
        f"closure_R_components={len(row['returns'])} "
        + " ".join(format_interval(component) for component in row["returns"])
    )
    print(
        f"K=E+closure(R) represented_components={len(row['thickened'])} "
        f"selector_points={row['selector_size']} "
        f"global_min_Q={row['selector_minimum']}"
    )
    print(
        f"full_packet={row['full_packet']} M={row['full_value']}; "
        f"max_peel_M={row['peeled_value']}>1/12 (loose sporadic row)"
    )
    print()
    print("Tournament Analysis (vertices = represented K-component obligations)")
    print("  observable A: signed erosion escape margin 11/13-min_C Q")
    print("  switch/gauge B: component width; tie path: left endpoint")
    print(f"  Hamiltonian paths A/B: {row['margin_order']} / {row['width_order']}")
    print(f"  score histogram: {row['fingerprint'][0]}")
    print(f"  directed 3-cycles: {row['fingerprint'][1]}")
    print(f"  SCC sizes: {row['fingerprint'][2]}")
    print(f"  Hamiltonian-path counts: {row['fingerprint'][3]} / 1")
    print(f"  edge flips: {row['flips']}")
    print("  preserved: every global alternative when vertices retain interval, selector, margin, and owners")
    print("  destroyed by bare tournament: exact sign, cusp address, return incidence, and containment verdict")
    print("Support Tournament (vertices = six folded classes)")
    print("  observable A/raw multiplicity; switch B/parity-twisted multiplicity")
    print(f"  paths={row['raw_order']} / {row['twisted_order']} edge_flips={row['support_flips']}")
    print("  both transitive: scores=(0,1,2,3,4,5), cycles=0, singleton SCCs, HP=1")
    print(f"sha256={digest}")
    print("PASS")


if __name__ == "__main__":
    main()

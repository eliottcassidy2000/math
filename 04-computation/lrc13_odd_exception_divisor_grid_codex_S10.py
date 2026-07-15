#!/usr/bin/env python3
"""Exact replay for the odd-grid obstruction inside THM-789's residual.

At an odd rational grid p/q, nearest-integer parity can be read directly
from balanced residues.  When q divides one odd exception this becomes an
exception-divisor sieve.  The mandatory q=13 grid in THM-772 then gives a
six-vertex folded-class obstruction; its one-sided walls recover the signed
residue multiset and a sharp metric pin.  All arithmetic below is integral or
uses Fraction; there is no floating-point or sampled-circle step.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import gcd


ALPHA = F(1, 13)
BETA = F(1, 11)
GAMMA = BETA - ALPHA
U0 = (1, 2, 3, 4, 7, 9, 10, 11, 12, 16)
X0, Y0 = 13, 5
U1 = (1, 2, 4, 6, 7, 9, 10, 11, 12, 16)
X1, Y1 = 13, 5


def norm(z: F) -> F:
    residue = z % 1
    return min(residue, 1 - residue)


def phi(speeds: tuple[int, ...], t: F) -> F:
    return min(norm(speed * t) for speed in speeds)


def folded_q(x: int, y: int, t: F) -> F:
    a = (x + y) // 2
    b = abs(x - y) // 2
    return norm(a * t) + norm(b * t)


def balanced_residue(value: int, modulus: int) -> int:
    """The unique residue in [-(q-1)/2,(q-1)/2] for odd q."""
    assert modulus >= 3 and modulus % 2 == 1
    residue = value % modulus
    if 2 * residue > modulus:
        residue -= modulus
    return residue


def folded_class(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, (-residue) % modulus)


def nearest_integer_data(speed: int, p: int, q: int) -> tuple[int, int]:
    residue = balanced_residue(speed * p, q)
    numerator = speed * p - residue
    assert numerator % q == 0
    return numerator // q, residue


def direct_sheet_predicate(p: int, q: int, x: int, y: int) -> bool:
    nx, rx = nearest_integer_data(x, p, q)
    ny, ry = nearest_integer_data(y, p, q)
    eligible = 13 * abs(rx) <= 2 * q and 13 * abs(ry) <= 2 * q
    return eligible and (nx - ny) % 2 == 1


def residue_sheet_predicate(p: int, q: int, x: int, y: int) -> bool:
    """Odd-grid form: small balanced residues of opposite parity."""
    rx = balanced_residue(x * p, q)
    ry = balanced_residue(y * p, q)
    shell = (2 * q) // 13
    return (
        abs(rx) <= shell
        and abs(ry) <= shell
        and (rx - ry) % 2 == 1
    )


def folded_diamond_predicate(p: int, q: int, x: int, y: int) -> bool:
    a = (x + y) // 2
    b = abs(x - y) // 2
    ra = abs(balanced_residue(a * p, q))
    rb = abs(balanced_residue(b * p, q))
    return 13 * (ra + rb) >= 11 * q


def deep_unit_predicate(
    speeds: tuple[int, ...], p: int, q: int
) -> bool:
    threshold = (q + 10) // 11
    return all(
        abs(balanced_residue(speed * p, q)) >= threshold
        for speed in speeds
    )


def audit_general_odd_grid(limit: int = 101) -> tuple[int, int, int]:
    """Exhaust all odd speed classes mod 2q and all unit numerators."""
    modulus_count = 0
    unit_count = 0
    tests = 0
    for q in range(3, limit + 1, 2):
        modulus_count += 1
        odd_speed_classes = range(1, 2 * q, 2)
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            unit_count += 1
            for x in odd_speed_classes:
                for y in odd_speed_classes:
                    direct = direct_sheet_predicate(p, q, x, y)
                    residue = residue_sheet_predicate(p, q, x, y)
                    diamond = folded_diamond_predicate(p, q, x, y)
                    assert direct == residue == diamond
                    tests += 1
    return modulus_count, unit_count, tests


def divisor_shell_predicate(p: int, q: int, y: int) -> bool:
    """Acceptance shell when the other odd exception is x=q."""
    if y % q == 0:
        return False
    residue = balanced_residue(y * p, q)
    return 1 <= abs(residue) <= (2 * q) // 13 and residue % 2 != 0


def audit_divisor_specialization(limit: int = 101) -> int:
    tests = 0
    for q in range(3, limit + 1, 2):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            for y in range(1, 2 * q, 2):
                assert direct_sheet_predicate(p, q, q, y) == (
                    divisor_shell_predicate(p, q, y)
                )
                tests += 1
    return tests


def audit_q13_corollary() -> tuple[int, int, int, dict[int, int]]:
    """Exhaust all folded-support subsets and odd y classes mod 26."""
    classes = set(range(1, 7))
    inverse = {
        c: folded_class(pow(c, -1, 13), 13)
        for c in classes
    }
    assertions = 0
    off_divisor_successes = 0
    double_divisor_successes = 0
    for mask in range(1 << 6):
        support = {c for c in classes if mask & (1 << (c - 1))}
        deep_classes = {
            p
            for p in classes
            if all(folded_class(u * p, 13) != 1 for u in support)
        }
        assert deep_classes == classes - {inverse[u] for u in support}
        for y in range(1, 26, 2):
            accepted = {
                p
                for p in classes
                if direct_sheet_predicate(p, 13, 13, y)
            }
            actual = deep_classes <= accepted
            if y % 13 == 0:
                expected = support == classes
                double_divisor_successes += int(actual)
            else:
                expected = (classes - support) <= {
                    folded_class(y, 13)
                }
                off_divisor_successes += int(actual)
            assert actual == expected
            assertions += 1
    return (
        assertions,
        off_divisor_successes,
        double_divisor_successes,
        inverse,
    )


def audit_q13_signed_wall_census() -> dict[str, object]:
    """Exhaust signed residue multisets after the folded support gate.

    A rejected q=13 grid wall can be locally tight only when its owner class
    occurs with both residue signs.  There are C(21,10)=352716 multisets of
    ten nonzero residues modulo 13; this census applies that exact wall rule
    after the folded support gate for every possible missing class.
    """
    classes = set(range(1, 7))
    profiles = 0
    off_divisor_gate_pairs = 0
    off_divisor_wall_survivors: list[tuple[int, tuple[int, ...]]] = []
    double_divisor_gate_profiles = 0
    double_divisor_wall_survivors = 0
    for residues in combinations_with_replacement(range(1, 13), 10):
        profiles += 1
        counts = {residue: residues.count(residue) for residue in range(1, 13)}
        support = {folded_class(residue, 13) for residue in residues}

        if support == classes:
            double_divisor_gate_profiles += 1
            double_wall = all(counts[c] and counts[13 - c] for c in classes)
            double_divisor_wall_survivors += int(double_wall)

        for missing in classes:
            gate = classes - support <= {missing}
            if not gate:
                continue
            off_divisor_gate_pairs += 1
            wall = all(
                counts[c] and counts[13 - c]
                for c in support
                if c != missing
            )
            if wall:
                off_divisor_wall_survivors.append((missing, residues))

    expected = [
        (
            missing,
            tuple(
                residue
                for residue in range(1, 13)
                if folded_class(residue, 13) != missing
            ),
        )
        for missing in range(1, 7)
    ]
    off_divisor_wall_survivors.sort()
    assert off_divisor_wall_survivors == expected
    assert double_divisor_wall_survivors == 0
    digest = sha256(repr(tuple(off_divisor_wall_survivors)).encode()).hexdigest()
    return {
        "profiles": profiles,
        "off_gate_pairs": off_divisor_gate_pairs,
        "off_survivors": tuple(off_divisor_wall_survivors),
        "double_gate_profiles": double_divisor_gate_profiles,
        "double_survivors": double_divisor_wall_survivors,
        "digest": digest,
    }


def audit_narrow_endpoint_owner(limit: int = 200) -> int:
    """Verify the equality case of the aligned metric pin.

    Normalize the accepted grid by yp=+1.  At displacement 1/(13B), among
    starting residues +/-2,...,+/-6 and speeds u<=B, clearance can reach
    1/13 only for the max speed B starting at residue -2.
    """
    tests = 0
    for B in range(1, limit + 1):
        for speed in range(1, B + 1):
            for residue in (-6, -5, -4, -3, -2, 2, 3, 4, 5, 6):
                clearance = norm(F(residue, 13) + F(speed, 13 * B))
                assert clearance >= ALPHA
                assert (clearance == ALPHA) == (speed == B and residue == -2)
                tests += 1
    return tests


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
    maximizers = tuple(
        sorted(t for t in candidates if phi(speeds, t) == value)
    )
    return value, maximizers


def intersect_closed(
    left: list[tuple[F, F]], right: list[tuple[F, F]]
) -> list[tuple[F, F]]:
    answer = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo <= hi:
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
        current = intersect_closed(current, safe)
    return tuple(current)


def component_q_minimum(
    component: tuple[F, F], x: int, y: int
) -> tuple[F, F]:
    left, right = component
    a = (x + y) // 2
    b = abs(x - y) // 2
    candidates = {left, right}
    for frequency in (a, b):
        low = int(2 * frequency * left) - 2
        high = int(2 * frequency * right) + 3
        for k in range(low, high + 1):
            cusp = F(k, 2 * frequency)
            if left <= cusp <= right:
                candidates.add(cusp)
    return min((folded_q(x, y, t), t) for t in candidates)


def total_order_edges(order: list[int]) -> set[tuple[int, int]]:
    rank = {vertex: index for index, vertex in enumerate(order)}
    edges = set()
    for left, right in combinations(order, 2):
        if rank[left] < rank[right]:
            edges.add((left, right))
        else:
            edges.add((right, left))
    return edges


def tournament_fingerprint(
    vertices: tuple[int, ...], edges: set[tuple[int, int]]
) -> tuple[tuple[int, ...], int]:
    scores = {vertex: 0 for vertex in vertices}
    for source, _target in edges:
        scores[source] += 1
    cycles = 0
    for triple in combinations(vertices, 3):
        induced = {
            edge for edge in edges if edge[0] in triple and edge[1] in triple
        }
        triple_scores = {vertex: 0 for vertex in triple}
        for source, _target in induced:
            triple_scores[source] += 1
        cycles += int(sorted(triple_scores.values()) == [1, 1, 1])
    return tuple(sorted(scores.values())), cycles


def audit_sharp_survivor() -> dict[str, object]:
    assert len(U0) == 10 and len(set(U0)) == 10 and gcd(*U0) == 1
    assert all(any(speed % modulus == 0 for speed in U0) for modulus in range(2, 13))
    assert all(speed % 13 != 0 for speed in U0)

    B = max(U0)
    value, maximizers = exact_loneliness(U0)
    assert value == F(2, 13)
    assert maximizers == (F(5, 13), F(8, 13))
    assert X0 % 2 == Y0 % 2 == 1 and X0 != Y0
    assert X0 <= 11 * B and Y0 <= 11 * B
    assert min(X0, Y0) <= 11 * B - 36

    rho = (value - ALPHA) / B
    astar_left = F(1, X0 * Y0) + 2 * rho
    astar_right = F(2, 13 * X0) + F(2, 13 * Y0)
    assert rho == F(1, 208)
    assert astar_left == F(1, 40) <= astar_right == F(36, 845)

    support_multiplicity = {
        c: sum(folded_class(speed, 13) == c for speed in U0)
        for c in range(1, 7)
    }
    assert tuple(support_multiplicity.values()) == (2, 2, 3, 2, 0, 1)
    support = {c for c, count in support_multiplicity.items() if count}
    assert support == {1, 2, 3, 4, 6}
    assert folded_class(Y0, 13) == 5

    deep_thirteen = tuple(
        p for p in range(1, 13) if deep_unit_predicate(U0, p, 13)
    )
    assert deep_thirteen == (5, 8)
    for p in deep_thirteen:
        t = F(p, 13)
        assert phi(U0, t) == F(2, 13)
        assert folded_q(X0, Y0, t) == F(12, 13)
        assert direct_sheet_predicate(p, 13, X0, Y0)

    escape = F(7, 33)
    assert phi(U0, escape) == BETA
    assert folded_q(X0, Y0, escape) == F(8, 33) < F(11, 13)
    assert max(norm(speed * F(0)) for speed in U0) == 0 < GAMMA

    full_packet = tuple(sorted({2 * speed for speed in U0} | {X0, Y0}))
    assert full_packet == (2, 4, 5, 6, 8, 13, 14, 18, 20, 22, 24, 32)
    full_value, full_maximizers = exact_loneliness(full_packet)
    peeled_packet = tuple(speed for speed in full_packet if speed != 32)
    peeled_value, peeled_maximizers = exact_loneliness(peeled_packet)
    assert full_value == F(2, 19) > ALPHA
    assert full_maximizers == (F(2, 19), F(8, 19), F(11, 19), F(17, 19))
    assert peeled_value == F(1, 8) > F(1, 12)
    assert peeled_maximizers == (F(1, 16), F(7, 16), F(9, 16), F(15, 16))

    components = deep_components(U0)
    assert components == (
        (F(12, 77), F(7, 44)),
        (F(23, 110), F(7, 33)),
        (F(67, 176), F(43, 110)),
        (F(9, 22), F(9, 22)),
        (F(13, 22), F(13, 22)),
        (F(67, 110), F(109, 176)),
        (F(26, 33), F(87, 110)),
        (F(37, 44), F(65, 77)),
    )
    endpoint_owners = tuple(
        (
            tuple(speed for speed in U0 if norm(speed * left) == BETA),
            tuple(speed for speed in U0 if norm(speed * right) == BETA),
        )
        for left, right in components
    )
    assert endpoint_owners == (
        ((7,), (12,)),
        ((10,), (9,)),
        ((16,), (10,)),
        ((10, 12), (10, 12)),
        ((10, 12), (10, 12)),
        ((10,), (16,)),
        ((9,), (10,)),
        ((12,), (7,)),
    )
    reflection_orbits = ((0, 7), (1, 6), (2, 5), (3, 4))
    for left_index, right_index in reflection_orbits:
        left_component = components[left_index]
        right_component = components[right_index]
        assert right_component == (
            1 - left_component[1],
            1 - left_component[0],
        )
    component_minima = tuple(
        component_q_minimum(components[index], X0, Y0)
        for index, _mate in reflection_orbits
    )
    assert component_minima == (
        (F(60, 77), F(12, 77)),
        (F(8, 33), F(7, 33)),
        (F(159, 176), F(67, 176)),
        (F(15, 22), F(9, 22)),
    )
    component_margins = {
        index: F(11, 13) - component_minima[index][0]
        for index in range(4)
    }
    component_widths = {
        index: components[orbit[0]][1] - components[orbit[0]][0]
        for index, orbit in enumerate(reflection_orbits)
    }
    assert component_margins == {
        0: F(67, 1001),
        1: F(259, 429),
        2: -F(131, 2288),
        3: F(47, 286),
    }
    assert component_widths == {
        0: F(1, 308),
        1: F(1, 330),
        2: F(9, 880),
        3: F(0),
    }
    component_margin_order = sorted(
        range(4), key=lambda index: (component_margins[index], index)
    )
    component_width_order = sorted(
        range(4), key=lambda index: (component_widths[index], index)
    )
    margin_edges = total_order_edges(component_margin_order)
    width_edges = total_order_edges(component_width_order)
    component_flips = len(margin_edges ^ width_edges) // 2
    component_scores, component_cycles = tournament_fingerprint(
        tuple(range(4)), margin_edges
    )
    assert component_margin_order == [2, 0, 3, 1]
    assert component_width_order == [3, 1, 0, 2]
    assert component_flips == 5
    assert component_scores == (0, 1, 2, 3) and component_cycles == 0

    vertices = tuple(range(1, 7))
    inverse = {
        c: folded_class(pow(c, -1, 13), 13)
        for c in vertices
    }
    multiplicity_order = sorted(
        vertices, key=lambda c: (support_multiplicity[c], c)
    )
    inverted_order = sorted(
        vertices, key=lambda c: (support_multiplicity[inverse[c]], c)
    )
    multiplicity_edges = total_order_edges(multiplicity_order)
    inverted_edges = total_order_edges(inverted_order)
    edge_flips = len(multiplicity_edges ^ inverted_edges) // 2
    scores, cycles = tournament_fingerprint(vertices, multiplicity_edges)
    assert multiplicity_order == [5, 6, 1, 2, 4, 3]
    assert inverted_order == [5, 2, 1, 3, 6, 4]
    assert edge_flips == 5
    assert scores == (0, 1, 2, 3, 4, 5) and cycles == 0

    return {
        "B": B,
        "value": value,
        "maximizers": maximizers,
        "rho": rho,
        "astar_left": astar_left,
        "astar_right": astar_right,
        "support_multiplicity": tuple(support_multiplicity.values()),
        "deep_thirteen": deep_thirteen,
        "escape": escape,
        "escape_phi": phi(U0, escape),
        "escape_q": folded_q(X0, Y0, escape),
        "full_packet": full_packet,
        "full_value": full_value,
        "full_maximizers": full_maximizers,
        "peeled_value": peeled_value,
        "peeled_maximizers": peeled_maximizers,
        "components": components,
        "endpoint_owners": endpoint_owners,
        "component_minima": component_minima,
        "component_margins": tuple(component_margins.values()),
        "component_widths": tuple(component_widths.values()),
        "component_margin_order": tuple(component_margin_order),
        "component_width_order": tuple(component_width_order),
        "component_flips": component_flips,
        "multiplicity_order": tuple(multiplicity_order),
        "inverted_order": tuple(inverted_order),
        "edge_flips": edge_flips,
        "scores": scores,
        "cycles": cycles,
    }


def odd_divisors(value: int) -> tuple[int, ...]:
    return tuple(
        divisor
        for divisor in range(3, value + 1, 2)
        if value % divisor == 0
    )


def folded_deep_and_shell(
    speeds: tuple[int, ...], q: int, x: int, y: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    units = range(1, q)
    deep = tuple(
        sorted(
            {
                folded_class(p, q)
                for p in units
                if gcd(p, q) == 1 and deep_unit_predicate(speeds, p, q)
            }
        )
    )
    shell = tuple(
        sorted(
            {
                folded_class(p, q)
                for p in units
                if gcd(p, q) == 1 and direct_sheet_predicate(p, q, x, y)
            }
        )
    )
    return deep, shell


def audit_signed_wall_survivor() -> dict[str, object]:
    """A row passing the signed q=13 wall but escaping at q=17."""
    assert len(U1) == 10 and len(set(U1)) == 10 and gcd(*U1) == 1
    assert all(any(speed % modulus == 0 for speed in U1) for modulus in range(2, 13))
    assert all(speed % 13 != 0 for speed in U1)
    B = max(U1)

    residues = tuple(sorted(speed % 13 for speed in U1))
    expected_residues = tuple(residue for residue in range(1, 13) if residue not in (5, 8))
    assert residues == expected_residues
    support_multiplicity = tuple(
        sum(folded_class(speed, 13) == c for speed in U1)
        for c in range(1, 7)
    )
    assert support_multiplicity == (2, 2, 2, 2, 0, 2)
    for c in (1, 2, 3, 4, 6):
        assert c in residues and 13 - c in residues

    value, maximizers = exact_loneliness(U1)
    assert value == F(2, 13)
    assert maximizers == (F(5, 13), F(8, 13))
    assert X1 % 13 == 0 and Y1 % 13 != 0
    assert X1 <= 2 * B - 1 and Y1 <= B - 1
    assert 13 * B + 2 * X1 * Y1 <= 2 * B * (X1 + Y1)
    rho = (value - ALPHA) / B
    astar_left = F(1, X1 * Y1) + 2 * rho
    astar_right = F(2, 13 * X1) + F(2, 13 * Y1)
    assert rho == F(1, 208)
    assert astar_left == F(1, 40) <= astar_right == F(36, 845)

    divisors = tuple(sorted(set(odd_divisors(X1)) | set(odd_divisors(Y1))))
    assert divisors == (5, 13)
    divisor_rows = {
        q: folded_deep_and_shell(U1, q, X1, Y1)
        for q in divisors
    }
    assert divisor_rows == {
        5: ((), ()),
        13: ((5,), (5,)),
    }
    for deep, shell in divisor_rows.values():
        assert set(deep) <= set(shell)

    witness = F(6, 17)
    witness_phi = phi(U1, witness)
    witness_q = folded_q(X1, Y1, witness)
    assert witness_phi == F(2, 17) > BETA
    assert witness_q == F(10, 17) < F(11, 13)
    deep17, shell17 = folded_deep_and_shell(U1, 17, X1, Y1)
    assert folded_class(6, 17) in deep17
    assert folded_class(6, 17) not in shell17

    payload = (
        U1,
        residues,
        value,
        maximizers,
        tuple(divisor_rows.items()),
        witness,
        witness_phi,
        witness_q,
    )
    digest = sha256(repr(payload).encode()).hexdigest()
    return {
        "B": B,
        "residues": residues,
        "support_multiplicity": support_multiplicity,
        "value": value,
        "maximizers": maximizers,
        "rho": rho,
        "astar_left": astar_left,
        "astar_right": astar_right,
        "divisor_rows": divisor_rows,
        "witness": witness,
        "witness_phi": witness_phi,
        "witness_q": witness_q,
        "deep17": deep17,
        "shell17": shell17,
        "digest": digest,
    }


def main() -> None:
    moduli, units, grid_tests = audit_general_odd_grid()
    divisor_tests = audit_divisor_specialization()
    q13_assertions, off_successes, double_successes, inverse = (
        audit_q13_corollary()
    )
    wall = audit_q13_signed_wall_census()
    endpoint_tests = audit_narrow_endpoint_owner()
    sharp = audit_sharp_survivor()
    signed = audit_signed_wall_survivor()

    print("Odd exception-divisor grid / global deep-component replay")
    print(
        "general_odd_grid: q=3..101 odd "
        f"moduli={moduli} unit_numerators={units} "
        f"(p,x,y)_tests={grid_tests} PASS"
    )
    print(
        "equivalence: lifted eligibility+opposite nearest parity "
        "<=> balanced residues in the 2q/13 shell with opposite parity "
        "<=> folded diamond"
    )
    print(
        f"divisor_specialization_tests={divisor_tests} PASS; "
        "q|x accepts exactly odd nonzero y-residues in "
        "|s|<=floor(2q/13), and accepts none when q|y"
    )
    print(
        "q13_folded_support: "
        f"subset/y_assertions={q13_assertions} "
        f"off_divisor_containment_profiles={off_successes} "
        f"double_divisor_containment_profiles={double_successes} PASS"
    )
    print(f"q13_inversion={inverse}")
    print(
        "q13_corollary: E_U subset H minus R_U with exactly one "
        "13-divisible exception forces "
        "C\\S(U) subset {fold(y)}; two 13-divisible exceptions force S(U)=C"
    )
    print(
        "q13_signed_wall_census: "
        f"residue_multisets={wall['profiles']} "
        f"off_divisor_gate_pairs={wall['off_gate_pairs']} "
        f"double_divisor_gate_profiles={wall['double_gate_profiles']} PASS"
    )
    print(
        "q13_signed_wall_survivors: "
        f"off_divisor={len(wall['off_survivors'])} "
        f"double_divisor={wall['double_survivors']} "
        f"digest={wall['digest']}"
    )
    print(
        "q13_signed_wall_conclusion: full support and double-13 are empty; "
        "for each missing class c the unique residue multiset is "
        "(Z/13Z)^*\\{+/-c}"
    )
    print(
        f"narrow_endpoint_owner_tests={endpoint_tests} PASS; "
        "normalized y=B equality forces the max-speed residue -2"
    )
    print()
    print(f"sharp_survivor_U={U0}")
    print(
        f"primitive=True divisor_complete_2_12=True no_13_multiple=True "
        f"B={sharp['B']} M(U)={sharp['value']} "
        f"maximizers={sharp['maximizers']}"
    )
    print(
        f"exceptions={(X0, Y0)} bounds=(<=11B, min<=11B-36)=PASS "
        f"rho={sharp['rho']} Astar={sharp['astar_left']}"
        f"<={sharp['astar_right']}"
    )
    print(
        f"folded_support_multiplicity={sharp['support_multiplicity']} "
        "missing_class=5=fold(y)"
    )
    print(
        f"q13_deep_numerators={sharp['deep_thirteen']} "
        "phi=2/13 Q=12/13 (all trapped)"
    )
    print(
        f"erosion_escape_t={sharp['escape']} "
        f"phi={sharp['escape_phi']} Q={sharp['escape_q']}<11/13; "
        "0 in R_U, so t notin H minus R_U"
    )
    print(
        f"full_packet={sharp['full_packet']} M={sharp['full_value']} "
        f"maximizers={sharp['full_maximizers']} (>1/13, loose)"
    )
    print(
        f"max_peel_M={sharp['peeled_value']} "
        f"maximizers={sharp['peeled_maximizers']} (>1/12, sporadic inequality)"
    )
    print(f"E_U_components={len(sharp['components'])} {sharp['components']}")
    print(f"component_endpoint_owners={sharp['endpoint_owners']}")
    print(
        "reflection_orbit_(minQ,argmin)="
        f"{sharp['component_minima']} margins={sharp['component_margins']}"
    )
    print()
    print(f"signed_wall_survivor_U={U1}")
    print(
        f"signed_residues_mod13={signed['residues']} "
        f"folded_multiplicity={signed['support_multiplicity']} "
        f"B={signed['B']} M(U)={signed['value']} "
        f"maximizers={signed['maximizers']}"
    )
    print(
        f"metric_pin=(x<=2B-1,y<=B-1,13B+2xy<=2B(x+y))=PASS "
        f"rho={signed['rho']} Astar={signed['astar_left']}"
        f"<={signed['astar_right']}"
    )
    print(
        f"all_exception_divisor_grids={signed['divisor_rows']} "
        "have D_q subset A_q"
    )
    print(
        f"nondivisor_escape_t={signed['witness']} "
        f"phi={signed['witness_phi']} Q={signed['witness_q']}<11/13; "
        f"q17_deep={signed['deep17']} q17_shell={signed['shell17']}"
    )
    print(f"signed_wall_survivor_digest={signed['digest']}")
    print()
    print("Tournament Analysis (vertices = six folded mod-13 obligations)")
    print("  pairwise observable A: core residue multiplicity")
    print("  switch/gauge B: inversion c -> c^(-1)")
    print(f"  tie Hamiltonian paths: {sharp['multiplicity_order']} / {sharp['inverted_order']}")
    print(f"  score histogram: {sharp['scores']}")
    print(f"  directed 3-cycles: {sharp['cycles']} / 0")
    print("  SCC sizes: (1,1,1,1,1,1)")
    print("  Hamiltonian-path counts: 1 / 1")
    print(f"  edge flips between gauges: {sharp['edge_flips']}")
    print("  preserved by decorated carrier: signed-wall/missing-class/deep-grid/exception-shell incidence")
    print("  destroyed by bare tournament: sign occupancy, class identities, inversion, widths, and E_U subset H verdict")
    print("Component Tournament (vertices = four reflection-paired E_U components)")
    print("  pairwise observable A: escape margin 11/13-min_C Q")
    print("  switch/gauge B: component width")
    print(
        "  tie Hamiltonian paths: "
        f"{sharp['component_margin_order']} / {sharp['component_width_order']}"
    )
    print("  score histograms: (0,1,2,3) / (0,1,2,3)")
    print("  directed 3-cycles: 0 / 0; SCC sizes: (1,1,1,1); HP counts: 1 / 1")
    print(f"  edge flips between gauges: {sharp['component_flips']}")
    print("  preserved by decorated carrier: exact intervals, endpoint owners, and margins")
    print("  destroyed by bare tournament: metric values and the sign of escape")
    print("PASS")


if __name__ == "__main__":
    main()

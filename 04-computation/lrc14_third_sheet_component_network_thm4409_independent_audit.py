#!/usr/bin/env python3
"""Clean-room exact audit of the third-sheet component-network certificate.

The producer implementation is neither imported nor consulted.  This file
uses only Python's standard library and exact Fraction arithmetic.  Every
proof gate is an explicit RuntimeError check, live under ``python -O``.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd
from hashlib import sha256
import sys


Q = Fraction
LAMBDA = Q(1, 14)
RAW_RADIUS = Q(3, 14)
SHEETS = tuple(permutations((0, 1, 2)))
PAIRS = ((0, 1), (0, 2), (1, 2))
TARGET = Q(6, 77)
CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("CHECK FAILED: " + label)


def floor_q(x):
    return x.numerator // x.denominator


def ceil_q(x):
    return -floor_q(-x)


def gcd_all(values):
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def integer_open_interval(lo, hi):
    return range(floor_q(lo) + 1, ceil_q(hi)) if lo < hi else range(0)


@lru_cache(maxsize=None)
def danger_components(speed, sheet):
    """Open danger arcs split into rational [0,1] chart components.

    The integer label is the nearest integer N for the raw variable y=3x:
    if m is nearest to speed*(x+sheet/3), then N=3m-speed*sheet.
    """
    half_width = LAMBDA / speed
    lower_m = Q(speed * sheet, 3) - LAMBDA
    upper_m = speed + Q(speed * sheet, 3) + LAMBDA
    answer = []
    for m in integer_open_interval(lower_m, upper_m):
        center = Q(m, speed) - Q(sheet, 3)
        lo = max(Q(0), center - half_width)
        hi = min(Q(1), center + half_width)
        if lo < hi:
            answer.append((lo, hi, 3 * m - speed * sheet))
    answer.sort()
    require(sum(hi - lo for lo, hi, _ in answer) == Q(1, 7), "one-sheet measure")
    require(all(answer[i][1] <= answer[i + 1][0] for i in range(len(answer) - 1)),
            "one-sheet components disjoint")
    return tuple(answer)


def pair_components(left, right):
    """Intersect two disjoint labeled interval lists by a two-pointer sweep."""
    answer = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            answer.append((lo, hi, (left[i][2], right[j][2])))
        end = min(left[i][1], right[j][1])
        if left[i][1] == end:
            i += 1
        if right[j][1] == end:
            j += 1
    require(all(answer[i][1] <= answer[i + 1][0] for i in range(len(answer) - 1)),
            "pair components disjoint")
    return tuple(answer)


def interval_edges(left, right):
    """Positive-length contacts, including their actual intersection interval."""
    answer = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            answer.append((i, j, lo, hi))
        end = min(left[i][1], right[j][1])
        if left[i][1] == end:
            i += 1
        if right[j][1] == end:
            j += 1
    return tuple(answer)


def sheet_geometry(w, pair, pi):
    i, j = pair
    k = next(index for index in range(3) if index not in pair)
    first = danger_components(w[i], pi[i])
    second = danger_components(w[j], pi[j])
    pair_list = pair_components(first, second)
    third = danger_components(w[k], pi[k])
    edges = interval_edges(pair_list, third)
    return pair_list, third, edges, k


def physical_sheet(w, pi):
    """Definition-level triple intersection, retaining the raw carrier label."""
    first = danger_components(w[0], pi[0])
    second = danger_components(w[1], pi[1])
    pair_list = pair_components(first, second)
    third = danger_components(w[2], pi[2])
    answer = defaultdict(Fraction)
    pieces = 0
    for index, third_index, lo, hi in interval_edges(pair_list, third):
        n0, n1 = pair_list[index][2]
        n2 = third[third_index][2]
        C = cross(w, (n0, n1, n2))
        answer[C] += hi - lo
        pieces += 1
    return dict(answer), pieces


def physical_comb(w):
    answer = defaultdict(Fraction)
    pieces = 0
    for pi in SHEETS:
        sheet, count = physical_sheet(w, pi)
        for C, length in sheet.items():
            answer[C] += length
        pieces += count
    return dict(answer), pieces


def roof(w, C):
    a, b, c = w
    A, B, Cc = C
    return max(Q(0), min(
        2 * RAW_RADIUS / a,
        2 * RAW_RADIUS / b,
        2 * RAW_RADIUS / c,
        RAW_RADIUS / a + RAW_RADIUS / b - Q(abs(Cc), a * b),
        RAW_RADIUS / a + RAW_RADIUS / c - Q(abs(B), a * c),
        RAW_RADIUS / b + RAW_RADIUS / c - Q(abs(A), b * c),
    ))


def raw_carrier_comb(w):
    """Independent THM-4392-style raw-lattice enumeration."""
    a, b, c = w
    bound_A = RAW_RADIUS * (b + c)
    bound_B = RAW_RADIUS * (a + c)
    bound_C = RAW_RADIUS * (a + b)
    answer = {}
    for A in integer_open_interval(-bound_A, bound_A):
        for B in integer_open_interval(-bound_B, bound_B):
            numerator = -(a * A + b * B)
            if numerator % c:
                continue
            Cc = numerator // c
            C = (A, B, Cc)
            if abs(Cc) >= bound_C or any(x % 3 == 0 for x in C):
                continue
            value = roof(w, C)
            if value > 0:
                answer[C] = value
    return answer


class Dinic:
    """Exact rational max flow with a residual min-cut certificate."""

    def __init__(self, size):
        self.graph = [[] for _ in range(size)]

    def add(self, source, target, capacity):
        forward = [target, len(self.graph[target]), capacity, capacity]
        reverse = [source, len(self.graph[source]), Q(0), Q(0)]
        self.graph[source].append(forward)
        self.graph[target].append(reverse)

    def solve(self, source, sink):
        total = Q(0)
        size = len(self.graph)
        while True:
            level = [-1] * size
            level[source] = 0
            queue = deque((source,))
            while queue:
                u = queue.popleft()
                for edge in self.graph[u]:
                    if edge[2] > 0 and level[edge[0]] < 0:
                        level[edge[0]] = level[u] + 1
                        queue.append(edge[0])
            if level[sink] < 0:
                break
            cursor = [0] * size

            def send(u, amount):
                if u == sink:
                    return amount
                while cursor[u] < len(self.graph[u]):
                    edge = self.graph[u][cursor[u]]
                    if edge[2] > 0 and level[edge[0]] == level[u] + 1:
                        pushed = send(edge[0], min(amount, edge[2]))
                        if pushed:
                            edge[2] -= pushed
                            self.graph[edge[0]][edge[1]][2] += pushed
                            return pushed
                    cursor[u] += 1
                return Q(0)

            while True:
                pushed = send(source, Q(1))
                if not pushed:
                    break
                total += pushed

        reachable = {source}
        queue = deque((source,))
        while queue:
            u = queue.popleft()
            for edge in self.graph[u]:
                if edge[2] > 0 and edge[0] not in reachable:
                    reachable.add(edge[0])
                    queue.append(edge[0])
        balance = [Q(0)] * size
        for u in range(size):
            for edge in self.graph[u]:
                if edge[3] <= 0:
                    continue
                used = edge[3] - edge[2]
                require(Q(0) <= used <= edge[3], "flow respects edge capacity")
                balance[u] += used
                balance[edge[0]] -= used
        require(balance[source] == total and balance[sink] == -total,
                "flow value equals terminal imbalance")
        require(all(balance[node] == 0 for node in range(size) if node not in (source, sink)),
                "flow conservation")
        cut = sum(
            edge[3]
            for u in reachable for edge in self.graph[u]
            if edge[0] not in reachable
        )
        require(total == cut, "fractional max-flow/min-cut certificate")
        return total


def capacity_flow(left_caps, right_caps, edges):
    source = 0
    left_base = 1
    right_base = left_base + len(left_caps)
    sink = right_base + len(right_caps)
    flow = Dinic(sink + 1)
    for i, capacity in enumerate(left_caps):
        flow.add(source, left_base + i, capacity)
    for j, capacity in enumerate(right_caps):
        flow.add(right_base + j, sink, capacity)
    for i, j, _, _ in edges:
        flow.add(left_base + i, right_base + j, Q(1))
    return flow.solve(source, sink)


def sheet_capacity_zero(w, pair, pi):
    pair_list, third, edges, _ = sheet_geometry(w, pair, pi)
    left_caps = tuple(hi - lo for lo, hi, _ in pair_list)
    right_caps = tuple(hi - lo for lo, hi, _ in third)
    kappa = capacity_flow(left_caps, right_caps, edges)
    edge_sum = sum(min(left_caps[i], right_caps[j]) for i, j, _, _ in edges)
    return kappa, edge_sum, pair_list, third, edges


def network_zero(w, pair, audit_edge_sum=False):
    total = Q(0)
    for pi in SHEETS:
        kappa, edge_sum, _, _, _ = sheet_capacity_zero(w, pair, pi)
        if audit_edge_sum:
            require(kappa == edge_sum, "edgewise minima form a feasible max flow")
        require(kappa <= Q(1, 7), "third-sheet source capacity")

        # Direct H=K=0 substitution into the signed network formula.
        G = Q(1, 343)
        cap_g = kappa / 49
        cap_one_minus_g = 48 * kappa / 49
        reduced = G + cap_one_minus_g - max(Q(0), G - cap_g)
        require(reduced == kappa, "degree-zero cancellation")
        total += kappa
    return total


PROFILE_DENOMINATOR = 17
PROFILES = (
    tuple(Q((7 * cell + 3) % 17, 16) for cell in range(PROFILE_DENOMINATOR)),
    tuple(Q((11 * cell + 5) % 17, 16) for cell in range(PROFILE_DENOMINATOR)),
    tuple(Q(1, 49) for _ in range(PROFILE_DENOMINATOR)),
)


def integrate_interval(lo, hi, profile):
    total = Q(0)
    for cell, value in enumerate(profile):
        left = max(lo, Q(cell, PROFILE_DENOMINATOR))
        right = min(hi, Q(cell + 1, PROFILE_DENOMINATOR))
        if left < right:
            total += (right - left) * value
    return total


def weighted_sheet_audit(w, pair, pi, profile):
    """Audit signs and both fractional capacity inequalities for arbitrary g."""
    require(all(Q(0) <= value <= Q(1) for value in profile), "0<=g<=1")
    pair_list, third, edges, _ = sheet_geometry(w, pair, pi)
    left_g = tuple(integrate_interval(lo, hi, profile) for lo, hi, _ in pair_list)
    right_g = tuple(integrate_interval(lo, hi, profile) for lo, hi, _ in third)
    left_one = tuple((hi - lo) - left_g[i] for i, (lo, hi, _) in enumerate(pair_list))
    right_one = tuple((hi - lo) - right_g[i] for i, (lo, hi, _) in enumerate(third))

    edge_g = tuple(integrate_interval(lo, hi, profile) for _, _, lo, hi in edges)
    edge_one = tuple((hi - lo) - edge_g[e] for e, (_, _, lo, hi) in enumerate(edges))
    overlap_length = sum(hi - lo for _, _, lo, hi in edges)
    overlap_g = sum(edge_g, Q(0))
    G = sum(right_g, Q(0))
    R = overlap_length - overlap_g
    L = G - overlap_g
    require(overlap_length == G + R - L, "signed third-sheet identity")
    require(R >= 0 and L >= 0, "positive deficit and leakage")

    for left_index in range(len(pair_list)):
        hits = tuple(e for e, row in enumerate(edges) if row[0] == left_index)
        require(sum((edge_g[e] for e in hits), Q(0)) <= left_g[left_index], "g left feasibility")
        require(sum((edge_one[e] for e in hits), Q(0)) <= left_one[left_index], "1-g left feasibility")
    for right_index in range(len(third)):
        hits = tuple(e for e, row in enumerate(edges) if row[1] == right_index)
        require(sum((edge_g[e] for e in hits), Q(0)) <= right_g[right_index], "g right feasibility")
        require(sum((edge_one[e] for e in hits), Q(0)) <= right_one[right_index], "1-g right feasibility")

    cap_g = capacity_flow(left_g, right_g, edges)
    cap_one = capacity_flow(left_one, right_one, edges)
    require(overlap_g <= cap_g, "fractional g capacity inequality")
    require(R <= cap_one, "fractional 1-g capacity inequality")
    leakage_floor = max(Q(0), G - cap_g)
    require(L >= leakage_floor, "negative leakage lower bound")
    upper = G + cap_one - leakage_floor
    require(overlap_length <= upper, "network upper direction")
    return overlap_length, upper, (G, R, L, cap_g, cap_one, leakage_floor)


def equality_geometry_audit():
    w = (1, 5, 11)
    expected = {
        (0, 1): (6, 6, 68, 62, 6, 6, 0, TARGET),
        (0, 2): (12, 6, 32, 26, 6, 6, 0, TARGET),
        (1, 2): (12, 6, 8, 2, 6, 6, 0, TARGET),
    }
    rows = []
    direct, pieces = physical_comb(w)
    raw = raw_carrier_comb(w)
    require(direct == raw, "equality physical/raw carrier dictionary")
    require(sum(direct.values(), Q(0)) == TARGET and len(direct) == 2 and pieces == 6,
            "equality comb anatomy")

    for pair in PAIRS:
        pair_pieces = touched_pair = third_pieces = isolated_third = 0
        edge_count = nested = crossing = 0
        total_capacity = Q(0)
        for pi_index, pi in enumerate(SHEETS):
            kappa, edge_sum, pair_list, third, edges = sheet_capacity_zero(w, pair, pi)
            require(kappa == edge_sum, "equality edge sum")
            left_degree = Counter(i for i, _, _, _ in edges)
            right_degree = Counter(j for _, j, _, _ in edges)
            require(all(degree <= 1 for degree in left_degree.values()), "matching on pair side")
            require(all(degree <= 1 for degree in right_degree.values()), "matching on third side")
            pair_pieces += len(pair_list)
            touched_pair += len(left_degree)
            third_pieces += len(third)
            isolated_third += len(third) - len(right_degree)
            edge_count += len(edges)
            total_capacity += kappa
            for i, j, lo, hi in edges:
                I = pair_list[i]
                J = third[j]
                is_nested = ((I[0] >= J[0] and I[1] <= J[1])
                             or (J[0] >= I[0] and J[1] <= I[1]))
                nested += is_nested
                crossing += not is_nested

            # Matching + nesting makes both capacity relaxations exact for
            # every nonnegative weight; three jagged rational controls test it.
            for profile in PROFILES:
                physical, upper, _ = weighted_sheet_audit(w, pair, pi, profile)
                require(physical == upper, "nested matching gives weighted exactness")

        row = (pair_pieces, touched_pair, third_pieces, isolated_third,
               edge_count, nested, crossing, total_capacity)
        require(row == expected[pair], "all equality graph totals " + repr(pair))
        require(network_zero(w, pair, audit_edge_sum=True) == TARGET,
                "all equality pair choices exact at degree zero")
        rows.append((pair, row))
    return tuple(rows), direct


def bounded_universe(height):
    speeds = tuple(x for x in range(1, height + 1, 2) if x % 3)
    return tuple(w for w in combinations(speeds, 3) if gcd_all(w) == 1)


def digest(rows):
    return sha256((repr(rows) + "\n").encode("ascii")).hexdigest()


EXPECTED_CENSUS_DIGEST = "fdc37d329434215b5442fd53e91af5d5f716f7165cba75530248f6e6b5c589ce"


def main():
    if "--tripwire" in sys.argv:
        require(False, "optimization-live tripwire")

    equality_rows, equality_dictionary = equality_geometry_audit()
    universe = bounded_universe(79)
    require(len(universe) == 2910, "height-79 hostile universe size")

    fixed_success = [0, 0, 0]
    selected = [0, 0, 0]
    physical_strict = physical_equal = physical_above = 0
    network_exact = 0
    network_target = []
    raw_counts = []
    records = []
    strict_gaps = []
    first_relaxed = None

    for w in universe:
        physical_dictionary, physical_pieces = physical_comb(w)
        raw_dictionary = raw_carrier_comb(w)
        require(physical_dictionary == raw_dictionary, "physical/raw hostile dictionary " + repr(w))
        mass = sum(physical_dictionary.values(), Q(0))
        raw_counts.append(len(raw_dictionary))
        physical_strict += mass < TARGET
        physical_equal += mass == TARGET
        physical_above += mass > TARGET

        bounds = []
        for pair_index, pair in enumerate(PAIRS):
            bound = network_zero(w, pair, audit_edge_sum=True)
            require(bound >= mass, "degree-zero network upper bound " + repr((w, pair)))
            fixed_success[pair_index] += bound <= TARGET
            bounds.append(bound)
        best_index, best = min(enumerate(bounds), key=lambda row: (row[1], row[0]))
        selected[best_index] += 1
        require(best <= TARGET, "minimized degree-zero certificate closes hostile row " + repr(w))
        network_exact += best == mass
        if best == TARGET:
            network_target.append(w)
        else:
            strict_gaps.append((TARGET - best, w, PAIRS[best_index], best))
        if first_relaxed is None and best > mass:
            first_relaxed = (w, PAIRS[best_index], mass, best)
        records.append((w, mass, tuple(bounds), best_index, best,
                        len(raw_dictionary), physical_pieces))

    require((physical_strict, physical_equal, physical_above) == (2909, 1, 0),
            "physical target trichotomy")
    require(tuple(fixed_success) == (2818, 2855, 2859), "fixed-pair success counts")
    require(tuple(selected) == (400, 533, 1977), "tie-broken minimizing pair counts")
    require(network_exact == 1747, "network exact-row count")
    require(tuple(network_target) == ((1, 5, 11),), "sole network target row")
    require((min(raw_counts), max(raw_counts)) == (0, 22), "raw carrier-count range")

    weakest = min(strict_gaps)
    require(weakest == (Q(6, 1771), (1, 11, 23), (0, 1), Q(12, 161)),
            "weakest strict certificate")
    census_digest = digest(tuple(records))
    require(census_digest == EXPECTED_CENSUS_DIGEST, "height-79 census digest")

    # A genuinely relaxed row is the hostile for signs/capacities.  Recheck
    # every sheet, all pairs, and jagged weights independently of degree zero.
    require(first_relaxed is not None, "strict graph-relaxation hostile exists")
    require(first_relaxed == ((1, 19, 79), (0, 2), Q(108, 10507), Q(8, 553)),
            "first strict graph-relaxation hostile")
    weighted_rows = []
    hostile_w = first_relaxed[0]
    for pair in PAIRS:
        for pi_index, pi in enumerate(SHEETS):
            for profile_index, profile in enumerate(PROFILES):
                physical, upper, scalars = weighted_sheet_audit(hostile_w, pair, pi, profile)
                weighted_rows.append((pair, pi_index, profile_index, physical, upper, scalars))

    print("status=PASS")
    print("scope=third_sheet_component_network_local_three_speed_certificate; LRC(14)=OPEN")
    print("measure_identity=mu_sheet=G+R-L; R=intersection_(1-g); L=outside_g")
    print("fractional_capacity=actual_edge_integrals_feasible; R<=Cap(1-g); L>=max(0,G-Cap(g))")
    print("degree_zero=G=1/343; Cap(g)=kappa/49; Cap(1-g)=48*kappa/49; U=kappa")
    print("equality_rows=" + repr(equality_rows))
    print("equality_dictionary=" + repr(equality_dictionary))
    print("height79_summary=" + repr((len(universe), physical_strict, physical_equal,
                                             physical_above, tuple(fixed_success), tuple(selected),
                                             network_exact, tuple(network_target),
                                             (min(raw_counts), max(raw_counts)), weakest)))
    print("first_relaxed_hostile=" + repr(first_relaxed))
    print("weighted_hostile_rows=" + repr(tuple(weighted_rows)))
    print("census_digest=" + census_digest)
    print("danger_cache=" + repr(danger_components.cache_info()))
    print("optimization_safe_checks=yes")
    print("explicit_checks=" + str(CHECKS))


if __name__ == "__main__":
    main()

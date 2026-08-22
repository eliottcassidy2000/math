#!/usr/bin/env python3
"""Independent hostile audit of THM-3366's refined k=3 composition.

The disputed composition has two inputs, and this audit rebuilds both without
importing either S174's recursive complement solver or S176's composition
engine.

* For the pool-14 terminal set it imports only the separately audited fixed
  arrangement/subset-enumeration geometry, then reconstructs every k=3
  body/divisor row from the pinned THM-2928 safe-cell word.
* For the final one-spike ledger it uses a literal denominator-symbol/current-
  lcm dynamic program, not divisor-Mobius grouped feature extraction.  The
  small-mask status bound is recomputed from vertices of the primal coupling
  polytope, rather than from the canonical dual fractional-cover routine.
* It intersects literal ``(body,D)`` keys and reconstructs both occurrence and
  feature-union summaries after deletion.

Hostile controls deliberately try body-only, divisor-only, and ``(L,D)`` row
keys; raw-support subtraction and double subtraction; and an overbroad reading
of finite pool-14 failure.  Every guard is a RuntimeError and remains active
under ``python -O``.
"""

from bisect import bisect_left
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = ROOT / "04-computation/lrc14_two_drift_body_projection_support_thm2928.py"
FIXED_GEOMETRY_PATH = (
    ROOT / "04-computation/lrc14_allk_complement_clock_independent_audit_20260814.py"
)
FIXED_GEOMETRY_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_allk_complement_clock_independent_audit_20260814.out"
)
TARGET_COMPOSITION_PATH = (
    ROOT / "04-computation/lrc14_k3_refined_complement_clock_composition_kps_s176.py"
)
TARGET_COMPOSITION_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_k3_refined_complement_clock_composition_kps_s176.out"
)
TARGET_ONE_SPIKE_PATH = (
    ROOT / "04-computation/lrc14_k3_four_drift_expected_spike_gf_thm2928.py"
)

EXPECTED_SUPPORT_LF_SHA256 = (
    "778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a"
)
EXPECTED_FIXED_GEOMETRY_LF_SHA256 = (
    "a43610abcb1d591f4b606677e3e6a9a32897b38fa2fd591a281ee8c4dc07d823"
)
EXPECTED_FIXED_GEOMETRY_OUTPUT_LF_SHA256 = (
    "31c096e528a148cd717f4d9c7485c93cdec6175b33b83bae9d2be627249d743c"
)
EXPECTED_TARGET_COMPOSITION_LF_SHA256 = (
    "27e4bff52705189bf8ff73db42d76d4e2fc94c44330d295f166f8d4217cb1804"
)
ADVERTISED_TARGET_COMPOSITION_SHA256 = (
    "127ef53b27f10a5c61ac273a49b13a5ae56ea4fa98809df8f6fd9accdce89d97"
)
EXPECTED_TARGET_COMPOSITION_OUTPUT_LF_SHA256 = (
    "4cb8f95113123007af9fb5a1f58b3b5373dd4637615b50be7168c6bf578b696b"
)
EXPECTED_TARGET_ONE_SPIKE_LF_SHA256 = (
    "05e365a654b32e66b814dcbce9385a2d13c22a2c84a5474e0855dcab6262b055"
)

POOL = tuple(range(1, 15))
K = 3
ARITY = 4
COVER_BUDGET = 5
SUPPORT_CUTOFF = Q(125, 143)
SAFE_MASS = Q(55, 91)

EXPECTED_BODY_COUNT = 3003
EXPECTED_BODY_DIVISOR_ROWS = 251_536
EXPECTED_SUPPORT_ROWS = 26_970
EXPECTED_SUPPORT_DIVISORS = 217
EXPECTED_TERMINAL_ROWS = 19_053
EXPECTED_TERMINAL_HISTOGRAM = ((1, 12_659), (2, 2_764), (3, 976), (4, 1_052), (5, 1_602))
EXPECTED_RAW_TERMINAL_OCCURRENCES = 5_887_257_171
EXPECTED_PRE = (398_241_574, 2_548_901_482, 1_904, 1_823, 107)
EXPECTED_PRE_BY_C = ((1, 71_619_386), (2, 1_351_841_956), (3, 1_065_317_472), (4, 60_122_668))
EXPECTED_PRE_SHAPES_BY_C = ((1, 30_050_140), (2, 255_882_564), (3, 103_263_188), (4, 9_045_682))
EXPECTED_PRE_SEMANTIC = "06b3e5a7a05d5c9a2f1633d74061f1605c2aab856b315026ad697245bac84964"
EXPECTED_DELETION = (7, 7_648)
EXPECTED_POST = (398_241_574, 2_548_893_834, 1_897, 1_823, 107)
EXPECTED_POST_BY_C = ((1, 71_619_386), (2, 1_351_841_956), (3, 1_065_309_824), (4, 60_122_668))
EXPECTED_POST_SHAPES_BY_C = EXPECTED_PRE_SHAPES_BY_C
EXPECTED_POST_SEMANTIC = "1dd557b3629ed4c48f2b85412b48476b29e300990f38bdfe65cef2372fbab8d4"
EXPECTED_HIT_BODIES = tuple(
    sorted(tuple(sorted((2, 6, 8, 10, 14, odd))) for odd in range(1, 14, 2))
)
EXPECTED_HIT_DIVISORS = (5_880, 5_880, 5_880, 5_880, 17_640, 64_680, 76_440)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def escaped_literal_sha256(path):
    """Hostile control for the historical hash-basis transcription error."""
    return sha256(path.read_bytes().replace(b"\\r\\n", b"\\n")).hexdigest()


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("module import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=None)
def divisors(number):
    low = []
    high = []
    candidate = 1
    while candidate * candidate <= number:
        if number % candidate == 0:
            low.append(candidate)
            quotient = number // candidate
            if quotient != candidate:
                high.append(quotient)
        candidate += 1
    return tuple(low + high[::-1])


def mobius(number):
    result = 1
    remaining = number
    prime = 2
    while prime * prime <= remaining:
        if remaining % prime:
            prime += 1
            continue
        remaining //= prime
        if remaining % prime == 0:
            return 0
        result = -result
        while remaining % prime == 0:
            remaining //= prime
        prime += 1
    return -result if remaining > 1 else result


@lru_cache(maxsize=None)
def exact_lcm_multiset_count(modulus, arity):
    total = 0
    for downward in divisors(modulus):
        alphabet = len(divisors(downward)) - 1
        total += mobius(modulus // downward) * comb(alphabet + arity - 1, arity)
    require(total >= 0, ("negative exact-lcm count", modulus, arity, total))
    return total


def projected_arcs(modulus, safe_ranges):
    """Independent merged cyclic projection of the exact safe-cell ranges."""
    pieces = []
    for left, right in safe_ranges:
        length = right - left
        if length >= modulus:
            return ((0, modulus),)
        start = left % modulus
        end = start + length
        if end <= modulus:
            pieces.append((start, end))
        else:
            pieces.append((start, modulus))
            pieces.append((0, end - modulus))
    pieces.sort()
    merged = []
    for left, right in pieces:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return tuple(merged)


def residue_load_histogram(arcs, modulus):
    """Direct residue loop, intentionally not the canonical event sweep."""
    loads = [0] * modulus
    for left, right in arcs:
        for residue in range(modulus):
            first = left + ((residue - left) % modulus)
            if first < right:
                loads[residue] += 1 + (right - 1 - first) // modulus
    histogram = tuple(sorted(Counter(loads).items()))
    require(sum(count for _load, count in histogram) == modulus, "histogram width")
    require(
        sum(load * count for load, count in histogram)
        == sum(right - left for left, right in arcs),
        "histogram mass",
    )
    return histogram


def largest_class_load(histogram):
    return histogram[-1][0]


def solve_square(rows, rhs):
    size = len(rows)
    matrix = [
        [Q(value) for value in row] + [Q(rhs[index])]
        for index, row in enumerate(rows)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if matrix[row][column]),
            None,
        )
        if pivot is None:
            return None
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        scale = matrix[column][column]
        matrix[column] = [value / scale for value in matrix[column]]
        for row in range(size):
            if row == column or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(matrix[row], matrix[column])
            ]
    return tuple(matrix[index][-1] for index in range(size))


@lru_cache(maxsize=None)
def coupling_vertices(marginals):
    """Vertices of the primal Boolean coupling polytope, exactly over Q."""
    marginals = tuple(marginals)
    count = len(marginals)
    if count == 0:
        return ((Q(1),),)
    states = tuple(range(1 << count))
    rhs = (Q(1),) + marginals
    vertices = set()
    for basis in combinations(states, count + 1):
        rows = [tuple(Q(1) for _state in basis)]
        rows.extend(
            tuple(Q((state >> coordinate) & 1) for state in basis)
            for coordinate in range(count)
        )
        solution = solve_square(tuple(rows), rhs)
        if solution is None or any(value < 0 for value in solution):
            continue
        vector = [Q(0)] * len(states)
        for state, value in zip(basis, solution):
            vector[state] = value
        require(sum(vector) == 1, ("coupling normalization", marginals, basis))
        for coordinate, marginal in enumerate(marginals):
            require(
                sum(vector[state] for state in states if (state >> coordinate) & 1)
                == marginal,
                ("coupling marginal", marginals, basis, coordinate),
            )
        vertices.add(tuple(vector))
    require(vertices, ("coupling polytope has no vertices", marginals))
    return tuple(sorted(vertices))


@lru_cache(maxsize=None)
def maximum_upward_mass(weights, marginals, need):
    weights = tuple(weights)
    marginals = tuple(marginals)
    require(len(weights) == len(marginals), "status arity")
    if need <= 0:
        return Q(1)
    if need > sum(weights):
        return Q(0)
    loads = tuple(
        sum(weight for index, weight in enumerate(weights) if (state >> index) & 1)
        for state in range(1 << len(weights))
    )
    return max(
        sum(probability for probability, load in zip(vertex, loads) if load >= need)
        for vertex in coupling_vertices(marginals)
    )


@lru_cache(maxsize=None)
def status_need_limit(weights, marginals):
    """Primal-vertex reconstruction of the strict 55/91 status threshold."""
    weights = tuple(weights)
    marginals = tuple(marginals)
    if not weights:
        return 0
    thresholds = sorted(
        {
            sum(weights[index] for index in range(len(weights)) if (mask >> index) & 1)
            for mask in range(1, 1 << len(weights))
        }
    )
    result = 0
    for threshold in thresholds:
        if maximum_upward_mass(weights, marginals, threshold) > SAFE_MASS:
            result = threshold
        else:
            break
    require(maximum_upward_mass(weights, marginals, result) > SAFE_MASS, "status lower")
    if result < sum(weights):
        require(
            maximum_upward_mass(weights, marginals, result + 1) <= SAFE_MASS,
            "status maximality",
        )
    return result


def small_status_limit(pattern, modulus, small_loads):
    weights = []
    marginals = []
    for denominator, copies in zip((2, 3, 4), pattern):
        if not copies:
            continue
        require(modulus % denominator == 0, ("small denominator", modulus, pattern))
        weight = small_loads[denominator]
        require(0 < weight <= modulus // denominator, ("small load", modulus, pattern))
        weights.extend([weight] * copies)
        marginals.extend([Q(denominator, 7)] * copies)
    return status_need_limit(tuple(weights), tuple(marginals))


def symbol_feature(modulus, denominator):
    q = modulus // 7
    pattern = (
        int(denominator == 2),
        int(denominator == 3),
        int(denominator == 4),
    )
    quotient_type = denominator // gcd(denominator, q)
    require(quotient_type in (1, 7), ("q-fibre dichotomy", modulus, denominator))
    if q % denominator:
        require(denominator % 7 == 0 and not any(pattern), ("uniform symbol", modulus, denominator))
        return pattern, 1, q, False
    capacity = 0
    if not any(pattern):
        capacity = (modulus // denominator) * ((denominator + 6) // 7)
    return pattern, 0, capacity, True


def one_spike_allowance(q, denominator):
    require(q % denominator == 0, ("spike denominator", q, denominator))
    quotient, remainder = divmod(denominator, 7)
    width = q // denominator
    return quotient * width + (width if remainder >= 5 else 0)


@lru_cache(maxsize=None)
def literal_feature_distributions(modulus):
    """Current-lcm multiset DP; no grouped Mobius feature extraction."""
    require(modulus % 7 == 0, ("nonseptimal modulus", modulus))
    q = modulus // 7
    # state=(used,current_lcm,m2,m3,m4,c,capacity,spike_marker)
    # spike_marker is 0 for none, d for exactly one, and -1 for at least two.
    states = {(0, 1, 0, 0, 0, 0, 0, 0): 1}
    for denominator in divisors(modulus)[1:]:
        pattern, unit_c, unit_capacity, is_spike = symbol_feature(modulus, denominator)
        additions = Counter()
        for state, multiplicity in tuple(states.items()):
            used, current, m2, m3, m4, c, capacity, spike_marker = state
            for copies in range(1, ARITY - used + 1):
                new_marker = spike_marker
                if is_spike:
                    if copies == 1 and spike_marker == 0:
                        new_marker = denominator
                    else:
                        new_marker = -1
                additions[
                    (
                        used + copies,
                        lcm(current, denominator),
                        m2 + copies * pattern[0],
                        m3 + copies * pattern[1],
                        m4 + copies * pattern[2],
                        c + copies * unit_c,
                        capacity + copies * unit_capacity,
                        new_marker,
                    )
                ] += multiplicity
        for new_state, multiplicity in additions.items():
            states[new_state] = states.get(new_state, 0) + multiplicity

    distribution = Counter()
    c3_distribution = Counter()
    for state, multiplicity in states.items():
        used, current, m2, m3, m4, c, capacity, spike_marker = state
        if used != ARITY or current != modulus:
            continue
        pattern = (m2, m3, m4)
        distribution[(pattern, c, capacity)] += multiplicity
        if c == 3:
            require(spike_marker > 0, ("c3 spike marker", modulus, state))
            feature = (
                pattern,
                spike_marker,
                capacity,
                one_spike_allowance(q, spike_marker),
            )
            c3_distribution[feature] += multiplicity

    require(all(value > 0 for value in distribution.values()), ("distribution sign", modulus))
    require(
        sum(distribution.values()) == exact_lcm_multiset_count(modulus, ARITY),
        ("current-lcm/Mobius total", modulus),
    )
    collapsed_c3 = Counter()
    for (pattern, _spike, capacity, _allowance), multiplicity in c3_distribution.items():
        collapsed_c3[(pattern, 3, capacity)] += multiplicity
    require(
        collapsed_c3
        == Counter({feature: count for feature, count in distribution.items() if feature[1] == 3}),
        ("c3 collapse", modulus),
    )
    return distribution, c3_distribution


def literal_distribution_controls():
    cases = 0
    shapes = 0
    for modulus in (14, 28, 42, 56, 70, 84):
        q = modulus // 7
        alphabet = divisors(modulus)[1:]
        brute = Counter()
        brute_c3 = Counter()
        for values in combinations_with_replacement(alphabet, ARITY):
            if lcm(*values) != modulus:
                continue
            pattern = tuple(values.count(denominator) for denominator in (2, 3, 4))
            c = 0
            capacity = 0
            spikes = []
            for denominator in values:
                _pattern, unit_c, unit_capacity, is_spike = symbol_feature(modulus, denominator)
                c += unit_c
                capacity += unit_capacity
                if is_spike:
                    spikes.append(denominator)
            brute[(pattern, c, capacity)] += 1
            if c == 3:
                require(len(spikes) == 1, ("literal c3 spike", modulus, values))
                denominator = spikes[0]
                brute_c3[(pattern, denominator, capacity, one_spike_allowance(q, denominator))] += 1
            shapes += 1
        distribution, c3_distribution = literal_feature_distributions(modulus)
        require(brute == distribution, ("literal distribution", modulus))
        require(brute_c3 == c3_distribution, ("literal c3 distribution", modulus))
        cases += 1
    return cases, shapes


def mean_passes(demand, c, q):
    return demand == 0 or 55 * demand < 13 * (4 - c) * q


def suffix(counter):
    keys = sorted(counter)
    totals = [0] * (len(keys) + 1)
    for index in range(len(keys) - 1, -1, -1):
        totals[index] = totals[index + 1] + counter[keys[index]]
    return keys, totals


def at_least(data, threshold):
    keys, totals = data
    return totals[bisect_left(keys, threshold)]


def build_support_and_terminal_rows(support, fixed):
    points = fixed.arrangement_points(POOL)
    atoms = tuple(zip(points, points[1:]))
    require((len(points), len(atoms)) == (206, 205), "pool-14 arrangement")
    clock_masks = fixed.make_clock_masks(points, atoms)
    subsets = fixed.enumerate_subset_masks(clock_masks)
    cover_cache = {}
    target_cover = {}
    by_divisor = defaultdict(list)
    terminal_keys = set()
    terminal_completions = {}
    terminal_histogram = Counter()
    terminal_semantic = sha256()
    body_count = 0
    body_divisor_rows = 0

    for body in combinations(range(1, 15), 6):
        body_count += 1
        L, safe_ranges = support.safe_cell_ranges(body)
        for modulus in divisors(L):
            body_divisor_rows += 1
            arcs = projected_arcs(modulus, safe_ranges)
            support_count = sum(right - left for left, right in arcs)
            require(
                support_count == support.support_size_bitset(modulus, safe_ranges),
                ("support projection", body, modulus),
            )
            if Q(support_count, modulus) > SUPPORT_CUTOFF:
                continue
            gaps = fixed.unsupported_gaps(modulus, arcs)
            target, _atom_target, _wrong_target, _endpoint_mask = fixed.target_masks(
                modulus, gaps, points, atoms
            )
            if target not in target_cover:
                target_cover[target] = fixed.least_cover(
                    target, subsets, COVER_BUDGET, cover_cache
                )
            completion = target_cover[target]
            record = (support_count, body, L, arcs, completion)
            by_divisor[modulus].append(record)
            if completion is not None:
                key = (body, modulus)
                terminal_keys.add(key)
                terminal_completions[key] = completion
                terminal_histogram[len(completion)] += 1

    for modulus in sorted(by_divisor):
        by_divisor[modulus].sort(key=lambda row: (row[0], row[1], row[2]))
        for support_count, body, L, _arcs, completion in by_divisor[modulus]:
            if completion is not None:
                terminal_semantic.update(
                    f"{body}|{L}|{modulus}|{support_count}|{completion}\n".encode()
                )

    require(body_count == EXPECTED_BODY_COUNT, ("body count", body_count))
    require(body_divisor_rows == EXPECTED_BODY_DIVISOR_ROWS, ("body/divisor rows", body_divisor_rows))
    require(sum(map(len, by_divisor.values())) == EXPECTED_SUPPORT_ROWS, "k3 support rows")
    require(len(by_divisor) == EXPECTED_SUPPORT_DIVISORS, "k3 support divisors")
    require(len(terminal_keys) == EXPECTED_TERMINAL_ROWS, "terminal rows")
    require(tuple(sorted(terminal_histogram.items())) == EXPECTED_TERMINAL_HISTOGRAM, "terminal histogram")
    return {
        "by_divisor": by_divisor,
        "terminal_keys": terminal_keys,
        "terminal_completions": terminal_completions,
        "terminal_histogram": terminal_histogram,
        "terminal_semantic": terminal_semantic.hexdigest(),
        "distinct_targets": len(target_cover),
        "cover_cache": len(cover_cache),
    }


def reconstruct_one_spike_ledger(by_divisor):
    rows = {}
    feature_multiplicity = {}
    feature_c = {}
    distribution_semantic = sha256()
    primal_patterns = set()

    for modulus in sorted(by_divisor):
        q = modulus // 7
        distribution, c3_distribution = literal_feature_distributions(modulus)
        for feature, multiplicity in sorted(distribution.items()):
            distribution_semantic.update(f"{modulus}|all|{feature}|{multiplicity}\n".encode())
        for feature, multiplicity in sorted(c3_distribution.items()):
            distribution_semantic.update(f"{modulus}|c3|{feature}|{multiplicity}\n".encode())

        by_pattern_c = defaultdict(Counter)
        for (pattern, c, capacity), multiplicity in distribution.items():
            by_pattern_c[(pattern, c)][capacity] += multiplicity
        suffixes = {key: suffix(counter) for key, counter in by_pattern_c.items()}

        for support_count, body, L, arcs, _completion in by_divisor[modulus]:
            q_histogram = residue_load_histogram(arcs, q)
            demand = {
                c: sum(count for load, count in q_histogram if load > c)
                for c in range(5)
            }
            small_loads = {}
            for denominator in (2, 3, 4):
                if modulus % denominator == 0:
                    small_loads[denominator] = largest_class_load(
                        residue_load_histogram(arcs, denominator)
                    )

            limits = {}
            for pattern, _c in by_pattern_c:
                if pattern not in limits:
                    limits[pattern] = small_status_limit(pattern, modulus, small_loads)
                    primal_patterns.add((pattern, tuple(small_loads.get(d) for d in (2, 3, 4))))

            row_count = 0
            row_by_c = Counter()
            row_features = set()
            for (pattern, c), data in suffixes.items():
                if c == 3 or not mean_passes(demand[c], c, q):
                    continue
                threshold = support_count - limits[pattern]
                count = at_least(data, threshold)
                if not count:
                    continue
                row_count += count
                row_by_c[c] += count
                capacities = data[0][bisect_left(data[0], threshold):]
                for capacity in capacities:
                    feature_id = ("mean", modulus, pattern, c, capacity)
                    row_features.add(feature_id)
                    feature_multiplicity[feature_id] = distribution[(pattern, c, capacity)]
                    feature_c[feature_id] = c

            if mean_passes(demand[3], 3, q):
                for feature, multiplicity in c3_distribution.items():
                    pattern, spike_denominator, capacity, allowance = feature
                    if capacity + limits[pattern] < support_count:
                        continue
                    if allowance < demand[3]:
                        continue
                    row_count += multiplicity
                    row_by_c[3] += multiplicity
                    feature_id = (
                        "one-spike",
                        modulus,
                        pattern,
                        spike_denominator,
                        capacity,
                        allowance,
                    )
                    row_features.add(feature_id)
                    feature_multiplicity[feature_id] = multiplicity
                    feature_c[feature_id] = 3

            if row_count:
                key = (body, modulus)
                require(key not in rows, ("duplicate row key", key))
                rows[key] = {
                    "L": L,
                    "support": support_count,
                    "count": row_count,
                    "by_c": row_by_c,
                    "features": frozenset(row_features),
                }

    return {
        "rows": rows,
        "feature_multiplicity": feature_multiplicity,
        "feature_c": feature_c,
        "distribution_semantic": distribution_semantic.hexdigest(),
        "primal_status_profiles": len(primal_patterns),
    }


def summarize(ledger, allowed_keys=None):
    rows = ledger["rows"]
    selected = set(rows) if allowed_keys is None else set(allowed_keys)
    selected &= set(rows)
    occurrences = 0
    bodies = set()
    moduli = set()
    features = set()
    by_c = Counter()
    semantic = sha256()
    for body, modulus in sorted(selected, key=lambda key: (key[1], rows[key]["support"], key[0], rows[key]["L"])):
        record = rows[(body, modulus)]
        occurrences += record["count"]
        bodies.add(body)
        moduli.add(modulus)
        features.update(record["features"])
        by_c.update(record["by_c"])
        semantic.update(
            f"{body}|{record['L']}|{modulus}|{record['support']}|{record['count']}\n".encode()
        )
    shape_by_c = Counter()
    for feature in features:
        multiplicity = ledger["feature_multiplicity"][feature]
        shape_by_c[ledger["feature_c"][feature]] += multiplicity
    shapes = sum(shape_by_c.values())
    return {
        "summary": (shapes, occurrences, len(selected), len(bodies), len(moduli)),
        "keys": selected,
        "features": features,
        "by_c": by_c,
        "shape_by_c": shape_by_c,
        "semantic": semantic.hexdigest(),
    }


def occurrence_sum(rows, keys):
    return sum(rows[key]["count"] for key in keys)


def first_false_positive(candidate_keys, correct_keys):
    false = candidate_keys - correct_keys
    require(false, "hostile key test found no false positive")
    return min(false, key=lambda key: (key[1], key[0]))


def main():
    dependencies = (
        (SUPPORT_PATH, EXPECTED_SUPPORT_LF_SHA256),
        (FIXED_GEOMETRY_PATH, EXPECTED_FIXED_GEOMETRY_LF_SHA256),
        (FIXED_GEOMETRY_OUTPUT_PATH, EXPECTED_FIXED_GEOMETRY_OUTPUT_LF_SHA256),
        (TARGET_COMPOSITION_PATH, EXPECTED_TARGET_COMPOSITION_LF_SHA256),
        (TARGET_COMPOSITION_OUTPUT_PATH, EXPECTED_TARGET_COMPOSITION_OUTPUT_LF_SHA256),
        (TARGET_ONE_SPIKE_PATH, EXPECTED_TARGET_ONE_SPIKE_LF_SHA256),
    )
    for path, expected in dependencies:
        require(lf_sha256(path) == expected, ("dependency changed", path, lf_sha256(path)))
    require(
        escaped_literal_sha256(TARGET_COMPOSITION_PATH)
        == ADVERTISED_TARGET_COMPOSITION_SHA256,
        "advertised S176 hash no longer has the diagnosed escaped-literal mechanism",
    )

    support = load_module("k3_refined_independent_support", SUPPORT_PATH)
    fixed = load_module("k3_refined_independent_fixed_geometry", FIXED_GEOMETRY_PATH)
    literal_control_cases, literal_control_shapes = literal_distribution_controls()
    geometry = build_support_and_terminal_rows(support, fixed)
    ledger = reconstruct_one_spike_ledger(geometry["by_divisor"])

    before = summarize(ledger)
    require(before["summary"] == EXPECTED_PRE, ("pre summary", before["summary"]))
    require(tuple(sorted(before["by_c"].items())) == EXPECTED_PRE_BY_C, ("pre by-c", before["by_c"]))
    require(tuple(sorted(before["shape_by_c"].items())) == EXPECTED_PRE_SHAPES_BY_C, ("pre shapes by-c", before["shape_by_c"]))
    require(before["semantic"] == EXPECTED_PRE_SEMANTIC, ("pre semantic", before["semantic"]))

    correct_hits = before["keys"] & geometry["terminal_keys"]
    killed_occurrences = occurrence_sum(ledger["rows"], correct_hits)
    require((len(correct_hits), killed_occurrences) == EXPECTED_DELETION, ("deletion", len(correct_hits), killed_occurrences))
    after = summarize(ledger, before["keys"] - correct_hits)
    require(after["summary"] == EXPECTED_POST, ("post summary", after["summary"]))
    require(tuple(sorted(after["by_c"].items())) == EXPECTED_POST_BY_C, ("post by-c", after["by_c"]))
    require(tuple(sorted(after["shape_by_c"].items())) == EXPECTED_POST_SHAPES_BY_C, ("post shapes by-c", after["shape_by_c"]))
    require(after["semantic"] == EXPECTED_POST_SEMANTIC, ("post semantic", after["semantic"]))
    require(before["summary"][1] - after["summary"][1] == killed_occurrences, "occurrence difference")

    hit_details = tuple(
        sorted(
            (
                modulus,
                body,
                ledger["rows"][(body, modulus)]["L"],
                geometry["terminal_completions"][(body, modulus)],
                ledger["rows"][(body, modulus)]["count"],
                tuple(sorted(ledger["rows"][(body, modulus)]["by_c"].items())),
            )
            for body, modulus in correct_hits
        )
    )
    require(tuple(sorted(body for _D, body, _L, _cover, _count, _by_c in hit_details)) == EXPECTED_HIT_BODIES, "hit family")
    require(tuple(D for D, _body, _L, _cover, _count, _by_c in hit_details) == EXPECTED_HIT_DIVISORS, "hit divisors")
    require(all(D * 2 == L for D, _body, L, _cover, _count, _by_c in hit_details), "D=L/2 family")
    require(all(by_c == ((3, count),) for _D, _body, _L, _cover, count, by_c in hit_details), "hit c scope")

    raw_terminal_occurrences = sum(
        exact_lcm_multiset_count(modulus, ARITY)
        for _body, modulus in geometry["terminal_keys"]
    )
    require(raw_terminal_occurrences == EXPECTED_RAW_TERMINAL_OCCURRENCES, "raw terminal occurrences")
    raw_hit_occurrences = sum(exact_lcm_multiset_count(modulus, ARITY) for _body, modulus in correct_hits)
    naive_global_subtraction = before["summary"][1] - raw_terminal_occurrences
    naive_hit_subtraction = before["summary"][1] - raw_hit_occurrences
    require(naive_global_subtraction != after["summary"][1], "raw global subtraction accidentally matched")
    require(naive_hit_subtraction != after["summary"][1], "raw hit subtraction accidentally matched")
    require(before["features"] == after["features"], "shape union should survive row deletion")
    exclusive_features = before["features"] - after["features"]
    require(not exclusive_features, "killed rows owned an exclusive feature")
    double_subtraction = before["summary"][1] - 2 * killed_occurrences
    require(double_subtraction != after["summary"][1], "double subtraction accidentally matched")

    terminal_bodies = {body for body, _modulus in geometry["terminal_keys"]}
    terminal_moduli = {modulus for _body, modulus in geometry["terminal_keys"]}
    terminal_LD = {
        (record["L"], modulus)
        for (body, modulus), record in ledger["rows"].items()
        if (body, modulus) in geometry["terminal_keys"]
    }
    wrong_body_hits = {key for key in before["keys"] if key[0] in terminal_bodies}
    wrong_modulus_hits = {key for key in before["keys"] if key[1] in terminal_moduli}
    wrong_LD_hits = {
        key
        for key in before["keys"]
        if (ledger["rows"][key]["L"], key[1]) in terminal_LD
    }
    wrong_key_controls = {
        "body_only": (
            len(wrong_body_hits),
            occurrence_sum(ledger["rows"], wrong_body_hits),
            first_false_positive(wrong_body_hits, correct_hits),
        ),
        "divisor_only": (
            len(wrong_modulus_hits),
            occurrence_sum(ledger["rows"], wrong_modulus_hits),
            first_false_positive(wrong_modulus_hits, correct_hits),
        ),
        "L_and_D": (
            len(wrong_LD_hits),
            occurrence_sum(ledger["rows"], wrong_LD_hits),
            first_false_positive(wrong_LD_hits, correct_hits),
        ),
    }
    require(all(value[:2] != EXPECTED_DELETION for value in wrong_key_controls.values()), "wrong key matched deletion")

    scope_survivors = before["keys"] - geometry["terminal_keys"]
    require(len(scope_survivors) == EXPECTED_POST[2], "scope survivor count")
    scope_key = min(scope_survivors, key=lambda key: (key[1], key[0]))
    scope_record = ledger["rows"][scope_key]
    require(scope_key not in geometry["terminal_completions"], "scope survivor was terminal")
    scope_witness = (
        scope_key[1],
        scope_key[0],
        scope_record["L"],
        scope_record["support"],
        scope_record["count"],
    )

    hit_semantic = sha256()
    for detail in hit_details:
        hit_semantic.update(f"{detail}\n".encode())

    print("LRC14 K=3 REFINED COMPLEMENT-CLOCK INDEPENDENT HOSTILE AUDIT")
    print("status=FINITE-EXACT independent current-lcm/primal-coupling/keywise replay of THM-3366 section 7")
    print(f"self_script_lf_sha256={lf_sha256(Path(__file__))}")
    print(
        "dependency_lf_sha256="
        f"support:{lf_sha256(SUPPORT_PATH)},"
        f"fixed_geometry:{lf_sha256(FIXED_GEOMETRY_PATH)},"
        f"fixed_geometry_output:{lf_sha256(FIXED_GEOMETRY_OUTPUT_PATH)},"
        f"target_composition:{lf_sha256(TARGET_COMPOSITION_PATH)},"
        f"target_composition_output:{lf_sha256(TARGET_COMPOSITION_OUTPUT_PATH)},"
        f"target_one_spike:{lf_sha256(TARGET_ONE_SPIKE_PATH)}"
    )
    print(
        "hostile_hash_basis_control="
        f"actual_target_LF:{lf_sha256(TARGET_COMPOSITION_PATH)},"
        f"advertised:{ADVERTISED_TARGET_COMPOSITION_SHA256},"
        f"escaped_literal_rewrite:{escaped_literal_sha256(TARGET_COMPOSITION_PATH)};"
        "verdict=advertised_value_rewrites_literal_backslash_r_backslash_n_not_line_endings"
    )
    print(
        f"fixed_geometry_scope=clocks_1_through_14,max_cover_size={COVER_BUDGET},strict_endpoints=True;"
        f"support_rows={sum(map(len, geometry['by_divisor'].values()))};"
        f"support_divisors={len(geometry['by_divisor'])};"
        f"distinct_targets={geometry['distinct_targets']};"
        f"terminal_rows={len(geometry['terminal_keys'])};"
        f"terminal_histogram={tuple(sorted(geometry['terminal_histogram'].items()))};"
        f"terminal_semantic_sha256={geometry['terminal_semantic']}"
    )
    print(
        f"independent_ledger_engine=literal_denominator_symbol_current_lcm_DP_plus_primal_Boolean_coupling_vertices;"
        f"literal_distribution_controls=({literal_control_cases},{literal_control_shapes});"
        f"coupling_marginal_types={coupling_vertices.cache_info().currsize};"
        f"status_profiles={ledger['primal_status_profiles']};"
        f"distribution_semantic_sha256={ledger['distribution_semantic']}"
    )
    print(f"pre_summary={before['summary']};pre_occurrence_by_c={tuple(sorted(before['by_c'].items()))};pre_shape_by_c={tuple(sorted(before['shape_by_c'].items()))};pre_semantic_sha256={before['semantic']}")
    print(f"exact_key=(body,D);refined_closed_rows={len(correct_hits)};refined_closed_occurrences={killed_occurrences};hit_semantic_sha256={hit_semantic.hexdigest()}")
    print(f"refined_closed_detail={hit_details}")
    print(f"post_summary={after['summary']};post_occurrence_by_c={tuple(sorted(after['by_c'].items()))};post_shape_by_c={tuple(sorted(after['shape_by_c'].items()))};post_semantic_sha256={after['semantic']}")
    print(f"hostile_wrong_row_keys={wrong_key_controls};correct=(7, 7648);verdict=PASS_wrong_keys_do_not_reproduce_exact_intersection")
    print(
        f"hostile_overlap_subtraction=raw_terminal_occurrences:{raw_terminal_occurrences},"
        f"raw_hit_occurrences:{raw_hit_occurrences},"
        f"naive_global_post:{naive_global_subtraction},"
        f"naive_hit_post:{naive_hit_subtraction},"
        f"double_subtraction_post:{double_subtraction},"
        f"exclusive_shape_features:{len(exclusive_features)},"
        f"correct_post:{after['summary'][1]};"
        "verdict=PASS_only_refined_row_counts_on_exact_intersection_may_be_deleted"
    )
    print(
        f"hostile_scope_control=first_pool14_unresolved_refined_row:{scope_witness};"
        f"unresolved_refined_rows:{len(scope_survivors)};"
        "meaning=no_cover_by_at_most_5_clocks_from_1..14_only;"
        "not_a_claim_of_physical_realizability_or_absolute_noncoverability_or_LRC14;verdict=PASS"
    )
    print("discrepancy=NONE")
    print("verdict=PASS")


if __name__ == "__main__":
    main()

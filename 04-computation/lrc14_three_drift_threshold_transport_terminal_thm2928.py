#!/usr/bin/env python3
"""Exact terminal referee for the THM-2928 mixed three-drift branch.

The frozen activity/fibre referee leaves 29,364 occurrences.  This script
reconstructs that residual from the original 544,571 all-top occurrences
and then applies two further necessary conditions.

Threshold transport.  Fix an outer quotient t|D and put M=D/t.  For a
denominator d, let ell=ceil(d/7), g=gcd(d,t), and ell=A*g+r.  On every
t-fibre the enlarged d-needle is a lifted unit arithmetic progression in
Z/MZ of denominator d/g and length A or A+1.  The high length occurs on
exactly

    R=(t/g)r

fibres.  Thus the eight joint high/low status counts n_E obey

    sum_E n_E=t,             sum_{E contains i} n_E=R_i.

For a status E, c_E is the sum of the three individual trace sizes, capped
at M.  If lambda_b is the target load on fibre b, then for every integer k

    #{b: lambda_b >= k}
        <= max sum_{E: c_E >= k} n_E.                 (*)

The maximum is an eight-variable real LP.  Both its primal basic feasible
solutions and its dual vertices are enumerated with exact Fractions; their
optima must agree.  Every t|D and every distinct positive target load is
tested.  This leaves exactly 19 occurrences.

Terminal local cover.  One explicit fibre is then treated exactly in each
of the 19 rows.  Each actual trace is contained in a lifted unit AP of
width ceil(ell/gcd(d,q)), where q=D/M.  We relax further by granting every
needle an independent phase and unit direction on that fibre.  A complete
depth-three grouped set-cover recursion still finds no cover: 18 witnesses
use M=99 and the remaining witness uses M=98.

The terminal computation is numerator-free and strictly relaxed at both
stages.  Consequently rejection is a necessary obstruction, not a search
over assumed numerators.  All guards use RuntimeError and remain active
under ``python -O``.
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
COMBINED_PATH = (
    HERE / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
COMBINED_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_three_drift_body_projection_fiber_thm2928.out"
)
ACTIVITY_PATH = (
    HERE / "lrc14_three_drift_mixed_activity_fiber_thm2928.py"
)
ACTIVITY_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_three_drift_mixed_activity_fiber_thm2928.out"
)

EXPECTED_COMBINED_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
EXPECTED_COMBINED_OUTPUT_SHA256 = (
    "2e211620ad7064ea06f7544b5fbac709d6d52d9a0e261b464ae26b595f09b669"
)
EXPECTED_ACTIVITY_SHA256 = (
    "067424a0edb126ad8e15f9f56ad77489bfdbb9a5fd2dd42e8a3eb4599438dfd9"
)
EXPECTED_ACTIVITY_OUTPUT_SHA256 = (
    "cab8e7b4a63177e89f19979c3c8129360dc486786b5a594bc871b25c4bf5c723"
)
EXPECTED_ALL_TOP_SHA256 = (
    "88ac1b6c34c11eadd286d02f930f69a44a82bd429fe90248417fd35caa947960"
)
EXPECTED_Q_SURVIVOR_SHA256 = (
    "87fc8e09d7e69c9d91cf7ebc895703393496e9b1d4b5a1623ca6c112c6cbb4fc"
)
EXPECTED_ACTIVITY_RESIDUAL_SHA256 = (
    "6a0b57e9da265d9f5aab55b8dd6166e7439d759a0f4ba5f1ef6e43b4d212cab5"
)
EXPECTED_THRESHOLD_SEMANTIC_SHA256 = (
    "088b4e16a21a7cc0a177e0e1148834e2e097bca348739b9174cce46033eda3cb"
)
EXPECTED_THRESHOLD_SURVIVOR_SHA256 = (
    "ca19dda6a982440c56461267e7b4fc649d026dabeddf42ce9a2e2879a1a55543"
)
EXPECTED_TERMINAL_SEMANTIC_SHA256 = (
    "3aa5d8b3e267bdc4bf5b90f315d8456be98a8b2b87026e7f8713d140d2960a4c"
)

PATTERNS = tuple(range(8))
PATTERN_VECTORS = {
    pattern: (1,) + tuple((pattern >> index) & 1 for index in range(3))
    for pattern in PATTERNS
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


for path, expected, label in (
    (COMBINED_PATH, EXPECTED_COMBINED_SHA256, "combined script"),
    (
        COMBINED_OUTPUT_PATH,
        EXPECTED_COMBINED_OUTPUT_SHA256,
        "combined output",
    ),
    (ACTIVITY_PATH, EXPECTED_ACTIVITY_SHA256, "activity script"),
    (
        ACTIVITY_OUTPUT_PATH,
        EXPECTED_ACTIVITY_OUTPUT_SHA256,
        "activity output",
    ),
):
    require(file_sha256(path) == expected, f"{label} dependency changed")

spec = spec_from_file_location("lrc14_three_drift_combined", COMBINED_PATH)
combined = module_from_spec(spec)
spec.loader.exec_module(combined)
support = combined.support_module


def solve_matrix(matrix, rhs):
    """Solve a 4 by 4 rational system, returning None when singular."""
    augmented = [
        [Q(value) for value in row] + [Q(rhs[index])]
        for index, row in enumerate(matrix)
    ]
    for column in range(4):
        pivot = next(
            (
                row
                for row in range(column, 4)
                if augmented[row][column]
            ),
            None,
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        scale = augmented[column][column]
        augmented[column] = [
            value / scale for value in augmented[column]
        ]
        for row in range(4):
            if row == column or not augmented[row][column]:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                left - scale * right
                for left, right in zip(
                    augmented[row], augmented[column]
                )
            ]
    return tuple(augmented[row][4] for row in range(4))


def solve_columns(patterns, rhs):
    """Solve with the four selected pattern vectors as columns."""
    return solve_matrix(
        tuple(
            tuple(PATTERN_VECTORS[pattern][row] for pattern in patterns)
            for row in range(4)
        ),
        rhs,
    )


PRIMAL_BASES = tuple(
    basis
    for basis in combinations(PATTERNS, 4)
    if solve_columns(basis, (0, 0, 0, 0)) is not None
)
require(len(PRIMAL_BASES) == 58, "transport basis count changed")


@lru_cache(maxsize=None)
def dual_vertices(costs):
    """All vertices of A^T y >= costs for one binary objective."""
    vertices = set()
    for tight in PRIMAL_BASES:
        vertex = solve_matrix(
            tuple(PATTERN_VECTORS[pattern] for pattern in tight),
            tuple(costs[pattern] for pattern in tight),
        )
        require(vertex is not None, "nonsingular dual basis became singular")
        if all(
            sum(
                Q(coefficient) * value
                for coefficient, value in zip(
                    PATTERN_VECTORS[pattern], vertex
                )
            )
            >= costs[pattern]
            for pattern in PATTERNS
        ):
            vertices.add(vertex)
    require(vertices, ("dual has no enumerated vertex", costs))
    return tuple(sorted(vertices))


@lru_cache(maxsize=None)
def transport_optimum(total, marginals, costs):
    """Exact primal/dual optimum for the eight-status relaxation."""
    rhs = (total,) + marginals
    primal = None
    for basis in PRIMAL_BASES:
        values = solve_columns(basis, rhs)
        require(values is not None, "transport basis became singular")
        if any(value < 0 for value in values):
            continue
        objective = sum(
            value * costs[pattern]
            for pattern, value in zip(basis, values)
        )
        if primal is None or objective > primal:
            primal = objective
    require(primal is not None, ("empty transport primal", total, marginals))

    dual = min(
        sum(Q(value) * coefficient for value, coefficient in zip(rhs, vertex))
        for vertex in dual_vertices(costs)
    )
    require(
        primal == dual,
        ("transport primal/dual gap", total, marginals, costs, primal, dual),
    )
    return primal


@lru_cache(maxsize=None)
def status_data(D, ds, t):
    """Trace denominators, low lengths, marginals, and raw status caps."""
    require(D % t == 0, "outer quotient does not divide D")
    M = D // t
    inner_ds = []
    lows = []
    marginals = []
    heights = []
    for d in ds:
        ell = (d + 6) // 7
        common = gcd(d, t)
        low, remainder = divmod(ell, common)
        inner_d = d // common
        require(M % inner_d == 0, "trace denominator does not divide fibre")
        inner_ds.append(inner_d)
        lows.append(low)
        marginals.append((t // common) * remainder)
        heights.append(M // inner_d)
    capacities = tuple(
        min(
            M,
            sum(
                heights[index]
                * (lows[index] + ((pattern >> index) & 1))
                for index in range(3)
            ),
        )
        for pattern in PATTERNS
    )
    return M, tuple(inner_ds), tuple(lows), tuple(marginals), capacities


def fibre_cap(D, d, q):
    """Sharp individual maximum on a q-fibre, used in inheritance."""
    ell = (d + 6) // 7
    common = gcd(d, q)
    height = D // lcm(d, q)
    return height * ((ell + common - 1) // common)


def reconstruct_activity_residual():
    """Rebuild the frozen 544,571 -> 29,364 inheritance chain."""
    by_divisor = {}
    body_count = 0
    divisor_rows = 0
    support_hard_rows = 0
    for F in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support.safe_cell_ranges(F)
        for D in support.divisors(L):
            divisor_rows += 1
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > combined.SUPPORT_CUTOFF:
                continue
            support_hard_rows += 1
            by_divisor.setdefault(D, []).append(
                (support_count, F, L, tuple(ranges))
            )
    require(body_count == 3003, "body universe changed")
    require(divisor_rows == 251536, "body/divisor universe changed")
    require(support_hard_rows == 13778, "support-hard universe changed")
    require(len(by_divisor) == 206, "support-hard divisor alphabet changed")

    all_top_count = 0
    q_kills = 0
    q_survivors = 0
    activity_kills = 0
    double_two_kills = 0
    all_top_semantic = sha256()
    q_semantic = sha256()
    residual_semantic = sha256()
    records = []

    for D in sorted(by_divisor):
        divisors = tuple(d for d in support.divisors(D) if d > 1)
        capacities = tuple(
            (D // d) * ((d + 6) // 7) for d in divisors
        )
        rows = sorted(
            by_divisor[D],
            key=lambda record: (record[0], record[1], record[2]),
        )
        support_sizes = tuple(record[0] for record in rows)
        top_loads_by_row = []
        maxima_by_row = []
        parity_by_row = []
        for support_count, F, _L, ranges in rows:
            arcs = combined.projected_support_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("support arc mismatch", F, D),
            )
            tops = []
            maxima = []
            parity = None
            for q in divisors:
                histogram = combined.residue_load_histogram(arcs, q)
                tops.append(
                    combined.top_class_load(histogram, (q + 6) // 7)
                )
                maxima.append(max(load for load, count in histogram if count))
                if q == 2:
                    require(
                        histogram == ((support_count // 2, 2),),
                        ("reflection parity balance changed", F, D),
                    )
                    parity = support_count // 2
            top_loads_by_row.append(tuple(tops))
            maxima_by_row.append(tuple(maxima))
            parity_by_row.append(parity)

        for first_index, d1 in enumerate(divisors):
            for second_index in range(first_index, len(divisors)):
                d2 = divisors[second_index]
                pair_lcm = lcm(d1, d2)
                for third_index in range(second_index, len(divisors)):
                    d3 = divisors[third_index]
                    if lcm(pair_lcm, d3) != D:
                        continue
                    full_capacity = (
                        capacities[first_index]
                        + capacities[second_index]
                        + capacities[third_index]
                    )
                    admissible_rows = bisect_right(
                        support_sizes, full_capacity
                    )
                    for row_index in range(admissible_rows):
                        support_count, F, L, _ranges = rows[row_index]
                        top_loads = top_loads_by_row[row_index]
                        if support_count > (
                            top_loads[first_index]
                            + top_loads[second_index]
                            + top_loads[third_index]
                        ):
                            continue

                        all_top_count += 1
                        all_top_semantic.update(
                            (
                                f"{F}|{L}|{D}|{support_count}|"
                                f"{d1}|{d2}|{d3}\n"
                            ).encode()
                        )
                        ds = (d1, d2, d3)
                        fibre_caps = tuple(
                            sum(fibre_cap(D, d, q) for d in ds)
                            for q in divisors
                        )
                        margins = tuple(
                            maximum - cap
                            for maximum, cap in zip(
                                maxima_by_row[row_index], fibre_caps
                            )
                        )
                        best_margin = max(margins)
                        if best_margin > 0:
                            q_kills += 1
                            continue

                        q_survivors += 1
                        record = (
                            F,
                            L,
                            D,
                            support_count,
                            d1,
                            d2,
                            d3,
                            best_margin,
                        )
                        q_semantic.update(f"{record}\n".encode())
                        global_caps = (
                            capacities[first_index],
                            capacities[second_index],
                            capacities[third_index],
                        )
                        if any(
                            d in (2, 3)
                            and support_count
                            > sum(
                                global_caps[other]
                                for other in range(3)
                                if other != index
                            )
                            for index, d in enumerate(ds)
                        ):
                            activity_kills += 1
                            continue

                        if d1 == d2 == 2:
                            parity = parity_by_row[row_index]
                            require(
                                parity is not None,
                                ("double-two row lacks parity", F, D),
                            )
                            third_cap = fibre_cap(D, d3, 2)
                            require(
                                parity > third_cap,
                                ("double-two parity obstruction lost", record),
                            )
                            double_two_kills += 1
                            continue

                        records.append(record)
                        residual_semantic.update(f"{record}\n".encode())

    require(all_top_count == 544571, "all-top occurrence count changed")
    require(q_kills == 125060, "q-fibre kill count changed")
    require(q_survivors == 419511, "q-fibre survivor count changed")
    require(activity_kills == 383391, "activity kill count changed")
    require(double_two_kills == 6756, "double-two kill count changed")
    require(len(records) == 29364, "activity residual count changed")
    require(
        all_top_semantic.hexdigest() == EXPECTED_ALL_TOP_SHA256,
        "all-top semantic atlas changed",
    )
    require(
        q_semantic.hexdigest() == EXPECTED_Q_SURVIVOR_SHA256,
        "q-fibre survivor atlas changed",
    )
    require(
        residual_semantic.hexdigest() == EXPECTED_ACTIVITY_RESIDUAL_SHA256,
        "activity residual atlas changed",
    )
    return records, by_divisor


def threshold_transport(records, by_divisor):
    """Apply (*) on every divisor quotient and return the 19 survivors."""
    ranges_by_body = {}
    arcs_cache = {}
    histogram_cache = {}
    kills = []
    survivors = []
    decisive_t = Counter()
    decisive_patterns = Counter()
    semantic = sha256()
    survivor_semantic = sha256()

    def arcs_for(F, L, D):
        key = (F, D)
        if key not in arcs_cache:
            if F not in ranges_by_body:
                actual_L, ranges = support.safe_cell_ranges(F)
                require(actual_L == L, ("body ruler changed", F, L, actual_L))
                ranges_by_body[F] = tuple(ranges)
            arcs_cache[key] = combined.projected_support_arcs(
                D, ranges_by_body[F]
            )
        return arcs_cache[key]

    for record in records:
        F, L, D, _support_count, d1, d2, d3, _margin = record
        ds = (d1, d2, d3)
        witness = None
        for t in support.divisors(D):
            if t == 1:
                continue
            histogram_key = (F, D, t)
            if histogram_key not in histogram_cache:
                histogram_cache[histogram_key] = (
                    combined.residue_load_histogram(arcs_for(F, L, D), t)
                )
            histogram = histogram_cache[histogram_key]
            M, inner_ds, lows, marginals, capacities = status_data(
                D, ds, t
            )
            target_count = 0
            for threshold, count in reversed(histogram):
                target_count += count
                if threshold <= 0:
                    continue
                admissible = tuple(
                    int(capacity >= threshold) for capacity in capacities
                )
                maximum = transport_optimum(t, marginals, admissible)
                if Q(target_count) > maximum:
                    witness = (
                        t,
                        M,
                        threshold,
                        target_count,
                        (maximum.numerator, maximum.denominator),
                        inner_ds,
                        lows,
                        marginals,
                        capacities,
                        admissible,
                    )
                    break
            if witness is not None:
                break

        if witness is None:
            survivors.append(record)
            semantic.update(f"{record}|survive\n".encode())
            survivor_semantic.update(f"{record}\n".encode())
        else:
            kills.append((record, witness))
            semantic.update(f"{record}|{witness}\n".encode())
            decisive_t[witness[0]] += 1
            decisive_patterns[witness[-1]] += 1

    require(len(kills) == 29345, "threshold kill count changed")
    require(len(survivors) == 19, "threshold survivor count changed")
    require(
        semantic.hexdigest() == EXPECTED_THRESHOLD_SEMANTIC_SHA256,
        "threshold semantic atlas changed",
    )
    require(
        survivor_semantic.hexdigest() == EXPECTED_THRESHOLD_SURVIVOR_SHA256,
        "threshold survivor atlas changed",
    )
    return (
        kills,
        survivors,
        decisive_t,
        decisive_patterns,
        semantic.hexdigest(),
        survivor_semantic.hexdigest(),
        arcs_for,
    )


@lru_cache(maxsize=None)
def base_masks(n, width):
    """All cyclic unit-step AP masks of a given width in Z/nZ."""
    if width <= 0:
        return (0,)
    if width >= n:
        return ((1 << n) - 1,)
    masks = set()
    for step in range(1, n):
        if gcd(step, n) != 1:
            continue
        for start in range(n):
            masks.add(
                sum(
                    1 << ((start + index * step) % n)
                    for index in range(width)
                )
            )
    return tuple(sorted(masks))


@lru_cache(maxsize=None)
def family_masks(M, n, width):
    """Lift every base AP mask from Z/nZ to Z/MZ."""
    require(M % n == 0, ("local denominator does not divide M", M, n))
    repeat = ((1 << M) - 1) // ((1 << n) - 1)
    return tuple(mask * repeat for mask in base_masks(n, width))


@lru_cache(maxsize=None)
def masks_through(M, n, width, point):
    return tuple(
        mask
        for mask in family_masks(M, n, width)
        if (mask >> point) & 1
    )


@lru_cache(maxsize=None)
def coverable(target, family_keys):
    """Complete grouped set-cover recursion, one mask from each family."""
    if target == 0:
        return True
    if not family_keys:
        return False
    capacity = sum(
        family_masks(*key)[0].bit_count() for key in family_keys
    )
    if target.bit_count() > capacity:
        return False
    point = (target & -target).bit_length() - 1
    used_keys = set()
    for index, key in enumerate(family_keys):
        if key in used_keys:
            continue
        used_keys.add(key)
        remaining = family_keys[:index] + family_keys[index + 1 :]
        for mask in masks_through(*key, point):
            if coverable(target & ~mask, remaining):
                return True
    return False


def slice_mask(arcs, q, M, residue):
    """Bit mask of one q-fibre of a cyclic interval union in Z/(qM)Z."""
    starts = tuple(left for left, _right in arcs)
    result = 0
    for height in range(M):
        point = residue + q * height
        index = bisect_right(starts, point) - 1
        if index >= 0 and point < arcs[index][1]:
            result |= 1 << height
    return result


def terminal_local_cover(survivors, arcs_for):
    """Certify one noncoverable relaxed fibre in every surviving row."""
    terminal = []
    semantic = sha256()
    kill_moduli = Counter()
    for record in survivors:
        F, L, D, _support_count, d1, d2, d3, _margin = record
        M = 98 if D == 129360 else 99
        require(D % M == 0, ("terminal modulus does not divide D", record, M))
        q = D // M
        residue = 0
        target = slice_mask(arcs_for(F, L, D), q, M, residue)
        keys = []
        for d in (d1, d2, d3):
            common = gcd(d, q)
            n = d // common
            ell = (d + 6) // 7
            width = (ell + common - 1) // common
            keys.append((M, n, width))
        keys = tuple(sorted(keys))
        require(
            not coverable(target, keys),
            ("terminal relaxed fibre became coverable", record, M),
        )
        target_hash = sha256(
            target.to_bytes((M + 7) // 8, "little")
        ).hexdigest()
        family_summary = tuple(
            (n, width, len(family_masks(M, n, width)))
            for _M, n, width in keys
        )
        witness = (
            M,
            q,
            residue,
            target.bit_count(),
            tuple(point for point in range(M) if (target >> point) & 1),
            target_hash,
            family_summary,
        )
        terminal.append((record, witness))
        semantic.update(f"{record}|{witness}\n".encode())
        kill_moduli[M] += 1

    require(len(terminal) == 19, "terminal witness count changed")
    require(kill_moduli == Counter({99: 18, 98: 1}), "kill moduli changed")
    require(
        semantic.hexdigest() == EXPECTED_TERMINAL_SEMANTIC_SHA256,
        "terminal witness atlas changed",
    )
    return terminal, kill_moduli, semantic.hexdigest()


def controls():
    """Positive and hostile finite controls for both new mechanisms."""
    # Exhaust the trace-size and high-fibre marginal law at small moduli.
    trace_cases = 0
    for D in range(1, 25):
        for d in support.divisors(D):
            ell = (d + 6) // 7
            for t in support.divisors(D):
                common = gcd(d, t)
                low, remainder = divmod(ell, common)
                M = D // t
                n = d // common
                height = M // n
                high_count = (t // common) * remainder
                for step in range(1, d + 1):
                    if gcd(step, d) != 1:
                        continue
                    classes = {(step * index) % d for index in range(ell)}
                    for phase in range(d):
                        sizes = tuple(
                            sum(
                                ((b + t * y - phase) % d) in classes
                                for y in range(M)
                            )
                            for b in range(t)
                        )
                        allowed = {height * low}
                        if remainder:
                            allowed.add(height * (low + 1))
                        require(
                            all(size in allowed for size in sizes),
                            ("trace size law failed", D, d, t),
                        )
                        if remainder:
                            require(
                                sum(
                                    size == height * (low + 1)
                                    for size in sizes
                                )
                                == high_count,
                                ("trace marginal law failed", D, d, t),
                            )
                        trace_cases += 1

    # Four points on one outer fibre exceed three length-one traces.
    hostile_D = 28
    hostile_t = 4
    hostile_arcs = tuple((point, point + 1) for point in (0, 4, 8, 12))
    hostile_histogram = combined.residue_load_histogram(
        hostile_arcs, hostile_t
    )
    M, _inner_ds, _lows, marginals, capacities = status_data(
        hostile_D, (28, 28, 28), hostile_t
    )
    require(M == 7, "hostile inner modulus changed")
    hostile_target = sum(
        count for load, count in hostile_histogram if load >= 4
    )
    hostile_costs = tuple(int(capacity >= 4) for capacity in capacities)
    hostile_maximum = transport_optimum(
        hostile_t, marginals, hostile_costs
    )
    require(
        (hostile_target, hostile_maximum) == (1, 0),
        "hostile threshold control lost obstruction",
    )

    # A literal union from the three local families must be coverable.
    positive_keys = ((21, 3, 1), (21, 7, 1), (21, 21, 3))
    positive_target = (
        family_masks(*positive_keys[0])[0]
        | family_masks(*positive_keys[1])[1]
        | family_masks(*positive_keys[2])[2]
    )
    require(
        coverable(positive_target, positive_keys),
        "positive grouped-cover control rejected",
    )
    hostile_keys = ((7, 7, 1),) * 3
    hostile_target_mask = sum(1 << point for point in range(4))
    require(
        not coverable(hostile_target_mask, hostile_keys),
        "hostile grouped-cover control became coverable",
    )
    return (
        trace_cases,
        hostile_target,
        hostile_maximum,
        positive_target.bit_count(),
        hostile_target_mask.bit_count(),
    )


def main():
    (
        trace_cases,
        hostile_target,
        hostile_maximum,
        positive_local_size,
        hostile_local_size,
    ) = controls()
    records, by_divisor = reconstruct_activity_residual()
    (
        threshold_kills,
        threshold_survivors,
        decisive_t,
        decisive_patterns,
        threshold_semantic,
        survivor_semantic,
        arcs_for,
    ) = threshold_transport(records, by_divisor)
    terminal, kill_moduli, terminal_semantic = terminal_local_cover(
        threshold_survivors, arcs_for
    )

    print("LRC14 mixed three-drift threshold-transport terminal referee")
    print(f"combined_script_sha256={file_sha256(COMBINED_PATH)}")
    print(f"combined_output_sha256={file_sha256(COMBINED_OUTPUT_PATH)}")
    print(f"activity_script_sha256={file_sha256(ACTIVITY_PATH)}")
    print(f"activity_output_sha256={file_sha256(ACTIVITY_OUTPUT_PATH)}")
    print(f"transport_nonsingular_bases={len(PRIMAL_BASES)}")
    print(f"trace_control_cases={trace_cases}")
    print(
        "hostile_threshold_control="
        f"target_fibres={hostile_target},lp_maximum={hostile_maximum}"
    )
    print(f"positive_local_control_size={positive_local_size}")
    print(f"hostile_local_control_size={hostile_local_size}")
    print("inherited_all_top_occurrences=544571")
    print("inherited_activity_residual=29364")
    print(f"threshold_kills={len(threshold_kills)}")
    print(f"threshold_survivors={len(threshold_survivors)}")
    print(f"threshold_semantic_sha256={threshold_semantic}")
    print(f"threshold_survivor_sha256={survivor_semantic}")
    print(f"threshold_decisive_outer_top={decisive_t.most_common(80)}")
    print(
        "threshold_decisive_patterns="
        f"{decisive_patterns.most_common()}"
    )
    print(f"terminal_local_kills={len(terminal)}")
    print("terminal_local_survivors=0")
    print(f"terminal_kill_moduli={kill_moduli}")
    print(f"terminal_semantic_sha256={terminal_semantic}")
    print(f"terminal_witnesses={terminal}")
    print(f"transport_cache={transport_optimum.cache_info()}")
    print(f"dual_vertex_cache={dual_vertices.cache_info()}")
    print(f"base_mask_cache={base_masks.cache_info()}")
    print(f"family_mask_cache={family_masks.cache_info()}")
    print(f"coverable_cache={coverable.cache_info()}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()

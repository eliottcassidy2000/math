#!/usr/bin/env python3
"""Full k=2/p=5 expected-spike compression at q=D/7.

Every support-transfer row has 7|D.  Put q=D/7.  For a denominator d|D:

* if d does not divide q, its enlarged needle contributes one point to every
  q-fibre (a uniform-one mask);
* if d|q, let Y_d(u) be the number of q-fibres hit by the literal mask at
  common phase u.  Exact averaging gives integral Y_d du=q/7, independently
  of d and its coprime numerator.

If c of the five masks are uniform-one and m=5-c are spike-type, write
N_c=#{b mod q: lambda_q(b)>c}.  Pointwise coverage forces

    sum_i Y_i(u) >= N_c

throughout the compact aligned-safe carrier.  Markov and compact/open
strictness give the necessary condition

    N_c < 13*m*q/66

whenever N_c>0.

The exact-lcm denominator multiset is aggregated only by c.  A degree-five
divisor-Mobius generating function counts every shape without enumeration.
"""

from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
COMBINED_PATH = (
    ROOT / "04-computation" / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
EXPECTED_COMBINED_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
SUPPORT_CUTOFF = Q(887, 990)
TWO_SAFE_FLOOR = Q(66, 91)
EXPECTED_ROWS = 27163
EXPECTED_DIVISORS = 219
EXPECTED_SHAPES = 50874159718
EXPECTED_OCCURRENCES = 951545890235


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(
    file_sha256(COMBINED_PATH) == EXPECTED_COMBINED_SHA256,
    "frozen body-projection dependency changed",
)
spec = spec_from_file_location("lrc14_three_drift_combined", COMBINED_PATH)
combined = module_from_spec(spec)
spec.loader.exec_module(combined)
support = combined.support_module


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
    if remaining > 1:
        result = -result
    return result


def multichoose(types, copies):
    if copies == 0:
        return 1
    if types == 0:
        return 0
    return comb(types + copies - 1, copies)


def lcm_multiset_count(D):
    total = 0
    for E in support.divisors(D):
        alphabet = len(support.divisors(E)) - 1
        total += (
            mobius(D // E)
            * multichoose(alphabet, 5)
        )
    return total


def uniform_count_distribution(D):
    """Exact lcm-D five-multiset count by c=#{d not dividing D/7}."""
    require(D % 7 == 0, "q7 distribution requested off the septimal locus")
    q = D // 7
    result = Counter()
    for E in support.divisors(D):
        sign = mobius(D // E)
        if not sign:
            continue
        total_types = len(support.divisors(E)) - 1
        spike_types = len(support.divisors(gcd(E, q))) - 1
        uniform_types = total_types - spike_types
        require(
            uniform_types >= 0,
            ("negative uniform alphabet", D, E),
        )
        for c in range(6):
            result[c] += (
                sign
                * multichoose(uniform_types, c)
                * multichoose(spike_types, 5 - c)
            )
    require(
        all(value >= 0 for value in result.values()),
        ("negative exact-lcm c coefficient", D, result),
    )
    result += Counter()
    require(
        sum(result.values()) == lcm_multiset_count(D),
        ("c-distribution lost exact-lcm shapes", D),
    )
    return result


def brute_distribution_controls():
    cases = 0
    for D in range(7, 101, 7):
        q = D // 7
        alphabet = tuple(d for d in support.divisors(D) if d > 1)
        brute = Counter()
        for values in combinations_with_replacement(alphabet, 5):
            if lcm(*values) != D:
                continue
            c = sum(d not in support.divisors(q) for d in values)
            brute[c] += 1
        require(
            brute == uniform_count_distribution(D),
            ("literal c-distribution control failed", D, brute),
        )
        cases += 1
    return cases


def norm_fraction(value):
    residue = value % 1
    return min(residue, 1 - residue)


def hit_count(d, q, c, u):
    require(q % d == 0, "spike control denominator does not divide q")
    classes = sum(
        norm_fraction(Q(c) * (residue + u) / d) < Q(1, 14)
        for residue in range(d)
    )
    return (q // d) * classes


def exact_mean_controls():
    """Exact cell-midpoint integration for small spike denominators."""
    cases = 0
    for d in range(2, 21):
        q = lcm(7, d)
        for c in range(1, d + 1):
            if gcd(c, d) != 1:
                continue
            cell_count = 14 * c
            total = sum(
                hit_count(
                    d,
                    q,
                    c,
                    Q(2 * index + 1, 2 * cell_count),
                )
                for index in range(cell_count)
            )
            mean = Q(total, cell_count)
            require(
                mean == Q(q, 7),
                ("expected-spike law failed", d, q, c, mean),
            )
            cases += 1
    return cases


def main():
    distribution_controls = brute_distribution_controls()
    mean_controls = exact_mean_controls()

    rows_by_D = defaultdict(list)
    body_count = 0
    body_divisor_rows = 0
    for body in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support.safe_cell_ranges(body)
        for D in support.divisors(L):
            body_divisor_rows += 1
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > SUPPORT_CUTOFF:
                continue
            require(D % 7 == 0, ("k2 support row is not septimal", body, D))
            q = D // 7
            arcs = combined.projected_support_arcs(D, ranges)
            histogram = combined.residue_load_histogram(arcs, q)
            N = tuple(
                sum(count for load, count in histogram if load > c)
                for c in range(6)
            )
            rows_by_D[D].append((body, L, support_count, histogram, N))

    require(body_count == 3003, "body universe changed")
    require(body_divisor_rows == 251536, "body/divisor universe changed")
    require(sum(map(len, rows_by_D.values())) == EXPECTED_ROWS, "k2 row count changed")
    require(len(rows_by_D) == EXPECTED_DIVISORS, "k2 divisor alphabet changed")

    raw_shapes = 0
    raw_occurrences = 0
    surviving_shapes = 0
    surviving_occurrences = 0
    surviving_rows = set()
    surviving_bodies = set()
    surviving_divisors = set()
    shapes_by_c = Counter()
    occurrences_by_c = Counter()
    surviving_shapes_by_c = Counter()
    surviving_occurrences_by_c = Counter()
    killed_occurrences_by_c = Counter()
    strict_equalities = Counter()
    row_N_types = Counter()
    semantic = sha256()
    minimum_survivor = None
    minimum_kill = None

    for D in sorted(rows_by_D):
        q = D // 7
        distribution = uniform_count_distribution(D)
        raw_shapes += sum(distribution.values())
        for c, shape_count in distribution.items():
            shapes_by_c[c] += shape_count
            raw_occurrences += shape_count * len(rows_by_D[D])
            occurrences_by_c[c] += shape_count * len(rows_by_D[D])

            shape_survives = False
            for body, L, support_count, histogram, N in rows_by_D[D]:
                N_c = N[c]
                m = 5 - c
                if N_c == 0:
                    passes = True
                else:
                    left = 66 * N_c
                    right = 13 * m * q
                    passes = left < right
                    if left == right:
                        strict_equalities[(D, c, N_c)] += shape_count
                row_N_types[(D, c, q, N_c)] += 1
                record = (
                    D,
                    body,
                    L,
                    support_count,
                    c,
                    q,
                    N_c,
                    m,
                    shape_count,
                )
                if passes:
                    shape_survives = True
                    surviving_occurrences += shape_count
                    surviving_occurrences_by_c[c] += shape_count
                    surviving_rows.add((body, D))
                    surviving_bodies.add(body)
                    surviving_divisors.add(D)
                    semantic.update(f"{record}\n".encode())
                    if minimum_survivor is None or record < minimum_survivor:
                        minimum_survivor = record
                else:
                    killed_occurrences_by_c[c] += shape_count
                    if minimum_kill is None or record < minimum_kill:
                        minimum_kill = record
            if shape_survives:
                surviving_shapes += shape_count
                surviving_shapes_by_c[c] += shape_count

    require(raw_shapes == EXPECTED_SHAPES, "full raw shape universe changed")
    require(
        raw_occurrences == EXPECTED_OCCURRENCES,
        "full raw occurrence universe changed",
    )
    require(
        surviving_occurrences + sum(killed_occurrences_by_c.values())
        == raw_occurrences,
        "expected-spike screen lost occurrences",
    )

    print("LRC14 k=2/p=5 full q=D/7 expected-spike GF scout")
    print(f"combined_script_sha256={file_sha256(COMBINED_PATH)}")
    print(f"weighted_distribution_control_cases={distribution_controls}")
    print(f"exact_mean_control_cases={mean_controls}")
    print(f"support_rows={sum(map(len, rows_by_D.values()))}")
    print(f"support_divisors={len(rows_by_D)}")
    print("all_support_rows_have_7_dividing_D=PASS")
    print(f"raw_shapes={raw_shapes}")
    print(f"raw_occurrences={raw_occurrences}")
    print(f"shapes_by_c={shapes_by_c}")
    print(f"occurrences_by_c={occurrences_by_c}")
    print(f"surviving_shapes={surviving_shapes}")
    print(f"surviving_occurrences={surviving_occurrences}")
    print(f"surviving_rows={len(surviving_rows)}")
    print(f"surviving_bodies={len(surviving_bodies)}")
    print(f"surviving_divisors={len(surviving_divisors)}")
    print(f"surviving_shapes_by_c={surviving_shapes_by_c}")
    print(f"surviving_occurrences_by_c={surviving_occurrences_by_c}")
    print(f"killed_occurrences_by_c={killed_occurrences_by_c}")
    print(f"strict_equality_kills={strict_equalities}")
    print(f"row_N_type_count={len(row_N_types)}")
    print(f"row_N_types_top30={row_N_types.most_common(30)}")
    print(f"semantic_sha256={semantic.hexdigest()}")
    print(f"minimum_kill={minimum_kill}")
    print(f"minimum_survivor={minimum_survivor}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()

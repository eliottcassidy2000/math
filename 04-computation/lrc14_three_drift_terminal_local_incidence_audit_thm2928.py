#!/usr/bin/env python3
"""Independent audit of THM-2928's terminal 19 local obstructions.

The primary terminal referee reconstructs the full mixed three-drift
residual and decides its last 19 rows by a grouped first-uncovered-point
recursion on bit-replicated arithmetic-progression masks.  This audit uses
only the 19 record identities frozen in that referee's canonical output.
It independently

* reconstructs every literal body slice from the original integer ranges;
* constructs every lifted unit progression by literal residue membership;
* builds point-to-needle incidence bitsets; and
* decides three-family coverability through a separate one-/two-family
  pair-cover recursion.

The primary terminal masks and their bit-replication implementation are
never imported.  A bounded exhaustive control also checks directly that
literal global needle traces are contained in the claimed quotient
families.  All guards use RuntimeError and remain active under ``python -O``.
"""

from ast import literal_eval
from bisect import bisect_right
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
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
TERMINAL_PATH = (
    HERE / "lrc14_three_drift_threshold_transport_terminal_thm2928.py"
)
TERMINAL_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_three_drift_threshold_transport_terminal_thm2928.out"
)

EXPECTED_COMBINED_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
EXPECTED_COMBINED_OUTPUT_SHA256 = (
    "2e211620ad7064ea06f7544b5fbac709d6d52d9a0e261b464ae26b595f09b669"
)
EXPECTED_TERMINAL_SHA256 = (
    "13e524e728736480798d52acc736afe0cbd7b651a487aeef5d771d2b5dfa1338"
)
EXPECTED_TERMINAL_OUTPUT_SHA256 = (
    "435c34b249255c8659e62e89e3a22e42b2fa8dd67c4cbd4b188d5e9806206f52"
)
EXPECTED_PRIMARY_SURVIVOR_SHA256 = (
    "ca19dda6a982440c56461267e7b4fc649d026dabeddf42ce9a2e2879a1a55543"
)
EXPECTED_PRIMARY_WITNESS_SHA256 = (
    "3aa5d8b3e267bdc4bf5b90f315d8456be98a8b2b87026e7f8713d140d2960a4c"
)
EXPECTED_AUDIT_WITNESS_SHA256 = (
    "107548c4bc078c68394807524d15d519e663d6ebd712b097a15345fa3b3eb49e"
)


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
    (TERMINAL_PATH, EXPECTED_TERMINAL_SHA256, "terminal script"),
    (
        TERMINAL_OUTPUT_PATH,
        EXPECTED_TERMINAL_OUTPUT_SHA256,
        "terminal output",
    ),
):
    require(file_sha256(path) == expected, f"{label} dependency changed")

spec = spec_from_file_location("lrc14_three_drift_combined", COMBINED_PATH)
combined = module_from_spec(spec)
spec.loader.exec_module(combined)
support = combined.support_module


def divisors(number):
    result = []
    candidate = 1
    while candidate * candidate <= number:
        if number % candidate == 0:
            result.append(candidate)
            if candidate * candidate != number:
                result.append(number // candidate)
        candidate += 1
    return tuple(sorted(result))


def primary_terminal_records():
    """Parse and audit the terminal identities in the canonical output."""
    lines = TERMINAL_OUTPUT_PATH.read_text().splitlines()

    def unique_value(prefix):
        values = [
            line[len(prefix) :]
            for line in lines
            if line.startswith(prefix)
        ]
        require(len(values) == 1, ("terminal output field changed", prefix))
        return values[0]

    require(
        unique_value("threshold_survivor_sha256=")
        == EXPECTED_PRIMARY_SURVIVOR_SHA256,
        "primary threshold survivor hash changed",
    )
    require(
        unique_value("terminal_semantic_sha256=")
        == EXPECTED_PRIMARY_WITNESS_SHA256,
        "primary terminal witness hash changed",
    )
    require(
        unique_value("terminal_local_kills=") == "19",
        "primary terminal kill count changed",
    )
    require(
        unique_value("terminal_local_survivors=") == "0",
        "primary terminal survivor count changed",
    )
    payload = literal_eval(unique_value("terminal_witnesses="))
    require(len(payload) == 19, "primary terminal payload changed")
    records = tuple(record for record, _witness in payload)
    require(len(set(records)) == 19, "primary terminal identities repeated")
    semantic = sha256()
    for record in records:
        semantic.update(f"{record}\n".encode())
    require(
        semantic.hexdigest() == EXPECTED_PRIMARY_SURVIVOR_SHA256,
        "parsed terminal record ledger changed",
    )
    return tuple(payload)


@lru_cache(maxsize=None)
def literal_family(M, n, width):
    """Every lifted unit AP, built by direct residue membership."""
    require(M % n == 0, ("trace period does not divide fibre", M, n))
    require(1 <= width <= n, ("trace width outside period", M, n, width))
    if width == n:
        return ((1 << M) - 1,)
    masks = set()
    for step in range(1, n):
        if gcd(step, n) != 1:
            continue
        for phase in range(n):
            classes = {
                (phase + index * step) % n for index in range(width)
            }
            mask = 0
            for point in range(M):
                if point % n in classes:
                    mask |= 1 << point
            masks.add(mask)
    result = tuple(sorted(masks))
    require(result, ("empty literal family", M, n, width))
    expected_size = (M // n) * width
    require(
        {mask.bit_count() for mask in result} == {expected_size},
        ("literal family cardinality changed", M, n, width),
    )
    return result


@lru_cache(maxsize=None)
def incidence(M, n, width):
    """For each point, the bitset of family members containing it."""
    masks = literal_family(M, n, width)
    rows = []
    for point in range(M):
        word = 0
        for index, mask in enumerate(masks):
            if (mask >> point) & 1:
                word |= 1 << index
        rows.append(word)
    return tuple(rows)


@lru_cache(maxsize=None)
def restricted_family(target, key):
    """Distinct intersections of one needle family with a fixed target."""
    return tuple(sorted({target & mask for mask in literal_family(*key)}))


@lru_cache(maxsize=None)
def single_cover_maximum(target, key):
    """Exact maximum number of target points met by one family member."""
    return max(mask.bit_count() for mask in restricted_family(target, key))


@lru_cache(maxsize=None)
def pair_cover_maximum(target, key):
    """Exact maximum target coverage by two members of one family."""
    masks = restricted_family(target, key)
    return max(
        (left | right).bit_count()
        for left_index, left in enumerate(masks)
        for right in masks[left_index:]
    )


@lru_cache(maxsize=None)
def triple_cover_maximum(target, key):
    """Exact maximum target coverage by three members of one family.

    All distinct two-mask unions are enumerated.  Process them in decreasing
    size and scan a third mask only while the pair size plus the global
    one-mask maximum can improve the incumbent.  The stopping inequality is
    an exact upper bound, so this is exhaustive rather than heuristic.
    """
    masks = restricted_family(target, key)
    pair_unions = {
        left | right
        for left_index, left in enumerate(masks)
        for right in masks[left_index:]
    }
    one_maximum = single_cover_maximum(target, key)
    best = max(union.bit_count() for union in pair_unions)
    for union in sorted(
        pair_unions,
        key=lambda mask: (-mask.bit_count(), mask),
    ):
        pair_size = union.bit_count()
        if pair_size + one_maximum <= best:
            break
        candidate = max((union | mask).bit_count() for mask in masks)
        if candidate > best:
            best = candidate
    return best


def first_point(mask):
    return (mask & -mask).bit_length() - 1


def indices(word):
    while word:
        bit = word & -word
        yield bit.bit_length() - 1
        word -= bit


@lru_cache(maxsize=None)
def one_coverable(target, key):
    """Whether one member of ``key`` covers the target."""
    if target == 0:
        return True
    candidates = (1 << len(literal_family(*key))) - 1
    remaining = target
    rows = incidence(*key)
    while remaining and candidates:
        point = first_point(remaining)
        candidates &= rows[point]
        remaining &= remaining - 1
    return bool(candidates)


@lru_cache(maxsize=None)
def two_coverable(target, left_key, right_key):
    """Complete pair-cover test, independent of grouped depth recursion."""
    if target == 0:
        return True
    point = first_point(target)
    left_masks = literal_family(*left_key)
    for index in indices(incidence(*left_key)[point]):
        if one_coverable(target & ~left_masks[index], right_key):
            return True
    # The first uncovered point may instead be assigned to the right mask.
    right_masks = literal_family(*right_key)
    for index in indices(incidence(*right_key)[point]):
        if one_coverable(target & ~right_masks[index], left_key):
            return True
    return False


@lru_cache(maxsize=None)
def three_coverable(target, keys):
    """Complete three-family decision reduced to the pair-cover oracle."""
    if target == 0:
        return True
    point = first_point(target)
    used_keys = set()
    for family_index, key in enumerate(keys):
        if key in used_keys:
            continue
        used_keys.add(key)
        masks = literal_family(*key)
        other_keys = keys[:family_index] + keys[family_index + 1 :]
        for index in indices(incidence(*key)[point]):
            if two_coverable(target & ~masks[index], *other_keys):
                return True
    return False


def direct_slice(ranges, D, q, M, residue):
    """Reconstruct a literal full-ruler support slice point by point."""
    require(q * M == D, ("bad fibre factorization", D, q, M))
    starts = tuple(left for left, _right in ranges)
    result = 0
    for height in range(M):
        point = residue + q * height
        index = bisect_right(starts, point) - 1
        if index >= 0 and point < ranges[index][1]:
            result |= 1 << height
    return result


def local_keys(D, denominators, q, M):
    """Return the three maximal quotient-needle family keys."""
    require(q * M == D, ("bad local quotient", D, q, M))
    keys = []
    for d in denominators:
        common = gcd(d, q)
        n = d // common
        ell = (d + 6) // 7
        width = (ell + common - 1) // common
        require(M % n == 0, ("trace denominator escaped fibre", D, d, q))
        require(
            width == (n + 6) // 7,
            ("toothpick width identity failed", D, d, q, n, width),
        )
        keys.append((M, n, width))
    return tuple(sorted(keys))


def controls():
    """Positive/hostile solver checks and literal trace containment."""
    trace_cases = 0
    for D in range(2, 29):
        for d in divisors(D):
            if d < 2:
                continue
            ell = (d + 6) // 7
            for q in divisors(D):
                M = D // q
                common = gcd(d, q)
                n = d // common
                width = (ell + common - 1) // common
                family = literal_family(M, n, width)
                for step in range(1, d):
                    if gcd(step, d) != 1:
                        continue
                    for phase in range(d):
                        classes = {
                            (phase + index * step) % d
                            for index in range(ell)
                        }
                        for residue in range(q):
                            trace = sum(
                                1 << height
                                for height in range(M)
                                if (residue + q * height) % d in classes
                            )
                            require(
                                any(trace & ~mask == 0 for mask in family),
                                (
                                    "literal trace escaped local family",
                                    D,
                                    d,
                                    q,
                                    phase,
                                    step,
                                    residue,
                                ),
                            )
                            trace_cases += 1

    positive_keys = ((21, 3, 1), (21, 7, 1), (21, 21, 3))
    positive_target = (
        literal_family(*positive_keys[0])[0]
        | literal_family(*positive_keys[1])[1]
        | literal_family(*positive_keys[2])[2]
    )
    require(
        three_coverable(positive_target, positive_keys),
        "positive literal union rejected",
    )
    hostile_keys = ((7, 7, 1),) * 3
    hostile_target = sum(1 << point for point in range(4))
    require(
        not three_coverable(hostile_target, hostile_keys),
        "hostile four-point target became coverable",
    )
    return trace_cases, positive_target.bit_count(), hostile_target.bit_count()


def main():
    trace_cases, positive_size, hostile_size = controls()
    primary = primary_terminal_records()
    expected_family_counts = {
        (99, 11, 2): 55,
        (99, 33, 5): 330,
        (99, 99, 15): 2970,
        (98, 49, 7): 1029,
        (98, 98, 14): 2058,
    }
    for key, expected in expected_family_counts.items():
        require(
            len(literal_family(*key)) == expected,
            ("terminal family count changed", key, expected),
        )

    witnesses = []
    kill_moduli = Counter()
    target_sizes = Counter()
    profile_counts = Counter()
    mechanism_counts = Counter()
    full_period_triple_maxima = {}
    crt_grid_omission_count = 0
    for record, primary_witness in primary:
        F, L, D, support_count, d1, d2, d3, margin = record
        require(D == L, ("terminal row left full ruler", record))
        require(margin == 0, ("terminal row margin changed", record))
        actual_L, ranges = support.safe_cell_ranges(F)
        require(actual_L == L, ("literal body ruler changed", F, L))
        ranges = tuple(ranges)
        require(
            sum(right - left for left, right in ranges) == support_count,
            ("literal support count changed", record),
        )

        (
            M,
            q,
            residue,
            primary_size,
            primary_points,
            primary_hash,
            primary_families,
        ) = primary_witness
        target = direct_slice(ranges, D, q, M, residue)
        points = tuple(
            point for point in range(M) if (target >> point) & 1
        )
        target_hash = sha256(
            target.to_bytes((M + 7) // 8, "little")
        ).hexdigest()
        require(
            (target.bit_count(), points, target_hash)
            == (primary_size, primary_points, primary_hash),
            ("direct slice disagrees with primary witness", record),
        )

        keys = local_keys(D, (d1, d2, d3), q, M)
        family_summary = tuple(
            (n, width, len(literal_family(M, n, width)))
            for _M, n, width in keys
        )
        require(
            family_summary == primary_families,
            ("literal families disagree with primary witness", record),
        )
        require(
            not three_coverable(target, keys),
            ("terminal target became independently coverable", record),
        )

        # Quantitative orbit-cover certificate.  Every terminal profile has
        # two copies of the full-period family.  First use exact individual
        # target intersections; then retain the exact union capacity of the
        # repeated pair.  Only four all-full profiles survive both bounds.
        single_maxima = tuple(
            single_cover_maximum(target, key) for key in keys
        )
        single_sum = sum(single_maxima)
        full_key = (M, M, (M + 6) // 7)
        require(keys.count(full_key) >= 2, ("full pair disappeared", record))
        repeated_pair_maximum = pair_cover_maximum(target, full_key)
        third_index = next(
            (
                index
                for index, key in enumerate(keys)
                if key != full_key
            ),
            None,
        )
        third_maximum = (
            single_cover_maximum(target, keys[third_index])
            if third_index is not None
            else single_cover_maximum(target, full_key)
        )
        pair_plus_third = repeated_pair_maximum + third_maximum
        exact_triple_maximum = None
        if single_sum < target.bit_count():
            mechanism = "single-family-sum"
        elif pair_plus_third < target.bit_count():
            mechanism = "repeated-pair-plus-third"
        else:
            require(
                keys == (full_key,) * 3,
                ("unexpected quantitative residual profile", record),
            )
            exact_triple_maximum = triple_cover_maximum(target, full_key)
            require(
                exact_triple_maximum < target.bit_count(),
                ("all-full exact maximum lost obstruction", record),
            )
            mechanism = "all-full-exact"
            full_period_triple_maxima[F] = exact_triple_maximum
        mechanism_counts[mechanism] += 1

        if M == 99:
            require(
                9 in F and 11 in F,
                ("Z/99 terminal body lost CRT wall speeds", F),
            )
            require(
                all(point % 9 and point % 11 for point in points),
                ("Z/99 target met forbidden CRT row or column", record),
            )
            crt_grid_omission_count += 1

        quantitative_certificate = (
            mechanism,
            single_maxima,
            repeated_pair_maximum,
            pair_plus_third,
            exact_triple_maximum,
        )
        witness = (
            record,
            M,
            q,
            residue,
            target.bit_count(),
            keys,
            target_hash,
            quantitative_certificate,
        )
        witnesses.append(witness)
        kill_moduli[M] += 1
        target_sizes[target.bit_count()] += 1
        profile_counts[(F, keys)] += 1

    require(len(witnesses) == 19, "independent witness count changed")
    require(kill_moduli == Counter({99: 18, 98: 1}), "kill moduli changed")
    require(
        sorted(target_sizes.items()) == [(28, 10), (30, 1), (32, 4), (34, 4)],
        "target-size census changed",
    )
    require(len(profile_counts) == 9, "terminal profile count changed")
    require(
        mechanism_counts
        == Counter(
            {
                "single-family-sum": 5,
                "repeated-pair-plus-third": 10,
                "all-full-exact": 4,
            }
        ),
        "quantitative obstruction ledger changed",
    )
    require(
        full_period_triple_maxima
        == {
            (1, 6, 7, 8, 9, 11): 27,
            (1, 5, 7, 8, 9, 11): 24,
        },
        "all-full exact maxima changed",
    )
    require(
        crt_grid_omission_count == 18,
        "Z/99 CRT grid omission count changed",
    )

    semantic = sha256()
    for witness in witnesses:
        semantic.update(f"{witness}\n".encode())
    require(
        semantic.hexdigest() == EXPECTED_AUDIT_WITNESS_SHA256,
        "independent terminal witness atlas changed",
    )

    print("LRC14 terminal-19 independent literal-incidence audit")
    print(f"combined_script_sha256={file_sha256(COMBINED_PATH)}")
    print(f"combined_output_sha256={file_sha256(COMBINED_OUTPUT_PATH)}")
    print(f"terminal_script_sha256={file_sha256(TERMINAL_PATH)}")
    print(f"terminal_output_sha256={file_sha256(TERMINAL_OUTPUT_PATH)}")
    print(f"literal_trace_control_cases={trace_cases}")
    print(f"positive_pair_solver_control_size={positive_size}")
    print(f"hostile_pair_solver_control_size={hostile_size}")
    print(f"primary_terminal_records={len(primary)}")
    print(f"independent_local_kills={len(witnesses)}")
    print("independent_local_survivors=0")
    print(f"kill_moduli={kill_moduli}")
    print(f"target_size_counts={sorted(target_sizes.items())}")
    print(f"distinct_terminal_profiles={len(profile_counts)}")
    print(
        "terminal_profile_counts="
        f"{sorted(profile_counts.items(), key=repr)}"
    )
    print(f"quantitative_mechanisms={mechanism_counts}")
    print(
        "full_period_triple_maxima="
        f"{sorted(full_period_triple_maxima.items())}"
    )
    print(f"z99_row9_column11_omissions={crt_grid_omission_count}")
    print(f"family_counts={sorted(expected_family_counts.items())}")
    print(f"audit_witness_sha256={semantic.hexdigest()}")
    print(f"literal_family_cache={literal_family.cache_info()}")
    print(f"one_coverable_cache={one_coverable.cache_info()}")
    print(f"two_coverable_cache={two_coverable.cache_info()}")
    print(f"three_coverable_cache={three_coverable.cache_info()}")
    print(f"restricted_family_cache={restricted_family.cache_info()}")
    print(f"single_cover_maximum_cache={single_cover_maximum.cache_info()}")
    print(f"pair_cover_maximum_cache={pair_cover_maximum.cache_info()}")
    print(f"triple_cover_maximum_cache={triple_cover_maximum.cache_info()}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()

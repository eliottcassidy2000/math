#!/usr/bin/env python3
"""Independent audit candidate for the THM-2928 three-drift terminal.

This referee is deliberately smaller than the primary census.  It audits the
two generic implications used after the 29,364-row inheritance boundary and
then independently checks the terminal 19-row local obstructions.

It does not read a scratch pickle or the Lorenz residual.  The 19 identities
are embedded below and tied to the primary referee by their semantic SHA-256.
For local coverability, body slices are rebuilt directly from the integer
safe-cell ranges, lifted AP families are built by literal residue tests, and a
pair-incidence solver is used instead of the primary grouped recursion.

The generic controls are:

* literal high/low fibre loads and maximal-local-AP containment for every
  D <= 84, every d,t | D, every unit direction and phase, and every t-fibre;
* the primary exact real transport LP versus an independent labelled-fibre
  integer dynamic program on every monotone three-bit threshold word through
  t=4 and deterministic boundary/interior samples at t=5,7,9;
* literal positive unions, a hostile non-cover, and the inclusive
  ``capacity >= threshold`` endpoint convention.

Set ``LRC14_AUDIT_ROOT`` only while developing this scratch candidate against
another worktree.  A promoted copy should use its own repository root and pin
the final primary source/output hashes.
"""

from functools import lru_cache
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import gcd, lcm
from pathlib import Path
import os


HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = HERE.parents[1]
ROOT = Path(os.environ.get("LRC14_AUDIT_ROOT", DEFAULT_ROOT)).resolve()
PRIMARY_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_three_drift_threshold_transport_terminal_thm2928.py"
)
EXPECTED_SURVIVOR_SEMANTIC_SHA256 = (
    "ca19dda6a982440c56461267e7b4fc649d026dabeddf42ce9a2e2879a1a55543"
)


TERMINAL_RECORDS = (
    (
        (1, 4, 7, 9, 11, 12),
        38808,
        38808,
        13248,
        9702,
        12936,
        19404,
        0,
    ),
    (
        (1, 4, 7, 9, 11, 12),
        38808,
        38808,
        13248,
        9702,
        19404,
        38808,
        0,
    ),
    (
        (1, 6, 7, 8, 9, 11),
        77616,
        77616,
        25390,
        19404,
        25872,
        38808,
        0,
    ),
    (
        (1, 6, 7, 8, 9, 11),
        77616,
        77616,
        25390,
        19404,
        38808,
        77616,
        0,
    ),
    (
        (1, 6, 7, 8, 9, 11),
        77616,
        77616,
        25390,
        25872,
        38808,
        38808,
        0,
    ),
    (
        (1, 4, 7, 8, 9, 11),
        77616,
        77616,
        26370,
        25872,
        38808,
        38808,
        0,
    ),
    (
        (1, 6, 7, 8, 9, 11),
        77616,
        77616,
        25390,
        25872,
        38808,
        77616,
        0,
    ),
    (
        (1, 4, 7, 8, 9, 11),
        77616,
        77616,
        26370,
        25872,
        38808,
        77616,
        0,
    ),
    (
        (1, 3, 7, 8, 10, 11),
        129360,
        129360,
        39274,
        25872,
        32340,
        129360,
        0,
    ),
    (
        (1, 5, 7, 8, 9, 11),
        388080,
        388080,
        109044,
        25872,
        97020,
        388080,
        0,
    ),
    (
        (1, 5, 7, 8, 9, 11),
        388080,
        388080,
        109044,
        25872,
        194040,
        388080,
        0,
    ),
    (
        (1, 5, 7, 8, 9, 11),
        388080,
        388080,
        109044,
        25872,
        388080,
        388080,
        0,
    ),
    (
        (1, 5, 7, 8, 9, 11),
        388080,
        388080,
        109044,
        38808,
        97020,
        388080,
        0,
    ),
    (
        (1, 5, 7, 8, 9, 11),
        388080,
        388080,
        109044,
        38808,
        129360,
        388080,
        0,
    ),
    (
        (1, 5, 7, 8, 9, 11),
        388080,
        388080,
        109044,
        43120,
        77616,
        388080,
        0,
    ),
    (
        (1, 5, 7, 8, 9, 11),
        388080,
        388080,
        109044,
        77616,
        97020,
        194040,
        0,
    ),
    (
        (1, 5, 7, 8, 9, 11),
        388080,
        388080,
        109044,
        77616,
        97020,
        388080,
        0,
    ),
    (
        (1, 5, 7, 8, 9, 11),
        388080,
        388080,
        109044,
        77616,
        129360,
        194040,
        0,
    ),
    (
        (1, 5, 7, 8, 9, 11),
        388080,
        388080,
        109044,
        77616,
        129360,
        388080,
        0,
    ),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def semantic_sha256(rows):
    digest = sha256()
    for row in rows:
        digest.update(f"{row}\n".encode())
    return digest.hexdigest()


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(PRIMARY_PATH.exists(), ("primary referee is missing", PRIMARY_PATH))
primary = load_module("lrc14_three_drift_terminal_primary_audited", PRIMARY_PATH)
support = primary.support
require(
    primary.EXPECTED_THRESHOLD_SURVIVOR_SHA256
    == EXPECTED_SURVIVOR_SEMANTIC_SHA256,
    "primary survivor-hash contract changed",
)
require(len(TERMINAL_RECORDS) == 19, "embedded terminal count changed")
require(
    semantic_sha256(TERMINAL_RECORDS) == EXPECTED_SURVIVOR_SEMANTIC_SHA256,
    "embedded terminal semantic ledger changed",
)


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


def units(modulus):
    return tuple(value for value in range(modulus) if gcd(value, modulus) == 1)


@lru_cache(maxsize=None)
def literal_family(M, n, width):
    """Every lifted unit AP, built by literal residue tests."""
    require(M % n == 0, ("trace period does not divide fibre", M, n))
    require(1 <= width <= n, ("trace width outside period", M, n, width))
    if width == n:
        return ((1 << M) - 1,)
    masks = set()
    for step in units(n):
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


def trace_law_and_containment_audit():
    """Exhaust the high/low law and local-AP containment through D=84."""
    needle_quotient_cases = 0
    literal_fibre_traces = 0
    for D in range(2, 85):
        for d in divisors(D):
            if d == 1:
                continue
            ell = (d + 6) // 7
            for t in divisors(D):
                M = D // t
                common = gcd(d, t)
                low, remainder = divmod(ell, common)
                n = d // common
                height = M // n
                width = (ell + common - 1) // common
                require(
                    width == (n + 6) // 7,
                    ("canonical local-width identity failed", D, d, t),
                )
                maximal_family = literal_family(M, n, width)
                for step in units(d):
                    for phase in range(d):
                        classes = {
                            (phase + index * step) % d
                            for index in range(ell)
                        }
                        loads = []
                        for residue in range(t):
                            trace = 0
                            for point in range(M):
                                if (residue + t * point) % d in classes:
                                    trace |= 1 << point
                            loads.append(trace.bit_count())
                            require(
                                any(
                                    trace & ~candidate == 0
                                    for candidate in maximal_family
                                ),
                                (
                                    "literal trace escaped maximal local AP",
                                    D,
                                    d,
                                    t,
                                    phase,
                                    step,
                                    residue,
                                ),
                            )
                            literal_fibre_traces += 1
                        allowed = {height * low}
                        if remainder:
                            allowed.add(height * (low + 1))
                        require(
                            set(loads).issubset(allowed),
                            ("high/low trace law failed", D, d, t),
                        )
                        if remainder:
                            require(
                                sum(
                                    load == height * (low + 1)
                                    for load in loads
                                )
                                == (t // common) * remainder,
                                ("high-fibre marginal failed", D, d, t),
                            )
                        needle_quotient_cases += 1
    require(
        needle_quotient_cases == 692416,
        ("needle/quotient control count changed", needle_quotient_cases),
    )
    require(
        literal_fibre_traces == 14279768,
        ("literal trace control count changed", literal_fibre_traces),
    )
    return needle_quotient_cases, literal_fibre_traces


def monotone_word(word):
    return all(
        word[pattern] <= word[pattern | (1 << index)]
        for pattern in range(8)
        for index in range(3)
    )


MONOTONE_WORDS = tuple(
    word
    for word in product((0, 1), repeat=8)
    if monotone_word(word)
)
require(len(MONOTONE_WORDS) == 20, "monotone three-bit word count changed")


def integer_transport_table(t, rewards):
    """Independent exact optimum over labelled integral status tables."""
    table = {(0, 0, 0): 0}
    for _position in range(t):
        updated = {}
        for marginal, value in table.items():
            for pattern in range(8):
                next_marginal = tuple(
                    marginal[index] + ((pattern >> index) & 1)
                    for index in range(3)
                )
                candidate = value + rewards[pattern]
                if candidate > updated.get(next_marginal, -1):
                    updated[next_marginal] = candidate
        table = updated
    return table


def transport_lp_direction_audit():
    """Compare the real LP with exact integral transport controls."""
    comparisons = 0
    equalities = 0
    fractional = 0
    maximum_gap = Q(0)
    tests = [(t, None) for t in range(1, 5)]
    for t in (5, 7, 9):
        values = (0, 1, t // 2, t - 1, t)
        selected = {
            (a, b, c)
            for a in values
            for b in values
            for c in values
            if (
                a in (0, t)
                or b in (0, t)
                or c in (0, t)
                or (3 * a + 5 * b + 7 * c) % 11 == 0
            )
        }
        tests.append((t, selected))

    for t, selected in tests:
        for rewards in MONOTONE_WORDS:
            exact = integer_transport_table(t, rewards)
            items = exact.items()
            if selected is not None:
                items = ((marginal, exact[marginal]) for marginal in selected)
            for marginals, integer_maximum in items:
                relaxed = primary.transport_optimum(t, marginals, rewards)
                require(
                    relaxed >= integer_maximum,
                    (
                        "real LP is not an upper relaxation",
                        t,
                        marginals,
                        rewards,
                        integer_maximum,
                        relaxed,
                    ),
                )
                comparisons += 1
                equalities += relaxed == integer_maximum
                fractional += relaxed.denominator != 1
                maximum_gap = max(maximum_gap, relaxed - integer_maximum)
    require(
        (comparisons, equalities, fractional, maximum_gap)
        == (10520, 10503, 17, Q(1, 2)),
        (
            "transport comparison ledger changed",
            comparisons,
            equalities,
            fractional,
            maximum_gap,
        ),
    )
    return comparisons, equalities, fractional, maximum_gap


def literal_needle(D, d, phase, step):
    width = (d + 6) // 7
    classes = {
        (phase + index * step) % d for index in range(width)
    }
    return frozenset(x for x in range(D) if x % d in classes)


def literal_fibre_loads(points, t):
    return tuple(
        sum(x % t == residue for x in points) for residue in range(t)
    )


def threshold_witness(points, D, ds, inclusive):
    for t in divisors(D):
        if t == 1:
            continue
        loads = literal_fibre_loads(points, t)
        _M, _inner_ds, _lows, marginals, capacities = primary.status_data(
            D, ds, t
        )
        for threshold in sorted(set(loads), reverse=True):
            if threshold <= 0:
                continue
            target_count = sum(load >= threshold for load in loads)
            if inclusive:
                rewards = tuple(
                    int(capacity >= threshold) for capacity in capacities
                )
            else:
                rewards = tuple(
                    int(capacity > threshold) for capacity in capacities
                )
            maximum = primary.transport_optimum(t, marginals, rewards)
            if Q(target_count) > maximum:
                return t, threshold, target_count, maximum, rewards
    return None


def transport_endpoint_controls():
    """Positive unions, a hostile target, and the inclusive boundary."""
    D = 42
    ds = (2, 3, 42)
    families = tuple(
        tuple(
            sorted(
                {
                    literal_needle(D, d, phase, step)
                    for phase in range(d)
                    for step in units(d)
                },
                key=lambda points: tuple(sorted(points)),
            )
        )
        for d in ds
    )
    positive_unions = tuple(
        sorted(
            {
                frozenset().union(*masks)
                for masks in product(*families)
            },
            key=lambda points: tuple(sorted(points)),
        )
    )
    strict_false_rejection = None
    for union in positive_unions:
        require(
            threshold_witness(union, D, ds, inclusive=True) is None,
            ("inclusive threshold rejected a positive union", union),
        )
        if strict_false_rejection is None:
            strict_false_rejection = threshold_witness(
                union, D, ds, inclusive=False
            )
    hostile = frozenset(range(12))
    require(
        not any(
            hostile <= frozenset().union(*masks)
            for masks in product(*families)
        ),
        "hostile threshold target has an exact cover",
    )
    hostile_witness = threshold_witness(hostile, D, ds, inclusive=True)
    require(
        hostile_witness is not None
        and hostile_witness[:4] == (6, 2, 6, Q(5)),
        ("hostile threshold witness changed", hostile_witness),
    )
    require(
        strict_false_rejection is not None,
        "strict threshold convention no longer false-rejects a positive union",
    )
    require(
        len(positive_unions) == 252,
        ("positive threshold union count changed", len(positive_unions)),
    )
    return len(positive_unions), hostile_witness, strict_false_rejection


@lru_cache(maxsize=None)
def incidence(M, n, width):
    """For each point, the bitset of literal family members through it."""
    masks = literal_family(M, n, width)
    rows = []
    for point in range(M):
        word = 0
        for index, mask in enumerate(masks):
            if (mask >> point) & 1:
                word |= 1 << index
        rows.append(word)
    return tuple(rows)


def first_point(mask):
    return (mask & -mask).bit_length() - 1


def indices(word):
    while word:
        bit = word & -word
        yield bit.bit_length() - 1
        word -= bit


@lru_cache(maxsize=None)
def one_coverable(target, key):
    if target == 0:
        return True
    masks = literal_family(*key)
    candidates = (1 << len(masks)) - 1
    remaining = target
    rows = incidence(*key)
    while remaining and candidates:
        point = first_point(remaining)
        candidates &= rows[point]
        remaining &= remaining - 1
    return bool(candidates)


@lru_cache(maxsize=None)
def two_coverable(target, left_key, right_key):
    if target == 0:
        return True
    point = first_point(target)
    left_masks = literal_family(*left_key)
    for index in indices(incidence(*left_key)[point]):
        if one_coverable(target & ~left_masks[index], right_key):
            return True
    right_masks = literal_family(*right_key)
    for index in indices(incidence(*right_key)[point]):
        if one_coverable(target & ~right_masks[index], left_key):
            return True
    return False


@lru_cache(maxsize=None)
def three_coverable(target, keys):
    if target == 0:
        return True
    point = first_point(target)
    for family_index, key in enumerate(keys):
        masks = literal_family(*key)
        other_keys = keys[:family_index] + keys[family_index + 1 :]
        for index in indices(incidence(*key)[point]):
            if two_coverable(target & ~masks[index], *other_keys):
                return True
    return False


def direct_slice(ranges, D, q, M, residue):
    """Literal support slice; every embedded terminal record has D=L."""
    require(q * M == D, ("bad fibre factorization", D, q, M))
    result = 0
    for height in range(M):
        point = residue + q * height
        if any(left <= point < right for left, right in ranges):
            result |= 1 << height
    return result


def terminal_pair_cover_audit():
    """Close the embedded 19 rows with the independent pair solver."""
    modulus_by_D = {
        38808: 99,
        77616: 99,
        129360: 98,
        388080: 99,
    }
    expected_D_counts = {
        38808: 2,
        77616: 6,
        129360: 1,
        388080: 10,
    }
    actual_D_counts = {
        D: sum(record[2] == D for record in TERMINAL_RECORDS)
        for D in modulus_by_D
    }
    require(
        actual_D_counts == expected_D_counts,
        ("terminal D profile changed", actual_D_counts),
    )
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
            ("terminal family count changed", key),
        )

    witnesses = []
    kill_moduli = {}
    for record in TERMINAL_RECORDS:
        F, L, D, support_count, d1, d2, d3, _margin = record
        require(D == L, ("terminal row left its full ruler", record))
        actual_L, ranges = support.safe_cell_ranges(F)
        require(actual_L == L, ("body ruler changed", F, L, actual_L))
        require(
            sum(right - left for left, right in ranges) == support_count,
            ("body support cardinality changed", record),
        )
        M = modulus_by_D[D]
        require(D % M == 0, ("terminal witness modulus escaped D", D, M))
        q = D // M
        residue = 0
        target = direct_slice(tuple(ranges), D, q, M, residue)
        keys = []
        for d in (d1, d2, d3):
            common = gcd(d, q)
            n = d // common
            ell = (d + 6) // 7
            width = (ell + common - 1) // common
            require(
                width == (n + 6) // 7,
                ("terminal canonical width identity failed", record, d),
            )
            require(
                n <= M and M % n == 0,
                ("terminal trace period escaped fibre", record, M, d, n),
            )
            keys.append((M, n, width))
        keys = tuple(sorted(keys))
        require(
            not three_coverable(target, keys),
            ("terminal target became pair-coverable", record, M),
        )
        witness = (
            M,
            q,
            residue,
            target.bit_count(),
            keys,
        )
        witnesses.append((record, witness))
        kill_moduli[M] = kill_moduli.get(M, 0) + 1

    require(
        kill_moduli == {99: 18, 98: 1},
        ("terminal kill-modulus ledger changed", kill_moduli),
    )
    witness_semantic = semantic_sha256(witnesses)
    return (
        len(witnesses),
        kill_moduli,
        expected_family_counts,
        witness_semantic,
        witnesses[0],
        witnesses[-1],
    )


def local_solver_controls():
    keys = ((12, 3, 1), (12, 4, 1), (12, 6, 1))
    positive = (
        literal_family(*keys[0])[0]
        | literal_family(*keys[1])[0]
        | literal_family(*keys[2])[0]
    )
    require(three_coverable(positive, keys), "positive local union rejected")
    hostile = (1 << 12) - 1
    hostile_keys = ((12, 12, 1),) * 3
    require(
        not three_coverable(hostile, hostile_keys),
        "hostile local target became coverable",
    )
    return positive.bit_count(), hostile.bit_count()


def main():
    (
        needle_quotient_cases,
        literal_fibre_traces,
    ) = trace_law_and_containment_audit()
    (
        lp_comparisons,
        lp_equalities,
        lp_fractional,
        lp_maximum_gap,
    ) = transport_lp_direction_audit()
    (
        positive_threshold_unions,
        hostile_threshold_witness,
        strict_false_rejection,
    ) = transport_endpoint_controls()
    positive_local_size, hostile_local_size = local_solver_controls()
    (
        terminal_kills,
        terminal_kill_moduli,
        terminal_family_counts,
        terminal_witness_semantic,
        first_terminal_witness,
        last_terminal_witness,
    ) = terminal_pair_cover_audit()

    print("LRC14 three-drift terminal independent audit candidate")
    print(f"primary_script_sha256={file_sha256(PRIMARY_PATH)}")
    print(
        "terminal_record_semantic_sha256="
        f"{semantic_sha256(TERMINAL_RECORDS)}"
    )
    print(f"needle_quotient_cases={needle_quotient_cases}")
    print(f"literal_fibre_traces={literal_fibre_traces}")
    print(f"monotone_status_words={len(MONOTONE_WORDS)}")
    print(f"lp_integer_comparisons={lp_comparisons}")
    print(f"lp_equalities={lp_equalities}")
    print(f"lp_fractional_cases={lp_fractional}")
    print(f"lp_maximum_gap={lp_maximum_gap}")
    print(f"positive_threshold_unions={positive_threshold_unions}")
    print(f"hostile_threshold_witness={hostile_threshold_witness}")
    print(f"strict_false_rejection={strict_false_rejection}")
    print(f"positive_local_control_size={positive_local_size}")
    print(f"hostile_local_control_size={hostile_local_size}")
    print(f"terminal_kills={terminal_kills}")
    print("terminal_survivors=0")
    print(f"terminal_kill_moduli={terminal_kill_moduli}")
    print(f"terminal_family_counts={terminal_family_counts}")
    print(
        "terminal_witness_semantic_sha256="
        f"{terminal_witness_semantic}"
    )
    print(f"first_terminal_witness={first_terminal_witness}")
    print(f"last_terminal_witness={last_terminal_witness}")
    print(f"literal_family_cache={literal_family.cache_info()}")
    print(f"one_coverable_cache={one_coverable.cache_info()}")
    print(f"two_coverable_cache={two_coverable.cache_info()}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()

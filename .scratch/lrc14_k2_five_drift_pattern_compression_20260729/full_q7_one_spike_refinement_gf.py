#!/usr/bin/env python3
"""Global exact c=4 one-spike refinement of the k=2 q=D/7 screen.

The expected-spike screen remembers only

    c = #{denominators d not dividing q},  q=D/7.

When c=4 there is exactly one spike denominator d|q.  Write d=7a+r.  At
common phase u its literal enlarged mask hits exactly

    Y_d(u) = (q/d) * (a + X_d(u))

q-fibres, where X_d is {0,1}-valued and mu(X_d=1)=r/7.  Four uniform masks
pay four points in every q-fibre, so coverage forces Y_d >= N_4 on the
compact aligned-safe carrier, where

    N_4 = #{b mod q : lambda_q(b)>4}.

If the threshold needs the extra tooth, the containing event is proper open;
therefore its mass r/7 must be strictly greater than u_2=66/91.  The script
retains d in the exact-lcm Mobius generating function and applies this exact
three-case law over the complete 27,163-row support-transfer ledger.
"""

from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE_PATH = (
    ROOT
    / ".scratch"
    / "lrc14_k2_five_drift_pattern_compression_20260729"
    / "full_q7_expected_spike_gf.py"
)
EXPECTED_BASE_SHA256 = (
    "3e26985308834c4c0dcece6f17a214fa7622c981a9a0e9f085930512347f39b7"
)
EXPECTED_BASE_OUTPUT_SHA256 = (
    "7ce4080f899f6e319b24977fc3b4f13a30290452f2b8113bb8ddb1f60e136d30"
)
BASE_OUTPUT_PATH = BASE_PATH.with_suffix(".ordinary.out")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(file_sha256(BASE_PATH) == EXPECTED_BASE_SHA256, "frozen base script changed")
require(
    file_sha256(BASE_OUTPUT_PATH) == EXPECTED_BASE_OUTPUT_SHA256,
    "frozen base ordinary output changed",
)
spec = spec_from_file_location("lrc14_k2_full_q7_expected_spike", BASE_PATH)
base = module_from_spec(spec)
spec.loader.exec_module(base)
support = base.support
combined = base.combined


TWO_SAFE_FLOOR = Q(66, 91)
EXPECTED_INCOMING_C4_SHAPES = 3246262178
EXPECTED_INCOMING_C4_OCCURRENCES = 41980227729
EXPECTED_FULL_INCOMING_SHAPES = 36962285549
EXPECTED_FULL_INCOMING_OCCURRENCES = 320011786356
EXPECTED_RESULT = {
    "controls": (14, 219, 127, 12208, 283),
    "incoming": (3246262178, 41980227729, 3609, 2940, 134),
    "surviving": (3208397602, 39762283189, 3373, 2871, 133),
    "killed_occurrences": 2217944540,
    "outcomes": {"full": 39336054399, "extra_r6": 426228790, "killed": 2217944540},
    "removed_divisors": (588,),
    "refined_full": (36924420973, 317793841816),
    "semantic_sha256": (
        "14c9e777f98ca3d44ec3747eeb08153f4edffd23e0d74ff8f4abf4f5195440fa"
    ),
    "minimum_kill": (
        168,
        (1, 2, 3, 4, 6, 12),
        168,
        88,
        24,
        4,
        2,
        0,
        12,
        "2/7",
        184,
    ),
    "minimum_survivor": (
        168,
        (1, 2, 3, 4, 6, 12),
        168,
        88,
        24,
        4,
        6,
        0,
        4,
        "6/7",
        204,
    ),
}


def one_spike_distribution(D):
    """Exact lcm-D five-multiset count by the lone spike denominator d|q.

    For E|D and fixed d|gcd(E,q), choose d once and choose a four-multiset
    from the U_D(E) divisor types that do not divide q.  Mobius inversion in
    E makes the total lcm exactly D.
    """
    require(D % 7 == 0, "one-spike distribution requested off septimal locus")
    q = D // 7
    result = Counter()
    divisors_D = support.divisors(D)
    divisors_q = tuple(d for d in support.divisors(q) if d > 1)
    for E in divisors_D:
        sign = base.mobius(D // E)
        if not sign:
            continue
        uniform_types = sum(
            divisor > 1 and q % divisor != 0
            for divisor in support.divisors(E)
        )
        uniform_multisets = base.multichoose(uniform_types, 4)
        if not uniform_multisets:
            continue
        for d in divisors_q:
            if E % d == 0:
                result[d] += sign * uniform_multisets
    result += Counter()
    require(
        all(value >= 0 for value in result.values()),
        ("negative exact-lcm one-spike coefficient", D, result),
    )
    require(
        sum(result.values()) == base.uniform_count_distribution(D)[4],
        ("one-spike distribution does not recover c=4", D),
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
            spikes = tuple(d for d in values if q % d == 0)
            if len(spikes) == 1:
                brute[spikes[0]] += 1
        require(
            brute == one_spike_distribution(D),
            ("literal one-spike distribution control failed", D, brute),
        )
        cases += 1
    return cases


def recurrence_distribution(D):
    """Independent divisor-downset recurrence for every weighted d count."""
    require(D % 7 == 0, "one-spike recurrence requested off septimal locus")
    q = D // 7
    divisors_D = tuple(sorted(support.divisors(D)))
    divisors_q = tuple(d for d in support.divisors(q) if d > 1)
    exact = {}
    for E in divisors_D:
        uniform_types = sum(
            divisor > 1 and q % divisor != 0
            for divisor in support.divisors(E)
        )
        at_most_count = base.multichoose(uniform_types, 4)
        current = Counter(
            {d: at_most_count for d in divisors_q if E % d == 0 and at_most_count}
        )
        for F in support.divisors(E):
            if F == E:
                continue
            for d, value in exact[F].items():
                current[d] -= value
        require(
            all(value >= 0 for value in current.values()),
            ("negative downward-recurrence coefficient", D, E, current),
        )
        current += Counter()
        exact[E] = current
    return exact[D]


def all_D_recurrence_controls(divisors):
    cases = 0
    for D in sorted(divisors):
        require(
            recurrence_distribution(D) == one_spike_distribution(D),
            ("all-D weighted recurrence control failed", D),
        )
        cases += 1
    return cases


def class_count(d, numerator, u):
    require(gcd(d, numerator) == 1, "two-level control numerator is not a unit")
    return sum(
        base.norm_fraction(Q(numerator) * (residue + u) / d) < Q(1, 14)
        for residue in range(d)
    )


def two_level_law_controls():
    """Literal cell audit of #active classes=a+X and mu(X)=r/7."""
    cases = 0
    probes = 0
    for d in range(2, 21):
        low, remainder = divmod(d, 7)
        high = low + bool(remainder)
        for numerator in range(1, d + 1):
            if gcd(d, numerator) != 1:
                continue
            # Every boundary lies on the 1/(14*numerator) grid.  Midpoints
            # therefore sample all equal-measure open cells exactly.
            cell_count = 14 * numerator
            counts = Counter(
                class_count(
                    d,
                    numerator,
                    Q(2 * index + 1, 2 * cell_count),
                )
                for index in range(cell_count)
            )
            require(
                set(counts) <= {low, high},
                ("literal hit count is not two-level", d, numerator, counts),
            )
            high_cells = counts[high] if high != low else 0
            require(
                Q(high_cells, cell_count) == Q(remainder, 7),
                ("extra-tooth mass changed", d, numerator, counts),
            )
            cases += 1
            probes += cell_count
    return cases, probes


def exact_event_mass(d, q, threshold):
    """Mass of {u:Y_d(u)>=threshold}; None denotes the full circle."""
    require(q % d == 0, "spike denominator does not divide q")
    low, remainder = divmod(d, 7)
    scale = q // d
    low_load = scale * low
    high_load = scale * (low + bool(remainder))
    if threshold <= low_load:
        return None, low_load, high_load
    if threshold <= high_load and remainder:
        return Q(remainder, 7), low_load, high_load
    return Q(0), low_load, high_load


def exact_one_spike_passes(d, q, threshold):
    event_mass, _low, _high = exact_event_mass(d, q, threshold)
    return event_mass is None or event_mass > TWO_SAFE_FLOOR


def event_rule_controls():
    checks = 0
    for d in range(2, 101):
        q = lcm(d, 7)
        low, remainder = divmod(d, 7)
        scale = q // d
        low_load = scale * low
        high_load = scale * (low + bool(remainder))
        mass, _, _ = exact_event_mass(d, q, low_load)
        require(mass is None, ("base-covered event is not full", d))
        checks += 1
        if high_load > low_load:
            mass, _, _ = exact_event_mass(d, q, high_load)
            require(mass == Q(remainder, 7), ("extra event mass changed", d))
            require(
                exact_one_spike_passes(d, q, high_load) == (remainder == 6),
                ("strict aligned-floor residue rule changed", d),
            )
            checks += 1
        mass, _, _ = exact_event_mass(d, q, high_load + 1)
        require(mass == 0, ("impossible event is nonempty", d))
        checks += 1
    return checks


def main():
    distribution_controls = brute_distribution_controls()
    law_control_cases, law_control_probes = two_level_law_controls()
    rule_controls = event_rule_controls()

    rows_by_D = defaultdict(list)
    for body in combinations(range(1, 15), 6):
        L, ranges = support.safe_cell_ranges(body)
        for D in support.divisors(L):
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > base.SUPPORT_CUTOFF:
                continue
            require(D % 7 == 0, ("k2 support row is not septimal", body, D))
            q = D // 7
            arcs = combined.projected_support_arcs(D, ranges)
            histogram = combined.residue_load_histogram(arcs, q)
            N4 = sum(count for load, count in histogram if load > 4)
            rows_by_D[D].append((body, L, support_count, histogram, N4))

    require(sum(map(len, rows_by_D.values())) == base.EXPECTED_ROWS, "row count changed")
    require(len(rows_by_D) == base.EXPECTED_DIVISORS, "divisor count changed")
    recurrence_controls = all_D_recurrence_controls(rows_by_D)

    incoming_shapes = 0
    incoming_occurrences = 0
    incoming_rows = set()
    incoming_bodies = set()
    incoming_divisors = set()
    surviving_shapes = 0
    surviving_occurrences = 0
    surviving_rows = set()
    surviving_bodies = set()
    surviving_divisors = set()
    incoming_by_residue = Counter()
    surviving_by_residue = Counter()
    killed_by_residue = Counter()
    incoming_by_d = Counter()
    surviving_by_d = Counter()
    killed_by_d = Counter()
    event_types = Counter()
    outcomes = Counter()
    semantic = sha256()
    minimum_survivor = None
    minimum_kill = None

    for D in sorted(rows_by_D):
        q = D // 7
        distribution = one_spike_distribution(D)
        for d, shape_count in sorted(distribution.items()):
            incoming_shape = False
            surviving_shape = False
            for body, L, support_count, histogram, N4 in rows_by_D[D]:
                # Reconstruct precisely the incoming c=4 Markov ledger.
                expected_spike_passes = N4 == 0 or 66 * N4 < 13 * q
                if not expected_spike_passes:
                    continue
                incoming_shape = True
                incoming_occurrences += shape_count
                incoming_rows.add((body, D))
                incoming_bodies.add(body)
                incoming_divisors.add(D)
                incoming_by_residue[d % 7] += shape_count
                incoming_by_d[d] += shape_count

                event_mass, low_load, high_load = exact_event_mass(d, q, N4)
                event_key = (
                    "full" if event_mass is None else event_mass,
                    d % 7,
                    N4,
                    low_load,
                    high_load,
                )
                event_types[event_key] += shape_count
                record = (
                    D,
                    body,
                    L,
                    support_count,
                    q,
                    N4,
                    d,
                    low_load,
                    high_load,
                    "full" if event_mass is None else str(event_mass),
                    shape_count,
                )
                if event_mass is None or event_mass > TWO_SAFE_FLOOR:
                    surviving_shape = True
                    surviving_occurrences += shape_count
                    surviving_rows.add((body, D))
                    surviving_bodies.add(body)
                    surviving_divisors.add(D)
                    surviving_by_residue[d % 7] += shape_count
                    surviving_by_d[d] += shape_count
                    outcomes["full" if event_mass is None else "extra_r6"] += shape_count
                    semantic.update(f"{record}\n".encode())
                    if minimum_survivor is None or record < minimum_survivor:
                        minimum_survivor = record
                else:
                    outcomes["killed"] += shape_count
                    killed_by_residue[d % 7] += shape_count
                    killed_by_d[d] += shape_count
                    if minimum_kill is None or record < minimum_kill:
                        minimum_kill = record
            if incoming_shape:
                incoming_shapes += shape_count
            if surviving_shape:
                surviving_shapes += shape_count

    require(
        (incoming_shapes, incoming_occurrences)
        == (EXPECTED_INCOMING_C4_SHAPES, EXPECTED_INCOMING_C4_OCCURRENCES),
        "incoming c=4 ledger changed",
    )
    require(
        surviving_occurrences + sum(killed_by_residue.values())
        == incoming_occurrences,
        "one-spike refinement lost incoming occurrences",
    )
    require(
        all(
            residue == 6
            for (mass, residue, _N4, _low, _high), count in event_types.items()
            if mass != "full" and mass > TWO_SAFE_FLOOR and count
        ),
        "a non-full surviving event does not have residue six",
    )

    refined_full_shapes = (
        EXPECTED_FULL_INCOMING_SHAPES - EXPECTED_INCOMING_C4_SHAPES + surviving_shapes
    )
    refined_full_occurrences = (
        EXPECTED_FULL_INCOMING_OCCURRENCES
        - EXPECTED_INCOMING_C4_OCCURRENCES
        + surviving_occurrences
    )

    require(
        (
            distribution_controls,
            recurrence_controls,
            law_control_cases,
            law_control_probes,
            rule_controls,
        )
        == EXPECTED_RESULT["controls"],
        "control census changed",
    )
    require(
        (
            incoming_shapes,
            incoming_occurrences,
            len(incoming_rows),
            len(incoming_bodies),
            len(incoming_divisors),
        )
        == EXPECTED_RESULT["incoming"],
        "incoming detailed census changed",
    )
    require(
        (
            surviving_shapes,
            surviving_occurrences,
            len(surviving_rows),
            len(surviving_bodies),
            len(surviving_divisors),
        )
        == EXPECTED_RESULT["surviving"],
        "surviving detailed census changed",
    )
    require(
        incoming_occurrences - surviving_occurrences
        == EXPECTED_RESULT["killed_occurrences"],
        "killed occurrence census changed",
    )
    require(
        (refined_full_shapes, refined_full_occurrences)
        == EXPECTED_RESULT["refined_full"],
        "refined full census changed",
    )
    require(
        semantic.hexdigest() == EXPECTED_RESULT["semantic_sha256"],
        "survivor semantic digest changed",
    )
    require(minimum_kill == EXPECTED_RESULT["minimum_kill"], "minimum kill changed")
    require(
        minimum_survivor == EXPECTED_RESULT["minimum_survivor"],
        "minimum survivor changed",
    )
    require(dict(outcomes) == EXPECTED_RESULT["outcomes"], "outcome split changed")
    require(
        tuple(sorted(incoming_divisors - surviving_divisors))
        == EXPECTED_RESULT["removed_divisors"],
        "removed divisor set changed",
    )

    print("LRC14 k=2/p=5 global c=4 exact one-spike refinement")
    print(f"base_script_sha256={file_sha256(BASE_PATH)}")
    print(f"base_output_sha256={file_sha256(BASE_OUTPUT_PATH)}")
    print(f"weighted_distribution_control_cases={distribution_controls}")
    print(f"all_D_weighted_recurrence_control_cases={recurrence_controls}")
    print(f"two_level_law_control_cases={law_control_cases}")
    print(f"two_level_law_control_probes={law_control_probes}")
    print(f"event_rule_control_cases={rule_controls}")
    print(f"support_rows={sum(map(len, rows_by_D.values()))}")
    print(f"support_divisors={len(rows_by_D)}")
    print(f"incoming_c4_shapes={incoming_shapes}")
    print(f"incoming_c4_occurrences={incoming_occurrences}")
    print(f"incoming_c4_rows={len(incoming_rows)}")
    print(f"incoming_c4_bodies={len(incoming_bodies)}")
    print(f"incoming_c4_divisors={len(incoming_divisors)}")
    print(f"surviving_c4_shapes={surviving_shapes}")
    print(f"surviving_c4_occurrences={surviving_occurrences}")
    print(f"surviving_c4_rows={len(surviving_rows)}")
    print(f"surviving_c4_bodies={len(surviving_bodies)}")
    print(f"surviving_c4_divisors={len(surviving_divisors)}")
    print(f"killed_c4_occurrences={incoming_occurrences - surviving_occurrences}")
    print(f"incoming_by_residue={incoming_by_residue}")
    print(f"surviving_by_residue={surviving_by_residue}")
    print(f"killed_by_residue={killed_by_residue}")
    print(f"incoming_distinct_spike_denominators={len(incoming_by_d)}")
    print(f"surviving_distinct_spike_denominators={len(surviving_by_d)}")
    print(f"surviving_by_d_top30={surviving_by_d.most_common(30)}")
    print(f"killed_by_d_top30={killed_by_d.most_common(30)}")
    print(f"outcome_occurrences={outcomes}")
    print(f"removed_divisors={sorted(incoming_divisors - surviving_divisors)}")
    print(f"event_type_count={len(event_types)}")
    print(f"event_types_top30={event_types.most_common(30)}")
    print(f"refined_full_shapes={refined_full_shapes}")
    print(f"refined_full_occurrences={refined_full_occurrences}")
    print(f"semantic_sha256={semantic.hexdigest()}")
    print(f"minimum_kill={minimum_kill}")
    print(f"minimum_survivor={minimum_survivor}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()

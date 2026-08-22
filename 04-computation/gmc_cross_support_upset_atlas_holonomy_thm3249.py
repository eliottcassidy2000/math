#!/usr/bin/env python3
"""Exact cross-support atlas and constant-gauge obstruction for THM-3249."""

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from functools import reduce
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import gcd, lcm
from multiprocessing import Pool, freeze_support
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
THM3238_SCRIPT = ROOT / (
    "04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py"
)
THM3238_OUTPUT = ROOT / (
    "05-knowledge/results/gmc_complete_physical_bank_unique_reset_thm3238.out"
)
THM3244_SCRIPT = ROOT / (
    "04-computation/gmc_unique_reset_rips_nonmorse_thm3244.py"
)
THM3244_OUTPUT = ROOT / (
    "05-knowledge/results/gmc_unique_reset_rips_nonmorse_thm3244.out"
)
PINS = {
    THM3238_SCRIPT:
        "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa",
    THM3238_OUTPUT:
        "77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b",
    THM3244_SCRIPT:
        "3ff0babc41e35e6a185b0ff442cfb9284d9688360c0b96cd947c1128e16400ba",
    THM3244_OUTPUT:
        "27dcd7c68e628465a1f09a564be0be366ded6075ef009a3d85d029f8f18605c9",
}


def lf_hash(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for pinned_path, expected_hash in PINS.items():
    require(lf_hash(pinned_path) == expected_hash,
            "dependency pin drift: " + pinned_path.name)

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)),
        "optimization-sensitive assert found")
require(not any(isinstance(node, ast.Constant) and isinstance(node.value, float)
                for node in ast.walk(syntax)), "floating literal found")

spec = spec_from_file_location("thm3238_companion", THM3238_SCRIPT)
T = module_from_spec(spec)
spec.loader.exec_module(T)
BANK = T.UP["BANKS"][1]


SMALL_POLES = (5, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1)
SMALL_RESET = (1, 2, 2, 3, 4, 5)
FULL_POLES = T.POLES
FULL_RESET = tuple(sorted(T.RESET))
EXPECTED_DEPTHS = (5, 13, 24, 35, 42, 42, 35, 24, 13, 5, 1)


def physical_states(poles):
    multiplicities = Counter(poles)
    roots = tuple(sorted(multiplicities))
    states = tuple(sorted((
        tuple(root for root, count in zip(roots, choices)
              for _ in range(count))
        for choices in product(
            *(range(multiplicities[root] + 1) for root in roots)
        )
        if any(choices)
    ), key=lambda state: (len(state), state)))
    return states


SMALL_STATES = physical_states(SMALL_POLES)
require(len(SMALL_STATES) == 239, "small physical-bank census drift")
require(tuple(sum(len(state) == depth for state in SMALL_STATES)
              for depth in range(1, len(SMALL_POLES) + 1)) == EXPECTED_DEPTHS,
        "small depth census drift")
require(Counter(SMALL_RESET) <= Counter(SMALL_POLES),
        "small reset left physical bank")


# These are the states in the exact separating Farkas mixture.  The first
# three belong to support (1,2), the remaining sixteen to support (1,3).
WITNESS = (
    (1, 2, (2, 3), Fraction(110021555, 210299496)),
    (1, 2, (2, 2), Fraction(244638624, 855590729)),
    (1, 2, (2, 2, 2), Fraction(13822588, 86007735)),
    (1, 3, (1,), Fraction(981, 345013099)),
    (1, 3, (2, 2), Fraction(740, 175674781)),
    (1, 3, (1, 1, 1, 1), Fraction(1988, 947209091)),
    (1, 3, (1, 1, 1, 2), Fraction(94, 29209079)),
    (1, 3, (2, 2, 2, 3, 3), Fraction(2741, 436074315)),
    (1, 3, (1, 1, 1, 1, 2, 2, 2), Fraction(1886, 428306147)),
    (1, 3, (4, 4, 5, 5, 6, 7, 8), Fraction(4040735, 402653491)),
    (1, 3, (3, 3, 4, 4, 5, 5, 7, 8), Fraction(1681, 222909716)),
    (1, 3, (1, 1, 1, 1, 5, 5, 6, 7, 8),
     Fraction(495494, 846338163)),
    (1, 3, (1, 1, 1, 1, 2, 2, 2, 3, 3, 4),
     Fraction(3301, 978297006)),
    (1, 3, (1, 1, 1, 1, 4, 4, 5, 6, 7, 8),
     Fraction(498441, 407384384)),
    (1, 3, (2, 2, 3, 3, 4, 4, 5, 5, 6, 7),
     Fraction(1, 223456666)),
    (1, 3, (2, 3, 3, 4, 4, 5, 5, 6, 7, 8),
     Fraction(1338659, 314691983)),
    (1, 3, (2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7),
     Fraction(17, 863033027)),
    (1, 3, (1, 1, 1, 1, 2, 2, 2, 5, 5, 6, 7, 8),
     Fraction(621161, 731987841)),
    (1, 3, (1, 1, 1, 1, 2, 2, 2, 4, 5, 5, 6, 7, 8),
     Fraction(11562097, 875223752)),
)
require(all(weight > 0 for _, _, _, weight in WITNESS),
        "Farkas weight lost positivity")


def row_coordinates(vector):
    return tuple(
        sum(vector[entry[0]][shape] for shape in upset)
        for entry, upset in zip(T.CERTIFICATE, T.UPSETS)
    )


def compute_task(task):
    a_value, b_value, state = task
    vector = T.coefficient_vectors(1, BANK, a_value, b_value, state)
    return task, row_coordinates(vector)


ACTIVE_ROWS = (5, 8, 10, 15, 16, 21)
ACTIVE_FACTORS = (
    Fraction(213915797, 685946076),
    Fraction(81328689, 164203297),
    Fraction(6210895, 131921019),
    Fraction(255973, 944937251),
    Fraction(36321374, 498402635),
    Fraction(55150150, 759378773),
)
ACTIVE_NORMALIZERS = (
    1538222247466205184,
    10668668469065531392,
    28568269398144,
    1484288,
    283951084704,
    386978585480040,
)
EXPECTED_PRIMITIVE = (
    31540805104781016214807762327098361741908455758236214857697597793166077169575373070006205,
    7222551733565666973704472289597930205521464117990467362419953714461056410853911375394360,
    256386550655705646640018746246846057992771457973407640357382724455640369346859861323130803200,
    28393046540266691701145009302241688233214029700961546492881520901040914942725957247847268249416960,
    39927984611766798370362293023501771795759571792906966553729869583513560353170793810845831929856,
    29197138787760557148002850343923059397376740901033387386489565737312725443107887357993779200,
)
EXPECTED_STATE_DIGEST = (
    "7a2dda702d6fd61f300771600a6cae098a81ec6bf74edb3e24def667f8dfd15e"
)
EXPECTED_PRIMITIVE_DIGEST = (
    "dccdc51d7b695304b4abc937b9957ef79ef1748206596f116a3183f2a8e0ad8d"
)
EXPECTED_WITNESS_DIGEST = (
    "89c1b11710a1c1c0793dfa0b0e8b385b9fb1fb8a265d324b4bde463d436a38ae"
)
EXPECTED_EXPECTATION_DIGEST = (
    "1d26244164b903e6d506d61962a4c49b580a73f0ce99a3d241fd884309616286"
)
EXPECTED_SMALL_COVER_DIGEST = (
    "4040c5bbb997626a0d63dca1b9e0b431d87191802ca29874812b2af126cb50d0"
)
EXPECTED_SMALL_PAIR_DIGEST = (
    "4e6aac29e5bd52a789a18f9b53ae6b8538a357f3754d13fd70294c53cda92a8f"
)
EXPECTED_SMALL_COVERING_PAIRS = (
    (2, 3), (2, 4), (2, 6), (2, 10), (2, 13), (2, 19),
    (3, 9), (3, 11), (3, 14), (3, 16), (3, 18), (3, 22),
    (4, 9), (4, 11), (4, 14), (4, 16), (4, 18), (5, 8),
    (6, 9), (6, 11), (6, 14), (6, 16), (6, 18),
    (7, 16), (7, 18), (7, 22),
    (9, 10), (9, 13), (9, 17), (9, 19),
    (10, 11), (10, 14), (10, 16), (10, 18), (10, 22),
    (11, 13), (11, 17), (11, 19),
    (13, 14), (13, 16), (13, 18), (13, 22), (14, 19),
    (16, 17), (16, 19), (16, 21),
    (17, 18), (17, 22), (18, 19), (18, 21),
    (19, 22), (21, 22),
)
EXPECTED_TARGET_COVERING_PAIRS = (
    (2, 7), (2, 10), (2, 13), (2, 17), (2, 19), (2, 21),
    (3, 9), (3, 11), (3, 16), (3, 22),
    (7, 14), (7, 18), (7, 22),
    (10, 11), (10, 16), (10, 22),
    (11, 13), (11, 17), (11, 21),
    (12, 13), (12, 19),
    (13, 14), (13, 18), (13, 22),
    (14, 19),
    (16, 17), (16, 21),
    (17, 22), (18, 19), (19, 22), (21, 22),
)
EXPECTED_COMMON_COVERING_PAIRS = (
    (2, 10), (2, 13), (2, 19),
    (3, 9), (3, 11), (3, 16), (3, 22),
    (7, 18), (7, 22),
    (10, 11), (10, 16), (10, 22),
    (11, 13), (11, 17),
    (13, 14), (13, 18), (13, 22),
    (14, 19),
    (16, 17), (16, 21),
    (17, 22), (18, 19), (19, 22), (21, 22),
)
EXPECTED_COMMON_PAIR_DIGEST = (
    "e622d25bfc96d28269eb9ecb86ec43211754530ba819361d880989321640dd63"
)
SMALL_BLEND_TRAP_A = (1, 1, 2, 2, 3, 4, 5)
SMALL_BLEND_TRAP_B = (1, 2, 2)


def clear_positive_weights(weights):
    common = reduce(lcm, (weight.denominator for weight in weights), 1)
    integers = tuple(weight.numerator * (common // weight.denominator)
                     for weight in weights)
    divisor = reduce(gcd, integers)
    return tuple(value // divisor for value in integers)


def digest(values):
    return hashlib.sha256("\n".join(map(str, values)).encode("ascii")).hexdigest()


def sign_census(values):
    return (sum(value < 0 for value in values),
            sum(value == 0 for value in values),
            sum(value > 0 for value in values))


def toward_reset_neighbours(state, poles, reset):
    """Return one-pole moves decreasing multiset l1 distance to reset."""

    current = Counter(state)
    pole_counts = Counter(poles)
    reset_counts = Counter(reset)
    answer = []
    for root in sorted(set(pole_counts) | set(reset_counts)):
        changed = current.copy()
        if current[root] > reset_counts[root]:
            changed[root] -= 1
            if changed[root] == 0:
                del changed[root]
        elif current[root] < reset_counts[root]:
            changed[root] += 1
        else:
            continue
        candidate = tuple(sorted(changed.elements()))
        if candidate:
            answer.append(candidate)
    return tuple(answer)


def positive_blend_trap_interval(state, state_index, rows, poles, reset):
    """Closed local-trap interval for lambda*row_2 + row_10, lambda>0."""

    lower = Fraction(0)
    upper = None
    for target in toward_reset_neighbours(state, poles, reset):
        delta_two = rows[state_index[target]][1] - rows[state_index[state]][1]
        delta_ten = rows[state_index[target]][9] - rows[state_index[state]][9]
        if delta_two > 0:
            bound = Fraction(-delta_ten, delta_two)
            require(bound > 0, "empty small blend trap interval")
            upper = bound if upper is None else min(upper, bound)
        elif delta_two < 0:
            lower = max(lower, Fraction(-delta_ten, delta_two))
        else:
            require(delta_ten <= 0, "empty small blend trap interval")
    require(upper is None or lower <= upper, "inconsistent small blend interval")
    return lower, upper


def main():
    require(T.reduced_poles(1, BANK, 1, 2)[0] == SMALL_POLES,
            "small pole bank drift")
    require(tuple(sorted(T.residual_roots(
        1, T.dominant_row(1), 1, 2))) == SMALL_RESET,
        "small reset drift")

    tasks = tuple((1, 2, state) for state in SMALL_STATES)
    tasks += tuple((a_value, b_value, state)
                   for a_value, b_value, state, _ in WITNESS
                   if (a_value, b_value) == (1, 3))
    with Pool(processes=4) as pool:
        rows = dict(pool.imap_unordered(compute_task, tasks, chunksize=1))
    require(len(rows) == len(set(tasks)), "response task reconstruction lost states")

    small_rows = tuple(rows[(1, 2, state)] for state in SMALL_STATES)
    observed_normalizers = tuple(
        max(abs(row[index - 1]) for row in small_rows)
        for index in ACTIVE_ROWS
    )
    require(observed_normalizers == ACTIVE_NORMALIZERS,
            "active response normalizer drift")
    rational_weights = tuple(
        factor / normalizer
        for factor, normalizer in zip(ACTIVE_FACTORS, ACTIVE_NORMALIZERS)
    )
    primitive = clear_positive_weights(rational_weights)
    require(primitive == EXPECTED_PRIMITIVE, "small primitive dual drift")

    active_weights = [0] * len(T.CERTIFICATE)
    for index, value in zip(ACTIVE_ROWS, primitive):
        active_weights[index - 1] = value
    active_coordinates = tuple(
        sum(weight * coordinate for weight, coordinate in zip(active_weights, row))
        for row in small_rows
    )
    require(sign_census(active_coordinates) == (238, 1, 0),
            "small active certificate sign drift")
    require(tuple(state for state, value in zip(SMALL_STATES, active_coordinates)
                  if value == 0) == (SMALL_RESET,),
            "small active certificate lost unique reset")

    # Setting every inactive chart coefficient to one produces an explicitly
    # strictly positive twenty-two-chart section without disturbing the signs.
    positive_weights = tuple(
        active_weights[index] if active_weights[index] else 1
        for index in range(len(active_weights))
    )
    positive_coordinates = tuple(
        sum(weight * coordinate for weight, coordinate in zip(positive_weights, row))
        for row in small_rows
    )
    require(sign_census(positive_coordinates) == (238, 1, 0),
            "all-positive local section sign drift")
    require(tuple(state for state, value in zip(SMALL_STATES, positive_coordinates)
                  if value == 0) == (SMALL_RESET,),
            "all-positive section lost unique reset")

    raw_weights = T.clear_multipliers(T.CERTIFICATE)
    raw_coordinates = tuple(
        sum(weight * coordinate for weight, coordinate in zip(raw_weights, row))
        for row in small_rows
    )
    require(sign_census(raw_coordinates) == (115, 1, 123),
            "transported THM-3238 raw-weight census drift")
    require(tuple(state for state, value in zip(SMALL_STATES, raw_coordinates)
                  if value == 0) == (SMALL_RESET,),
            "raw transported weights gained an extra zero")

    # Local atlas test: a row covers a state if it rises along at least one
    # one-pole move that reduces multiset l1 distance to Q_12.  Any covering
    # pair therefore supplies a monotone chart-by-chart route to the reset.
    state_index = {state: index for index, state in enumerate(SMALL_STATES)}
    nonreset_states = frozenset(SMALL_STATES) - {SMALL_RESET}
    covers = []
    for chart in range(len(T.CERTIFICATE)):
        covered = set()
        for state in nonreset_states:
            source = small_rows[state_index[state]][chart]
            if any(small_rows[state_index[target]][chart] > source
                   for target in toward_reset_neighbours(
                       state, SMALL_POLES, SMALL_RESET)):
                covered.add(state)
        covers.append(frozenset(covered))
    covering_pairs = tuple(
        (left + 1, right + 1)
        for left in range(len(covers))
        for right in range(left + 1, len(covers))
        if covers[left] | covers[right] == nonreset_states
    )
    require(covering_pairs == EXPECTED_SMALL_COVERING_PAIRS,
            "small one-pole covering-pair census drift")
    single_cover_counts = tuple(map(len, covers))
    require(max(single_cover_counts) == 230
            and single_cover_counts[15] == 230,
            "small best single-chart coverage drift")
    require((len(covers[1]), len(covers[9]),
             len(covers[1] & covers[9]),
             len(covers[1] - covers[9]),
             len(covers[9] - covers[1])) == (219, 186, 167, 52, 19),
            "small distinguished pair cover split drift")

    # THM-3244 independently verifies the complete target-face family.  Read
    # its pinned exact constant rather than silently duplicating provenance.
    thm3244_syntax = ast.parse(THM3244_SCRIPT.read_text(encoding="utf-8"))
    target_pair_nodes = tuple(
        node.value for node in thm3244_syntax.body
        if isinstance(node, ast.Assign)
        and len(node.targets) == 1
        and isinstance(node.targets[0], ast.Name)
        and node.targets[0].id == "ROW_COVERING_PAIRS"
    )
    require(len(target_pair_nodes) == 1, "THM-3244 pair constant missing")
    require(ast.literal_eval(target_pair_nodes[0])
            == EXPECTED_TARGET_COVERING_PAIRS,
            "THM-3244 target covering-pair drift")

    # Their set-theoretic intersection is the exact cross-support atlas.
    common_covering_pairs = tuple(
        pair for pair in covering_pairs
        if pair in frozenset(EXPECTED_TARGET_COVERING_PAIRS)
    )
    require(common_covering_pairs == EXPECTED_COMMON_COVERING_PAIRS,
            "cross-support covering-pair intersection drift")

    # The adaptive choice is necessary even on the small face: no one fixed
    # positive blend of the distinguished two rows has a Q-directed ascent
    # everywhere.  Equality stays in the trap because ascent is strict.
    small_blend_interval_a = positive_blend_trap_interval(
        SMALL_BLEND_TRAP_A, state_index, small_rows, SMALL_POLES, SMALL_RESET
    )
    small_blend_interval_b = positive_blend_trap_interval(
        SMALL_BLEND_TRAP_B, state_index, small_rows, SMALL_POLES, SMALL_RESET
    )
    require(small_blend_interval_a == (
        Fraction(0), Fraction(23938, 39079645),
    ), "small-ratio small-face blend trap")
    require(small_blend_interval_b == (
        Fraction(66495115323, 16401877394431324), None,
    ), "large-ratio small-face blend trap")
    require(small_blend_interval_b[0] < small_blend_interval_a[1],
            "small-face blend traps lost overlap")
    require(SMALL_BLEND_TRAP_A in covers[1] - covers[9]
            and SMALL_BLEND_TRAP_B in covers[9] - covers[1],
            "small-face blend traps lost exclusive-chart typing")

    witness_rows = []
    for a_value, b_value, state, weight in WITNESS:
        poles = SMALL_POLES if b_value == 2 else FULL_POLES
        reset = SMALL_RESET if b_value == 2 else FULL_RESET
        require(Counter(state) <= Counter(poles) and state != reset,
                "Farkas witness left its nonreset physical bank")
        witness_rows.append((weight, rows[(a_value, b_value, state)]))
    expectations = tuple(
        sum(weight * row[index] for weight, row in witness_rows)
        for index in range(len(T.CERTIFICATE))
    )
    require(sign_census(expectations) == (0, 0, 22),
            "Farkas expectation lost strict coordinate positivity")

    closest_active = max(
        (value, state) for state, value in zip(SMALL_STATES, active_coordinates)
        if value < 0
    )
    farthest_active = min(
        (value, state) for state, value in zip(SMALL_STATES, active_coordinates)
    )
    closest_positive = max(
        (value, state) for state, value in zip(SMALL_STATES, positive_coordinates)
        if value < 0
    )

    state_payload = tuple(
        ",".join(map(str, state)) for state in SMALL_STATES
    )
    witness_payload = tuple(
        f"{a_value},{b_value}|{state}|{weight.numerator}/{weight.denominator}"
        for a_value, b_value, state, weight in WITNESS
    )
    expectation_payload = tuple(
        f"{value.numerator}/{value.denominator}" for value in expectations
    )
    cover_payload = tuple(
        f"{state}|" + "".join("1" if state in cover else "0"
                                for cover in covers)
        for state in SMALL_STATES if state != SMALL_RESET
    )
    pair_payload = tuple(f"{left},{right}" for left, right in covering_pairs)
    common_pair_payload = tuple(
        f"{left},{right}" for left, right in common_covering_pairs
    )
    state_digest = digest(state_payload)
    primitive_digest = digest((",".join(map(str, primitive)),))
    witness_digest = digest(witness_payload)
    expectation_digest = digest(expectation_payload)
    cover_digest = digest(cover_payload)
    pair_digest = digest(pair_payload)
    common_pair_digest = digest(common_pair_payload)
    require(state_digest == EXPECTED_STATE_DIGEST, "small state digest drift")
    require(primitive_digest == EXPECTED_PRIMITIVE_DIGEST,
            "small primitive digest drift")
    require(witness_digest == EXPECTED_WITNESS_DIGEST,
            "Farkas witness digest drift")
    require(expectation_digest == EXPECTED_EXPECTATION_DIGEST,
            "Farkas expectation digest drift")
    require(cover_digest == EXPECTED_SMALL_COVER_DIGEST,
            "small chart-cover digest drift")
    require(pair_digest == EXPECTED_SMALL_PAIR_DIGEST,
            "small covering-pair digest drift")
    require(common_pair_digest == EXPECTED_COMMON_PAIR_DIGEST,
            "cross-support covering-pair digest drift")

    print("dependency_pins=4:PASS")
    print("faces=support_(1,2)_I2_and_support_(1,3)_I2")
    print("small_poles=" + repr(SMALL_POLES))
    print("small_reset=" + repr(SMALL_RESET))
    print("small_states=239:depths=" + repr(EXPECTED_DEPTHS))
    print("small_state_digest=" + state_digest)
    print("transported_thm3238_signs=" + repr(sign_census(raw_coordinates)))
    print("transported_thm3238_top=" + repr(max(
        (value, state) for state, value in zip(SMALL_STATES, raw_coordinates))))
    print("active_rows=" + repr(ACTIVE_ROWS))
    print("active_normalizers=" + repr(ACTIVE_NORMALIZERS))
    print("active_factors=" + repr(ACTIVE_FACTORS))
    print("active_primitive=" + repr(primitive))
    print("active_primitive_bits=%d..%d" % (
        min(value.bit_length() for value in primitive),
        max(value.bit_length() for value in primitive)))
    print("active_primitive_digest=" + primitive_digest)
    print("active_signs=" + repr(sign_census(active_coordinates)))
    print("active_closest=" + repr(closest_active))
    print("active_farthest=" + repr(farthest_active))
    print("positive_22_weights=" + repr(positive_weights))
    print("positive_22_signs=" + repr(sign_census(positive_coordinates)))
    print("positive_22_closest=" + repr(closest_positive))
    print("small_single_chart_best=row16:230/238")
    print("small_covering_pair_count=52")
    print("small_covering_pairs=" + repr(covering_pairs))
    print("small_cover_digest=" + cover_digest)
    print("small_pair_digest=" + pair_digest)
    print("small_pair_(2,10)_split=219,186,167,52,19,0")
    print("small_pair_(2,10)=Q_monotone_one_pole_atlas:PASS")
    print("target_covering_pair_count=31")
    print("cross_support_covering_pair_count=24")
    print("cross_support_covering_pairs=" + repr(common_covering_pairs))
    print("cross_support_pair_digest=" + common_pair_digest)
    print("small_constant_blend_trap_A=(0,%s],state=%s" % (
        str(small_blend_interval_a[1]), repr(SMALL_BLEND_TRAP_A),
    ))
    print("small_constant_blend_trap_B=[%s,infinity),state=%s" % (
        str(small_blend_interval_b[0]), repr(SMALL_BLEND_TRAP_B),
    ))
    print("no_constant_positive_small_row2_row10_blend=PASS,two_states_sharp=PASS")
    print("farkas_witness_count=19:small=3:full=16")
    for index, record in enumerate(WITNESS, 1):
        print("farkas_%02d=" % index + repr(record))
    print("farkas_witness_digest=" + witness_digest)
    print("farkas_expectation_signs=" + repr(sign_census(expectations)))
    print("farkas_expectation_digest=" + expectation_digest)
    print("constant_nonnegative_common_gauge=IMPOSSIBLE_on_these_22_templates")
    print("scope=two_complete_I2_faces_adaptive_atlas_and_no_constant_gauge_only")
    print("status=PASS")


if __name__ == "__main__":
    freeze_support()
    main()

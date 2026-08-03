#!/usr/bin/env python3
"""Exact first-shell and all-covering-pair clutch audit for THM-3254."""

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
THM3238_SCRIPT = (
    ROOT / "04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py"
)
THM3238_OUTPUT = (
    ROOT / "05-knowledge/results/gmc_complete_physical_bank_unique_reset_thm3238.out"
)
THM3244_SCRIPT = (
    ROOT / "04-computation/gmc_unique_reset_rips_nonmorse_thm3244.py"
)
THM3244_OUTPUT = (
    ROOT / "05-knowledge/results/gmc_unique_reset_rips_nonmorse_thm3244.out"
)
DEPENDENCIES = {
    THM3238_SCRIPT:
        "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa",
    THM3238_OUTPUT:
        "77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b",
    THM3244_SCRIPT:
        "3ff0babc41e35e6a185b0ff442cfb9284d9688360c0b96cd947c1128e16400ba",
    THM3244_OUTPUT:
        "27dcd7c68e628465a1f09a564be0be366ded6075ef009a3d85d029f8f18605c9",
}


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


for dependency, expected in DEPENDENCIES.items():
    require(lf_sha256(dependency) == expected,
            ("dependency drift", dependency.name))

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assert_nodes == 0, "optimization-sensitive assert")
require(float_literals == 0, "floating literal")


# Load the immutable THM-3238 exact definitions but not its exhaustive main.
# This reconstructs every row value used below from the upstream coefficient
# formulas; no scratch cache or embedded response value is trusted.
tree = ast.parse(THM3238_SCRIPT.read_text(encoding="utf-8"))
prefix = []
for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name == "main":
        break
    prefix.append(node)
namespace = {
    "__file__": str(THM3238_SCRIPT),
    "__name__": "thm3238_exact_prefix",
}
module = ast.fix_missing_locations(ast.Module(body=prefix, type_ignores=[]))
exec(compile(module, str(THM3238_SCRIPT), "exec"), namespace)

Q = namespace["RESET"]
VALUES = namespace["VALUES"]
MULTIPLICITY = namespace["MULTIPLICITY"]
CERTIFICATE = namespace["CERTIFICATE"]
UPSETS = namespace["UPSETS"]
BANK = namespace["BANK"]
coefficient_vectors = namespace["coefficient_vectors"]
require(Q == (1, 3, 3, 4, 5, 6, 7, 8), "reset drift")
require(len(CERTIFICATE) == 22, "lawful response-row drift")
q_counts = Counter(Q)


ROW_COVERING_PAIRS = (
    (2, 7), (2, 10), (2, 13), (2, 17), (2, 19), (2, 21),
    (3, 9), (3, 11), (3, 16), (3, 22), (7, 14), (7, 18),
    (7, 22), (10, 11), (10, 16), (10, 22), (11, 13), (11, 17),
    (11, 21), (12, 13), (12, 19), (13, 14), (13, 18), (13, 22),
    (14, 19), (16, 17), (16, 21), (17, 22), (18, 19), (19, 22),
    (21, 22),
)

# These are exactly the eight covering pairs whose traps on the distance-one
# link leave a positive ratio gap.  Two explicit global states close each gap.
DELAYED_WITNESSES = {
    (2, 7): (
        (1, 3, 3, 4, 4, 5, 5, 6, 7, 8),
        (6,),
    ),
    (3, 9): (
        (1, 3, 4, 4, 5, 6, 7, 8),
        (8,),
    ),
    (7, 14): (
        (1, 3, 4, 4, 5, 6, 7, 8),
        (8,),
    ),
    (11, 17): (
        (8,),
        (1,),
    ),
    (11, 21): (
        (8,),
        (6,),
    ),
    (12, 13): (
        (3, 3, 4, 4, 5, 5, 6, 7, 8),
        (6,),
    ),
    (12, 19): (
        (3, 3, 4, 4, 5, 5, 6, 7, 8),
        (6,),
    ),
    (14, 19): (
        (1, 3, 3, 4, 4, 5, 5, 6, 7, 8),
        (4,),
    ),
}


def physical(state):
    counts = Counter(state)
    return bool(state) and all(
        counts[value] <= MULTIPLICITY[value] for value in VALUES
    )


def reset_distance(state):
    counts = Counter(state)
    return sum(abs(counts[value] - q_counts[value]) for value in VALUES)


def q_neighbors(state):
    answer = []
    counts = Counter(state)
    for value in VALUES:
        if counts[value] == q_counts[value]:
            continue
        changed = list(state)
        if counts[value] < q_counts[value]:
            changed.append(value)
        else:
            changed.remove(value)
        if changed:
            target = tuple(sorted(changed))
            require(physical(target), "Q-neighbor left physical bank")
            require(reset_distance(target) + 1 == reset_distance(state),
                    "Q-neighbor lost strict distance descent")
            answer.append((value, target))
    return tuple(answer)


# The distance-one link consists of the eleven physical add/delete states.
link_states = tuple(sorted(
    state for state in namespace["STATES"] if reset_distance(state) == 1
))
require(len(link_states) == 11, "distance-one link census")
require(all(q_neighbors(state) == ((
    next(value for value in VALUES
         if Counter(state)[value] != q_counts[value]), Q),)
            for state in link_states), "distance-one edge not unique")

needed_states = {Q, *link_states}
for witnesses in DELAYED_WITNESSES.values():
    for state in witnesses:
        require(physical(state), "nonphysical delayed witness")
        needed_states.add(state)
        needed_states.update(target for _, target in q_neighbors(state))

# Include the particularly transparent row-(2,10) first-link witness chosen
# for the displayed two-inequality Farkas circuit.
PAIR_TWO_TEN_STATES = (
    (3, 3, 4, 5, 6, 7, 8),
    (1, 3, 3, 4, 4, 5, 6, 7, 8),
)
needed_states.update(PAIR_TWO_TEN_STATES)

vectors = {
    state: coefficient_vectors(1, BANK, 1, 3, state)
    for state in sorted(needed_states, key=lambda item: (len(item), item))
}


def response_values(state):
    vector = vectors[state]
    return tuple(
        sum(vector[degree][shape] for shape in upset).numerator
        for (degree, _, _, _), upset in zip(CERTIFICATE, UPSETS)
    )


responses = {state: response_values(state) for state in vectors}
require(responses[Q] == (0,) * 22, "reset response not zero")


def trap_interval(state, left_row, right_row):
    """Positive lambda for which lambda*f_left+f_right has no Q ascent."""
    lower = Fraction(0)
    upper = None
    source = responses[state]
    for _, target in q_neighbors(state):
        target_values = responses[target]
        delta_left = target_values[left_row - 1] - source[left_row - 1]
        delta_right = target_values[right_row - 1] - source[right_row - 1]
        if delta_left > 0:
            bound = Fraction(-delta_right, delta_left)
            upper = bound if upper is None else min(upper, bound)
        elif delta_left < 0:
            lower = max(lower, Fraction(-delta_right, delta_left))
        elif delta_right > 0:
            return None
    lower = max(lower, Fraction(0))
    if upper is not None and (upper < 0 or lower > upper):
        return None
    return lower, upper


def farther(left, right):
    """Whether endpoint left is farther right than endpoint right."""
    if left is None:
        return right is not None
    if right is None:
        return False
    return left > right


def minimum_interval_cover(intervals):
    """Greedy minimum cover of the positive half-line by closed intervals."""
    current = Fraction(0)
    chosen = []
    while True:
        candidates = sorted(
            (state, interval) for state, interval in intervals.items()
            if interval[0] <= current
        )
        if not candidates:
            return None, current
        best_state, best_interval = candidates[0]
        for state, interval in candidates[1:]:
            if farther(interval[1], best_interval[1]):
                best_state, best_interval = state, interval
        if best_interval[1] is None:
            chosen.append(best_state)
            return tuple(chosen), None
        if best_interval[1] <= current:
            return None, current
        chosen.append(best_state)
        current = best_interval[1]


link_covers = {}
link_gap_points = {}
for pair in ROW_COVERING_PAIRS:
    intervals = {
        state: interval for state in link_states
        if (interval := trap_interval(state, *pair)) is not None
    }
    cover, endpoint = minimum_interval_cover(intervals)
    if cover is not None:
        require(len(cover) == 2, ("first-link cover not sharp two", pair))
        link_covers[pair] = cover
        continue
    next_starts = sorted(
        interval[0] for interval in intervals.values()
        if interval[0] > endpoint
    )
    gap_point = (
        (endpoint + next_starts[0]) / 2
        if next_starts else endpoint + 1
    )
    require(gap_point > 0, "nonpositive link gap point")
    require(all(not (interval[0] <= gap_point
                    and (interval[1] is None or gap_point <= interval[1]))
                for interval in intervals.values()), "link gap not open")
    link_gap_points[pair] = gap_point

require(tuple(pair for pair in ROW_COVERING_PAIRS if pair in link_covers)
        == tuple(pair for pair in ROW_COVERING_PAIRS
                 if pair not in DELAYED_WITNESSES),
        "distance-one blocked-pair classification")
require(len(link_covers) == 23 and len(link_gap_points) == 8,
        "distance-one 23/8 split")

# Close the eight link gaps with two exact global witnesses.  In each ordered
# pair, the first interval starts at zero and the second reaches infinity.
delayed_thresholds = {}
for pair, (small_state, large_state) in DELAYED_WITNESSES.items():
    small_interval = trap_interval(small_state, *pair)
    large_interval = trap_interval(large_state, *pair)
    require(small_interval is not None and small_interval[0] == 0
            and small_interval[1] is not None,
            ("bad small-ratio delayed witness", pair))
    require(large_interval is not None and large_interval[1] is None,
            ("bad large-ratio delayed witness", pair))
    require(large_interval[0] <= small_interval[1],
            ("delayed witness intervals lost overlap", pair))
    delayed_thresholds[pair] = (large_interval[0], small_interval[1])

global_witnesses = {
    pair: (link_covers[pair] if pair in link_covers
           else DELAYED_WITNESSES[pair])
    for pair in ROW_COVERING_PAIRS
}
require(all(len(states) == 2 for states in global_witnesses.values()),
        "global witness count not two")

# The row-(2,10) obstruction already occurs on the unique-edge reset link.
state_s, state_t = PAIR_TWO_TEN_STATES
values_s = (responses[state_s][1], responses[state_s][9])
values_t = (responses[state_t][1], responses[state_t][9])
require(values_s == (-93114226186560, 4806299808),
        "row-(2,10) first-shell state S drift")
require(values_t == (16680033356234096640, -48264627879936),
        "row-(2,10) first-shell state T drift")
A, B = -values_s[0], values_s[1]
C, D = values_t[0], -values_t[1]
threshold_s = Fraction(B, A)
threshold_t = Fraction(D, C)
require(threshold_s == Fraction(16688541, 323313285370),
        "state-S threshold drift")
require(threshold_t == Fraction(239733, 82850621920),
        "state-T threshold drift")
threshold_gap = threshold_s - threshold_t
require(threshold_gap
        == Fraction(130514713694581251, 2678670676790293731040) > 0,
        "first-shell threshold order")
farkas_determinant = A * D - C * B
require(farkas_determinant
        == -75675117640279023737068584960 < 0,
        "first-shell Farkas determinant")

witness_digest = hashlib.sha256(
    repr(tuple(global_witnesses.items())).encode("ascii")
).hexdigest()
link_gap_digest = hashlib.sha256(
    repr(tuple(link_gap_points.items())).encode("ascii")
).hexdigest()
delayed_digest = hashlib.sha256(
    repr(tuple(delayed_thresholds.items())).encode("ascii")
).hexdigest()

print("THM-3254 first-shell and all-covering-pair clutch exact control")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("direct_exact_states_reconstructed=%d" % len(vectors))
print("reset_link_states=%d,unique_Q_edges=PASS" % len(link_states))
print("covering_pairs_total=%d" % len(ROW_COVERING_PAIRS))
print("first_link_blocked_delayed=%s" % repr((
    len(link_covers), len(link_gap_points))))
print("first_link_delayed_pairs=%s" % repr(tuple(link_gap_points)))
print("all_31_pairs_two_state_trapped=PASS,one_state_impossible_by_cover=PASS")
print("global_witness_map_sha256=%s" % witness_digest)
print("first_link_gap_points_sha256=%s" % link_gap_digest)
print("delayed_thresholds_sha256=%s" % delayed_digest)
print("row2_row10_first_shell_values=%s" % repr((
    (state_s, values_s), (state_t, values_t))))
print("row2_row10_ascent_thresholds=%s" % repr((
    threshold_s, threshold_t, threshold_gap)))
print("row2_row10_Farkas_determinant=%d" % farkas_determinant)
print("reset_distance_graded_positive_row2_row10_gauge=IMPOSSIBLE")
print("scope=support_(1,3)_bank_I2_lawful_rows_and_Q_monotone_edges_only")
print("all_exact_checks=PASS")

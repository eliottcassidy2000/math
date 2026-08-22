#!/usr/bin/env python3
"""Exact companion for THM-3258's depth-two graded clutch theorem."""

import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
THM3238_SCRIPT = (
    ROOT / "04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py"
)
DEPENDENCIES = {
    ROOT / "01-canon/theorems/THM-3238-complete-physical-product-gamma-bank-unique-reset-stitch.md":
        "3c60ad9f7b74a10df5e6bba5b999e2dffbeb08a8b9d0886bc88a6a923dce1ca1",
    THM3238_SCRIPT:
        "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa",
    ROOT / "05-knowledge/results/gmc_complete_physical_bank_unique_reset_thm3238.out":
        "77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b",
    ROOT / "01-canon/theorems/THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go.md":
        "c7c5948d5181dc845da8f40683517c06e3eecfe837e17cb84d356b91eae1a081",
    ROOT / "04-computation/gmc_first_shell_pair_clutch_thm3254.py":
        "05efd37eeedeca7e3be581977a894592a7873d94a966f06d9533482cc8498fee",
    ROOT / "05-knowledge/results/gmc_first_shell_pair_clutch_thm3254.out":
        "dd415c8ce6e2e196c115421d3508addabb724305a16843509264a8b3205beee9",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected in DEPENDENCIES.items():
    require(sha256(lf_bytes(dependency)).hexdigest() == expected,
            ("dependency hash drift", dependency.name))

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assert_nodes == 0, "optimization-sensitive assert")
require(float_literals == 0, "floating literal")


# Load the immutable THM-3238 definitions without its exhaustive main.  Every
# response used below is recomputed from the original coefficient formulas.
tree = ast.parse(THM3238_SCRIPT.read_text(encoding="utf-8"))
prefix = []
for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name == "main":
        break
    prefix.append(node)
namespace = {
    "__file__": str(THM3238_SCRIPT),
    "__name__": "thm3238_depth_two_exact_prefix",
}
module = ast.fix_missing_locations(ast.Module(body=prefix, type_ignores=[]))
exec(compile(module, str(THM3238_SCRIPT), "exec"), namespace)

RESET = namespace["RESET"]
VALUES = namespace["VALUES"]
MULTIPLICITY = namespace["MULTIPLICITY"]
CERTIFICATE = namespace["CERTIFICATE"]
UPSETS = namespace["UPSETS"]
BANK = namespace["BANK"]
coefficient_vectors = namespace["coefficient_vectors"]
RESET_COUNTS = Counter(RESET)
require(RESET == (1, 3, 3, 4, 5, 6, 7, 8), "reset drift")
require(len(CERTIFICATE) == 22, "response-row drift")


def physical(state):
    counts = Counter(state)
    return bool(state) and all(
        counts[value] <= MULTIPLICITY[value] for value in VALUES
    )


def reset_distance(state):
    counts = Counter(state)
    return sum(abs(counts[value] - RESET_COUNTS[value]) for value in VALUES)


def q_neighbors(state):
    counts = Counter(state)
    answer = []
    for value in VALUES:
        if counts[value] == RESET_COUNTS[value]:
            continue
        changed = list(state)
        if counts[value] < RESET_COUNTS[value]:
            changed.append(value)
        else:
            changed.remove(value)
        target = tuple(sorted(changed))
        if target and physical(target) and reset_distance(target) + 1 == reset_distance(state):
            answer.append(target)
    return tuple(answer)


LINK = tuple(sorted(
    state for state in namespace["STATES"] if reset_distance(state) == 1
))
SHELL_TWO = tuple(sorted(
    state for state in namespace["STATES"] if reset_distance(state) == 2
))
NEIGHBORS = {state: q_neighbors(state) for state in SHELL_TWO}
require(len(LINK) == 11 and len(SHELL_TWO) == 55, "shell census")
require(all(len(targets) in (1, 2) for targets in NEIGHBORS.values()),
        "depth-two target count")
require(sum(len(targets) == 1 for targets in NEIGHBORS.values()) == 3,
        "unique depth-two target census")

needed_states = {RESET} | set(LINK) | set(SHELL_TWO)
for targets in NEIGHBORS.values():
    needed_states.update(targets)


def response_values(state):
    vector = coefficient_vectors(1, BANK, 1, 3, state)
    require(all(vector[degree][shape].denominator == 1
                for degree in range(5, 15)
                for shape in namespace["partitions"](degree)),
            ("nonintegral response vector", state))
    return tuple(
        sum(vector[degree][shape] for shape in upset).numerator
        for (degree, _, _, _), upset in zip(CERTIFICATE, UPSETS)
    )


RESPONSES = {
    state: response_values(state)
    for state in sorted(needed_states, key=lambda item: (len(item), item))
}
require(RESPONSES[RESET] == (0,) * 22, "reset response")


SURVIVING_PAIRS = (
    (2, 7), (3, 9), (7, 14), (11, 17),
    (11, 21), (12, 13), (12, 19), (14, 19),
)


def first_shell_gap(pair):
    """Open ratios lambda for which every distance-one state has an ascent."""
    left_row, right_row = pair
    lower = Fraction(0)
    upper = None
    for state in LINK:
        left = RESPONSES[state][left_row - 1]
        right = RESPONSES[state][right_row - 1]
        # The unique target is RESET, so strict ascent is left*lambda+right<0.
        if left > 0:
            bound = Fraction(-right, left)
            upper = bound if upper is None else min(upper, bound)
        elif left < 0:
            lower = max(lower, Fraction(-right, left))
        else:
            require(right < 0, ("empty link gap", pair, state))
    require(upper is None or lower < upper, ("nonpositive link gap", pair))
    return lower, upper


EXPECTED_GAPS = {
    (2, 7): (Fraction(221630501005, 672431853909), Fraction(1)),
    (3, 9): (Fraction(1, 8), Fraction(5576212, 26946933)),
    (7, 14): (Fraction(1, 55), Fraction(21657540264, 221630501005)),
    (11, 17): (Fraction(554766960, 270430481), Fraction(28, 9)),
    (11, 21): (Fraction(453787039, 270430481), Fraction(3)),
    (12, 13): (Fraction(1420574831, 3464404479), None),
    (12, 19): (Fraction(38075529020, 3464404479), None),
    (14, 19): (Fraction(351232551587, 161618731130),
               Fraction(133113695, 47074217)),
}
require({pair: first_shell_gap(pair) for pair in SURVIVING_PAIRS}
        == EXPECTED_GAPS, "first-shell gap bank")


def ratio_cells(pair):
    """Exact target-envelope cells inside one first-shell ratio gap."""
    left_row, right_row = pair
    lower, upper = EXPECTED_GAPS[pair]
    boundaries = {lower}
    if upper is not None:
        boundaries.add(upper)
    for state in SHELL_TWO:
        targets = NEIGHBORS[state]
        if len(targets) != 2:
            continue
        first, second = targets
        slope = (RESPONSES[first][left_row - 1]
                 - RESPONSES[second][left_row - 1])
        intercept = (RESPONSES[first][right_row - 1]
                     - RESPONSES[second][right_row - 1])
        if slope:
            crossing = Fraction(-intercept, slope)
            if lower < crossing and (upper is None or crossing < upper):
                boundaries.add(crossing)
    ordered = sorted(boundaries)
    cells = list(zip(ordered, ordered[1:]))
    if upper is None:
        cells.append((ordered[-1], None))
    return tuple(cells)


EXPECTED_CELL_COUNTS = {
    (2, 7): 14, (3, 9): 15, (7, 14): 14, (11, 17): 15,
    (11, 21): 15, (12, 13): 18, (12, 19): 14, (14, 19): 4,
}
require({pair: len(ratio_cells(pair)) for pair in SURVIVING_PAIRS}
        == EXPECTED_CELL_COUNTS, "raw ratio-cell census")
require(sum(EXPECTED_CELL_COUNTS.values()) == 109, "raw cell total")


def shell_label(state, target):
    return "shell2", state, target


# These fourteen small supports were found by discovery, but no discovered
# coefficient is trusted: exact response rows, positive null weights, target
# dominance intervals and negative affine constants are all rebuilt below.
CIRCUIT_LABELS = {
    (2, 7): (
        (("ratio-upper", Fraction(121588653639453734947671,
                                  356989167296855203779847)),
         ("positive-b2",),
         shell_label((1, 1, 1, 3, 3, 4, 5, 6, 7, 8),
                     (1, 1, 3, 3, 4, 5, 6, 7, 8)),
         shell_label((1, 3, 3, 4, 4, 5, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
        (shell_label((1, 3, 3, 4, 5, 6), (1, 3, 3, 4, 5, 6, 8)),
         shell_label((1, 3, 3, 5, 5, 6, 7, 8),
                     (1, 3, 3, 4, 5, 5, 6, 7, 8)),
         shell_label((2, 3, 3, 4, 5, 6, 7, 8),
                     (3, 3, 4, 5, 6, 7, 8)),
         shell_label((1, 3, 3, 4, 4, 5, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
        (shell_label((1, 3, 3, 4, 5, 6), (1, 3, 3, 4, 5, 6, 8)),
         shell_label((1, 3, 3, 5, 5, 6, 7, 8),
                     (1, 3, 3, 4, 5, 5, 6, 7, 8)),
         shell_label((2, 3, 3, 4, 5, 6, 7, 8),
                     (1, 2, 3, 3, 4, 5, 6, 7, 8)),
         shell_label((1, 1, 1, 3, 3, 4, 5, 6, 7, 8),
                     (1, 1, 3, 3, 4, 5, 6, 7, 8))),
    ),
    (3, 9): (
        (("positive-a2",),
         shell_label((1, 3, 3, 4, 5, 6), (1, 3, 3, 4, 5, 6, 8)),
         shell_label((1, 3, 3, 5, 5, 6, 7, 8),
                     (1, 3, 3, 4, 5, 5, 6, 7, 8)),
         shell_label((1, 3, 4, 4, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
        (("positive-b2",),
         shell_label((1, 3, 3, 4, 5, 6), (1, 3, 3, 4, 5, 6, 7)),
         shell_label((1, 3, 3, 4, 5, 8), (1, 3, 3, 4, 5, 7, 8)),
         shell_label((1, 3, 4, 4, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
    ),
    (7, 14): (
        (("positive-a2",), ("positive-b2",),
         shell_label((1, 3, 3, 4, 5, 6), (1, 3, 3, 4, 5, 6, 8)),
         shell_label((1, 3, 4, 4, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
        (("ratio-lower", Fraction(33500762305293296528429,
                                  364765960918361204843013)),
         ("positive-b2",),
         shell_label((1, 4, 5, 6, 7, 8), (1, 3, 4, 5, 6, 7, 8)),
         shell_label((1, 3, 4, 4, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
    ),
    (11, 17): (
        (("positive-a2",),
         shell_label((1, 3, 3, 4, 5, 6), (1, 3, 3, 4, 5, 6, 7)),
         shell_label((1, 4, 5, 6, 7, 8), (1, 3, 4, 5, 6, 7, 8)),
         shell_label((1, 3, 4, 4, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
        (("positive-a2",), ("positive-b2",),
         shell_label((1, 3, 3, 4, 5, 6), (1, 3, 3, 4, 5, 6, 8)),
         shell_label((1, 3, 4, 4, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
    ),
    (11, 21): (
        (("positive-a2",),
         shell_label((1, 3, 3, 4, 5, 6), (1, 3, 3, 4, 5, 6, 7)),
         shell_label((1, 4, 5, 6, 7, 8), (1, 3, 4, 5, 6, 7, 8)),
         shell_label((1, 3, 4, 4, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
        (("positive-b2",),
         shell_label((1, 3, 3, 4, 5, 6), (1, 3, 3, 4, 5, 6, 8)),
         shell_label((1, 4, 5, 6, 7, 8), (1, 3, 4, 5, 6, 7, 8)),
         shell_label((1, 3, 4, 4, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
    ),
    (12, 13): (
        (("positive-a2",),
         shell_label((1, 1, 1, 3, 3, 4, 5, 6, 7, 8),
                     (1, 1, 3, 3, 4, 5, 6, 7, 8)),
         shell_label((1, 3, 3, 4, 4, 5, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
    ),
    (12, 19): (
        (("positive-a2",),
         shell_label((1, 3, 4, 4, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8)),
         shell_label((1, 3, 3, 4, 4, 5, 5, 6, 7, 8),
                     (1, 3, 3, 4, 4, 5, 6, 7, 8))),
    ),
    (14, 19): (
        (("ratio-lower", Fraction(351232551587, 161618731130)),
         shell_label((1, 3, 3, 4, 5, 6), (1, 3, 3, 4, 5, 6, 8)),
         shell_label((1, 3, 3, 5, 5, 6, 7, 8),
                     (1, 3, 3, 4, 5, 5, 6, 7, 8)),
         shell_label((1, 3, 4, 5, 5, 6, 7, 8),
                     (1, 3, 3, 4, 5, 5, 6, 7, 8))),
    ),
}
require(sum(map(len, CIRCUIT_LABELS.values())) == 14, "circuit census")
require(Counter(len(labels) for bank in CIRCUIT_LABELS.values()
                for labels in bank) == {3: 2, 4: 12},
        "circuit support sizes")


def row_from_label(pair, label):
    left_row, right_row = pair
    if label[0] == "positive-a2":
        return (0, 1, 0), Fraction(0)
    if label[0] == "positive-b2":
        return (0, 0, 1), Fraction(0)
    if label[0] == "ratio-lower":
        return (1, 0, 0), -label[1]
    if label[0] == "ratio-upper":
        return (-1, 0, 0), label[1]
    require(label[0] == "shell2", ("unknown label", label))
    state, target = label[1], label[2]
    require(state in NEIGHBORS and target in NEIGHBORS[state],
            ("illegal directed shell edge", pair, state, target))
    return (
        (RESPONSES[target][left_row - 1],
         -RESPONSES[state][left_row - 1],
         -RESPONSES[state][right_row - 1]),
        Fraction(RESPONSES[target][right_row - 1]),
    )


def positive_null_relation(vectors):
    """Primitive positive relation among three or four vectors in Q^3."""
    width = len(vectors)
    matrix = [[Fraction(vectors[column][row]) for column in range(width)]
              for row in range(3)]
    pivot_columns = []
    pivot_row = 0
    for column in range(width):
        selected = next((row for row in range(pivot_row, 3)
                         if matrix[row][column]), None)
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        pivot = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot for value in matrix[pivot_row]]
        for row in range(3):
            if row == pivot_row or not matrix[row][column]:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [left - multiplier * right
                           for left, right in zip(matrix[row], matrix[pivot_row])]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == 3:
            break
    free_columns = [column for column in range(width)
                    if column not in pivot_columns]
    require(len(free_columns) == 1, ("circuit nullity", vectors, free_columns))
    free = free_columns[0]
    weights = [Fraction(0) for _ in range(width)]
    weights[free] = 1
    for row, column in enumerate(pivot_columns):
        weights[column] = -matrix[row][free]
    denominator = 1
    for value in weights:
        denominator = lcm(denominator, value.denominator)
    integers = [int(value * denominator) for value in weights]
    divisor = 0
    for value in integers:
        divisor = gcd(divisor, abs(value))
    integers = [value // divisor for value in integers]
    if all(value < 0 for value in integers):
        integers = [-value for value in integers]
    require(all(value > 0 for value in integers),
            ("nonpositive circuit weights", vectors, integers))
    require(all(sum(integers[index] * vectors[index][coordinate]
                    for index in range(width)) == 0
                for coordinate in range(3)), "circuit vector cancellation")
    return tuple(integers)


def intersect_interval(lower, upper, new_lower, new_upper):
    lower = max(lower, new_lower)
    if upper is None:
        upper = new_upper
    elif new_upper is not None:
        upper = min(upper, new_upper)
    return lower, upper


def dominance_interval(pair, state, target):
    """Ratios where target attains the maximum over Q-neighbors of state."""
    left_row, right_row = pair
    lower, upper = EXPECTED_GAPS[pair]
    for other in NEIGHBORS[state]:
        if other == target:
            continue
        slope = (RESPONSES[target][left_row - 1]
                 - RESPONSES[other][left_row - 1])
        intercept = (RESPONSES[target][right_row - 1]
                     - RESPONSES[other][right_row - 1])
        if slope > 0:
            lower = max(lower, Fraction(-intercept, slope))
        elif slope < 0:
            bound = Fraction(-intercept, slope)
            upper = bound if upper is None else min(upper, bound)
        else:
            require(intercept >= 0, ("never-maximal target", pair, state, target))
    return lower, upper


def circuit_domain(pair, labels):
    lower, upper = EXPECTED_GAPS[pair]
    for label in labels:
        if label[0] == "ratio-lower":
            lower = max(lower, label[1])
        elif label[0] == "ratio-upper":
            upper = label[1] if upper is None else min(upper, label[1])
        elif label[0] == "shell2":
            lower, upper = intersect_interval(
                lower, upper, *dominance_interval(pair, label[1], label[2])
            )
    require(upper is None or lower <= upper,
            ("empty circuit domain", pair, labels, lower, upper))
    return lower, upper


certificates = []
for pair in SURVIVING_PAIRS:
    domains = []
    for labels in CIRCUIT_LABELS[pair]:
        rows = tuple(row_from_label(pair, label) for label in labels)
        weights = positive_null_relation(tuple(row[0] for row in rows))
        constant = sum(weight * row[1]
                       for weight, row in zip(weights, rows))
        require(constant < 0, ("nonnegative Farkas constant", pair, labels))
        require(any(label[0] == "shell2" and weight > 0
                    for label, weight in zip(labels, weights)),
                "circuit has no strict shell inequality")
        domain = circuit_domain(pair, labels)
        domains.append(domain)
        certificates.append((pair, domain, labels, weights, constant))

    # The declared circuits form an exact interval chain across the complete
    # open first-shell gap.  Shared endpoints are safe because shell ascent
    # inequalities remain strict even when a ratio-bound row is zero.
    gap_lower, gap_upper = EXPECTED_GAPS[pair]
    require(domains[0][0] == gap_lower, ("circuit chain start", pair, domains))
    current = domains[0][1]
    for lower, upper in domains[1:]:
        require(current is not None and lower <= current,
                ("circuit chain gap", pair, current, lower))
        if current is None or upper is None:
            current = None
        else:
            current = max(current, upper)
    require((gap_upper is None and current is None)
            or (gap_upper is not None and current is not None
                and current >= gap_upper),
            ("circuit chain end", pair, current, gap_upper))

    # Independent raw-envelope coverage: each of the 109 exact open cells is
    # contained in at least one certified domain.
    for cell_lower, cell_upper in ratio_cells(pair):
        require(any(domain_lower <= cell_lower
                    and (domain_upper is None
                         or (cell_upper is not None and cell_upper <= domain_upper))
                    for domain_lower, domain_upper in domains),
                ("uncovered raw ratio cell", pair, cell_lower, cell_upper))

require(len(certificates) == 14, "final certificate count")
certificate_digest = sha256("\n".join(map(repr, certificates)).encode("ascii")).hexdigest()
require(certificate_digest ==
        "6e0cbf932f99cea022dcf9fa81bbaa2adec0f6251e27b5a72fcbd32aca6a9263",
        "circuit certificate digest")


# A compact human-readable control: the whole unbounded (12,13) gap is killed
# by one three-row circuit, whose exact weights and constant are small enough
# to display in the theorem.
simple = next(item for item in certificates if item[0] == (12, 13))
require(simple[3] == (4452809708057982010202880,
                      4277107928796831, 371631008),
        "simple circuit weights")
require(simple[4] == -9167765286556241877685416,
        "simple circuit constant")

print("THM-3258 DEPTH-TWO GRADED PAIR CLUTCH EXACT AUDIT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("direct_exact_states_reconstructed=%d" % len(RESPONSES))
print("reset_link_states=%d,depth_two_states=%d,target_counts=(3*1,52*2)" %
      (len(LINK), len(SHELL_TWO)))
print("first_link_surviving_pairs=%s" % (SURVIVING_PAIRS,))
print("exact_first_link_gap_bank=%s" %
      (tuple((pair, EXPECTED_GAPS[pair]) for pair in SURVIVING_PAIRS),))
print("raw_target_envelope_cells=%s,total=%d" %
      (tuple(EXPECTED_CELL_COUNTS[pair] for pair in SURVIVING_PAIRS),
       sum(EXPECTED_CELL_COUNTS.values())))
print("compressed_affine_Farkas_circuits=(3,2,2,2,2,1,1,1),total=14")
print("circuit_support_histogram=((3,2),(4,12)),all_weights_positive,all_constants_negative")
print("circuit_certificate_sha256=%s" % certificate_digest)
print("simple_pair_(12,13)_weights=%s,constant=%d" %
      (simple[3], simple[4]))
print("all_31_covering_pairs_reset_distance_graded_positive_gauge=IMPOSSIBLE_by_depth<=2")
print("scope=support_(1,3)_bank_I2_lawful_rows_Q_monotone_edges_two-row_grading_only")
print("all_exact_checks=PASS")

#!/usr/bin/env python3
"""Exact threshold-layer multicover anatomy of the z216 status prefix.

The maintained projected-k3 screen assigns one common four-bit status table
to every load threshold.  This scout ignores the noncanonical Farkas vectors
returned by the inherited floating search.  It reconstructs each deterministic
status instance and looks instead for an integral modular majorant

    sum_(t in T) 1[c(P) >= t] <= sum_(i in P) w_i                 (1)

on all sixteen patterns P.  Summing (1) against a common status table gives

    sum_(t in T) H_t <= sum_i w_i m_i,                            (2)

where H_t is the required load tail and m_i is the i-th marginal.  A strict
reverse inequality is therefore an exact, solver-free contradiction.

The script also proves tail-support minimality.  For every smaller threshold
support it constructs an exact feasible common table by deterministic rational
vertex enumeration.  This rules out *all* Farkas certificates on that smaller
tail support, not merely binary or nonnegative-coordinate certificates.

Everything here remains inside the projected-k3, z1=216 necessary atlas.  No
physical speed entry, owner, phase, current, arbitrary-k, rung, or LRC(14)
consequence is asserted.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import itertools
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py"
)
AUDIT_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.out"
)
AUDIT_SOURCE_SHA256 = (
    "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3"
)
AUDIT_OUTPUT_SHA256 = (
    "e035531d0af15998a8a07d083042ff2486c80c7fdd729d902922179409d1725a"
)
AUDIT_SEMANTIC_SHA256 = (
    "4bb4eca9949eb4f07b618f9e49bbb5de63492364a29a4a5abfc84a6acf105e11"
)

LEVEL = 216
SELECTED_INDICES = (
    134, 121, 17, 27, 46, 56, 79, 89, 126, 138, 205, 248, 258, 286,
    298, 347, 376, 425,
)
EXPECTED_STATUS_BY_ROW = (
    (134, 1), (121, 16), (17, 13), (27, 18), (46, 0), (56, 24),
    (79, 0), (89, 0), (126, 0), (138, 98), (205, 3), (248, 0),
    (258, 0), (286, 2), (298, 89), (347, 3), (376, 0), (425, 4),
)
EXPECTED_LAYER_COUNTS = ((1, 260), (2, 10), (3, 1))
EXPECTED_OLD_TO_LAYER = (
    ("coordinate_union", ((1, 248),)),
    ("two_fan", ((1, 11),)),
    ("weighted_core", ((2, 10), (3, 1))),
    ("zero_reduced_union", ((1, 1),)),
)
EXPECTED_MULTI_COUNT = 11
EXPECTED_UNIQUE_MULTI_COUNT = 8
EXPECTED_PROPER_SUPPORT_ADDRESSES = 59
EXPECTED_DISTINCT_PROPER_SUPPORTS = 47
EXPECTED_SUBMODULAR_COUNTS = ((False, 3), (True, 8))
EXPECTED_SUPERMODULAR_COUNTS = ((False, 11),)

EXPECTED_CERTIFICATE_SHA256 = (
    "9636ec6677f20c0cd317bf8d4f9ad9455f6f237334e83faa5c139f155fff398c"
)
EXPECTED_MULTI_SHA256 = (
    "fe4627c932651bca1e2f6776db42d472fa1536a467e92863bde9cecfdab0c37f"
)
EXPECTED_MINIMALITY_SHA256 = (
    "e3830599d2987136ac501be6cc435751d6f396b51b9b9db8a7a9e3080fa098a2"
)
EXPECTED_HOSTILE_SHA256 = (
    "b5531108a2deca1e4042c7ae25c4b5d9f0de52904ffc759c255af176eda6f3cd"
)
EXPECTED_SEMANTIC_SHA256 = (
    "bd7ac011d960afe7c4566ac9c17e9df3d8c2cf710e6023a715cb6a6d1ce1408d"
)

HOSTILE_INDEX = 64
HOSTILE_COUNTS = (115, 50, 57, 8)
HOSTILE_DIVISORS = (80, 2156, 4312, 5390)
HOSTILE_Q = 43_120
HOSTILE_MARGINALS = (6468, 6160, 6160, 6160)
HOSTILE_CAPACITIES = (
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
)
HOSTILE_HISTOGRAM = ((0, 28_046), (1, 15_074))
HOSTILE_TABLE = (
    28_046, 0, 0, 3234, 0, 3234, 2446, 0,
    5680, 0, 0, 0, 0, 0, 480, 0,
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
    return hashlib.sha256(payload).hexdigest()


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tail_data(capacities, histogram):
    thresholds = []
    demands = {}
    good_rows = {}
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        good = tuple(int(capacity >= threshold) for capacity in capacities)
        if all(good):
            continue
        thresholds.append(threshold)
        demands[threshold] = sum(
            count for load_value, count in histogram
            if load_value >= threshold
        )
        good_rows[threshold] = good
    return tuple(thresholds), demands, good_rows


def layer_values(capacities, thresholds):
    return tuple(
        sum(capacities[pattern] >= threshold for threshold in thresholds)
        for pattern in range(16)
    )


def dominates_layers(weights, values):
    return all(
        sum(
            weights[coordinate] * ((pattern >> coordinate) & 1)
            for coordinate in range(4)
        )
        >= values[pattern]
        for pattern in range(16)
    )


def minimal_binary_multicovers(capacities, marginals, histogram):
    """All strict certificates at the least binary tail-support size.

    With k selected layers the layer count is at most k.  Hence every
    nonnegative integral coordinate weight may canonically be truncated to
    [0,k], making the search finite without loss inside this template.
    """
    thresholds, demands, _good_rows = tail_data(capacities, histogram)
    for support_size in range(1, len(thresholds) + 1):
        certificates = []
        for selected in itertools.combinations(thresholds, support_size):
            values = layer_values(capacities, selected)
            tail_sum = sum(demands[threshold] for threshold in selected)
            for weights in itertools.product(
                range(support_size + 1), repeat=4
            ):
                marginal_cost = sum(
                    weight * marginal
                    for weight, marginal in zip(weights, marginals)
                )
                if tail_sum <= marginal_cost:
                    continue
                if not dominates_layers(weights, values):
                    continue
                certificates.append(
                    (selected, weights, tail_sum - marginal_cost)
                )
        if certificates:
            return support_size, tuple(certificates)
    return None, ()


def tight_patterns(capacities, thresholds, weights):
    values = layer_values(capacities, thresholds)
    return tuple(
        pattern for pattern in range(16)
        if sum(
            weights[coordinate] * ((pattern >> coordinate) & 1)
            for coordinate in range(4)
        )
        == values[pattern]
    )


def coordinate_minimal(capacities, thresholds, weights):
    values = layer_values(capacities, thresholds)
    for coordinate, weight in enumerate(weights):
        if not weight:
            continue
        lowered = tuple(
            value - int(index == coordinate)
            for index, value in enumerate(weights)
        )
        if dominates_layers(lowered, values):
            return False
    return True


def modularity_status(capacities, thresholds):
    values = layer_values(capacities, thresholds)
    sub_witness = None
    super_witness = None
    for left in range(16):
        for right in range(16):
            defect = (
                values[left] + values[right]
                - values[left | right] - values[left & right]
            )
            if defect < 0 and sub_witness is None:
                sub_witness = (left, right, defect)
            if defect > 0 and super_witness is None:
                super_witness = (left, right, defect)
    return sub_witness is None, super_witness is None, sub_witness, super_witness


def solve_square(rows, rhs):
    size = len(rows)
    require(all(len(row) == size for row in rows), "nonsquare system")
    matrix = [
        [Fraction(value) for value in row] + [Fraction(rhs[index])]
        for index, row in enumerate(rows)
    ]
    for column in range(size):
        pivot = next(
            (
                row for row in range(column, size)
                if matrix[row][column]
            ),
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


def equality_rows():
    return (
        (1,) * 16,
        *tuple(
            tuple((pattern >> coordinate) & 1 for pattern in range(16))
            for coordinate in range(4)
        ),
    )


@lru_cache(maxsize=None)
def marginal_vertices(q, marginals):
    rows = equality_rows()
    rhs = (q, *marginals)
    vertices = []
    for columns in itertools.combinations(range(16), 5):
        square = tuple(
            tuple(row[column] for column in columns) for row in rows
        )
        solution = solve_square(square, rhs)
        if solution is None or any(value < 0 for value in solution):
            continue
        table = [Fraction(0)] * 16
        for column, value in zip(columns, solution):
            table[column] = value
        vertices.append(tuple(table))
    require(vertices, (q, marginals, "no marginal vertex"))
    return tuple(vertices)


def verify_table(q, marginals, good_rows, demands, selected, table):
    if len(table) != 16 or any(value < 0 for value in table):
        return False
    rows = equality_rows()
    rhs = (q, *marginals)
    if any(
        sum(coefficient * value for coefficient, value in zip(row, table))
        != target
        for row, target in zip(rows, rhs)
    ):
        return False
    return all(
        sum(
            coefficient * value
            for coefficient, value in zip(good_rows[threshold], table)
        )
        >= demands[threshold]
        for threshold in selected
    )


def exact_feasible_table(q, marginals, capacities, histogram, selected):
    """Find a lexicographically first exact vertex satisfying selected tails."""
    _thresholds, demands, good_rows = tail_data(capacities, histogram)
    for table in marginal_vertices(q, tuple(marginals)):
        if verify_table(
            q, marginals, good_rows, demands, selected, table
        ):
            return table

    base_rows = equality_rows()
    base_rhs = (q, *marginals)
    for active_count in range(1, len(selected) + 1):
        for active in itertools.combinations(selected, active_count):
            rows = base_rows + tuple(good_rows[threshold] for threshold in active)
            rhs = base_rhs + tuple(demands[threshold] for threshold in active)
            size = len(rows)
            for columns in itertools.combinations(range(16), size):
                square = tuple(
                    tuple(row[column] for column in columns) for row in rows
                )
                solution = solve_square(square, rhs)
                if solution is None or any(value < 0 for value in solution):
                    continue
                table = [Fraction(0)] * 16
                for column, value in zip(columns, solution):
                    table[column] = value
                table = tuple(table)
                if verify_table(
                    q, marginals, good_rows, demands, selected, table
                ):
                    return table
    return None


def old_branch(audit, capacities, marginals, histogram):
    unions = audit.union_cuts(capacities, marginals, histogram)
    reduced = () if unions else audit.zero_reduced_union_cuts(
        capacities, marginals, histogram
    )
    fans = () if (unions or reduced) else audit.two_fan_cuts(
        capacities, marginals, histogram
    )
    if unions:
        return "coordinate_union"
    if reduced:
        return "zero_reduced_union"
    if fans:
        return "two_fan"
    return "weighted_core"


def format_multi(record):
    return (
        f"row:{record[0]};divisors:{record[1]};q:{record[2]};"
        f"thresholds:{record[8]};weights:{record[9]};gap:{record[10]};"
        f"tight:{record[11]};submodular:{int(record[12])};"
        f"submodular_witness:{record[14]}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    require(lf_sha(AUDIT_SOURCE) == AUDIT_SOURCE_SHA256, "audit source changed")
    require(lf_sha(AUDIT_OUTPUT) == AUDIT_OUTPUT_SHA256, "audit output changed")
    require(
        f"semantic_sha256={AUDIT_SEMANTIC_SHA256}"
        in AUDIT_OUTPUT.read_text(encoding="utf-8"),
        "audit semantic changed",
    )
    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assert_nodes, float_nodes) == (0, 0), "truth-gate syntax")

    audit = load("threshold_layer_audit", AUDIT_SOURCE)
    natural = audit.load("threshold_layer_natural", audit.NATURAL_SOURCE)
    rows, _components = natural.atlas_rows()
    require(len(rows) == 480, len(rows))

    certificate_records = []
    multi_records = []
    old_to_layer = defaultdict(Counter)
    status_by_row = []
    for index in SELECTED_INDICES:
        eng = audit.status_engine(natural, str(index))
        eng.FIRST = LEVEL
        eng.ray.FIRST = LEVEL
        body, ruler, high, wall = rows[index]
        stream = eng.ray.Stream(body)
        require((stream.L, stream.high_floor) == (ruler, high), (index, "atlas"))
        _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
        _crude, status, residual = eng.exact_common_status_screen(stream, states)
        require(wall and not residual, (index, "selected row did not close"))
        status_by_row.append((index, len(status)))

        for divisors, witness in sorted(status.items()):
            q, cofactor, marginals, capacity_set, histogram, _raw_dual = witness
            divisor_lcm = lcm(*divisors)
            require(divisor_lcm == q * cofactor, (index, divisors, "cofactor"))
            rebuilt_marginals, capacities = eng.ray.local.hunter_status_data(
                divisor_lcm, divisors, q
            )
            require(rebuilt_marginals == marginals, (index, divisors, "marginals"))
            require(
                tuple(sorted(set(capacities))) == capacity_set,
                (index, divisors, "capacities"),
            )
            support_size, certificates = minimal_binary_multicovers(
                capacities, marginals, histogram
            )
            require(support_size is not None and certificates, (index, divisors))
            selected, weights, gap = certificates[0]
            require(
                coordinate_minimal(capacities, selected, weights),
                (index, divisors, "nonminimal coordinate weights"),
            )
            require(
                gcd(*tuple(value for value in (*weights, *([1] * len(selected))) if value))
                == 1,
                (index, divisors, "nonprimitive certificate"),
            )
            tight = tight_patterns(capacities, selected, weights)
            submodular, supermodular, sub_witness, super_witness = (
                modularity_status(capacities, selected)
            )
            branch = old_branch(audit, capacities, marginals, histogram)
            old_to_layer[branch][support_size] += 1
            record = (
                index,
                divisors,
                q,
                cofactor,
                marginals,
                histogram,
                tuple(capacities),
                branch,
                selected,
                weights,
                gap,
                tight,
                submodular,
                supermodular,
                sub_witness,
                super_witness,
            )
            certificate_records.append(record)
            if support_size > 1:
                require(len(certificates) == 1, (index, divisors, certificates))
                multi_records.append(record)

    status_by_row = tuple(status_by_row)
    certificate_records = tuple(certificate_records)
    multi_records = tuple(multi_records)
    require(status_by_row == EXPECTED_STATUS_BY_ROW, status_by_row)
    require(len(certificate_records) == 271, len(certificate_records))
    layer_counts = tuple(sorted(Counter(len(row[8]) for row in certificate_records).items()))
    require(layer_counts == EXPECTED_LAYER_COUNTS, layer_counts)
    old_to_layer_packet = tuple(
        (branch, tuple(sorted(counts.items())))
        for branch, counts in sorted(old_to_layer.items())
    )
    require(old_to_layer_packet == EXPECTED_OLD_TO_LAYER, old_to_layer_packet)
    require(len(multi_records) == EXPECTED_MULTI_COUNT, len(multi_records))

    unique = {}
    address_multiplicity = Counter()
    for record in multi_records:
        key = (record[2], record[4], record[5], record[6])
        unique.setdefault(key, record)
        address_multiplicity[key] += 1
    require(len(unique) == EXPECTED_UNIQUE_MULTI_COUNT, len(unique))

    minimality_records = []
    for key, record in sorted(unique.items(), key=lambda item: repr(item[0])):
        q, marginals, histogram, capacities = key
        thresholds, _demands, _good_rows = tail_data(capacities, histogram)
        minimum = len(record[8])
        for support_size in range(1, minimum):
            for selected in itertools.combinations(thresholds, support_size):
                table = exact_feasible_table(
                    q, marginals, capacities, histogram, selected
                )
                require(table is not None, (record[0], selected, "infeasible proper support"))
                minimality_records.append((key, selected, table))
    minimality_records = tuple(minimality_records)
    require(
        len(minimality_records) == EXPECTED_DISTINCT_PROPER_SUPPORTS,
        len(minimality_records),
    )
    addressed_proper_supports = sum(
        (
            sum(
                len(tuple(itertools.combinations(
                    tail_data(record[6], record[5])[0], support_size
                )))
                for support_size in range(1, len(record[8]))
            )
        )
        for record in multi_records
    )
    require(
        addressed_proper_supports == EXPECTED_PROPER_SUPPORT_ADDRESSES,
        addressed_proper_supports,
    )

    submodular_counts = tuple(sorted(Counter(row[12] for row in multi_records).items()))
    supermodular_counts = tuple(sorted(Counter(row[13] for row in multi_records).items()))
    require(submodular_counts == EXPECTED_SUBMODULAR_COUNTS, submodular_counts)
    require(supermodular_counts == EXPECTED_SUPERMODULAR_COUNTS, supermodular_counts)

    # Hostile control: the old row 64 has 57 genuine status kills, all caught
    # by one-layer covers, but its first screen survivor admits an explicit
    # common table at q=43120 and must not be killed by this template.
    index = HOSTILE_INDEX
    eng = audit.status_engine(natural, "hostile64")
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(rows[index][0])
    _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
    crude, status, residual = eng.exact_common_status_screen(stream, states)
    hostile_counts = (len(states), len(crude), len(status), len(residual))
    require(hostile_counts == HOSTILE_COUNTS, hostile_counts)
    hostile_status_packet = []
    for divisors, witness in sorted(status.items()):
        q, _cofactor, marginals, _capacity_set, histogram, _raw_dual = witness
        divisor_lcm = lcm(*divisors)
        _rebuilt, capacities = eng.ray.local.hunter_status_data(
            divisor_lcm, divisors, q
        )
        support_size, certificates = minimal_binary_multicovers(
            capacities, marginals, histogram
        )
        require(support_size == 1 and certificates, (index, divisors))
        hostile_status_packet.append((divisors, certificates[0]))
    require(HOSTILE_DIVISORS in residual, residual)
    divisor_lcm = lcm(*HOSTILE_DIVISORS)
    require(divisor_lcm == HOSTILE_Q, divisor_lcm)
    arcs = eng.ray.fibre.projected_support_arcs(divisor_lcm, stream.ranges)
    hostile_marginals, hostile_capacities = eng.ray.local.hunter_status_data(
        divisor_lcm, HOSTILE_DIVISORS, HOSTILE_Q
    )
    hostile_histogram = eng.ray.fibre.residue_load_histogram(arcs, HOSTILE_Q)
    require(hostile_marginals == HOSTILE_MARGINALS, hostile_marginals)
    require(tuple(hostile_capacities) == HOSTILE_CAPACITIES, hostile_capacities)
    require(hostile_histogram == HOSTILE_HISTOGRAM, hostile_histogram)
    hostile_thresholds, hostile_demands, hostile_good = tail_data(
        hostile_capacities, hostile_histogram
    )
    require(hostile_thresholds == (1,), hostile_thresholds)
    require(
        verify_table(
            HOSTILE_Q,
            HOSTILE_MARGINALS,
            hostile_good,
            hostile_demands,
            hostile_thresholds,
            tuple(map(Fraction, HOSTILE_TABLE)),
        ),
        "hostile common table invalid",
    )
    hostile_support, hostile_certificates = minimal_binary_multicovers(
        hostile_capacities, hostile_marginals, hostile_histogram
    )
    require(
        hostile_support is None and not hostile_certificates,
        "hostile survivor falsely killed",
    )

    certificate_sha = hashlib.sha256(repr(certificate_records).encode()).hexdigest()
    multi_sha = hashlib.sha256(repr(multi_records).encode()).hexdigest()
    minimality_sha = hashlib.sha256(repr(minimality_records).encode()).hexdigest()
    hostile_packet = (
        hostile_counts,
        tuple(hostile_status_packet),
        HOSTILE_DIVISORS,
        HOSTILE_Q,
        HOSTILE_MARGINALS,
        HOSTILE_CAPACITIES,
        HOSTILE_HISTOGRAM,
        HOSTILE_TABLE,
    )
    hostile_sha = hashlib.sha256(repr(hostile_packet).encode()).hexdigest()
    if EXPECTED_CERTIFICATE_SHA256 is not None:
        require(certificate_sha == EXPECTED_CERTIFICATE_SHA256, certificate_sha)
    if EXPECTED_MULTI_SHA256 is not None:
        require(multi_sha == EXPECTED_MULTI_SHA256, multi_sha)
    if EXPECTED_MINIMALITY_SHA256 is not None:
        require(minimality_sha == EXPECTED_MINIMALITY_SHA256, minimality_sha)
    if EXPECTED_HOSTILE_SHA256 is not None:
        require(hostile_sha == EXPECTED_HOSTILE_SHA256, hostile_sha)

    multi_gap_counts = tuple(sorted(Counter(row[10] for row in multi_records).items()))
    multiplicity_counts = tuple(sorted(Counter(address_multiplicity.values()).items()))
    semantic_packet = (
        AUDIT_SOURCE_SHA256,
        AUDIT_OUTPUT_SHA256,
        AUDIT_SEMANTIC_SHA256,
        status_by_row,
        layer_counts,
        old_to_layer_packet,
        certificate_sha,
        multi_records,
        multi_sha,
        multiplicity_counts,
        addressed_proper_supports,
        minimality_sha,
        submodular_counts,
        supermodular_counts,
        hostile_sha,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected-k3 z216 threshold-layer multicover dual anatomy",
        (
            f"dependency=cost_prefix_audit_source:{AUDIT_SOURCE_SHA256};"
            f"output:{AUDIT_OUTPUT_SHA256};semantic:{AUDIT_SEMANTIC_SHA256}"
        ),
        (
            f"universe=z1:{LEVEL};selected_rows:{len(SELECTED_INDICES)};"
            f"status_instances:{len(certificate_records)};"
            f"status_by_row:{status_by_row}"
        ),
        (
            "template=sum_selected_tail_indicators_is_pointwise_at_most_"
            "a_zero_constant_nonnegative_integral_modular_coordinate_cover;"
            "summing_against_one_common_status_table_bounds_the_sum_of_tail_"
            "demands_by_the_weighted_marginals"
        ),
        (
            f"layer_taxonomy={layer_counts};old_to_layer:{old_to_layer_packet};"
            f"certificate_sha256:{certificate_sha}"
        ),
        (
            f"multi_layer=addresses:{len(multi_records)};unique_instances:"
            f"{len(unique)};unique_multiplicity:{multiplicity_counts};"
            f"gap_counts:{multi_gap_counts};multi_sha256:{multi_sha}"
        ),
    ]
    lines.extend(f"multi={format_multi(record)}" for record in multi_records)
    lines.extend(
        [
            (
                f"tail_support_minimality=addressed_proper_supports:"
                f"{addressed_proper_supports};distinct_proper_supports:"
                f"{len(minimality_records)};all_have_exact_feasible_common_"
                f"tables:PASS;minimality_sha256:{minimality_sha}"
            ),
            (
                "primitive_normalization=every_multi_layer_certificate_has_"
                "binary_tail_weights_zero_constant_coordinatewise_minimal_"
                "nonnegative_integer_coordinate_weights_and_gcd_one;the_"
                "inherited_solver_vector_and_its_scale_are_not_used_or_hashed"
            ),
            (
                f"submodularity_test=selected_layer_functions_submodular:"
                f"{submodular_counts};supermodular:{supermodular_counts};"
                "row27_and_row56_packets_supply_exact_submodularity_violations_"
                "while_the_modular_majorants_remain_valid"
            ),
            (
                "first_failed_bold_implication=nested_tail_events_do_not_force_"
                "a_two_layer_certificate;the_row138_divisors_18_49_196_882_"
                "instance_has_exact_feasible_common_tables_for_every_single_"
                "and_pair_of_its_six_tail_rows_but_the_three_layer_support_"
                "1_4_5_has_gap_11"
            ),
            (
                f"hostile=row:{HOSTILE_INDEX};screen_counts:{hostile_counts};"
                "all_57_status_kills_have_one_layer_covers;first_residual_"
                f"divisors:{HOSTILE_DIVISORS};q:{HOSTILE_Q};explicit_common_"
                f"table:PASS;no_binary_multicover_contradiction:PASS;"
                f"hostile_sha256:{hostile_sha}"
            ),
            (
                "mechanism=the_obstruction_is_cross_threshold_compatibility_"
                "of_one_joint_status_table;each_proper_tail_subsystem_can_pass_"
                "even_when_the_selected_layer_sum_cannot;this_is_a_threshold_"
                "chain_multicover_not_a_coordinate_union_identity_and_not_"
                "a_submodularity_argument"
            ),
            (
                "direction=each_strict_modular_cover_gap_is_an_exact_Farkas_"
                "contradiction_for_the_deterministic_common_status_instance;"
                "it_recertifies_only_the_271_already_selected_projected_status_"
                "kills_and_changes_no_row_ledger"
            ),
            (
                "scope=projected_k3_z216_necessary_atlas_only;no_physical_"
                "speed_entry_endpoint_owner_phase_current_arbitrary_k_le_1_"
                "rung_or_LRC14_claim"
            ),
            f"truth_gates=assert_nodes:{assert_nodes};float_literals:{float_nodes}",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    payload = "\n".join(lines) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()

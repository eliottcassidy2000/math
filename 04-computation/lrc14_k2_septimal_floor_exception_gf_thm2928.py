#!/usr/bin/env python3
"""Complete k=2/p=5 q=D/7 floor/exception fractional-cover GF.

For every support-hard row put q=D/7.  A denominator not dividing q gives
one point uniformly on every q-fibre.  For each remaining denominator d|q,
write d=7a+r and w=q/d.  Its exact number of hit fibres is

    a*w + w*X_d(u),       Pr(X_d=1)=r/7.

This program retains every nonconstant bit (w,r), grants the deterministic
floor sum, and applies the exact upward-event fractional-cover maximum with
the compact/open strict endpoint >66/91.  A divisor-Mobius completion count
groups all uniform symbols and all r=0 spike symbols, so no five-tuple is
materialized.  The finite arity-four weighted-state engine is exact rational
arithmetic; this Python referee reconstructs the literal body rows, baseline
expected-spike stage, bounded literal controls, and all survivor ledgers.
"""

from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, lcm
from pathlib import Path
from subprocess import PIPE, run
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
COMBINED_PATH = (
    ROOT / "04-computation" / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
EXPECTED_COMBINED_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
ENGINE_PATH = Path(__file__).with_name(
    "lrc14_k2_septimal_floor_exception_engine_thm2928.cpp"
)
EXPECTED_ENGINE_SHA256 = (
    "664d0df36d104d959279605c8ea8539d61ab595b155e5157fa7d0433f1b7944c"
)
EXPECTED_ENGINE_OUTPUT_SHA256 = (
    "49cbf8fb160a78125d426cf7348c3c35e4859f110f493bb4b03acfc2f0c92125"
)
SUPPORT_CUTOFF = Q(887, 990)
TWO_SAFE_FLOOR = Q(66, 91)
EXPECTED_ROWS = 27163
EXPECTED_DIVISORS = 219
EXPECTED_RAW_SHAPES = 50874159718
EXPECTED_RAW_OCCURRENCES = 951545890235
EXPECTED_NONSEPTIMAL_FULL_ROWS = 96235
LITERAL_CONTROL_MAX_D = 300
EXPECTED_STAGE = {
    "raw_shapes_by_c": {
        1: 17780752220,
        2: 18817476992,
        3: 10575488581,
        4: 3250143848,
        5: 450298077,
    },
    "raw_occurrences_by_c": {
        1: 300597339936,
        2: 342425871732,
        3: 217083433126,
        4: 78454285495,
        5: 12984959946,
    },
    "expected_shapes": 36962285549,
    "expected_occurrences": 320011786356,
    "expected_shapes_by_c": {
        1: 4129900738,
        2: 18792244265,
        3: 10568748763,
        4: 3246262178,
        5: 225129605,
    },
    "expected_occurrences_by_c": {
        1: 13975689268,
        2: 143128657438,
        3: 119436903764,
        4: 41980227729,
        5: 1490308157,
    },
    "exact_c4": (3208397602, 39762283189),
    "expected_support": (4592, 2999, 149),
    "expected_rows_by_c": {1: 575, 2: 1151, 3: 2712, 4: 3609, 5: 2409},
    "expected_bodies_by_c": {1: 575, 2: 1148, 3: 2572, 4: 2940, 5: 1244},
    "expected_divisors_by_c": {1: 85, 2: 87, 3: 118, 4: 134, 5: 91},
}
EXPECTED_EXACT = {
    "literal_controls": (18, 52925, 1880),
    "shapes": 26908162790,
    "occurrences": 200389247292,
    "shapes_by_c": {
        1: 3457885,
        2: 12908211669,
        3: 10562966029,
        4: 3208397602,
        5: 225129605,
    },
    "occurrences_by_c": {
        1: 3524756,
        2: 46320209782,
        3: 112812921408,
        4: 39762283189,
        5: 1490308157,
    },
    "support": (4414, 2977, 149),
    "rows_by_c": {1: 45, 2: 779, 3: 2425, 4: 3373, 5: 2409},
    "bodies_by_c": {1: 45, 2: 779, 3: 2337, 4: 2871, 5: 1244},
    "divisors_by_c": {1: 15, 2: 71, 3: 115, 4: 133, 5: 91},
    "minimum_D_by_c": {1: 13860, 2: 840, 3: 336, 4: 168, 5: 392},
    "semantic_sha256": (
        "2eec9a97f02a7b8f8e36e50f747d53186ff5b84234a14fc3ac64818b54033675"
    ),
    "minimum_survivor": (
        168,
        (1, 2, 3, 4, 6, 12),
        168,
        88,
        ((3, 14), (4, 6), (5, 2), (6, 2)),
        4,
        4,
        4,
        204,
    ),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(file_sha256(COMBINED_PATH) == EXPECTED_COMBINED_SHA256, "tracked dependency changed")
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


@lru_cache(maxsize=None)
def c_shape_distribution(D):
    """Exact lcm-D shape count by c=#{d not dividing D/7}."""
    q = D // 7
    result = Counter()
    for E in support.divisors(D):
        sign = mobius(D // E)
        if not sign:
            continue
        alphabet = tuple(d for d in support.divisors(E) if d > 1)
        uniform = sum(q % d != 0 for d in alphabet)
        spike = len(alphabet) - uniform
        for c in range(1, 6):
            result[c] += (
                sign
                * multichoose(uniform, c)
                * multichoose(spike, 5 - c)
            )
    require(all(value >= 0 for value in result.values()), ("negative c GF", D))
    result += Counter()
    return result


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
def maximum_upward_mass(weights, marginals, need):
    weights = tuple(weights)
    marginals = tuple(marginals)
    count = len(weights)
    if need <= 0:
        return Q(1)
    if need > sum(weights):
        return Q(0)
    minimal = []
    for mask in range(1 << count):
        weight = sum(
            weights[index]
            for index in range(count)
            if (mask >> index) & 1
        )
        if weight < need:
            continue
        if any(
            weight - weights[index] >= need
            for index in range(count)
            if (mask >> index) & 1
        ):
            continue
        minimal.append(tuple((mask >> index) & 1 for index in range(count)))
    require(minimal, "nonempty literal-control event has no minimal state")
    constraints = [(row, Q(1)) for row in minimal]
    constraints.extend(
        (
            tuple(int(index == coordinate) for index in range(count)),
            Q(0),
        )
        for coordinate in range(count)
    )
    optimum = None
    for chosen in combinations(range(len(constraints)), count):
        rows = tuple(constraints[index][0] for index in chosen)
        rhs = tuple(constraints[index][1] for index in chosen)
        point = solve_square(rows, rhs)
        if point is None or any(value < 0 for value in point):
            continue
        if any(
            sum(value * coefficient for value, coefficient in zip(point, row)) < 1
            for row in minimal
        ):
            continue
        objective = sum(
            marginal * value
            for marginal, value in zip(marginals, point)
        )
        if optimum is None or objective < optimum:
            optimum = objective
    require(optimum is not None, "literal-control fractional cover has no vertex")
    return min(Q(1), optimum)


@lru_cache(maxsize=None)
def status_allowance(weights, marginals):
    weights = tuple(weights)
    marginals = tuple(marginals)
    if not weights:
        return 0
    thresholds = sorted(
        {
            sum(
                weights[index]
                for index in range(len(weights))
                if (mask >> index) & 1
            )
            for mask in range(1 << len(weights))
        }
    )
    answer = 0
    for threshold in thresholds[1:]:
        if maximum_upward_mass(weights, marginals, threshold) > TWO_SAFE_FLOOR:
            answer = threshold
        else:
            break
    return answer


def literal_allowance(D, values):
    q = D // 7
    uniform = tuple(d for d in values if q % d)
    spikes = tuple(d for d in values if q % d == 0)
    deterministic = sum((d // 7) * (q // d) for d in spikes)
    weights = tuple(q // d for d in spikes if d % 7)
    marginals = tuple(Q(d % 7, 7) for d in spikes if d % 7)
    return len(uniform), deterministic + status_allowance(weights, marginals)


def build_rows():
    rows_by_D = defaultdict(list)
    body_count = 0
    body_divisor_rows = 0
    nonseptimal_rows = 0
    nonseptimal_full = 0
    for body in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support.safe_cell_ranges(body)
        for D in support.divisors(L):
            body_divisor_rows += 1
            support_count = support.support_size_bitset(D, ranges)
            if D % 7:
                nonseptimal_rows += 1
                nonseptimal_full += support_count == D
            if Q(support_count, D) > SUPPORT_CUTOFF:
                continue
            require(D % 7 == 0, ("support-hard row is not septimal", body, D))
            q = D // 7
            arcs = combined.projected_support_arcs(D, ranges)
            histogram = combined.residue_load_histogram(arcs, q)
            N = tuple(
                sum(count for load, count in histogram if load > c)
                for c in range(1, 6)
            )
            rows_by_D[D].append((body, L, support_count, histogram, N))
    require(body_count == 3003, "body count changed")
    require(body_divisor_rows == 251536, "body/divisor universe changed")
    require(
        (nonseptimal_rows, nonseptimal_full)
        == (EXPECTED_NONSEPTIMAL_FULL_ROWS,) * 2,
        "nonseptimal full-support census changed",
    )
    require(sum(map(len, rows_by_D.values())) == EXPECTED_ROWS, "support row count changed")
    require(len(rows_by_D) == EXPECTED_DIVISORS, "support divisor count changed")
    return rows_by_D, nonseptimal_rows


def engine_queries(rows_by_D):
    result = {}
    for D, rows in rows_by_D.items():
        q = D // 7
        queries = {
            (c, N[c - 1])
            for _body, _L, _support_count, _histogram, N in rows
            for c in range(1, 6)
        }
        if D <= LITERAL_CONTROL_MAX_D:
            queries.update((c, need) for c in range(1, 6) for need in range(q + 1))
        result[D] = tuple(sorted(queries))
    return result


def run_engine(queries):
    input_lines = [str(len(queries))]
    for D in sorted(queries):
        input_lines.append(f"{D} {len(queries[D])}")
        input_lines.extend(f"{c} {need}" for c, need in queries[D])
    input_text = "\n".join(input_lines) + "\n"
    engine_outputs = []
    with TemporaryDirectory(prefix="lrc14-k2-allc-") as temporary:
        for optimization in ("-O2", "-O3"):
            executable = Path(temporary) / f"engine-{optimization[1:]}"
            compile_result = run(
                [
                    "/usr/bin/clang++",
                    "-std=c++17",
                    optimization,
                    "-DNDEBUG",
                    str(ENGINE_PATH),
                    "-o",
                    str(executable),
                ],
                stdout=PIPE,
                stderr=PIPE,
                text=True,
                check=False,
            )
            require(
                compile_result.returncode == 0,
                ("engine compilation failed", optimization, compile_result.stderr),
            )
            engine_result = run(
                [str(executable)],
                input=input_text,
                stdout=PIPE,
                stderr=PIPE,
                text=True,
                check=False,
            )
            require(
                engine_result.returncode == 0,
                ("engine execution failed", optimization, engine_result.stderr),
            )
            engine_outputs.append(engine_result.stdout)
    require(engine_outputs[0] == engine_outputs[1], "O2/O3 engine outputs differ")
    engine_output = engine_outputs[0]
    raw = {}
    answers = {}
    engine_lines = engine_output.splitlines()
    require(engine_lines[0] == "ENGINE all-c exact floor/exception GF", "engine banner changed")
    require(engine_lines[-1] == "all_exact_engine_controls=PASS", "engine did not pass")
    for line in engine_lines:
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "C":
            _tag, D, c, count, state_count, minimum, maximum = fields
            raw[(int(D), int(c))] = (
                int(count),
                int(state_count),
                int(minimum),
                int(maximum),
            )
        elif fields[0] == "Q":
            _tag, D, c, need, count, first = fields
            answers[(int(D), int(c), int(need))] = (int(count), int(first))
    require(len(raw) == 5 * len(queries), "engine raw record count changed")
    require(
        len(answers) == sum(map(len, queries.values())),
        "engine query record count changed",
    )
    return raw, answers, engine_output


def bounded_literal_controls(rows_by_D, raw, answers):
    cases = 0
    shapes = 0
    queries_checked = 0
    for D in sorted(D for D in rows_by_D if D <= LITERAL_CONTROL_MAX_D):
        q = D // 7
        alphabet = tuple(d for d in support.divisors(D) if d > 1)
        distribution = Counter()
        for values in combinations_with_replacement(alphabet, 5):
            if lcm(*values) != D:
                continue
            distribution[literal_allowance(D, values)] += 1
            shapes += 1
        for c in range(1, 6):
            raw_count = sum(
                count for (local_c, _allowance), count in distribution.items()
                if local_c == c
            )
            require(raw[(D, c)][0] == raw_count, ("bounded raw mismatch", D, c))
            for need in range(q + 1):
                expected = sum(
                    count
                    for (local_c, allowance), count in distribution.items()
                    if local_c == c and allowance >= need
                )
                possible = sorted(
                    allowance
                    for (local_c, allowance), count in distribution.items()
                    if local_c == c and count and allowance >= need
                )
                expected_first = possible[0] if possible else -1
                require(
                    answers[(D, c, need)] == (expected, expected_first),
                    ("bounded exact query mismatch", D, c, need),
                )
                queries_checked += 1
        cases += 1
    return cases, shapes, queries_checked


def main():
    rows_by_D, nonseptimal_rows = build_rows()
    queries = engine_queries(rows_by_D)
    raw, answers, engine_output = run_engine(queries)
    literal_cases, literal_shapes, literal_queries = bounded_literal_controls(
        rows_by_D, raw, answers
    )

    raw_shapes = 0
    raw_occurrences = 0
    raw_shapes_by_c = Counter()
    raw_occurrences_by_c = Counter()
    expected_shapes = 0
    expected_occurrences = 0
    expected_shapes_by_c = Counter()
    expected_occurrences_by_c = Counter()
    expected_rows = set()
    expected_bodies = set()
    expected_divisors = set()
    expected_rows_by_c = {c: set() for c in range(1, 6)}
    expected_bodies_by_c = {c: set() for c in range(1, 6)}
    expected_divisors_by_c = {c: set() for c in range(1, 6)}
    exact_shapes = 0
    exact_occurrences = 0
    exact_shapes_by_c = Counter()
    exact_occurrences_by_c = Counter()
    exact_rows = set()
    exact_bodies = set()
    exact_divisors = set()
    exact_rows_by_c = {c: set() for c in range(1, 6)}
    exact_bodies_by_c = {c: set() for c in range(1, 6)}
    exact_divisors_by_c = {c: set() for c in range(1, 6)}
    semantic = sha256()
    minimum_survivor = None
    minimum_survivor_by_c = {c: None for c in range(1, 6)}
    exact_occurrences_by_D_c = Counter()
    expected_violation = None

    for D in sorted(rows_by_D):
        q = D // 7
        c_distribution = c_shape_distribution(D)
        require(
            all(raw[(D, c)][0] == c_distribution[c] for c in range(1, 6)),
            ("engine c coefficient mismatch", D),
        )
        shapes_D = sum(c_distribution.values())
        raw_shapes += shapes_D
        raw_occurrences += shapes_D * len(rows_by_D[D])
        for c, count in c_distribution.items():
            raw_shapes_by_c[c] += count
            raw_occurrences_by_c[c] += count * len(rows_by_D[D])

        minimum_N = {
            c: min(row[4][c - 1] for row in rows_by_D[D])
            for c in range(1, 6)
        }
        for c in range(1, 6):
            # Expected-spike shape existence is governed by the smallest N_c.
            N = minimum_N[c]
            expected_pass = N == 0 or 66 * N < 13 * (5 - c) * q
            if expected_pass:
                expected_shapes += c_distribution[c]
                expected_shapes_by_c[c] += c_distribution[c]
            exact_count, _first = answers[(D, c, N)]
            exact_shapes += exact_count
            exact_shapes_by_c[c] += exact_count

        for body, L, support_count, histogram, N_values in rows_by_D[D]:
            row_has_exact = False
            row_has_expected = False
            for c, N in enumerate(N_values, 1):
                expected_pass = N == 0 or 66 * N < 13 * (5 - c) * q
                if expected_pass:
                    row_has_expected = True
                    expected_occurrences += c_distribution[c]
                    expected_occurrences_by_c[c] += c_distribution[c]
                    expected_rows_by_c[c].add((body, D))
                    expected_bodies_by_c[c].add(body)
                    expected_divisors_by_c[c].add(D)
                count, first = answers[(D, c, N)]
                if count and not expected_pass and expected_violation is None:
                    expected_violation = (D, body, c, N, count)
                if not count:
                    continue
                row_has_exact = True
                exact_occurrences += count
                exact_occurrences_by_c[c] += count
                exact_occurrences_by_D_c[(D, c)] += count
                exact_rows_by_c[c].add((body, D))
                exact_bodies_by_c[c].add(body)
                exact_divisors_by_c[c].add(D)
                record = (
                    D,
                    body,
                    L,
                    support_count,
                    histogram,
                    c,
                    N,
                    first,
                    count,
                )
                semantic.update(f"{record}\n".encode())
                if minimum_survivor is None or record < minimum_survivor:
                    minimum_survivor = record
                if (
                    minimum_survivor_by_c[c] is None
                    or record < minimum_survivor_by_c[c]
                ):
                    minimum_survivor_by_c[c] = record
            if row_has_expected:
                expected_rows.add((body, D))
                expected_bodies.add(body)
                expected_divisors.add(D)
            if row_has_exact:
                exact_rows.add((body, D))
                exact_bodies.add(body)
                exact_divisors.add(D)

    require(raw_shapes == EXPECTED_RAW_SHAPES, "raw shape total changed")
    require(raw_occurrences == EXPECTED_RAW_OCCURRENCES, "raw occurrence total changed")
    require(dict(raw_shapes_by_c) == EXPECTED_STAGE["raw_shapes_by_c"], "raw shapes-by-c changed")
    require(
        dict(raw_occurrences_by_c) == EXPECTED_STAGE["raw_occurrences_by_c"],
        "raw occurrences-by-c changed",
    )
    require(expected_shapes == EXPECTED_STAGE["expected_shapes"], "expected shape stage changed")
    require(
        expected_occurrences == EXPECTED_STAGE["expected_occurrences"],
        "expected occurrence stage changed",
    )
    require(
        dict(expected_shapes_by_c) == EXPECTED_STAGE["expected_shapes_by_c"],
        "expected shapes-by-c changed",
    )
    require(
        dict(expected_occurrences_by_c)
        == EXPECTED_STAGE["expected_occurrences_by_c"],
        "expected occurrences-by-c changed",
    )
    require(expected_violation is None, ("exact screen failed to dominate mean", expected_violation))
    require(
        (exact_shapes_by_c[4], exact_occurrences_by_c[4])
        == EXPECTED_STAGE["exact_c4"],
        "all-c engine does not reproduce exact one-spike c4 stage",
    )

    engine_sha = sha256(engine_output.encode()).hexdigest()
    require(file_sha256(ENGINE_PATH) == EXPECTED_ENGINE_SHA256, "engine source changed")
    require(engine_sha == EXPECTED_ENGINE_OUTPUT_SHA256, "engine output changed")
    require(
        (literal_cases, literal_shapes, literal_queries)
        == EXPECTED_EXACT["literal_controls"],
        "bounded literal control census changed",
    )
    require(
        (len(expected_rows), len(expected_bodies), len(expected_divisors))
        == EXPECTED_STAGE["expected_support"],
        "expected support census changed",
    )
    require(
        {c: len(value) for c, value in expected_rows_by_c.items()}
        == EXPECTED_STAGE["expected_rows_by_c"],
        "expected rows-by-c changed",
    )
    require(
        {c: len(value) for c, value in expected_bodies_by_c.items()}
        == EXPECTED_STAGE["expected_bodies_by_c"],
        "expected bodies-by-c changed",
    )
    require(
        {c: len(value) for c, value in expected_divisors_by_c.items()}
        == EXPECTED_STAGE["expected_divisors_by_c"],
        "expected divisors-by-c changed",
    )
    require(exact_shapes == EXPECTED_EXACT["shapes"], "exact shape total changed")
    require(exact_occurrences == EXPECTED_EXACT["occurrences"], "exact occurrence total changed")
    require(dict(exact_shapes_by_c) == EXPECTED_EXACT["shapes_by_c"], "exact shapes-by-c changed")
    require(
        dict(exact_occurrences_by_c) == EXPECTED_EXACT["occurrences_by_c"],
        "exact occurrences-by-c changed",
    )
    require(
        (len(exact_rows), len(exact_bodies), len(exact_divisors))
        == EXPECTED_EXACT["support"],
        "exact support census changed",
    )
    require(
        {c: len(value) for c, value in exact_rows_by_c.items()}
        == EXPECTED_EXACT["rows_by_c"],
        "exact rows-by-c changed",
    )
    require(
        {c: len(value) for c, value in exact_bodies_by_c.items()}
        == EXPECTED_EXACT["bodies_by_c"],
        "exact bodies-by-c changed",
    )
    require(
        {c: len(value) for c, value in exact_divisors_by_c.items()}
        == EXPECTED_EXACT["divisors_by_c"],
        "exact divisors-by-c changed",
    )
    require(
        {c: value[0] for c, value in minimum_survivor_by_c.items()}
        == EXPECTED_EXACT["minimum_D_by_c"],
        "minimum exact divisor by c changed",
    )
    require(
        semantic.hexdigest() == EXPECTED_EXACT["semantic_sha256"],
        "exact semantic digest changed",
    )
    require(
        minimum_survivor == EXPECTED_EXACT["minimum_survivor"],
        "minimum exact survivor changed",
    )
    print("LRC14 k=2/p=5 complete q=D/7 floor/exception GF")
    print(f"combined_script_sha256={file_sha256(COMBINED_PATH)}")
    print(f"engine_source_sha256={file_sha256(ENGINE_PATH)}")
    print(f"engine_output_sha256={engine_sha}")
    print(f"nonseptimal_full_support_rows={nonseptimal_rows}/{nonseptimal_rows}")
    print(f"support_rows={sum(map(len, rows_by_D.values()))}")
    print(f"support_divisors={len(rows_by_D)}")
    print(f"bounded_literal_D_cases={literal_cases}")
    print(f"bounded_literal_shapes={literal_shapes}")
    print(f"bounded_literal_suffix_queries={literal_queries}")
    print("engine_O2_O3_byte_identical=PASS")
    print(f"raw_shapes={raw_shapes}")
    print(f"raw_occurrences={raw_occurrences}")
    print(f"raw_shapes_by_c={raw_shapes_by_c}")
    print(f"raw_occurrences_by_c={raw_occurrences_by_c}")
    print(f"expected_shapes={expected_shapes}")
    print(f"expected_occurrences={expected_occurrences}")
    print(f"expected_shapes_by_c={expected_shapes_by_c}")
    print(f"expected_occurrences_by_c={expected_occurrences_by_c}")
    print(f"expected_rows={len(expected_rows)}")
    print(f"expected_bodies={len(expected_bodies)}")
    print(f"expected_divisors={len(expected_divisors)}")
    print(f"expected_rows_by_c={Counter({c: len(v) for c, v in expected_rows_by_c.items()})}")
    print(f"expected_bodies_by_c={Counter({c: len(v) for c, v in expected_bodies_by_c.items()})}")
    print(f"expected_divisors_by_c={Counter({c: len(v) for c, v in expected_divisors_by_c.items()})}")
    print(f"exact_shapes={exact_shapes}")
    print(f"exact_occurrences={exact_occurrences}")
    print(f"exact_shapes_by_c={exact_shapes_by_c}")
    print(f"exact_occurrences_by_c={exact_occurrences_by_c}")
    print(f"exact_rows={len(exact_rows)}")
    print(f"exact_bodies={len(exact_bodies)}")
    print(f"exact_divisors={len(exact_divisors)}")
    print(f"exact_rows_by_c={Counter({c: len(v) for c, v in exact_rows_by_c.items()})}")
    print(f"exact_bodies_by_c={Counter({c: len(v) for c, v in exact_bodies_by_c.items()})}")
    print(f"exact_divisors_by_c={Counter({c: len(v) for c, v in exact_divisors_by_c.items()})}")
    print(
        "exact_c1_occurrences_by_D="
        f"{Counter({D: count for (D, c), count in exact_occurrences_by_D_c.items() if c == 1})}"
    )
    print(f"killed_from_expected_shapes={expected_shapes - exact_shapes}")
    print(f"killed_from_expected_occurrences={expected_occurrences - exact_occurrences}")
    print(f"semantic_sha256={semantic.hexdigest()}")
    print(f"minimum_survivor={minimum_survivor}")
    print(f"minimum_survivor_by_c={minimum_survivor_by_c}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()

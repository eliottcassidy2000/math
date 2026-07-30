#!/usr/bin/env python3
"""Located-phase closure of the k=2,c=4 sole-d=6 septimal subfamily.

Universe: the exact support-hard body/modulus rows and exact-lcm unordered
five-denominator multisets of THM-2928's k=2 floor/exception ledger.  Here
c=4 means four transverse sections and one vertical denominator.  If that
vertical denominator is d=6 and N_4>0, coverage would force its high bit on
throughout the compact two-aligned safe carrier.

For reduced numerator s in {1,5}, the high event is

    V_s = {u : ||s*u|| < 3/7}.

Its low set C_s is the radius-1/14 comb translated by half a tooth.  An exact
Bernoulli/Fourier overlap law proves that two distinct aligned danger combs
cover at most 3/28 of C_s, whereas |C_s|=1/7.  Thus the forced containment is
impossible.  This referee verifies the overlap formula by exact interval
arithmetic, counts fixed-d=6 completions by divisor-Mobius inversion, checks
bounded literal multisets, and recomputes the full residual-union ledger.

The imported floor/exception referee compiles its required C++ engine under
-O2 and -O3 and requires byte-identical engine output.  This script itself is
also intended to be replayed under ordinary and optimized Python.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations_with_replacement
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/lrc14_k2_septimal_floor_exception_gf_thm2928.py"
EXPECTED_BASE_SHA256 = (
    "085b4e2747a48bdbc1125e894af7d4f647dfdd7be86a00cf02dea2a8667e26dc"
)
EXPECTED_H_TABLES = {
    1: (0, -1, 5, -3, 3, -5, 1, 0, -1, 5, -3, 3, -5, 1),
    5: (0, -5, 4, -8, 8, -4, 5, 0, -5, 4, -8, 8, -4, 5),
}
EXPECTED = {
    "literal_overlap_controls": 1000,
    "literal_lcm_controls": (18, 1029),
    "participating_shapes": 18064772,
    "deleted_shapes": 8998004,
    "deleted_occurrences": 248154771,
    "touched": (2406, 2148, 105),
    "killed_semantic": (
        "a14c202ee92936958eab161b5aad857c5ce8f243b27a0db0778c5d62ef437d89"
    ),
    "residual": (26899164786, 200141092521, 4354, 2966, 147),
    "residual_shapes_by_c": {
        1: 3457885,
        2: 12908211669,
        3: 10562966029,
        4: 3199399598,
        5: 225129605,
    },
    "residual_occurrences_by_c": {
        1: 3524756,
        2: 46320209782,
        3: 112812921408,
        4: 39514128418,
        5: 1490308157,
    },
    "residual_rows_by_c": {1: 45, 2: 779, 3: 2425, 4: 3185, 5: 2409},
    "residual_bodies_by_c": {1: 45, 2: 779, 3: 2337, 4: 2812, 5: 1244},
    "residual_divisors_by_c": {1: 15, 2: 71, 3: 115, 4: 130, 5: 91},
    "residual_minimum": (
        336,
        (1, 2, 3, 4, 8, 12),
        336,
        160,
        ((2, 14), (3, 14), (4, 14), (5, 2), (6, 4)),
        4,
        6,
        2095,
    ),
    "residual_minimum_D_by_c": {1: 13860, 2: 840, 3: 336, 4: 336, 5: 392},
    "residual_semantic": (
        "c4eca31da58ee8db462970519e45edb9c514fd43a9642c2b8e925a8d9987b7e9"
    ),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(file_sha256(BASE_PATH) == EXPECTED_BASE_SHA256, "base referee changed")
spec = spec_from_file_location("lrc14_k2_septimal_base", BASE_PATH)
base = module_from_spec(spec)
spec.loader.exec_module(base)


def bernoulli2(value):
    value %= 1
    return value * value - value + Q(1, 6)


def h_table(S):
    return tuple(
        int(
            49
            * (
                bernoulli2(Q((S + 6 * A) % 14, 14))
                - bernoulli2(Q((S + 8 * A) % 14, 14))
            )
        )
        for A in range(14)
    )


def overlap_formula(a, s):
    """Haar measure of C_s intersect D_a, with null endpoints immaterial."""
    g = gcd(a, s)
    A = a // g
    S = s // g
    require(S in (1, 5), ("unexpected reduced spike numerator", a, s, A, S))
    h = EXPECTED_H_TABLES[S][A % 14]
    return Q(1, 49) + Q(h, 49 * A * S)


def literal_overlap(a, s):
    """Independent exact interval intersection of C_s with D_a."""
    result = Q(0)
    danger_radius = Q(1, 14 * a)
    for cell in range(s):
        low_left = (Q(cell) + Q(3, 7)) / s
        low_right = (Q(cell) + Q(4, 7)) / s
        for tooth in range(a + 1):
            center = Q(tooth, a)
            overlap = min(low_right, center + danger_radius) - max(
                low_left, center - danger_radius
            )
            if overlap > 0:
                result += overlap
    return result


def overlap_controls():
    require(
        {S: h_table(S) for S in (1, 5)} == EXPECTED_H_TABLES,
        "Bernoulli residue table changed",
    )
    controls = 0
    for s in (1, 5):
        for a in range(1, 501):
            require(
                literal_overlap(a, s) == overlap_formula(a, s),
                ("literal/Fourier overlap mismatch", a, s),
            )
            controls += 1

    # The residue tables make these bounds uniform, not merely bounded tests.
    # For h>0 the correction h/(49*A*S) decreases within each residue class.
    for s in (1, 5):
        candidates = []
        for a in range(1, 351):  # one complete lcm(14,5) tail plus slack
            value = overlap_formula(a, s)
            candidates.append((value, a))
        maximum = max(candidates)
        require(maximum == (Q(1, 14), 2 * s), ("sharp maximum changed", s, maximum))
        second = max(value for value, a in candidates if a != 2 * s)
        require(second == Q(1, 28), ("second overlap bound changed", s, second))
    require(Q(1, 14) + Q(1, 28) < Q(1, 7), "located two-comb margin vanished")
    return controls


def fixed_d6_shape_count(D):
    """Exact-lcm count with four transverse symbols and one fixed d=6."""
    q = D // 7
    if q % 6:
        return 0
    answer = 0
    for E in base.support.divisors(D):
        if E % 6:
            continue
        sign = base.mobius(D // E)
        if not sign:
            continue
        transverse = sum(
            q % d != 0 for d in base.support.divisors(E) if d > 1
        )
        answer += sign * base.multichoose(transverse, 4)
    require(answer >= 0, ("negative fixed-d6 coefficient", D, answer))
    return answer


def literal_lcm_controls(rows_by_D):
    cases = 0
    shapes = 0
    for D in sorted(D for D in rows_by_D if D <= base.LITERAL_CONTROL_MAX_D):
        q = D // 7
        literal = 0
        alphabet = tuple(d for d in base.support.divisors(D) if d > 1)
        for values in combinations_with_replacement(alphabet, 5):
            if lcm(*values) != D:
                continue
            vertical = tuple(d for d in values if q % d == 0)
            if vertical == (6,):
                literal += 1
        require(
            literal == fixed_d6_shape_count(D),
            ("literal fixed-d6 completion mismatch", D, literal),
        )
        cases += 1
        shapes += literal
    return cases, shapes


def main():
    overlap_count = overlap_controls()
    require(
        overlap_count == EXPECTED["literal_overlap_controls"],
        "literal overlap control count changed",
    )
    rows_by_D, _nonseptimal = base.build_rows()
    literal_lcm = literal_lcm_controls(rows_by_D)
    queries = base.engine_queries(rows_by_D)
    _raw, answers, engine_output = base.run_engine(queries)

    participating_shapes = set()
    deleted_shapes = 0
    deleted_occurrences = 0
    touched_rows = set()
    touched_bodies = set()
    touched_divisors = set()
    killed_semantic = sha256()
    residual_semantic = sha256()
    residual_shapes_by_c = Counter()
    residual_occurrences_by_c = Counter()
    residual_rows = set()
    residual_bodies = set()
    residual_divisors = set()
    residual_rows_by_c = {c: set() for c in range(1, 6)}
    residual_bodies_by_c = {c: set() for c in range(1, 6)}
    residual_divisors_by_c = {c: set() for c in range(1, 6)}
    residual_minimum = None
    residual_minimum_by_c = {c: None for c in range(1, 6)}

    for D in sorted(rows_by_D):
        q = D // 7
        coefficient = fixed_d6_shape_count(D)
        minimum_N = {
            c: min(row[4][c - 1] for row in rows_by_D[D])
            for c in range(1, 6)
        }
        for c in range(1, 6):
            count = answers[(D, c, minimum_N[c])][0]
            if c == 4 and coefficient and 1 <= minimum_N[c] <= q // 6:
                count -= coefficient
                deleted_shapes += coefficient
            require(count >= 0, ("negative residual shape count", D, c, count))
            residual_shapes_by_c[c] += count

        for body, L, support_count, histogram, N_values in rows_by_D[D]:
            N4 = N_values[3]
            decrement = coefficient if coefficient and 1 <= N4 <= q // 6 else 0
            if decrement:
                record = (D, body, L, support_count, histogram, N4, coefficient)
                killed_semantic.update(f"{record}\n".encode())
                participating_shapes.add((D, coefficient))
                deleted_occurrences += coefficient
                touched_rows.add((body, D))
                touched_bodies.add(body)
                touched_divisors.add(D)

            row_alive = False
            for c, N in enumerate(N_values, 1):
                count = answers[(D, c, N)][0]
                if c == 4:
                    count -= decrement
                require(count >= 0, ("negative residual occurrence", D, body, c))
                if not count:
                    continue
                row_alive = True
                record = (
                    D,
                    body,
                    L,
                    support_count,
                    histogram,
                    c,
                    N,
                    count,
                )
                residual_semantic.update(f"{record}\n".encode())
                residual_occurrences_by_c[c] += count
                residual_rows_by_c[c].add((body, D))
                residual_bodies_by_c[c].add(body)
                residual_divisors_by_c[c].add(D)
                if residual_minimum is None or record < residual_minimum:
                    residual_minimum = record
                if (
                    residual_minimum_by_c[c] is None
                    or record < residual_minimum_by_c[c]
                ):
                    residual_minimum_by_c[c] = record
            if row_alive:
                residual_rows.add((body, D))
                residual_bodies.add(body)
                residual_divisors.add(D)

    participating_total = sum(value for _D, value in participating_shapes)
    residual = (
        sum(residual_shapes_by_c.values()),
        sum(residual_occurrences_by_c.values()),
        len(residual_rows),
        len(residual_bodies),
        len(residual_divisors),
    )
    require(participating_total == EXPECTED["participating_shapes"], "shape support changed")
    require(deleted_shapes == EXPECTED["deleted_shapes"], "deleted shape count changed")
    require(
        deleted_occurrences == EXPECTED["deleted_occurrences"],
        "deleted occurrence count changed",
    )
    require(
        (len(touched_rows), len(touched_bodies), len(touched_divisors))
        == EXPECTED["touched"],
        "touched support changed",
    )
    require(
        killed_semantic.hexdigest() == EXPECTED["killed_semantic"],
        "killed semantic digest changed",
    )
    require(residual == EXPECTED["residual"], "residual union changed")
    require(
        dict(residual_shapes_by_c) == EXPECTED["residual_shapes_by_c"],
        "residual shapes-by-c changed",
    )
    require(
        dict(residual_occurrences_by_c) == EXPECTED["residual_occurrences_by_c"],
        "residual occurrences-by-c changed",
    )
    require(
        {c: len(value) for c, value in residual_rows_by_c.items()}
        == EXPECTED["residual_rows_by_c"],
        "residual rows-by-c changed",
    )
    require(
        {c: len(value) for c, value in residual_bodies_by_c.items()}
        == EXPECTED["residual_bodies_by_c"],
        "residual bodies-by-c changed",
    )
    require(
        {c: len(value) for c, value in residual_divisors_by_c.items()}
        == EXPECTED["residual_divisors_by_c"],
        "residual divisors-by-c changed",
    )
    require(residual_minimum == EXPECTED["residual_minimum"], "minimum changed")
    require(
        {c: value[0] for c, value in residual_minimum_by_c.items()}
        == EXPECTED["residual_minimum_D_by_c"],
        "minimum D-by-c changed",
    )
    if EXPECTED["literal_lcm_controls"] is not None:
        require(literal_lcm == EXPECTED["literal_lcm_controls"], "lcm controls changed")
    if EXPECTED["residual_semantic"] is not None:
        require(
            residual_semantic.hexdigest() == EXPECTED["residual_semantic"],
            "residual semantic digest changed",
        )

    print("LRC14 k2 c4 sole-d6 located-phase closure")
    print(f"base_script_sha256={file_sha256(BASE_PATH)}")
    print(f"base_engine_output_sha256={sha256(engine_output.encode()).hexdigest()}")
    print(f"bernoulli_h_tables={EXPECTED_H_TABLES}")
    print("sharp_one_comb_overlap=1/14 equality_only_at_a=2s")
    print("second_distinct_overlap_bound=1/28")
    print("two_comb_low_cover_bound=3/28<1/7")
    print(f"literal_overlap_controls={overlap_count}")
    print(f"literal_lcm_controls={literal_lcm}")
    print(f"participating_d6_shapes={participating_total}")
    print(f"deleted_global_shapes={deleted_shapes}")
    print(f"deleted_occurrences={deleted_occurrences}")
    print(f"touched_rows_bodies_divisors={EXPECTED['touched']}")
    print(f"killed_semantic_sha256={killed_semantic.hexdigest()}")
    print(f"residual_shapes={residual[0]}")
    print(f"residual_occurrences={residual[1]}")
    print(f"residual_rows={residual[2]}")
    print(f"residual_bodies={residual[3]}")
    print(f"residual_divisors={residual[4]}")
    print(f"residual_shapes_by_c={residual_shapes_by_c}")
    print(f"residual_occurrences_by_c={residual_occurrences_by_c}")
    print(
        "residual_rows_by_c="
        f"{Counter({c: len(v) for c, v in residual_rows_by_c.items()})}"
    )
    print(
        "residual_bodies_by_c="
        f"{Counter({c: len(v) for c, v in residual_bodies_by_c.items()})}"
    )
    print(
        "residual_divisors_by_c="
        f"{Counter({c: len(v) for c, v in residual_divisors_by_c.items()})}"
    )
    print(f"residual_minimum={residual_minimum}")
    print(f"residual_minimum_by_c={residual_minimum_by_c}")
    print(f"residual_semantic_sha256={residual_semantic.hexdigest()}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()

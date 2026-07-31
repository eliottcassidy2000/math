#!/usr/bin/env python3
"""Count the k=2,c=4 survivors killed by the located d=6 no-containment lemma.

The lemma applies when the sole vertical denominator is 6 and N_4>0.
The canonical floor/exception screen then necessarily has 1<=N_4<=4.
For fixed D, Mobius inversion counts the four transverse completions around
the fixed d=6 symbol without materializing denominator multisets.
"""

from collections import Counter
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "04-computation/lrc14_k2_septimal_floor_exception_gf_thm2928.py"

spec = spec_from_file_location("k2_septimal", SOURCE)
k2 = module_from_spec(spec)
spec.loader.exec_module(k2)


def fixed_d6_shape_count(D):
    q = D // 7
    if q % 6:
        return 0
    answer = 0
    for E in k2.support.divisors(D):
        if E % 6:
            continue
        sign = k2.mobius(D // E)
        if not sign:
            continue
        transverse = sum(
            q % d != 0
            for d in k2.support.divisors(E)
            if d > 1
        )
        answer += sign * k2.multichoose(transverse, 4)
    if answer < 0:
        raise RuntimeError(("negative fixed-d6 coefficient", D, answer))
    return answer


def main():
    rows_by_D, _nonseptimal = k2.build_rows()
    queries = k2.engine_queries(rows_by_D)
    _raw, answers, _engine_output = k2.run_engine(queries)
    killed_shapes = set()
    deleted_shapes = 0
    killed_occurrences = 0
    killed_rows = set()
    killed_bodies = set()
    killed_divisors = set()
    by_D = Counter()
    by_N = Counter()
    semantic = sha256()
    minimum = None
    residual_rows = set()
    residual_bodies = set()
    residual_divisors = set()
    residual_rows_by_c = {c: set() for c in range(1, 6)}
    residual_bodies_by_c = {c: set() for c in range(1, 6)}
    residual_divisors_by_c = {c: set() for c in range(1, 6)}
    residual_shapes_by_c = Counter()
    residual_occurrences_by_c = Counter()
    residual_minimum = None
    residual_minimum_by_c = {c: None for c in range(1, 6)}
    for D in sorted(rows_by_D):
        coefficient = fixed_d6_shape_count(D)
        q = D // 7
        minimum_N4 = min(row[4][3] for row in rows_by_D[D])
        original_shape_count, _first = answers[(D, 4, minimum_N4)]
        shape_decrement = (
            coefficient if coefficient and 1 <= minimum_N4 <= q // 6 else 0
        )
        deleted_shapes += shape_decrement
        residual_shapes_by_c[4] += original_shape_count - shape_decrement
        for c in (1, 2, 3, 5):
            minimum_N = min(row[4][c - 1] for row in rows_by_D[D])
            residual_shapes_by_c[c] += answers[(D, c, minimum_N)][0]
        if not coefficient:
            coefficient = 0
        for body, L, support_count, histogram, N_values in rows_by_D[D]:
            N4 = N_values[3]
            decrement = coefficient if 1 <= N4 <= q // 6 else 0
            if decrement:
                record = (D, body, L, support_count, histogram, N4, coefficient)
                semantic.update(f"{record}\n".encode())
                killed_shapes.add((D, coefficient))
                killed_occurrences += coefficient
                killed_rows.add((body, D))
                killed_bodies.add(body)
                killed_divisors.add(D)
                by_D[D] += coefficient
                by_N[N4] += coefficient
                if minimum is None or record < minimum:
                    minimum = record
            row_alive = False
            for c, N in enumerate(N_values, 1):
                count = answers[(D, c, N)][0]
                if c == 4:
                    count -= decrement
                if not count:
                    continue
                row_alive = True
                residual_record = (
                    D,
                    body,
                    L,
                    support_count,
                    histogram,
                    c,
                    N,
                    count,
                )
                if residual_minimum is None or residual_record < residual_minimum:
                    residual_minimum = residual_record
                if (
                    residual_minimum_by_c[c] is None
                    or residual_record < residual_minimum_by_c[c]
                ):
                    residual_minimum_by_c[c] = residual_record
                residual_occurrences_by_c[c] += count
                residual_rows_by_c[c].add((body, D))
                residual_bodies_by_c[c].add(body)
                residual_divisors_by_c[c].add(D)
            if row_alive:
                residual_rows.add((body, D))
                residual_bodies.add(body)
                residual_divisors.add(D)
    print("k2 c4 sole-d6 located-phase closure ledger")
    print(f"killed_denominator_shapes={sum(value for _D, value in killed_shapes)}")
    print(f"deleted_global_shapes={deleted_shapes}")
    print(f"killed_occurrences={killed_occurrences}")
    print(f"killed_rows={len(killed_rows)}")
    print(f"killed_bodies={len(killed_bodies)}")
    print(f"killed_divisors={len(killed_divisors)}")
    print(f"by_N4={by_N}")
    print(f"by_D_first20={by_D.most_common(20)}")
    print(f"minimum={minimum}")
    print(f"semantic_sha256={semantic.hexdigest()}")
    print(f"residual_shapes={sum(residual_shapes_by_c.values())}")
    print(f"residual_occurrences={sum(residual_occurrences_by_c.values())}")
    print(f"residual_rows={len(residual_rows)}")
    print(f"residual_bodies={len(residual_bodies)}")
    print(f"residual_divisors={len(residual_divisors)}")
    print(f"residual_shapes_by_c={residual_shapes_by_c}")
    print(f"residual_occurrences_by_c={residual_occurrences_by_c}")
    print(f"residual_minimum={residual_minimum}")
    print(f"residual_minimum_by_c={residual_minimum_by_c}")
    print(
        "residual_rows_by_c=",
        Counter({c: len(value) for c, value in residual_rows_by_c.items()}),
    )
    print(
        "residual_bodies_by_c=",
        Counter({c: len(value) for c, value in residual_bodies_by_c.items()}),
    )
    print(
        "residual_divisors_by_c=",
        Counter({c: len(value) for c, value in residual_divisors_by_c.items()}),
    )


if __name__ == "__main__":
    main()

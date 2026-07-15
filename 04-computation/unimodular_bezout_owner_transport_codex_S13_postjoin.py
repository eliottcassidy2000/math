#!/usr/bin/env python3
"""Exact replay for THM-829's unimodular inverse-owner transport.

The theorem-bearing object is a primitive integer column v=(a,q)^T together
with a Bezout row b=(j,k), b v=1.  A unimodular column action v'=A v has the
contragredient row action b'=b A^{-1}.  The script checks this on a broad
finite matrix/column bank, verifies reflection conjugacy and its tiny GL2(Z)
centralizer, emits q=13 tables for the THM-812/813 convergent matrices, and
replays THM-808's five prime-sheet root translations using Bezout owners.

Tournament Analysis compares candidate arithmetic carriers.  Its pairwise
observable is the number of target-distinct finite-bank record pairs separated;
retention and retention per partition cell are the switches.  These planning
vertices are arithmetic packets, not runners or unmarked tournament classes.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from tournament_tiling_metagraph_address_codex_S4 import carrier_tournament


Matrix = tuple[tuple[int, int], tuple[int, int]]
Column = tuple[int, int]
Row = tuple[int, int]

ROOT_INPUT = Path("05-knowledge/results/continued_fraction_redundancy_root_transport_codex_S13.json")
OUTPUT = Path("05-knowledge/results/unimodular_bezout_owner_transport_codex_S13_postjoin.out")
JSON_OUTPUT = Path("05-knowledge/results/unimodular_bezout_owner_transport_codex_S13_postjoin.json")

IDENTITY: Matrix = ((1, 0), (0, 1))
REFLECTION: Matrix = ((-1, 1), (0, 1))


def mat_mul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )  # type: ignore[return-value]


def mat_vec(matrix: Matrix, vector: Column) -> Column:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def row_mat(row: Row, matrix: Matrix) -> Row:
    return (
        row[0] * matrix[0][0] + row[1] * matrix[1][0],
        row[0] * matrix[0][1] + row[1] * matrix[1][1],
    )


def determinant(matrix: Matrix) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def inverse(matrix: Matrix) -> Matrix:
    det = determinant(matrix)
    assert abs(det) == 1
    return (
        (matrix[1][1] // det, -matrix[0][1] // det),
        (-matrix[1][0] // det, matrix[0][0] // det),
    )


def conjugate_reflection(matrix: Matrix) -> Matrix:
    return mat_mul(mat_mul(REFLECTION, matrix), REFLECTION)


def canonical_bezout(a: int, q: int) -> Row:
    assert q > 1 and 0 < a < q and math.gcd(a, q) == 1
    j = pow(a, -1, q)
    k = (1 - j * a) // q
    assert j * a + k * q == 1
    return j, k


def normalize(vector: Column, row: Row) -> tuple[Column, Row]:
    """Normalize q>1 and 0<a<q while retaining row*column=1."""
    a_raw, q_raw = vector
    j_raw, k_raw = row
    if q_raw < 0:
        a_raw, q_raw = -a_raw, -q_raw
        j_raw, k_raw = -j_raw, -k_raw
    assert q_raw > 1
    a = a_raw % q_raw
    assert 0 < a < q_raw
    quotient = (a_raw - a) // q_raw
    # (a_raw,q) = ((1,quotient),(0,1)) (a,q), so the row is
    # multiplied on the right by that normalization shear.
    j = j_raw
    k = k_raw + quotient * j_raw
    assert j * a + k * q_raw == 1
    j_canonical = j % q_raw
    row_shift = (j - j_canonical) // q_raw
    k_canonical = k + row_shift * a
    assert j_canonical * a + k_canonical * q_raw == 1
    return (a, q_raw), (j_canonical, k_canonical)


def act(matrix: Matrix, a: int, q: int) -> dict[str, object]:
    source_column = (a, q)
    source_row = canonical_bezout(a, q)
    raw_column = mat_vec(matrix, source_column)
    if abs(raw_column[1]) <= 1:
        raise ValueError("action reaches a cusp or unit denominator")
    raw_row = row_mat(source_row, inverse(matrix))
    assert raw_row[0] * raw_column[0] + raw_row[1] * raw_column[1] == 1
    target_column, target_row = normalize(raw_column, raw_row)
    assert target_row == canonical_bezout(*target_column)
    return {
        "source_column": source_column,
        "source_row": source_row,
        "matrix": matrix,
        "raw_target_column": raw_column,
        "raw_target_row": raw_row,
        "target_column": target_column,
        "target_row": target_row,
        "target_inverse_pair": (target_row[0], target_column[1] - target_row[0]),
    }


def digit_matrix(digit: int) -> Matrix:
    assert digit >= 0
    return ((digit, 1), (1, 0))


def continued_fraction_matrix(digits: tuple[int, ...]) -> Matrix:
    matrix = IDENTITY
    for digit in digits:
        matrix = mat_mul(matrix, digit_matrix(digit))
    return matrix


def transport_bank(matrices: dict[str, Matrix]) -> tuple[list[dict], int]:
    records = []
    exact_failures = 0
    for q in range(2, 31):
        for a in range(1, q):
            if math.gcd(a, q) != 1:
                continue
            j, k = canonical_bezout(a, q)
            for name, matrix in matrices.items():
                raw_q = mat_vec(matrix, (a, q))[1]
                if abs(raw_q) <= 1:
                    continue
                result = act(matrix, a, q)
                target_a, target_q = result["target_column"]
                target_j, target_k = result["target_row"]
                exact_failures += target_j * target_a + target_k * target_q != 1
                records.append(
                    {
                        "a": a,
                        "q": q,
                        "j": j,
                        "k": k,
                        "matrix_name": name,
                        "matrix": matrix,
                        "target": (target_a, target_q, target_j, target_k),
                    }
                )
    assert exact_failures == 0
    return records, exact_failures


def carrier_stats(records: list[dict], fields: tuple[str, ...]) -> dict[str, int]:
    groups: dict[tuple, list[int]] = defaultdict(list)
    targets = [record["target"] for record in records]
    for index, record in enumerate(records):
        groups[tuple(record[field] for field in fields)].append(index)
    conflict_pairs = 0
    for indices in groups.values():
        counts = Counter(targets[index] for index in indices)
        total = len(indices)
        conflict_pairs += total * (total - 1) // 2
        conflict_pairs -= sum(size * (size - 1) // 2 for size in counts.values())
    all_counts = Counter(targets)
    total = len(records)
    target_distinct_pairs = total * (total - 1) // 2 - sum(
        size * (size - 1) // 2 for size in all_counts.values()
    )
    return {
        "separated_pairs": target_distinct_pairs - conflict_pairs,
        "cells": len(groups),
        "conflict_pairs": conflict_pairs,
        "target_distinct_pairs": target_distinct_pairs,
    }


def audit(root_input: Path) -> dict:
    a32 = continued_fraction_matrix((1, 2))
    a43 = continued_fraction_matrix((1, 3))
    composition = mat_mul(a43, a32)
    assert a32 == ((3, 1), (2, 1))
    assert a43 == ((4, 1), (3, 1))
    assert composition == ((14, 5), (11, 4))
    matrices = {"A32_[1;2]": a32, "A43_[1;3]": a43, "A43_A32": composition}

    records, exact_failures = transport_bank(matrices)

    # The reflection square uses R A R on the reflected source, never A
    # itself unless A lies in the GL2(Z) centralizer of R.
    reflection_failures = owner_swap_failures = 0
    for q in range(2, 31):
        for a in range(1, q):
            if math.gcd(a, q) != 1:
                continue
            for matrix in matrices.values():
                if abs(mat_vec(matrix, (a, q))[1]) <= 1:
                    continue
                forward = act(matrix, a, q)
                reflected = act(conjugate_reflection(matrix), q - a, q)
                a1, q1 = forward["target_column"]
                ar, qr = reflected["target_column"]
                j1, _ = forward["target_row"]
                jr, _ = reflected["target_row"]
                reflection_failures += (ar, qr) != (q1 - a1, q1)
                owner_swap_failures += jr != q1 - j1
    assert reflection_failures == owner_swap_failures == 0

    centralizer = []
    for alpha in range(-5, 6):
        for beta in range(-5, 6):
            for gamma in range(-5, 6):
                for delta in range(-5, 6):
                    matrix = ((alpha, beta), (gamma, delta))
                    if abs(determinant(matrix)) != 1:
                        continue
                    if mat_mul(matrix, REFLECTION) == mat_mul(REFLECTION, matrix):
                        centralizer.append(matrix)
    expected_centralizer = {IDENTITY, ((-1, 0), (0, -1)), REFLECTION, ((1, -1), (0, -1))}
    assert set(centralizer) == expected_centralizer

    q13_table = []
    for a in range(1, 13):
        row = {"a": a, "old_inverse_pair": (pow(a, -1, 13), 13 - pow(a, -1, 13))}
        for name, matrix in matrices.items():
            result = act(matrix, a, 13)
            row[name] = {
                "target_column": result["target_column"],
                "target_inverse_pair": result["target_inverse_pair"],
            }
        q13_table.append(row)

    # THM-808 fixed-prime-sheet compatibility.  An integer lift shear
    # T_m:(a,q)->(a+m q,q) fixes the inverse owner j.  Hence each published
    # root translation is recovered from the Bezout first coordinates.
    root_data = json.loads(root_input.read_text())
    prime = root_data["prime"]
    speeds = root_data["speeds"]
    root_blocks = []
    root_translation_failures = shear_owner_failures = 0
    for block in root_data["block_types"]:
        counts = block["owner_counts"]
        owners = [canonical_bezout(speed % prime, prime)[0] for speed in speeds]
        translation = -sum(count * owner for count, owner in zip(counts, owners)) % prime
        expected = block["root_translation_mod_prime"][0]
        root_translation_failures += translation != expected
        for speed in speeds:
            a = speed % prime
            for lift in range(-3, 4):
                shear = ((1, lift), (0, 1))
                transformed = act(shear, a, prime)
                shear_owner_failures += transformed["target_row"][0] != pow(a, -1, prime)
        root_blocks.append(
            {
                "rank": block["rank"],
                "occurrences": block["occurrences"],
                "translation": translation,
                "expected": expected,
            }
        )
    assert root_translation_failures == shear_owner_failures == 0

    candidates = {
        "inverse_only": ("j",),
        "inverse_matrix": ("j", "matrix_name"),
        "bezout_matrix": ("j", "k", "matrix_name"),
        "column_matrix": ("a", "q", "matrix_name"),
        "full_stalk": ("a", "q", "j", "k", "matrix_name"),
    }
    stats = {name: carrier_stats(records, fields) for name, fields in candidates.items()}
    assert stats["inverse_matrix"]["conflict_pairs"] > 0
    assert stats["bezout_matrix"]["conflict_pairs"] > 0
    assert stats["column_matrix"]["conflict_pairs"] == 0
    assert stats["full_stalk"]["conflict_pairs"] == 0
    tournament_input = {
        name: {"separated_pairs": row["separated_pairs"], "cells": row["cells"]}
        for name, row in stats.items()
    }
    retention = carrier_tournament(tournament_input, "retention")
    economy = carrier_tournament(tournament_input, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(candidates))
        for j in range(i + 1, len(candidates))
    )

    return {
        "schema_version": 1,
        "theorem": "THM-829",
        "matrices": {name: matrix for name, matrix in matrices.items()},
        "reflected_matrices": {
            name: conjugate_reflection(matrix) for name, matrix in matrices.items()
        },
        "finite_transport_records": len(records),
        "exact_failures": exact_failures,
        "reflection_failures": reflection_failures,
        "inverse_owner_swap_failures": owner_swap_failures,
        "centralizer_GL2Z_of_reflection": centralizer,
        "q13_table": q13_table,
        "thm808_root_blocks": root_blocks,
        "thm808_root_translation_failures": root_translation_failures,
        "fixed_denominator_shear_owner_failures": shear_owner_failures,
        "carrier_partitions": stats,
        "tournament_analysis": {
            "vertices": list(candidates),
            "pairwise_observable": "number of target-distinct finite-bank transport pairs separated",
            "switches": ["separation retention", "separation retention per partition cell"],
            "tie_hamiltonian_path": list(candidates),
            "retention": retention,
            "economy": economy,
            "edge_flips": flips,
        },
    }


def matrix_text(matrix: Matrix) -> str:
    return f"(({matrix[0][0]},{matrix[0][1]}),({matrix[1][0]},{matrix[1][1]}))"


def render(result: dict) -> str:
    lines = [
        "THM-829 UNIMODULAR BEZOUT / INVERSE-OWNER TRANSPORT",
        "=" * 76,
        "law: v'=A v, b'=b A^(-1), with b v=b' v'=1",
        f"finite records={result['finite_transport_records']} exact failures={result['exact_failures']}",
        f"reflection/owner-swap failures={result['reflection_failures']}/{result['inverse_owner_swap_failures']}",
        "",
        "CENTERED-RATIO CONVERGENT MATRICES",
    ]
    for name, matrix in result["matrices"].items():
        lines.append(
            f"  {name}: A={matrix_text(matrix)} RAR={matrix_text(result['reflected_matrices'][name])}"
        )
    lines.extend(
        [
            "",
            f"GL2(Z) centralizer of R in search box={tuple(matrix_text(m) for m in result['centralizer_GL2Z_of_reflection'])}",
            "proved centralizer={+I,-I,+R,-R}",
            "",
            "q=13 PRIMITIVE WITNESS TABLE",
            "  a old(j+,j-) | A32:(a',q';j+',j-') | A43 | A43*A32",
        ]
    )
    for row in result["q13_table"]:
        cells = []
        for name in result["matrices"]:
            cell = row[name]
            cells.append(f"{tuple(cell['target_column'])};{tuple(cell['target_inverse_pair'])}")
        lines.append(f"  {row['a']:2d} {tuple(row['old_inverse_pair'])} | " + " | ".join(cells))
    lines.extend(["", "THM-808 ROOT BLOCKS VIA BEZOUT OWNERS"])
    for block in result["thm808_root_blocks"]:
        lines.append(
            f"  rank={block['rank']} occurrences={tuple(block['occurrences'])} translation={block['translation']} expected={block['expected']}"
        )
    lines.extend(
        [
            f"  translation/shear-owner failures={result['thm808_root_translation_failures']}/{result['fixed_denominator_shear_owner_failures']}",
            "",
            "CARRIER PURITY",
        ]
    )
    for name, row in result["carrier_partitions"].items():
        lines.append(
            f"  {name}: cells={row['cells']} separated={row['separated_pairs']} conflict_pairs={row['conflict_pairs']}"
        )
    ta = result["tournament_analysis"]
    lines.extend(
        [
            "",
            "TOURNAMENT ANALYSIS (arithmetic carriers as planning vertices)",
            f"  vertices={tuple(ta['vertices'])}",
            f"  observable={ta['pairwise_observable']}",
            f"  switches={tuple(ta['switches'])}; flips={ta['edge_flips']}",
            "",
            "Verdict: inverse owners form a contragredient Bezout-row stalk.",
            "Reflection requires the conjugate action RAR except for +/-I,+/-R.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root-input", type=Path, default=ROOT_INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--json", type=Path, default=JSON_OUTPUT)
    args = parser.parse_args()
    result = audit(args.root_input)
    text = render(result)
    print(text, end="")
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()

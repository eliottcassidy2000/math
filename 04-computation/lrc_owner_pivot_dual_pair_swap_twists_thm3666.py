#!/usr/bin/env python3
"""Exact owner-pivot dual-basis and pair-swap certificate for THM-3666."""

from __future__ import annotations

import ast
from hashlib import sha256
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENTS = (
    ROOT / "04-computation/lrc14_owner_aligned_pivot_thm2309.py",
    ROOT / "05-knowledge/results/lrc14_owner_aligned_pivot_thm2309.out",
    ROOT / "04-computation/lrc14_relation_residue_pushforward_thm2334.py",
    ROOT / "05-knowledge/results/lrc14_relation_residue_pushforward_thm2334.out",
    ROOT / "04-computation/lrc_minimal_three_twist_target_detector_thm3665.py",
    ROOT / "05-knowledge/results/lrc_minimal_three_twist_target_detector_thm3665.out",
)
EXPECTED_PARENT_HASHES = (
    "a17a015f9187b59a36cef541100ae29a83395690294148ddee6ec5a45e0ea889",
    "f487850adfe0e0c5cc1f71efabf85f5e8c150bc42453b08f6d87afc41a322b4b",
    "ce220ede12175b6851810782c880f95048fe9f4643cc4f52f47a7f4d8dcb7b0c",
    "d2d9b49db9ef3eabf7e3ae17cea247da554a9f8df2abfc3907243d317b21fec1",
    "a5ef75f038d80c5d91308bb5379303970b44e9e538323eb49cf8779386356938",
    "172fe3e32fc27bb2abb21f4c7a7af59e71cdfa4c604586b2d4f8725a5ac6211a",
)
EXPECTED_SEMANTIC_SHA256 = "b40649847f46b03b54cc48cb0c6968c861936e4c81441eaaf33adc2f6636ffa0"

P = 13
N = P * P
UNITS = tuple(range(6))
SOURCE = 6
TARGET_A = 7
TARGET_B = 8
LABEL_COUNT = 9


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(left, right)) % P


def owner_rows(
    speeds: tuple[int, ...], omitted: int, graft_a: int, graft_b: int
) -> tuple[tuple[int, ...], ...]:
    pivot_labels = (SOURCE,) + tuple(unit for unit in UNITS if unit != omitted)
    rows = []
    for label in pivot_labels:
        row = [0] * LABEL_COUNT
        row[omitted] = speeds[label] % P
        row[label] = (row[label] - speeds[omitted]) % P
        if label == graft_a:
            row[omitted] = (row[omitted] + speeds[TARGET_A]) % P
            row[TARGET_A] = (row[TARGET_A] - speeds[omitted]) % P
        if label == graft_b:
            row[omitted] = (row[omitted] + speeds[TARGET_B]) % P
            row[TARGET_B] = (row[TARGET_B] - speeds[omitted]) % P
        rows.append(tuple(row))
    return tuple(rows)


def dual_basis(graft_a: int, graft_b: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    alpha = [0] * LABEL_COUNT
    beta = [0] * LABEL_COUNT
    alpha[TARGET_A] = 1
    alpha[graft_a] = -1 % P
    beta[TARGET_B] = 1
    beta[graft_b] = -1 % P
    return tuple(alpha), tuple(beta)


def main() -> None:
    parent_hashes = tuple(lf_sha256(path) for path in PARENTS)
    require(parent_hashes == EXPECTED_PARENT_HASHES,
            ("parent hashes", parent_hashes))

    omitted = 0
    graft_a = 1
    graft_b = 2
    alpha, beta = dual_basis(graft_a, graft_b)
    require(dot(alpha, beta) == 0, "disjoint pair-swap support")

    # Scale the omitted unit to one.  This exhausts all 12^5 projective unit
    # scalar types for the fixed label choice; the formulas are label-
    # symmetric for the other 119 owner-packet choices.
    scalar_digest = sha256()
    scalar_count = 0
    for tail in itertools.product(range(1, P), repeat=5):
        speeds = (1,) + tail + (0, 0, 0)
        rows = owner_rows(speeds, omitted, graft_a, graft_b)
        pivot_labels = (SOURCE,) + tuple(unit for unit in UNITS if unit != omitted)
        require(all(dot(row, speeds) == 0 for row in rows),
                ("relation row", speeds, rows))
        require(all(rows[index][label] == -speeds[omitted] % P
                    for index, label in enumerate(pivot_labels)),
                ("pivot diagonal", speeds, rows))
        require(all(rows[index][other] == 0
                    for index, _label in enumerate(pivot_labels)
                    for other in pivot_labels
                    if other != pivot_labels[index]),
                ("pivot off diagonal", speeds, rows))
        require(all(dot(row, alpha) == 0 and dot(row, beta) == 0
                    for row in rows),
                ("dual basis orthogonality", speeds, rows))
        require(dot(alpha, tuple(1 if index == TARGET_A else 0
                                 for index in range(LABEL_COUNT))) == 1
                and dot(alpha, tuple(1 if index == TARGET_B else 0
                                     for index in range(LABEL_COUNT))) == 0
                and dot(beta, tuple(1 if index == TARGET_A else 0
                                    for index in range(LABEL_COUNT))) == 0
                and dot(beta, tuple(1 if index == TARGET_B else 0
                                    for index in range(LABEL_COUNT))) == 1,
                "target dual pairing")
        scalar_digest.update(json.dumps((speeds, rows), separators=(",", ":")).encode("ascii"))
        scalar_digest.update(b"\n")
        scalar_count += 1
    require(scalar_count == 12 ** 5, scalar_count)
    scalar_sha256 = scalar_digest.hexdigest()

    # Check all six omitted-unit and twenty ordered graft choices on a hostile
    # nonconstant unit residue profile.  Rank is supplied by the explicit
    # diagonal pivot, so orthogonality plus dimension gives
    # L^perp=<w,alpha,beta> in every case.
    control_speeds = (1, 2, 3, 4, 5, 6, 0, 0, 0)
    choice_ledger = []
    for omitted in UNITS:
        available = tuple(unit for unit in UNITS if unit != omitted)
        for graft_a in available:
            for graft_b in available:
                if graft_b == graft_a:
                    continue
                rows = owner_rows(control_speeds, omitted, graft_a, graft_b)
                alpha, beta = dual_basis(graft_a, graft_b)
                require(all(dot(row, control_speeds) == 0
                            and dot(row, alpha) == 0
                            and dot(row, beta) == 0 for row in rows),
                        ("choice orthogonality", omitted, graft_a, graft_b))
                pivot_labels = (SOURCE,) + tuple(
                    unit for unit in UNITS if unit != omitted
                )
                require(all(rows[index][label] == -control_speeds[omitted] % P
                            for index, label in enumerate(pivot_labels)),
                        ("choice pivot", omitted, graft_a, graft_b))
                choice_ledger.append((omitted, graft_a, graft_b, alpha, beta))
    require(len(choice_ledger) == 6 * 5 * 4, len(choice_ledger))
    choice_digest = digest_json(tuple(choice_ledger))

    # Pair-swap coordinate action and word neutrality.  A twist x*alpha+y*beta
    # is supported on exactly the four advertised coordinates (unless x or y
    # vanishes), and R*ell/13 is integral for every LRC depth exponent.
    alpha, beta = dual_basis(1, 2)
    twist_ledger = []
    for x in range(P):
        for y in range(P):
            twist = tuple((x * alpha[index] + y * beta[index]) % P
                          for index in range(LABEL_COUNT))
            require(twist[TARGET_A] == x
                    and twist[1] == -x % P
                    and twist[TARGET_B] == y
                    and twist[2] == -y % P,
                    ("pair-swap coordinates", x, y, twist))
            require(all(twist[index] == 0
                        for index in (0, 3, 4, 5, SOURCE)),
                    ("pair-swap leakage", x, y, twist))
            for exponent in range(1, 9):
                dilation = P ** exponent
                require(all((dilation * coordinate) % P == 0
                            for coordinate in twist),
                        ("word neutrality", exponent, x, y, twist))
            twist_ledger.append(((x, y), twist))
    require(len(twist_ledger) == N, len(twist_ledger))
    twist_digest = digest_json(tuple(twist_ledger))

    # If the deepest coordinate is a unit, the graft labels can always avoid
    # it; source and all ungrafted labels then have zero twist coordinate.
    phase_avoidance = []
    for deepest in UNITS:
        choices_list = []
        for packet in choice_ledger:
            packet_omitted, packet_graft_a, packet_graft_b, _alpha, _beta = packet
            if deepest != packet_graft_a and deepest != packet_graft_b:
                choices_list.append(
                    (packet_omitted, packet_graft_a, packet_graft_b)
                )
        choices = tuple(choices_list)
        require(len(choices) > 0, ("deep phase avoidance", deepest))
        phase_avoidance.append((deepest, len(choices)))
    require(tuple(count for _deepest, count in phase_avoidance) == (80,) * 6,
            phase_avoidance)

    semantic = digest_json((
        parent_hashes,
        P, UNITS, SOURCE, TARGET_A, TARGET_B,
        scalar_count, scalar_sha256,
        tuple(choice_ledger), choice_digest,
        tuple(twist_ledger), twist_digest,
        tuple(phase_avoidance),
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3666 LRC owner-pivot dual pair-swap twists ==")
    print(f"parents_sha256_lf={parent_hashes}")
    print(f"normalized_scalar_types={scalar_count};sha256={scalar_sha256}")
    print("dual_basis=alpha=e_a-e_ka,beta=e_b-e_kb;L_perp=<w,alpha,beta>")
    print("target_pairing=<alpha,e_a>=<beta,e_b>=1;cross_pairings=0")
    print(f"owner_packet_choices={len(choice_ledger)};sha256={choice_digest}")
    print(f"twist_table=169;support_labels=(a,ka,b,kb);sha256={twist_digest}")
    print("word_action=neutral because 13 divides R;present_action=opposite 1/13 pair shifts")
    print(f"unit_deep_phase_avoidance_counts={tuple(phase_avoidance)}")
    print("three_twist_test=H(x,y)+H(x-1,y)-2H(x,y-1)")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("scope=typed relation-twist coordinate action;not ancestry-digit alignment or nonvanishing")
    print("PASS")


if __name__ == "__main__":
    main()

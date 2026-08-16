#!/usr/bin/env python3
"""Finite-exact K4/V4/XOR model of the nine THM-3479 relation slots."""

from __future__ import annotations

import ast
from functools import reduce
from hashlib import sha256
from itertools import permutations
from json import dumps
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PINS = (
    (
        "THM-3479",
        ROOT / "01-canon/theorems/THM-3479-literal-half-twist-relation-current-two-transplant-certificate.md",
        "025998551e3cdf3c6e4db5c0a4f208dd32f6845970fd4729d4a276035e0fdfeb",
        "PROVED-STRUCTURAL-AUDITED",
    ),
    (
        "ROLE-CHART-SIDECAR",
        ROOT / "04-computation/lrc_relation_role_chart_weighted_closure_probe_20260815.py",
        "207c65ca235ea5647e346027d424264e8abbcf27c5f574b5901cca13611d7e03",
        "FINITE-EXACT",
    ),
)

EXPECTED_SEMANTIC_SHA256 = (
    "904764346382212e617f2bc15bb95e235af42c0cf1e79131842ced17e39538c8"
)

LABELS = ("c1", "c2", "c3", "H", "q1", "q2", "q3", "q4", "q5")
RELATION = dict(zip(
    LABELS,
    (-27, -27, -27, 20110798, -41, -27, -27, -27, 38),
))

TUPLES = {
    "U_full": dict(zip(
        LABELS,
        (13, 2197, 742586, 1, 183, 27, 131, 53, 313),
    )),
    "U_clock": dict(zip(
        LABELS,
        (65, 2197, 742586, 5, 661549, 655231, 658533,
         661445, 291),
    )),
    "U_q27": dict(zip(
        LABELS,
        (28405, 7599423, 18279868269, 3459, 2016, 2757,
         1041, 3693, 11163142875),
    )),
    "U_q51": dict(zip(
        LABELS,
        (70993, 7199569, 30105550319, 5825, 4200, 5214,
         7684, 4421, 18313194875),
    )),
}

DELTAS = (1, 2, 3)
C_LABELS = ("c1", "c2", "c3")
Q_LABELS = ("q2", "q3", "q4")
AXIS_LABELS = ("H", "q1", "q5")


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def source_hashes() -> tuple[tuple[str, str, str], ...]:
    rows = tuple(
        (label, status, lf_sha256(path))
        for label, path, _expected, status in SOURCE_PINS
    )
    expected = tuple(
        (label, status, digest)
        for label, _path, digest, status in SOURCE_PINS
    )
    require(rows == expected, ("source hash drift", rows, expected))
    return rows


def security_certificate() -> tuple[object, ...]:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert forbidden")
    forbidden = {"eval", "exec", "compile", "__import__"}
    calls = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    require(not (calls & forbidden), ("dynamic execution", calls & forbidden))
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add(node.module or "")
    allowed = {
        "__future__", "ast", "functools", "hashlib", "itertools", "json",
        "math", "pathlib",
    }
    require(imports <= allowed, ("unexpected imports", imports - allowed))
    return ("NO_ASSERT_AST", tuple(sorted(imports)), tuple(sorted(forbidden)))


def edge(left: int, right: int) -> frozenset[int]:
    require(left != right, (left, right))
    return frozenset((left, right))


MATCHINGS = {
    delta: frozenset(
        edge(vertex, vertex ^ delta) for vertex in range(4)
    )
    for delta in DELTAS
}

# Choose the c-triple as the star at 00 and the q-triple as its complementary
# triangle.  Each (c_i,q_(i+1)) pair is one XOR-difference matching.
LABEL_TO_EDGE = {
    "c1": edge(0, 1),
    "q2": edge(2, 3),
    "c2": edge(0, 2),
    "q3": edge(1, 3),
    "c3": edge(0, 3),
    "q4": edge(1, 2),
}
EDGE_TO_LABEL = {value: label for label, value in LABEL_TO_EDGE.items()}
DELTA_TO_AXIS = dict(zip(DELTAS, AXIS_LABELS))
AXIS_TO_DELTA = {label: delta for delta, label in DELTA_TO_AXIS.items()}


def matching_image(
    permutation: tuple[int, ...], delta: int
) -> int:
    image = frozenset(
        edge(permutation[min(pair)], permutation[max(pair)])
        for pair in MATCHINGS[delta]
    )
    return next(candidate for candidate in DELTAS
                if MATCHINGS[candidate] == image)


def induced_label_permutation(
    permutation: tuple[int, ...]
) -> dict[str, str]:
    image: dict[str, str] = {}
    for label, pair in LABEL_TO_EDGE.items():
        left, right = tuple(pair)
        image[label] = EDGE_TO_LABEL[edge(permutation[left], permutation[right])]
    for label in AXIS_LABELS:
        image[label] = DELTA_TO_AXIS[
            matching_image(permutation, AXIS_TO_DELTA[label])
        ]
    return image


def relation_stabilizer_certificate() -> tuple[object, ...]:
    group = tuple(tuple(row) for row in permutations(range(4)))
    stabilizer = []
    axis_actions = set()
    for permutation in group:
        label_image = induced_label_permutation(permutation)
        axis_actions.add(tuple(label_image[label] for label in AXIS_LABELS))
        if all(RELATION[label] == RELATION[label_image[label]]
               for label in LABELS):
            stabilizer.append(permutation)
    stabilizer = tuple(stabilizer)
    translations = tuple(
        tuple(vertex ^ translation for vertex in range(4))
        for translation in range(4)
    )
    require(set(stabilizer) == set(translations),
            ("stabilizer is not V4", stabilizer, translations))
    require(len(axis_actions) == 6, axis_actions)
    return (len(group), len(stabilizer), len(axis_actions), stabilizer)


def character_table() -> tuple[tuple[int, ...], ...]:
    rows = []
    for delta in DELTAS:
        signs = []
        c_edge = LABEL_TO_EDGE[C_LABELS[delta - 1]]
        for translation in range(4):
            moved = edge(*(vertex ^ translation for vertex in c_edge))
            signs.append(1 if moved == c_edge else -1)
        rows.append(tuple(signs))
    result = tuple(rows)
    require(result == (
        (1, 1, -1, -1),
        (1, -1, 1, -1),
        (1, -1, -1, 1),
    ), result)
    return result


def selected_type(bits: tuple[int, int, int]) -> str:
    return "STAR" if sum(bits) % 2 == 1 else "TRIANGLE"


def canonical_graph_type(
    values: dict[str, int]
) -> tuple[tuple[int, ...], tuple[int, ...], str]:
    differences = tuple(
        values[c_label] - values[q_label]
        for c_label, q_label in zip(C_LABELS, Q_LABELS)
    )
    require(all(differences), ("tie in canonical matching", differences))
    bits = tuple(int(value > 0) for value in differences)
    selected_labels = tuple(
        c_label if bit else q_label
        for bit, c_label, q_label in zip(bits, C_LABELS, Q_LABELS)
    )
    degrees = [0] * 4
    for label in selected_labels:
        for vertex in LABEL_TO_EDGE[label]:
            degrees[vertex] += 1
    degree_profile = tuple(sorted(degrees))
    kind = selected_type(bits)
    require(
        degree_profile == ((1, 1, 1, 3) if kind == "STAR"
                           else (0, 2, 2, 2)),
        (bits, kind, degree_profile),
    )
    return differences, bits, kind


def pairing_census(name: str, values: dict[str, int]) -> tuple[object, ...]:
    rows = []
    for q_order in permutations(Q_LABELS):
        differences = tuple(
            values[c_label] - values[q_label]
            for c_label, q_label in zip(C_LABELS, q_order)
        )
        require(all(differences), (name, q_order, differences))
        product = reduce(lambda left, right: left * right, differences, 1)
        kind = "STAR" if product > 0 else "TRIANGLE"
        rows.append((q_order, differences, product, kind))
    kinds = tuple(sorted({row[3] for row in rows}))
    canonical = canonical_graph_type(values)
    require(kinds == (canonical[2],), (name, kinds, canonical))
    products = tuple(row[2] for row in rows)
    residues = tuple(value % 13 for value in canonical[0])
    require(len(set(residues)) == 1 and residues[0] != 0,
            (name, "non-diagonal mod13 current", residues))
    return (
        name,
        canonical,
        len(rows),
        kinds,
        min(abs(value) for value in products),
        max(abs(value) for value in products),
        reduce(gcd, (abs(value) for value in products)),
        residues,
        tuple(sorted(values[label] for label in C_LABELS)),
        tuple(sorted(values[label] for label in Q_LABELS)),
    )


def hostile_census() -> tuple[object, ...]:
    values = {
        "c1": 0, "c2": 4, "c3": 8,
        "q2": 2, "q3": 6, "q4": 10,
    }
    kinds = []
    products = []
    for q_order in permutations(Q_LABELS):
        product = reduce(
            lambda left, right: left * right,
            (values[c_label] - values[q_label]
             for c_label, q_label in zip(C_LABELS, q_order)),
            1,
        )
        products.append(product)
        kinds.append("STAR" if product > 0 else "TRIANGLE")
    require(set(kinds) == {"STAR", "TRIANGLE"}, (kinds, products))
    return tuple(kinds), tuple(products)


def main() -> None:
    hashes = source_hashes()
    security = security_certificate()

    require(tuple(label for label in LABELS if RELATION[label] == -27)
            == ("c1", "c2", "c3", "q2", "q3", "q4"),
            "six-edge coefficient packet drift")
    require(tuple(len(MATCHINGS[delta]) for delta in DELTAS) == (2, 2, 2),
            MATCHINGS)
    require(set().union(*MATCHINGS.values()) == set(LABEL_TO_EDGE.values()),
            "matching partition misses an edge")

    stabilizer = relation_stabilizer_certificate()
    characters = character_table()
    censuses = tuple(
        pairing_census(name, values) for name, values in TUPLES.items()
    )
    require(tuple(census[1][2] for census in censuses)
            == ("TRIANGLE", "STAR", "STAR", "STAR"), censuses)
    hostile = hostile_census()

    semantic_payload = {
        "source_hashes": hashes,
        "security": security,
        "relation": tuple(RELATION[label] for label in LABELS),
        "edge_model": tuple(
            (label, tuple(sorted(pair)))
            for label, pair in sorted(LABEL_TO_EDGE.items())
        ),
        "matching_axes": tuple(zip(DELTAS, AXIS_LABELS)),
        "stabilizer": stabilizer,
        "characters": characters,
        "censuses": censuses,
        "hostile": hostile,
        "representation_multiplicities": ("trivial^6", "chi1", "chi2", "chi3"),
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":"),
              default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
                (semantic_hash, EXPECTED_SEMANTIC_SHA256))

    print("LRC NINE-SLOT K4/V4 XOR STAR-TRIANGLE FINITE-EXACT SIDECAR")
    print("STATUS: FINITE-EXACT STRUCTURAL COMPONENT OF PROVED THM-3479; INDEPENDENTLY AUDITED")
    print(f"SOURCE_HASHES: {hashes}")
    print(f"SECURITY: {security}")
    print("K4_EDGE_MODEL: c1/q2, c2/q3, c3/q4 are the three opposite-edge pairs; c-triple is a star and q234-triple its complementary triangle")
    print("AXIS_MODEL: H,q1,q5 label the three perfect-matching/XOR-difference classes")
    print(f"S4_RELATION_STABILIZER: {stabilizer[:3]}; stabilizer is the translation V4")
    print(f"NONTRIVIAL_V4_CHARACTER_TABLE: {characters}")
    print("REPRESENTATION_SPLIT: six edge slots = three invariant pair sums plus one copy of each nontrivial V4 character; H/q1/q5 add three invariant axes")
    for census in censuses:
        print(f"PAIRING_CENSUS: {census}")
    print(f"INTERLACING_HOSTILE: {hostile}")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print("VERDICT: U_full is pairing-robust TRIANGLE type while U_clock and both CRT lifts are pairing-robust STAR type; this is an exact K4/XOR selector invariant correlated with, but not proved to cause, the two-transplant split and it supplies no tournament orientation or physical endpoint current")


if __name__ == "__main__":
    main()

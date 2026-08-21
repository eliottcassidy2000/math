#!/usr/bin/env python3
"""Exact address atlas for THM-3593's rank-two LRC correction quotient."""

from __future__ import annotations

from hashlib import sha256
import ast
import importlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PARENT_PATH = ROOT / "04-computation/lrc_common_a4_anova_graph_flag_thm3593.py"
EXPECTED_PARENT_SHA256 = "3ee8486833d539599a4c8add304172b72929c69f3859e9cc1ce22e3018199516"
EXPECTED_SEMANTIC_SHA256 = "792059590ad928cdb096763c6fd429b8f02b991593efef68964e26c57f7b7059"

PRIME = 755373809845391722745761
P = 13
B = 371578917865089240854253
B_INVERSE = 555904782330327598552358
REVERSAL_ACTION = ((0, B), (B_INVERSE, 0))
EXCEPTIONAL_ORBITS = (
    ((12, 0), (0, 12)),
    ((12, 1), (0, 11)),
    ((6, 5), (6, 7)),
    ((3, 9), (9, 3)),
)


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


require(lf_sha256(PARENT_PATH) == EXPECTED_PARENT_SHA256, "THM-3593 parent hash")
sys.path.insert(0, str(PARENT_PATH.parent))
G = importlib.import_module(PARENT_PATH.stem)
M = G.M
require((G.PRIME, G.P) == (PRIME, P), "parent field drift")


def add_rows(*rows):
    require(rows and all(len(row) == len(rows[0]) for row in rows), "row widths")
    return tuple(sum(values) % PRIME for values in zip(*rows))


def state_side(row):
    constant, state, _relation, _interaction = G.anova_components(row)
    return add_rows(constant, state)


def coordinates_in_rref(basis, row):
    pivots = tuple(next(index for index, value in enumerate(vector) if value)
                   for vector in basis)
    coordinates = tuple(row[index] % PRIME for index in pivots)
    rebuilt = tuple(
        sum(coordinates[index] * basis[index][column]
            for index in range(len(basis))) % PRIME
        for column in range(len(row))
    )
    require(rebuilt == tuple(value % PRIME for value in row),
            "row outside correction plane")
    return coordinates


def reversal(label):
    return tuple((-value - 1) % P for value in label)


def act(vector, matrix=REVERSAL_ACTION):
    return tuple(
        sum(vector[row] * matrix[row][column] for row in range(2)) % PRIME
        for column in range(2)
    )


def multiply(left, right):
    return tuple(tuple(
        sum(left[row][middle] * right[middle][column] for middle in range(2))
        % PRIME
        for column in range(2)
    ) for row in range(2))


def slope(vector):
    if vector == (0, 0):
        return None
    require(vector[0] != 0, ("unexpected vertical quotient line", vector))
    return vector[1] * pow(vector[0], -1, PRIME) % PRIME


def main() -> None:
    source_tensor, source_reconstruction = M.reconstruct_source_current()
    two_tensor, two_reconstruction = M.reconstruct_two_current()
    require(source_reconstruction == M.EXPECTED_PARENT_RECONSTRUCTION_SHA256[2:],
            "source reconstruction drift")
    require(two_reconstruction == M.EXPECTED_PARENT_RECONSTRUCTION_SHA256[:2],
            "two-current reconstruction drift")

    source_raw = M.flatten_source_current(source_tensor)
    two_raw = M.flatten_two_current(two_tensor)
    source_e = tuple(state_side(row) for row in source_raw)
    two_e = tuple(state_side(row) for row in two_raw)
    basis = M.canonical_row_basis(source_e)
    require(len(basis) == 2, "source correction rank")
    require(M.canonical_row_basis(two_e) == basis, "two-current correction plane")
    basis_digest = M.rowspace_digest(basis)
    require(basis_digest == G.EXPECTED_DIGESTS["state_side"],
            ("correction basis digest", basis_digest))

    source_coordinates = tuple(coordinates_in_rref(basis, row) for row in source_e)
    two_coordinates = tuple(coordinates_in_rref(basis, row) for row in two_e)
    source_labels = tuple(M.SOURCE.CHARACTER_PAIRS)
    two_labels = tuple((r0, r1) for r0 in range(P) for r1 in range(P))
    require((len(source_labels), len(two_labels)) == (17, 169), "address counts")
    source_by_label = dict(zip(source_labels, source_coordinates))
    two_by_label = dict(zip(two_labels, two_coordinates))

    source_slopes = tuple(slope(source_by_label[label]) for label in source_labels)
    two_slopes = {label: slope(two_by_label[label]) for label in two_labels}
    require(None not in source_slopes and len(set(source_slopes)) == 17,
            "source lines are not 17 distinct nonzero lines")

    zero_labels = tuple(
        label for label in two_labels
        if (label[1] == 6
            or (label[1] == 0 and label[0] != 12)
            or (label[1] == 12 and label[0] != 0))
    )
    require(tuple(label for label in two_labels if two_slopes[label] is None)
            == zero_labels, "three-row kernel law")
    require(len(zero_labels) == 37, "kernel address count")

    exceptional_labels = frozenset(
        label for orbit in EXCEPTIONAL_ORBITS for label in orbit)
    generic_labels = tuple(
        label for label in two_labels
        if label not in exceptional_labels and label not in zero_labels
    )
    require(len(generic_labels) == 124, "generic address count")
    require(set(two_slopes[label] for label in generic_labels) == {B},
            "generic plus-line law")
    exceptional_slopes = tuple(two_slopes[label]
                               for orbit in EXCEPTIONAL_ORBITS for label in orbit)
    require(None not in exceptional_slopes and len(set(exceptional_slopes)) == 8,
            "exceptional line multiplicity")
    require(B not in exceptional_slopes, "exceptional label on generic line")
    require(len(set(two_slopes.values()) - {None}) == 9,
            "two-current projective class count")

    identity = ((1, 0), (0, 1))
    require(B * B_INVERSE % PRIME == 1, "reversal scalar inverse")
    require(multiply(REVERSAL_ACTION, REVERSAL_ACTION) == identity,
            "reversal action is not involutive")
    require(all(act(two_by_label[label]) == two_by_label[reversal(label)]
                for label in two_labels), "exact reversal covariance")
    require(all(reversal(left) == right and reversal(right) == left
                for left, right in EXCEPTIONAL_ORBITS),
            "exceptional reversal orbits")
    require(set(map(reversal, zero_labels)) == set(zero_labels),
            "kernel set reversal")
    require(set(map(reversal, generic_labels)) == set(generic_labels),
            "generic set reversal")

    plus_line = B
    minus_line = (-B) % PRIME
    require(act((1, plus_line)) == (1, plus_line), "plus eigenline")
    require(act((1, minus_line))
            == ((-1) % PRIME, B), "minus eigenline")
    require(minus_line not in set(two_slopes.values()),
            "hostile minus line unexpectedly realized")

    two_nonzero_lines = set(two_slopes.values()) - {None}
    require(set(source_slopes).isdisjoint(two_nonzero_lines),
            "source/two-current projective collision")
    require(all(M.rank_mod((source_by_label[label], (1, B))) == 2
                for label in source_labels), "multi-stranger gate rank")

    source_atlas = tuple((label, source_slopes[index])
                         for index, label in enumerate(source_labels))
    exceptional_atlas = tuple(
        (label, two_slopes[label])
        for orbit in EXCEPTIONAL_ORBITS for label in orbit
    )
    semantic = digest_json((
        PRIME, P, EXPECTED_PARENT_SHA256,
        source_reconstruction, two_reconstruction, basis_digest,
        source_atlas, zero_labels, generic_labels, exceptional_atlas,
        REVERSAL_ACTION, EXCEPTIONAL_ORBITS, plus_line, minus_line,
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3657 LRC two-current quotient address/reversal gate ==")
    print(f"field=p:{PRIME};labels=source:17,two_current:169")
    print(f"parent_sha256_lf={EXPECTED_PARENT_SHA256}")
    print(f"reconstruction_sha256=(source:{source_reconstruction},two:{two_reconstruction})")
    print(f"correction_plane=rank:2;rref_sha256:{basis_digest}")
    print(f"source_atlas={source_atlas};nonzero_distinct_lines=17")
    print("two_kernel_law=r1=6 OR (r1=0 AND r0!=12) OR (r1=12 AND r0!=0)")
    print(f"two_kernel=count:{len(zero_labels)};labels_sha256:{digest_json(zero_labels)}")
    print(f"two_generic=count:{len(generic_labels)};slope:{B};labels_sha256:{digest_json(generic_labels)}")
    print(f"two_exceptional_atlas={exceptional_atlas};orbits={EXCEPTIONAL_ORBITS}")
    print(f"two_nonzero_projective_classes=9;sizes=(124,1,1,1,1,1,1,1,1)")
    print(f"reversal=label:(r0,r1)->(12-r0,12-r1);action:{REVERSAL_ACTION};square:I")
    print(f"reversal_eigenlines=plus:{plus_line},minus:{minus_line};minus_realized=False")
    print("cross_source_two_projective_intersection=empty")
    print("multi_stranger_gate=any source cancellation by two-current rows uses an exceptional address")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("status=FINITE-EXACT+VERIFIED-EXACT over the pinned prime")
    print("scope=static quotient atlas;not chronology/current/entry/characteristic-zero/LRC14")
    print("PASS")


if __name__ == "__main__":
    main()

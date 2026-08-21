#!/usr/bin/env python3
"""Exact ANOVA graph flag for the common A4 plane of THM-3585."""

from __future__ import annotations

from hashlib import sha256
import importlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PARENT_PATH = ROOT / "04-computation/lrc_common_a4_channel_plane_centering_complement_thm3585.py"
EXPECTED_PARENT_SHA256 = "dad78dd317a25958f42b3df44f395d4da344330b07ea4f5041a7e838e46dd0e8"

EXPECTED_DIGESTS = {
    "raw": "1d9293d05fa3551b785a1537e78bc8be585fcc43dbb5172036c9b32546ca8560",
    "constant": "a29ec8f53318e8d984e6f2ab4d7bfe6d26c6a4523dea113576347d350cf4b9f4",
    "state": "8ade1771d2e173068be16e186ea588bea3afbeccdb37b6e15da0133bf41f21c8",
    "relation": "c92171956f95eeb8cbfd53bad43299c55b3278270a33bea4b3baf26f236a3412",
    "interaction": "0cfa2e3330f92ab59fd183e5664715c490a702df2ad74491a8180793cae4a21e",
    "additive": "3e7a63217755ab17e3dc5cfbbb9d8cc6e4efdf721cecd532814574764777152a",
    "state_side": "2d82d5adbb3c5a71035ab622e36e51e2a6a701a941a9f6f31ca308b8de13d17c",
    "relation_side": "0c9bccb7f9a354cf316d7ea8be0bf12e6d7007ce25ac69926f95c1527a6ea0c3",
    "graph_stack": "9f4ec33d31337b0100a55871f6284443c6e5cfc0f5133a1493f807d563670821",
}
EXPECTED_SEMANTIC_SHA256 = "44d2d062a447eb68fc58b33eee8d23fb8092e1829534fb7c4f9e9f082e037f76"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


require(lf_sha256(PARENT_PATH) == EXPECTED_PARENT_SHA256, "THM-3585 parent hash")
sys.path.insert(0, str(PARENT_PATH.parent))
M = importlib.import_module(PARENT_PATH.stem)
P, V, PRIME = M.P, M.V, M.PRIME


def add_rows(*families):
    require(families and all(len(family) == len(families[0]) for family in families),
            "row-family arity")
    return tuple(
        tuple(sum(values) % PRIME for values in zip(*(family[i] for family in families)))
        for i in range(len(families[0]))
    )


def anova_components(row):
    require(len(row) == V * P, "ANOVA row width")
    inverse_p = pow(P, -1, PRIME)
    inverse_v = pow(V, -1, PRIME)
    inverse_vp = pow(V * P, -1, PRIME)
    matrix = [row[state * P:(state + 1) * P] for state in range(V)]
    grand = sum(map(sum, matrix)) * inverse_vp % PRIME
    state_means = [sum(matrix[state]) * inverse_p % PRIME for state in range(V)]
    relation_means = [
        sum(matrix[state][relation] for state in range(V)) * inverse_v % PRIME
        for relation in range(P)
    ]
    constant = tuple(grand for _state in range(V) for _relation in range(P))
    state = tuple(
        (state_means[s] - grand) % PRIME
        for s in range(V) for _relation in range(P)
    )
    relation = tuple(
        (relation_means[r] - grand) % PRIME
        for _state in range(V) for r in range(P)
    )
    interaction = tuple(
        (matrix[s][r] - state_means[s] - relation_means[r] + grand) % PRIME
        for s in range(V) for r in range(P)
    )
    rebuilt = tuple(
        sum(values) % PRIME
        for values in zip(constant, state, relation, interaction)
    )
    require(rebuilt == tuple(value % PRIME for value in row), "ANOVA reconstruction")
    return constant, state, relation, interaction


def component_families(raw):
    rows = tuple(anova_components(row) for row in raw)
    families = tuple(tuple(row[index] for row in rows) for index in range(4))
    zero = (0,) * (V * P)
    # Each component is an actual idempotent projection, not merely a formula
    # whose four outputs happen to add back to the row.
    for index, family in enumerate(families):
        for row in family:
            projected = anova_components(row)
            for other, value in enumerate(projected):
                require(value == (row if other == index else zero),
                        ("ANOVA idempotent", index, other))
    return families


def standard_spaces():
    pure_state = tuple(
        tuple(
            (int(state == target) - int(state == V - 1)) % PRIME
            for state in range(V) for _relation in range(P)
        )
        for target in range(V - 1)
    )
    state_only = tuple(
        tuple(int(state == target) for state in range(V) for _relation in range(P))
        for target in range(V)
    )
    pure_relation = tuple(
        tuple(
            (int(relation == target) - int(relation == P - 1)) % PRIME
            for _state in range(V) for relation in range(P)
        )
        for target in range(P - 1)
    )
    relation_only = tuple(
        tuple(int(relation == target) for _state in range(V) for relation in range(P))
        for target in range(P)
    )
    require(tuple(map(M.rank_mod, (pure_state, state_only, pure_relation, relation_only)))
            == (3, 4, 12, 13), "standard ANOVA dimensions")
    return pure_state, state_only, pure_relation, relation_only


def main() -> None:
    tensor, reconstruction = M.reconstruct_source_current()
    require(reconstruction == M.EXPECTED_PARENT_RECONSTRUCTION_SHA256[2:],
            "source reconstruction digests")
    raw = M.flatten_source_current(tensor)
    constant, state, relation, interaction = component_families(raw)
    additive = add_rows(constant, state, relation)
    state_side = add_rows(constant, state)
    relation_side = add_rows(constant, relation)
    pure_state, state_only, pure_relation, relation_only = standard_spaces()

    named = {
        "raw": raw,
        "constant": constant,
        "state": state,
        "relation": relation,
        "interaction": interaction,
        "additive": additive,
        "state_side": state_side,
        "relation_side": relation_side,
    }
    ranks = {name: M.rank_mod(rows) for name, rows in named.items()}
    digests = {name: M.rowspace_digest(rows) for name, rows in named.items()}
    require(ranks == {
        "raw": 4, "constant": 1, "state": 2, "relation": 4,
        "interaction": 4, "additive": 4, "state_side": 2,
        "relation_side": 4,
    }, ("component ranks", ranks))
    require(digests == {name: EXPECTED_DIGESTS[name] for name in named},
            ("component digests", digests))

    stack_ranks = {
        "raw_additive": M.rank_mod(raw + additive),
        "raw_interaction": M.rank_mod(raw + interaction),
        "additive_interaction": M.rank_mod(additive + interaction),
        "additive_pure_state": M.rank_mod(additive + pure_state),
        "additive_state_only": M.rank_mod(additive + state_only),
        "additive_pure_relation": M.rank_mod(additive + pure_relation),
        "additive_relation_only": M.rank_mod(additive + relation_only),
    }
    require(stack_ranks == {
        "raw_additive": 8,
        "raw_interaction": 8,
        "additive_interaction": 8,
        "additive_pure_state": 7,
        "additive_state_only": 8,
        "additive_pure_relation": 14,
        "additive_relation_only": 15,
    }, ("stack ranks", stack_ranks))

    raw_stack = M.canonical_row_basis(raw + interaction)
    graph_stack = M.canonical_row_basis(additive + interaction)
    require(raw_stack == graph_stack, "raw/additive graph spans differ")
    graph_digest = digest_json(graph_stack)
    require(graph_digest == EXPECTED_DIGESTS["graph_stack"],
            ("graph stack digest", graph_digest))

    intersection_dimensions = {
        "raw_additive": 4 + 4 - stack_ranks["raw_additive"],
        "raw_interaction": 4 + 4 - stack_ranks["raw_interaction"],
        "additive_pure_state": 4 + 3 - stack_ranks["additive_pure_state"],
        "additive_state_only": 4 + 4 - stack_ranks["additive_state_only"],
        "additive_pure_relation": 4 + 12 - stack_ranks["additive_pure_relation"],
        "additive_relation_only": 4 + 13 - stack_ranks["additive_relation_only"],
    }
    require(intersection_dimensions == {
        "raw_additive": 0, "raw_interaction": 0,
        "additive_pure_state": 0, "additive_state_only": 0,
        "additive_pure_relation": 2, "additive_relation_only": 2,
    }, ("intersection dimensions", intersection_dimensions))

    semantic = digest_json((
        PRIME, V, P, EXPECTED_PARENT_SHA256, reconstruction,
        tuple(sorted(ranks.items())), tuple(sorted(digests.items())),
        tuple(sorted(stack_ranks.items())),
        tuple(sorted(intersection_dimensions.items())), graph_digest,
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== THM-3593 LRC common A4 ANOVA graph flag ==")
    print(f"field=(prime={PRIME},state={V},relation={P},ambient={V*P})")
    print(f"parent_sha256_lf={EXPECTED_PARENT_SHA256}")
    print(f"component_ranks={ranks}")
    print(f"component_digests={digests}")
    print(f"stack_ranks={stack_ranks}")
    print(f"intersection_dimensions={intersection_dimensions}")
    print(f"graph_stack_equal=True;graph_stack_sha256={graph_digest}")
    print("flag=raw4_graph(additive4,interaction4);additive4_graph(relation4,state_correction2)")
    print("pure_relation_kernel_dimension=2;pure_state_intersection_dimension=0")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("status=FINITE-EXACT+VERIFIED-EXACT over the pinned prime")
    print("scope=static ANOVA flag;not chronology/current/entry/characteristic-zero/LRC14")
    print("PASS")


if __name__ == "__main__":
    main()

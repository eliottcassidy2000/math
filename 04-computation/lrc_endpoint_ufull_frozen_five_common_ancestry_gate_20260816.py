#!/usr/bin/env python3
"""Exact gate audit for the five frozen refined U_full role residues.

The incoming 2,197-twist computation proves that five endpoint-aggregate
values have nonzero graph factors after reduction at a certified split prime.
This sidecar asks the next, differently typed question: do those five frozen
residues by themselves realize a lawful nonnegative Boolean/common-ancestry
current?

Two exact controls localize the answer.

* A five-atom weighted Boolean base realizes the residue representatives as
  vertex potentials and preserves all 72 nonzero graph products.  Thus there
  is no scalar or graph-polynomial obstruction.
* Reduction modulo the split prime forgets positivity, and separated endpoint
  marginals forget their coupling.  A two-atom Boolean hostile has identical
  left/right marginals and identical product of marginals but common-base
  intersection 1/2 versus 0.  Therefore neither the frozen residues nor the
  endpoint engine's post-marginalization product certifies lawful ancestry.

This is a data/API obstruction, not a proof that U_full has no lawful
realization.  The missing object is an atom-level coupling supported on a
typed U_full ancestry relation before endpoint marginalization and reduction.
"""

from __future__ import annotations

import ast
import hashlib
import json
from fractions import Fraction
from itertools import permutations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_endpoint_ufull_frozen_five_common_ancestry_gate_20260816.py"
OUTPUT = "05-knowledge/results/lrc_endpoint_ufull_frozen_five_common_ancestry_gate_20260816.out"
EXPECTED_SEMANTIC_SHA256 = "674f52e86cabcc72d97a193528434cc9524ae3111439e53e49a43a811d4bf41a"

PINS = (
    (
        "PRIMARY-2197-SCRIPT",
        ROOT / "04-computation/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.py",
        "ee2105742abee578a9c41ff7ec954a07ada324fccc2c643429e7ac6e6e6f8fc2",
        None,
    ),
    (
        "PRIMARY-2197-OUTPUT",
        ROOT / "05-knowledge/results/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.out",
        "10a98351cc59615a5b6d2b8f555e0936d1a39566d9906127edc2b0fbc3918e73",
        "STATUS=PASS",
    ),
    (
        "INDEPENDENT-GRAPH-AUDIT-SCRIPT",
        ROOT / "04-computation/lrc14_guard_deleted_refined_endpoint_role_graph_audit_20260816.py",
        "b75f5ef933c18c07b4fd2c4812fea468a9faaa6c5847a5c6a11190cc04676261",
        None,
    ),
    (
        "INDEPENDENT-GRAPH-AUDIT-OUTPUT",
        ROOT / "05-knowledge/results/lrc14_guard_deleted_refined_endpoint_role_graph_audit_20260816.out",
        "d4f5e4b12854bd5802842df2784d580abdd75d9a4fc58dd77db3a185a10403f1",
        "STATUS=PASS",
    ),
    (
        "THM-2471-BOOLEAN-STALK",
        ROOT / "01-canon/theorems/THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary.md",
        "423882a4e3b06e2020fd7f3d0748fc9e904ced4d3214ef8cd77bfd5c6b9597eb",
        "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED",
    ),
    (
        "THM-2538-COMMON-ANCESTRY-GATE",
        ROOT / "01-canon/theorems/THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary.md",
        "b3179e0849f581c6f5f5a44dab9b7d57d7c87f96d71259863315b23a241ae614",
        "PROVED + VERIFIED-EXACT + INDEPENDENTLY SPLIT-AUDITED",
    ),
    (
        "THM-2594-REALIZED-CANONICAL-INSTANCE",
        ROOT / "01-canon/theorems/THM-2594-realized-theta-slaved-contraction-at-the-r5-window.md",
        "b46e955c4dfc7fb1019791db82295643042c319eed52c02ed8d89573d4725cdc",
        "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED",
    ),
)

PRIME = 572252886246508880869
CLASS_VALUES = {
    (0, 0, 0): 405336876493642499425,
    (0, 1, 0): 518539850465495448196,
    (1, 0, 0): 503604956476841920373,
    (1, 0, 1): 320618948602619577408,
    (1, 12, 0): 15703541686881447885,
}
ROLE_CLASSES = {
    "c1": (0, 0, 0),
    "c2": (1, 0, 0),
    "c3": (0, 1, 0),
    "H": (1, 0, 1),
    "q2": (1, 12, 0),
    "q3": (1, 0, 0),
    "q4": (1, 0, 0),
    "q5": (1, 0, 0),
}
ROLE_VALUES = {label: CLASS_VALUES[address] for label, address in ROLE_CLASSES.items()}

EDGES = (
    (0, 3), (0, 4), (0, 5),
    (1, 2), (1, 4), (1, 7),
    (2, 4), (2, 7),
    (3, 4), (3, 5),
    (4, 5), (4, 6), (4, 7),
)
HUB = 4
LEAF = 6
WINGS = ((0, 3, 5), (1, 2, 7))
BLOCKERS = ("c1", "c2", "c3")
UNITS = ("q2", "q3", "q4")
EXPECTED_CHART_SHA256 = "b7d8c2c9860e4f1aa542b1c85fdb7b65cf4985aba5a81a84ff3a324834d51c51"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(
        path.read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(body).hexdigest()


def role_charts() -> tuple[tuple[str, ...], ...]:
    charts = []
    for swap in (0, 1):
        for blockers in permutations(BLOCKERS):
            for units in permutations(UNITS):
                chart = {HUB: "H", LEAF: "q5"}
                chart.update(zip(WINGS[swap], blockers))
                chart.update(zip(WINGS[1 - swap], units))
                charts.append(tuple(chart[index] for index in range(8)))
    answer = tuple(sorted(set(charts)))
    require(len(answer) == 72, ("role chart count", len(answer)))
    return answer


def bareiss_determinant(matrix: list[list[int]]) -> int:
    size = len(matrix)
    if size == 0:
        return 1
    work = [list(row) for row in matrix]
    sign = 1
    denominator = 1
    for column in range(size - 1):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = (
                    work[row][target] * pivot_value
                    - work[row][column] * work[column][target]
                )
                require(numerator % denominator == 0,
                        ("Bareiss divisibility", column, row, target))
                work[row][target] = numerator // denominator
            work[row][column] = 0
        denominator = pivot_value
    return sign * work[-1][-1]


def k4_tree_sum(
    values: dict[str, int], chart: tuple[str, ...], vertices: tuple[int, ...]
) -> int:
    positions = {vertex: index for index, vertex in enumerate(vertices)}
    laplacian = [[0] * 4 for _ in range(4)]
    for left, right in EDGES:
        if left not in positions or right not in positions:
            continue
        weight = values[chart[left]] - values[chart[right]]
        i, j = positions[left], positions[right]
        laplacian[i][i] += weight
        laplacian[j][j] += weight
        laplacian[i][j] -= weight
        laplacian[j][i] -= weight
    return bareiss_determinant([row[:-1] for row in laplacian[:-1]])


def chart_rows(values: dict[str, int]) -> tuple[tuple[object, ...], ...]:
    rows = []
    for chart in role_charts():
        bridge = values["H"] - values["q5"]
        left = k4_tree_sum(values, chart, WINGS[0] + (HUB,))
        right = k4_tree_sum(values, chart, WINGS[1] + (HUB,))
        rows.append((chart, bridge, left, right, bridge * left * right))
    return tuple(rows)


def integral(weight: tuple[int, ...], event: tuple[int, ...]) -> int:
    return sum(value * flag for value, flag in zip(weight, event))


def boolean_vertex_realization() -> tuple[object, ...]:
    atoms = tuple(sorted(CLASS_VALUES))
    weight = tuple(CLASS_VALUES[atom] for atom in atoms)
    events = {
        label: tuple(int(atom == ROLE_CLASSES[label]) for atom in atoms)
        for label in ROLE_CLASSES
    }
    masses = {label: integral(weight, event) for label, event in events.items()}
    require(all(value > 0 for value in weight), "nonpositive synthetic weight")
    require(masses == ROLE_VALUES, ("synthetic role masses", masses))
    rows = chart_rows(masses)
    require(all(row[1] and row[2] and row[3] and row[4] for row in rows),
            "synthetic integer graph factor vanished")
    reduced_rows = tuple(
        (row[0],) + tuple(value % PRIME for value in row[1:]) for row in rows
    )
    require(digest_json(reduced_rows) == EXPECTED_CHART_SHA256,
            ("mod-p chart digest", digest_json(reduced_rows)))
    products = tuple(abs(row[4]) for row in rows)
    return (
        atoms,
        weight,
        tuple((label, events[label]) for label in ROLE_CLASSES),
        tuple((label, masses[label]) for label in ROLE_CLASSES),
        min(products),
        max(products),
        digest_json(rows),
        digest_json(reduced_rows),
    )


def reduction_lift_hostile() -> tuple[object, ...]:
    positive = tuple((address, CLASS_VALUES[address]) for address in sorted(CLASS_VALUES))
    negative_values = dict(CLASS_VALUES)
    negative_values[(1, 0, 1)] -= PRIME
    negative = tuple((address, negative_values[address]) for address in sorted(negative_values))
    require(all(value >= 0 for _address, value in positive), "positive lift")
    require(negative_values[(1, 0, 1)] < 0, "negative lift did not turn negative")
    require(all(
        positive_value % PRIME == negative_values[address] % PRIME
        for address, positive_value in positive
    ), "lift residues differ")
    return (
        positive,
        negative,
        (1, 0, 1),
        negative_values[(1, 0, 1)],
        "NONNEGATIVE_EVENT_MASS_DOES_NOT_DESCEND_MOD_P",
    )


def coupling_hostile() -> tuple[object, ...]:
    measure = (Fraction(1, 2), Fraction(1, 2))
    left = (1, 0)
    right_aligned = (1, 0)
    right_crossed = (0, 1)

    def mass(field: tuple[int, int]) -> Fraction:
        return sum(mu * value for mu, value in zip(measure, field))

    def joint(first: tuple[int, int], second: tuple[int, int]) -> Fraction:
        return sum(
            mu * left_value * right_value
            for mu, left_value, right_value in zip(measure, first, second)
        )

    left_mass = mass(left)
    aligned_mass = mass(right_aligned)
    crossed_mass = mass(right_crossed)
    require((left_mass, aligned_mass, crossed_mass) == (
        Fraction(1, 2), Fraction(1, 2), Fraction(1, 2)
    ), "coupling marginals")
    separate_aligned = left_mass * aligned_mass
    separate_crossed = left_mass * crossed_mass
    joint_aligned = joint(left, right_aligned)
    joint_crossed = joint(left, right_crossed)
    require(separate_aligned == separate_crossed == Fraction(1, 4),
            "separate products")
    require((joint_aligned, joint_crossed) == (Fraction(1, 2), Fraction(0)),
            "joint hostile")
    return (
        measure,
        left,
        right_aligned,
        right_crossed,
        (left_mass, aligned_mass, crossed_mass),
        (separate_aligned, separate_crossed),
        (joint_aligned, joint_crossed),
        "SEPARATE_MARGINALS_DO_NOT_DETERMINE_COMMON_ANCESTRY",
    )


def endpoint_api_audit() -> tuple[object, ...]:
    source_path = PINS[0][1]
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    function_names = tuple(sorted(
        node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)
    ))
    require("worker" in function_names, "primary worker absent")
    require("ax, _ = M.fast_x_sweep(" in source, "left marginal call drift")
    require("by = M.fast_endpoint_sum(" in source, "right marginal call drift")
    require("rows.append(phase * ax[0] % p * by[0] % p)" in source,
            "post-marginalization product drift")
    atom_api = tuple(
        name for name in function_names
        if any(token in name.lower() for token in ("ancestry", "stalk", "joint", "atom"))
    )
    require(atom_api == (), ("unexpected atom-level API", atom_api))
    return (
        ("left", "fast_x_sweep -> one field scalar ax"),
        ("right", "fast_endpoint_sum -> one field scalar by"),
        ("composition", "phase*ax*by after both marginalizations"),
        ("worker_output", "alpha, interval_count, tuple(gamma scalars)"),
        ("atom_or_ancestry_functions", atom_api),
    )


def security_certificate(path: Path) -> tuple[object, ...]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    forbidden_calls = {"eval", "exec", "compile", "__import__"}
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            bad.append("Assert")
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id in forbidden_calls:
                bad.append(node.func.id)
    require(not bad, ("security", bad))
    return len(tuple(ast.walk(tree))), tuple(bad)


def main() -> None:
    pin_rows = []
    for label, path, expected_hash, status_text in PINS:
        actual_hash = lf_sha256(path)
        require(actual_hash == expected_hash, (label, "hash drift", actual_hash))
        if status_text is not None:
            require(status_text in path.read_text(encoding="utf-8"),
                    (label, "status drift"))
        pin_rows.append((label, actual_hash))

    primary_output = PINS[1][1].read_text(encoding="utf-8")
    expected_value_line = (
        "refined_role_values=(((0, 0, 0), 405336876493642499425), "
        "((0, 1, 0), 518539850465495448196), "
        "((1, 0, 0), 503604956476841920373), "
        "((1, 0, 1), 320618948602619577408), "
        "((1, 12, 0), 15703541686881447885))"
    )
    require(expected_value_line in primary_output, "frozen role values drift")
    require("factor_zero_counts=(bridge,left_K4,right_K4,product)=(0, 0, 0, 0)"
            in primary_output, "primary factor census drift")

    synthetic = boolean_vertex_realization()
    lift_hostile = reduction_lift_hostile()
    coupling = coupling_hostile()
    api = endpoint_api_audit()
    security = security_certificate(ROOT / SCRIPT)

    consequence = (
        tuple(pin_rows),
        PRIME,
        tuple((address, CLASS_VALUES[address]) for address in sorted(CLASS_VALUES)),
        tuple((label, ROLE_CLASSES[label]) for label in ROLE_CLASSES),
        synthetic,
        lift_hostile,
        coupling,
        api,
    )
    semantic_hash = hashlib.sha256(repr(consequence).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic_hash))

    print("LRC U_FULL FROZEN-FIVE COMMON-ANCESTRY GATE AUDIT")
    print("status=FINITE-EXACT scope gate for proved THM-3479; U_full endpoint nonvanishing inherited; physical current and LRC(14) OPEN")
    print(f"script={SCRIPT}")
    print(f"stored_output={OUTPUT}")
    print(f"dependency_hashes={tuple(pin_rows)}")
    print(f"frozen_field_data=(prime={PRIME},classes={tuple((address, CLASS_VALUES[address]) for address in sorted(CLASS_VALUES))})")
    print(f"synthetic_boolean_base=(atoms={synthetic[0]},positive_weights={synthetic[1]},role_masses={synthetic[3]})")
    print(f"synthetic_integer_graph=(zero_products=0/72,min_abs_product={synthetic[4]},max_abs_product={synthetic[5]},integer_chart_sha256={synthetic[6]},mod_p_chart_sha256={synthetic[7]})")
    print("SCALAR_GATE=PASS: five disjoint weighted Boolean atoms realize the residue representatives and all 72 graph products")
    print(f"reduction_lift_hostile=(negative_class={lift_hostile[2]},negative_lift={lift_hostile[3]},verdict={lift_hostile[4]})")
    print(f"coupling_hostile=(marginals={coupling[4]},separate_products={coupling[5]},common_base_intersections={coupling[6]},verdict={coupling[7]})")
    print(f"endpoint_worker_api={api}")
    print("FIRST_UNAVAILABLE_COORDINATE=a shared U_full ancestry key linking one left fast_x_sweep atom to one right fast_endpoint_sum atom before either sum")
    print("MISSING_API=an atom-level coupling J_ell(omega_L,omega_R) supported on the lawful U_full ancestry relation, with both endpoint factors inserted before marginalization and with inverse DFT/reduction recovering the frozen bridge classes")
    print("WHY_THM2594_DOES_NOT_FILL_IT=THM-2594 supplies one linked-node contraction on the different canonical row and explicitly is not the THM-2334 endpoint bank or a generic transplant")
    print("CHEAPEST_NEXT_TEST=construct only the two-class bridge coupling q_H=(1,0,1), q_q5=(1,0,0); require its atom-level inverse-DFT difference to reduce to 389266878372286537904 mod p before testing either K4")
    print("PRESERVED=five labels/classes, positive synthetic vertex masses, graph orbit, finite-field nonvanishing")
    print("LOST_BY_REDUCTION=order, sign, positivity, characteristic-zero lift; LOST_BY_MARGINALIZATION=left/right pairing, ancestry support, chronology, owner/deep semantics")
    print("NONCONSEQUENCES=no obstruction to existence of a U_full lawful coupling, no physical current, no grouped C(a;X,m), no B(q), no scalar-row exclusion, no U_clock statement, no LRC(14)")
    print(f"semantic_sha256={semantic_hash}")
    print(f"security_ast_nodes_forbidden={security}")
    print("VERDICT=SCALAR_BOOLEAN_REALIZATION_EXISTS_BUT_THE_LAWFUL_COMMON_ANCESTRY_GATE_IS_NOT_DECIDABLE_FROM_THE_FROZEN_FIVE_OR_POST_MARGINALIZATION_ENDPOINT_API")


if __name__ == "__main__":
    main()

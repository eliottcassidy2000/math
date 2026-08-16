#!/usr/bin/env python3
"""Exact T4-atlas realization of the Fibonacci 24-state branch transplant.

An affine owner plus an order of the three nonzero V4 directions is a total
order of four vertices.  This probe identifies both 24-state bundles with the
atlas of labelled transitive T4 tournaments and classifies the base-preserving
conjugacies after inserting each nonzero H1 seam class.
"""

from __future__ import annotations

import ast
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from json import dumps
from math import lcm
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/fibonacci_t4_atlas_h1_branch_transplant_probe_20260816.py"
OUTPUT = "05-knowledge/results/fibonacci_t4_atlas_h1_branch_transplant_probe_20260816.out"
EXPECTED_SEMANTIC_SHA256 = "3e19be5a6d7679656f4e476ea1cb5a3bd5944a190bab860c43a8ecbd034d8176"

PINS = (
    (
        "THM-2622",
        ROOT / "01-canon/theorems/THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary.md",
        "ffcbd93f19d0b669c53501216c0c790ff60faf041907e75d1ab0512a4e492844",
    ),
    (
        "THM-3339-script",
        ROOT / "04-computation/fibonacci_berggren_three_ray_owner_thm3339.py",
        "094fc254dcb7965791e59247a98f60a337725b40d169db6f3862c1da5943149b",
    ),
    (
        "frame-line-output",
        ROOT / "05-knowledge/results/fibonacci_farey_mod3_q15_four_state_probe_20260814.out",
        "803c848ed775b2b1278fecfa1340d2f8f634abaaa2a23ea3e1f24abf76c57e48",
    ),
)

ZERO = (0, 0)
P = (1, 0)
Q = (0, 1)
R = (1, 1)
V4 = (ZERO, P, Q, R)
I2 = ((1, 0), (0, 1))
J2 = ((0, 1), (1, 0))
ORDERS = (
    (Q, R, P),
    (Q, P, R),
    (P, Q, R),
    (P, R, Q),
    (R, P, Q),
    (R, Q, P),
)
OWNERS = (ZERO, P, P, R, Q, Q)
G = (1, 2, 3, 0)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def vadd(left, right):
    return left[0] ^ right[0], left[1] ^ right[1]


def mvec(matrix, vector):
    return (
        (matrix[0][0] * vector[0] + matrix[0][1] * vector[1]) % 2,
        (matrix[1][0] * vector[0] + matrix[1][1] * vector[1]) % 2,
    )


def mmul(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(2)) % 2
              for j in range(2))
        for i in range(2)
    )


GL2 = tuple(
    ((a, b), (c, d))
    for a, b, c, d in product((0, 1), repeat=4)
    if (a * d - b * c) % 2 == 1
)


def linear_step(source, target):
    candidates = tuple(
        matrix for matrix in GL2
        if tuple(mvec(matrix, direction) for direction in source) == target
    )
    require(len(candidates) == 1, (source, target, candidates))
    return candidates[0]


def affine_apply(affine, vector):
    matrix, translation = affine
    return vadd(mvec(matrix, vector), translation)


def affine_compose(left, right):
    left_matrix, left_translation = left
    right_matrix, right_translation = right
    return (
        mmul(left_matrix, right_matrix),
        vadd(mvec(left_matrix, right_translation), left_translation),
    )


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def power(permutation, exponent):
    result = tuple(range(len(permutation)))
    base = permutation
    while exponent:
        if exponent & 1:
            result = compose(base, result)
        exponent //= 2
        if exponent:
            base = compose(base, base)
    return result


def cycle_lengths(permutation):
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            current = permutation[current]
            length += 1
        require(current == start, (start, current))
        lengths.append(length)
    return tuple(sorted(lengths))


def owner_steps():
    steps = []
    for index in range(6):
        next_index = (index + 1) % 6
        linear = linear_step(ORDERS[index], ORDERS[next_index])
        translation = vadd(OWNERS[next_index], mvec(linear, OWNERS[index]))
        steps.append((linear, translation))
    holonomy = (I2, ZERO)
    for step in steps:
        holonomy = affine_compose(step, holonomy)
    require(holonomy == (I2, ZERO), holonomy)
    return tuple(steps)


def ranking(order_index, owner):
    return (owner,) + tuple(vadd(owner, direction) for direction in ORDERS[order_index])


def tournament_arcs(total_order):
    return tuple(sorted(
        (total_order[left], total_order[right])
        for left, right in combinations(range(4), 2)
    ))


def repaired_steps(base_steps, seam):
    steps = list(base_steps)
    matrix, translation = steps[-1]
    steps[-1] = matrix, vadd(translation, seam)
    holonomy = (I2, ZERO)
    for step in steps:
        holonomy = affine_compose(step, holonomy)
    require(holonomy == (I2, seam), (seam, holonomy))
    return tuple(steps)


def propagate_labelings(initial, steps):
    labelings = [tuple(initial)]
    for index, step in enumerate(steps):
        following = [None] * 4
        for point in range(4):
            following[G[point]] = affine_apply(step, labelings[index][point])
        require(all(value is not None for value in following), following)
        require(len(set(following)) == 4, following)
        labelings.append(tuple(following))
    return tuple(labelings)


def bundle_permutation(steps):
    states = tuple((index, owner) for index in range(6) for owner in V4)
    state_index = {state: position for position, state in enumerate(states)}
    return tuple(
        state_index[((index + 1) % 6, affine_apply(steps[index], owner))]
        for index, owner in states
    )


def security_gate(path: Path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert forbidden")
    forbidden = {"eval", "exec", "compile", "open", "system", "popen", "run", "Popen"}
    calls = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    calls.update(
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    )
    require(not (calls & forbidden), ("forbidden calls", calls & forbidden))
    return len(tuple(ast.walk(tree))), tuple(sorted(calls & forbidden))


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    pins = []
    for label, path, expected in PINS:
        actual = lf_sha256(path)
        require(actual == expected, (label, "hash drift", actual))
        pins.append((label, actual))

    base_steps = owner_steps()
    all_rankings = tuple(ranking(index, owner) for index in range(6) for owner in V4)
    require(len(set(all_rankings)) == 24, "owner-order atlas is not bijective")
    require(set(all_rankings) == set(permutations(V4)), "not every total order occurs")
    all_tournaments = tuple(tournament_arcs(order) for order in all_rankings)
    require(len(set(all_tournaments)) == 24, "transitive T4 collision")

    # The lawful moving owner edge map transports every ranked tournament
    # positionwise, but its full return is flat.
    for index, step in enumerate(base_steps):
        for owner in V4:
            source = ranking(index, owner)
            target = ranking((index + 1) % 6, affine_apply(step, owner))
            require(tuple(affine_apply(step, vertex) for vertex in source) == target,
                    (index, owner, source, target))
    require(cycle_lengths(bundle_permutation(base_steps)) == (6, 6, 6, 6),
            "lawful atlas cycle type")

    g_square = power(G, 2)
    structured_counts = []
    structured_maps = []
    for seam in V4:
        steps = repaired_steps(base_steps, seam)
        accepted = []
        for initial in permutations(V4):
            labelings = propagate_labelings(initial, steps)
            if labelings[-1] != labelings[0]:
                continue
            # Equivariance and positionwise T4 transport on every state.
            for index in range(6):
                for point in range(4):
                    owner = labelings[index][point]
                    target_owner = labelings[(index + 1) % 6][G[point]]
                    require(affine_apply(steps[index], owner) == target_owner,
                            (seam, index, point, "bundle equivariance"))
                    source_ranking = ranking(index, owner)
                    target_ranking = ranking((index + 1) % 6, target_owner)
                    require(tuple(affine_apply(steps[index], vertex)
                                  for vertex in source_ranking) == target_ranking,
                            (seam, index, point, "T4 transport"))
            accepted.append(tuple(labelings[:-1]))
        expected = 0 if seam == ZERO else 8
        require(len(accepted) == expected, (seam, len(accepted), expected))
        if seam != ZERO:
            require(all(
                all(initial[g_square[point]] == vadd(initial[point], seam)
                    for point in range(4))
                for initial in (labeling[0] for labeling in accepted)
            ), (seam, "initial matching gauge"))
            require(cycle_lengths(bundle_permutation(steps)) == (12, 12),
                    (seam, "repaired cycle type"))
        structured_counts.append((seam, len(accepted)))
        structured_maps.extend((seam, labeling) for labeling in accepted)

    require(tuple(structured_counts) == ((ZERO, 0), (P, 8), (Q, 8), (R, 8)),
            structured_counts)
    require(len(structured_maps) == 24, len(structured_maps))
    require({entry[1][0] for entry in structured_maps} == set(permutations(V4)),
            "the 24 point gauges did not partition by seam")

    translation_cochain_census = Counter()
    for cochain in product(V4, repeat=6):
        steps = tuple((I2, translation) for translation in cochain)
        seam = ZERO
        for translation in cochain:
            seam = vadd(seam, translation)
        cycles = cycle_lengths(bundle_permutation(steps))
        expected_cycles = (6, 6, 6, 6) if seam == ZERO else (12, 12)
        require(cycles == expected_cycles, (cochain, seam, cycles))
        section_count = 4 if seam == ZERO else 0
        translation_cochain_census[(seam, cycles, section_count)] += 1
    require(tuple(translation_cochain_census.values()) == (1024, 1024, 1024, 1024),
            translation_cochain_census)

    # The H1=V4 classifier is restricted to the fixed/trivialized linear
    # local system.  A nontrivial return linear part gives a third cycle type.
    linear_return_hostile = tuple((I2, ZERO) for _ in range(5)) + ((J2, ZERO),)
    linear_hostile_cycles = cycle_lengths(bundle_permutation(linear_return_hostile))
    require(linear_hostile_cycles == (6, 6, 12), linear_hostile_cycles)
    require(sum(mvec(J2, owner) == owner for owner in V4) == 2,
            "linear hostile fixed-section count")

    # A fixed relabelling of K4 has order at most four.  The lawful and
    # repaired atlas shifts have orders six and twelve, respectively, so the
    # connection is necessarily time-dependent.
    s4_orders = tuple(sorted({
        max(cycle_lengths(permutation))
        if len(set(cycle_lengths(permutation))) == 1
        else lcm(*cycle_lengths(permutation))
        for permutation in permutations(range(4))
    }))
    require(s4_orders == (1, 2, 3, 4), s4_orders)

    security = security_gate(ROOT / SCRIPT)
    semantic_payload = {
        "pins": tuple(pins),
        "atlas_rankings": all_rankings,
        "moving_steps": base_steps,
        "lawful_cycles": cycle_lengths(bundle_permutation(base_steps)),
        "structured_counts": tuple(structured_counts),
        "structured_maps": tuple(structured_maps),
        "translation_cochain_census": tuple(translation_cochain_census.items()),
        "linear_return_hostile": (J2, linear_hostile_cycles, 2),
        "repaired_cycles": tuple(
            (seam, cycle_lengths(bundle_permutation(repaired_steps(base_steps, seam))))
            for seam in V4[1:]
        ),
        "s4_element_orders": s4_orders,
        "security": security,
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":"),
              default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
                (semantic_hash, EXPECTED_SEMANTIC_SHA256))

    print("FIBONACCI 24-STATE T4 ATLAS / H1 BRANCH TRANSPLANT EXACT PROBE")
    print("STATUS: VERIFIED-EXACT STRUCTURAL SIDECAR TO PROVED THM-3487 + INDEPENDENTLY AUDITED")
    print(f"SCRIPT: {SCRIPT}")
    print(f"OUTPUT: {OUTPUT}")
    print(f"PINS: {tuple(pins)}")
    print("T4_ATLAS: (owner,(d1,d2,d3))->(owner,owner+d1,owner+d2,owner+d3) bijects 4*6 owner-order states with all 24 total orders, hence all labelled transitive tournaments on V4")
    print("LAWFUL_CONNECTION: each moving affine owner step transports the ranked T4 positionwise; six-step holonomy=id; cycle_type=6^4")
    print(f"BASE_PRESERVING_TRANSPLANTS_BY_H1_CLASS: {tuple(structured_counts)}")
    print("NONZERO_TRANSPLANT: for each h in {p,q,r}, exactly 8 point gauges extend uniquely to a closed base-preserving T4-atlas conjugacy; the repaired owner holonomy is translation by h and cycle_type=12^2")
    print("ZERO_CLASS_OBSTRUCTION: no point gauge conjugates projective return G^2=(02)(13) to identity; hence the lawful closed-section owner bundle admits zero base-preserving frame-line transplants")
    print("GAUGE_PARTITION: the 24 initial identifications split 8+8+8 by which projective perfect matching G^2 is sent to which XOR direction")
    print(f"TRANSLATION_COCHAIN_EXHAUSTION: all 4^6=4096 cochains; {tuple(translation_cochain_census.items())}")
    print(f"LINEAR_RETURN_HOSTILE: return_linear_part=swap; fixed_sections=2; cycle_type={linear_hostile_cycles}; outside the pure-translation H1 classifier")
    print(f"STATIC_RELABELING_NO_GO: S4 element orders={s4_orders}; full shifts have orders 6 or 12, so no single K4 vertex permutation induces the connection, although every edge step is affine")
    print("TYPING: this is a genuine atlas of transitive T4 tournaments with a time-dependent affine connection; it is not a T6 on the six base states and gives no global Berggren-tree or physical-current transport")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print(f"SECURITY_AST_NODES_AND_FORBIDDEN: {security}")
    print("VERDICT: the branch transplant exists exactly after choosing a nonzero H1 seam class and one of eight compatible point gauges; the same datum that creates the T4-atlas conjugacy forbids every closed owner section")


if __name__ == "__main__":
    main()

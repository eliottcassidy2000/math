#!/usr/bin/env python3
"""Exact finite controls for the K4 hafnian initial-face sidecar.

The mathematical statement is the ordinary associated-graded addition law
for a three-term hafnian.  The finite controls exhaust truncated p-adic
rings, audit the sharp lift ambiguity when the initial face sum vanishes,
and check the S4/V4 matching action and common vertex-unit gauge.

Truth-bearing checks use explicit raises, so ``python -O`` performs the same
work.
"""

from __future__ import annotations

import hashlib
from itertools import permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}
MATCHINGS = ((0, 5), (1, 4), (2, 3))


def vp_integer(x, p):
    require(x != 0, "integer valuation requested at zero")
    x = abs(x)
    ans = 0
    while x % p == 0:
        ans += 1
        x //= p
    return ans


def truncated_datum(x, p, depth):
    modulus = p**depth
    x %= modulus
    if x == 0:
        return depth, 0
    value = 0
    while x % p == 0:
        value += 1
        x //= p
    return value, x % p


def face_sum(values, residues, p):
    minimum = min(values)
    face = tuple(i for i, value in enumerate(values) if value == minimum)
    sigma = sum(residues[i] for i in face) % p
    return minimum, face, sigma


def matching_products(edge_values):
    return tuple(edge_values[i] * edge_values[j] for i, j in MATCHINGS)


def matching_action(g):
    pair_sets = [frozenset((EDGES[i], EDGES[j])) for i, j in MATCHINGS]
    action = []
    for pair_set in pair_sets:
        moved = frozenset(tuple(sorted((g[a], g[b]))) for a, b in pair_set)
        action.append(pair_sets.index(moved))
    return tuple(action)


def clutch(values):
    return tuple(value % 2 for value in values), sum(values) % 3


def finite_field_zero_counts(q):
    counts = {}
    for size in (1, 2, 3):
        count = sum(
            1
            for entries in product(range(1, q), repeat=size)
            if sum(entries) % q == 0
        )
        expected = {1: 0, 2: q - 1, 3: (q - 1) * (q - 2)}[size]
        require(count == expected, f"F_{q} zero-sum count changed at size {size}")
        counts[size] = count
    return counts


def audit_truncated_ring(p, depth):
    modulus = p**depth
    elements = range(1, modulus)
    outcomes = {}
    triple_count = 0
    sigma_zero_count = 0
    sigma_nonzero_count = 0
    valuation_exact_count = 0

    for entries in product(elements, repeat=3):
        data = tuple(truncated_datum(entry, p, depth) for entry in entries)
        values = tuple(item[0] for item in data)
        residues = tuple(item[1] for item in data)
        minimum, _, sigma = face_sum(values, residues, p)
        total = sum(entries) % modulus
        total_value = truncated_datum(total, p, depth)[0]

        if sigma:
            require(total_value == minimum, "nonzero initial face did not fix valuation")
            sigma_nonzero_count += 1
            valuation_exact_count += 1
        else:
            require(total_value > minimum, "zero initial face failed to cancel leading order")
            sigma_zero_count += 1

        signature = values, residues
        mask = outcomes.get(signature, 0)
        mask |= 1 if total == 0 else 2
        outcomes[signature] = mask
        triple_count += 1

    robust_groups = 0
    ambiguous_groups = 0
    top_boundary_groups = 0
    for (values, residues), mask in outcomes.items():
        minimum, _, sigma = face_sum(values, residues, p)
        if sigma:
            require(mask == 2, "nonzero sigma signature admitted a zero lift")
            robust_groups += 1
        elif minimum < depth - 1:
            require(mask == 3, "zero sigma signature did not admit zero and nonzero lifts")
            ambiguous_groups += 1
        else:
            require(mask == 1, "top truncated layer should have only zero sums")
            top_boundary_groups += 1

    require(triple_count == (modulus - 1) ** 3, "truncated universe size changed")
    require(sigma_zero_count + sigma_nonzero_count == triple_count, "sigma partition failed")
    return {
        "p": p,
        "depth": depth,
        "triples": triple_count,
        "signatures": len(outcomes),
        "sigma_nonzero": sigma_nonzero_count,
        "sigma_zero": sigma_zero_count,
        "valuation_exact": valuation_exact_count,
        "robust_groups": robust_groups,
        "ambiguous_groups": ambiguous_groups,
        "top_boundary_groups": top_boundary_groups,
    }


def audit_symmetry():
    actions = [matching_action(g) for g in permutations(range(4))]
    require(len(set(actions)) == 6, "S4 matching image is not S3")
    require(actions.count((0, 1, 2)) == 4, "S4 matching kernel is not V4")

    checks = 0
    for p in (2, 3, 5, 7):
        for values in product(range(3), repeat=3):
            for residues in product(range(1, p), repeat=3):
                minimum, face, sigma = face_sum(values, residues, p)
                for action in set(actions):
                    moved_values = [0, 0, 0]
                    moved_residues = [0, 0, 0]
                    for old, new in enumerate(action):
                        moved_values[new] = values[old]
                        moved_residues[new] = residues[old]
                    moved = face_sum(tuple(moved_values), tuple(moved_residues), p)
                    require(moved[0] == minimum, "matching permutation changed minimum")
                    require(moved[2] == sigma, "matching permutation changed face sum")
                    require(
                        tuple(sorted(moved[1]))
                        == tuple(sorted(action[i] for i in face)),
                        "matching permutation changed the face incorrectly",
                    )
                    checks += 1
    return len(set(actions)), actions.count((0, 1, 2)), checks


def audit_vertex_gauge():
    checks = 0
    for p in (3, 5):
        units = range(1, p)
        for edge_residues in product(units, repeat=6):
            before = matching_products(edge_residues)
            before_sigma = sum(before) % p
            for gauge in product(units, repeat=4):
                moved_edges = tuple(
                    edge_residues[i] * gauge[a] * gauge[b] % p
                    for i, (a, b) in enumerate(EDGES)
                )
                after_sigma = sum(matching_products(moved_edges)) % p
                scale = 1
                for entry in gauge:
                    scale = scale * entry % p
                require(after_sigma == scale * before_sigma % p, "vertex gauge law failed")
                require((after_sigma == 0) == (before_sigma == 0), "gauge changed vanishing")
                checks += 1
    return checks


lines = []


def emit(line):
    lines.append(line)


# The hostile and its arbitrarily deep live perturbations.
hostile_checks = 0
for p in (2, 3, 5, 7):
    for depth in range(1, 7):
        cancel = (1, 1, -2)
        live = (1 + p**depth, 1, -2)
        require(sum(cancel) == 0, "base hostile stopped cancelling")
        require(sum(live) == p**depth, "deep perturbation has wrong amplitude")
        cancel_values = tuple(vp_integer(entry, p) for entry in cancel)
        live_values = tuple(vp_integer(entry, p) for entry in live)
        require(cancel_values == live_values, "deep hostile changed valuation data")
        cancel_residues = tuple(
            (entry // p**value) % p for entry, value in zip(cancel, cancel_values)
        )
        live_residues = tuple(
            (entry // p**value) % p for entry, value in zip(live, live_values)
        )
        require(cancel_residues == live_residues, "deep hostile changed leading residues")
        require(
            all((live[i] - cancel[i]) % p**depth == 0 for i in range(3)),
            "deep hostile changed an earlier jet",
        )
        hostile_checks += 1

# Same clutch and all three leading residues, but different tropical faces.
p = 5
clutch_live = (1, p**2, -2 * p**4)
clutch_zero = (p**2, p**2, -2 * p**2)
live_values = tuple(vp_integer(entry, p) for entry in clutch_live)
zero_values = tuple(vp_integer(entry, p) for entry in clutch_zero)
live_residues = tuple((entry // p**value) % p for entry, value in zip(clutch_live, live_values))
zero_residues = tuple((entry // p**value) % p for entry, value in zip(clutch_zero, zero_values))
require(clutch(live_values) == clutch(zero_values) == ((0, 0, 0), 0), "clutch hostile changed")
require(live_residues == zero_residues == (1, 1, 3), "clutch residue hostile changed")
require(sum(clutch_live) != 0 and sum(clutch_zero) == 0, "clutch face hostile failed")
require(matching_products((clutch_live[0], clutch_live[1], clutch_live[2], 1, 1, 1)) == clutch_live,
        "live triple is not realized by K4 edge weights")
require(matching_products((clutch_zero[0], clutch_zero[1], clutch_zero[2], 1, 1, 1)) == clutch_zero,
        "zero triple is not realized by K4 edge weights")

field_counts = {q: finite_field_zero_counts(q) for q in (2, 3, 5, 7, 11)}
ring_records = [
    audit_truncated_ring(2, 6),
    audit_truncated_ring(3, 4),
    audit_truncated_ring(5, 3),
    audit_truncated_ring(7, 2),
]
matching_image, matching_kernel, symmetry_checks = audit_symmetry()
gauge_checks = audit_vertex_gauge()

emit("K4 HAFNIAN INITIAL-FACE SIDECAR: FINITE-EXACT PASS")
emit("hafnian=A01*A23+A02*A13+A03*A12")
emit("face=argmin(lambda);sigma=sum_face(initial_matching_units) in residue_field")
emit("law=sigma!=0 iff valuation(hafnian)=min(lambda)")
emit("robust_lift_law=sigma!=0 iff every lift with fixed valuations/residues is nonzero")
emit("zero_sigma_boundary=both exact-zero and higher-valuation nonzero lifts exist")
emit("actual_nonvanishing=requires the unbounded filtered contraction jet when sigma=0")
for q, counts in field_counts.items():
    emit(
        f"F{q}_zero_face_counts=size1:{counts[1]};size2:{counts[2]};size3:{counts[3]}"
    )
for record in ring_records:
    emit(
        "truncated="
        + ";".join(f"{key}:{value}" for key, value in record.items())
    )
emit(f"deep_hostile_checks={hostile_checks}")
emit("clutch_hostile=lambda024/live versus lambda222/zero;clutch=(000,0);residues=(1,1,3)")
emit(f"S4_matching_image={matching_image};V4_kernel={matching_kernel};symmetry_checks={symmetry_checks}")
emit(f"vertex_unit_gauge_checks={gauge_checks}")
emit("selector_scope=apply only after THM2290 endpoint-colour selection and pair aggregation")
emit("torsor_scope=initial face is S3-equivariant and V4-blind;no quartic origin is recovered")

semantic = hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()
for line in lines:
    print(line)
print(f"semantic_sha256={semantic}")

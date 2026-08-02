#!/usr/bin/env python3
"""Exact controls for THM-3058's finite-fibre initial-face theorem.

The proof is the associated-graded addition law over a DVR.  These finite
controls exhaust the q=5 residue and truncated-ring cases used by the K4
hafnian specialization, verify the general zero-sum count and lift
construction in bounded banks, and audit the S4/V4 action, vertex-unit gauge,
same-clutch hostile, deep cancellation hostile, and oriented Pluecker wall.

Truth-bearing checks use explicit raises, so ``python -O`` performs the same
work.
"""

from __future__ import annotations

import hashlib
from itertools import permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 5
EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
MATCHINGS = ((0, 5), (1, 4), (2, 3))


def vp_integer(x, p):
    require(x != 0, "integer valuation requested at zero")
    x = abs(x)
    answer = 0
    while x % p == 0:
        answer += 1
        x //= p
    return answer


def leading_residue(x, p):
    value = vp_integer(x, p)
    return value, (x // p**value) % p


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


def initial_face(values, residues, p):
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


def zero_sum_count(q, arity):
    numerator = (q - 1) ** arity + (-1) ** arity * (q - 1)
    require(numerator % q == 0, "closed zero-sum count is not integral")
    return numerator // q


def audit_q5_zero_sum_formula(max_arity=8):
    counts = {}
    tuple_count = 0
    for arity in range(1, max_arity + 1):
        count = sum(
            1
            for entries in product(range(1, P), repeat=arity)
            if sum(entries) % P == 0
        )
        require(
            count == zero_sum_count(P, arity),
            f"F_5 zero-sum formula changed at arity {arity}",
        )
        counts[arity] = count
        tuple_count += (P - 1) ** arity
    return counts, tuple_count


def audit_zero_sigma_lift_constructor(max_arity=5):
    checks = 0
    for arity in range(1, max_arity + 1):
        for values in product(range(3), repeat=arity):
            for residues in product(range(1, P), repeat=arity):
                minimum, face, sigma = initial_face(values, residues, P)
                if sigma != 0:
                    continue

                pivot = face[0]
                terms = [residues[i] * P ** values[i] for i in range(arity)]
                terms[pivot] = -sum(terms[i] for i in range(arity) if i != pivot)
                rebuilt = tuple(leading_residue(term, P) for term in terms)
                require(
                    tuple(value for value, _ in rebuilt) == values,
                    "zero-lift construction changed a valuation",
                )
                require(
                    tuple(residue for _, residue in rebuilt) == residues,
                    "zero-lift construction changed a leading residue",
                )
                require(sum(terms) == 0, "zero-lift construction did not cancel")

                live = list(terms)
                live[pivot] += P ** (minimum + 1)
                live_data = tuple(leading_residue(term, P) for term in live)
                require(live_data == rebuilt, "live perturbation changed first-order data")
                require(
                    sum(live) == P ** (minimum + 1),
                    "live perturbation has the wrong higher valuation",
                )
                checks += 1
    return checks


def audit_q5_truncated_triples(depth=3):
    modulus = P**depth
    elements = range(1, modulus)
    outcomes = {}
    triple_count = 0
    sigma_zero_count = 0
    sigma_nonzero_count = 0
    valuation_exact_count = 0

    for entries in product(elements, repeat=3):
        data = tuple(truncated_datum(entry, P, depth) for entry in entries)
        values = tuple(item[0] for item in data)
        residues = tuple(item[1] for item in data)
        minimum, _, sigma = initial_face(values, residues, P)
        total = sum(entries) % modulus
        total_value = truncated_datum(total, P, depth)[0]

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
        minimum, _, sigma = initial_face(values, residues, P)
        if sigma:
            require(mask == 2, "nonzero sigma signature admitted a zero lift")
            robust_groups += 1
        elif minimum < depth - 1:
            require(mask == 3, "zero sigma signature lacked zero/nonzero lifts")
            ambiguous_groups += 1
        else:
            require(mask == 1, "top truncated layer should contain only truncated zeros")
            top_boundary_groups += 1

    require(triple_count == (modulus - 1) ** 3, "truncated universe size changed")
    require(
        sigma_zero_count + sigma_nonzero_count == triple_count,
        "sigma partition failed",
    )
    return {
        "q": P,
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


def audit_deep_hostile(max_depth=12):
    checks = 0
    for depth in range(1, max_depth + 1):
        cancel = (1, 1, -2)
        live = (1 + P**depth, 1, -2)
        require(sum(cancel) == 0, "base hostile stopped cancelling")
        require(sum(live) == P**depth, "deep perturbation has wrong amplitude")
        cancel_data = tuple(leading_residue(entry, P) for entry in cancel)
        live_data = tuple(leading_residue(entry, P) for entry in live)
        require(cancel_data == live_data, "deep hostile changed first-order data")
        require(
            all((live[i] - cancel[i]) % P**depth == 0 for i in range(3)),
            "deep hostile changed a bounded jet",
        )
        checks += 1
    return checks


def audit_same_clutch_hostile():
    live = (1, P**2, -2 * P**4)
    zero = (P**2, P**2, -2 * P**2)
    live_data = tuple(leading_residue(entry, P) for entry in live)
    zero_data = tuple(leading_residue(entry, P) for entry in zero)
    live_values = tuple(value for value, _ in live_data)
    zero_values = tuple(value for value, _ in zero_data)
    live_residues = tuple(residue for _, residue in live_data)
    zero_residues = tuple(residue for _, residue in zero_data)
    require(
        clutch(live_values) == clutch(zero_values) == ((0, 0, 0), 0),
        "same-clutch hostile changed",
    )
    require(
        live_residues == zero_residues == (1, 1, 3),
        "same-clutch leading residues changed",
    )
    require(sum(live) != 0 and sum(zero) == 0, "same-clutch hostile failed")
    require(
        matching_products((live[0], live[1], live[2], 1, 1, 1)) == live,
        "live triple is not realized by K4 edge weights",
    )
    require(
        matching_products((zero[0], zero[1], zero[2], 1, 1, 1)) == zero,
        "zero triple is not realized by K4 edge weights",
    )


def audit_pluecker_boundary():
    checks = 0
    for roots in permutations(range(10), 4):
        z0, z1, z2, z3 = roots
        p1 = (z0 - z1) * (z2 - z3)
        p2 = (z0 - z2) * (z1 - z3)
        p3 = (z0 - z3) * (z1 - z2)
        signed = (p1, -p2, p3)
        require(sum(signed) == 0, "oriented Pluecker identity failed")
        data = tuple(leading_residue(entry, P) for entry in signed)
        values = tuple(value for value, _ in data)
        residues = tuple(residue for _, residue in data)
        _, face, sigma = initial_face(values, residues, P)
        require(len(face) >= 2, "Pluecker minimum was unique")
        require(sigma == 0, "Pluecker initial face did not cancel")
        checks += 1
    return checks


def audit_symmetry():
    actions = [matching_action(g) for g in permutations(range(4))]
    distinct_actions = set(actions)
    require(len(distinct_actions) == 6, "S4 matching image is not S3")
    require(actions.count((0, 1, 2)) == 4, "S4 matching kernel is not V4")

    checks = 0
    for values in product(range(3), repeat=3):
        for residues in product(range(1, P), repeat=3):
            minimum, face, sigma = initial_face(values, residues, P)
            for action in distinct_actions:
                moved_values = [0, 0, 0]
                moved_residues = [0, 0, 0]
                for old, new in enumerate(action):
                    moved_values[new] = values[old]
                    moved_residues[new] = residues[old]
                moved = initial_face(tuple(moved_values), tuple(moved_residues), P)
                require(moved[0] == minimum, "matching permutation changed the minimum")
                require(moved[2] == sigma, "matching permutation changed sigma")
                require(
                    tuple(sorted(moved[1]))
                    == tuple(sorted(action[i] for i in face)),
                    "matching permutation changed the face incorrectly",
                )
                checks += 1
    return len(distinct_actions), actions.count((0, 1, 2)), checks


def audit_vertex_unit_gauge():
    checks = 0
    units = range(1, P)
    for edge_residues in product(units, repeat=6):
        before_sigma = sum(matching_products(edge_residues)) % P
        for gauge in product(units, repeat=4):
            moved_edges = tuple(
                edge_residues[i] * gauge[a] * gauge[b] % P
                for i, (a, b) in enumerate(EDGES)
            )
            after_sigma = sum(matching_products(moved_edges)) % P
            scale = 1
            for entry in gauge:
                scale = scale * entry % P
            require(after_sigma == scale * before_sigma % P, "vertex gauge law failed")
            require(
                (after_sigma == 0) == (before_sigma == 0),
                "vertex gauge changed vanishing",
            )
            checks += 1
    return checks


lines = []


def emit(line):
    lines.append(line)


counts, residue_tuple_count = audit_q5_zero_sum_formula()
constructor_checks = audit_zero_sigma_lift_constructor()
ring_record = audit_q5_truncated_triples()
deep_hostile_checks = audit_deep_hostile()
audit_same_clutch_hostile()
pluecker_checks = audit_pluecker_boundary()
matching_image, matching_kernel, symmetry_checks = audit_symmetry()
gauge_checks = audit_vertex_unit_gauge()

emit("FINITE-FIBRE INITIAL-FACE AUGMENTATION / K4 HAFNIAN: VERIFIED-EXACT PASS")
emit("general_law=sigma!=0 iff valuation(sum_i A_i)=minimum valuation")
emit("robust_lift_law=sigma!=0 iff every lift with fixed valuations/residues is nonzero")
emit("zero_sigma_boundary=same first-order data admit exact-zero and higher-valuation nonzero lifts")
emit("actual_nonvanishing=requires the unbounded filtered contraction jet when sigma=0")
emit("zero_sum_formula=((q-1)^r+(-1)^r*(q-1))/q")
emit(
    "F5_zero_sum_counts="
    + ";".join(f"r{arity}:{count}" for arity, count in counts.items())
    + f";exhaustive_tuples:{residue_tuple_count}"
)
emit(f"F5_zero_sigma_lift_constructor_checks={constructor_checks};arities=1..5;values=0..2")
emit("truncated=" + ";".join(f"{key}:{value}" for key, value in ring_record.items()))
emit(f"deep_hostile_checks={deep_hostile_checks};depths=1..12")
emit("clutch_hostile=lambda024/live versus lambda222/zero;clutch=(000,0);residues=(1,1,3)")
emit(f"oriented_Pluecker_sigma_zero_checks={pluecker_checks};q=5")
emit(
    f"S4_matching_image={matching_image};V4_kernel={matching_kernel};"
    f"symmetry_checks={symmetry_checks}"
)
emit(f"vertex_unit_gauge_checks={gauge_checks};q=5")
emit("selector_scope=apply only after THM2290 endpoint-colour selection and pair aggregation")
emit("torsor_scope=initial face is S3-equivariant and V4-blind;no quartic origin is recovered")
emit("scratch_origin_commit=631ae23d4c76c7cbba40a6aef9a1b3c02d878c8d")

semantic = hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()
for line in lines:
    print(line)
print(f"semantic_sha256={semantic}")

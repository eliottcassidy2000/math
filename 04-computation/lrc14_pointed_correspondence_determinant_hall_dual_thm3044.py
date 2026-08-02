#!/usr/bin/env python3
"""Exact referee for THM-3044.

The program uses only the Python standard library and integer/Fraction
arithmetic.  Truth-bearing checks raise explicit exceptions, so ``-O`` runs
the same audit.
"""

from fractions import Fraction
from itertools import permutations, product
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def jdump(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def matpow(a, exponent):
    out = eye(len(a))
    base = a
    e = exponent
    while e:
        if e & 1:
            out = matmul(out, base)
        base = matmul(base, base)
        e //= 2
    return out


def permutation_matrix(perm):
    n = len(perm)
    # Column j is sent to row perm[j].
    return [[int(i == perm[j]) for j in range(n)] for i in range(n)]


def poly_add(a, b):
    out = [Fraction(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return out


def poly_mul(a, b, cutoff=None):
    size = len(a) + len(b) - 1
    if cutoff is not None:
        size = min(size, cutoff + 1)
    out = [Fraction(0)] * size
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < size:
                out[i + j] += x * y
    return out


def permutation_sign(perm):
    inversions = sum(perm[i] > perm[j] for i in range(len(perm)) for j in range(i + 1, len(perm)))
    return -1 if inversions % 2 else 1


def det_i_plus_ta(a):
    n = len(a)
    total = [Fraction(0)]
    for perm in permutations(range(n)):
        term = [Fraction(1)]
        for i in range(n):
            term = poly_mul(term, [Fraction(int(i == perm[i])), Fraction(a[i][perm[i]])])
        term = [permutation_sign(perm) * x for x in term]
        total = poly_add(total, term)
    while len(total) > 1 and total[-1] == 0:
        total.pop()
    return total


def log_series(poly, degree):
    require(poly[0] == 1, "formal logarithm needs constant one")
    x = list(poly[: degree + 1]) + [Fraction(0)] * max(0, degree + 1 - len(poly))
    x[0] -= 1
    out = [Fraction(0)] * (degree + 1)
    power = [Fraction(1)]
    for k in range(1, degree + 1):
        power = poly_mul(power, x, degree)
        sign = 1 if k % 2 else -1
        for j, value in enumerate(power):
            out[j] += sign * value / k
    return out


def compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, length - 1):
            yield (first,) + rest


def bounded_alloc(total, caps, allowed):
    allowed = tuple(allowed)
    if not allowed:
        if total == 0:
            yield (0,) * len(caps)
        return

    out = [0] * len(caps)

    def rec(position, remaining):
        if position == len(allowed):
            if remaining == 0:
                yield tuple(out)
            return
        j = allowed[position]
        for value in range(min(caps[j], remaining) + 1):
            out[j] = value
            yield from rec(position + 1, remaining - value)
        out[j] = 0

    yield from rec(0, total)


def transport_min_cost(p, q, edges, pointed_columns):
    rows = len(p)
    cols = len(q)
    remaining = list(q)
    best = None

    def rec(i, cost):
        nonlocal best
        if best is not None and cost >= best:
            return
        if i == rows:
            if all(value == 0 for value in remaining):
                best = cost
            return
        allowed = [j for j in range(cols) if (i, j) in edges]
        for allocation in bounded_alloc(p[i], remaining, allowed):
            added = allocation[pointed_columns[i]]
            for j, value in enumerate(allocation):
                remaining[j] -= value
            rec(i + 1, cost + added)
            for j, value in enumerate(allocation):
                remaining[j] += value

    rec(0, 0)
    return best


def hall_deficiency(p, q, edges, pointed_columns):
    rows = len(p)
    best = 0
    for mask in range(1 << rows):
        subset = [i for i in range(rows) if mask & (1 << i)]
        neighbours = {
            j
            for i in subset
            for j in range(len(q))
            if (i, j) in edges and j != pointed_columns[i]
        }
        deficit = sum(p[i] for i in subset) - sum(q[j] for j in neighbours)
        best = max(best, deficit)
    return best


def dual_max_n2(p, q, edges, pointed_columns):
    # Total row and column mass agree, so a common potential translation is
    # harmless.  Fix alpha_0=0 and exhaust the integral network potentials.
    best = None
    for alpha_1 in range(-4, 5):
        alpha = (0, alpha_1)
        beta = []
        lawful = True
        for j in range(len(q)):
            incident = [i for i in range(2) if (i, j) in edges]
            if not incident:
                if q[j]:
                    lawful = False
                    break
                beta.append(0)
                continue
            beta.append(min(int(j == pointed_columns[i]) - alpha[i] for i in incident))
        if not lawful:
            continue
        value = sum(alpha[i] * p[i] for i in range(2)) + sum(beta[j] * q[j] for j in range(len(q)))
        if best is None or value > best:
            best = value
    return best


print("THM-3044 POINTED CORRESPONDENCE DETERMINANT AND HALL DUAL")

# Gauge covariance and the formal logarithmic trace law.
base_matrices = [
    [[2, 1, 0], [0, 3, 4], [5, 0, 1]],
    [[0, 2, 1], [3, 0, 1], [1, 4, 0]],
]
perms3 = list(permutations(range(3)))
gauge_cells = 0
log_cells = 0
for c in base_matrices:
    for iota in perms3:
        jmat = permutation_matrix(iota)
        a = matmul(c, jmat)
        poly = det_i_plus_ta(a)
        logs = log_series(poly, 3)
        for m in range(1, 4):
            expected = Fraction((1 if m % 2 else -1) * trace(matpow(a, m)), m)
            require(logs[m] == expected, "log-determinant trace identity failed")
            log_cells += 1
        for ph in perms3:
            pmat = permutation_matrix(ph)
            for qb in perms3:
                qmat = permutation_matrix(qb)
                cprime = matmul(matmul(pmat, c), transpose(qmat))
                jprime = matmul(matmul(qmat, jmat), transpose(pmat))
                aprime = matmul(cprime, jprime)
                require(aprime == matmul(matmul(pmat, a), transpose(pmat)), "pointed gauge covariance failed")
                require(det_i_plus_ta(aprime) == poly, "pointed determinant gauge invariance failed")
                gauge_cells += 1

print(jdump({"pointed_gauge":{"covariance_cells":gauge_cells,"log_trace_cells":log_cells}}))

# Exhaust all two-head/two-root-plus-cemetery support graphs and small margins.
all_edges = [(i, j) for i in range(2) for j in range(3)]
pointed = (0, 1)
feasible = hall_checks = dual_checks = complete_checks = 0
strict_hall_cases = 0
for edge_mask in range(1 << len(all_edges)):
    edges = {edge for k, edge in enumerate(all_edges) if edge_mask & (1 << k)}
    for mass in range(1, 5):
        for p in compositions(mass, 2):
            for q in compositions(mass, 3):
                minimum = transport_min_cost(p, q, edges, pointed)
                if minimum is None:
                    continue
                feasible += 1
                delta = hall_deficiency(p, q, edges, pointed)
                require((minimum == 0) == (delta == 0), "Hall zero gate failed")
                require(minimum >= delta, "Hall lower bound failed")
                hall_checks += 1
                dual = dual_max_n2(p, q, edges, pointed)
                require(dual == minimum, "transport dual mismatch")
                dual_checks += 1
                if minimum > delta:
                    strict_hall_cases += 1
                if len(edges) == len(all_edges):
                    overload = max(0, p[0] + q[0] - mass, p[1] + q[1] - mass)
                    require(minimum == overload, "complete-graph overload formula failed")
                    complete_checks += 1

print(jdump({"transport_exhaustion":{"complete_formula_cases":complete_checks,"dual_cases":dual_checks,"feasible_cases":feasible,"hall_gate_cases":hall_checks,"strict_hall_lower_bound_cases":strict_hall_cases,"universe":"2 heads x (2 roots + cemetery), all supports, total mass 1..4"}}))

# The smallest transparent strict quantitative hostile.
strict_edges = {(0, 0), (0, 1), (1, 1)}
strict_p = (1, 1)
strict_q = (1, 1)
strict_min = transport_min_cost(strict_p, strict_q, strict_edges, (0, 1))
strict_delta = hall_deficiency(strict_p, strict_q, strict_edges, (0, 1))
strict_dual = dual_max_n2(strict_p, strict_q, strict_edges, (0, 1))
require((strict_min, strict_delta, strict_dual) == (2, 1, 2), "strict Hall hostile changed")
print(jdump({"restricted_graph_hostile":{"exact_minimum":strict_min,"max_single_hall_deficiency":strict_delta,"dual_witness_value":strict_dual,"unique_table":[[1,0],[0,1]]}}))

# A robust cascade: every row and column has an off-diagonal incident edge,
# but one Hall cut still misses half of the exact cost.
robust_off = {(0, 1), (0, 2), (1, 4), (2, 3), (3, 0), (4, 0)}
robust_edges = robust_off | {(i, i) for i in range(5)}
robust_p = (1,) * 5
robust_q = (1,) * 5
robust_pointed = tuple(range(5))
robust_min = transport_min_cost(robust_p, robust_q, robust_edges, robust_pointed)
robust_delta = hall_deficiency(robust_p, robust_q, robust_edges, robust_pointed)
robust_u = (-2, -1, -1, 0, 0)
robust_v = (0, 2, 2, 1, 1)
robust_perm = (1, 4, 2, 3, 0)
for i, j in robust_edges:
    require(robust_u[i] + robust_v[j] <= int(i == j), "robust dual witness is infeasible")
robust_dual = sum(robust_u) + sum(robust_v)
require(sorted(robust_perm) == list(range(5)), "robust primal is not a permutation")
require(all((i, robust_perm[i]) in robust_edges for i in range(5)), "robust primal uses a forbidden edge")
require(sum(int(i == robust_perm[i]) for i in range(5)) == 2, "robust primal cost changed")
require((robust_min, robust_delta, robust_dual) == (2, 1, 2), "robust Hall hostile changed")
require(all(any((i, j) in robust_off for j in range(5)) for i in range(5)), "robust row incidence failed")
require(all(any((i, j) in robust_off for i in range(5)) for j in range(5)), "robust column incidence failed")
print(jdump({"robust_multilevel_hostile":{"dual_value":robust_dual,"every_column_offdiagonal":True,"every_row_offdiagonal":True,"exact_minimum":robust_min,"max_single_hall_deficiency":robust_delta,"primal_permutation":list(robust_perm)}}))

# Balanced nonnegative flows: zero first ghost forces a higher cycle ghost.
balanced = balanced_loopless = 0
cycle_lengths = {2: 0, 3: 0}
for entries in product(range(3), repeat=9):
    a = [list(entries[3 * i : 3 * i + 3]) for i in range(3)]
    if not any(entries):
        continue
    if [sum(row) for row in a] != [sum(a[i][j] for i in range(3)) for j in range(3)]:
        continue
    balanced += 1
    if trace(a):
        continue
    balanced_loopless += 1
    witnesses = [m for m in (2, 3) if trace(matpow(a, m)) > 0]
    require(witnesses, "balanced loopless flow has no cycle ghost")
    cycle_lengths[min(witnesses)] += 1

nilpotent_path = [[0, 1, 0], [0, 0, 0], [0, 0, 0]]
require(all(trace(matpow(nilpotent_path, m)) == 0 for m in range(1, 5)), "path hostile should have no ghosts")
print(jdump({"cycle_boundary":{"balanced_nonzero_tables":balanced,"balanced_zero_first_ghost":balanced_loopless,"first_witness_length_counts":cycle_lengths,"unbalanced_path_all_ghosts_zero_through":4}}))

# Independent gauges identify all permutation correspondences; a pointing does not.
identity3 = eye(3)
cycle3 = permutation_matrix((1, 2, 0))
aligned_poly = det_i_plus_ta(identity3)
cycle_poly = det_i_plus_ta(cycle3)
require(aligned_poly == [1, 3, 3, 1], "aligned determinant control failed")
require(cycle_poly == [1, 0, 0, 1], "three-cycle determinant control failed")
require(matmul(identity3, transpose(identity3)) == matmul(cycle3, transpose(cycle3)), "unpointed Gram hostile failed")
same_double_gauge_orbit = any(
    matmul(matmul(permutation_matrix(ph), identity3), transpose(permutation_matrix(qb))) == cycle3
    for ph in perms3
    for qb in perms3
)
require(same_double_gauge_orbit, "aligned/cycle tables left the independent double-gauge orbit")
swap2 = permutation_matrix((1, 0))
require(trace(eye(2)) == 2 and trace(swap2) == 0, "minimal two-root hostile failed")
print(jdump({"unpointed_hostiles":{"aligned_det":"(1+t)^3","cycle_det":"1+t^3","same_double_gauge_orbit":same_double_gauge_orbit,"same_gram":True,"two_root_arrivals":[2,0]}}))

# THM-2549's conditional chart is diagonal once the semantic identification is supplied.
weights = (2, 0, 3, 1)
diagonal = [[int(i == j) * weights[i] for j in range(4)] for i in range(4)]
require(trace(diagonal) == sum(weights), "conditional ancestry diagonal failed")
require(det_i_plus_ta(diagonal) == [Fraction(1), Fraction(6), Fraction(11), Fraction(6)], "conditional determinant control failed")
print(jdump({"conditional_thm2549":{"b_w_equals_head":True,"first_ghost":sum(weights),"packet_mass":sum(weights),"semantic_identification_required":True}}))

print("all_exact_checks=PASS")

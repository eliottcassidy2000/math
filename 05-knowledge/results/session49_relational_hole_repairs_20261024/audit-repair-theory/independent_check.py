#!/usr/bin/env python3
"""Independent Session49 checks. No producer imports, asserts, or nonstdlib modules."""
from collections import Counter, deque
from itertools import combinations, combinations_with_replacement, permutations, product
from pathlib import Path
import argparse
import hashlib
import json
from math import factorial

HERE = Path(__file__).resolve().parent
CHECKS = Counter()


def need(truth, name):
    CHECKS[name] += 1
    if not truth:
        raise RuntimeError(name)


def powerset(xs):
    xs = tuple(xs)
    return [frozenset(c) for k in range(len(xs) + 1) for c in combinations(xs, k)]


def holes():
    cases = 0
    for n in range(1, 6):
        group = set(range(n))
        subsets = powerset(group)
        for C, D in product(subsets, repeat=2):
            for a, b in product(range(n), repeat=2):
                direct = {(x + a) % n for x in C} | {(x + b) % n for x in D}
                branch_holes = {(x + a) % n for x in group - C} & {(x + b) % n for x in group - D}
                need(group - direct == branch_holes, 'arbitrary branch-conditioned hole identity')
                if C == D:
                    H = group - C
                    differences = {(x - y) % n for x in H for y in H}
                    need((direct == group) == ((b - a) % n not in differences), 'equal-bank two-branch difference iff')
                cases += 1
    # Both unconditioned banks are complete. A joint relation removes different rows.
    accepted = [(repair, receiver) for repair, receiver in product(range(2), repeat=2) if repair == receiver]
    actual = {(repair + receiver) % 2 for repair, receiver in accepted}
    need(actual == {0}, 'unconditioned completeness can fail after joint relation')
    joint_cases = 0
    group = set(range(3))
    for branches in product(powerset(group), repeat=3):
        for labels in product(range(3), repeat=3):
            missing = {q for q in group if all((q - c) % 3 not in C for c, C in zip(labels, branches))}
            H = [group - set(C) for C in branches]
            joint = {((a - b) % 3, (a - c) % 3) for a, b, c in product(*H)}
            forbidden = ((labels[1] - labels[0]) % 3, (labels[2] - labels[0]) % 3)
            need(bool(missing) == (forbidden in joint), 'joint multibranch difference obstruction retains one shared hole')
            joint_cases += 1
    H = [{0, 1}, {1, 2}, {0, 2}]
    need(all(a & b for a, b in combinations(H, 2)) and not set.intersection(*H), 'pairwise hole overlaps need not give joint obstruction')
    return {'finite_cyclic_branch_pairs': cases, 'three_branch_joint_difference_cases': joint_cases, 'joint_relation_hostile': {'group_order': 2, 'relation': accepted, 'support': sorted(actual)}}


def cut_number(family, universe):
    for ban in powerset(universe):
        if all(ban & edge for edge in family):
            return len(ban)
    return None  # An empty witness cannot be destroyed by any ban.


def cuts():
    universe = range(3)
    edges = powerset(universe)
    cases = 0
    for family in powerset(range(len(edges))):
        F = [edges[i] for i in family]
        tau = cut_number(F, universe)
        for budget in range(4):
            robust = all(any(not (ban & edge) for edge in F) for ban in powerset(universe) if len(ban) <= budget)
            need(robust == (tau is None or tau > budget), 'universal ban robustness iff transversal exceeds budget')
            cases += 1
    record_to_unit = {'a0': 'a', 'a1': 'a', 'b': 'b', 'c': 'c'}
    admissible = [e for e in powerset(record_to_unit) if len({record_to_unit[v] for v in e}) == len(e)]
    compared = 0
    for indices in powerset(range(len(admissible))):
        F = [admissible[i] for i in indices]
        units = [{record_to_unit[v] for v in e} for e in F]
        tau_record = cut_number(F, record_to_unit)
        tau_unit = cut_number(units, set(record_to_unit.values()))
        need(tau_record is None or (tau_unit is not None and tau_unit <= tau_record), 'resource identification cannot increase transversal number')
        compared += 1
    need(cut_number([{'a0'}, {'a1'}], record_to_unit) == 2, 'two lifts need two full-record deletions')
    need(cut_number([{'a'}, {'a'}], {'a', 'b', 'c'}) == 1, 'same two lifts need one resource ban')
    graphs = 0
    edges = list(product(range(3), repeat=2))
    left_weight, right_weight = (1, 2, 3), (3, 2, 1)
    for selected in powerset(range(9)):
        graph = {edges[i] for i in selected}
        covers = [(U, V) for U in powerset(range(3)) for V in powerset(range(3)) if all(i in U or j in V for i, j in graph)]
        min_cover = min(sum(left_weight[i] for i in U) + sum(right_weight[j] for j in V) for U, V in covers)
        M = sum(left_weight) + sum(right_weight) + 1
        min_cut = min(sum(left_weight[i] for i in range(3) if i not in U) + M * sum(i in U and j not in V for i, j in graph) + sum(right_weight[j] for j in V) for U in powerset(range(3)) for V in powerset(range(3)))
        need(min_cover == min_cut, 'whole-assignment weighted vertex-cover equals explicit network cut')
        for a, b in product(range(4), repeat=2):
            actual = any(len(U) <= a and len(V) <= b for U, V in covers)
            rectangle = any(len(U) <= a and len({j for i, j in graph if i not in U}) <= b for U in powerset(range(3)))
            need(actual == rectangle, 'separate-side deletion budgets retain rectangular quantifier')
        graphs += 1
    return {'robustness_budget_cases': cases, 'capacity_valid_lift_families': compared, 'three_by_three_whole_assignment_graphs': graphs, 'record_vs_resource_hostile': {'record_transversal': 2, 'resource_transversal': 1}, 'empty_family_transversal': 0, 'empty_witness_transversal': 'infinity'}


def assignments(families):
    answer = Counter()
    for choice in product(*families):
        flattened = [resource for edge in choice for resource in edge]
        if len(flattened) == len(set(flattened)):
            answer[frozenset(flattened)] += 1
    return answer


def squarezero(families):
    polynomial = Counter({frozenset(): 1})
    for family in families:
        nxt = Counter()
        for used, multiplicity in polynomial.items():
            for edge in family:
                if not used & edge:
                    nxt[used | edge] += multiplicity
        polynomial = nxt
    return polynomial


def allocations():
    edges = [frozenset(e) for e in combinations(range(4), 2)]
    families = [[edges[i] for i in indices] for indices in powerset(range(6))]
    cases = 0
    for A, B in product(families, repeat=2):
        need(squarezero([A, B]) == assignments([A, B]), 'square-zero coefficient equals full disjoint witness count')
        cases += 1
    A = [frozenset((0, 1)), frozenset((2, 3))]
    B = [frozenset((0, 2)), frozenset((1, 3))]
    need(not assignments([A, B]), 'cross-matching hyperedge hostile has no allocation')
    for rows in ([A], [B], [A, B]):
        need(len(set().union(*(e for row in rows for e in row))) >= 2 * len(rows), 'hyperedge hostile passes even weighted union cardinalities')
    need(all(cut_number(row, range(4)) == 2 for row in (A, B)), 'hyperedge hostile has transversal two in each family')
    singleton = [[frozenset((0,)), frozenset((1,))]] * 2
    need(squarezero(singleton) == {frozenset((0, 1)): 2}, 'nonnegative characteristic-zero multiplicity retains two assignments')
    # Exterior signs give e0*e1 + e1*e0 = 1 - 1 = 0. Modulo 2 also loses 2.
    need(1 + (-1) == 0 and 2 % 2 == 0, 'signed or positive-characteristic coefficients can cancel valid assignments')
    greedy_cases = 0
    for A, B in product(families, repeat=2):
        k = max(map(len, A), default=0)
        tau = cut_number(B, range(4))
        if A and (tau is None or tau > k):
            need(all(any(not a & b for b in B) for a in A), 'transversal depletion condition extends every greedy prefix')
            greedy_cases += 1
    return {'two_rank_two_family_cases': cases, 'greedy_extension_cases': greedy_cases, 'hostile': {'A': [[0, 1], [2, 3]], 'B': [[0, 2], [1, 3]], 'allocation_count': 0}, 'commutative_two_singleton_coefficient': 2, 'exterior_and_F2_coefficient': 0}


def quota_hall():
    resources = range(3)
    domains = powerset(resources)
    cases = 0
    for neighborhoods in product(domains, repeat=3):
        for quotas in product(range(3), repeat=3):
            # Enumerate complete disjoint assignments, independently of flow algorithms.
            menus = [[frozenset(c) for c in combinations(neighbors, quota)] for neighbors, quota in zip(neighborhoods, quotas)]
            actual = bool(assignments(menus))
            cut = all(sum(quotas[i] for i in I) <= len(set().union(*(neighborhoods[i] for i in I))) for I in powerset(range(3)))
            need(actual == cut, 'unit-resource quota allocation iff every demand-side cut passes')
            cases += 1
    return {'three_by_three_incidence_quota_cases': cases, 'quotas_per_row': [0, 1, 2], 'scope': 'fixed boundary, single independent resource per unit of quota'}


def relational_quotas():
    # A relation table is a subset of full compatible pairs, not a simultaneous
    # allocation of records. Enumerate entire table families instead of flows.
    nx, nz = 2, 3
    full = list(product(range(nx), range(nz)))
    cases = 0; infeasible = 0
    all_demands = [(d, e) for d in product(range(nz + 1), repeat=nx) for e in product(range(nx + 1), repeat=nz) if sum(d) == sum(e)]
    for indices in powerset(range(len(full))):
        E = {full[i] for i in indices}
        tables = {}
        for table in powerset(E):
            d = tuple(sum(x == i for x, z in table) for i in range(nx))
            e = tuple(sum(z == j for x, z in table) for j in range(nz))
            tables.setdefault((d, e), []).append(table)
        for d, e in all_demands:
            # Verify nonemptiness first, independently against the analogous cut.
            def deficit(edges, U, V):
                return sum(d[i] for i in U) - sum(e[j] for j in V) - sum(i in U and j not in V for i, j in edges)
            cuts = [(U, V) for U in powerset(range(nx)) for V in powerset(range(nz))]
            feasible = (d, e) in tables
            need(feasible == (max(deficit(E, U, V) for U, V in cuts) <= 0), 'relation-quota feasibility against all complete-table cut deficits')
            if not feasible:
                infeasible += 1
                continue
            for q in range(2):
                W = {(i, j) for i, j in E if (i + j) % 2 == q}
                universal = all(bool(table & W) for table in tables[d, e])
                obstruction = max(deficit(E - W, U, V) for U, V in cuts) > 0
                need(universal == obstruction, 'universal charge retention iff positive relation-table quota deficit')
                cases += 1
    left = (0, 1, 1, 1); right = (0, 0, 1, 1)
    retained = [sorted({(left[i] + right[p[i]]) % 2 for i in range(4)}) for p in permutations(range(4))]
    need(all(q == [0, 1] for q in retained), 'all twenty-four perfect-table relations retain both charges')
    # Every branch in each chosen table has one receiver, so it is not complete.
    need(all(len({(left[i] + right[p[i]]) % 2}) == 1 for p in permutations(range(4)) for i in range(4)), 'quota-universal support needs no individually complete conditioned branch')
    profile_cases = 0; permutation_cases = 0
    for n in range(1, 7):
        for a, b in product(range(n + 1), repeat=2):
            distribution = Counter()
            for p in permutations(range(n)):
                even = sum((i < a) == (p[i] < b) for i in range(n))
                distribution[even] += 1
                permutation_cases += 1
            expected = Counter()
            for x in range(max(0, a + b - n), min(a, b) + 1):
                numerator = factorial(a) * factorial(n - a) * factorial(b) * factorial(n - b)
                denominator = factorial(x) * factorial(a - x) * factorial(b - x) * factorial(n - a - b + x)
                need(numerator % denominator == 0, 'two-charge contingency coefficient is integral')
                expected[n - a - b + 2 * x] += numerator // denominator
            need(distribution == expected and sum(expected.values()) == factorial(n), 'exact C2 matching-table polynomial coefficients')
            need(min(distribution) == abs(n - a - b) and n - max(distribution) == abs(a - b), 'sharp universal two-charge margins')
            need(all(0 < q < n for q in distribution) == (a != b and a + b != n), 'universal two-charge coverage iff unequal and noncomplementary counts')
            profile_cases += 1
    return {'feasible_two_by_three_quota_charge_cases': cases, 'infeasible_quota_families_excluded': infeasible, 'four_by_four_unavoidable_two_charge_tables': len(retained), 'binary_charge_matching_profiles_through_n6': profile_cases, 'binary_charge_full_permutations_counted': permutation_cases, 'quota_semantics': 'number of allowed full pairs, never physical unit capacity'}


def hamming(x, y):
    return sum(a != b for a, b in zip(x, y))


def strict_history_route(states, start, finish, budget):
    resources = {x: frozenset(2 * i + b for i, b in enumerate(x)) for x in states}
    initial = (start, frozenset())
    queue = deque([initial]); seen = {initial}
    while queue:
        x, retired = queue.popleft()
        if x == finish:
            return True
        for y in states:
            if y == x or hamming(x, y) > budget or resources[y] & retired:
                continue
            z = (y, retired | (resources[x] - resources[y]))
            if z not in seen:
                seen.add(z); queue.append(z)
    return False


def endpoint_dag_route(states, start, finish, budget):
    differing = {i for i, (a, b) in enumerate(zip(start, finish)) if a != b}
    admissible = {x: frozenset(i for i in differing if x[i] == finish[i]) for x in states if all(x[i] == start[i] for i in range(len(x)) if i not in differing)}
    reachable = {start}
    for x in sorted(admissible, key=lambda x: len(admissible[x])):
        if x in reachable:
            reachable.update(y for y in admissible if admissible[x] < admissible[y] and len(admissible[y] - admissible[x]) <= budget)
    return finish in reachable


def nonreturn():
    cases = 0
    for n in range(1, 4):
        vertices = list(product(range(2), repeat=n))
        for indices in powerset(range(len(vertices))):
            states = [vertices[i] for i in indices]
            for x, y in combinations(states, 2):
                for budget in range(1, n + 1):
                    need(strict_history_route(states, x, y, budget) == endpoint_dag_route(states, x, y, budget), 'binary endpoint DAG equals global retired-resource search')
                    cases += 1
    hostile = [(0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (1, 1, 0)]
    edges = [(i, j) for i, j in combinations(range(5), 2) if hamming(hostile[i], hostile[j]) == 1]
    need(edges == [(0, 1), (1, 2), (2, 3), (3, 4)], 'native hostile table is precisely five-state one-slot path')
    need(not strict_history_route(hostile, hostile[0], hostile[-1], 1), 'native hostile cannot avoid resource reacquisition at one slot')
    need(strict_history_route(hostile, hostile[0], hostile[-1], 2), 'native hostile direct two-slot strict route')
    return {'binary_accepted_table_endpoint_budget_cases': cases, 'hostile_table': hostile, 'one_slot_edges': edges, 'strict_threshold': {'slots': 2, 'owner_clocks': 4}}


def antimatroids():
    n = 3; ground = frozenset(range(n)); subsets = powerset(ground)
    count = 0; endpoints = 0
    slot = [1, 2, 3]; clock = [7, 2, 5]
    for indices in powerset(range(len(subsets))):
        family = {subsets[i] for i in indices}
        if frozenset() not in family or ground not in family:
            continue
        if any(A | B not in family for A, B in product(family, repeat=2)):
            continue
        if any(A and not any(A - {i} in family for i in A) for A in family):
            continue
        count += 1
        for A, B in product(family, repeat=2):
            difference = A ^ B
            budget_s = max((slot[i] for i in difference), default=0)
            budget_c = max((clock[i] for i in difference), default=0)
            # Only exact endpoint-aligned changes, and arbitrary retained subsets.
            allowed = {C for C in family if A & B <= C <= A | B}
            reached = {A}
            queue = [A]
            while queue:
                C = queue.pop()
                for D in allowed - reached:
                    changes = C ^ D
                    if len(D ^ B) >= len(C ^ B) or not changes <= C ^ B:
                        continue
                    if sum(slot[i] for i in changes) <= budget_s and sum(clock[i] for i in changes) <= budget_c:
                        reached.add(D); queue.append(D)
            need(B in reached, 'antimatroid exact endpoint simultaneous bottlenecks have monotone route')
            endpoints += 1
    return {'full_ground_three_element_antimatroids': count, 'ordered_endpoint_cases': endpoints, 'block_slot_weights': slot, 'block_clock_weights': clock}


def paid_boundaries():
    # Coordinate minima must not be combined across distinct endpoints.
    options = [(1, 10), (2, 2)]
    need(not any(s <= 1 and c <= 2 for s, c in options), 'separate slot-clock minima are not one event')
    # Any additional directed transition predicate needs directed reachability.
    directed = {('a', 'b')}
    need(('b', 'a') not in directed, 'weak component contraction loses directed reachability')
    # Reusing x in two endpoint alphabets is legal; x*x in one alphabet vanishes.
    before = [('before', 0)]; after = [('after', 0)]
    need(len(set(before + after)) == 2 and len({0} | {0}) == 1, 'two capacity alphabets preserve continuous occupancy')
    return {'incompatible_pareto_minima': options, 'directed_hostile': 'a to b only; weak contraction invents b to a', 'same_resource_before_after': 'legal in two alphabets, annihilated in one'}


def davenport():
    # Auxiliary F2 group algebra is appropriate for zero-sum obstruction,
    # not for the cancellation-free occupancy/existence compiler above.
    def plus(a, b):
        return ((a % 128 + b % 128) % 128) + 128 * ((a // 128) ^ (b // 128))
    def multiply_one_plus(poly, g):
        return poly ^ {plus(x, g) for x in poly}
    generators = [1, 128, 256, 512, 1024]
    orders = [128, 2, 2, 2, 2]
    for generator, order in zip(generators, orders):
        poly = {0}
        for _ in range(order):
            poly = multiply_one_plus(poly, generator)
        need(not poly, 'augmentation generators have exact advertised nilpotence upper powers')
    poly = {0}
    for generator, order in zip(generators, orders):
        for _ in range(order - 1):
            poly = multiply_one_plus(poly, generator)
    need(poly == set(range(2048)), 'degree131 sharp augmentation product is full group sum')
    for g in range(2048):
        need(not multiply_one_plus(poly, g), 'every further augmentation factor kills sharp degree131 product')
    need(sum(order - 1 for order in orders) == 131, 'truncated polynomial maximum nonzero total degree131')
    need(all((cyclic + 128 * mask != 0) for cyclic, mask in product(range(128), range(16)) if cyclic or mask), 'sharp length131 generator multiset has no nonempty zero subset')
    zero_patterns = []
    for cyclic, mask, final in product(range(128), range(16), range(2)):
        value = plus(cyclic + 128 * mask, 1 + 128 * 15 if final else 0)
        if value == 0:
            zero_patterns.append((cyclic, mask, final))
    need(zero_patterns == [(0, 0, 0), (127, 15, 1)], 'length132 generator-plus-completion sequence is minimally zero-sum')
    finite_cases = 0
    for moduli in ((2,), (4,), (2, 2), (4, 2), (2, 2, 2)):
        group = list(product(*(range(m) for m in moduli)))
        bound = 1 + sum(m - 1 for m in moduli)
        for sequence in combinations_with_replacement(group, bound):
            has_zero = any(all(sum(sequence[i][j] for i in I) % moduli[j] == 0 for j in range(len(moduli))) for I in powerset(range(bound)) if I)
            need(has_zero, 'small two-group Davenport bound by literal subset sums')
            finite_cases += 1
    # Two independent population cohorts, binary endpoints and C3 charges.
    # Enumerate neutral endpoint differences, pair plus/minus inside cohorts,
    # then execute an inclusion-minimal zero-sum partition of the pair sequence.
    pair_cases = 0
    for increments in product(range(3), repeat=4):
        if sum(increments) % 3:
            continue
        # Each pair consists of one +population and one -population switch.
        pairs = [(0, 1), (2, 3), (4, 5), (6, 7)]
        cohort = [0, 0, 1, 1, 0, 0, 1, 1]
        old = [0, 1, 0, 1, 0, 1, 0, 1]
        new = [1 - v for v in old]
        state = old[:]; left = set(range(4)); retired = set()
        while left:
            neutral = next(I for I in powerset(sorted(left)) if I and sum(increments[i] for i in I) % 3 == 0)
            need(len(neutral) <= 3, 'minimal cyclic-three paired trade has at most D pairs')
            changed = [slot for i in neutral for slot in pairs[i]]
            for slot in changed:
                resource = 2 * slot + state[slot]
                retired.add(resource); state[slot] = new[slot]
                need(2 * slot + state[slot] not in retired, 'disjoint endpoint pair trades do not reacquire resources')
            need(all(sum(state[i] for i in range(8) if cohort[i] == c) == sum(old[i] for i in range(8) if cohort[i] == c) for c in (0, 1)), 'paired trade preserves each population margin separately')
            left -= neutral
        need(state == new, 'paired partition restores exact requested endpoint')
        pair_cases += 1
    return {'guard_group': 'C128 x C2^4', 'group_order': 2048, 'maximum_augmentation_degree': 131, 'exact_Davenport_constant': 132, 'conditional_endpoint_packet_bound_slots': 264, 'abstract_sharpness': {'separate_two_slot_population_classes': 132, 'accepted_full_states': 2, 'minimal_zero_sum_subset_patterns': zero_patterns, 'not_sharpness_for_only_core_helper_margins': True, 'not_a_native_electrical_bank': True}, 'sharp_degree131_coefficient_support': len(poly), 'small_two_group_multisets_at_bound': finite_cases, 'two_cohort_cyclic_three_pair_sequences': pair_cases, 'scope': 'complete additive guard; matched opposite population increments; no unverified relational prefix constraints'}


P0 = (3, 0, 6, 7, 1, 5, 4, 2)
P1 = (3, 1, 6, 7, 0, 5, 4, 2)
G = (4, 0, 2, 3, 1, 5, 6, 7)


def seed_arcs():
    arcs = set()
    for i, (a, b) in enumerate(combinations(range(8), 2)):
        arcs.add((a, b) if 29558785 & (1 << (27 - i)) else (b, a))
    return arcs


def source_arcs(r):
    seed = seed_arcs(); q = 8 * r
    arcs = set()
    for a, b in combinations(range(q + 3), 2):
        if b < q and a // 8 == b // 8:
            forward = (a % 8, b % 8) in seed
        elif a >= q:
            forward = (b - a) % 3 == 1
        else:
            forward = True
        arcs.add((a, b) if forward else (b, a))
    return arcs


def native_path(r, key, phase):
    path = tuple(8 * block + u for block in range(r) for u in (P1 if key & (1 << block) else P0))
    return path + tuple(8 * r + (phase + j) % 3 for j in range(3))


def partner(arcs, path):
    adjacent = {frozenset((a, b)) for a, b in zip(path, path[1:])}
    return {(a, b) if frozenset((a, b)) in adjacent else (b, a) for a, b in arcs}


def word(arcs, path):
    return tuple(int((path[i], path[j]) in arcs) for i, j in combinations(range(len(path)), 2))


def reservoir(arcs, path):
    target = partner(arcs, path)
    return frozenset((word(arcs, path), word(target, path), word({(b, a) for a, b in arcs}, path[::-1]), word({(b, a) for a, b in target}, path[::-1])))


def native():
    seed = seed_arcs()
    automorphisms = [p for p in permutations(range(8)) if {(p[a], p[b]) for a, b in seed} == seed]
    need(automorphisms == [tuple(range(8))], 'seed rigidity by all literal permutations')
    rows = []
    for r in range(1, 5):
        q = 8 * r; S = source_arcs(r)
        records = {(key, phase): native_path(r, key, phase) for key in range(1 << r) for phase in range(3)}
        tiles = {}
        for (key, phase), path in records.items():
            need(set(path) == set(range(q + 3)) and all(edge in S for edge in zip(path, path[1:])), 'actual native Hamiltonian path')
            tiles[key, phase] = reservoir(S, path)
            need(len(tiles[key, phase]) == 4, 'four distinct inverse-converse words')
            target = partner(S, path)
            for d in (1, 2):
                rotation = {u: u if u < q else q + (u - q + d) % 3 for u in range(q + 3)}
                changed = records[key, (phase + d) % 3]
                need(tuple(rotation[u] for u in path) == changed, 'source kernel lift literal phase')
                need({(rotation[a], rotation[b]) for a, b in S} == S, 'source kernel automorphism')
                need({(rotation[a], rotation[b]) for a, b in target} == partner(S, changed), 'actual partner kernel transport')
                need(hamming(path, changed) == 3, 'phase change exact three-clock cost')
            need(path[q] - q == phase, 'canonical suffix phase subtracts raw owner offset')
            need((path[::-1][0] - q) % 3 == (phase + 2) % 3, 'converse-head chart adds two')
            need(path[::-1][2] == path[q], 'converse reflected clock two preserves phase')
            need(partner({(b, a) for a, b in S}, path[::-1]) == {(b, a) for a, b in target}, 'literal converse partner identity')
            if key % 2 == 0:
                other = records[key + 1, phase]
                need(hamming(path, other) == 2, 'paired key event exact two-clock cost')
                perm = tuple(G[u] if u < 8 else u for u in range(q + 3))
                need({(perm[a], perm[b]) for a, b in target} == partner(S, other), 'paired key endpoint-marked partner transport')
            for block in range(r):
                if key & (1 << block):
                    continue
                perm = {u: 8 * block + G[u % 8] if 8 * block <= u < 8 * block + 8 else u for u in range(q + 3)}
                changed = records[key | (1 << block), phase]
                need({(perm[a], perm[b]) for a, b in target} == partner(S, changed), 'every prefix-bit change authenticates actual target class')
        need(all(tiles[key, phase] == tiles[key, 0] for key, phase in records), 'all three lifts use one normalized unit')
        need(len(set().union(*(tiles[key, 0] for key in range(1 << r)))) == 4 * (1 << r), 'distinct keys literal four-word reservoirs are disjoint')
        for (b, k), (c, ell) in combinations(records, 2):
            need(hamming(records[b, k], records[c, ell]) == 2 * (b ^ c).bit_count() + 3 * (k != ell), 'all native record pairs have exact key-phase cost formula')
        for key in range(0, 1 << r, 2):
            fixed = records[key, 0]
            completions = set()
            for first_free in permutations((0, 1)):
                for tail in permutations(range(q, q + 3)):
                    candidate = list(fixed)
                    candidate[1], candidate[4] = first_free
                    candidate[q:] = tail
                    candidate = tuple(candidate)
                    if all(edge in S for edge in zip(candidate, candidate[1:])):
                        completions.add(candidate)
            offered = {records[key + bit, phase] for bit in range(2) for phase in range(3)}
            need(completions == offered, 'complete six-record paired-key menu under literal fixed-clock pins')
            phase_zero = {p for p in completions if p[q] == q}
            need(phase_zero == {records[key, 0], records[key + 1, 0]}, 'literal suffix-owner pin leaves exactly two binary records')
        rows.append({'r': r, 'order': q + 3, 'records': len(records), 'resource_units': 1 << r})
    # Explicit raw-coordinate shift hostile, not hidden by r divisible by 3.
    need(native_path(1, 0, 0)[8] % 3 == 2 and native_path(1, 0, 1)[8] % 3 == 0, 'raw-owner pin zero is not canonical phase zero')
    # Reconstruct all six actual full boundary records and count all changed owners.
    r = 3; S = source_arcs(r)
    boundary = {}
    for mode, keys in enumerate(((0, 1, 4, 5), (2, 3, 4, 5))):
        for phase in range(3):
            boundary[mode, phase] = tuple(native_path(r, key, phase) for key in keys) + (native_path(r, 6, mode),)
    costs = Counter()
    adjacency = {7: {}, 12: {}}
    for a, b in combinations(boundary, 2):
        cost = sum(hamming(x, y) for x, y in zip(boundary[a], boundary[b]))
        costs[cost] += 1
        expected = 12 if a[0] == b[0] else 7 if a[1] == b[1] else 19
        need(cost == expected, 'actual phase/key boundary literal owner-clock cost')
    def sizes(budget):
        unseen = set(boundary); out = []
        while unseen:
            queue = [unseen.pop()]; count = 0
            while queue:
                x = queue.pop(); count += 1
                adjacent = {y for y in unseen if sum(hamming(a, b) for a, b in zip(boundary[x], boundary[y])) <= budget}
                unseen -= adjacent; queue.extend(adjacent)
            out.append(count)
        return sorted(out)
    need(sizes(7) == [2, 2, 2] and sizes(12) == [6], 'paid boundary requires budget-component retention')
    return {'literal_native_scales': rows, 'phase_zero_raw_owner_at_r1': 8, 'phase_one_raw_owner_residue_at_r1': 0, 'six_state_boundary': {'pair_cost_multiplicities': dict(sorted(costs.items())), 'clock7_components': sizes(7), 'clock12_components': sizes(12)}}


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', type=Path); args = parser.parse_args()
    receipt = {'schema': 'session49-independent-relational-repair-audit-v1', 'holes': holes(), 'cuts': cuts(), 'allocations': allocations(), 'quota_hall': quota_hall(), 'relational_quotas': relational_quotas(), 'nonreturn': nonreturn(), 'antimatroids': antimatroids(), 'paid_boundaries': paid_boundaries(), 'davenport': davenport(), 'native': native(), 'checks': dict(sorted(CHECKS.items()))}
    data = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(data)
    print('PASS independent relational repair audit; ' + str(sum(CHECKS.values())) + ' checks')


if __name__ == '__main__':
    main()

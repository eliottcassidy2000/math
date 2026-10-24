#!/usr/bin/env python3
"""Exact abstract controls. Explicit checks remain active under python -O."""
import argparse
import hashlib
import itertools as it
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
COUNTS = {}


def need(condition, name):
    COUNTS[name] = COUNTS.get(name, 0) + 1
    if not condition:
        raise RuntimeError(name)


def submasks(mask):
    sub = mask
    while True:
        yield sub
        if not sub:
            break
        sub = (sub - 1) & mask


def members(mask, n):
    return tuple(i for i in range(n) if mask >> i & 1)


def subset(a, b):
    return not (a & ~b)


def reachability(n, edges):
    r = [(1 << i) for i in range(n)]
    for i, j in edges:
        r[i] |= 1 << j
    for k in range(n):
        for i in range(n):
            if r[i] >> k & 1:
                r[i] |= r[k]
    return r


def components(n, edges):
    reach = reachability(n, edges)
    unseen = set(range(n))
    comps = []
    while unseen:
        i = min(unseen)
        comp = tuple(j for j in sorted(unseen) if reach[i] >> j & 1 and reach[j] >> i & 1)
        comps.append(comp)
        unseen.difference_update(comp)
    return comps


def is_ideal(mask, edges):
    return all(not (mask >> i & 1) or mask >> j & 1 for i, j in edges)


def charge(mask, q):
    total = 0
    for i, value in enumerate(q):
        if mask >> i & 1:
            total ^= value
    return total


def grade(mask, d):
    return sum(value for i, value in enumerate(d) if mask >> i & 1)


def cost(mask, s, w):
    return (sum(a for i, a in enumerate(s) if mask >> i & 1),
            sum(a for i, a in enumerate(w) if mask >> i & 1))


def pareto(pairs):
    return sorted(p for p in set(pairs) if not any(q != p and q[0] <= p[0] and q[1] <= p[1] for q in pairs))


def bottlenecks(states, s, w, covers_only=False):
    states = sorted(states, key=lambda z: (z.bit_count(), z))
    dp = {0: [(0, 0)]}
    for target in states[1:]:
        pairs = []
        for source in states:
            if source == target or not subset(source, target):
                continue
            if covers_only and any(source != mid != target and subset(source, mid) and subset(mid, target) for mid in states):
                continue
            c = cost(target ^ source, s, w)
            pairs.extend((max(p[0], c[0]), max(p[1], c[1])) for p in dp[source])
        dp[target] = pareto(pairs)
    return dp[states[-1]]


def partitions(n):
    """Every set partition once, represented by bitmask blocks."""
    def rec(i, parts):
        if i == n:
            yield tuple(parts)
            return
        for j in range(len(parts)):
            changed = parts.copy()
            changed[j] |= 1 << i
            yield from rec(i + 1, changed)
        yield from rec(i + 1, parts + [1 << i])
    return rec(0, [])


def contracted_edges(parts, edges):
    owner = {i: k for k, p in enumerate(parts) for i in members(p, p.bit_length())}
    return sorted({(owner[i], owner[j]) for i, j in edges if owner[i] != owner[j]})


def requirement_order(k, edges):
    done = set()
    order = []
    while len(done) < k:
        available = [i for i in range(k) if i not in done and all(j in done for a, j in edges if a == i)]
        if not available:
            return None
        chosen = min(available)
        done.add(chosen)
        order.append(chosen)
    return order


def compile_interval(n, clauses, x, y):
    active = members(x ^ y, n)
    local = {i: j for j, i in enumerate(active)}
    edges = []
    for clause in clauses:
        residual = []
        tautology = False
        for i, sign in clause:
            if i not in local:
                if (x >> i & 1) == sign:
                    tautology = True
                    break
            else:
                residual.append((local[i], sign ^ (x >> i & 1)))
        if tautology:
            continue
        need(len(residual) == 2, 'residual nontrivial clause has two differing literals')
        need(residual[0][1] != residual[1][1], 'residual nontrivial clause is implication')
        source = next(i for i, sign in residual if sign == 0)
        target = next(i for i, sign in residual if sign == 1)
        edges.append((source, target))
    return active, sorted(set(edges))


def cnf_families(n):
    clauses = [((i, sign),) for i in range(n) for sign in (0, 1)]
    clauses += [((i, a), (j, b)) for i in range(n) for j in range(i + 1, n) for a in (0, 1) for b in (0, 1)]
    clause_masks = [sum(1 << x for x in range(1 << n) if any((x >> i & 1) == sign for i, sign in clause)) for clause in clauses]
    full = (1 << (1 << n)) - 1
    families = {full}
    for cm in clause_masks:
        families |= {f & cm for f in families}
    return clauses, clause_masks, sorted(families)


def classification_checks():
    rows = []
    for n in range(1, 5):
        clauses, masks, families = cnf_families(n)
        endpoint_pairs = 0
        for family in families:
            formula = [c for c, cm in zip(clauses, masks) if family & ~cm == 0]
            reconstructed = sum(1 << x for x in range(1 << n) if all(any((x >> i & 1) == sign for i, sign in c) for c in formula))
            need(reconstructed == family, 'canonical clauses reconstruct every 2-CNF family')
            accepted = members(family, 1 << n)
            for x, y in it.product(accepted, repeat=2):
                endpoint_pairs += 1
                active, edges = compile_interval(n, formula, x, y)
                comps = components(len(active), edges)
                local_expected = {z for z in range(1 << len(active)) if family >> (x ^ sum(1 << active[i] for i in members(z, len(active)))) & 1}
                direct = {z for z in range(1 << len(active)) if is_ideal(z, edges)}
                need(local_expected == direct, 'all endpoint interval states equal implication closed sets')
                owner = {v: c for c, comp in enumerate(comps) for v in comp}
                quotient = sorted({(owner[i], owner[j]) for i, j in edges if owner[i] != owner[j]})
                from_ideals = {sum(1 << v for c in members(bits, len(comps)) for v in comps[c])
                               for bits in range(1 << len(comps)) if is_ideal(bits, quotient)}
                need(from_ideals == local_expected, 'SCC ideals reconstruct every literal endpoint state')
                need(requirement_order(len(comps), quotient) is not None, 'SCC condensation is acyclic')
                for comp in comps:
                    need(all(all((z >> comp[0] & 1) == (z >> v & 1) for v in comp) for z in local_expected), 'SCC member slots synchronize')
        rows.append({'slots': n, 'all_2cnf_solution_families_including_empty': len(families), 'accepted_ordered_endpoint_pairs': endpoint_pairs})
    return rows


def poset_packet_checks():
    rows = []
    for n in range(1, 5):
        possible = [(i, j) for i in range(n) for j in range(i)]
        ps = list(partitions(n))
        s = tuple(1 + i % 2 for i in range(n))
        w = tuple(1 + (2 * i + 1) % 5 for i in range(n))
        instances = 0
        neutral_partitions = 0
        grades = {tuple(0 for _ in range(n))}
        if n % 2 == 0:
            grades.add(tuple(1 if i % 2 == 0 else -1 for i in range(n)))
        if n > 1:
            grades.add((1,) + (0,) * (n - 2) + (-1,))
        for edge_bits in range(1 << len(possible)):
            edges = [edge for i, edge in enumerate(possible) if edge_bits >> i & 1]
            ideals = [z for z in range(1 << n) if is_ideal(z, edges)]
            for q in it.product(range(4), repeat=n):
                if charge((1 << n) - 1, q):
                    continue
                for d in sorted(grades):
                    instances += 1
                    states = [z for z in ideals if charge(z, q) == 0 and grade(z, d) == 0]
                    direct = bottlenecks(states, s, w)
                    covers = bottlenecks(states, s, w, True)
                    need(direct == covers, 'cover refinement preserves exact joint bottleneck frontier')
                    candidates = []
                    merged_candidates = []
                    for parts in ps:
                        if any(charge(p, q) or grade(p, d) for p in parts):
                            continue
                        neutral_partitions += 1
                        deps = contracted_edges(parts, edges)
                        order = requirement_order(len(parts), deps)
                        if order is not None:
                            running = 0
                            for i in order:
                                running |= parts[i]
                                need(running in states, 'acyclic neutral partition has legal accepted prefixes')
                            candidates.append(tuple(max(cost(p, s, w)[axis] for p in parts) for axis in (0, 1)))
                        ccs = components(len(parts), deps)
                        merged = [sum(parts[j] for j in comp) for comp in ccs]
                        merged_deps = contracted_edges(merged, edges)
                        merged_order = requirement_order(len(merged), merged_deps)
                        need(merged_order is not None, 'neutral packet SCC merging removes every dependency cycle')
                        running = 0
                        for i in merged_order:
                            running |= merged[i]
                            need(running in states, 'cycle merged neutral packets give accepted strict route')
                        merged_candidates.append(tuple(max(cost(p, s, w)[axis] for p in merged) for axis in (0, 1)))
                    need(pareto(candidates) == direct, 'neutral acyclic packetizations equal exact route Pareto frontier')
                    need(pareto(merged_candidates) == direct, 'cycle merging over all neutral partitions gives exact frontier')
                    accessible = all(any((z ^ (1 << i)) in states for i in members(z, n)) for z in states if z)
                    all_neutral = all(q[i] == 0 and d[i] == 0 for i in range(n))
                    need(accessible == all_neutral, 'neutral ideal accessibility iff each SCC defect vanishes')
                    for a, b in it.product(states, repeat=2):
                        if a == b or not subset(a, b):
                            continue
                        diff = b ^ a
                        relative = [k for k in submasks(diff) if all(not (k >> i & 1) or not (diff >> j & 1) or k >> j & 1 for i, j in edges)]
                        mid_expected = {a | k for k in relative if charge(k, q) == 0 and grade(k, d) == 0}
                        mid_actual = {z for z in states if subset(a, z) and subset(z, b)}
                        need(mid_expected == mid_actual, 'relative zero-defect ideals equal accepted intermediate interval')
        rows.append({'SCC_vertices': n, 'forward_DAG_presentations': 1 << len(possible), 'charge_and_integer_margin_instances': instances, 'neutral_partitions_tested': neutral_partitions})
    return rows


def path_charges(n):
    prefixes = [0] + [1 if k % 2 else 2 for k in range(1, n)] + [0]
    return [prefixes[i] ^ prefixes[i + 1] for i in range(n)], prefixes


def hostile_checks():
    certificates = []
    for n in range(2, 65):
        q, prefix = path_charges(n)
        edges = [(i, i - 1) for i in range(1, n)]
        need(all(q), 'chain every physical slot charge increment nonzero')
        need(len(components(n, edges)) == n, 'chain original SCC size one')
        need([k for k in range(n + 1) if charge((1 << k) - 1, q) == 0] == [0, n], 'chain only empty full prefixes are neutral')
        if n <= 12:
            states = [z for z in range(1 << n) if is_ideal(z, edges)]
            need(states == [(1 << k) - 1 for k in range(n + 1)], 'literal chain states exhaust exactly prefixes')
        if n % 2 == 0:
            x = sum(1 << i for i in range(n) if i % 2)
            y = ((1 << n) - 1) ^ x
            need(x.bit_count() == y.bit_count() == n // 2, 'even chain has exact equal population endpoints')
            balanced_prefixes = [k for k in range(n + 1) if (x ^ ((1 << k) - 1)).bit_count() == n // 2]
            need(balanced_prefixes == list(range(0, n + 1, 2)), 'even chain exact integer population not reduced modulo group')
        if n in (4, 6, 10, 64):
            certificates.append({'physical_slots': n, 'dependencies': [list(e) for e in edges], 'charge_increments_C2_squared': q, 'proper_prefix_charges': prefix[1:-1], 'exact_population': n // 2, 'accepted_switch_prefix_lengths': [0, n], 'all_increment_nonzero': True, 'minimum_slot_event': n, 'clock_costs': list(range(1, n + 1)), 'minimum_clock_event': n * (n + 1) // 2, 'scope': 'ABSTRACT'})
    # Four-slot contraction is an actual two-cycle, despite two neutral pairs.
    q, _ = path_charges(4)
    parts = [0b1001, 0b0110]
    edges = [(1, 0), (2, 1), (3, 2)]
    need(all(charge(p, q) == 0 for p in parts), 'four slot chosen packets both neutral')
    need(contracted_edges(parts, edges) == [(0, 1), (1, 0)], 'four slot neutral contraction creates a directed two cycle')
    need(requirement_order(2, contracted_edges(parts, edges)) is None, 'four slot small neutral packets have no legal ordering')
    for n in range(2, 11):
        q, _ = path_charges(n)
        physical = 2 * n
        edges = [(2 * i, 2 * i + 1) for i in range(n)] + [(2 * i + 1, 2 * i) for i in range(n)] + [(2 * i, 2 * i - 2) for i in range(1, n)]
        charges = []
        for d in q:
            e, f = [a for a in (1, 2, 3) if a != d]
            need(e ^ f == d and e != 0 and f != 0, 'paired lift has nonzero physical increments summing to packet defect')
            charges.extend((e, f))
        comps = components(physical, edges)
        need(all(len(c) == 2 for c in comps), 'paired lift exact SCC size two')
        need(len({tuple(sorted(e)) for e in edges}) == physical - 1, 'paired lift relation graph has tree edge count')
        full = (1 << physical) - 1
        source = sum(1 << (2 * i + 1) for i in range(n))
        target = full ^ source
        need(source.bit_count() == target.bit_count() == n, 'paired lift exact total population')
        states = [sum((3 << (2 * i)) for i in range(k)) for k in range(n + 1)]
        for z in states:
            need(is_ideal(z, edges), 'paired literal prefix satisfies binary relations')
            need((source ^ z).bit_count() == n, 'paired all relation prefixes preserve population')
        need([z for z in states if charge(z, charges) == 0] == [0, full], 'paired lift only endpoints charge neutral')
        if n <= 6:
            literal_states = [z for z in range(1 << physical) if is_ideal(z, edges)]
            need(literal_states == states, 'paired full binary cube enumeration reconstructs all legal states')
    # D(C2^2)=3 and explicit nonzero two-term zero-sumfree witness.
    for q in it.product(range(4), repeat=3):
        need(any(charge(z, q) == 0 for z in range(1, 8)), 'every three term C2 squared sequence has a zero sum subset')
    need(all(charge(z, (1, 2)) != 0 for z in range(1, 4)), 'C2 squared two distinct nonzero terms prove D lower bound')
    return {'chain_lengths': [2, 64], 'literal_chain_cube_exhaustion_through': 12, 'paired_logical_lengths': [2, 10], 'paired_full_cube_exhaustion_through': 6, 'D_C2_squared': 3, 'four_slot_packet_dependency_edges': [[0, 1], [1, 0]], 'certificates': certificates}


def negative_controls():
    # A ternary table containing 000 and 111 need not be an ideal lattice.
    # This specific table is not closed under intersection.
    family = {0, 3, 5, 7}
    need((3 & 5) not in family, 'arbitrary ternary table is not automatically an ideal lattice')
    # A charge fibre inside a Boolean cube need not be union closed.
    q = (1, 1, 1, 1)
    states = {z for z in range(16) if charge(z, q) == 0}
    need(3 in states and 5 in states and (3 | 5) not in states, 'neutral ideal subfamily need not inherit relation-only union closure')
    # Zero-charge subset can fail ideal closure and be an illegal first event.
    need(charge(0b0110, (1, 3, 3, 1)) == 0 and not is_ideal(0b0110, [(1, 0), (2, 1), (3, 2)]), 'neutral subset alone is not a legal first event')
    return {'non_binary_relation_table': sorted(family), 'neutral_union_counterexample': [3, 5, 7], 'illegal_zero_sum_first_packet': [2, 3]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'schema': 'session51-implication-intervals-v1', 'scope': 'ABSTRACT: no native electrical charge authentication', 'classification': classification_checks(), 'packet_and_cover_controls': poset_packet_checks(), 'hostiles': hostile_checks(), 'negative_controls': negative_controls(), 'checks': dict(sorted(COUNTS.items()))}
    out = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(out)
    print(out, end='')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Independent bounded controls, using affine orbit-cube coordinates.

No producer imports, repository atlas, graph library, or full S10 census.
All checks remain active under python -O.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
import json


def check(condition, explanation):
    if not condition:
        raise RuntimeError(explanation)


def compose(g, h):
    return tuple(g[h[i]] for i in range(len(g)))


def inverse(g):
    return tuple(g.index(i) for i in range(len(g)))


def generated(generators, n):
    identity = tuple(range(n))
    answer, fringe = {identity}, [identity]
    while fringe:
        h = fringe.pop()
        for g in generators:
            gh = compose(g, h)
            if gh not in answer:
                answer.add(gh)
                fringe.append(gh)
    return answer


def conjugate(H, g):
    gi = inverse(g)
    return frozenset(compose(compose(g, h), gi) for h in H)


class Cube:
    def __init__(self, n, K):
        self.n, self.K = n, K
        unseen = set(combinations(range(n), 2))
        self.orbits = []
        self.base_arcs = []
        while unseen:
            i, j = min(unseen)
            arcs = frozenset((g[i], g[j]) for g in K)
            edges = frozenset(tuple(sorted(e)) for e in arcs)
            check(len(edges) % 2 == 1, 'odd pair orbit')
            check(len(edges) == len(arcs), 'no edge inversion by K')
            self.orbits.append(edges)
            self.base_arcs.append(arcs)
            unseen -= edges
        self.d = len(self.orbits)
        self.edge_orbit = {e: i for i, orb in enumerate(self.orbits) for e in orb}

    def action(self, g):
        dest, flips = [], 0
        for arcs in self.base_arcs:
            i, j = min(arcs)
            arc = (g[i], g[j])
            target = self.edge_orbit[tuple(sorted(arc))]
            dest.append(target)
            if arc not in self.base_arcs[target]:
                flips |= 1 << target
        check(len(set(dest)) == self.d, 'normalizer permutes pair orbits')
        return tuple(dest), flips

    def materialize(self, x, tournament):
        if tournament:
            return frozenset(edge if not x >> i & 1 else edge[::-1]
                             for i, arcs in enumerate(self.base_arcs) for edge in arcs)
        return frozenset(edge for i, edges in enumerate(self.orbits) if x >> i & 1 for edge in edges)


def permute_bits(x, dest):
    return sum(((x >> i) & 1) << j for i, j in enumerate(dest))


def act(x, action, tournament):
    dest, flips = action
    return permute_bits(x, dest) ^ (flips if tournament else 0)


def signed_image(x, action):
    dest, flips = action
    y = permute_bits(x, dest)
    return y, (-1 if (y & flips).bit_count() % 2 else 1)


def relative_records(cube, normalizer, tournament):
    actions = [cube.action(g) for g in normalizer]
    unseen, records = set(range(1 << cube.d)), []
    while unseen:
        x = min(unseen)
        images = [act(x, a, tournament) for a in actions]
        unseen.difference_update(images)
        signed = {}
        consistent = True
        for y, a in zip(images, actions):
            value = 1 if tournament else signed_image(x, a)[1]
            if y in signed and signed[y] != value:
                consistent = False
            signed[y] = value
        if consistent:
            stabilizer = frozenset(g for g, y in zip(normalizer, images) if y == x)
            records.append((x, signed, stabilizer))
    return records


def relabel(edges, g, tournament):
    if tournament:
        return tuple(sorted((g[i], g[j]) for i, j in edges))
    return tuple(sorted(tuple(sorted((g[i], g[j]))) for i, j in edges))


def full_type(edges, perms, tournament):
    return min(relabel(edges, g, tournament) for g in perms)


def global_even(edges, perms):
    edge_tuple = tuple(sorted(edges))
    return all(sum(g[i] > g[j] for i, j in edges) % 2 == 0
               for g in perms if relabel(edges, g, False) == edge_tuple)


def n6_control():
    n = 6
    a = (1, 2, 0, 3, 4, 5)
    b = (0, 1, 2, 4, 5, 3)
    full = list(permutations(range(n)))
    out, all_t, all_e = [], [], []
    for label, gens, expected in [
            ('single C3', [a], (128, 12, 11)),
            ('diagonal C3', [compose(a, b)], (32, 4, 3)),
            ('C3 x C3', [a, b], (8, 1, 1))]:
        K = generated(gens, n)
        normalizer = [g for g in full if conjugate(K, g) == K]
        cube = Cube(n, K)
        T = relative_records(cube, normalizer, True)
        E = relative_records(cube, normalizer, False)
        t_counts = Counter(full_type(cube.materialize(x, True), full, True) for x, _, _ in T)
        e_counts, bad_edges = Counter(), []
        for x, _, _ in E:
            graph = cube.materialize(x, False)
            if global_even(graph, full):
                e_counts[full_type(graph, full, False)] += 1
            else:
                bad_edges.append(len(graph))
        check((1 << cube.d, len(T), len(t_counts)) == expected, label + ' tournament counts')
        check(len(T) == len(E) and len(t_counts) == len(e_counts), label + ' graph counts')
        all_t.append(set(t_counts)); all_e.append(set(e_counts))
        out.append({'action': label, 'fixed_objects': 1 << cube.d,
                    'relative_types': len(T), 'global_types': len(t_counts),
                    'tournament_fusion_multiplicities': sorted(t_counts.values()),
                    'globally_odd_relative_graph_edge_counts': bad_edges})
    check(all_t[0] & all_t[1] == all_t[2], 'joint tournament intersection')
    check(all_e[0] & all_e[1] == all_e[2], 'joint graph intersection')
    check(out[0]['globally_odd_relative_graph_edge_counts'] == [9], 'K3,3 edge count')
    check(out[1]['globally_odd_relative_graph_edge_counts'] == [6], 'two triangles edge count')
    return out


def determinant(A):
    A = [[Fraction(v) for v in row] for row in A]
    result = Fraction(1)
    for i in range(len(A)):
        pivot = next(j for j in range(i, len(A)) if A[j][i])
        if pivot != i:
            A[i], A[pivot] = A[pivot], A[i]
            result = -result
        v = A[i][i]
        result *= v
        for j in range(i + 1, len(A)):
            factor = A[j][i] / v
            A[j] = [b - factor * a for a, b in zip(A[i], A[j])]
    check(result.denominator == 1, 'integral determinant')
    return result.numerator


def c7_control():
    K = generated([tuple(list(range(1, 7)) + [0, 7, 8, 9])], 10)
    gamma = [tuple([(a * i) % 7 for i in range(7)] + list(r))
             for a, r in product(range(1, 7), permutations((7, 8, 9)))]
    cube = Cube(10, K)
    actions = dict(zip(gamma, map(cube.action, gamma)))
    T = relative_records(cube, gamma, True)
    E = relative_records(cube, gamma, False)
    check(cube.d == 9 and len(T) == len(E) == 24, 'C7 cube size and orbit count')
    M = [[sum(eps * (-1 if (t & y).bit_count() % 2 else 1) for y, eps in copies.items())
          for _, copies, _ in E] for t, _, _ in T]
    for i, (_, _, tstab) in enumerate(T):
        for j in range(len(T)):
            check(sum(M[i][k] * len(E[k][2]) * M[j][k] for k in range(len(E)))
                  == ((1 << cube.d) * len(tstab) if i == j else 0), 'weighted identity')
    det = determinant(M)
    check(abs(det) == 2 ** 88, 'C7 determinant exponent')
    a = tuple([(2 * i) % 7 for i in range(7)] + [7, 8, 9])
    b = tuple(list(range(7)) + [8, 9, 7])
    Q = generated([a, b], 10)
    canon = lambda H: min(tuple(sorted(conjugate(H, g))) for g in gamma)
    tp = [stab & Q for _, _, stab in T]
    ep = [stab & Q for _, _, stab in E]
    for records, ps in [(T, tp), (E, ep)]:
        for (_, _, stab), p in zip(records, ps):
            weight = len(stab)
            while weight % 3 == 0:
                weight //= 3
            check(len(p) * weight == len(stab), 'quotient Sylow orders')
    checked_units = 0
    for i in range(24):
        for j in range(24):
            if M[i][j] % 3:
                checked_units += 1
                copies = E[j][1]
                fixed = [y for y in copies if all(act(y, actions[g], False) == y for g in tp[i])]
                check(fixed, 'localization witness for every unit coefficient')
                localized = sum(copies[y] * (-1 if (T[i][0] & y).bit_count() % 2 else 1) for y in fixed)
                check((M[i][j] - localized) % 3 == 0, 'signed localization congruence')
    # Direct Hall augmenting search in equal Sylow-action blocks.
    owners = {}
    def augment(i, visited):
        for j in range(24):
            if j in visited or M[i][j] % 3 == 0 or canon(tp[i]) != canon(ep[j]):
                continue
            visited.add(j)
            if j not in owners or augment(owners[j], visited):
                owners[j] = i
                return True
        return False
    check(all(augment(i, set()) for i in range(24)), 'full action-preserving unit matching')
    predicates = [[], [a], [b], [compose(a, b)], [a, b]]
    def signature(record, tournament):
        x, copies, _ = record
        return tuple(any(all(act(y, actions[g], tournament) == y for g in gens)
                         for y in copies) for gens in predicates)
    ts = list(map(lambda x: signature(x, True), T))
    es = list(map(lambda x: signature(x, False), E))
    counts_t = [sum(s[k] for s in ts) for k in range(5)]
    counts_e = [sum(s[k] for s in es) for k in range(5)]
    check(counts_t == counts_e == [24, 12, 4, 2, 2], 'all extension counts')
    check(all(ts[i] == es[j] for j, i in owners.items()), 'joint predicates preserved')
    hist = Counter(str(s[1:4]) for s in ts)
    check(hist == Counter(str(s[1:4]) for s in es), 'joint predicate histogram')
    H = generated([tuple(list(range(1, 7)) + [0, 7, 8, 9]), compose(a, b)], 10)
    root_swap = tuple(list(range(7)) + [8, 7, 9])
    H2 = conjugate(H, root_swap)
    L = generated(list(H | H2), 10)
    check(len(H) == 21 and len(L) == 63, 'diagonal nonpronormal witness orders')
    check(all(conjugate(H, g) != H2 for g in L), 'diagonal F21 nonpronormality')
    return {'dimension': cube.d, 'relative_types_each': 24, 'determinant_abs_exponent': 88,
            'all_q_unit_entries_localized': checked_units, 'extension_counts_each': counts_t,
            'joint_histogram': dict(sorted(hist.items())), 'diagonal_F21_nonpronormal': True}


def c15_boundary():
    n = 10
    first = tuple([1, 2, 3, 4, 0, 5, 6, 7, 8, 9])
    swap = tuple(list(range(5, 10)) + list(range(5)))
    K = generated([first], n)
    K2 = conjugate(K, swap)
    L = generated(list(K | K2), n)
    check(len(L) == 25 and all(conjugate(K, g) != K2 for g in L), 'C5 nonpronormal witness')
    complete_bipartite = {(i, j) for i in range(5) for j in range(5, 10)}
    check(relabel(complete_bipartite, swap, False) == tuple(sorted(complete_bipartite)), 'part swap automorphism')
    check(sum(swap[i] > swap[j] for i, j in complete_bipartite) == 25, 'odd K5,5 part swap')
    graph_counts = []
    for roots in range(3):
        weights = [3, 5] + [1] * roots
        quotient_pairs = list(combinations(range(len(weights)), 2))
        block_of = [i for i, w in enumerate(weights) for _ in range(w)]
        nv = len(block_of)
        types = set()
        for mask in range(1 << len(quotient_pairs)):
            quotient_edges = {e for k, e in enumerate(quotient_pairs) if mask >> k & 1}
            adj = [[int(tuple(sorted((block_of[i], block_of[j]))) in quotient_edges)
                    for j in range(nv)] for i in range(nv)]
            unseen, classes = set(range(nv)), []
            while unseen:
                i = min(unseen)
                twin_class = tuple(j for j in sorted(unseen)
                                   if all(adj[i][k] == adj[j][k] for k in range(nv) if k not in (i, j)))
                classes.append(twin_class)
                unseen.difference_update(twin_class)
            size = tuple(map(len, classes))
            internal = tuple(adj[c[0]][c[1]] if len(c) > 1 else 0 for c in classes)
            if any(internal):
                continue  # A transposition within a clique-twin class is odd.
            m = len(classes)
            qa = [[adj[classes[i][0]][classes[j][0]] if i != j else 0 for j in range(m)] for i in range(m)]
            orders = list(permutations(range(m)))
            even = True
            for order in orders:
                if any(size[i] != size[order[i]] for i in range(m)):
                    continue
                if any(qa[i][j] != qa[order[i]][order[j]] for i, j in combinations(range(m), 2)):
                    continue
                reversals = sum(size[i] * size[j] * qa[i][j]
                                for i, j in combinations(range(m), 2) if order[i] > order[j])
                if reversals % 2:
                    even = False
            if even:
                types.add(min((tuple(size[i] for i in order),
                               tuple(qa[order[i]][order[j]] for i, j in combinations(range(m), 2)))
                              for order in orders))
        graph_counts.append(len(types))
    check(graph_counts == [2, 8, 31], 'C15 intrinsic twin-quotient graph counts')
    # Tournament counts follow from the unique prime C3 and C5 modules; only
    # the r <= 2 remaining singleton vertices can be permuted.
    tournament_counts = [2, 8, 32]
    return {'C5_generated_group_order': 25, 'C5_is_nonpronormal': True,
            'part_swap_edge_reversals': 25, 'even_graph_counts_n8_n9_n10': graph_counts,
            'tournament_counts_from_prime_module_proof': tournament_counts,
            'graph_quotient_presentations_checked': 2 + 8 + 64}


if __name__ == '__main__':
    print(json.dumps({'status': 'PASS', 'n6': n6_control(), 'C7_quotient': c7_control(),
                      'C15_boundary': c15_boundary()}, indent=2, sort_keys=True))

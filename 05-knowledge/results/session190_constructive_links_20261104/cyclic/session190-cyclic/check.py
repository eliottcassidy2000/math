#!/usr/bin/env python3
"""Independent standard-library C15 orbit-defect audit at degree ten.

No repository imports, nauty, external graph library, or canonical labels.
Generate every H-fixed tournament and graph directly from H-orbits on pairs;
quotient by the complete normalizer; compare representatives under unrestricted
vertex isomorphism by exhaustive backtracking; test graph evenness using an
exact stabilizer-chain transversal, not a heuristic sample of automorphisms.
"""
import itertools
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


N = 10
PAIRS = [(i, j) for i in range(N) for j in range(i + 1, N)]
H = (1, 2, 0, 4, 5, 6, 7, 3, 8, 9)


def pair_orbits():
    todo = set(PAIRS)
    out = []
    while todo:
        e = min(todo)
        orbit = []
        current = e
        while current not in orbit:
            orbit.append(current)
            current = (H[current[0]], H[current[1]])
        require(current == e, 'Oriented orbit did not close')
        undirected = {tuple(sorted(x)) for x in orbit}
        require(len(undirected) == len(orbit), 'H reverses a pair')
        require(undirected <= todo, 'Pair orbit overlap')
        todo -= undirected
        out.append(orbit)
    return out


def normalizer():
    group = []
    for a, b, c, d, swap in itertools.product(range(1, 3), range(3), range(1, 5), range(5), range(2)):
        g = tuple([(a * i + b) % 3 for i in range(3)]
                  + [3 + (c * j + d) % 5 for j in range(5)]
                  + ([9, 8] if swap else [8, 9]))
        group.append(g)
    require(len(set(group)) == 240, 'Normalizer size')
    return group


def action(edges, g, directed):
    if directed:
        return frozenset((g[u], g[v]) for u, v in edges)
    return frozenset(tuple(sorted((g[u], g[v]))) for u, v in edges)


def sign(g, edges):
    return -1 if sum(g[u] > g[v] for u, v in edges) % 2 else 1


def all_fixed(orbits, directed):
    objects = []
    for bits in itertools.product(range(2), repeat=len(orbits)):
        edges = []
        for orbit, bit in zip(orbits, bits):
            if directed:
                edges.extend(orbit if bit else [(v, u) for u, v in orbit])
            elif bit:
                edges.extend(tuple(sorted(e)) for e in orbit)
        objects.append(frozenset(edges))
    require(len(set(objects)) == 512, 'H-fixed cube size')
    return objects


def n_representatives(objects, group, directed):
    seen, reps = set(), []
    for obj in objects:
        if obj in seen:
            continue
        orbit = {action(obj, g, directed) for g in group}
        require(orbit <= set(objects), 'Normalizer leaves H-fixed cube')
        seen |= orbit
        if directed or all(sign(g, obj) == 1 for g in group if action(obj, g, False) == obj):
            reps.append(obj)
    require(len(seen) == len(objects), 'Normalizer coverage')
    return reps


def matrix(edges, directed):
    adj = [[0] * N for _ in range(N)]
    for u, v in edges:
        adj[u][v] = 1
        if not directed:
            adj[v][u] = 1
    return adj


def invariants(adj):
    # Refinement is only a safe pruning invariant, never the isomorphism test.
    colours = [sum(row) for row in adj]
    for _ in range(N):
        profiles = [(colours[v], tuple(sorted((adj[v][w], adj[w][v], colours[w])
                                              for w in range(N) if w != v))) for v in range(N)]
        unique = {x: i for i, x in enumerate(sorted(set(profiles)))}
        new = [unique[x] for x in profiles]
        if new == colours:
            break
        colours = new
    return colours


def find_isomorphism(a, b, prescribed=None):
    """Return a checked isomorphism or exhaust every compatible permutation."""
    ca, cb = invariants(a), invariants(b)
    if sorted(ca) != sorted(cb) or sorted(map(sum, a)) != sorted(map(sum, b)):
        return None
    assigned = {} if prescribed is None else dict(prescribed)
    if len(set(assigned.values())) != len(assigned):
        return None
    for x, y in assigned.items():
        if ca[x] != cb[y] or sum(a[x]) != sum(b[y]):
            return None
        for u, v in assigned.items():
            if a[x][u] != b[y][v]:
                return None

    def recurse():
        if len(assigned) == N:
            g = tuple(assigned[i] for i in range(N))
            require(all(a[i][j] == b[g[i]][g[j]] for i in range(N) for j in range(N)), 'Invalid isomorphism')
            return g
        used = set(assigned.values())
        best = None
        for x in range(N):
            if x in assigned:
                continue
            candidates = [y for y in range(N) if y not in used and ca[x] == cb[y]
                          and sum(a[x]) == sum(b[y])
                          and all(a[x][u] == b[y][v] and a[u][x] == b[v][y]
                                  for u, v in assigned.items())]
            if not candidates:
                return None
            if best is None or len(candidates) < len(best[1]):
                best = (x, candidates)
        x, candidates = best
        for y in candidates:
            assigned[x] = y
            found = recurse()
            del assigned[x]
            if found is not None:
                return found
        return None

    return recurse()


def automorphism_even(edges):
    adj = matrix(edges, False)
    # At level i, collect one coset representative taking i to each point of
    # its orbit under the pointwise stabilizer of 0,...,i-1. Together over
    # all levels these generate Aut(adj), by induction down the stabilizer
    # chain. Thus checking the sign character on these elements is exact.
    tested = 0
    for i in range(N):
        for j in range(i + 1, N):
            fixed = {k: k for k in range(i)}
            fixed[i] = j
            g = find_isomorphism(adj, adj, fixed)
            if g is not None:
                tested += 1
                if sign(g, edges) == -1:
                    return False, g, tested
    return True, None, tested


def type_classes(reps, directed):
    types, matrices = [], []
    for obj in reps:
        a = matrix(obj, directed)
        for i, b in enumerate(matrices):
            if find_isomorphism(a, b) is not None:
                types[i].append(obj)
                break
        else:
            matrices.append(a)
            types.append([obj])
    return types


def twin_quotient(edges):
    """Independent graph classification and sign test on <=4 false-twin cells.

    Equal open-neighbourhood classes are intrinsic independent sets. Every
    automorphism induces a size-preserving quotient automorphism, and every
    such quotient automorphism lifts. Permutations within a cell reverse no
    edges when the cells are ordered consecutively. A cell inversion across
    an edge contributes the product of the two cell sizes to reversal parity.
    """
    adj = matrix(edges, False)
    classes = {}
    for v, row in enumerate(adj):
        classes.setdefault(tuple(row), []).append(v)
    cells = list(classes.values())
    m = len(cells)
    require(m <= 4, 'Unexpected false-twin quotient size')
    q = [[adj[a[0]][b[0]] for b in cells] for a in cells]
    sizes = [len(c) for c in cells]
    canonical = None
    even = True
    for perm in itertools.permutations(range(m)):
        code = (tuple(sizes[i] for i in perm),
                tuple(q[perm[i]][perm[j]] for i in range(m) for j in range(i + 1, m)))
        canonical = code if canonical is None else min(canonical, code)
        if all(sizes[i] == sizes[perm[i]] for i in range(m)) and all(
                q[i][j] == q[perm[i]][perm[j]] for i in range(m) for j in range(m)):
            parity = sum(sizes[i] * sizes[j] for i in range(m) for j in range(i + 1, m)
                         if q[i][j] and perm[i] > perm[j]) % 2
            if parity:
                even = False
    return canonical, even


def main():
    orbits = pair_orbits()
    group = normalizer()
    t_reps = n_representatives(all_fixed(orbits, True), group, True)
    g_reps = n_representatives(all_fixed(orbits, False), group, False)
    t_types = type_classes(t_reps, True)
    g_types = type_classes(g_reps, False)
    even, odd, transversals = [], [], 0
    for objects in g_types:
        good, witness, tested = automorphism_even(objects[0])
        transversals += tested
        if good:
            even.append(objects)
        else:
            odd.append((objects, witness))
    require(len(t_reps) == len(g_reps) == 32, 'Normalizer orbit counts')
    require(len(t_types) == 32 and len(even) == 31 and len(odd) == 1, 'Defect counts')
    witness_edges = odd[0][0][0]
    part_a, part_b = {0, 1, 2, 8, 9}, {3, 4, 5, 6, 7}
    k55 = frozenset(tuple(sorted((a, b))) for a in part_a for b in part_b)
    require(witness_edges == k55, 'Unique odd representative must be K5,5')
    side_swap = (3, 4, 5, 0, 1, 2, 8, 9, 6, 7)
    require(action(k55, side_swap, False) == k55 and sign(side_swap, k55) == -1, 'Explicit odd side swap')
    require(all(sign(g, k55) == 1 for g in group), 'K5,5 is normalizer-even')
    quotient_types = {}
    for edges in g_reps:
        code, good = twin_quotient(edges)
        old = quotient_types.setdefault(code, good)
        require(old == good, 'Twin quotient sign inconsistency')
        require(good == automorphism_even(edges)[0], 'Independent twin quotient sign check')
    require(len(quotient_types) == 32 and sum(quotient_types.values()) == 31,
            'Independent twin quotient count check')
    result = dict(n=N, H_generator=H, normalizer_order=len(group), pair_orbits=len(orbits),
                  H_fixed_tournaments=512, H_fixed_graphs=512,
                  normalizer_tournament_orbits=len(t_reps), normalizer_even_graph_orbits=len(g_reps),
                  tau_T=len(t_types), tau_E=len(even), defect=len(t_types)-len(even),
                  W=len(odd), E_T=len(t_reps)-len(t_types),
                  E_even=sum(len(x)-1 for x in even), E_wit=sum(len(x)-1 for x, _ in odd),
                  witness='K_{5,5}', witness_parts=[sorted(part_a), sorted(part_b)],
                  explicit_odd_side_swap=side_swap, tested_automorphism_transversals=transversals,
                  independent_twin_quotient_types=len(quotient_types),
                  independent_twin_quotient_even_types=sum(quotient_types.values()),
                  all_checks_passed=True)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

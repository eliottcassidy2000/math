#!/usr/bin/env python3
"""Exact controls for octonion permutation transports. Standard library only."""
import hashlib
import json
from collections import deque
from itertools import combinations, product


def require(value, message):
    if not value:
        raise RuntimeError(message)


def qmul(a, b):
    if a == 0:
        return 1, b
    if b == 0:
        return 1, a
    if a == b:
        return -1, 0
    return (1 if (a, b) in ((1, 2), (2, 3), (3, 1)) else -1), a ^ b


def basis_mul(a, b):
    """Explicit quaternion-pair Cayley-Dickson rule from the note."""
    A, B = a % 4, b % 4
    ah, bh = a // 4, b // 4
    if ah == bh == 0:
        return qmul(A, B)
    if ah == 0 and bh == 1:
        s, c = qmul(B, A)
        return s, c + 4
    if ah == 1 and bh == 0:
        s, c = qmul(A, B)
        return s * (-1 if B else 1), c + 4
    s, c = qmul(B, A)
    return -s * (-1 if B else 1), c


def recursive_basis_mul(level, a, b):
    """Independent recursion, beginning at the real algebra."""
    if level == 0:
        require(a == b == 0, 'real basis')
        return 1, 0
    half = 1 << (level - 1)
    A, B = a % half, b % half
    ah, bh = a // half, b // half
    if ah == 0 and bh == 0:
        return recursive_basis_mul(level - 1, A, B)
    if ah == 0:
        sign, index = recursive_basis_mul(level - 1, B, A)
        return sign, index + half
    if bh == 0:
        sign, index = recursive_basis_mul(level - 1, A, B)
        return sign * (-1 if B else 1), index + half
    sign, index = recursive_basis_mul(level - 1, B, A)
    return -sign * (-1 if B else 1), index


def mul(u, v):
    s, a = u
    t, b = v
    sign, c = basis_mul(a, b)
    return s * t * sign, c


U = tuple(product((1, -1), range(8)))
INDEX = {u: i for i, u in enumerate(U)}
I = tuple(range(16))
Z = tuple(INDEX[(-s, a)] for s, a in U)
L = [tuple(INDEX[mul((1, a), u)] for u in U) for a in range(8)]


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    out = [None] * len(p)
    for i, v in enumerate(p):
        out[v] = i
    return tuple(out)


def generated(gens, degree=16):
    identity = tuple(range(degree))
    known, todo = {identity}, [identity]
    while todo:
        p = todo.pop()
        for q in gens:
            r = compose(q, p)
            if r not in known:
                known.add(r)
                todo.append(r)
    return known


def dot(a, b):
    return (a & b).bit_count() % 2


def determinant(a, b, c):
    av, bv, cv = [[(x >> i) & 1 for i in range(3)] for x in (a, b, c)]
    return (av[0]*bv[1]*cv[2] + av[0]*bv[2]*cv[1]
            + av[1]*bv[0]*cv[2] + av[1]*bv[2]*cv[0]
            + av[2]*bv[0]*cv[1] + av[2]*bv[1]*cv[0]) % 2


def diagonal(c, b):
    return tuple(INDEX[(s * (-1 if c ^ dot(b, a) else 1), a)] for s, a in U)


def centralizer_from_root(gens, degree):
    """Enumerate all commuting permutations for a transitive action.

    Independently propagates each root image. It does not assume the group
    centralizer formula or that a commuting permutation is in the group.
    """
    transitions = list(gens) + [inverse(p) for p in gens]
    answers = []
    for root_image in range(degree):
        values, queue, valid = {0: root_image}, deque([0]), True
        while queue and valid:
            x = queue.popleft()
            for p in transitions:
                next_x, next_y = p[x], p[values[x]]
                if next_x in values:
                    if values[next_x] != next_y:
                        valid = False
                        break
                else:
                    values[next_x] = next_y
                    queue.append(next_x)
        if valid:
            require(len(values) == degree, 'centralizer action must be transitive')
            candidate = tuple(values[i] for i in range(degree))
            require(len(set(candidate)) == degree, 'commuting map permutation')
            require(all(compose(candidate, p) == compose(p, candidate) for p in gens),
                    'centralizer all-generator check')
            answers.append(candidate)
    return answers


def adjacency(vertex_count, edges):
    neighbors = [set() for _ in range(vertex_count)]
    require(len(edges) == len(set(edges)), 'simple lifted edge set')
    for v, w in edges:
        require(v != w, 'simple graph has no loops')
        neighbors[v].add(w)
        neighbors[w].add(v)
    return neighbors


def component_sizes(neighbors):
    unseen, sizes = set(range(len(neighbors))), []
    while unseen:
        root = min(unseen)
        seen, todo = {root}, [root]
        while todo:
            for w in neighbors[todo.pop()]:
                if w not in seen:
                    seen.add(w)
                    todo.append(w)
        unseen.difference_update(seen)
        sizes.append(len(seen))
    return sorted(sizes)


def make_cover(base_edges, voltages, degree):
    edges = []
    for v, w in base_edges:
        p = voltages[(v, w)]
        require(sorted(p) == list(range(degree)), 'edge voltage permutation')
        for x in range(degree):
            edges.append((v * degree + x, w * degree + p[x]))
    neighbors = adjacency(5 * degree, edges)
    for v in range(5):
        expected = set(range(5)) - {v}
        for x in range(degree):
            actual = [w // degree for w in neighbors[v * degree + x]]
            require(len(actual) == 4 and set(actual) == expected, 'local cover bijection')
    return edges, neighbors


def main():
    for a, b in product(range(8), repeat=2):
        require(basis_mul(a, b) == recursive_basis_mul(3, a, b), 'independent product recursion')
        require(basis_mul(a, b)[1] == a ^ b, 'index XOR')
    nonassociative_triples = 0
    for a, b, c in product(range(8), repeat=3):
        left = mul((1, a), mul((1, b), (1, c)))
        right = mul(mul((1, a), (1, b)), (1, c))
        expected = right[0] * (-1 if determinant(a, b, c) else 1), right[1]
        require(left == expected, 'associator determinant')
        nonassociative_triples += left != right
    require(nonassociative_triples == 168, 'associator triple count')
    D = {(c, b): diagonal(c, b) for c in (0, 1) for b in range(8)}
    require(len(set(D.values())) == 16, 'all characters distinct')
    for a, b in product(range(8), repeat=2):
        s, c = basis_mul(a, b)
        lab = L[c] if s == 1 else compose(Z, L[c])
        defect = compose(inverse(lab), compose(L[a], L[b]))
        expected = tuple(INDEX[(sign * (-1 if determinant(a, b, x) else 1), x)]
                         for sign, x in U)
        require(defect == expected, 'full permutation defect')
        for character in range(8):
            for parity in (0, 1):
                conjugated = compose(L[a], compose(D[(parity, character)], inverse(L[a])))
                require(conjugated == D[(parity ^ dot(character, a), character)], 'character conjugation')
    require(compose(inverse(L[3]), compose(L[1], L[2])) == D[(0, 4)], 'named defect')
    for a in range(1, 8):
        require(compose(L[a], L[a]) == Z, 'unit squares')
        for b in range(1, 8):
            if a != b:
                require(compose(L[a], L[b]) == compose(Z, compose(L[b], L[a])), 'Clifford relation')
    volume = I
    for a in range(1, 8):
        volume = compose(volume, L[a])
    require(volume == Z, 'volume relation')
    G = generated(L[1:])
    normal_forms = {compose(d, l) for d in D.values() for l in L}
    require(len(G) == 128 and G == normal_forms, '128 distinct complete normal forms')
    require(generated(L[1:7]) == G, 'six generators suffice')
    for p in G:
        shift = U[p[0]][1]
        require(all(U[p[i]][1] == U[i][1] ^ shift for i in range(16)), 'translation quotient')
    kernel = {p for p in G if U[p[0]][1] == 0}
    require(kernel == set(D.values()), 'kernel character group')
    stabilizer = {p for p in G if p[0] == 0}
    require(stabilizer == {D[(0, b)] for b in range(8)}, 'root stabilizer')
    require(len({p[0] for p in G}) == 16, 'transitive action')
    center = {p for p in G if all(compose(p, l) == compose(l, p) for l in L[1:])}
    require(center == {I, Z}, 'center')
    commutators = {compose(p, compose(q, compose(inverse(p), inverse(q)))) for p in G for q in G}
    require(commutators == {I, Z}, 'all commutators')
    require(all(compose(p, p) in {I, Z} for p in G), 'elementary abelian quotient')
    centralizer = centralizer_from_root(L[1:7], 16)
    require(set(centralizer) == {I, Z}, 'full permutation centralizer')
    subset_counts = {}
    for size in range(1, 8):
        group_sizes, transitive_count = set(), 0
        for subset in combinations(range(1, 8), size):
            H = generated([L[a] for a in subset])
            group_sizes.add(len(H))
            transitive_count += len({h[0] for h in H}) == 16
        require(group_sizes == {2 ** (min(size, 6) + 1)}, 'generator subset sizes')
        subset_counts[str(size)] = {'subsets': len(list(combinations(range(1,8),size))),
                                   'group_order': next(iter(group_sizes)),
                                   'transitive_subsets': transitive_count}
    base_edges = list(combinations(range(5), 2))
    chords = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    voltage16 = {edge: I for edge in base_edges}
    voltage8 = {edge: tuple(range(8)) for edge in base_edges}
    for label, edge in enumerate(chords, 1):
        voltage16[edge] = L[label]
        voltage8[edge] = tuple(x ^ label for x in range(8))
    edges16, neighbors16 = make_cover(base_edges, voltage16, 16)
    edges8, neighbors8 = make_cover(base_edges, voltage8, 8)
    require(component_sizes(neighbors16) == [80], '16-sheet cover connected')
    require(component_sizes(neighbors8) == [40], '8-sheet cover connected')
    require(len(edges16) == 160 and len(edges8) == 80, 'edge counts')
    project = [8 * (v // 16) + U[v % 16][1] for v in range(80)]
    for v in range(80):
        images = [project[w] for w in neighbors16[v]]
        require(len(images) == len(set(images)) and set(images) == neighbors8[project[v]], 'double cover local bijection')
    require(all(project.count(v) == 2 for v in range(40)), 'double cover fibers')
    centralizer8 = centralizer_from_root([voltage8[e] for e in chords], 8)
    require(len(centralizer8) == 8, '8-sheet cover regular deck count')
    defect = D[(0, 4)]
    fixed = [i for i in range(16) if defect[i] == i]
    transpositions = [[i, defect[i]] for i in range(16) if i < defect[i]]
    require(len(fixed) == 8 and len(transpositions) == 4, 'defect cycle structure')
    # Order-eight phase action on the four odd eighths of a turn.
    rotation = (1, 2, 3, 0)
    reflection = (3, 2, 1, 0)
    require(len(generated([rotation, reflection], 4)) == 8, 'odd-eighth dihedral group')
    out = {
        'status': 'PASS',
        'product_convention': '(a,b)(c,d)=(ac-conj(d)b, da+b conj(c))',
        'basis_products_compared_independently': 64,
        'associator_triples_checked': 512,
        'nonassociative_basis_triples': nonassociative_triples,
        'transport_group_order': len(G),
        'transport_group_sha256': hashlib.sha256(json.dumps(sorted(G)).encode()).hexdigest(),
        'translation_quotient_order': 8,
        'character_kernel_order': len(kernel),
        'root_stabilizer_order': len(stabilizer),
        'center_order': len(center),
        'centralizer_in_S16_order': len(centralizer),
        'abelianization_order': 64,
        'minimum_monodromy_generators': 6,
        'generator_subset_controls': subset_counts,
        'defect_L3_inverse_L1_L2': {'fixed_points': len(fixed), 'transpositions': len(transpositions),
                                  'fixed_indices': [0,1,2,3], 'sign_swapped_indices': [4,5,6,7]},
        'K5_cover': {'base_vertices': 5, 'base_edges': 10, 'base_cycle_rank': 6,
                     'E_vertices': len(neighbors16), 'E_edges': len(edges16), 'E_components': 1,
                     'C_vertices': len(neighbors8), 'C_edges': len(edges8), 'C_components': 1,
                     'E_over_C_degree': 2, 'E_over_C_regular': True,
                     'C_over_base_degree': 8, 'C_over_base_deck_order': len(centralizer8),
                     'C_over_base_regular': True, 'E_over_base_degree': 16,
                     'E_over_base_deck_order': len(centralizer), 'E_over_base_regular': False},
        'cube_cycle_rank': 5,
        'odd_eighth_phase_dihedral_order': 8,
        'scope': 'Formal permutation cover. No tournament realization or knot/type-theory equivalence claimed.'
    }
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

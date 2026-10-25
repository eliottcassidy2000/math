"""Literal finite checks of the clone-cut sign-reversing involution."""
from itertools import combinations, permutations
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def graph(n, mask):
    arcs = set()
    for i, (a, b) in enumerate(combinations(range(n), 2)):
        arcs.add((a, b) if (mask >> i) & 1 else (b, a))
    return frozenset(arcs)


def code(n, arcs, colors=None):
    colors = (0,) * n if colors is None else tuple(colors)
    pairs = list(combinations(range(n), 2))
    return min((tuple(colors[p[i]] for i in range(n)),
                sum(((p[a], p[b]) in arcs) << j for j, (a, b) in enumerate(pairs)))
               for p in permutations(range(n)))


def automorphisms(n, arcs):
    return [p for p in permutations(range(n)) if all((p[a], p[b]) in arcs for a, b in arcs)]


def modules(n, arcs):
    neighbors = {v: set() for v in range(n)}
    for a, b in combinations(range(n), 2):
        if all(((a, c) in arcs) == ((b, c) in arcs) for c in range(n) if c not in (a, b)):
            neighbors[a].add(b)
            neighbors[b].add(a)
    result = []
    unseen = set(range(n))
    while unseen:
        component = {min(unseen)}
        frontier = list(component)
        while frontier:
            for v in neighbors[frontier.pop()]:
                if v not in component:
                    component.add(v)
                    frontier.append(v)
        unseen -= component
        chain = tuple(sorted(component, key=lambda v: sum((u, v) in arcs for u in component)))
        need(all((chain[i], chain[j]) in arcs for i in range(len(chain)) for j in range(i + 1, len(chain))),
             'component is transitive in its intrinsic order')
        need(all(len({(v, u) in arcs for u in component}) == 1 for v in range(n) if v not in component),
             'component is a module')
        need(sum(len(neighbors[v]) for v in component) == 2 * (len(chain) - 1), 'homogeneous-pair component is a path')
        result.append(chain)
    return result


def quotient_key(arcs, chains, gaps, cut):
    index = {g: i for i, g in enumerate(gaps)}
    blocks = []
    for chain in chains:
        block = [chain[0]]
        for a, b in zip(chain, chain[1:]):
            if cut >> index[(a, b)] & 1:
                blocks.append(tuple(block))
                block = []
            block.append(b)
        blocks.append(tuple(block))
    q_arcs = set()
    for i, j in combinations(range(len(blocks)), 2):
        direction = (blocks[i][0], blocks[j][0]) in arcs
        need(all(((a, b) in arcs) == direction for a in blocks[i] for b in blocks[j]), 'uniform block pair')
        q_arcs.add((i, j) if direction else (j, i))
    need(sum(len(b) for b in blocks) == sum(len(c) for c in chains), 'same literal inflated size')
    return code(len(blocks), q_arcs, [len(block) - 1 for block in blocks])


def inspect(n, arcs):
    chains = modules(n, arcs)
    gaps = [g for chain in chains for g in zip(chain, chain[1:])]
    index = {g: i for i, g in enumerate(gaps)}
    group = automorphisms(n, arcs)
    need(len(group) % 2 == 1, 'odd automorphism order')
    actions = [tuple(index[(p[a], p[b])] for a, b in gaps) for p in group]

    def act(cut, perm):
        return sum(1 << perm[i] for i in range(len(gaps)) if cut >> i & 1)

    def orbit_code(cut):
        return min(act(cut, perm) for perm in actions)

    nontrivial = [chain for chain in chains if len(chain) > 1]
    if nontrivial:
        rooted = {chain: code(n, arcs, [int(v in chain) for v in range(n)]) for chain in nontrivial}
        minimum = min(rooted.values())
        chosen = [chain for chain in nontrivial if rooted[chain] == minimum]
        seed = chosen[0]
        literal_orbit = {frozenset(p[v] for v in seed) for p in group}
        need(literal_orbit == {frozenset(chain) for chain in chosen}, 'canonical rooted code chooses one orbit')
        toggle = sum(1 << index[(chain[0], chain[1])] for chain in chosen)
        need(toggle.bit_count() % 2 == 1, 'odd gap set')
        need(all(act(toggle, perm) == toggle for perm in actions), 'invariant gap set')
    else:
        toggle = 0

    q = len(chains)
    key_to_orbit = {}
    orbit_to_key = {}
    signed = 0
    for cut in range(1 << len(gaps)):
        canonical = orbit_code(cut)
        key = quotient_key(arcs, chains, gaps, cut)
        need(key_to_orbit.setdefault(key, canonical) == canonical, 'weighted quotient equivalence implies cut orbit')
        need(orbit_to_key.setdefault(canonical, key) == key, 'cut orbit implies weighted quotient equivalence')
        need(len(key[0]) == q + cut.bit_count(), 'block count')
        if cut == canonical:
            signed += (-1) ** (n - q - cut.bit_count())
        if toggle:
            paired = cut ^ toggle
            need((paired ^ toggle) == cut, 'involution')
            need((paired.bit_count() - cut.bit_count()) % 2 == 1, 'sign reversal')
            need(orbit_code(paired) != canonical, 'no fixed orbit')
            need(all(act(paired, perm) == (act(cut, perm) ^ toggle) for perm in actions), 'equivariance')
    need(signed == int(not gaps), 'fiber signed sum')
    return {'modules': [len(c) for c in chains], 'gap_count': len(gaps),
            'automorphisms': len(group), 'cut_orbits': len(orbit_to_key),
            'signed_sum': signed, 'toggle_size': toggle.bit_count()}, orbit_code, gaps


def cyclic_inflation(sizes):
    blocks = []
    start = 0
    for size in sizes:
        blocks.append(tuple(range(start, start + size)))
        start += size
    arcs = set()
    for block in blocks:
        arcs.update(combinations(block, 2))
    for i, j in ((0, 1), (1, 2), (2, 0)):
        arcs.update((a, b) for a in blocks[i] for b in blocks[j])
    return start, frozenset(arcs)


def main():
    counts = []
    cut_orbits = 0
    for n in range(6):
        seen = set()
        for mask in range(1 << (n * (n - 1) // 2)):
            arcs = graph(n, mask)
            canonical = code(n, arcs)
            if canonical in seen:
                continue
            seen.add(canonical)
            result, _, _ = inspect(n, arcs)
            cut_orbits += result['cut_orbits']
        counts.append(len(seen))
    need(counts == [1, 1, 1, 2, 4, 12], 'complete small isomorphism universe')
    n, arcs = cyclic_inflation((2, 2, 2))
    symmetric, orbit, gaps = inspect(n, arcs)
    need(symmetric == {'modules': [2, 2, 2], 'gap_count': 3, 'automorphisms': 3,
                       'cut_orbits': 4, 'signed_sum': 0, 'toggle_size': 3}, 'symmetric hostile data')
    need(orbit(1) == orbit(2) and orbit(1 ^ 1) != orbit(2 ^ 1), 'single labelled gap fails to descend')
    n, arcs = cyclic_inflation((2, 2, 1))
    asymmetric, _, _ = inspect(n, arcs)
    need(asymmetric['gap_count'] == 2 and asymmetric['toggle_size'] == 1,
         'all nontrivial modules would give an even toggle')
    print(json.dumps({'status': 'PASS', 'complete_tournament_class_counts_n0_to_5': counts,
                      'small_fiber_cut_orbits': cut_orbits, 'symmetric_C3_TT2_hostile': symmetric,
                      'two_nontrivial_orbits_control': asymmetric,
                      'scope': 'literal inflation fibers and cut-orbit involution; no enumeration speed claim'}, indent=2))


if __name__ == '__main__':
    main()

"""Explicit type bijection via SCC words, Lyndon factors, and connected atoms.

The theorem's finite algorithm applies at every order; this deliberately bounded
implementation authenticates orders <=6. Atom rank choices are retained, not
misrepresented as a relabelling-equivariant or rank-free natural construction.
"""
from functools import lru_cache
from itertools import combinations, permutations, product
from collections import Counter
import json

LIMIT = 6
CHECKS = 0


def need(ok, why):
    if not ok:
        raise ValueError(why)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def source(n, mask):
    need(type(n) is int and 0 <= n <= LIMIT, 'exact order in declared range 0..6')
    need(type(mask) is int and 0 <= mask < 1 << (n*(n-1)//2), 'exact pair-bit mask')


@lru_cache(None)
def pairs(n):
    return tuple(combinations(range(n), 2))


@lru_cache(None)
def actions(n):
    edges = pairs(n); index = {e: i for i, e in enumerate(edges)}
    out = []
    for p in permutations(range(n)):
        dst, flip, input_flip = [], 0, 0
        for i, (u, v) in enumerate(edges):
            a, b = p[u], p[v]
            j = index[tuple(sorted((a, b)))]; dst.append(j)
            if a > b:
                flip |= 1 << j; input_flip |= 1 << i
        out.append((p, tuple(dst), flip, input_flip))
    return tuple(out)


def act(mask, dst, flip=0):
    result = flip
    for i, j in enumerate(dst):
        if mask >> i & 1:
            result ^= 1 << j
    return result


@lru_cache(None)
def _catalogue(n, tournament):
    canon = [-1]*(1 << len(pairs(n))); representatives = []
    for mask in range(len(canon)):
        if canon[mask] >= 0:
            continue
        orbit = {act(mask, dst, flip if tournament else 0)
                 for _, dst, flip, _ in actions(n)}
        need(mask == min(orbit), 'increasing representative enumeration')
        representatives.append(mask)
        for item in orbit:
            canon[item] = mask
    return tuple(representatives), tuple(canon)


def canonical(n, mask, *, tournament):
    source(n, mask); need(type(tournament) is bool, 'typed object species')
    return _catalogue(n, tournament)[1][mask]


def automorphism_even(n, mask):
    source(n, mask)
    return _automorphism_even(n, canonical(n, mask, tournament=False))


@lru_cache(None)
def _automorphism_even(n, mask):
    return all((mask & input_flip).bit_count() % 2 == 0
               for _, dst, _, input_flip in actions(n) if act(mask, dst) == mask)


def adjacency(n, mask, tournament):
    result = [[False]*n for _ in range(n)]
    for i, (u, v) in enumerate(pairs(n)):
        if tournament:
            a, b = (v, u) if mask >> i & 1 else (u, v)
            result[a][b] = True
        elif mask >> i & 1:
            result[u][v] = result[v][u] = True
    return result


def induced(n, mask, vertices, tournament):
    matrix = adjacency(n, mask, tournament)
    result = 0
    for i, (a, b) in enumerate(pairs(len(vertices))):
        u, v = vertices[a], vertices[b]
        if matrix[v][u] if tournament else matrix[u][v]:
            result |= 1 << i
    return canonical(len(vertices), result, tournament=tournament)


def components(n, mask, *, tournament):
    source(n, mask); need(type(tournament) is bool, 'typed component species')
    reach = adjacency(n, mask, tournament)
    for i in range(n):
        reach[i][i] = True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] |= reach[k][j]
    unused, parts = set(range(n)), []
    while unused:
        i = min(unused)
        part = tuple(j for j in sorted(unused) if reach[i][j] and reach[j][i])
        parts.append(part); unused.difference_update(part)
    if tournament:
        parts.sort(key=lambda p: -sum(reach[p[0]][j] for j in range(n) if j not in p))
        need(all(reach[a[0]][b[0]] and not reach[b[0]][a[0]]
                 for i, a in enumerate(parts) for b in parts[i+1:]), 'total SCC order')
    return tuple((len(p), induced(n, mask, p, tournament)) for p in parts)


def assemble(parts, *, tournament):
    need(type(parts) is tuple and type(tournament) is bool, 'typed ordered parts')
    for part in parts:
        need(type(part) is tuple and len(part) == 2, 'typed component')
        source(*part); need(part[0] > 0, 'nonempty components')
    n = sum(p[0] for p in parts); source(n, 0)
    matrix = [[False]*n for _ in range(n)]; offset = 0
    for size, mask in parts:
        local = adjacency(size, mask, tournament)
        for u in range(size):
            for v in range(size):
                matrix[offset+u][offset+v] = local[u][v]
        if tournament:
            for u in range(offset):
                for v in range(offset, offset+size):
                    matrix[u][v] = True
        offset += size
    out = sum(1 << i for i, (u, v) in enumerate(pairs(n))
              if (matrix[v][u] if tournament else matrix[u][v]))
    return n, canonical(n, out, tournament=tournament)


def word_input(word):
    need(type(word) is tuple, 'ordered SCC word')
    total = 0
    for letter in word:
        need(type(letter) is tuple and len(letter) == 2, 'exact strong letter')
        n, mask = letter; source(n, mask); need(n > 0, 'positive letter weight')
        need(canonical(n, mask, tournament=True) == mask and
             len(components(n, mask, tournament=True)) == 1, 'canonical strong letter')
        total += n
    need(total <= LIMIT, 'word weight in declared range')


def lyndon(word):
    return bool(word) and all(word < word[i:] for i in range(1, len(word)))


def factor(word):
    """Smallest-suffix recursion: finite, exact, and deliberately not optimized."""
    word_input(word)
    if not word:
        return ()
    i = min(range(len(word)), key=lambda j: word[j:])
    answer = factor(word[:i]) + (word[i:],)
    need(all(lyndon(w) for w in answer) and
         all(a >= b for a, b in zip(answer, answer[1:])), 'CFL invariants')
    return answer


def weight(word):
    return sum(n for n, _ in word)


@lru_cache(None)
def _atlas(bound):
    strong, connected = {}, {}
    for n in range(1, bound+1):
        strong[n] = tuple(m for m in _catalogue(n, True)[0]
                          if len(components(n, m, tournament=True)) == 1)
        connected[n] = tuple(m for m in _catalogue(n, False)[0]
                             if automorphism_even(n, m) and
                             len(components(n, m, tournament=False)) == 1)
    alphabet = tuple((n, m) for n in range(1, bound+1) for m in strong[n])
    words = {n: [] for n in range(1, bound+1)}
    def generate(prefix, size):
        if prefix and lyndon(prefix):
            words[size].append(prefix)
        for letter in alphabet:
            if size+letter[0] <= bound:
                generate(prefix+(letter,), size+letter[0])
    generate((), 0)
    forward, inverse = {}, {}
    for n in range(1, bound+1):
        ordered = tuple(sorted(words[n]))
        need(len(ordered) == len(connected[n]), 'independent atom-count equality')
        for w, m in zip(ordered, connected[n]):
            forward[w] = (n, m); inverse[n, m] = w
    return forward, inverse, strong, connected


def tournament_to_even(n, mask):
    source(n, mask)
    word = components(n, mask, tournament=True)
    table = _atlas(n)[0]
    return assemble(tuple(table[w] for w in factor(word)), tournament=False)


def even_to_tournament(n, mask):
    source(n, mask)
    need(automorphism_even(n, mask), 'automorphism-even graph required, not just Eulerian')
    table = _atlas(n)[1]
    words = tuple(sorted((table[p] for p in components(n, mask, tournament=False)), reverse=True))
    return assemble(tuple(letter for w in words for letter in w), tournament=True)


def decode_factors(parts):
    """A marked factor packet must actually be nonincreasing Lyndon factors."""
    need(type(parts) is tuple, 'exact factor packet')
    for w in parts:
        word_input(w); need(lyndon(w), 'nonempty Lyndon factor')
    need(all(a >= b for a, b in zip(parts, parts[1:])), 'retained nonincreasing factor order')
    word = tuple(letter for w in parts for letter in w); word_input(word)
    need(factor(word) == parts, 'complete factor packet')
    return assemble(word, tournament=True)


def section_to_graph(a, b):
    """The fixed path 0->1->2->3, flips02/13, with flip03 fixed at zero."""
    need(type(a) is int and a in (0, 1) and type(b) is int and b in (0, 1),
         'two exact section bits')
    slots = {e: i for i, e in enumerate(pairs(4))}
    return tournament_to_even(4, (a << slots[0, 2]) | (b << slots[1, 3]))


def graph_to_section(n, mask):
    source(n, mask); need(n == 4, 'the marked four-vertex section only')
    need(automorphism_even(n, mask), 'automorphism-even graph required')
    wanted = canonical(n, mask, tournament=False)
    solutions = tuple((a, b) for a, b in product((0, 1), repeat=2)
                      if section_to_graph(a, b)[1] == wanted)
    need(len(solutions) == 1, 'complete four-type section')
    return solutions[0]


def all_factorizations(word):
    if not word:
        return ((),)
    out = []
    for mask in range(1 << (len(word)-1)):
        start, blocks = 0, []
        for i in range(1, len(word)):
            if mask >> (i-1) & 1:
                blocks.append(word[start:i]); start = i
        blocks.append(word[start:]); blocks = tuple(blocks)
        if all(lyndon(w) for w in blocks) and all(a >= b for a, b in zip(blocks, blocks[1:])):
            out.append(blocks)
    return tuple(out)


def reject(call):
    try:
        call()
    except ValueError:
        check(True, 'malformed or false-domain input rejected')
    else:
        check(False, 'rejection required')


def main():
    table = []; relabels = 0
    for n in range(LIMIT+1):
        tours = _catalogue(n, True)[0]
        even = tuple(m for m in _catalogue(n, False)[0] if automorphism_even(n, m))
        forward, inverse, strong, connected = _atlas(n)
        check(len(tours) == len(even), 'complete type counts')
        image = set()
        for t in tours:
            g = tournament_to_even(n, t)
            check(automorphism_even(*g), 'image has the correct even predicate')
            check(even_to_tournament(*g) == (n, t), 'tournament inverse')
            word = components(n, t, tournament=True); factors = factor(word)
            check(all_factorizations(word) == (factors,), 'independent all-cut factor uniqueness')
            check(decode_factors(factors) == (n, t), 'marked factor inverse')
            check(sorted(weight(w) for w in factors) ==
                  sorted(s for s, _ in components(*g, tournament=False)), 'exact preserved weight multiset')
            image.add(g[1])
            # Every relabelling of each representative, through six vertices.
            for _, dst, flip, _ in actions(n):
                tt = act(t, dst, flip)
                check(tournament_to_even(n, tt) == g, 'forward factors through isomorphism type')
                relabels += 1
        check(image == set(even), 'onto every automorphism-even type')
        for g in even:
            t = even_to_tournament(n, g)
            check(tournament_to_even(*t) == (n, g), 'graph inverse')
            for _, dst, _, _ in actions(n):
                check(even_to_tournament(n, act(g, dst)) == t, 'inverse factors through graph type')
                relabels += 1
        table.append({'n': n, 'tournament_types': len(tours), 'even_types': len(even),
                      'strong_types': len(strong.get(n, ())), 'connected_even_types': len(connected.get(n, ()))})
    # The native fixed-Hamiltonian-path four-type transversal.
    slots = {e: i for i, e in enumerate(pairs(4))}
    four = []
    for a, b in product((0, 1), repeat=2):
        mask = (a << slots[0, 2]) | (b << slots[1, 3])
        g = section_to_graph(a, b)
        check(graph_to_section(*g) == (a, b), 'explicit marked two-bit inverse')
        parity = g[1].bit_count() % 2
        check(parity == (b & (1 ^ a))
              and (parity ^ (section_to_graph(1 ^ a, b)[1].bit_count() % 2)) == b
              and (parity ^ (section_to_graph(a, 1 ^ b)[1].bit_count() % 2)) == (1 ^ a),
              'quadratic interaction coordinate and context-dependent port derivatives')
        four.append({'a': a, 'b': b, 'tournament_mask': mask, 'graph_mask': g[1]})
    check([r['graph_mask'] for r in four] == [0, 7, 3, 30], 'explicit four-type graph table')
    check(len({canonical(4, r['tournament_mask'], tournament=True) for r in four}) == 4,
          'the partial-tournament section is a transversal')
    check(sum(r['graph_mask'].bit_count() % 2 for r in four) == 1,
          'four graph types have odd parity count1, impossible for an affine edge-toggle square')
    even4 = tuple(m for m in range(64) if automorphism_even(4, m))
    affine_transversals = 0
    for x, y, z in combinations(even4, 3):
        fourth = x ^ y ^ z
        if fourth in even4:
            types = {canonical(4, m, tournament=False) for m in (x, y, z, fourth)}
            affine_transversals += len(types) == 4
    check(affine_transversals == 0, 'all labelled affine-square candidates independently fail')
    # Losing chronology by sorting SCC letters is not our factorization map.
    a = (1, 0); b = (3, _atlas(3)[2][3][0])
    ab, ba = assemble((a, b), tournament=True), assemble((b, a), tournament=True)
    check(ab != ba and tournament_to_even(*ab) == (4, 7) and
          tournament_to_even(*ba) == (4, 3), 'ordered join is not sent to disjoint union')
    check(assemble((tournament_to_even(*a), tournament_to_even(*b)), tournament=False) == (4, 3),
          'disjoint union is the ba outcome, not the ab outcome')
    for r in range(1, 4):
        front, back = (a,) * r + (b,), (b,) + (a,) * r
        gf = tournament_to_even(*assemble(front, tournament=True))
        gb = tournament_to_even(*assemble(back, tournament=True))
        check(lyndon(front) and factor(back) == ((b,),) + ((a,),) * r
              and len(components(*gf, tournament=False)) == 1
              and gb == assemble(((3, 3),) + ((1, 0),) * r, tournament=False),
              'moving one C3 across r singleton letters fuses r+1 graph components')
    # Actual equivariance obstruction already for this chosen map at n=3.
    cyclic = b[1]; path = tournament_to_even(*b)[1]
    bad = None
    for p, dst, flip, _ in actions(3):
        if act(cyclic, dst, flip) == cyclic and act(path, dst) != path:
            bad = p; break
    check(bad is not None, 'chosen type codec is not a labelled equivariant map')
    # All conjugates of a C3 rotation likewise cannot act on a path of order3.
    check(all(sum(act(path, dst) == path for _, dst, _, _ in actions(3)) == 2 for _ in (0,)),
          'path automorphism order2 versus cycle order3')
    # One Eulerian triangle is sign-odd and outside the graph-side domain.
    reject(lambda: even_to_tournament(3, 7))
    # Primitive necklace and periodic word differ; one may not collapse aa.
    check(not lyndon((a, a)) and factor((a, a)) == ((a,), (a,)), 'periodicity retains multiplicity')
    reject(lambda: decode_factors(((a, a),)))
    reject(lambda: decode_factors(((a,), (b,))))
    # Type checks must run before cache lookups with bool/float aliases.
    tournament_to_even(1, 0)
    reject(lambda: tournament_to_even(True, 0))
    reject(lambda: tournament_to_even(1, False))
    reject(lambda: even_to_tournament(3.0, 3))
    reject(lambda: canonical(3, 0, tournament=1))
    reject(lambda: tournament_to_even(7, 0))
    reject(lambda: decode_factors((((1, False),),)))
    reject(lambda: section_to_graph(True, 0))
    reject(lambda: graph_to_section(4, 1))
    print(json.dumps({'scope': 'explicit inverse pair on types; atom-ranking is non-natural',
                      'checks': CHECKS, 'complete_orders': '0..6', 'type_counts': table,
                      'relabelled_roundtrip_controls': relabels, 'four_type_transversal': four,
                      'edge_parity_polynomial_F2': 'b + a*b',
                      'order_join_hostile': {'ab_graph_mask': 7, 'ba_graph_mask': 3},
                      'equivariance_hostile_permutation': bad}, indent=2))


if __name__ == '__main__':
    main()

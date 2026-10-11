"""Deletion-context defect antichains and exact finite-prime escape guards.

No seven-chromatic unit graph is constructed. Profiles retain the actual
source graph, palette and marked boundary. Number-field splitting tests
are sufficient screen-removal certificates, not chromatic lower bounds.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, product
from math import isqrt, prod
import boundary_palette_compression_20261011 as p
import boundary_context_weights_20261011d as w

CHECKS = 0


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def minimal(sets):
    ordered = sorted(set(sets), key=lambda s: (len(s), tuple(sorted(s))))
    answer = []
    for s in ordered:
        if not any(t <= s for t in answer):
            answer.append(s)
    return tuple(answer)


def defects(g, colours):
    p.audit_graph(g)
    p.need(type(colours) is tuple and len(colours) == len(g.vertices)
           and all(type(c) is int and c >= 0 for c in colours), 'literal colours')
    a = dict(zip(g.vertices, colours))
    return frozenset(e for e in g.edges if a[e[0]] == a[e[1]])


@dataclass(frozen=True)
class Profile:
    graph: p.Graph
    boundary: tuple
    palette: int
    entries: tuple


def profile(g, boundary, k):
    p.check_domain(g, boundary, k)
    positions = tuple(g.vertices.index(v) for v in boundary)
    result = {}
    for cs in p.patterns(len(g.vertices), k):
        pi = p.canonical(tuple(cs[i] for i in positions))
        result.setdefault(pi, set()).add(defects(g, cs))
    return Profile(g, boundary, k, tuple((pi, minimal(ss)) for pi, ss in sorted(result.items())))


def audit_profile(t):
    p.need(type(t) is Profile, 'typed profile')
    p.check_domain(t.graph, t.boundary, t.palette)
    p.need(type(t.entries) is tuple, 'immutable profile entries')
    for row in t.entries:
        p.need(type(row) is tuple and len(row) == 2, 'typed profile row')
        pi, bank = row
        p.need(type(pi) is tuple and len(pi) == len(t.boundary)
               and all(type(c) is int and c >= 0 for c in pi)
               and p.canonical(pi) == pi and len(set(pi)) <= t.palette,
               'exact canonical profile pattern')
        p.need(type(bank) is tuple, 'immutable defect bank')
        for d in bank:
            p.need(type(d) is frozenset
                   and all(type(e) is tuple and len(e) == 2
                           and all(type(v) is int for v in e) for e in d)
                   and d <= t.graph.edges, 'exact immutable source defects')
    # Complete reauthentication, not a checksum of the supplied entries.
    p.need(t == profile(t.graph, t.boundary, t.palette), 'complete source profile')


def allows(t, removed, boundary_pattern):
    audit_profile(t)
    p.need(type(removed) is frozenset and removed <= t.graph.edges
           and all(type(e) is tuple and len(e) == 2
                   and all(type(v) is int for v in e) for e in removed),
           'actual exact removed edge set')
    p.need(type(boundary_pattern) is tuple
           and len(boundary_pattern) == len(t.boundary)
           and all(type(c) is int and c >= 0 for c in boundary_pattern)
           and p.canonical(boundary_pattern) == boundary_pattern
           and len(set(boundary_pattern)) <= t.palette, 'canonical marked pattern')
    return any(d <= removed for d in dict(t.entries).get(boundary_pattern, ()))


def glue(t, u):
    audit_profile(t)
    audit_profile(u)
    p.need(t.palette == u.palette and t.boundary == u.boundary,
           'same ordered interface and palette')
    p.need(set(t.graph.vertices) & set(u.graph.vertices) == set(t.boundary),
           'private interiors; retained complete overlap')
    g = p.Graph(tuple(sorted(set(t.graph.vertices) | set(u.graph.vertices))),
                t.graph.edges | u.graph.edges)
    left, right = dict(t.entries), dict(u.entries)
    rows = tuple((pi, minimal(a | b for a, b in product(left[pi], right[pi])))
                 for pi in sorted(left.keys() & right.keys()))
    return Profile(g, t.boundary, t.palette, rows)


def mycielski(g):
    p.audit_graph(g)
    n = len(g.vertices)
    p.need(g.vertices == tuple(range(n)), 'canonical labels for clone construction')
    edges = set(g.edges)
    for a, b in g.edges:
        edges.add((a, b+n))
        edges.add((b, a+n))
    edges.update((v+n, 2*n) for v in range(n))
    return p.Graph(tuple(range(2*n+1)), frozenset(edges))


def prime(n):
    return type(n) is int and n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def prime_escape(primes):
    p.need(type(primes) is tuple and len(set(primes)) == len(primes)
           and all(prime(q) for q in primes), 'distinct actual rational primes')
    product_odd = prod(q for q in primes if q != 2)
    t = 2*product_odd
    return t, 4*t-1


def split_witness(d, q):
    p.need(type(d) is int and d > 0 and prime(q), 'imaginary radicand and prime')
    if q == 2:
        return 'one mod eight' if (-d) % 8 == 1 else None
    return next((a for a in range(1, q) if (a*a+d) % q == 0), None)


def split_profile(t, primes):
    p.need(type(t) is int and t >= 1, 'positive exact rotation parameter')
    p.need(type(primes) is tuple and all(prime(q) for q in primes), 'prime bank')
    ds = (3, 11, 4*t-1)
    return tuple((q, next(((d, split_witness(d, q)) for d in ds
                          if split_witness(d, q) is not None), None)) for q in primes)


def v2(n):
    p.need(type(n) is int and n > 0, 'positive integer')
    return (n & -n).bit_length()-1


def rejects(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


def main():
    print('Universe: labelled graphs on0..4 vertices, palettes2..4, every edge-deletion context.')
    cases = 0
    for n in range(5):
        for mask in range(1 << (n*(n-1)//2)):
            g = w.graph_from_mask(n, mask)
            edges = tuple(sorted(g.edges))
            for k in range(2, 5):
                t = profile(g, (), k)
                bank = dict(t.entries)[()]
                for deletion in range(1 << len(edges)):
                    removed = frozenset(e for i, e in enumerate(edges) if deletion >> i & 1)
                    fragment = p.Graph(g.vertices, g.edges-removed)
                    literal = any(p.valid_colouring(fragment, k, cs)
                                  for cs in product(range(k), repeat=n))
                    check(literal == any(d <= removed for d in bank), 'all deletion contexts')
                    cases += 1
                for removed in bank:
                    fragment = p.Graph(g.vertices, g.edges-removed)
                    for edge in removed:
                        check(w.conditional_count(fragment, edge, k, (0, 0)) > 0,
                              'minimal nonempty defect retains equal witness')
                        check(w.conditional_count(fragment, edge, k, (0, 1)) == 0,
                              'every minimal defect edge is forced equal')
    print('Independent literal deletion comparisons:', cases)

    # The marked interface is essential: a two-edge path and an edge are
    # individually2-colourable but their gluing is the odd triangle.
    path = p.Graph((0, 1, 2), frozenset({(0, 2), (1, 2)}))
    edge = p.Graph((0, 1), frozenset({(0, 1)}))
    joined = glue(profile(path, (0, 1), 2), profile(edge, (0, 1), 2))
    check(joined == profile(joined.graph, (0, 1), 2), 'literal marked antichain gluing')
    check(not allows(joined, frozenset(), (0, 0)), 'triangle equal boundary cannot extend')
    check(not allows(joined, frozenset(), (0, 1)), 'triangle different boundary cannot extend')
    print('Forgetting shared port patterns before gluing loses the triangle obstruction.')
    for k in range(2, 7):
        for a, b in product(range(8), repeat=2):
            left0, right0 = w.graph_from_mask(3, a), w.graph_from_mask(3, b)
            right = p.Graph((0, 1, 3), frozenset(tuple(sorted(3 if v == 2 else v for v in e))
                                               for e in right0.edges))
            out = glue(profile(left0, (0, 1), k), profile(right, (0, 1), k))
            check(out == profile(out.graph, (0, 1), k), 'all small marked gluing contexts')

    # A fixed colour assignment determines its strongest compatible fragment.
    spindle = p.spindle()
    for k in (3, 4, 6):
        bank = dict(profile(spindle, (), k).entries)[()]
        check((frozenset() in bank) == (k >= 4), 'Moser minimal-defect support')
        if k == 3:
            check(set(bank) == {frozenset({e}) for e in spindle.edges}, 'all11 minimal Moser repairs')
        print(f'Moser palette{k}: minimal defect supports={len(bank)}, sizes={sorted(set(map(len, bank)))}')
    complete = p.Graph(tuple(range(7)), frozenset(combinations(range(7), 2)))
    bank = dict(profile(complete, (), 6).entries)[()]
    check(set(bank) == {frozenset({e}) for e in complete.edges}, 'abstract six-palette positive control')
    print('K7 palette6:21 singleton minimal repairs; abstract, not planar unit geometry.')

    # Standard chromatic lifting is not automatically a unit-distance recursion.
    claw = p.Graph((0, 1, 2, 3), frozenset({(0, 1), (0, 2), (0, 3)}))
    mc = mycielski(claw)
    adj = {v: set() for v in mc.vertices}
    for a, b in mc.edges:
        adj[a].add(b)
        adj[b].add(a)
    check(adj[0] & adj[4] == {1, 2, 3}, 'clone keeps three common neighbours: forbidden unit geometry')
    print('Mycielski lift of any source with degree>=3 contains the unit-forbidden K2,3.')

    bank_primes = tuple(q for q in range(2, 20) if prime(q))
    for primes in ((2,), (2, 3, 11), bank_primes, tuple(q for q in range(2, 32) if prime(q))):
        t, d = prime_escape(primes)
        check(Q(2*t-1, 2*t)**2+d*Q(1, 2*t)**2 == 1, 'exact norm-one rotation')
        check(all(split_witness(d, q) is not None for q in primes), 'all requested primes split')
        # A compatible squared Eisenstein displacement uses m=2P diamonds.
        m = 2*prod(q for q in primes if q != 2)
        compatible_t, compatible_d = 3*m*m, 12*m*m-1
        check(all(split_witness(compatible_d, q) is not None for q in primes),
              'norm-compatible prime escape parameter')
        print('Prime escape:', primes, 'rotation t,D=', (t, d),
              'diamond family m,t,D=', (m, compatible_t, compatible_d))
    passing = [t for t in range(2, 121, 2) if all(cert is not None for _, cert in split_profile(t, bank_primes))]
    compatible = [t for t in passing if v2(t) % 2 == 0]
    check(passing[0] == 22 and compatible[0] == 60, 'bounded first screens; norm guard only necessary')
    print('Even t<=120 with a split-subfield certificate at every prime<=19:', passing)
    print('First passing t22 has odd2-adic norm valuation; first surviving that necessary guard t60.')
    print('t60 split certificates:', split_profile(60, bank_primes))
    check(33+3*3**2 == 60, 'explicit Moser displacement sqrt33+3sqrt-3')
    check(Q(119, 120)**2+239*Q(1, 120)**2 == 1, 'omega60 norm')
    check(60*((1-Q(119, 120))**2+239*Q(1, 120)**2) == 1, 'omega60 closes certified displacement')
    print('Certified geometric parameter: z=sqrt33+3sqrt(-3), norm60; omega=(119+sqrt(-239))/120.')
    print('No equality rigidity or high chromatic conclusion follows from the displacement.')

    t = profile(path, (0, 1), 2)
    forged = Profile(t.graph, t.boundary, t.palette, ())
    check(rejects(lambda: audit_profile(forged)), 'forged complete table')
    check(rejects(lambda: profile(path, (False, 1), 2)), 'bool vertex alias')
    check(rejects(lambda: profile(path, (0, 1), True)), 'bool palette alias')
    check(rejects(lambda: prime_escape((2, 9))), 'nonprime input')
    check(rejects(lambda: prime_escape((2, True))), 'bool prime alias')
    check(rejects(lambda: allows(t, frozenset({(0, 1)}), (0, 0))), 'foreign removed edge')
    te = profile(edge, (0, 1), 2)
    malformed = (
        (((False, False), (frozenset({(0, 1)}),)), ((0, 1), (frozenset(),))),
        (((0, 0), (frozenset({(0.0, 1.0)}),)), ((0, 1), (frozenset(),))),
        (((0, 0), ({(0, 1)},)), ((0, 1), (frozenset(),))),
        list(te.entries),
        (list(te.entries[0]), te.entries[1]),
        (((0, 0), list(te.entries[0][1])), te.entries[1]),
        (((0, 2), te.entries[0][1]), te.entries[1]),
        (((0, 0), (frozenset({(0, 2)}),)), te.entries[1]),
    )
    for rows in malformed:
        check(rejects(lambda: audit_profile(Profile(te.graph, te.boundary, te.palette, rows))),
              'exact packet types and source edges precede Python equality')
    check(rejects(lambda: allows(te, frozenset(), (False, True))), 'bool query pattern alias')
    print(f'PASS: {CHECKS} exact checks. Six-colour equality remains OPEN.')


if __name__ == '__main__':
    main()

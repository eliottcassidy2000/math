"""Exact odd-cover descent, signed fibres, and prime-power odd symbols.

Run with Python (also -B -O); JSON on stdout.  No graph discovery or metric
embedding is asserted.  Every cover retains its base, sheets and edge matchings.
"""
from dataclasses import dataclass
from itertools import combinations, permutations, product
from math import gcd, isqrt
import json


def integer(x, low=None):
    if type(x) is not int or (low is not None and x < low):
        raise ValueError("exact integer outside domain")
    return x


def prime(p):
    integer(p, 3)
    if p % 2 == 0 or any(p % d == 0 for d in range(3, isqrt(p) + 1, 2)):
        raise ValueError("odd prime required")
    return p


@dataclass(frozen=True)
class Graph:
    n: int
    edges: tuple


def audit_graph(g):
    if type(g) is not Graph:
        raise ValueError("Graph required")
    integer(g.n, 1)
    if type(g.edges) is not tuple:
        raise ValueError("tuple edges required")
    for e in g.edges:
        if type(e) is not tuple or len(e) != 2:
            raise ValueError("edge type")
        u, v = e
        integer(u, 0); integer(v, 0)
        if not u < v < g.n:
            raise ValueError("simple labelled graph required")
    if g.edges != tuple(sorted(set(g.edges))):
        raise ValueError("canonical edge list required")
    return g


def graph(n, edges):
    return audit_graph(Graph(n, tuple(sorted(tuple(sorted(e)) for e in edges))))


def audit_pins(g, pins):
    audit_graph(g)
    if type(pins) is not tuple:
        raise ValueError("tuple pins")
    seen = set()
    for pin in pins:
        if type(pin) is not tuple or len(pin) != 2:
            raise ValueError("pin type")
        v, c = pin
        integer(v, 0); integer(c, 0)
        if v >= g.n or c >= 6 or v in seen:
            raise ValueError("unique source vertex and six-colour pin required")
        seen.add(v)


def labels_and_pins(g, labels, pins):
    audit_pins(g, pins)
    if type(labels) is not tuple or len(labels) != g.n:
        raise ValueError("one retained ternary label per base vertex")
    for a in labels:
        integer(a, 0)
        if a >= 3:
            raise ValueError("ternary label")


def audit_colouring(g, colours, labels=None, pins=()):
    audit_pins(g, pins)
    if type(colours) is not tuple or len(colours) != g.n:
        raise ValueError("one exact colour per vertex")
    for c in colours:
        integer(c, 0)
        if c >= 6:
            raise ValueError("six-colour range")
    if any(colours[u] == colours[v] for u, v in g.edges):
        raise ValueError("monochromatic actual edge")
    if labels is not None:
        labels_and_pins(g, labels, pins)
        if any(c % 3 != a for c, a in zip(colours, labels)):
            raise ValueError("fixed ternary source labels lost")
    for v, c in pins:
        if colours[v] != c:
            raise ValueError("source pin lost")
    return colours


def fixed_fibre(g, labels, pins=()):
    """Exact count and one witness for this fixed-label, source-pinned fibre."""
    labels_and_pins(g, labels, pins)
    if any(c % 3 != labels[v] for v, c in pins):
        return 0, None
    adj = [[] for _ in range(g.n)]
    for u, v in g.edges:
        if labels[u] == labels[v]:
            adj[u].append(v); adj[v].append(u)
    component = [-1] * g.n
    parity = [0] * g.n
    roots = []
    for root in range(g.n):
        if component[root] >= 0:
            continue
        ci = len(roots); roots.append(root)
        component[root] = ci
        todo = [root]
        for u in todo:
            for v in adj[u]:
                if component[v] < 0:
                    component[v] = ci; parity[v] = 1 - parity[u]; todo.append(v)
                elif parity[u] == parity[v]:
                    return 0, None
    fixed = {}
    for v, c in pins:
        b = (c % 2) ^ parity[v]
        ci = component[v]
        if ci in fixed and fixed[ci] != b:
            return 0, None
        fixed[ci] = b
    colours = tuple((3 * (parity[v] ^ fixed.get(component[v], 0)) + 4 * labels[v]) % 6
                    for v in range(g.n))
    audit_colouring(g, colours, labels, pins)
    return 2 ** (len(roots) - len(fixed)), colours


@dataclass(frozen=True)
class Cover:
    base: Graph
    sheets: int
    matchings: tuple


def audit_cover(c):
    if type(c) is not Cover:
        raise ValueError("source-owned Cover required")
    audit_graph(c.base); integer(c.sheets, 1)
    if type(c.matchings) is not tuple or len(c.matchings) != len(c.base.edges):
        raise ValueError("one matching per source edge")
    for p in c.matchings:
        if type(p) is not tuple or len(p) != c.sheets:
            raise ValueError("matching arity")
        for x in p:
            integer(x, 0)
        if sorted(p) != list(range(c.sheets)):
            raise ValueError("cover needs bijective edge matchings")
    return c


def lift_graph(c):
    audit_cover(c)
    return graph(c.base.n * c.sheets,
                 [(u * c.sheets + i, v * c.sheets + p[i])
                  for (u, v), p in zip(c.base.edges, c.matchings)
                  for i in range(c.sheets)])


def lift_context(c, labels, pins=()):
    audit_cover(c); labels_and_pins(c.base, labels, pins)
    return (tuple(a for a in labels for _ in range(c.sheets)),
            tuple((v * c.sheets + i, a) for v, a in pins for i in range(c.sheets)))


def pullback(c, labels, pins, colours):
    audit_cover(c); labels_and_pins(c.base, labels, pins)
    audit_colouring(c.base, colours, labels, pins)
    result = tuple(a for a in colours for _ in range(c.sheets))
    ll, pp = lift_context(c, labels, pins)
    audit_colouring(lift_graph(c), result, ll, pp)
    return result


def majority_descend(c, labels, pins, lifted_colours):
    """Odd-cover support descent. Every pin is replicated to its entire fibre."""
    audit_cover(c); labels_and_pins(c.base, labels, pins)
    if c.sheets % 2 == 0:
        raise ValueError("odd sheet count is required")
    ll, pp = lift_context(c, labels, pins)
    audit_colouring(lift_graph(c), lifted_colours, ll, pp)
    sums = tuple(sum(1 - 2 * (lifted_colours[v * c.sheets + i] % 2)
                     for i in range(c.sheets)) for v in range(c.base.n))
    if any(s == 0 for s in sums):
        raise ArithmeticError("impossible odd fibre tie")
    result = tuple((3 * int(s < 0) + 4 * labels[v]) % 6 for v, s in enumerate(sums))
    audit_colouring(c.base, result, labels, pins)
    return result, sums


def audit_thickening(c, host):
    """A matching skeleton plus only edges over genuine, nonloop base edges."""
    audit_cover(c); audit_graph(host)
    if host.n != c.base.n * c.sheets:
        raise ValueError("host changed the source fibres")
    old = set(c.base.edges); actual = set(host.edges)
    if not set(lift_graph(c).edges) <= actual:
        raise ValueError("a retained perfect matching edge is missing")
    if any(tuple(sorted((u // c.sheets, v // c.sheets))) not in old for u, v in host.edges):
        raise ValueError("new projected pair or collapsed edge")
    return host


def thickening_descend(c, host, labels, pins, colours):
    audit_thickening(c, host)
    ll, pp = lift_context(c, labels, pins)
    audit_colouring(host, colours, ll, pp)
    return majority_descend(c, labels, pins, colours)


def thickening_pullback(c, host, labels, pins, colours):
    audit_thickening(c, host)
    result = pullback(c, labels, pins, colours)
    ll, pp = lift_context(c, labels, pins)
    return audit_colouring(host, result, ll, pp)


def odd_cycle_lift(c, cycle):
    """Authenticate an actual odd simple cycle in the lift using monodromy."""
    audit_cover(c)
    if c.sheets % 2 == 0 or type(cycle) is not tuple or len(cycle) < 3 or len(cycle) % 2 == 0:
        raise ValueError("odd sheets and an odd simple base cycle required")
    for v in cycle:
        integer(v, 0)
        if v >= c.base.n:
            raise ValueError("cycle source")
    if len(set(cycle)) != len(cycle):
        raise ValueError("base cycle must be simple")
    moves = {}
    for (u, v), p in zip(c.base.edges, c.matchings):
        inv = tuple(p.index(i) for i in range(c.sheets))
        moves[u, v] = p; moves[v, u] = inv
    steps = tuple(zip(cycle, cycle[1:] + cycle[:1]))
    if any(e not in moves for e in steps):
        raise ValueError("cycle uses a nonexistent source edge")
    def tour(i):
        for e in steps:
            i = moves[e][i]
        return i
    seen = set(); chosen = None
    for start in range(c.sheets):
        if start in seen:
            continue
        orbit = []; i = start
        while i not in seen:
            seen.add(i); orbit.append(i); i = tour(i)
        if len(orbit) % 2 == 1 and chosen is None:
            chosen = orbit
    if chosen is None:
        raise ArithmeticError("odd permutation has an odd orbit")
    result = []; i = chosen[0]
    for _ in chosen:
        for u, v in steps:
            result.append(u * c.sheets + i); i = moves[u, v][i]
    lifted = lift_graph(c)
    es = set(lifted.edges)
    if (i != chosen[0] or len(set(result)) != len(result) or len(result) % 2 == 0
            or any(tuple(sorted(e)) not in es for e in zip(result, result[1:] + result[:1]))):
        raise ArithmeticError("invalid lifted cycle")
    return tuple(result)


def translation_count(domain, step, modulus):
    integer(domain, 1); integer(step); integer(modulus, 2)
    orbits = gcd(domain, step)
    length = domain // orbits
    return (modulus if length % 2 == 0 else gcd(modulus, 2)) ** orbits


def inversion_count(domain, modulus):
    integer(domain, 1); integer(modulus, 2)
    fixed = gcd(domain, 2)
    return gcd(modulus, 2) ** fixed * modulus ** ((domain - fixed) // 2)


def legendre(u, p):
    prime(p); integer(u)
    a = pow(u % p, (p - 1) // 2, p)
    return -1 if a == p - 1 else a


def square_lift(p, r, u, root):
    """The unique square root extending the supplied nonzero residue root."""
    prime(p); integer(r, 1); integer(u); integer(root, 1)
    if root >= p or u % p == 0 or (root * root - u) % p:
        raise ValueError("native nonzero square-root branch required")
    x = root; modulus = p
    for _ in range(1, r):
        digit = ((u - x * x) // modulus * pow(2 * x, -1, p)) % p
        x += modulus * digit; modulus *= p
    if (x * x - u) % modulus:
        raise ArithmeticError("square lifting identity")
    return x


def quadratic_symbol(p, r, layer_signs):
    """All quadratic-unit-covariant signs on the nonzero residues of Z/p^r."""
    prime(p); integer(r, 1)
    if type(layer_signs) is not tuple or len(layer_signs) != r:
        raise ValueError("one sign per valuation layer required")
    if any(type(s) is not int or s not in (-1, 1) for s in layer_signs):
        raise ValueError("exact signed layer frame required")
    result = [0]
    for d in range(1, p ** r):
        v = 0; a = d
        while a % p == 0:
            v += 1; a //= p
        result.append(layer_signs[v] * legendre(a, p))
    return tuple(result)


def unsigned_rows(g):
    audit_graph(g)
    return tuple(tuple(int(i in e) for i in range(g.n)) for e in g.edges)


def rank_mod(rows, n, p):
    prime(p); integer(n, 1)
    if type(rows) is not tuple:
        raise ValueError("matrix tuple required")
    for row in rows:
        if type(row) is not tuple or len(row) != n:
            raise ValueError("matrix shape")
        for x in row:
            integer(x)
    a = [[x % p for x in row] for row in rows]
    rank = 0
    for col in range(n):
        pivot = next((j for j in range(rank, len(a)) if a[j][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, p)
        a[rank] = [x * inv % p for x in a[rank]]
        for j in range(len(a)):
            if j != rank and a[j][col]:
                factor = a[j][col]
                a[j] = [(x - factor * y) % p for x, y in zip(a[j], a[rank])]
        rank += 1
    return rank


def component_types(g):
    audit_graph(g)
    adj = [[] for _ in range(g.n)]
    for u, v in g.edges:
        adj[u].append(v); adj[v].append(u)
    colours = {}; good = bad = 0
    for r in range(g.n):
        if r in colours:
            continue
        colours[r] = 0; todo = [r]; bip = True
        for u in todo:
            for v in adj[u]:
                if v not in colours:
                    colours[v] = 1 - colours[u]; todo.append(v)
                elif colours[v] == colours[u]:
                    bip = False
        good += bip; bad += not bip
    return good, bad


def kernel_count(g, modulus):
    integer(modulus, 2)
    good, bad = component_types(g)
    return modulus ** good * gcd(modulus, 2) ** bad


def cubic_frame(coefficients):
    """Exact square and norm in Q[t]/(t^3-2), retaining the basis 1,t,t^2."""
    if type(coefficients) is not tuple or len(coefficients) != 3:
        raise ValueError("three exact source coefficients")
    for x in coefficients:
        integer(x)
    x, y, z = coefficients
    return ((x*x + 4*y*z, 2*x*y + 2*z*z, y*y + 2*x*z),
            x**3 + 2*y**3 + 4*z**3 - 6*x*y*z)


def verify():
    checks = 0; rejected = 0
    def check(condition):
        nonlocal checks
        if not condition:
            raise ArithmeticError("exact control failed")
        checks += 1
    def reject(f):
        nonlocal rejected
        try:
            f()
        except ValueError:
            rejected += 1; check(True); return
        raise ArithmeticError("malformed input was accepted")

    # Direct enumeration checks the all-parameter permutation formulas.
    for m in range(1, 7):
        for q in (2, 3, 4):
            vectors = tuple(product(range(q), repeat=m))
            check(sum(all((v[i] + v[-i % m]) % q == 0 for i in range(m))
                      for v in vectors) == inversion_count(m, q))
            for step in range(m):
                check(sum(all((v[i] + v[(i + step) % m]) % q == 0 for i in range(m))
                          for v in vectors) == translation_count(m, step, q))

    graph_count = 0
    for n in range(1, 6):
        possible = tuple(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            g = Graph(n, tuple(e for i, e in enumerate(possible) if mask >> i & 1))
            b, _ = component_types(g); graph_count += 1
            for p in (3, 5, 7):
                check(rank_mod(unsigned_rows(g), n, p) == n - b)
            if n <= 4:
                for q in (2, 3, 4, 6):
                    check(sum(all((x[u] + x[v]) % q == 0 for u, v in g.edges)
                              for x in product(range(q), repeat=n)) == kernel_count(g, q))
                # Every ternary label, with no pins or a fixed colour at 0.
                for labels in product(range(3), repeat=n):
                    for pins in ((), ((0, 0),), ((0, 1),)):
                        count, witness = fixed_fibre(g, labels, pins)
                        brute = 0
                        for bits in product(range(2), repeat=n):
                            cs = tuple((3 * b + 4 * a) % 6 for b, a in zip(bits, labels))
                            brute += (all(cs[u] != cs[v] for u, v in g.edges)
                                      and all(cs[v] == c for v, c in pins))
                        check(count == brute)
                        if witness is not None:
                            check(audit_colouring(g, witness, labels, pins) == witness)

    c3 = graph(3, ((0, 1), (1, 2), (0, 2)))
    c4 = graph(4, ((0, 1), (1, 2), (2, 3), (0, 3)))
    perms3 = tuple(permutations(range(3)))
    cycle_lengths = {}
    for pp in product(perms3, repeat=3):
        c = Cover(c3, 3, pp)
        cy = odd_cycle_lift(c, (0, 1, 2))
        cycle_lengths[len(cy)] = cycle_lengths.get(len(cy), 0) + 1
        check(fixed_fibre(lift_graph(c), (0,) * 9)[0] == 0)
        check(len(cy) in (3, 9))
        for labels in product(range(3), repeat=3):
            for pins in ((), ((0, 0),), ((0, 0), (1, 3))):
                bc, _ = fixed_fibre(c3, labels, pins)
                ll, lp = lift_context(c, labels, pins)
                hc, hw = fixed_fibre(lift_graph(c), ll, lp)
                check(bool(bc) == bool(hc))
                if hw is not None:
                    down, _ = majority_descend(c, labels, pins, hw)
                    check(audit_colouring(c3, down, labels, pins) == down)
    cover_tests = 0
    # All 6^4 edge-matchings of the three-sheet square cover.
    for pp in product(perms3, repeat=4):
        c = Cover(c4, 3, pp); gg = lift_graph(c)
        labels = (0, 0, 0, 0); pins = ((0, 0),)
        ll, lp = lift_context(c, labels, pins)
        count, witness = fixed_fibre(gg, ll, lp)
        check(count > 0)
        down, sums = majority_descend(c, labels, pins, witness)
        check(down == (0, 3, 0, 3) and all(abs(s) in (1, 3) for s in sums))
        check(pullback(c, labels, pins, down) == tuple(a for a in down for _ in range(3)))
        cover_tests += 1
    # Every shade colouring of a disconnected positive cover, including nonconstant fibres.
    trivial = Cover(graph(2, ((0, 1),)), 3, ((0, 1, 2),))
    for gauges in product(range(2), repeat=3):
        cc = tuple(3 * x for x in gauges) + tuple(3 * (1 - x) for x in gauges)
        down, sums = majority_descend(trivial, (0, 0), (), cc)
        check(down[0] != down[1] and sums[0] == -sums[1])
    # Boolean support is preserved, counts are not: 2 versus 2^3.
    check(fixed_fibre(trivial.base, (0, 0))[0] == 2)
    check(fixed_fibre(lift_graph(trivial), (0,) * 6)[0] == 8)
    dense = graph(6, ((i, 3 + j) for i in range(3) for j in range(3)))
    dc = thickening_pullback(trivial, dense, (0, 0), (), (0, 3))
    check(thickening_descend(trivial, dense, (0, 0), (), dc)[0] == (0, 3))
    check(fixed_fibre(dense, (0,) * 6)[0] == 2)
    reject(lambda: audit_thickening(trivial, graph(6, dense.edges + ((0, 1),))))
    reject(lambda: audit_thickening(trivial, graph(6, dense.edges[1:])))
    # One pin per fibre is insufficient; uniformly replicated pins are necessary.
    cc = (0, 3, 0, 3, 0, 3)
    check(cc[0] == cc[4] == 0 and all(cc[u] != cc[v] for u, v in lift_graph(trivial).edges))
    reject(lambda: majority_descend(trivial, (0, 0), ((0, 0), (1, 0)), cc))
    double = Cover(c3, 2, ((1, 0),) * 3)
    double_colours = (0, 3) * 3
    check(audit_colouring(lift_graph(double), double_colours) == double_colours)
    reject(lambda: majority_descend(double, (0, 0, 0), (), double_colours))
    # Unrestricted ternary labels destroy odd-cover descent.
    k7 = graph(7, combinations(range(7), 2))
    k7cover = Cover(k7, 3, ((1, 2, 0),) * 21)
    sheet_colours = (0, 1, 2) * 7
    check(audit_colouring(lift_graph(k7cover), sheet_colours) == sheet_colours)
    reject(lambda: majority_descend(k7cover, (0,) * 7, (), sheet_colours))

    symbol_rows = []
    for p, r in ((3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 1), (7, 2)):
        m = p ** r; good = 0
        for signs in product((-1, 1), repeat=r):
            sig = quadratic_symbol(p, r, signs)
            check(all(sig[u * d % m] == legendre(u, p) * sig[d]
                      for u in range(1, m) if u % p for d in range(1, m)))
            odd = all(sig[-d % m] == -sig[d] for d in range(1, m))
            check(odd == (p % 4 == 3)); good += odd
        check(good == (2 ** r if p % 4 == 3 else 0))
        symbol_rows.append([p, r, good])
        for u in range(1, m):
            if u % p == 0:
                continue
            roots = [a for a in range(1, p) if (a * a - u) % p == 0]
            for a in roots:
                z = square_lift(p, r, u, a)
                check(z % p == a and z * z % m == u)
            check(len(roots) == (2 if legendre(u, p) == 1 else 0))
        if m <= 9:
            actual = set()
            for half in product((-1, 1), repeat=(m - 1) // 2):
                ss = [0] * m
                for d, val in enumerate(half, 1):
                    ss[d] = val; ss[m - d] = -val
                if all(ss[u * d % m] == legendre(u, p) * ss[d]
                       for u in range(1, m) if u % p for d in range(1, m)):
                    actual.add(tuple(ss))
            predicted = {quadratic_symbol(p, r, s) for s in product((-1, 1), repeat=r)} if p % 4 == 3 else set()
            check(actual == predicted)
    check(legendre(3, 3) == 0 and quadratic_symbol(3, 2, (1, 1))[3] == 1)
    # Modulo 3, Q=B^T B kills 1 on K4, but B does not.
    k4 = graph(4, combinations(range(4), 2)); rows = unsigned_rows(k4)
    gram = tuple(tuple(sum(r[i] * r[j] for r in rows) for j in range(4)) for i in range(4))
    check(rank_mod(rows, 4, 3) == 4 and rank_mod(gram, 4, 3) == 3)
    check(all(sum(r) % 3 == 0 for r in gram) and all(sum(r) % 3 == 2 for r in rows))
    check(kernel_count(c3, 6) == 2 and all((3 + 3) % 6 == 0 for _ in c3.edges))
    # Norm descent in a marked cubic field, with determinant checked independently.
    for x, y, z in product(range(-3, 4), repeat=3):
        square, nb = cubic_frame((x, y, z))
        matrix = ((x, 2*z, 2*y), (y, x, 2*z), (z, y, x))
        det = sum((1 if i in (0, 3, 4) else -1) *
                  matrix[0][p[0]] * matrix[1][p[1]] * matrix[2][p[2]]
                  for i, p in enumerate(permutations(range(3))))
        check(nb == det)
        check((square[1:] == (0, 0)) == (y == z == 0))
        if square[1:] == (0, 0) and square[0] != 0:
            a = square[0]
            check(nb * nb == a ** 3 and nb == x * a)
    # Exact types, wrong source matchings, pin labels and native root branches.
    reject(lambda: audit_graph(Graph(True, ())))
    reject(lambda: audit_graph(Graph(2, ((0, 1.0),))))
    reject(lambda: audit_cover(Cover(c3, 3, ((0, 1, 1),) * 3)))
    reject(lambda: audit_cover(Cover(c3, 3, ((False, 1, 2),) * 3)))
    reject(lambda: audit_colouring(Graph(1, ()), (0,), pins=((False, 0),)))
    reject(lambda: audit_colouring(Graph(1, ()), (0,), pins=((0, 0.0),)))
    reject(lambda: quadratic_symbol(3, 2, (True, 1)))
    reject(lambda: square_lift(3, 2, 2, 1))
    reject(lambda: majority_descend(trivial, (0.0, 0), (), (0, 0, 0, 3, 3, 3)))
    reject(lambda: odd_cycle_lift(Cover(c3, 3, ((0, 1, 2),) * 3), (0, 1, 0)))
    check(fixed_fibre(Graph(1, ()), (1,), ((0, 0),)) == (0, None))
    return {
        "status": "FINITE-EXACT controls for stated all-parameter elementary theorems",
        "checks": checks,
        "universe": {"all_simple_graphs_through_5": graph_count,
                     "odd_characteristics": [3, 5, 7],
                     "fixed_label_bruteforce": "all graphs n<=4, every ternary label, no pins or vertex 0 pinned 0/1",
                     "three_sheet_triangle_covers": 216,
                     "triangle_cover_contexts": "every ternary label and three prescribed pin contexts per cover",
                     "lifted_odd_cycle_length_histogram": cycle_lengths,
                     "three_sheet_square_covers": cover_tests,
                     "literal_permutation_kernels": "domain 1..6, codomain Z/2,Z/3,Z/4"},
        "norm_descent_control": "all343 integer coefficient vectors in[-3,3]^3 in the marked basis of Q(cuberoot2)",
        "quadratic_covariant_odd_symbol_counts": symbol_rows,
        "positive": {"base_square": c4.edges, "uniform_pin": [0, 0],
                     "base_colouring": [0, 3, 0, 3], "every_three_sheet_edge_matching": True,
                     "mod9_symbol": quadratic_symbol(3, 2, (1, 1))},
        "hostiles": {"unrestricted_odd_cover": {"base": "K7", "sheets": 3,
                     "every_edge_matching": [1, 2, 0], "colour": "sheet index"},
                     "even_cover": "C3 double cover is C6",
                     "support_not_counts": [2, 8],
                     "modular_Gram": {"graph": "K4", "p": 3, "rank_B": 4, "rank_Q": 3},
                     "nonunit_triangle_mod6": [3, 3, 3],
                     "one_pin_per_fibre": cc, "type_or_guard_rejections": rejected},
        "scope": "Fixed ternary labels and uniformly replicated pins; odd covers and authenticated matching thickenings; no metric realization or new HN colouring bound."
    }


if __name__ == "__main__":
    print(json.dumps(verify(), sort_keys=True, indent=2))

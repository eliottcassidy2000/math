"""Exact marked digit, tournament and colour-difference interfaces.

No floating point, discovery of a geometric equality gadget, or reliance on
idoneal completeness. Running this file emits its deterministic JSON receipt.
"""
from itertools import combinations, permutations, product
from math import gcd, isqrt
import json


def integer(x, minimum=None):
    if type(x) is not int or (minimum is not None and x < minimum):
        raise ValueError("exact integer required")
    return x


def residue(x, q):
    integer(x)
    if not 0 <= x < q:
        raise ValueError("noncanonical residue")
    return x


def graph(n, edges):
    integer(n, 1)
    if type(edges) is not tuple:
        raise ValueError("ordered tuple of edges required")
    for e in edges:
        if type(e) is not tuple or len(e) != 2:
            raise ValueError("bad edge")
        u, v = e
        integer(u); integer(v)
        if not 0 <= u < v < n:
            raise ValueError("canonical simple edge required")
    if tuple(sorted(set(edges))) != edges:
        raise ValueError("duplicate/unordered edges")
    adj = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i, 1)); adj[v].append((u, i, -1))
    seen = {0}; queue = [0]
    for u in queue:
        for v, _, _ in adj[u]:
            if v not in seen:
                seen.add(v); queue.append(v)
    if len(seen) != n:
        raise ValueError("connected marked graph required")
    return adj


def tree_code(n, edges):
    """Return tree-edge indices, vertex path rows, and fundamental cycle rows."""
    adj = graph(n, edges)
    m = len(edges); paths = [None] * n; paths[0] = (0,) * m
    queue = [0]; tree = []
    for u in queue:
        for v, i, sign in adj[u]:
            if paths[v] is None:
                row = list(paths[u]); row[i] += sign
                paths[v] = tuple(row); tree.append(i); queue.append(v)
    cycles = []
    for i, (u, v) in enumerate(edges):
        if i not in tree:
            row = [paths[v][j] - paths[u][j] for j in range(m)]
            row[i] -= 1
            cycles.append(tuple(row))
    return tuple(tree), tuple(paths), tuple(cycles)


def decode(n, edges, q, gains, root=0):
    """Exact source-bound potential, or None for inconsistent cycle gains."""
    integer(q, 2); residue(root, q)
    if type(gains) is not tuple or len(gains) != len(edges):
        raise ValueError("one marked gain per edge required")
    for x in gains:
        residue(x, q)
    _, paths, cycles = tree_code(n, edges)
    if any(sum(a*b for a, b in zip(row, gains)) % q for row in cycles):
        return None
    return tuple((root + sum(a*b for a, b in zip(row, gains))) % q
                 for row in paths)


def encode(n, edges, q, colours):
    graph(n, edges); integer(q, 2)
    if type(colours) is not tuple or len(colours) != n:
        raise ValueError("marked colours required")
    for c in colours:
        residue(c, q)
    return tuple((colours[v] - colours[u]) % q for u, v in edges)


def matrix_box(q, rows, lists):
    """Iff EVERY point of the nonempty Cartesian box lies in ker(rows)."""
    integer(q, 2)
    if type(rows) is not tuple or type(lists) is not tuple:
        raise ValueError("literal matrix and lists required")
    for ss in lists:
        if type(ss) is not tuple or not ss or tuple(sorted(set(ss))) != ss:
            raise ValueError("nonempty canonical list required")
        for x in ss:
            residue(x, q)
    for row in rows:
        if type(row) is not tuple or len(row) != len(lists):
            raise ValueError("matrix dimensions")
        for x in row:
            integer(x)
    base = tuple(ss[0] for ss in lists)
    base_ok = all(sum(a*b for a, b in zip(row, base)) % q == 0
                  for row in rows)
    independent_ok = all(a*(x-ss[0]) % q == 0
                         for row in rows for a, ss in zip(row, lists)
                         for x in ss)
    return base_ok and independent_ok


def digit_word(res, depth=2):
    """Descending ordered three-lift fibre; residues, not decimal concatenation."""
    integer(depth, 1); modulus = 3**depth; residue(res, modulus)
    step = modulus // 3
    return tuple((res-j*step) % modulus for j in range(3))


def lift_step(depth, lower, phase):
    integer(depth, 1); modulus = 3**(depth-1)
    residue(lower, modulus); residue(phase, 3)
    carry, newlower = divmod(4*lower+1, modulus)
    return newlower, (phase+carry) % 3, carry


def lift_inverse(depth, lower, phase):
    integer(depth, 1); modulus = 3**(depth-1)
    residue(lower, modulus); residue(phase, 3)
    previous = ((lower-1)*pow(4, -1, modulus)) % modulus if modulus > 1 else 0
    carry = (4*previous+1) // modulus
    return previous, (phase-carry) % 3


CONNECTION = (1, 3, 4, 7)


def arc(u, v):
    residue(u, 9); residue(v, 9)
    return (v-u) % 9 in CONNECTION


def reduced_forms(discriminant):
    integer(discriminant)
    if discriminant >= 0 or discriminant % 4 not in (0, 1):
        raise ValueError("negative quadratic discriminant required")
    forms = []
    for a in range(1, isqrt((-discriminant)//3)+1):
        for b in range(-a, a+1):
            if (b*b-discriminant) % (4*a):
                continue
            c = (b*b-discriminant)//(4*a)
            if a > c or gcd(gcd(a, b), c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            forms.append((a, b, c))
    return tuple(forms)


def genus120(m):
    integer(m, 1)
    if gcd(m, 120) != 1:
        raise ValueError("genus evaluation requires a represented 120-unit")
    return (int(m % 8 not in (1, 7)), int(m % 3 != 1))


def rotate_genus(pair):
    if type(pair) is not tuple or len(pair) != 2:
        raise ValueError("marked binary pair required")
    for x in pair:
        residue(x, 2)
    a, b = pair
    return b, a ^ b


def edge_colour(edge):
    """K4 on F2^2, origin0 fixed: direction pair then incidence-with0 shade."""
    if type(edge) is not tuple or len(edge) != 2:
        raise ValueError("literal edge pair")
    u, v = edge
    residue(u, 4); residue(v, 4)
    if not u < v:
        raise ValueError("canonical unordered K4 edge")
    return 2*((u ^ v)-1)+int(u != 0)


def vertex_palette_action(perm):
    if type(perm) is not tuple or len(perm) != 4:
        raise ValueError("four marked vertex images")
    for x in perm:
        residue(x, 4)
    if len(set(perm)) != 4:
        raise ValueError("vertex bijection required")
    action = [None]*6
    for edge in combinations(range(4), 2):
        image = tuple(sorted(perm[x] for x in edge))
        action[edge_colour(edge)] = edge_colour(image)
    return tuple(action)


def verify():
    checks = 0

    def check(condition):
        nonlocal checks
        if not condition:
            raise RuntimeError("exact control failed")
        checks += 1

    def rejects(fn):
        try:
            fn()
        except ValueError:
            check(True)
        else:
            check(False)

    # All digit residues through depth seven, with exact inverse and carry.
    for depth in range(1, 8):
        modulus = 3**depth; step = modulus//3
        for x in range(modulus):
            lower, phase = x % step, x // step
            rr, tt, carry = lift_step(depth, lower, phase)
            check(rr + step*tt == (4*x+1) % modulus)
            check(lift_inverse(depth, rr, tt) == (lower, phase))
            k = 3**(depth-1)
            # Compute the divided geometric sum at precision 3*modulus.
            factor = pow(4, k, 3*modulus)
            check((factor*x+(factor-1)//3) % modulus == (x+step) % modulus)
    words = tuple(int(''.join(str(x or 9) for x in digit_word(a % 9))) for a in range(1, 10))
    check(set(words) == {396, 417, 528, 639, 741, 852, 963, 174, 285})
    for u in range(9):
        for v in range(9):
            r, t = u % 3, u // 3; s, z = v % 3, v // 3
            check(arc(u, v) == ((s-r) % 3 == 1 if r != s else (z-t) % 3 == 1))
            check(arc(u, v) == arc((4*u+1) % 9, (4*v+1) % 9))
    automorphisms = []
    # Intrinsic 3-cycle counts identify exactly the three fibre blocks.
    for triple in combinations(range(9), 3):
        cyclic = all(sum(arc(u, v) for v in triple if v != u) == 1 for u in triple)
        if cyclic:
            check(len({x % 3 for x in triple}) in (1, 3))
    for u, v in combinations(range(9), 2):
        cyclic_count = 0
        for w in range(9):
            if w in (u, v):
                continue
            triple = (u, v, w)
            cyclic_count += all(sum(arc(a, b) for b in triple if b != a) == 1
                                for a in triple)
        check(cyclic_count == (1 if u % 3 == v % 3 else 3))
    # Independent complete 9! permutation audit: 81 actual automorphisms.
    for perm in permutations(range(9)):
        if all(arc(u, v) == arc(perm[u], perm[v]) for u in range(9) for v in range(u+1, 9)):
            automorphisms.append(perm)
    check(len(automorphisms) == 81)
    for perm in automorphisms:
        check(all(len({perm[r+3*t] % 3 for t in range(3)}) == 1 for r in range(3)))

    # Exact inverse for all connected graphs <=4 and all rooted q=2,3,4 colourings.
    graph_count = 0
    for n in range(1, 5):
        pairs = tuple(combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            edges = tuple(e for i, e in enumerate(pairs) if mask >> i & 1)
            try:
                tree, paths, cycles = tree_code(n, edges)
            except ValueError:
                continue
            graph_count += 1
            check(len(tree) == n-1 and len(cycles) == len(edges)-n+1)
            for q in (2, 3, 4):
                for tail in product(range(q), repeat=n-1):
                    cols = (0,)+tail; gains = encode(n, edges, q, cols)
                    check(decode(n, edges, q, gains) == cols)
                    check(all(gains) == all(cols[u] != cols[v] for u, v in edges))

    # Exhaustive Cartesian box theorem on triangle plus bridge, q=3.
    edges = ((0, 1), (0, 2), (1, 2), (2, 3))
    rows = tree_code(4, edges)[2]
    subsets = tuple(tuple(x for x in range(3) if m >> x & 1) for m in range(1, 8))
    for lists in product(subsets, repeat=4):
        literal = all(decode(4, edges, 3, tuple(gains)) is not None for gains in product(*lists))
        check(matrix_box(3, rows, tuple(lists)) == literal)
        simple = all(len(lists[i]) == 1 for i in range(3)) and decode(
            4, edges, 3, tuple(ss[0] for ss in lists)) is not None
        check(literal == simple)
    # Repeated-label annihilator: varying a gain can be killed by cycle coefficient3.
    check(matrix_box(6, ((3,),), ((2, 4),)))
    check(not matrix_box(6, ((3,),), ((1, 2),)))
    # Marginally allowed nonzero gains can have an impossible cycle sum.
    tri = ((0, 1), (0, 2), (1, 2))
    check(decode(3, tri, 6, (1, 5, 1)) is None)  # cyclic 0->1->2->0 all +1
    check(decode(3, tri, 6, (1, 2, 1)) == (0, 1, 2))
    for gains in product(range(6), repeat=3):
        six = decode(3, tri, 6, gains)
        two = decode(3, tri, 2, tuple(x % 2 for x in gains))
        three = decode(3, tri, 3, tuple(x % 3 for x in gains))
        check((six is not None) == (two is not None and three is not None))
        if six is not None:
            check(six == tuple((3*a+4*b) % 6 for a, b in zip(two, three)))

    # Original eight-port graph, root fixed: tree enumeration, independent cycle formula.
    old = ((0, 1), (0, 4), (0, 5), (0, 7), (1, 2), (2, 3), (3, 4), (5, 6), (6, 7))
    tree, paths, rows = tree_code(8, old)
    count = 0; margins = [set() for _ in old]
    patterns = set()
    for values in product(range(1, 6), repeat=7):
        partial = [0]*9
        for i, value in zip(tree, values):
            partial[i] = value
        cols = tuple(sum(a*b for a, b in zip(row, partial)) % 6 for row in paths)
        gains = encode(8, old, 6, cols)
        if not all(gains):
            continue
        count += 1
        for ss, x in zip(margins, gains):
            ss.add(x)
        names = {}; pattern = []
        for c in cols:
            if c not in names:
                names[c] = len(names)
            pattern.append(names[c])
        patterns.add(tuple(pattern))
    check(count == ((5**5-5)*(5**4+5))//36 == 54600)
    check(len(patterns) == 570)
    check(all(ss == set(range(1, 6)) for ss in margins))
    check(not matrix_box(6, rows, ((1, 2, 3, 4, 5),)*9))
    # A permutation of named colours need not preserve their arithmetic difference.
    colour_swap = (0, 2, 1, 3, 4, 5)
    check((1-0) % 6 == (2-1) % 6)
    check((colour_swap[1]-colour_swap[0]) % 6 !=
          (colour_swap[2]-colour_swap[1]) % 6)

    # Complete reduced forms, with explicit represented-unit genus labels.
    forms = reduced_forms(-120)
    check(forms == ((1, 0, 30), (2, 0, 15), (3, 0, 10), (5, 0, 6)))
    labels = []
    for a, b, c in forms:
        observed = {genus120(a*x*x+b*x*y+c*y*y)
                    for x in range(-8, 9) for y in range(-8, 9)
                    if (x, y) != (0, 0) and gcd(a*x*x+b*x*y+c*y*y, 120) == 1}
        check(len(observed) == 1)
        labels.append(next(iter(observed)))
    check(len(set(labels)) == 4)
    for pair in product(range(2), repeat=2):
        check(rotate_genus(rotate_genus(rotate_genus(pair))) == pair)
        check((pair == (0, 0)) == (rotate_genus(pair) == pair))
    # No nontrivial homomorphism C3 -> C2^r; order-three automorphism is different.
    for rank in range(1, 6):
        for x in product(range(2), repeat=rank):
            check(tuple((3*t) % 2 for t in x) == x)
    check(reduced_forms(-44) == ((1, 0, 11), (3, -2, 4), (3, 2, 4)))
    check(reduced_forms(-76) == ((1, 0, 19), (4, -2, 5), (4, 2, 5)))

    # Exact marked K4-edge palette bridge; a separate odd duality doubles S4.
    vertex_actions = {vertex_palette_action(p) for p in permutations(range(4))}
    pair_actions = set()
    even_actions = set()
    for sigma in permutations(range(3)):
        for flips in product(range(2), repeat=3):
            action = tuple(2*sigma[d]+(s ^ flips[d]) for d in range(3) for s in range(2))
            pair_actions.add(action)
            if sum(flips) % 2 == 0:
                even_actions.add(action)
    check(len(vertex_actions) == 24 and vertex_actions == even_actions)
    check(len(pair_actions) == 48)
    duality = tuple(c ^ 1 for c in range(6))
    check(duality not in vertex_actions)
    doubled = vertex_actions | {tuple(c ^ 1 for c in action) for action in vertex_actions}
    check(doubled == pair_actions)
    for action in pair_actions:
        check(tuple(action[c ^ 1] for c in range(6)) == tuple(c ^ 1 for c in action))
    for edge in combinations(range(4), 2):
        opposite = tuple(x for x in range(4) if x not in edge)
        check(edge_colour(opposite) == (edge_colour(edge) ^ 1))
    for t in range(4):
        action = vertex_palette_action(tuple(x ^ t for x in range(4)))
        expected = tuple(2*(d-1)+(s ^ int(t not in (0, d))) for d in (1, 2, 3) for s in (0, 1))
        check(action == expected)

    for bad in (True, 1.0):
        rejects(lambda bad=bad: digit_word(bad))
        rejects(lambda bad=bad: decode(2, ((0, 1),), 6, (bad,)))
        rejects(lambda bad=bad: rotate_genus((bad, 0)))
        rejects(lambda bad=bad: matrix_box(6, ((bad,),), ((1,),)))
        rejects(lambda bad=bad: vertex_palette_action((0, bad, 2, 3)))
    rejects(lambda: decode(2, ((0, 1),), 6, (6,)))
    rejects(lambda: decode(3, ((0, 1),), 6, (1,)))
    rejects(lambda: matrix_box(6, ((1,),), ((),)))
    rejects(lambda: genus120(3))
    return {
        "checks": checks,
        "digit_words_first_digit_order": words,
        "digit_carry_depths": [1, 7],
        "tournament": {"source": "C9{1,3,4,7}", "automorphisms": len(automorphisms),
                       "fibres": [[r+3*t for t in range(3)] for r in range(3)],
                       "listed_word_orientation": "reverse of +3"},
        "small_connected_graphs": graph_count,
        "old_port_graph": {"n": 8, "edges": old, "tree_edges": tree, "cycle_rows": rows,
                           "root_fixed_colourings": count, "labelled_colourings": 6*count,
                           "equality_patterns": len(patterns), "independent_edge_choices": 5**9,
                           "tree_choices": 5**7},
        "idoneal_marked_control": {"n": 30, "discriminant": -120, "forms": forms,
                                   "character_labels_chi8_chiminus3": labels,
                                   "completeness_of_65_used": False},
        "marked_K4_palette": {"vertex_actions": len(vertex_actions),
                              "pair_preserving_actions": len(pair_actions),
                              "vertex_actions_have_even_flip_parity": True,
                              "opposite_edge_duality": duality,
                              "duality_induced_by_vertex_permutation": False},
        "boundaries": ["No plane equality gadget or universal colouring theorem",
                       "Named Z6 palette and exact cycle relations retained",
                       "Genus rotation is a chosen automorphism, not a field automorphism",
                       "No canonical numerical bijection of the 65-element lists"]}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))

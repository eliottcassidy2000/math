#!/usr/bin/env python3
"""Stdlib exact checks for the odd-order orientation-anchor theorem.
Universe: all labelled binary n*n matrices; ordered tournament pair on common
V=[n], with a distinguished additional root of the second tournament.
No filtering of matrices or tournaments. No repository writes.
"""
from collections import Counter
from itertools import permutations
import json


def require(condition, message="verification failed"):
    if not condition:
        raise AssertionError(message)


def pairs(n):
    return tuple((i, j) for i in range(n) for j in range(i + 1, n))


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def closure(n, gens):
    ident = tuple(range(n))
    known, todo = {ident}, [ident]
    while todo:
        p = todo.pop()
        for q in gens:
            r = compose(p, q)
            if r not in known:
                known.add(r)
                todo.append(r)
    return frozenset(known)


def subgroups(n):
    ps = tuple(permutations(range(n)))
    ident = frozenset([tuple(range(n))])
    seen, todo = {ident}, [ident]
    while todo:
        h = todo.pop()
        for p in ps:
            if p not in h:
                k = closure(n, tuple(h) + (p,))
                if k not in seen:
                    seen.add(k)
                    todo.append(k)
    return sorted(seen, key=lambda h: (len(h), sorted(h)))


def orbit_count(items, h, action):
    todo, out = set(items), 0
    while todo:
        x = min(todo)
        todo -= {action(p, x) for p in h}
        out += 1
    return out


def adjacency(n, mask):
    rows = [0] * n
    for k, (i, j) in enumerate(pairs(n)):
        u, v = (i, j) if (mask >> k) & 1 else (j, i)
        rows[u] |= 1 << v
    return tuple(rows)


def relabel_adj(rows, p):
    out = [0] * len(rows)
    for i, row in enumerate(rows):
        for j in range(len(rows)):
            out[p[i]] |= ((row >> j) & 1) << p[j]
    return tuple(out)


def ordinary(m, p):
    n = len(p)
    return sum(((m >> (i * n + j)) & 1) << (p[i] * n + p[j])
               for i in range(n) for j in range(n))


def encode(m, a):
    n = len(a)
    t, u = [0] * n, [0] * (n + 1)
    for i, j in pairs(n):
        s, r = (i, j) if (a[i] >> j) & 1 else (j, i)
        x, y = (m >> (s * n + r)) & 1, (m >> (r * n + s)) & 1
        t[s if x else r] |= 1 << (r if x else s)
        u[s if y else r] |= 1 << (r if y else s)
    for i in range(n):
        d = (m >> (i * n + i)) & 1
        u[n if d else i] |= 1 << (i if d else n)
    return tuple(t), tuple(u)


def decode(t, u, a):
    n, m = len(a), 0
    for i, j in pairs(n):
        s, r = (i, j) if (a[i] >> j) & 1 else (j, i)
        m |= ((t[s] >> r) & 1) << (s * n + r)
        m |= ((u[s] >> r) & 1) << (r * n + s)
    for i in range(n):
        m |= ((u[n] >> i) & 1) << (i * n + i)
    return m


def anchors_by_orbits(n, h):
    """Construct every invariant anchor; a reversed directed orbit forbids one."""
    remaining, blocks = set(pairs(n)), []
    while remaining:
        i, j = min(remaining)
        directed = {(p[i], p[j]) for p in h}
        if (j, i) in directed:
            return []
        undirected = {tuple(sorted(e)) for e in directed}
        require(len(directed) == len(undirected))
        remaining -= undirected
        blocks.append(directed)
    out = []
    for mask in range(1 << len(blocks)):
        rows = [0] * n
        for k, directed in enumerate(blocks):
            for i, j in directed:
                s, r = (i, j) if (mask >> k) & 1 else (j, i)
                rows[s] |= 1 << r
        out.append(tuple(rows))
    return out


def solve_anchor_constraints(n, generators):
    """Compact affine solution or reversal-word certificate; never lists H.

    Elementary graph work is O((s+1)*n*n); this function does not enumerate
    the 2**e solution family. Signed word letters are 1-based generator
    indices, negative for inverses, applied from left to right.
    """
    es = pairs(n)
    graph = {e: [] for e in es}
    for k, p in enumerate(generators, 1):
        for i, j in es:
            image = tuple(sorted((p[i], p[j])))
            flip = int(p[i] > p[j])
            graph[(i, j)].append((image, flip, k))
            graph[image].append(((i, j), flip, -k))
    values, components, parent = {}, [], {}

    def path(e):
        letters = []
        while parent[e] is not None:
            e, letter = parent[e]
            letters.append(letter)
        return list(reversed(letters))

    for root in es:
        if root in values:
            continue
        values[root], component, todo, parent[root] = 0, [], [root], None
        while todo:
            e = todo.pop()
            component.append(e)
            for f, flip, letter in graph[e]:
                expected = values[e] ^ flip
                if f in values:
                    if values[f] != expected:
                        # Root->e->f followed by inverse(root->f).
                        word = path(e) + [letter] + [-k for k in reversed(path(f))]
                        return {"consistent": False, "edge": root, "word": word}
                else:
                    values[f] = expected
                    parent[f] = (e, letter)
                    todo.append(f)
        components.append(tuple(component))
    mask = sum(values[e] << k for k, e in enumerate(es))
    return {"consistent": True, "values": values,
            "particular_anchor": adjacency(n, mask), "components": tuple(components)}


def inverse_permutation(p):
    inv = [0] * len(p)
    for i, j in enumerate(p):
        inv[j] = i
    return tuple(inv)


def verify_conflict_word(n, generators, certificate):
    """Independent negative verifier: evaluates only the supplied word."""
    require(not certificate["consistent"])
    g = tuple(range(n))
    for letter in certificate["word"]:
        require(1 <= abs(letter) <= len(generators))
        p = generators[abs(letter)-1]
        g = compose(p if letter > 0 else inverse_permutation(p), g)
    i, j = certificate["edge"]
    require(i != j and g[i] == j and g[j] == i)
    require(len(certificate["word"]) <= 2 * len(pairs(n)) - 1)
    return g


def anchors_by_generators(n, generators):
    """Exponential output wrapper around the compact parity solver."""
    result = solve_anchor_constraints(n, generators)
    if not result["consistent"]:
        return []
    es, components, out = pairs(n), result["components"], []
    for mask in range(1 << len(components)):
        bits = dict(result["values"])
        for k, component in enumerate(components):
            for e in component:
                bits[e] ^= (mask >> k) & 1
        out.append(adjacency(n, sum(bits[e] << k for k, e in enumerate(es))))
    return out


def fixed_dimensions(n, h):
    vertices = orbit_count(range(n), h, lambda p, i: p[i])
    edges = orbit_count(pairs(n), h, lambda p, e: tuple(sorted((p[e[0]], p[e[1]]))))
    ordered = orbit_count(((i, j) for i in range(n) for j in range(n)), h,
                          lambda p, e: (p[e[0]], p[e[1]]))
    return vertices, edges, ordered


def test_codec(n, h, anchors, exhaustive):
    probes = range(1 << (n * n)) if exhaustive else [0] + [1 << k for k in range(n*n)]
    checks = 0
    for a in anchors:
        for m in probes:
            t, u = encode(m, a)
            require(decode(t, u, a) == m)
            for p in h:
                t2, u2 = encode(ordinary(m, p), a)
                require(t2 == relabel_adj(t, p))
                require(u2 == relabel_adj(u, p + (n,)))
                checks += 1
    return checks


def test_gauge_and_cocycle(n):
    """Independent physical adjacency check, all anchors/perms, affine basis."""
    es = pairs(n)
    probes = [0] + [1 << k for k in range(n * n)]
    count = 0
    for a_mask in range(1 << len(es)):
        a = adjacency(n, a_mask)
        for flip in range(len(es)):
            b = adjacency(n, a_mask ^ (1 << flip))
            i, j = es[flip]
            s, r = (i, j) if (a[i] >> j) & 1 else (j, i)
            for m in probes:
                ta, ua = encode(m, a)
                tb, ub = encode(m, b)
                require(((tb[s] >> r) & 1) == (1 ^ ((ua[s] >> r) & 1)))
                require(((ub[s] >> r) & 1) == (1 ^ ((ta[s] >> r) & 1)))
                for k, l in es:
                    if (k, l) != (i, j):
                        require(((ta[k] >> l) & 1) == ((tb[k] >> l) & 1))
                        require(((ua[k] >> l) & 1) == ((ub[k] >> l) & 1))
                count += 1
        # c(g,e) records reversal of the source anchor arc at its image edge.
        def c(p, e):
            i, j = e
            s, r = (i, j) if (a[i] >> j) & 1 else (j, i)
            return 1 ^ ((a[p[s]] >> p[r]) & 1)
        ps = tuple(permutations(range(n)))
        for g in ps:
            for h in ps:
                gh = compose(g, h)
                for e in es:
                    he = tuple(sorted((h[e[0]], h[e[1]])))
                    require(c(gh, e) == (c(h, e) ^ c(g, he)))
    return count


def aut_order_transitive(rows):
    """Full exact aut order, using proven transitivity to fix vertex zero."""
    n, stabilizer = len(rows), 0
    # Every regular C3^2 translation has already been checked as an automorphism.
    for tail in permutations(range(1, n)):
        p = (0,) + tail
        good = True
        for i, j in pairs(n):
            if ((rows[i] >> j) & 1) != ((rows[p[i]] >> p[j]) & 1):
                good = False
                break
        stabilizer += good
    return n * stabilizer


def run():
    report = {"schema": 1, "status": "PASS", "source_pin": "03408f84ac43fcdcfb0887e5d6bc1fc8792415cf",
              "universe": "all binary matrices; all natural ordered aligned consecutive tournament pairs",
              "small_subgroups": [], "gauge_checks": 0, "hostiles": {},
              "check_style": "explicit require; all tests survive python -O"}
    for n in range(5):
        hs = subgroups(n)
        codec_checks = 0
        odd_count = 0
        for h in hs:
            anchors = anchors_by_orbits(n, h)
            brute = [a for mask in range(1 << len(pairs(n)))
                     if all(relabel_adj(a := adjacency(n, mask), p) == a for p in h)]
            require(set(anchors) == set(brute))
            require(set(anchors) == set(anchors_by_generators(n, h)))
            v, e, o = fixed_dimensions(n, h)
            compact = solve_anchor_constraints(n, tuple(h))
            if compact["consistent"]:
                require(len(compact["components"]) == e)
                require(all(relabel_adj(compact["particular_anchor"], p) ==
                            compact["particular_anchor"] for p in h))
            else:
                verify_conflict_word(n, tuple(h), compact)
            require(len(anchors) == (1 << e if len(h) % 2 else 0))
            if len(h) % 2:
                odd_count += 1
                require(o == 2 * e + v)
                codec_checks += test_codec(n, h, anchors, exhaustive=(n <= 3))
            else:
                require(not anchors)
        report["small_subgroups"].append({"n": n, "all_subgroups": len(hs),
                                           "odd_subgroups": odd_count,
                                           "codec_equivariance_checks": codec_checks})
        if n <= 3:
            report["gauge_checks"] += test_gauge_and_cocycle(n)
    c2 = closure(2, [(1, 0)])
    c3 = closure(3, [(1, 2, 0)])
    c3x3 = frozenset(tuple(3*((i//3+a)%3)+(i%3+b)%3 for i in range(9))
                     for a in range(3) for b in range(3))
    a4_generators = [(1, 2, 0, 3), (1, 3, 2, 0)]
    a4 = closure(4, a4_generators)
    require(len(a4) == 12)
    require(not anchors_by_generators(4, a4_generators))
    a4_certificate = solve_anchor_constraints(4, a4_generators)
    a4_reverser = verify_conflict_word(4, a4_generators, a4_certificate)
    c3x3_gens = [tuple(3*((i//3+1)%3)+i%3 for i in range(9)),
                 tuple(3*(i//3)+(i%3+1)%3 for i in range(9))]
    require(closure(9, c3x3_gens) == c3x3)
    require(set(anchors_by_generators(9, c3x3_gens)) == set(anchors_by_orbits(9, c3x3)))
    for name, n, h in [("C2_on_two", 2, c2), ("C3_on_three", 3, c3),
                       ("C3xC3_regular_on_nine", 9, c3x3), ("A4_odd_order_generators", 4, a4)]:
        anchors = anchors_by_orbits(n, h)
        v, e, o = fixed_dimensions(n, h)
        row = {"n": n, "group_order": len(h), "vertex_orbits": v,
               "unordered_pair_orbits": e, "ordered_pair_orbits_including_diagonal": o,
               "invariant_anchor_count": len(anchors), "fixed_matrix_count": 1 << o,
               "fixed_natural_pair_count": len(anchors)**2 * (1 << v)}
        if name == "A4_odd_order_generators":
            row["generator_orders"] = [3, 3]
            row["conflict_certificate"] = {
                "generators": a4_generators,
                "word_convention": "signed 1-based indices, left-to-right; negative means inverse",
                "edge": a4_certificate["edge"], "word": a4_certificate["word"],
                "evaluated_reversing_permutation": a4_reverser}
            row["generator_constraint_solution_count"] = len(anchors_by_generators(n, a4_generators))
        if anchors:
            row["affine_basis_codec_checks"] = test_codec(n, h, anchors, exhaustive=False)
        if name == "C3xC3_regular_on_nine":
            # Direct full automorphism count for all sixteen labelled anchors.
            row["full_anchor_automorphism_order_histogram"] = dict(Counter(
                aut_order_transitive(a) for a in anchors))
        report["hostiles"][name] = row
    # Nonfaithful warning: abstract C2 acting trivially has odd image order 1.
    report["nonfaithful_control"] = {
        "abstract_group_order": 2, "image_order": 1, "n": 2,
        "invariant_anchor_count": len(anchors_by_orbits(2, closure(2, []))),
        "codec_exists": True}
    # Pointwise fixed names leave S_(n-k), whose order is odd exactly n-k <= 1.
    report["point_individualization"] = [
        {"n": n, "minimum_distinct_ordered_names": max(0, n-1)} for n in range(7)]
    return report


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))

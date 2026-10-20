"""Exact target certificates, physical Kempe coordinates, and golden controls.

Run with python or python -O. Dependencies: networkx, sympy.
No optimizer is needed: the fourteen-probe cover is a fixed verified receipt.
"""
from collections import Counter
from itertools import combinations, permutations, product
from math import lcm
import networkx as nx
import sympy as sp


def need(value, message):
    if not value:
        raise RuntimeError(message)


def canon(word):
    seen = {}
    return tuple(seen.setdefault(c, len(seen)) for c in word)


def matchings(items):
    if not items:
        yield ()
        return
    a = items[0]
    for j in range(1, len(items)):
        for rest in matchings(items[1:j] + items[j + 1:]):
            yield ((a, items[j]),) + rest


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def rank_mod(rows, p):
    a = [list(row) for row in rows]
    r = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][col] % p), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][col] % p, -1, p)
        a[r] = [x * inv % p for x in a[r]]
        for i in range(r + 1, len(a)):
            c = a[i][col] % p
            if c:
                a[i] = [(x - c * y) % p for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def junction(word, label):
    triple, pairs = label
    return int(len({word[v] for v in triple}) == 3
               and all(word[a] == word[b] for a, b in pairs))


def pole_profile(vertices, edges, ports):
    vertices, edges, ports = tuple(vertices), tuple(edges), tuple(ports)
    variables = edges + tuple((v,) for v in ports)
    degree = Counter(v for e in variables for v in e)
    need(set(degree) == set(vertices) and all(degree[v] == 3 for v in vertices),
         "cubic interior")
    used = {v: 0 for v in vertices}
    values = [-1] * len(variables)
    result = Counter()

    def visit(left):
        if not left:
            result[tuple(values[len(edges):])] += 1
            return
        best, choices = None, 7
        for j in left:
            mask = 7
            for v in variables[j]:
                mask &= 7 ^ used[v]
            if not mask:
                return
            if best is None or mask.bit_count() < choices.bit_count():
                best, choices = j, mask
        for c in range(3):
            bit = 1 << c
            if not choices & bit:
                continue
            values[best] = c
            for v in variables[best]:
                used[v] |= bit
            visit(tuple(j for j in left if j != best))
            for v in variables[best]:
                used[v] ^= bit

    visit(tuple(range(len(variables))))
    return result


def monomer_defects(vertices, edges, ports):
    """No colouring enumeration: internal near-perfect matchings and path parity."""
    vertices = set(vertices)
    adj = {v: set() for v in vertices}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)

    def pms(left):
        if not left:
            yield set()
            return
        v = min(left, key=lambda x: (len(adj[x] & left), x))
        for u in sorted(adj[v] & left):
            for rest in pms(left - {u, v}):
                yield rest | {frozenset((u, v))}

    answer = []
    for i, monomer in enumerate(ports):
        count = 0
        for matching in pms(vertices - {monomer}):
            residual = {("v", v): [] for v in vertices}
            for u, v in edges:
                if frozenset((u, v)) not in matching:
                    residual["v", u].append(("v", v))
                    residual["v", v].append(("v", u))
            for j, v in enumerate(ports):
                if j == i:
                    continue
                residual["p", j] = [("v", v)]
                residual["v", v].append(("p", j))
            need(all(len(residual["v", v]) == 2 for v in vertices),
                 "two-regular coloured complement")
            unseen = set(residual)
            cycles, paths, valid = 0, 0, True
            while unseen:
                start = next(iter(unseen))
                component, todo = {start}, [start]
                unseen.remove(start)
                while todo:
                    for v in residual[todo.pop()]:
                        if v in unseen:
                            unseen.remove(v)
                            component.add(v)
                            todo.append(v)
                terminals = sum(v[0] == "p" for v in component)
                if terminals == 0:
                    cycles += 1
                    valid &= len(component) % 2 == 0
                else:
                    need(terminals == 2, "not a boundary path")
                    paths += 1
                    valid &= (len(component) - 2) % 2 == 1
            need(paths == 3, "three boundary paths")
            if valid:
                count += 2 ** cycles
        answer.append(24 * count)
    return tuple(answer)


P = ((0, 1), (0, 4), (0, 5), (1, 2), (1, 6), (2, 3), (2, 7),
     (3, 4), (3, 8), (4, 9), (5, 7), (5, 8), (6, 8), (6, 9), (7, 9))


def connected_j():
    edges = {frozenset(e) for e in ((0, -1), (0, -2), (0, -3), (-4, -5), (-6, -7))}
    for fresh, first, second in ((1, (0, -1), (-4, -5)),
                                  (9, (0, -2), (-6, -7))):
        ren = {v: fresh + j for j, v in enumerate(v for v in range(10) if v not in (0, 5))}
        edges.remove(frozenset(first))
        edges.remove(frozenset(second))
        edges.update(frozenset((ren[u], ren[v])) for u, v in P if 0 not in (u, v) and 5 not in (u, v))
        edges.update(frozenset((v, ren[p])) for v, p in zip(first + second, (1, 4, 7, 8)))
    ports = tuple(next(next(v for v in e if v != -i - 1) for e in edges if -i - 1 in e)
                  for i in range(7))
    internal = tuple(sorted(tuple(sorted(e)) for e in edges if min(e) >= 0))
    return tuple(range(17)), internal, ports


def phi_powers(n):
    values = [(1, 0)]
    for _ in range(n):
        a, b = values[-1]
        values.append((b, a + b))
    return values


def phi_moment(coefficients):
    powers = phi_powers(len(coefficients))
    return tuple(sum(c * powers[i][j] for i, c in enumerate(coefficients)) for j in (0, 1))


def loops(a, b):
    parent = list(range(2 * len(a)))

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x

    for u, v in a + b:
        parent[find(u)] = find(v)
    return len({find(i) for i in range(len(parent))})


def noncrossing(matching):
    return not any(a < c < b < d or c < a < d < b
                   for (a, b), (c, d) in combinations(matching, 2))


def main():
    all_words = tuple(product(range(3), repeat=7))
    physical = tuple(w for w in all_words if all(w.count(c) % 2 for c in range(3)))
    words = sorted({canon(w) for w in physical})
    need((len(all_words), len(physical), len(words)) == (2187, 546, 91), "word bank")
    labels = [(t, m) for t in combinations(range(7), 3)
              for m in matchings(tuple(v for v in range(7) if v not in t))]
    A = [tuple(junction(w, label) for w in words) for label in labels]
    W = [tuple((1 if sorted(Counter(w).values()) == [1, 1, 5] else -1)
               if w.count(w[i]) == 1 else 0 for w in words) for i in range(7)]
    need(all(dot(a, w) == 0 for a in A for w in W), "seven blind directions")
    need(rank_mod(A, 101) == 84 and rank_mod(W, 101) == 7, "exact-rank bounds")
    target = ((0, 1, 2), ((3, 4), (5, 6)))
    need(labels[0] == target, "target index")
    cover = (27, 28, 29, 42, 57, 58, 59, 72, 75, 90, 91, 92, 93, 96)
    zeros = [j for j, row in enumerate(A) if dot(row, A[0]) == 0]
    need(len(zeros) == 42 and set(cover) <= set(zeros), "zero-overlap bank")
    for word in all_words:
        if not all(word.count(c) % 2 for c in range(3)):
            continue
        rejected = any(junction(word, labels[j]) for j in cover)
        need(rejected == (not junction(word, target)), "all-word fourteen-probe coverage")
    for new_target in labels:
        order = new_target[0] + new_target[1][0] + new_target[1][1]
        for word in physical:
            moved = tuple(word[order[i]] for i in range(7))
            need(any(junction(moved, labels[j]) for j in cover)
                 == (not junction(word, new_target)), "all 105 relabelled targets")
    print("J certificate: 14 explicit zero probes cover all 492 forbidden parity words / 82 orbits.")
    print("All 105 target relabellings pass; normalization is self-closure 54m.")
    print("Fourteen labels:", tuple(labels[j] for j in cover))
    fake = [0] * len(words)
    fake[next(i for i, v in enumerate(A[0]) if v)] = 9
    need(all(dot(fake, A[j]) == 0 for j in cover) and dot(fake, A[0]) == 9,
         "nonphysical positive fifteen-probe hostile")
    need(tuple(fake) != A[0] and any(dot(fake, row) != dot(A[0], row) for row in A),
         "physical predicate is necessary; all 105 detect positive hostile")
    print("Physical hypothesis is essential: a weight-nine supported orbit passes the 15 probes but is not J.")
    print("All seven generic signed blind vectors remain invisible to all 105 probes.")

    # Enumerate the 630 physically distinct labelled P5 tree closures.
    arrangements = set()
    for p in permutations(range(7)):
        a = tuple(sorted(p[:2])) + p[2:5] + tuple(sorted(p[5:]))
        b = tuple(sorted(p[5:])) + p[2:5][::-1] + tuple(sorted(p[:2]))
        arrangements.add(min(a, b))
    need(len(arrangements) == 630, "path tree closure bank")
    Ksum = [0] * len(words)
    for p in arrangements:
        for j, word in enumerate(words):
            flow, good = 0, True
            for t, i in enumerate(p, 1):
                flow ^= word[i] + 1
                if t in (2, 3, 4, 5):
                    good &= flow != 0
            Ksum[j] += int(good)
    Jsum = [sum(row[j] for row in A) for j in range(len(words))]
    Wsum = [sum(row[j] for row in W) for j in range(len(words))]
    need(all(11 * k - 112 * j == -180 * w for k, j, w in zip(Ksum, Jsum, Wsum)),
         "positive aggregate identity")
    need(set(zip(Jsum, Ksum)) == {(15, 120), (9, 108)}, "aggregate orbit values")
    need(all(dot(row, [11 * x for x in Ksum]) == dot(row, [112 * x for x in Jsum])
             for row in A), "physical-ensemble hostile")
    print("Aggregate: 11*C_path630 - 112*C_junction105 = 180*sum_i N_i.")
    print("Hostile physical ensembles 11*all630 path trees and 112*all105 junctions have identical 105 readouts.")

    def audit_pole(name, vertices, edges, ports):
        need(nx.is_connected(nx.Graph(edges)), "connected physical control")
        prof = pole_profile(vertices, edges, ports)
        need(all(all(w.count(c) % 2 for c in range(3)) for w in prof), "physical parity")
        f = tuple(prof[w] for w in words)
        N = tuple(-6 * dot(f, w) for w in W)
        need(N == monomer_defects(vertices, edges, ports), "independent monomer-dimer identity")
        need(all(v >= 0 and v % 24 == 0 for v in N), "Kempe divisibility")
        ck, cj = 6 * dot(f, Ksum), 6 * dot(f, Jsum)
        need(11 * ck - 112 * cj == 180 * sum(N), "physical aggregate")
        if name:
            print(name, "boundary mass=", sum(prof.values()), "N=", N)
        return prof, N

    tree_ports = (0, 0, 2, 1, 3, 4, 4)  # singleton middle port now labelled 2
    tree_profile, tree_N = audit_pole("P5 middle-at-2", range(5),
                                      tuple((i, i + 1) for i in range(4)), tree_ports)
    removed = {0, 1, 3}
    pet_edges = tuple(e for e in P if not removed & set(e))
    pet_ports = (5, 6, 4, 2, 2, 4, 8)  # two nonzero defects now at ports 0,1
    _, pet_N = audit_pole("Petersen minus 0,1,3", sorted(set(range(10)) - removed),
                         pet_edges, pet_ports)
    need(tree_N == (0, 0, 24, 0, 0, 0, 0), "tree defect")
    need(pet_N == (24, 24, 0, 0, 0, 0, 0), "Petersen defect")
    need(phi_moment(tree_N) == phi_moment(pet_N) == (24, 24), "physical golden collision")
    print("Physical golden collision: 24*phi^2 = 24*(1+phi), distinct N vectors of two connected single poles.")
    q = connected_j()
    qprofile, qN = audit_pole("Connected 17-vertex Q", *q)
    need(qprofile == Counter({w: 4 for w in physical if junction(w, target)}), "Q=4J all words")
    need(all(sum(qprofile[w] * junction(w, labels[j]) for w in physical) == 0 for j in cover),
         "Q zero-probe certificate")
    need(sum(qprofile[w] * junction(w, target) for w in physical) == 216, "Q normalization")
    atlas = 0
    for g in nx.graph_atlas_g():
        n = len(g)
        if n not in (5, 7) or 3 * n - 2 * g.number_of_edges() != 7:
            continue
        if not nx.is_connected(g) or max(dict(g.degree()).values()) > 3:
            continue
        ports = tuple(v for v in g for _ in range(3 - g.degree(v)))
        audit_pole(None, tuple(g), tuple(g.edges()), ports)
        atlas += 1
    print("Physical monomer-dimer and aggregate controls: all", atlas,
          "connected simple 5/7-interior-vertex seven-poles from NetworkX atlas.")

    # Golden loop parameter: exact field relation and a crossing detector.
    bank = tuple(matchings(tuple(range(8))))
    planar = tuple(m for m in bank if noncrossing(m))
    need((len(bank), len(planar)) == (105, 14), "eight-port carriers")
    powers = phi_powers(4)
    pair_gram = [[powers[loops(a, b)] for b in planar] for a in planar]
    A0 = sp.Matrix([[v[0] for v in row] for row in pair_gram])
    A1 = sp.Matrix([[v[1] for v in row] for row in pair_gram])
    block = A0.row_join(A1).col_join(A1.row_join(A0 + A1))
    null = block.nullspace()
    need(len(null) == 2, "one-dimensional golden field radical")
    scale = lcm(*(int(x.q) for x in null[0]))
    vec = [(int(null[0][i] * scale), int(null[0][i + 14] * scale)) for i in range(14)]
    need(rank_mod([[pow(4, loops(a, b), 11) for b in planar] for a in planar], 11) == 13,
         "golden planar lower rank")

    def pair_eval(test, value):
        out = [0, 0]
        for m, (x, y) in zip(planar, vec):
            a, b = powers[loops(test, m)] if value == "phi" else (3 ** loops(test, m), 0)
            out[0] += a * x + b * y
            out[1] += a * y + b * x + b * y
        return tuple(out)

    need(all(pair_eval(m, "phi") == (0, 0) for m in planar), "exact golden null vector")
    crossing = next(m for m in bank if pair_eval(m, "phi") != (0, 0))
    ternary = next(m for m in planar if pair_eval(m, 3) != (0, 0))
    need(rank_mod([[pow(4, loops(a, b), 11) for b in bank] for a in bank], 11) == 105,
         "full golden rank")
    need(rank_mod([[pow(3, loops(a, b), 101) for b in planar] for a in planar], 101) == 14,
         "ternary planar rank")
    print("Golden planar 14-diagram Gram rank=13; full 105-diagram golden rank=105; ternary planar rank=14.")
    print("Golden null-vector coefficients (a,b) for a+b*phi:", tuple(vec))
    print("Crossing detector:", crossing, "value=", pair_eval(crossing, "phi"))
    print("Ternary detector:", ternary, "value=", pair_eval(ternary, 3))
    # Exact stack-height-at-most-three noncrossing matching counts.
    state = [1, 0, 0, 0]
    heights = []
    for steps in range(1, 15):
        state = [sum(state[j] for j in (i - 1, i + 1) if 0 <= j < 4) for i in range(4)]
        if steps % 2 == 0:
            heights.append(state[0])
    need(heights == [1, 2, 5, 13, 34, 89, 233], "bounded-stack Fibonacci recurrence")
    print("Noncrossing matching nesting height<=3, 2..14 ports:", heights)
    print("PASS exact covers, physical controls, positive aggregates, and golden quotient hostiles.")


if __name__ == "__main__":
    main()

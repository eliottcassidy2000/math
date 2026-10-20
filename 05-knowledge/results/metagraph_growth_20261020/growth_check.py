"""Exact controls for canonical tournament inflation and fixed-path growth.

Run from repository root:
  python 05-knowledge/results/metagraph_growth_20261020/growth_check.py
All checks remain active under python -O; no third-party dependencies.
"""
from collections import Counter
from functools import reduce
from itertools import combinations, permutations, product
from math import factorial, gcd


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def pairs(n):
    return list(combinations(range(n), 2))


def free_pairs(n):
    return [(i, j) for i, j in pairs(n) if j > i + 1]


def tournament(n, code, fixed_path=False):
    a = [0] * n
    slots = free_pairs(n) if fixed_path else pairs(n)
    if fixed_path:
        for i in range(n - 1):
            a[i] |= 1 << (i + 1)
    for r, (i, j) in enumerate(slots):
        u, v = (i, j) if code >> r & 1 else (j, i)
        a[u] |= 1 << v
    return tuple(a)


def arc(a, i, j):
    return (a[i] >> j) & 1


def induced(a, vertices):
    return tuple(sum(arc(a, u, v) << j for j, v in enumerate(vertices))
                 for u in vertices)


def packed(a, order):
    return sum(arc(a, order[i], order[j]) << r
               for r, (i, j) in enumerate(pairs(len(order))))


def canonical(a):
    # Isomorphisms preserve outdegree; sorting the degree blocks loses nothing.
    groups = {}
    for v, row in enumerate(a):
        groups.setdefault(row.bit_count(), []).append(v)
    return min(packed(a, tuple(v for block in blocks for v in block))
               for blocks in product(*(permutations(groups[d])
                                       for d in sorted(groups))))


def clone_chains(a):
    n = len(a)
    successor, predecessor = {}, {}
    for u, v in pairs(n):
        if all(arc(a, u, z) == arc(a, v, z)
               for z in range(n) if z not in (u, v)):
            if not arc(a, u, v):
                u, v = v, u
            need(u not in successor and v not in predecessor, "clone branching")
            successor[u], predecessor[v] = v, u
    chains = []
    for start in range(n):
        if start in predecessor:
            continue
        chain = [start]
        while chain[-1] in successor:
            chain.append(successor[chain[-1]])
        chains.append(tuple(chain))
    need(sum(map(len, chains)) == n, "clone cycle")
    for c in chains:
        need(all(arc(a, u, v) for u, v in combinations(c, 2)), "nontransitive")
        need(all(len({arc(a, u, z) for u in c}) == 1
                 for z in range(n) if z not in c), "nonmodule")
    return tuple(chains)


def maximal_transitive_modules(a):
    # Independent all-subset test, separate from the homogeneous-pair algorithm.
    n = len(a)
    modules = []
    for mask in range(1, 1 << n):
        vertices = tuple(v for v in range(n) if mask >> v & 1)
        if any(len({arc(a, u, z) for u in vertices}) != 1
               for z in range(n) if not mask >> z & 1):
            continue
        degrees = [sum(arc(a, u, v) for v in vertices) for u in vertices]
        if sorted(degrees) == list(range(len(vertices))):
            modules.append(mask)
    return {m for m in modules if not any(m != q and m & q == m for q in modules)}


def inflate(a, k):
    n = len(a)
    return tuple(sum((int(u // k == v // k and u < v)
                      if u // k == v // k else arc(a, u // k, v // k)) << v
                     for v in range(k * n))
                 for u in range(k * n))


def uniform_root(a, k):
    chains = clone_chains(a)
    need(all(len(c) % k == 0 for c in chains), "nondivisible root")
    order = tuple(v for c in chains for v in c)
    root_vertices = tuple(c[j] for c in chains for j in range(0, len(c), k))
    root = induced(a, root_vertices)
    need(induced(a, order) == inflate(root, k), "root reconstruction")
    return root


def transpose(a):
    n = len(a)
    return tuple(sum(arc(a, n - 1 - v, n - 1 - u) << v
                     for v in range(n)) for u in range(n))


def tau(a):
    n = len(a)
    b = list(a)
    for u, v in free_pairs(n):
        b[u] ^= 1 << v
        b[v] ^= 1 << u
    return tuple(b)


def diff_edges(a, b):
    return {(u, v) for u, v in pairs(len(a)) if arc(a, u, v) != arc(b, u, v)}


def defect_edges(n, k):
    return {(u, v) for u, v in free_pairs(k * n)
            if u // k == v // k or v // k == u // k + 1}


def inflate_checks(a, k):
    b = inflate(a, k)
    base = sorted(map(len, clone_chains(a)))
    need(sorted(map(len, clone_chains(b))) == [k * s for s in base], "length law")
    need(canonical(uniform_root(b, k)) == canonical(a), "root class")
    need(transpose(b) == inflate(transpose(a), k), "transpose inflation")
    actual = diff_edges(tau(b), inflate(tau(a), k))
    need(actual == defect_edges(len(a), k), "tau defect support")
    need({(len(b) - 1 - v, len(b) - 1 - u) for u, v in actual} == actual,
         "transpose-invariant defect")
    expected = len(a) * (k - 1) * (k - 2) // 2 + (len(a) - 1) * (k * k - 1)
    need(len(actual) == expected, "tau defect cardinality")


def hp_subset_dp(a):
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            if not dp[mask][v]:
                continue
            for u in range(n):
                if not mask >> u & 1 and arc(a, v, u):
                    dp[mask | 1 << u][u] += dp[mask][v]
    return sum(dp[-1])


def block_word_count(a, k):
    n = len(a)
    counts = [0] * n
    total = 0
    numerator = factorial(k) ** n

    def visit(word):
        nonlocal total
        if len(word) == n * k:
            denominator = 1
            run = 1
            for u, v in zip(word, word[1:]):
                if u == v:
                    run += 1
                else:
                    denominator *= factorial(run)
                    run = 1
            denominator *= factorial(run)
            need(numerator % denominator == 0, "run multinomial")
            total += numerator // denominator
            return
        for v in range(n):
            if counts[v] == k:
                continue
            if word and word[-1] != v and not arc(a, word[-1], v):
                continue
            counts[v] += 1
            visit(word + (v,))
            counts[v] -= 1

    visit(())
    return total


def main():
    print("FINITE-EXACT: canonical transitive inflation and path growth.")
    module_cases = 0
    raw_classes = {}
    for n in range(1, 6):
        raw_classes[n] = set()
        for code in range(1 << (n * (n - 1) // 2)):
            a = tournament(n, code)
            raw_classes[n].add(canonical(a))
            chains = clone_chains(a)
            masks = {sum(1 << v for v in c) for c in chains}
            need(masks == maximal_transitive_modules(a), "independent module check")
            quotient = induced(a, tuple(c[0] for c in chains))
            need(all(len(c) == 1 for c in clone_chains(quotient)), "quotient twins")
            module_cases += 1
    print("All labelled tournaments through n=5: module controls =", module_cases)

    classes = {}
    inflation_cases = 0
    for n in range(1, 7):
        classes[n] = {}
        for code in range(1 << len(free_pairs(n))):
            a = tournament(n, code, fixed_path=True)
            key = canonical(a)
            classes[n].setdefault(key, a)
            if n <= 5:
                for k in (2, 3):
                    inflate_checks(a, k)
                    inflation_cases += 1
        hist = Counter(reduce(gcd, map(len, clone_chains(a)))
                       for a in classes[n].values())
        if n <= 5:
            need(set(classes[n]) == raw_classes[n], "all-labelled versus path census")
        print("n =", n, "classes =", len(classes[n]), "chain-gcd histogram =", dict(sorted(hist.items())))
        for k in range(1, n + 1):
            if n % k:
                continue
            roots = [canonical(uniform_root(a, k)) for a in classes[n].values()
                     if all(len(c) % k == 0 for c in clone_chains(a))]
            need(len(roots) == len(set(roots)) == len(classes[n // k]), "inflation image census")
            need(set(roots) == set(classes[n // k]), "inflation image identities")
    need([len(classes[n]) for n in range(1, 7)] == [1, 1, 2, 4, 12, 56], "known small census")
    print("Path presentations with k=2,3: exact inflation checks =", inflation_cases)

    # Independent class-equivalence control: no degree filtering.
    independent = 0
    for n in range(1, 7):
        seen = {}
        for key, a in classes[n].items():
            full = min(packed(a, order) for order in permutations(range(n)))
            need(full not in seen, "independent class collision")
            seen[full] = key
            independent += 1
    print("Independent full-permutation class representatives checked =", independent)

    # The black-line defect already appears at the first nontrivial doubling.
    transitive2 = tournament(2, 0, fixed_path=True)
    doubled = inflate(transitive2, 2)
    flipped = tau(doubled)
    need(sorted(map(len, clone_chains(doubled))) == [4], "positive doubling control")
    need(sorted(map(len, clone_chains(flipped))) == [1, 1, 2], "hostile doubling control")
    print("Hostile tau(D2(T2)): chain lengths (4) -> (1,1,2); tau exits doubling image.")

    c3 = tournament(3, 0, fixed_path=True)
    c3_double = inflate(c3, 2)
    fibre_data = []
    for a in (c3, c3_double):
        n = len(a)
        orders = list(permutations(range(n)))
        hp_count = sum(all(arc(a, u, v) for u, v in zip(order, order[1:]))
                       for order in orders)
        automorphisms = sum(induced(a, order) == a for order in orders)
        fibre_data.append((hp_count, automorphisms, hp_count // automorphisms))
    need(fibre_data == [(3, 3, 1), (45, 3, 15)], "non-block HP fibre")
    print("C3 and D2(C3) (Hamilton paths, automorphisms, fibre):", fibre_data)
    word_twins = ((0, 2, 3, 4, 5, 1), (1, 2, 3, 4, 5, 0))
    need(tuple(v // 2 for v in word_twins[0]) == tuple(v // 2 for v in word_twins[1]),
         "same block word")
    need(all(arc(c3_double, u, v) for hp in word_twins for u, v in zip(hp, hp[1:])),
         "block-word hostile HPs")
    partner_scores = [sorted(row.bit_count() for row in tau(induced(c3_double, hp)))
                      for hp in word_twins]
    need(partner_scores == [[2, 2, 2, 3, 3, 3], [1, 2, 2, 3, 3, 4]],
         "block word loses clone allocation")
    print("Same block word ABBCCA, tau-partner scores:", partner_scores)
    walk_cases = 0
    for n in range(1, 5):
        for code in range(1 << len(free_pairs(n))):
            a = tournament(n, code, fixed_path=True)
            for k in ((2, 3) if n <= 3 else (2,)):
                need(block_word_count(a, k) == hp_subset_dp(inflate(a, k)),
                     "independent weighted-walk count")
                walk_cases += 1
    print("Block-run multinomial versus independent HP subset-DP controls =", walk_cases)

    # The fixed example has two HP endpoints with nonisomorphic deletion children.
    arcs = [(0, 1), (1, 2), (2, 3), (0, 2), (1, 3), (3, 0)]
    a = tuple(sum(1 << v for u, v in arcs if u == x) for x in range(4))
    hp1, hp2 = (0, 1, 2, 3), (2, 3, 0, 1)
    need(all(arc(a, u, v) for hp in (hp1, hp2) for u, v in zip(hp, hp[1:])), "HP witness")
    children = [induced(a, hp[:-1]) for hp in (hp1, hp2)]
    need(sorted(row.bit_count() for row in children[0]) == [0, 1, 2], "transitive child")
    need(sorted(row.bit_count() for row in children[1]) == [1, 1, 1], "cyclic child")
    print("Hostile endpoint deletion: one order-four class -> transitive and cyclic triples.")

    # Check per-parent Mode-B colour law in the nondegenerate line range n=3..5.
    growth_cases = 0
    for n in range(3, 6):
        counts = {}
        for code in range(1 << len(free_pairs(n + 2))):
            a = tournament(n + 2, code, fixed_path=True)
            child = induced(a, tuple(range(1, n + 1)))
            key = packed(child, tuple(range(n)))
            counts.setdefault(key, [0, 0])[int(a == transpose(a))] += 1
            need(tau(child) == induced(tau(a), tuple(range(1, n + 1))), "growth tau")
            need(transpose(child) == induced(transpose(a), tuple(range(1, n + 1))), "growth sigma")
            growth_cases += 1
        for key, (black, blue) in counts.items():
            child = tournament(n, key)
            expected_blue = (1 << n) if child == transpose(child) else 0
            need(blue == expected_blue and black + blue == 1 << (2 * n - 1), "growth colour")
    print("Inherited Mode-B per-parent colour controls (target n=5..7) =", growth_cases)
    print("PASS: all checks active under normal and optimized Python.")


if __name__ == "__main__":
    main()

"""Minimum arc repair for a pair-module halving of Q[C3,C3,C3,1].

Exact finite universe: all 64 labeled four-vertex cores and all 945
perfect matchings of ten labeled vertices. No external dependencies.
"""
from collections import Counter
from itertools import combinations, permutations


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def matchings(vertices):
    if not vertices:
        yield ()
        return
    a = vertices[0]
    for i in range(1, len(vertices)):
        for rest in matchings(vertices[1:i] + vertices[i + 1:]):
            yield ((a, vertices[i]),) + rest


PAIRS4 = tuple(combinations(range(4), 2))
PAIRS5 = tuple(combinations(range(5), 2))
MATCHINGS10 = tuple(matchings(tuple(range(10))))
MATCHINGS4 = tuple(matchings(tuple(range(4))))


def core(mask):
    q = [[0] * 4 for _ in range(4)]
    for bit, (a, b) in enumerate(PAIRS4):
        q[a][b] = (mask >> bit) & 1
        q[b][a] = 1 - q[a][b]
    return q


def group(a):
    return a // 3 if a < 9 else 3


def expand(q):
    a = [[0] * 10 for _ in range(10)]
    for x, y in combinations(range(10), 2):
        bit = int((y - x) % 3 == 1) if group(x) == group(y) else q[group(x)][group(y)]
        a[x][y], a[y][x] = bit, 1 - bit
    return a


def cost(a, matching):
    answer = 0
    for (u, v), (x, y) in combinations(matching, 2):
        k = a[u][x] + a[u][y] + a[v][x] + a[v][y]
        answer += min(k, 4 - k)
    return answer


def repair(a, matching, target=None):
    b = [row[:] for row in a]
    flips = []
    for index, (left, right) in enumerate(combinations(matching, 2)):
        k = sum(a[x][y] for x in left for y in right)
        # With no prescribed quotient, use a fixed majority/tie gauge.
        direction = int(k >= 2) if target is None else target[index]
        for x in left:
            for y in right:
                if a[x][y] != direction:
                    flips.append(tuple(sorted((x, y))))
                b[x][y], b[y][x] = direction, 1 - direction
    # Verify modules from individual outside vertices, not 2x2 sums.
    for x, y in matching:
        for z in range(10):
            if z not in (x, y):
                check(b[x][z] == b[y][z], ("not module", matching, x, y, z))
    actual = sum(a[x][y] != b[x][y] for x, y in combinations(range(10), 2))
    expected = cost(a, matching) if target is None else sum(
        4 - k if bit else k for bit, k in zip(target, block_counts(a, matching)))
    check(actual == len(flips) == expected, "cost realization")
    check(all(b[x][y] == a[x][y] for x, y in matching), "internal arcs")
    return b, sorted(flips)


def block_counts(a, matching):
    return tuple(sum(a[x][y] for x in left for y in right)
                 for left, right in combinations(matching, 2))


def regular_quotients():
    regular = []
    for mask in range(1 << 10):
        word = tuple((mask >> i) & 1 for i in range(10))
        degree = [0] * 5
        for bit, (x, y) in zip(word, PAIRS5):
            degree[x if bit else y] += 1
        if degree == [2] * 5:
            regular.append(word)
    check(len(regular) == 24, "regular5 count")
    # Independent orbit of the circulant regular tournament H5.
    orbit = {tuple(int((p[y] - p[x]) % 5 in (1, 2)) for x, y in PAIRS5)
             for p in permutations(range(5))}
    check(set(regular) == orbit, "regular5 unique isomorphism type")
    return tuple(regular)


def internal_pairs(matching):
    return sum(group(x) == group(y) for x, y in matching)


def macro_cost(q, matching):
    defect = 0
    for i in range(3):
        opposite = next(pair for pair in matching if i not in pair)
        defect += q[i][opposite[0]] != q[i][opposite[1]]
    left, right = matching
    k = sum(q[x][y] for x in left for y in right)
    return 3 + 2 * defect + min(k, 4 - k)


def natural_lift(matching):
    rep = (2, 5, 8, 9)
    return ((0, 1), (3, 4), (6, 7)) + tuple(
        (rep[x], rep[y]) for x, y in matching)


def cyclic_triple(q, vertices):
    return all(sum(q[x][y] for y in vertices) == 1 for x in vertices)


def rooted_key(q):
    # Root3 fixed; permute three equal C3 copies. Include global reversal.
    words = []
    for p in permutations(range(3)):
        order = p + (3,)
        word = tuple(q[order[x]][order[y]] for x, y in PAIRS4)
        words.extend((word, tuple(1 - x for x in word)))
    return min(words)


def main():
    check(len(MATCHINGS10) == 945 and len(set(MATCHINGS10)) == 945,
          "matching universe")
    print("Universe:64 labeled cores,945 perfect matchings per core; total60480")
    distribution, rows, gauge = Counter(), [], {}
    macro_total = 0
    for mask in range(64):
        q = core(mask)
        a = expand(q)
        values = [(cost(a, m), m) for m in MATCHINGS10]
        best = min(c for c, _ in values)
        opt = [m for c, m in values if c == best]
        check(all(internal_pairs(m) == 3 for m in opt),
              ("nonnatural optimum", mask))
        macros = [(macro_cost(q, m), m) for m in MATCHINGS4]
        for c, m in macros:
            check(c == cost(a, natural_lift(m)), ("macro formula", mask, m))
        check(min(c for c, _ in macros) == best, ("macro minimum", mask))
        macro_opt = [m for c, m in macros if c == best]
        check(len(opt) == 27 * len(macro_opt), ("rotational count", mask))
        macro_total += len(macro_opt)
        # Verify every optimal matching has an actual repaired tournament.
        for matching in opt:
            repair(a, matching)
        scores = tuple(sorted(map(sum, q)))
        deleted_cyclic = cyclic_triple(q, (0, 1, 2))
        key = rooted_key(q)
        invariant = (best, len(opt))
        check(key not in gauge or gauge[key] == invariant,
              ("rooted relabel/dual gauge", mask))
        gauge[key] = invariant
        distribution[best] += 1
        rows.append((mask, scores, sum(q[3]), deleted_cyclic,
                     best, len(opt), macro_opt, opt[0]))
    check(distribution == {3: 24, 6: 24, 7: 12, 8: 4}, "minimum census")
    print("minimum reversal distribution:", sorted(distribution.items()))
    print("root-preserving relabeling + global reversal orbits:", len(gauge))
    print("ALL optimum matchings have exactly3 internal C3 pairs")
    print("ALL optimum counts =27 times optimum macro matchings")
    print("macro formula and realized repairs: PASS")
    print("columns: mask,sorted_core_outdegrees,root_outdegree,root_deleted_C3,")
    print("         minimum,#optimal_matchings,optimal_macro_matchings,example")
    for row in rows:
        print(row)
    for mask in (0, 4, 5, 6, 2):
        row = rows[mask]
        q = core(mask)
        repaired, flips = repair(expand(q), row[-1])
        matching = row[-1]
        quotient = [[0] * 5 for _ in range(5)]
        for i, j in combinations(range(5), 2):
            quotient[i][j] = repaired[matching[i][0]][matching[j][0]]
            quotient[j][i] = 1 - quotient[i][j]
        print("representative", mask, "flip_edges", flips,
              "quotient_outdegrees", sorted(map(sum, quotient)))
    # Positive control: a transitive tournament of even order already halves.
    tt10 = [[int(x < y) for y in range(10)] for x in range(10)]
    adjacent = tuple((i, i + 1) for i in range(0, 10, 2))
    check(cost(tt10, adjacent) == 0, "TT10 positive control")
    check(min(cost(tt10, m) for m in MATCHINGS10) == 0, "TT10 full census")
    print("TT10 positive control:minimum0; F_Q(C3) all64 hostile:min>=3")
    targets = regular_quotients()
    print("CANONICAL QUOTIENT:24 labeled regular5 tournaments,one iso type")
    print("Universe64*945*24=1451520 matching/quotient candidates")
    canonical_distribution = Counter()
    canonical_gauge = {}
    canonical_examples = {}
    for mask in range(64):
        q = core(mask)
        a = expand(q)
        best, natural_best, opt = 99, 99, []
        for matching in MATCHINGS10:
            counts = block_counts(a, matching)
            natural = internal_pairs(matching) == 3
            for target in targets:
                value = sum(4 - k if bit else k for bit, k in zip(target, counts))
                if natural:
                    natural_best = min(natural_best, value)
                if value < best:
                    best, opt = value, [(matching, target)]
                elif value == best:
                    opt.append((matching, target))
        unique_matchings = {matching for matching, _ in opt}
        profile = tuple(sorted(Counter(internal_pairs(m) for m in unique_matchings).items()))
        # Realize every optimum and independently inspect quotient outdegrees.
        for matching, target in opt:
            repaired, flips = repair(a, matching, target)
            degree = [sum(repaired[pair[0]][other[0]] for other in matching
                          if other != pair) for pair in matching]
            check(degree == [2] * 5, ("canonical quotient degree", mask))
            check(len(flips) == best, "canonical minimum realization")
        key = rooted_key(q)
        invariant = (best, natural_best, len(unique_matchings), len(opt), profile)
        check(key not in canonical_gauge or canonical_gauge[key] == invariant,
              ("canonical rooted gauge", mask))
        canonical_gauge[key] = invariant
        canonical_distribution[(rows[mask][4], best)] += 1
        print("canonical", mask, "unrestricted", rows[mask][4], "minimum", best,
              "natural_minimum", natural_best, "optimal_matchings", len(unique_matchings),
              "optimal_matching_quotient_pairs", len(opt), "internal_pair_profile", profile)
        if mask in (0, 2, 4, 5, 6):
            canonical_examples[mask] = (opt[0], repair(a, *opt[0])[1])
    check(canonical_distribution == {(3, 15): 24, (6, 13): 12, (6, 12): 12,
                                      (7, 7): 12, (8, 10): 4}, "canonical census")
    print("unrestricted->canonical distribution:", sorted(canonical_distribution.items()))
    for mask, ((matching, target), flips) in canonical_examples.items():
        print("canonical example", mask, "matching", matching,
              "quotient_bits", target, "flip_edges", flips)
    e10 = [[0] * 10 for _ in range(10)]
    target = targets[0]
    for bit, (i, j) in zip(target, PAIRS5):
        for x in (2 * i, 2 * i + 1):
            for y in (2 * j, 2 * j + 1):
                e10[x][y], e10[y][x] = bit, 1 - bit
    for i in range(5):
        e10[2 * i][2 * i + 1] = 1
    check(sum(4-k if bit else k for bit, k in
              zip(target, block_counts(e10, adjacent))) == 0, "E10 positive control")
    print("canonical E10 positive control:cost0")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()

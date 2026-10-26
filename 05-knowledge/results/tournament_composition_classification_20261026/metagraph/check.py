"""Independent finite controls for the marked composition proofs in RESULT.md.

No external packages, producer data, precomputed classes, or assert statements.
Free input bits use i-major nonconsecutive pairs. Canonical ordinary-class
integers use j-major upper-triangle bits; these are different encodings.
"""

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations, product
import json
import sys

sys.stdout.reconfigure(newline="\n")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def presentation(n, mask):
    a = [[0] * n for _ in range(n)]
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            forward = 1 if j == i + 1 else (mask >> bit) & 1
            bit += j > i + 1
            a[i][j], a[j][i] = forward, 1 - forward
    return tuple(map(tuple, a))


@lru_cache(None)
def universe(n):
    return tuple(presentation(n, m) for m in range(1 << max(0, (n - 1) * (n - 2) // 2))) if n else ((),)


def tau(a):
    n = len(a)
    return tuple(tuple(a[i][j] if abs(i - j) <= 1 else 1 - a[i][j]
                       for j in range(n)) for i in range(n))


def sigma(a):
    n = len(a)
    return tuple(tuple(a[n - 1 - j][n - 1 - i] for j in range(n)) for i in range(n))


def reindex(a, p):
    return tuple(tuple(a[i][j] for j in p) for i in p)


def join(a, b):
    m, n = len(a), len(b)
    return tuple(tuple(a[i][j] if i < m and j < m else
                       b[i - m][j - m] if i >= m and j >= m else int(i < m)
                       for j in range(m + n)) for i in range(m + n))


def hinge(a, b):
    """Direct cross-arc definition, independent of complement conjugation."""
    m, n = len(a), len(b)
    if not m:
        return b
    if not n:
        return a
    c = [list(row) for row in join(a, b)]
    for i in range(m):
        for j in range(m, m + n):
            c[i][j] = int(i == m - 1 and j == m)
            c[j][i] = 1 - c[i][j]
    return tuple(map(tuple, c))


def rooted_glue(a, u, b, v):
    """Arbitrary root labels, with no assumption that either root is an end."""
    m, n = len(a), len(b)
    c = [list(row) for row in join(a, b)]
    for i in range(m):
        for j in range(n):
            c[i][m + j] = int(i == u and j == v)
            c[m + j][i] = 1 - c[i][m + j]
    return tuple(map(tuple, c))


def relocations(a):
    """Full word reindexing, retaining move names and any loops."""
    n = len(a)
    ans = {}
    for i in range(n):
        for j in range(i + 2, n):
            p = tuple(range(i)) + tuple(range(i + 1, j + 1)) + (i,) + tuple(range(j + 1, n))
            if all(a[p[k]][p[k + 1]] for k in range(n - 1)):
                ans[i, j] = reindex(a, p)
    return ans


def c3(a):
    return sum(a[i][j] == a[j][k] == a[k][i]
               for i, j, k in combinations(range(len(a)), 3))


def strong_intervals(a):
    n = len(a)
    reach = [sum(a[i][j] << j for j in range(n)) | (1 << i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i] & (1 << k):
                reach[i] |= reach[k]
    blocks = []
    unused = set(range(n))
    while unused:
        i = min(unused)
        block = tuple(j for j in sorted(unused) if reach[i] & (1 << j) and reach[j] & (1 << i))
        require(block == tuple(range(block[0], block[-1] + 1)), "SCC is not a contiguous path interval")
        blocks.append(block)
        unused.difference_update(block)
    return blocks


def fold(operation, factors):
    answer = ()
    for factor in factors:
        answer = operation(answer, factor)
    return answer


def word_code(a, p):
    return sum(a[p[i]][p[j]] << (j * (j - 1) // 2 + i)
               for j in range(len(a)) for i in range(j))


@lru_cache(None)
def canon(a):
    return min(word_code(a, p) for p in permutations(range(len(a))))


def rooted_canon(a, root):
    others = tuple(i for i in range(len(a)) if i != root)
    return len(a), min(word_code(a, (root,) + p) for p in permutations(others))


def h_dp(a):
    n = len(a)
    table = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        table[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            if table[mask][v]:
                for w in range(n):
                    if not mask & (1 << w) and a[v][w]:
                        table[mask | (1 << w)][w] += table[mask][v]
    return sum(table[-1])


def h_words(a):
    return sum(all(a[p[i]][p[i + 1]] for i in range(len(a) - 1))
               for p in permutations(range(len(a))))


def aut(a):
    return sum(reindex(a, p) == a for p in permutations(range(len(a))))


def c4(a):
    return sum(all(a[p[i]][p[(i + 1) % 4]] for i in range(4))
               for subset in combinations(range(len(a)), 4)
               for tail in permutations(subset[1:]) for p in [(subset[0],) + tail])


def factor_controls():
    pairs = move_instances = triples = 0
    for m in range(1, 7):
        for n in range(1, 7):
            if m + n > 8:
                continue
            for a, b in product(universe(m), universe(n)):
                z = join(a, b)
                expected = {(i, j): join(q, b) for (i, j), q in relocations(a).items()}
                expected.update({(i + m, j + m): join(a, q) for (i, j), q in relocations(b).items()})
                actual = relocations(z)
                require(actual == expected, "join relocation product failed")
                w = hinge(a, b)
                require(w == tau(join(tau(a), tau(b))), "direct hinge/conjugate mismatch")
                require(sigma(z) == join(sigma(b), sigma(a)), "join reflection failed")
                require(sigma(w) == hinge(sigma(b), sigma(a)), "hinge reflection failed")
                require(c3(w) == c3(a) + c3(b) + sum(row[-1] for row in a) + sum(b[0]), "two-factor cyclic formula failed")
                require(m + n == 2 or len(strong_intervals(w)) == 1, "hinge is not strong")
                pairs += 1
                move_instances += len(actual)
    for sizes in product(range(1, 5), repeat=3):
        if sum(sizes) > 7:
            continue
        for factors in product(*(universe(n) for n in sizes)):
            a, b, c = factors
            require(join(join(a, b), c) == join(a, join(b, c)), "join associativity failed")
            w = hinge(hinge(a, b), c)
            require(w == hinge(a, hinge(b, c)), "hinge associativity failed")
            count = sum(map(c3, factors))
            count += sum(sum(row[-1] for row in x) + sum(y[0]) for x, y in zip(factors, factors[1:]))
            count += int(sizes[1] == 1)
            require(c3(w) == count, "three-factor cyclic formula failed")
            triples += 1
    for n in range(2, 16):
        require(c3(fold(hinge, [universe(1)[0]] * n)) == n - 2, "singleton boundary failed")
    return {"pair_size_bounds": "1<=a,b<=6; a+b<=8", "pair_presentations": pairs,
            "named_relocation_instances": move_instances,
            "triple_size_bounds": "1<=a,b,c<=4; a+b+c<=7", "triple_presentations": triples,
            "singleton_products_through_order": 15}


def decomposition_controls():
    strong_counts, total_counts = [], [1]
    tested = 0
    for n in range(1, 8):
        strong = 0
        for a in universe(n):
            blocks = strong_intervals(a)
            factors = [reindex(a, block) for block in blocks]
            require(fold(join, factors) == a, "SCC join reconstruction failed")
            require(all(len(strong_intervals(f)) == 1 for f in factors), "SCC factor not strong")
            b = tau(a)
            dual_factors = [tau(reindex(b, block)) for block in strong_intervals(b)]
            require(fold(hinge, dual_factors) == a, "star factor reconstruction failed")
            require(tau(tau(a)) == a and sigma(tau(a)) == tau(sigma(a)), "involution control failed")
            strong += len(blocks) == 1
            tested += 1
        strong_counts.append(strong)
        total_counts.append(len(universe(n)))
        require(total_counts[n] == sum(strong_counts[k - 1] * total_counts[n - k] for k in range(1, n + 1)),
                "free join monoid count failed")
    require(strong_counts == [1, 0, 1, 5, 50, 903, 30773], "unexpected strong presentation control")
    return {"all_X_n_through_order": 7, "presentations": tested, "strong_presentation_counts": strong_counts}


def minimality_and_hostile():
    counts = {}
    for total in range(2, 5):
        count = 0
        for m in range(1, total):
            n = total - m
            seen = {}
            for a, b in product(universe(m), universe(n)):
                key = canon(a), canon(tau(a)), canon(b), canon(tau(b))
                target = canon(tau(join(a, b)))
                require(key not in seen or seen[key] == target, "smaller coarse composition obstruction")
                seen[key] = target
                count += 1
        counts[total] = count
    x, y = presentation(4, 2), presentation(4, 5)
    e = presentation(1, 0)
    require(tau(x) == y and sigma(x) == x and sigma(y) == y, "hostile is not one blue line")
    require(reindex(y, (1, 3, 0, 2)) == x, "explicit hostile isomorphism failed")
    require(canon(x) == canon(y) and aut(x) == aut(y) == 1, "factor class or automorphism control failed")
    require(h_dp(x) == h_words(x) == h_dp(y) == h_words(y) == 5, "factor Hamiltonian control failed")
    require(sum(tau(x)[0]) == sum(tau(y)[0]) == 2, "root degree control failed")
    require(rooted_canon(tau(x), 0) != rooted_canon(tau(y), 0), "rooted ports failed to separate hostile")
    outputs = [tau(join(e, x)), tau(join(e, y))]
    facts = [{"canonical_class": canon(a), "H_dp": h_dp(a), "H_words": h_words(a),
              "automorphism_order": aut(a), "cyclic_triples": c3(a), "directed_four_cycles": c4(a),
              "sorted_outdegrees": sorted(map(sum, a))} for a in outputs]
    require([f["H_dp"] for f in facts] == [11, 13], "hostile H values failed")
    require(all(f["H_dp"] == f["H_words"] for f in facts), "independent H paths disagree")
    require([f["directed_four_cycles"] for f in facts] == [4, 3], "hostile four-cycle values failed")
    require(all(f["cyclic_triples"] == 4 and f["sorted_outdegrees"] == [1, 2, 2, 2, 3]
                and f["automorphism_order"] == 1 for f in facts), "hostile coarse sidecars failed")
    return {"all_splits_below_order_five": counts, "witness_total_order": 5,
            "factor_masks": [2, 5], "one_literal_blue_loop": True, "outputs": facts}


def histogram_controls():
    class_fibres = {}
    for n in range(1, 5):
        fibres = defaultdict(list)
        for a in universe(n):
            fibres[canon(a)].append(a)
        for members in fibres.values():
            require(len(members) * aut(members[0]) == h_dp(members[0]), "presentation mass H/Aut failed")
        class_fibres[n] = fibres
    class_pairs = presentation_pairs = 0
    for m in range(1, 5):
        for n in range(1, 5):
            if m + n > 6:
                continue
            for left in class_fibres[m].values():
                ends, end_reps = Counter(), {}
                for a in left:
                    c = tau(a)
                    key = rooted_canon(c, m - 1)
                    ends[key] += 1
                    end_reps[key] = c, m - 1
                for right in class_fibres[n].values():
                    starts, start_reps = Counter(), {}
                    for b in right:
                        d = tau(b)
                        key = rooted_canon(d, 0)
                        starts[key] += 1
                        start_reps[key] = d, 0
                    convolution = Counter()
                    for ckey, ecount in ends.items():
                        c, u = end_reps[ckey]
                        for dkey, bcount in starts.items():
                            d, v = start_reps[dkey]
                            convolution[canon(rooted_glue(c, u, d, v))] += ecount * bcount
                    direct = Counter(canon(tau(join(a, b))) for a, b in product(left, right))
                    require(convolution == direct, "rooted multiplicity convolution failed")
                    require(sum(convolution.values()) == len(left) * len(right), "histogram total mass failed")
                    class_pairs += 1
                    presentation_pairs += len(left) * len(right)
    return {"factor_orders": "1..4; total<=6", "ordinary_class_pairs": class_pairs,
            "presentation_pairs": presentation_pairs, "target_multiplicities_equal": True}


def main():
    result = {"status": "FINITE-EXACT",
              "scope": "bounded controls of all-order proofs; ordinary isomorphism, no converse quotient",
              "carrier": "fixed forward Hamiltonian path; free bits i-major; all presentations retained",
              "operations": factor_controls(),
              "factorizations": decomposition_controls(),
              "coarse_composition_hostile": minimality_and_hostile(),
              "rooted_histogram_convolution": histogram_controls()}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

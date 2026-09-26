"""Exact finite audit of order laws and the Collatz height-conditioning loss.

No third-party dependencies. Run with Python 3.11+; all checks survive -O.
The paper's asymptotic theorem is not established by this finite census.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import permutations, product
from math import factorial


def check(predicate, message):
    if not predicate:
        raise RuntimeError(message)


def natural_posets(n):
    """Every distinct transitive relation with i<j in its fixed topological order."""
    pairs = [(i, j) for j in range(n) for i in range(j)]
    seen = set()
    for mask in range(1 << len(pairs)):
        pred = [0] * n
        for k, (i, j) in enumerate(pairs):
            if mask >> k & 1:
                pred[j] |= 1 << i
        for j in range(n):
            for i in range(j):
                if pred[j] >> i & 1:
                    pred[j] |= pred[i]
        seen.add(tuple(pred))
    return sorted(seen)


def extensions(pred):
    n = len(pred)
    out = []

    def visit(mask, order):
        if len(order) == n:
            out.append(tuple(order))
            return
        for i in range(n):
            if not (mask >> i & 1) and pred[i] & mask == pred[i]:
                visit(mask | (1 << i), order + [i])

    visit(0, [])
    return out


def rank_moments(orders):
    """Numerators for Cov(Z_i,Z_j), Z=(n+1)*uniform order coordinates."""
    n, m = len(orders[0]), len(orders)
    heights = [0] * n
    moments = [[0] * n for _ in range(n)]
    rank_rows = []
    for order in orders:
        ranks = [0] * n
        for k, v in enumerate(order, 1):
            ranks[v] = k
        rank_rows.append(ranks)
        for i in range(n):
            heights[i] += ranks[i]
            for j in range(n):
                moments[i][j] += min(ranks[i], ranks[j]) * (max(ranks[i], ranks[j]) + 1)
    denominator = (n + 2) * m * m
    cov = [[(n + 1) * m * moments[i][j] - (n + 2) * heights[i] * heights[j]
            for j in range(n)] for i in range(n)]
    return heights, cov, denominator, rank_rows


def xyz(orders, x, y, z):
    _, cov, denominator, _ = rank_moments(orders)
    return F(cov[y][z] - cov[y][y] - cov[x][z] + cov[x][y], denominator)


def census_posets():
    totals = Counter()
    expected = [1, 2, 7, 40, 357, 4824]
    for n in range(1, 7):
        posets = natural_posets(n)
        check(len(posets) == expected[n - 1], "natural-poset count")
        for pred in posets:
            orders = extensions(pred)
            m = len(orders)
            h, cov, denominator, ranks = rank_moments(orders)
            gap = [sum(row[i] - max([0] + [row[j] for j in range(n) if pred[i] >> j & 1])
                       for row in ranks) for i in range(n)]
            for x, y, z in permutations(range(n), 3):
                check(cov[y][z] - cov[y][y] - cov[x][z] + cov[x][y] <= 0, "XYZ covariance")
                totals['xyz_triples'] += 1
            for x in range(n):
                for y in range(x):
                    variance = cov[x][x] + cov[y][y] - 2 * cov[x][y]
                    check(3 * variance >= (n + 2) * max(gap[x], gap[y]) ** 2,
                          "pair variance/window bound")
                    totals['pair_window_bounds'] += 1
            for mask in range(1, 1 << n):
                vertices = [i for i in range(n) if mask >> i & 1]
                if all(pred[i] & mask == pred[i] for i in vertices):
                    maxima = [i for i in vertices if not any(pred[j] >> i & 1 for j in vertices)]
                    check(max(h[i] for i in vertices) >= m * (len(vertices) - len(maxima) + 1),
                          "Haqi ideal mean-height bound")
                    totals['ideal_height_bounds'] += 1
                if all(not (pred[i] & mask) for i in vertices):
                    check(2 * n * sum(gap[i] for i in vertices) >=
                          m * len(vertices) ** 2 * (n + 1), "antichain window sum")
                    totals['antichain_bounds'] += 1
            totals['posets'] += 1
            totals['linear_extensions'] += m
        print(f"naturally_labelled_posets n={n}: {len(posets)}")
    print("finite_poset_lemma_census:", dict(sorted(totals.items())))


def insertion_controls():
    # An isolated x=0 and a t-chain 1<...<t give exactly the insertion law.
    for t in range(1, 7):
        pred = (0,) + tuple(sum(1 << j for j in range(1, i)) for i in range(1, t + 1))
        orders = extensions(pred)
        check(len(orders) == t + 1, "insertion extension count")
        for i in range(1, t + 1):
            check(sum(o.index(0) < o.index(i) for o in orders) == i,
                  "insertion comparison probability")
        distance = F(1) - F(t + 1, factorial(t + 1))
        print(f"exact_insertor t={t}: extensions={t + 1}, TV_from_uniform_orders={distance}")
    # A genuine antichain restriction need not have a uniform order law.
    v_orders = extensions((0, 0, 1))  # a<c, b isolated
    p = F(sum(o.index(0) < o.index(1) for o in v_orders), len(v_orders))
    check(p == F(2, 3), "antichain restriction hostile")
    print("antichain_restriction hostile: P(a before b | a<c,b isolated)=", p)


def parity_word(n, length):
    word = []
    for _ in range(length):
        odd = n & 1
        word.append(odd)
        n = (3 * n + 1) // 2 if odd else n // 2
    return tuple(word)


def residue(word):
    exponent, carry = 0, 0
    for j, bit in enumerate(word):
        if bit:
            exponent += 1
            carry = 3 * carry + (1 << j)
    modulus = 1 << len(word)
    return (-carry * pow(3 ** exponent, -1, modulus)) % modulus


def positive_prefix(word):
    exponent = 0
    for j, bit in enumerate(word, 1):
        exponent += bit
        if 3 ** exponent <= 2 ** j:
            return False
    return True


def event_order(word):
    a = b = 0
    ones = sum(word)
    order = []
    for bit in word:
        if bit:
            order.append(a)
            a += 1
        else:
            order.append(ones + b)
            b += 1
    return tuple(order)


def unanimous_predecessors(orders):
    n = len(orders[0])
    rank_rows = [{v: i for i, v in enumerate(o)} for o in orders]
    return tuple(sum(1 << a for a in range(n) if a != b and
                     all(r[a] < r[b] for r in rank_rows)) for b in range(n))


def law_tv(first, second):
    m, n = sum(first.values()), sum(second.values())
    return sum((abs(F(first[k], m) - F(second[k], n)) for k in first.keys() | second.keys()), F(0)) / 2


def parity_laws():
    good = {}
    checked = 0
    for length in range(1, 13):
        good[length] = []
        for word in product((0, 1), repeat=length):
            r = residue(word)
            check(parity_word(r, length) == word, "carry inverse disagrees with direct orbit")
            checked += 1
            if positive_prefix(word):
                good[length].append(word)
    print("parity_inverse/direct_checks:", checked)
    for length in range(1, 12):
        current = Counter(good[length])
        projected = Counter(w[:-1] for w in good[length + 1])
        print(f"conditional_good_words L={length}: count={sum(current.values())}, "
              f"TV_from_next_level_projection={law_tv(current, projected)}")
    check(set(good[3]) == {(1, 1, 0), (1, 1, 1)}, "level-three witness")
    check(set(good[4]) == {(1, 1, 0, 1), (1, 1, 1, 0), (1, 1, 1, 1)}, "level-four witness")
    check(law_tv(Counter(good[3]), Counter(w[:-1] for w in good[4])) == F(1, 6),
          "projective-law hostile")

    # Exhaustive minimality in this explicitly bounded universe, not all posets.
    first = None
    tested = 0
    for length in range(1, 7):
        groups = defaultdict(list)
        for w in good[length]:
            groups[sum(w)].append((residue(w), w, event_order(w)))
        for ones, rows in sorted(groups.items()):
            rows.sort()
            for m in range(1, len(rows)):
                chosen = [row[2] for row in rows[:m]]
                restored = extensions(unanimous_predecessors(chosen))
                tested += 1
                if set(restored) != set(chosen) and first is None:
                    first = length, ones, rows[m - 1][0], rows[:m], restored
    check(first is not None and first[:3] == (6, 5, 31), "height-fibre minimal witness")
    length, ones, cutoff, selected, restored = first
    rows = sorted((residue(w), w, event_order(w)) for w in good[6] if sum(w) == 5)
    chosen = [row[2] for row in selected]
    check({r for r, _, _ in rows} == {27, 31, 39, 47}, "height-fibre residue values")
    check({r for r, _, _ in selected} == {27, 31}, "height-fibre selected starts")
    check(set(restored) == {row[2] for row in rows}, "unanimous completion has all four orders")
    check(xyz(chosen, 2, 1, 5) == F(1, 4), "positive XYZ covariance hostile")
    check(xyz(restored, 2, 1, 5) <= 0, "unconditioned XYZ positive control")
    p = tuple(F(v, 10) for v in (1, 2, 4, 5, 6, 3))
    qpoint = tuple(F(v, 10) for v in (1, 2, 3, 4, 5, 7))
    midpoint = tuple((u + v) / 2 for u, v in zip(p, qpoint))
    coordinate_order = lambda values: tuple(sorted(range(6), key=values.__getitem__))
    check(coordinate_order(p) in chosen and coordinate_order(qpoint) in chosen,
          "nonconvex support endpoint controls")
    check(coordinate_order(midpoint) not in chosen and coordinate_order(midpoint) in restored,
          "nonconvex support midpoint hostile")
    print("height_fibre_minimality_universe: L<=6, all endpoint counts, least-residue initial cuts; tests=", tested)
    print("height_fibre_rows:", [(r, ''.join(map(str, w))) for r, w, _ in rows])
    print("height_fibre_selected: n<=31 gives27,31; unanimous_poset_completion_count=", len(restored))
    print("height_fibre_XYZ rank units:", xyz(chosen, 2, 1, 5),
          "; full-four-order control:", xyz(restored, 2, 1, 5))
    print("height_fibre_XYZ original unit coordinates:", xyz(chosen, 2, 1, 5) / 49)
    print("height_fibre_nonconvex_support: two strict chamber points have midpoint in missing chamber")

    # Positive control: restore the full law over many complete residue periods.
    # Compare the closed formula with actual integer trajectories, not with a
    # second rearrangement of the residue-count formula.
    cutoff_checks = 0
    for length in range(1, 11):
        modulus = 1 << length
        groups = defaultdict(list)
        for w in good[length]:
            groups[sum(w)].append(w)
        cutoffs = sorted({modulus // 3, modulus // 2, modulus - 1, modulus,
                          modulus + 1, 2 * modulus - 1, 3 * modulus + modulus // 3})
        for cutoff in cutoffs:
            if cutoff < 1:
                continue
            actual = Counter(parity_word(n, length) for n in range(1, cutoff + 1))
            q, s = divmod(cutoff, modulus)
            for words in groups.values():
                filtered = Counter({w: actual[w] for w in words})
                if not sum(filtered.values()):
                    continue
                m = len(words)
                t = sum(1 <= residue(w) <= s for w in words)
                exact_tv = F(t * (m - t), m * (q * m + t))
                check(law_tv(filtered, Counter(words)) == exact_tv, "periodic TV recovery")
                if q:
                    check(exact_tv <= F(1, 4 * q), "uniform many-period TV bound")
                cutoff_checks += 1
    print("many_period_uniformity_direct_integer_checks:", cutoff_checks)


def main():
    print("STATUS FINITE-EXACT; the attached asymptotic claim remains UNDER REVIEW")
    census_posets()
    insertion_controls()
    parity_laws()
    print("PASS: all exact checks completed; no assertion-disabled checks")


if __name__ == '__main__':
    main()

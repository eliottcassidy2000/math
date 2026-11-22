"""Exact local branching, inverse fibres, and Paley stabilizer counts."""
from collections import Counter
from math import comb

from grammar import (ROOT, canon, capped_masks, cut, endpoint_cut, strong_classes,
                     strongly_connected)


def divisor_list(a):
    return [j for j in range(1, a + 1) if a % j == 0]


def cyclic_sizes(p, k):
    """Counts of Paley children at both complementary mask weights by |Aut|."""
    m = (p - 1) // 2
    div = divisor_list(m)
    exact = {}
    for d in reversed(div[1:]):
        q = (p - 1) // d
        fixed = sum(comb(q, (k - e) // d) for e in (0, 1)
                    if (k - e) % d == 0 and 0 <= (k - e) // d <= q)
        exact[d] = p * fixed - sum(exact[t] for t in exact if t % d == 0)
    exact[1] = comb(p, k) - sum(exact.values())
    out = {d: 2 * masks // (p * m // d) for d, masks in exact.items()}
    if k == m:
        out[m] -= 1
    return {d: n for d, n in sorted(out.items()) if n}


def paley(p):
    squares = {j*j % p for j in range(1, p)}
    return tuple(sum(1 << ((x+s) % p) for s in squares) for x in range(p))


def local_branching(n):
    directory = ROOT / 'session212/cycle-levels'
    strong = 0
    fewest = 1 << n
    largest = strict = 0
    hostile = None
    for row in (directory / f'classes{n}.rows').read_text().splitlines():
        rows = tuple(map(int, row.split()))
        full = (1 << n) - 1
        if not strongly_connected(rows, full):
            continue
        strong += 1
        h, weights = endpoint_cut(rows)
        values = [h + 2*cut(weights, mask) for mask in range(full+1)]
        fewest = min(fewest, len(set(values)))
        largest = max(largest, max(Counter(values).values()))
        best = min(values[1:full])
        singleton = min(min(values[1 << i], values[full ^ (1 << i)]) for i in range(n))
        if best < singleton:
            strict += 1
            if hostile is None:
                hostile = (rows, h, best, singleton)
        if n <= 5:
            for ceiling in range(min(values), max(values)+2, 2):
                filtered = dict(capped_masks(weights, (ceiling-h)//2))
                assert set(filtered) == {m for m in range(1, full) if values[m] <= ceiling}
    print('order/strong/min distinct extension counts/max masks at one count/'
          'parents with all minimum strong cuts nonsingleton:', n, strong, fewest, largest, strict)
    if hostile:
        print('first hostile rows, H, minimum, best singleton:', hostile)


def small_fibres():
    wanted = (3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 33)
    inventory = Counter()
    nine = []
    for n, code, h in strong_classes(max(wanted)):
        inventory[h, n] += 1
        if h == 9:
            rows = canon.rows(n, code)
            deck = sorted(canon.hamilton_count(n - 1, canon.encode(rows,
                          [i for i in range(n) if i != v])) for v in range(n))
            reverse = code ^ ((1 << (n * (n - 1) // 2)) - 1)
            assert canon.canonical(n, reverse).code == code
            assert len(canon.canonical(n, code).orders) == 1
            nine.append(deck)
    assert sorted(nine) == [[1, 1, 5, 5, 5], [1, 3, 3, 5, 5]]
    for h in wanted:
        print('all strong orders, H', h, [(n, inventory[h,n]) for n in range(3, 8)
                                           if inventory[h, n]])


def first_paley():
    from csv import DictReader
    path = ROOT / '04-computation/experiments/tiling_line_metagraph_n9_20261020_classes_n8.tsv'
    with path.open() as src:
        all_eight = list(DictReader(src, delimiter='\t'))
    for p in (3, 7, 11, 19, 23, 31):
        m = (p-1)//2
        print('Paley central H-fibre via parent, p / classes by automorphism order',
              p, cyclic_sizes(p, m))
    p = 7
    base = paley(p)
    values = {}
    for mask in range(1, (1<<p)-1):
        n = p+1
        code = canon.extend(p, canon.encode(base), ((1<<p)-1)^mask)
        cf = canon.canonical(n, code)
        h = canon.hamilton_count(n, cf.code)
        values.setdefault(h, {})[cf.code] = len(cf.orders)
    for k, (h, group) in enumerate(sorted(values.items()), start=1):
        local = dict(sorted(Counter(group.values()).items()))
        assert local == cyclic_sizes(p, k)
        complete = [r for r in all_eight if int(r['H']) == h]
        print('first-order-8 H / Paley-card automorphisms / all classes',
              h, local, len(complete))


if __name__ == '__main__':
    small_fibres()
    for n in range(3, 8):
        local_branching(n)
    first_paley()

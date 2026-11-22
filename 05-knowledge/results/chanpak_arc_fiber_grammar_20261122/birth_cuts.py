"""Find pairs of strong component heights that cannot coexist at any first birth."""

from collections import defaultdict
from itertools import combinations_with_replacement, product

from fibers import census_rows


def main():
    birth, by_order = {}, defaultdict(list)
    for n in range(3, 12):
        for h, _, _, _ in census_rows(n):
            if h > 1 and h not in birth:
                birth[h] = n
                by_order[n].append(h)

    for size in range(6, 15):
        pairs = compressed = 0
        minimum_saving = size
        first_unknown = None
        for n in range(3, size // 2 + 1):
            m = size - n
            options = (combinations_with_replacement(by_order[n], 2) if n == m
                       else product(by_order[n], by_order[m]))
            for a, b in options:
                pairs += 1
                product_order = birth.get(a * b)
                if product_order is None or product_order >= size:
                    if first_unknown is None:
                        first_unknown = (a, b)
                    continue
                compressed += 1
                minimum_saving = min(minimum_saving, size - product_order)
        print(size, 'pairs', pairs, 'strictly compressed', compressed,
              'least certified saving', minimum_saving, 'first undecided', first_unknown)
        if size <= 13:
            assert pairs == compressed
        else:
            assert (pairs, compressed) == (51818, 22565)
    assert (birth[3], birth[9]) == (3, 5)


if __name__ == '__main__':
    main()

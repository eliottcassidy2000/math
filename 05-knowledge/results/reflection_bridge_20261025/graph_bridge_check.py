"""Exact controls for graph_bridge.md; explicit checks also run under -O."""
import argparse
from collections import Counter
from itertools import combinations, permutations, product
import json
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def divisors(k):
    return [x for x in product(range(k + 1), range(2), range(2))
            if x not in ((0, 0, 0), (k, 1, 1))]


def complement(k, x):
    a, b, c = x
    return k-a, 1-b, 1-c


def swap(x):
    a, b, c = x
    return a, c, b


def divides(x, y):
    return all(a <= b for a, b in zip(x, y))


def graph(k, x):
    a, b, c = x
    if k == 1:
        return {e for e, bit in zip(((0, 1), (0, 2), (1, 2)), (b, a, c)) if bit}
    require(k == 3, 'map defined only at k=1 or k=3')
    bits = (a+2*b >= 3, b, a >= 2, a+b+c >= 3, c, a+2*c >= 3)
    return {e for e, bit in zip(combinations(range(4), 2), bits) if bit}


def reflect(m, edges):
    return {tuple(sorted((m-1-i, m-1-j))) for i, j in edges}


def tournament(m, edges):
    return {(j, i) if (i, j-1) in edges else (i, j)
            for i, j in combinations(range(m+1), 2)}


def strong(n, arcs):
    for start in range(n):
        reached = {start}
        while True:
            grown = reached | {b for a, b in arcs if a in reached}
            if grown == reached:
                break
            reached = grown
        if len(reached) != n:
            return False
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    checks = Counter()
    rows = []
    for k, m in ((1, 3), (3, 4)):
        d = divisors(k)
        universe = set(combinations(range(m), 2))
        image = {x: graph(k, x) for x in d}
        require(len({tuple(sorted(v)) for v in image.values()}) == len(d), f'injection k={k}')
        for x in d:
            edges = image[x]
            require(image[complement(k, x)] == universe-edges, f'C law k={k}, x={x}')
            require(image[swap(x)] == reflect(m, edges), f'Q law k={k}, x={x}')
            checks['equivariance_states'] += 1
            arcs = tournament(m, edges)
            require(all((i, i+1) in arcs for i in range(m)), f'path lift k={k}, x={x}')
            require((x == swap(x)) == (edges == reflect(m, edges)), f'blue predicate {x}')
            if k == 3:
                a, b, c = x
                require(((m, 0) in arcs) == (a >= 2), f'closure predicate {x}')
                require(strong(m+1, arcs) == (a >= 2 or b == c == 1), f'strong predicate {x}')
                require((int((0, 2) in edges)-int((1, 3) in edges)) == b-c, f'g predicate {x}')
                require((not edges) == (x == (1, 0, 0)), f'transitivity predicate {x}')
                code = sum(int(e in edges) << i for i, e in enumerate(combinations(range(m), 2)))
                rows.append(dict(exponents=x, graph_code=code, graph_edges=sorted(edges),
                                 closed=a >= 2, strong=strong(m+1, arcs), blue=b == c))
        for x in d:
            for y in d:
                if divides(x, y):
                    require(image[x] <= image[y], f'monotonicity k={k}: {x}, {y}')
                    checks['divisibility_pairs'] += 1
        if k == 3:
            p, q = (1, 0, 0), (0, 1, 0)
            require(image[p] <= image[q] and not divides(p, q), 'order-reflection hostile')

    fixed_rows = []
    for k in range(1, 13):
        d = divisors(k)
        qfixed = [x for x in d if swap(x) == x]
        cqfixed = [x for x in d if complement(k, swap(x)) == x]
        forced_empty = [x for x in d if divides(x, complement(k, x))]
        require(len(qfixed) == 2*k, f'Q fixed count k={k}')
        require(len(cqfixed) == (0 if k % 2 else 2), f'CQ fixed count k={k}')
        require(forced_empty == [(a, 0, 0) for a in range(1, k//2+1)], f'forced empty set k={k}')
        fixed_rows.append(dict(k=k, size=len(d), Q_fixed=len(qfixed), CQ_fixed=len(cqfixed),
                               forced_empty=forced_empty))

    # Independent recovery of the prior D_2-minus-central-pair monotone count.
    d8 = [x for x in divisors(2) if x not in ((1, 1, 0), (1, 0, 1))]
    target = list(range(8))
    reverse = lambda w: ((w & 1) << 2) | (w & 2) | ((w & 4) >> 2)
    equivariant = monotone = 0
    for values in permutations(target):
        image = dict(zip(d8, values))
        if not all(image[complement(2, x)] == 7 ^ image[x]
                   and image[swap(x)] == reverse(image[x]) for x in d8):
            continue
        equivariant += 1
        if all(image[x] & image[y] == image[x] for x in d8 for y in d8 if divides(x, y)):
            monotone += 1
    require((equivariant, monotone) == (32, 2), 'inherited D8 counts')
    result = dict(status='FINITE-EXACT; all explicit checks passed, including under -O',
                  universe='D1 and D3 full maps; Dk fixed/witness sets k=1..12; all 8! D8 bijections',
                  checks=dict(checks), k3_table=rows, fixed_and_forced_rows=fixed_rows,
                  inherited_d8_counts=dict(equivariant=equivariant, monotone=monotone))
    rendered = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(rendered, encoding='utf-8')
    print(rendered)


if __name__ == '__main__':
    main()

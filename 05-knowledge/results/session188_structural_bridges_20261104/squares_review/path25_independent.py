#!/usr/bin/env python3
"""Independent endpoint-rooted memoized path count, anchored at forced leaf 18."""
from functools import lru_cache
from collections import Counter
from pathlib import Path
import json


def check(x, why):
    if not x:
        raise RuntimeError(why)


def main():
    n = 25
    square_values = [k*k for k in range(2,8)]
    adjacency = {u:tuple(sorted(s-u for s in square_values if 1 <= s-u <= n and s-u != u))
                 for u in range(1,n+1)}
    check([u for u in adjacency if len(adjacency[u]) == 1] == [18], 'forced leaf')
    all_edges = {(u,v) for u in adjacency for v in adjacency[u] if u < v}
    inventory = Counter(('oo' if u%2 and v%2 else 'ee' if u%2 == v%2 == 0 else 'cross')
                        for u,v in all_edges)
    check(inventory['oo'] == 9 and inventory['ee'] == 6, 'inventory')

    @lru_cache(None)
    def count(current, remaining):
        if remaining == 0:
            return {(current,0,0):1}
        total = Counter()
        for child in adjacency[current]:
            bit = 1 << (child-1)
            if remaining & bit:
                odd = int(current%2 == child%2 == 1)
                even = int(current%2 == child%2 == 0)
                for (end,oo,ee),weight in count(child, remaining ^ bit).items():
                    total[end,oo+odd,ee+even] += weight
        return dict(total)

    rooted = count(18, ((1 << n)-1) ^ (1 << 17))
    check(sum(rooted.values()) == 10, 'rooted count')
    histogram, endpoints = Counter(),Counter()
    for (end,oo,ee),weight in rooted.items():
        check(end != 18, 'same endpoints')
        histogram[oo,ee] += 2*weight
        endpoints[18,end] += weight
        endpoints[end,18] += weight
    published = json.loads((Path(__file__).resolve().parent.parent/'squares'/'path25.json').read_text())
    hist_rows = [[oo,ee,w] for (oo,ee),w in sorted(histogram.items())]
    endpoint_rows = [[u,v,w] for (u,v),w in sorted(endpoints.items())]
    check(hist_rows == published['odd_odd_even_even_histogram'], 'histogram mismatch')
    check(endpoint_rows == published['endpoint_pairs'], 'endpoint mismatch')
    check(all(oo-ee in (1,2) and 1 <= 9-oo <= 4 for oo,ee in histogram), 'parity application')
    print(json.dumps(dict(status='PASS', rooted_paths=10,directed_paths=20,
                          all_paths_accounted_for='Every Hamilton path has endpoint 18; reversing a rooted path gives its unique other orientation.',
                          inventory=dict(inventory),memoized_states=count.cache_info().currsize,
                          odd_odd_even_even_histogram=hist_rows, endpoint_pairs=endpoint_rows),
                     indent=2,sort_keys=True))


if __name__ == '__main__':
    main()

from itertools import combinations, permutations
from collections import Counter

def certificate(n):
    edges = list(combinations(range(n), 2))
    forced = {(i, i + 1): True for i in range(n - 1)}
    forced[(0, n - 1)] = False
    free = [e for e in edges if e not in forced]
    histogram = Counter()
    hs = Counter()
    for bitfield in range(1 << len(free)):
        orient = dict(forced)
        orient.update((e, bool(bitfield & (1 << i))) for i, e in enumerate(free))
        def arc(u, v):
            a, b = sorted((u, v))
            return orient[a, b] == (u == a)
        cycles = []
        for size in range(3, n + 1, 2):
            for vertices in combinations(range(n), size):
                for rest in permutations(vertices[1:]):
                    cycle = (vertices[0],) + rest
                    if all(arc(cycle[i], cycle[(i+1) % size]) for i in range(size)):
                        cycles.append(frozenset(cycle))
        C = len(cycles)
        D = sum(not (a & b) for a, b in combinations(cycles, 2))
        # No three vertex-disjoint cycles fit in n<=6.
        h = 1 + 2*C + 4*D
        # Independent direct counting of directed Hamiltonian paths.
        paths = sum(all(arc(order[i], order[i+1]) for i in range(n-1))
                    for order in permutations(range(n)))
        if h != paths:
            raise RuntimeError(f"identity mismatch: n={n}, bitfield={bitfield}, H={h}, paths={paths}")
        histogram[(C,D)] += 1
        hs[h] += 1
    return sorted(histogram.items()), sorted(hs.items())

for n in range(3, 7):
    pairs, hs = certificate(n)
    print('n=', n, 'H histogram=', hs)
    if n == 6: print('n=6, C<=10:', [(k,v) for k,v in pairs if k[0] <= 10])

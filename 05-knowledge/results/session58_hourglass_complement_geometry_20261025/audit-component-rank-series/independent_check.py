"""Independent adjacent-generator signed-orbit audit; stdlib, n <= 6 only."""
from collections import deque
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import hashlib
import json

OUT = Path(__file__).parent

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

@lru_cache(None)
def pairs(n):
    return tuple(combinations(range(n), 2))

@lru_cache(None)
def adjacent_maps(n):
    ps = pairs(n)
    index = {edge: k for k, edge in enumerate(ps)}
    result = []
    for j in range(n - 1):
        p = list(range(n))
        p[j], p[j + 1] = p[j + 1], p[j]
        destinations = tuple(index[tuple(sorted((p[a], p[b])))] for a, b in ps)
        # An adjacent transposition reverses only the pair (j,j+1).
        result.append((destinations, 1 << index[(j, j + 1)]))
    return tuple(result)

def transport(mask, destinations):
    result = 0
    while mask:
        low = mask & -mask
        result |= 1 << destinations[low.bit_length() - 1]
        mask -= low
    return result

def signed_orbits(n):
    seen = set()
    rows = []
    for representative in range(1 << len(pairs(n))):
        if representative in seen:
            continue
        signs = {representative: 1}
        queue = deque([representative])
        coherent = True
        while queue:
            mask = queue.popleft()
            for destinations, reversed_pair in adjacent_maps(n):
                image = transport(mask, destinations)
                sign = signs[mask] * (-1 if mask & reversed_pair else 1)
                if image not in signs:
                    signs[image] = sign
                    queue.append(image)
                elif signs[image] != sign:
                    coherent = False
        seen.update(signs)
        # Odd stabilizer forces cancellation, regardless of chosen BFS signs.
        alpha = Fraction(sum(signs.values()), len(signs)) if coherent else Fraction(0)
        rows.append((representative, coherent, alpha, len(signs)))
    require(len(seen) == 1 << len(pairs(n)), "orbit coverage")
    return rows

@lru_cache(None)
def direct(n, mask):
    edges = [edge for k, edge in enumerate(pairs(n)) if mask >> k & 1]
    edge_set = frozenset(edges)
    signed_sum = 0
    even = True
    automorphisms = 0
    for p in permutations(range(n)):
        sign = -1 if sum(p[a] > p[b] for a, b in edges) % 2 else 1
        signed_sum += sign
        image = frozenset(tuple(sorted((p[a], p[b]))) for a, b in edges)
        if image == edge_set:
            automorphisms += 1
            even = even and sign == 1
    return even, Fraction(signed_sum, factorial(n)), automorphisms

def components(n, mask):
    adjacency = [set() for _ in range(n)]
    for k, (a, b) in enumerate(pairs(n)):
        if mask >> k & 1:
            adjacency[a].add(b)
            adjacency[b].add(a)
    result = []
    unseen = set(range(n))
    while unseen:
        todo = [min(unseen)]
        vertices = set(todo)
        unseen.remove(todo[0])
        while todo:
            for other in adjacency[todo.pop()]:
                if other in unseen:
                    unseen.remove(other)
                    vertices.add(other)
                    todo.append(other)
        vertices = sorted(vertices)
        submask = sum(1 << k for k, (a, b) in enumerate(pairs(len(vertices)))
                      if vertices[b] in adjacency[vertices[a]])
        result.append((len(vertices), submask))
    return result

def product_series(connected, selector, signed=False, bound=6):
    values = [1] + [0] * bound
    for n, mask, alpha in connected:
        if not selector(n, mask, alpha):
            continue
        parity = -1 if signed and mask.bit_count() % 2 else 1
        old = values[:]
        values = [0] * (bound + 1)
        for base, coefficient in enumerate(old):
            for multiplicity in range((bound - base) // n + 1):
                values[base + multiplicity * n] += coefficient * parity ** multiplicity
    return values

def main():
    census = []
    connected = []
    detailed = []
    for n in range(7):
        rows = signed_orbits(n)
        totals = dict(n=n, graph_types=len(rows), royle_even=0, even_edge_royle_even=0,
                      rank=0, connected_royle_even=0, connected_odd_edge_royle_even=0,
                      connected_nonzero_alpha=0, connected_even_edge_zero_alpha=0)
        for mask, even, alpha, orbit_size in rows:
            direct_even, direct_alpha, aut = direct(n, mask)
            require((even, alpha) == (direct_even, direct_alpha), f"independent signed orbit {n},{mask}")
            comps = components(n, mask)
            component_data = [direct(k, submask) for k, submask in comps]
            require(even == all(row[0] for row in component_data), f"component Royle parity {n},{mask}")
            product = Fraction(1)
            for row in component_data:
                product *= row[1]
            require(alpha == product, f"component alpha factorization {n},{mask}")
            if mask.bit_count() % 2:
                require(alpha == 0, f"odd edge alpha {n},{mask}")
            if even:
                totals['royle_even'] += 1
                totals['even_edge_royle_even'] += mask.bit_count() % 2 == 0
                totals['rank'] += alpha != 0
                if len(comps) == 1:
                    connected.append((n, mask, alpha))
                    totals['connected_royle_even'] += 1
                    totals['connected_odd_edge_royle_even'] += mask.bit_count() % 2
                    totals['connected_nonzero_alpha'] += alpha != 0
                    totals['connected_even_edge_zero_alpha'] += mask.bit_count() % 2 == 0 and alpha == 0
                detailed.append(dict(n=n, mask=mask, edge_count=mask.bit_count(), connected=len(comps)==1,
                                     alpha=str(alpha), orbit_size=orbit_size, automorphisms=aut))
        census.append(totals)
    T = product_series(connected, lambda *args: True)
    D = product_series(connected, lambda *args: True, signed=True)
    R = product_series(connected, lambda n, mask, alpha: alpha != 0)
    C = [(t + d) // 2 for t, d in zip(T, D)]
    require(all((t + d) % 2 == 0 for t, d in zip(T, D)), 'parity projector integrality')
    require(T == [r['royle_even'] for r in census], 'T series')
    require(C == [r['even_edge_royle_even'] for r in census], 'C series')
    require(R == [r['rank'] for r in census], 'R series')
    require(C == [1,1,1,2,3,10,34] and R == [1,1,1,2,3,10,32], 'frozen dimension agreement')
    require(direct(6,255) == (True,Fraction(0),6), 'F3 witness')
    require(direct(6,2015) == (True,Fraction(0),2), 'second connected zero witness')
    require(len(components(6,255)) == len(components(6,2015)) == 1, 'connected zero witnesses')
    require((255).bit_count() == 8 and (2015).bit_count() == 10, 'distinct even-edge witnesses')
    result = {'status':'PASS', 'method':'adjacent-generator signed-orbit BFS plus independent permutation signs',
              'census':census, 'series':dict(T=T,D=D,C=C,R=R,K=[c-r for c,r in zip(C,R)]),
              'royle_even_types':detailed}
    (OUT/'independent_certificate.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status':result['status'],'census':census,'series':result['series']}, indent=2))

if __name__ == '__main__':
    main()
